// Per-molecule rule tables split by rule shape (issue #71).
//
// `RuleState` used to hold one 80-byte row per molecule per rule, carrying
// every field any rule shape could want.  It now holds a 12-byte row every
// rule carries and a 64-byte row only the shapes that read one allocate, in
// a second table (`mol_aux`) kept exactly as long as the first.  Three
// things can now go wrong that could not before:
//
//   * a shape that reads the wide half is not given one, or is given one
//     that has not grown with the narrow half — the rule then prices a row
//     that is missing, or that belongs to another molecule;
//   * the shape predicate is read as false for a rule that does need the
//     wide half, so its per-molecule rates and its shared-component split
//     are silently zero.
//
// (The third change in #71, folding the `cache_init` flag into the value
// the rescan fills the table with instead of a second pass over it, is not
// something a test can reach: `AddMolecule` marks every bit of a new
// molecule changed, so a row that came up flagged would be recomputed
// anyway.  The on-demand growth path still takes the unflagged default,
// which is what makes that true rather than merely lucky.)
//
// `rule_table_shape_model.bngl` puts one rule of each wide-half group and
// one of the no-wide-half shape in a single pool, and makes every reacting
// molecule with a maker rule rather than seeding it, so each arm runs on
// rows created after the session build sized the tables.  See that file for
// the arm-by-arm layout.
//
// The exact assertions are structural: the makers fire exactly once per
// maker molecule, A is conserved at 8 and paired two at a time, the E pair
// is one complex, and the three decay arms run to zero over a horizon 200
// times their mean lifetime.  The four means are the statistical ones, and
// each is there to PRICE its arm rather than merely observe that it fired:
// a rule reading a missing or misaligned wide-half row still runs, just at
// the wrong rate.  Each is a 4-sigma band on the standard error of the mean
// over 200 replicates, against an analytic reference — the CME steady state
// for the homodimer, and the plain exponential survivor count for the three
// decay arms.
//
// The catalysis arm's band is also what says the pure-context slot is
// engaged at all: its E arrives as a bonded pair, so a slot priced per
// complex runs S at kcat and a slot priced per molecule would run it at
// 2*kcat — 14.7 survivors against 5.4, some 40 sigma apart.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

double value_at(const rulemonkey::Result& r, const std::string& name, std::size_t row) {
  for (std::size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name) {
      if (row >= r.observable_data[i].size())
        throw std::runtime_error("sample row out of range for observable " + name);
      return r.observable_data[i][row];
    }
  throw std::runtime_error("observable not found: " + name);
}

// Mean and variance of the dimer count k at the CME steady state of
// A(a) + A(a) <-> A(a!1).A(a!1), from detailed balance over k:
//
//     pi(k+1)/pi(k) = kp*(A-2k)*(A-2k-1) / (2*km*(k+1))
//
// The same recurrence homodimer_rate_test uses; repeated here rather than
// shared because the two tests pin different things with it.
struct Moments {
  double mean = 0.0;
  double var = 0.0;
};

Moments cme_dimer_moments(int a_tot, double kp, double km) {
  int const kmax = a_tot / 2;
  std::vector<double> pi(kmax + 1, 0.0);
  pi[0] = 1.0;
  for (int k = 0; k < kmax; ++k) {
    int const free_a = a_tot - (2 * k);
    double const fwd = kp * free_a * (free_a - 1) / 2.0;
    double const rev = km * (k + 1);
    pi[k + 1] = pi[k] * fwd / rev;
  }
  double z = 0.0;
  for (double const p : pi)
    z += p;
  Moments m;
  for (int k = 0; k <= kmax; ++k)
    m.mean += k * pi[k] / z;
  for (int k = 0; k <= kmax; ++k)
    m.var += (k - m.mean) * (k - m.mean) * pi[k] / z;
  return m;
}

// Survivors of `n` independent molecules each decaying at unit rate, at
// time `t`: Binomial(n, exp(-t)).
Moments survivor_moments(int n, double t) {
  double const p = std::exp(-t);
  Moments m;
  m.mean = n * p;
  m.var = n * p * (1.0 - p);
  return m;
}

// Mean over replicates of one observable at one sample row.
struct Mean {
  double sum = 0.0;
  int n = 0;
  void add(double v) {
    sum += v;
    ++n;
  }
  double get() const { return (n > 0) ? sum / n : 0.0; }
};

constexpr int kReps = 200;
constexpr double kEarly = 1.0; // sample that prices the three decay arms
constexpr double kEnd = 200.0; // sample that reads the homodimer equilibrium
constexpr int kATot = 8;       // molecules the MkA rule makes
constexpr int kLTot = 20;      // ... MkL
constexpr int kDTot = 30;      // ... MkD
constexpr int kSTot = 40;      // S seeded, converted by the catalysis arm
constexpr double kKp = 0.5;
constexpr double kKm = 0.5;

// Price one arm against its analytic reference: report the comparison,
// then band it at 4 sigma on the standard error of the mean over kReps
// replicates, plus whatever absolute `pad` the call site's own bias needs.
void price_arm(const char* what, double empirical, const Moments& ref, double pad,
               const std::string& msg) {
  double const se = std::sqrt(ref.var / kReps);
  std::fprintf(stderr, "  %-28s empirical=%7.4f  analytic=%7.4f  |diff|=%6.4f  se=%6.4f\n", what,
               empirical, ref.mean, std::fabs(empirical - ref.mean), se);
  check(std::fabs(empirical - ref.mean) < (4.0 * se) + pad, msg);
}

void test_every_rule_shape_prices_rows_born_after_the_build(const std::string& xml) {
  Mean lu_early, su_early, dn_early, adim_end;

  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    rulemonkey::TimeSpec ts;
    ts.t_start = 0.0;
    ts.t_end = kEnd;
    ts.sample_times = {kEarly, kEnd};
    auto const r = sim.run(ts, /*seed=*/std::uint64_t{7000} + static_cast<std::uint64_t>(rep));

    for (std::size_t row = 0; row < 2; ++row) {
      // Every maker fires exactly once per maker molecule, so these are
      // exact counts of how many molecules each arm was given — and the
      // makers are themselves the shape that carries no wide half.
      check(value_at(r, "MkAd", row) == kATot && value_at(r, "MkLd", row) == kLTot &&
                value_at(r, "MkEd", row) == 1.0 && value_at(r, "MkDd", row) == kDTot,
            "rep " + std::to_string(rep) +
                ": every maker rule should have fired once per maker "
                "molecule by the first sample");
      // Nothing creates or destroys A once the maker is done, and the
      // homodimer rules move it two molecules at a time.
      check(value_at(r, "Atot", row) == kATot,
            "rep " + std::to_string(rep) + ": A is conserved at " + std::to_string(kATot));
      check(value_at(r, "Amon", row) + value_at(r, "Adim", row) == kATot,
            "rep " + std::to_string(rep) + ": every A is either free or in a dimer");
      check(std::fmod(value_at(r, "Adim", row), 2.0) == 0.0,
            "rep " + std::to_string(rep) + ": bonded A come in pairs");
      // The catalyst arrives as a bonded pair: two molecules, one complex.
      // The context slot prices the complex, which is what makes the S
      // band below able to tell the two readings apart.
      check(value_at(r, "Emol", row) == 2.0 && value_at(r, "Ecx", row) == 1.0,
            "rep " + std::to_string(rep) +
                ": the catalyst should be two E molecules in one "
                "complex");
      // The local-rate arm changes L's state, never its population.
      check(value_at(r, "Lall", row) == kLTot,
            "rep " + std::to_string(rep) + ": L is conserved at " + std::to_string(kLTot));
    }

    // Each decay arm has a mean lifetime of 1, so a horizon of 200 leaves
    // nothing.  A row the rule cannot price is a molecule that never
    // reacts, which is exactly what these would catch.
    check(value_at(r, "Lu", 1) == 0.0,
          "rep " + std::to_string(rep) +
              ": every L should have flipped by t=" + std::to_string(kEnd));
    check(value_at(r, "Su", 1) == 0.0,
          "rep " + std::to_string(rep) +
              ": every S should have been converted by t=" + std::to_string(kEnd));
    check(value_at(r, "Dn", 1) == 0.0,
          "rep " + std::to_string(rep) +
              ": every D should have decayed by t=" + std::to_string(kEnd));

    lu_early.add(value_at(r, "Lu", 0));
    su_early.add(value_at(r, "Su", 0));
    dn_early.add(value_at(r, "Dn", 0));
    adim_end.add(value_at(r, "Adim", 1));

    sim.destroy_session();
  }

  std::fprintf(stderr, "rule-shape arms over %d replicates:\n", kReps);

  // Local rate law: the per-molecule rate is kdec * (L in the tagged
  // molecule's own complex) = kdec, so survivors are Binomial(nL, e^-t).
  price_arm("local-rate L(t~u) at t=1", lu_early.get(), survivor_moments(kLTot, kEarly), 0.0,
            "the local-rate arm should run at kdec per molecule");
  // Pure-context slot: one E complex, so kcat per S.  A per-molecule
  // reading of the same pair would give 2*kcat and land near 5.4.
  price_arm("context-slot S(x~u) at t=1", su_early.get(), survivor_moments(kSTot, kEarly), 0.0,
            "the catalysis arm should run at kcat per S, i.e. per E COMPLEX");
  // No wide half at all — the shape issue #71 is about.
  price_arm("plain D() -> 0 at t=1", dn_early.get(), survivor_moments(kDTot, kEarly), 0.0,
            "the plain unimolecular arm should run at kdel per molecule");
  // Homodimer: the shared-component split and the ab_both_sq deflation.
  // The observable counts bonded A molecules, i.e. twice the dimer count.
  // t_end is 100 relaxation times but not infinite, so the band carries
  // the same 0.05 residual pad homodimer_rate_test uses.
  {
    Moments const k = cme_dimer_moments(kATot, kKp, kKm);
    price_arm("homodimer A(a!+) at t=200", adim_end.get(), Moments{2.0 * k.mean, 4.0 * k.var}, 0.05,
              "the homodimer arm should sit at its CME steady state");
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <rule_table_shape_model.xml>\n", argv[0]);
    return 2;
  }
  try {
    test_every_rule_shape_prices_rows_born_after_the_build(argv[1]);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: all four rule shapes price molecules born after the session build "
                       "sized their tables\n");
  return 0;
}
