// Observable scope inside a local function (issue #38).
//
// An observable applied to the function's own tag — `Mod(x)` — is evaluated
// at the tagged molecule.  One written bare — `Vol` — is the system-wide
// count.  BNG2's `<ListOfReferences>` marks both `type="Observable"` with
// nothing to separate them, so RM used to sweep every one of them into the
// local set and evaluate it at the tagged molecule.  A bare observable is
// almost never present in the tagged molecule's own complex, so it read 0,
// the propensity went to 0, and the rule never fired — no warning, mass
// conserved, trajectory flat.  The classification now comes from parsing
// the `<Expression>`, which is the only place the distinction survives.
//
// Three things are pinned here, on two models:
//
//   ft_local_fcn_global_obs   — the base fix, on both propensity paths.
//       RU  uflip(x) = kf*Vol*Uphos(x)   unimolecular
//       RB  ecat(x)  = kc*Vol*Mod(x)     bimolecular DOR1
//     Each arm carries a subpopulation whose TAGGED observable is 0 —
//     U(p~0) and E(e~0) — which must never react.  That separates the two
//     possible scope errors: localizing the bare observable stalls both
//     rules outright, while globalizing the tagged one would let the inert
//     subpopulations react and run everything orders of magnitude fast.
//     Nothing consumes V, U(p~*) or E, so `Vol` stays 50 and each arm is an
//     exact binomial death process at c = 0.02:
//         Uoff(t) = 200 + 400·exp(-c·t)      (only the U(p~1) half reacts)
//         Sub(t)  = 600·exp(-c·t)
//
//   ft_local_fcn_mixed_scope  — the two cases the base model cannot reach.
//       lfA(x) = ka*Aoff(x)         Aoff at local scope
//       lfB(x) = kb*Aoff*Boff(x)    the same Aoff at global scope
//     (a) RM holds one eval slot per observable, so the union of the
//         per-function local sets cannot stand in for them: the localized
//         `Aoff` would leak into lfB and stall RB at 200.
//     (b) `Aoff` decays as RA fires, so RB's rate moves for every B
//         instance while no B molecule is marked affected.  Holding the
//         seed value instead gives Boff(100) = 200·exp(-kb·200·100) = 27.1
//         against the true 56.5 — which is what the pinned NFsim does.
//     Analytic, with Aoff(t) = 200·exp(-ka·t):
//         Aoff(100) = 200·exp(-1)                      = 73.58
//         Boff(100) = 200·exp(-2·(1-exp(-1)))          = 56.49
//     matching BNG2's network -> ODE to the digits printed here.
//
//   local_fcn_mixed_ref_model.xml — one observable written BOTH ways in a
//     single function, `kc*(Aoff+Aoff(x))`.  One slot cannot hold two
//     values; RM resolves it at local scope and must say so rather than
//     mis-evaluate it quietly.
//
// Tolerances come from the replicates' own spread (4 SE of the sample
// mean), so no arm depends on a hand-derived variance.

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

int idx_of(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return static_cast<int>(i);
  return -1;
}

const std::vector<double>& series(const rulemonkey::Result& r, const std::string& name) {
  int const i = idx_of(r, name);
  if (i < 0)
    throw std::runtime_error("observable not found: " + name);
  return r.observable_data[i];
}

double final_value(const rulemonkey::Result& r, const std::string& name) {
  return series(r, name).back();
}

constexpr int kReps = 40;
constexpr double kSigmaTol = 4.0;

// Compare a set of replicate terminal values against an analytic mean,
// with the tolerance taken from the replicates' own standard error.
void check_mean(const std::string& label, const std::vector<double>& vals, double analytic) {
  auto const n = static_cast<double>(vals.size());
  double mean = 0.0;
  for (double const v : vals)
    mean += v;
  mean /= n;
  double var = 0.0;
  for (double const v : vals)
    var += (v - mean) * (v - mean);
  var /= (n - 1.0);
  double const se = std::sqrt(var / n);
  double const diff = std::fabs(mean - analytic);

  std::fprintf(stderr, "%-26s empirical=%8.3f analytic=%8.3f |diff|=%6.3f (%.2f SE)\n",
               label.c_str(), mean, analytic, diff, se > 0 ? diff / se : 0.0);
  check(diff < kSigmaTol * se, label + " terminal mean within 4 SE of the analytic value");
}

// ---------------------------------------------------------------------------
// ft_local_fcn_global_obs: a bare observable multiplying a tagged one, on the
// unimolecular and the bimolecular (DOR1) propensity paths.
// ---------------------------------------------------------------------------
void test_global_obs_arms(const std::string& xml) {
  constexpr double kTEnd = 60.0;
  double const p = std::exp(-0.02 * kTEnd);

  std::vector<double> uoff, sub;
  int stalled = 0;
  int inert_reacted = 0;

  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, kTEnd, 2},
                     /*seed=*/std::uint64_t{3800} + static_cast<std::uint64_t>(rep));

    // Context only — no rule touches V, the `p` site, or the enzymes, so
    // every arm's hazard is constant and the binomial oracle holds.  `Vol`
    // is also the observable under test: if it were being localized, this
    // globally-read value would still be 50 while the rate saw 0.
    check(final_value(r, "Vol") == 50.0, "Vol conserved at 50");
    check(final_value(r, "Uphos") == 400.0, "Uphos conserved at 400");
    check(final_value(r, "Mod") == 20.0, "Mod conserved at 20");

    double const u = final_value(r, "Uoff");
    double const s = final_value(r, "Sub");
    uoff.push_back(u);
    sub.push_back(s);

    // Pre-fix signature: `Vol` localized to 0, both rules inert, both
    // observables flat at their seed values of 600.
    if (u > 595.0 || s > 595.0)
      ++stalled;
    // Opposite error: `Uphos` read globally as 400 would give every U —
    // including the 200 in state p~0 — a rate of 8.0/s, emptying the pool.
    // Only 400 of the 600 U molecules can ever flip, so Uoff has a hard
    // floor at 200 that no amount of sampling luck can cross.
    if (u < 200.0)
      ++inert_reacted;
  }

  check(stalled == 0, "neither arm stalled at its seed value");
  check(inert_reacted == 0, "U(p~0) never reacts — the tagged observable stayed local");

  check_mean("Uoff (unimolecular)", uoff, 200.0 + (400.0 * p));
  check_mean("Sub  (bimolecular DOR1)", sub, 600.0 * p);
}

// ---------------------------------------------------------------------------
// ft_local_fcn_mixed_scope: one observable local in one function and global in
// another, and a global whose value moves during the run.
// ---------------------------------------------------------------------------
void test_mixed_scope(const std::string& xml) {
  constexpr double kTEnd = 100.0;
  double const a_final = 200.0 * std::exp(-0.01 * kTEnd);
  // Boff(t) = 200·exp(-kb·Integral(Aoff)), Integral = 20000·(1-exp(-ka·t)).
  double const b_final = 200.0 * std::exp(-1.0e-4 * 20000.0 * (1.0 - std::exp(-0.01 * kTEnd)));

  std::vector<double> aoff, boff;
  int stalled = 0;

  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, kTEnd, 2},
                     /*seed=*/std::uint64_t{3810} + static_cast<std::uint64_t>(rep));
    double const b = final_value(r, "Boff");
    aoff.push_back(final_value(r, "Aoff"));
    boff.push_back(b);
    // Union scoping would hand lfB the localized Aoff — 0 at a B molecule —
    // and RB would never fire at all.
    if (b > 195.0)
      ++stalled;
  }

  check(stalled == 0, "RB fired — lfB read Aoff at global scope, not lfA's local scope");
  check_mean("Aoff (tagged in lfA)", aoff, a_final);
  check_mean("Boff (bare in lfB)", boff, b_final);

  // Holding the bare observable at its seed value — the pinned NFsim's
  // behaviour — lands at 27.1 rather than 56.5.  Assert the separation
  // directly so a regression there reads as itself and not as noise.
  double mean_b = 0.0;
  for (double const v : boff)
    mean_b += v;
  mean_b /= static_cast<double>(boff.size());
  double const stale = 200.0 * std::exp(-1.0e-4 * 200.0 * kTEnd);
  std::fprintf(stderr, "Boff mean %.3f vs live-global %.3f vs frozen-global %.3f\n", mean_b,
               b_final, stale);
  check(std::fabs(mean_b - b_final) < std::fabs(mean_b - stale),
        "Boff tracks the live global observable, not its seed value");
}

// ---------------------------------------------------------------------------
// One observable, both scopes, one function: `kc*(Aoff+Aoff(x))`.  RM has a
// single eval slot per observable and cannot hold both readings, so it takes
// the tagged one and reports the model rather than mis-evaluating in silence.
// ---------------------------------------------------------------------------
void test_both_scopes_in_one_function_warns(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  bool warned = false;
  for (const auto& w : sim.unsupported_features()) {
    if (w.severity == rulemonkey::Severity::Warn && w.feature.find("'Aoff'") != std::string::npos &&
        w.feature.find("bare") != std::string::npos)
      warned = true;
    check(w.severity != rulemonkey::Severity::Error,
          "mixed-scope reference is a warning, not a refusal");
  }
  check(warned, "an observable used at both scopes in one function is reported");

  // And the model still runs, at the tagged reading: kc·(1+1) = 2e-3 per
  // molecule, so the pool depletes rather than stalling or exploding.
  auto r = sim.run({0.0, 500.0, 2}, /*seed=*/17);
  double const left = final_value(r, "Aoff");
  std::fprintf(stderr, "both-scopes fixture: Aoff 100 -> %.0f (expected ~37)\n", left);
  check(left > 10.0 && left < 70.0, "both-scopes fixture runs at the tagged reading");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::fprintf(stderr, "Usage: local_fcn_global_obs_test <ft_local_fcn_global_obs.xml> "
                         "<ft_local_fcn_mixed_scope.xml> <local_fcn_mixed_ref_model.xml>\n");
    return 2;
  }
  try {
    test_global_obs_arms(argv[1]);
    test_mixed_scope(argv[2]);
    test_both_scopes_in_one_function_warns(argv[3]);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }
  if (g_failures != 0) {
    std::fprintf(stderr, "%d failure(s)\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "PASS\n");
  return 0;
}
