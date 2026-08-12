// Rate-parity test for BNG2's symmetry_factor on an MM(kcat,Km) rate law
// (issue #37).
//
// `compute_propensity` applied the factor on every rate law it handles
// except the Michaelis-Menten branch, which returned
// `kcat*sFree*E/(Km+sFree)` off the raw substrate match count.  Nothing in
// any of the three corpora paired a non-unit symmetry_factor with MM, so
// the omission was invisible: a symmetric-reaction-center MM rule ran fast,
// by a factor that depends on where in the saturation curve the rule sits.
//
// The correction belongs on a reactant COUNT, inside the law:
//
//     S = a_used * symmetry_factor        (canonical shape)
//     a = kcat * sFree * E / (Km + sFree),  sFree from (S, E, Km)
//
// What the factor corrects is a match multiplicity.  A pattern with a
// non-trivial reaction-center automorphism — A(d!1).A(d!1), where either
// molecule can seed the match — matches each complex twice, so the count is
// 2x the number of complexes and the law is handed more of that reactant
// than exists.  MM is not linear in that count, so scaling the finished
// propensity instead is exact only below saturation, while dropping the
// factor is exact only deep inside it.  Scaling the count is exact
// everywhere and is what BNG2's network expansion integrates: it builds one
// reaction per rule application and folds the factor into the reaction
// multiplicity (2 maps x 0.5 = 1), leaving the law reading species counts.
// bngsim's patched NFsim (lanl/bngsim#195 -> #278) does the same.
//
// Which count is the one the rule TRANSFORMS.  BNG2 refuses MM on a rule
// whose two reactant patterns are isomorphic, so no automorphism carries one
// pattern onto the other and a pure-context pattern contributes 1 (all its
// automorphisms stabilize the empty reaction center) — it is also already
// counted per complex, issue #33.  For the canonical shape, where the enzyme
// is a catalyst and so context, that slot is the substrate.
//
// Model: tests/cpp/mm_symmetry_model.bngl, four arms, each sized so a 2 s
// window depletes its own pool by ~1-10%:
//
//   RA  symmetric dimer, saturated (S=200 dimers, Km=1)
//   RB  symmetric dimer, linear    (S=1000 dimers, Km=10000)
//   RC  asymmetric control, symmetry_factor="1", saturated
//   RD  symmetry in the enzyme slot: the rule transforms its second reactant
//       (a symmetric dimer) and its substrate is context.  Off-spec for MM,
//       but BNG2 writes it and puts symmetry_factor="0.5" there.
//
// No single arm separates all the candidate propensities, which is why there
// are four.  Expected firings over t in [0,2], from the mean-field
// integration below:
//
//                            RA (sat)  RB (linear)  RC (ctrl)  RD (enzyme)
//   factor on transformed      9.95       9.05        9.95       19.01  <- correct
//   no factor at all           9.97      16.55        9.95       36.20  <- pre-fix
//   factor on propensity       4.99       8.30        9.95       19.00
//   factor always on substrate 9.95       9.05        9.95       36.02
//
// Firings are read off the observables: each dimer-splitting firing frees
// two monomers (Afree/2, Bfree/2, Gfree/2), each control firing makes one
// C(x~P).
//
// The oracle is the mean-field integral rather than a closed form because
// each arm's pool does deplete slightly over the window; the residual error
// against the true stochastic mean is second order in that depletion (the
// law is concave, and Var(N) over the window is ~1e-2 of N), far inside the
// tolerance.  It is not self-referential: BNG2's own network -> ODE run of
// the RD arm alone gives Gfree(2) = 38.02, i.e. 19.01 firings, the number
// below.  The ensemble comparison against BNG2's ODE lives in
// tests/models/feature_coverage/ft_mm_ratelaw_sym.bngl.

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

double final_value(const rulemonkey::Result& r, const std::string& name) {
  int const i = idx_of(r, name);
  if (i < 0)
    throw std::runtime_error("observable not found: " + name);
  return r.observable_data[i].back();
}

// NFsim's MMRxnClass::update_a, on COMPLEX counts: the semantics the engine
// is expected to implement once symmetry_factor has deflated the transformed
// slot's match count back to a complex count.
double mm_propensity(double S, double E, double kcat, double Km) {
  double const diff = S - Km - E;
  double const sFree = 0.5 * (diff + std::sqrt((diff * diff) + (4.0 * Km * S)));
  return kcat * sFree * E / (Km + sFree);
}

// One arm: how its firings are read off the run, and the parameters that
// determine the expected count.  Counts are COMPLEXES — what the law reads
// once symmetry_factor has deflated the transformed slot's match count.
struct Arm {
  std::string name;
  std::string observable;
  double per_firing;    // observable units produced per firing
  double s0;            // substrate complexes at t=0
  double e0;            // enzyme complexes at t=0
  bool enzyme_depletes; // true iff the rule transforms its enzyme slot
  double kcat;
  double km;
};

// Mean-field firings over [0, T] for one arm, counting each consumed complex
// as one firing.  Whichever slot the rule transforms is the one that drains;
// the other is context and holds.  Step size is far below the arm's timescale
// (the fastest pool turns over at ~5e-2 /s), so the quadrature error is
// negligible next to the tolerance.
double expected_firings(const Arm& arm, double T) {
  constexpr int kSteps = 20000;
  double const dt = T / kSteps;
  double s = arm.s0;
  double e = arm.e0;
  double fired = 0.0;
  for (int i = 0; i < kSteps; ++i) {
    double const a = mm_propensity(s, e, arm.kcat, arm.km);
    fired += a * dt;
    (arm.enzyme_depletes ? e : s) -= a * dt;
  }
  return fired;
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <mm_symmetry_model.xml>\n", argv[0]);
    return 2;
  }
  std::string const xml = argv[1];

  constexpr double kTEnd = 2.0;
  constexpr int kReps = 200;

  // Mirrors the parameter and seed-species blocks of mm_symmetry_model.bngl.
  std::vector<Arm> const arms = {
      {"RA (symmetric, saturated)", "Afree", 2.0, 200.0, 5.0, false, 1.0, 1.0},
      {"RB (symmetric, linear)", "Bfree", 2.0, 1000.0, 5.0, false, 10.0, 10000.0},
      {"RC (asymmetric control)", "Cprod", 1.0, 200.0, 5.0, false, 1.0, 1.0},
      {"RD (symmetry in enzyme slot)", "Gfree", 2.0, 1000.0, 200.0, true, 0.05, 1.0},
  };

  rulemonkey::TimeSpec ts;
  ts.t_start = 0.0;
  ts.t_end = kTEnd;
  ts.n_points = 2;

  std::vector<double> sum(arms.size(), 0.0);
  std::vector<double> sumsq(arms.size(), 0.0);

  try {
    rulemonkey::RuleMonkeySimulator sim(xml);
    for (int rep = 0; rep < kReps; ++rep) {
      rulemonkey::Result const r =
          sim.run(ts, std::uint64_t{1000} + static_cast<std::uint64_t>(rep));
      for (size_t i = 0; i < arms.size(); ++i) {
        double const fired = final_value(r, arms[i].observable) / arms[i].per_firing;
        sum[i] += fired;
        sumsq[i] += fired * fired;
      }
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: simulation threw: %s\n", e.what());
    return 1;
  }

  for (size_t i = 0; i < arms.size(); ++i) {
    const Arm& arm = arms[i];
    double const mean = sum[i] / kReps;
    double const var = (sumsq[i] / kReps) - (mean * mean);
    double const se = std::sqrt((var > 0 ? var : 0.0) / kReps);
    double const want = expected_firings(arm, kTEnd);

    // 4 standard errors of the ensemble mean, plus 2% of the expectation for
    // the mean-field-vs-stochastic residue.  ~11% of the expected count with
    // 200 reps, against contrasts of 83% (RB, factor dropped), 100% (RA,
    // factor on the propensity) and 90% (RD, factor attributed to the
    // substrate rather than to the slot the rule transforms).
    double const tol = (4.0 * se) + (0.02 * want);
    std::printf("%-28s fired=%6.3f  expected=%6.3f  tol=%5.3f  (se=%.3f)\n", arm.name.c_str(), mean,
                want, tol, se);
    check(std::fabs(mean - want) <= tol, arm.name + ": mean firings " + std::to_string(mean) +
                                             " outside " + std::to_string(want) + " +/- " +
                                             std::to_string(tol));
  }

  if (g_failures) {
    std::fprintf(stderr, "%d check(s) failed\n", g_failures);
    return 1;
  }
  std::printf("mm_symmetry_test: all checks passed\n");
  return 0;
}
