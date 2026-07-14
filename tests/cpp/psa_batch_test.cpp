// Partial-scaling (opt-in) Phase-1 + Phase-2 behavioral test.
//
// Exercises the scaled batch-fire path that `set_critical_population(Nc)`
// enables (Lin, Feng & Hlavacek 2019).  The exact/large-Nc byte-identity
// contract is pinned in api_coverage_test; here we drive the *scaled* path and
// assert the properties that make it correct and useful:
//
//   1. Batching actually happens (ps_reaction_firings exceeds the SSA step
//      count) for both the unimolecular (Phase 1) and the bimolecular
//      association + dissociation (Phase 2) paths.
//   2. First moment stays unbiased: the ensemble mean of the scaled run agrees
//      with the analytic mean (Poisson for birth_death, binding equilibrium
//      for the reversible dimer) and with the exact-SSA ensemble.
//   3. Per-rule multiplier telemetry is sane: a batched rule has a time-
//      averaged multiplier > 1, a zeroth-order rule stays at 1.
//   4. The Nc-too-small guard throws on a model whose only channel is the
//      scaled rule when Nc is aggressive enough to collapse the run.
//
// Under Debug / ASan builds the partial-scaling batch self-check
// (kBatchInvariant, plan §5) runs inside every scaled batch here, so this test
// also exercises the epoch-cache reconciliation invariant in CI — for the
// bimolecular association + dissociation batches too.
//
// argv[1] = birth_death.xml   (0 -> A ; A -> 0, k_death=1, A0=100 => Poisson(100))
// argv[2] = psa_pure_death.xml (A -> 0 only; no source)
// argv[3] = binding.xml        (A(b)+B(a)<->A(b!1).B(a!1); equilibrium Abound=100)

#include "rulemonkey/simulator.hpp"

#include <algorithm>
#include <cmath>
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

// Final value of the (single) observable after one run.
double final_obs(rulemonkey::Result& r, int obs = 0) { return r.observable_data.at(obs).back(); }

// birth_death: A ~ Poisson(100).  Analytic mean = variance = 100.
void test_birth_death_scaled(const std::string& xml) {
  constexpr int kNc = 20;   // K ~ floor(N/20) ~ 5 at steady state
  constexpr int kReps = 80; // SE(mean) ~ std/sqrt(reps); std ~ 17 scaled => ~1.9

  // --- scaled ensemble ---
  double sum = 0.0;
  int64_t total_firings = 0, total_steps = 0;
  std::vector<double> mult; // last run's per-rule multipliers
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_critical_population(kNc);
    auto r = sim.run({0.0, 10.0, 100}, /*seed=*/1000 + s);
    sum += final_obs(r);
    total_firings += r.ps_reaction_firings;
    total_steps += r.event_count;
    mult = r.ps_multipliers;
  }
  double const scaled_mean = sum / kReps;

  // 1. Batching happened: more firings than SSA steps.
  check(total_firings > total_steps,
        "scaled run must fire more reactions than SSA steps (batching active)");

  // 2. First moment unbiased vs the analytic Poisson mean (100).  Loose ±8
  //    tolerance keeps this deterministic across platforms while still
  //    catching a real bias (a mis-scaled batch drifts by tens).
  check(std::fabs(scaled_mean - 100.0) < 8.0,
        "scaled ensemble mean must stay near the analytic Poisson mean 100 (got " +
            std::to_string(scaled_mean) + ")");

  // 3. Telemetry: exactly the unimolecular death rule is scaled.  birth is
  //    zeroth-order (multiplier 1); death carries the batch multiplier (>1).
  check(mult.size() == 2, "birth_death should expose 2 per-rule multipliers");
  if (mult.size() == 2) {
    double const lo = std::min(mult[0], mult[1]);
    double const hi = std::max(mult[0], mult[1]);
    check(std::fabs(lo - 1.0) < 1e-9,
          "the zeroth-order (birth) rule multiplier must be exactly 1 (got " + std::to_string(lo) +
              ")");
    check(hi > 1.5,
          "the batched (death) rule multiplier must exceed 1 (got " + std::to_string(hi) + ")");
  }

  // 4. Cross-check the scaled mean against an exact-SSA ensemble on the same
  //    model/seeds — they must agree (both ~100).
  double exact_sum = 0.0;
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, 10.0, 100}, /*seed=*/1000 + s);
    exact_sum += final_obs(r);
  }
  double const exact_mean = exact_sum / kReps;
  check(std::fabs(scaled_mean - exact_mean) < 8.0,
        "scaled and exact ensemble means must agree (scaled=" + std::to_string(scaled_mean) +
            " exact=" + std::to_string(exact_mean) + ")");

  // An exact (Nc unset) run must not populate the scaled telemetry.
  rulemonkey::RuleMonkeySimulator exact(xml);
  auto er = exact.run({0.0, 10.0, 100}, /*seed=*/7);
  check(er.ps_multipliers.empty(), "an exact run must not expose ps_multipliers");
  check(er.ps_reaction_firings == er.event_count,
        "an exact run's firing count must equal its step count");
}

// psa_pure_death: A -> 0 as the only channel.  With an aggressive Nc the scaled
// total propensity collapses and the guard must refuse to run.
void test_nc_too_small_guard(const std::string& xml) {
  // Exact run is fine.
  {
    rulemonkey::RuleMonkeySimulator sim(xml);
    bool threw = false;
    try {
      sim.run({0.0, 10.0, 10}, /*seed=*/1);
    } catch (const std::exception&) {
      threw = true;
    }
    check(!threw, "exact run of the pure-death model must not throw");
  }
  // Aggressive Nc trips the guard.
  {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_critical_population(1);
    bool threw = false;
    try {
      sim.run({0.0, 10.0, 10}, /*seed=*/1);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    check(threw, "Nc=1 on the pure-death model must trip the Nc-too-small guard");
  }
  // A large Nc leaves K=1 everywhere: no scaling, no guard, a normal run.
  {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_critical_population(1'000'000);
    bool threw = false;
    try {
      sim.run({0.0, 10.0, 10}, /*seed=*/1);
    } catch (const std::exception&) {
      threw = true;
    }
    check(!threw, "a large Nc must not trip the guard (K=1 everywhere => exact)");
  }
}

// binding: A(b) + B(a) <-> A(b!1).B(a!1).  Reversible heterodimer at its
// analytic equilibrium (Abound = 100), so the ensemble is stationary and any
// first-moment bias in the scaled bimolecular association / dissociation batch
// shows up as a drift of the mean.  Exercises BOTH new Phase-2 batch paths in
// one model: `bind` (bimolecular association, AddBond) and `unbind`
// (dissociation, DeleteBond).
void test_binding_scaled(const std::string& xml) {
  constexpr int kNc = 20; // K ~ floor(100/20) = 5 at equilibrium
  constexpr int kReps = 80;
  constexpr int kAboundObs = 1; // observables: {Afree, Abound}

  double sum = 0.0;
  int64_t total_firings = 0, total_steps = 0;
  std::vector<double> mult;
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_critical_population(kNc);
    auto r = sim.run({0.0, 10.0, 100}, /*seed=*/2000 + s);
    sum += final_obs(r, kAboundObs);
    total_firings += r.ps_reaction_firings;
    total_steps += r.event_count;
    mult = r.ps_multipliers;
  }
  double const scaled_mean = sum / kReps;

  // 1. Batching happened: more reactions fired than SSA steps.  Both bind and
  //    unbind batch, so the excess is substantial.
  check(total_firings > total_steps,
        "binding scaled run must fire more reactions than SSA steps (bimolecular "
        "batching active)");

  // 2. First moment unbiased vs the analytic equilibrium Abound = 100.
  check(std::fabs(scaled_mean - 100.0) < 8.0,
        "binding scaled ensemble mean (Abound) must stay near the analytic "
        "equilibrium 100 (got " +
            std::to_string(scaled_mean) + ")");

  // 3. Telemetry: BOTH rules batch (association + dissociation), so both
  //    multipliers exceed 1.
  check(mult.size() == 2, "binding should expose 2 per-rule multipliers");
  if (mult.size() == 2) {
    check(mult[0] > 1.2 && mult[1] > 1.2,
          "both binding rules (bind + unbind) must carry a batch multiplier > 1 "
          "(got " +
              std::to_string(mult[0]) + ", " + std::to_string(mult[1]) + ")");
  }

  // 4. Cross-check against an exact-SSA ensemble on the same seeds.
  double exact_sum = 0.0;
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, 10.0, 100}, /*seed=*/2000 + s);
    exact_sum += final_obs(r, kAboundObs);
  }
  double const exact_mean = exact_sum / kReps;
  check(std::fabs(scaled_mean - exact_mean) < 8.0,
        "binding scaled and exact ensemble means must agree (scaled=" +
            std::to_string(scaled_mean) + " exact=" + std::to_string(exact_mean) + ")");
}

// intra_ring: A(x,y) <-> A(x!1,y!1).  Unimolecular intramolecular ring closure
// (kind 1, AddBond) + ring opening (kind 4, DeleteBond n_product==1) on a single
// molecule.  Every molecule is its own singleton complex, so the K batched
// firings are molecule-disjoint — the cleanest Phase-3 test.  Reversible
// two-state system at its analytic equilibrium (Aopen = Aclosed = 500), so the
// ensemble is stationary and any first-moment bias in either scaled batch path
// drifts the mean.
void test_intra_ring_scaled(const std::string& xml) {
  // Kept lean: the debug kBatchInvariant runs a full O(rules*mols) recompute
  // after every batch, so a short window / modest ensemble keeps this a CI smoke
  // check (the tight N=300 z-tests live in dev/psa_harness/validate_psa.py).
  constexpr int kNc = 50; // K ~ floor(500/50) = 10 for both close and open
  constexpr int kReps = 40;
  constexpr int kAclosedObs = 1; // observables: {Aopen, Aclosed}

  double sum = 0.0;
  int64_t total_firings = 0, total_steps = 0;
  std::vector<double> mult;
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_critical_population(kNc);
    auto r = sim.run({0.0, 3.0, 30}, /*seed=*/9000 + s);
    sum += final_obs(r, kAclosedObs);
    total_firings += r.ps_reaction_firings;
    total_steps += r.event_count;
    mult = r.ps_multipliers;
  }
  double const scaled_mean = sum / kReps;

  check(total_firings > total_steps,
        "intra_ring scaled run must fire more reactions than SSA steps (ring "
        "close/open batching active)");
  check(std::fabs(scaled_mean - 500.0) < 30.0,
        "intra_ring scaled ensemble mean (Aclosed) must stay near the analytic "
        "equilibrium 500 (got " +
            std::to_string(scaled_mean) + ")");
  // Telemetry: BOTH rules batch (close = kind 1 AddBond, open = kind 4
  // DeleteBond) — positive proof the ring-close/open kinds are admitted.
  check(mult.size() == 2, "intra_ring should expose 2 per-rule multipliers");
  if (mult.size() == 2)
    check(mult[0] > 1.2 && mult[1] > 1.2,
          "both intra_ring rules (close + open) must carry a batch multiplier > 1 "
          "(got " +
              std::to_string(mult[0]) + ", " + std::to_string(mult[1]) + ")");

  double exact_sum = 0.0;
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, 3.0, 30}, /*seed=*/9000 + s);
    exact_sum += final_obs(r, kAclosedObs);
  }
  double const exact_mean = exact_sum / kReps;
  check(std::fabs(scaled_mean - exact_mean) < 30.0,
        "intra_ring scaled and exact ensemble means must agree (scaled=" +
            std::to_string(scaled_mean) + " exact=" + std::to_string(exact_mean) + ")");
}

// ab_ring: finite A-B double-bond network.  Exercises the DISJOINT
// multi-molecule ring closure (kind 3, close: A(q).B(q)->A(q!1).B(q!1), a
// needs_complex_expansion rule) + ring opening (kind 4, open) + bimolecular
// association (bind) + dissociation (unbind).  `close`/`open` fire inside an
// existing A-B complex, so a scaled batch exercises the same-complex within-batch
// interaction the Phase-2 gate used to forbid.  Analytic equilibrium: every
// observable = 200, ensemble stationary, so a biased scaled batch drifts the
// mean.
void test_ab_ring_scaled(const std::string& xml) {
  // Kept lean (see intra_ring): the kind-3 ncx complex-expansion makes each
  // scaled ab_ring step the most expensive in the debug kBatchInvariant, so use
  // a short window / modest ensemble.  Analytic-vs-scaled precision is checked
  // at N=300 in validate_psa.py; here we smoke-test batching + no gross bias.
  constexpr int kNc = 40; // K ~ floor(200/40) = 5 for every rule
  constexpr int kReps = 40;
  constexpr int kDringObs = 2; // observables: {Afree, Ssingle, Dring}

  double sum = 0.0;
  int64_t total_firings = 0, total_steps = 0;
  std::vector<double> mult;
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_critical_population(kNc);
    auto r = sim.run({0.0, 3.0, 30}, /*seed=*/11000 + s);
    sum += final_obs(r, kDringObs);
    total_firings += r.ps_reaction_firings;
    total_steps += r.event_count;
    mult = r.ps_multipliers;
  }
  double const scaled_mean = sum / kReps;

  check(total_firings > total_steps,
        "ab_ring scaled run must fire more reactions than SSA steps (disjoint "
        "ring-close batching active)");
  check(std::fabs(scaled_mean - 200.0) < 25.0,
        "ab_ring scaled ensemble mean (Dring) must stay near the analytic "
        "equilibrium 200 (got " +
            std::to_string(scaled_mean) + ")");
  // Telemetry: all four rules batch (bind assoc, unbind dissoc, close = kind 3
  // disjoint ring close, open = kind 4).  The `close` multiplier > 1 is the
  // positive proof the disjoint-pattern (needs_complex_expansion) ring-closure
  // rule is admitted to batching.
  check(mult.size() == 4, "ab_ring should expose 4 per-rule multipliers");
  if (mult.size() == 4)
    for (size_t i = 0; i < mult.size(); ++i)
      check(mult[i] > 1.1, "every ab_ring rule must carry a batch multiplier > 1 (rule " +
                               std::to_string(i) + " got " + std::to_string(mult[i]) + ")");

  double exact_sum = 0.0;
  for (int s = 0; s < kReps; ++s) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r = sim.run({0.0, 3.0, 30}, /*seed=*/11000 + s);
    exact_sum += final_obs(r, kDringObs);
  }
  double const exact_mean = exact_sum / kReps;
  check(std::fabs(scaled_mean - exact_mean) < 25.0,
        "ab_ring scaled and exact ensemble means must agree (scaled=" +
            std::to_string(scaled_mean) + " exact=" + std::to_string(exact_mean) + ")");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 6) {
    std::fprintf(stderr, "Usage: psa_batch_test <birth_death.xml> <psa_pure_death.xml> "
                         "<binding.xml> <intra_ring.xml> <ab_ring.xml>\n");
    return 2;
  }
  try {
    test_birth_death_scaled(argv[1]);
    test_nc_too_small_guard(argv[2]);
    test_binding_scaled(argv[3]);
    test_intra_ring_scaled(argv[4]);
    test_ab_ring_scaled(argv[5]);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }
  if (g_failures == 0)
    std::fprintf(stderr, "psa_batch_test: all checks passed\n");
  return g_failures == 0 ? 0 : 1;
}
