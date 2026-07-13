// Partial-scaling (opt-in) Phase-1 behavioral test.
//
// Exercises the scaled batch-fire path that `set_critical_population(Nc)`
// enables (Lin, Feng & Hlavacek 2019).  The exact/large-Nc byte-identity
// contract is pinned in api_coverage_test; here we drive the *scaled* path and
// assert the properties that make it correct and useful:
//
//   1. Unimolecular batching actually happens (ps_reaction_firings exceeds the
//      SSA step count on a large-population unimolecular model).
//   2. First moment stays unbiased: the ensemble mean of the scaled run agrees
//      with the analytic Poisson mean and with the exact-SSA ensemble.
//   3. Per-rule multiplier telemetry is sane: the batched unimolecular rule has
//      a time-averaged multiplier > 1, the zeroth-order rule stays at 1.
//   4. The Nc-too-small guard throws on a model whose only channel is the
//      scaled rule when Nc is aggressive enough to collapse the run.
//
// Under Debug / ASan builds the partial-scaling batch self-check
// (kBatchInvariant, plan §5) runs inside every scaled batch here, so this test
// also exercises the epoch-cache reconciliation invariant in CI.
//
// argv[1] = birth_death.xml   (0 -> A ; A -> 0, k_death=1, A0=100 => Poisson(100))
// argv[2] = psa_pure_death.xml (A -> 0 only; no source)

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

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: psa_batch_test <birth_death.xml> <psa_pure_death.xml>\n");
    return 2;
  }
  try {
    test_birth_death_scaled(argv[1]);
    test_nc_too_small_guard(argv[2]);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }
  if (g_failures == 0)
    std::fprintf(stderr, "psa_batch_test: all checks passed\n");
  return g_failures == 0 ? 0 : 1;
}
