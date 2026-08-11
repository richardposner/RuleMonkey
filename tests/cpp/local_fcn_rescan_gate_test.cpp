// A local rate over a global that never moves must cost nothing (issue #40).
//
// Issue #38 gave a local rate that reads a bare observable a full O(N)
// rescan after every event: the observable moves for every instance at
// once with no molecule marked affected, so the affected-molecule delta
// path cannot see it.  That is the right answer to "can this rule read a
// moving global?" and the wrong question to gate an O(N) rescan on.  The
// question that gates the rescan is "did that global actually move?", and
// a bare observable is usually a volume proxy, a total or a clock — the
// answer is no at every event of the run, and every one of those rescans
// recomputes rates that cannot have changed.  On bngsim's
// context_symmetry fixture it was 20x; on the pair below, ~200x.
//
// Two models, identical in every respect but the middle factor of one
// local function:
//
//   ss_local_fcn_const_global        aflip(x) = kf*Vol*Aphos(x)
//   local_fcn_const_global_folded    aflip(x) = kf*VolConst*Aphos(x)
//
// `Vol` counts V(), nothing in either model creates or destroys a V, and
// `VolConst` is the parameter 50 that `Vol` equals throughout.  So the two
// are the same arithmetic on the same RNG stream — but only the first
// reads an observable, so only the first is marked as tracking a global.
// The folded twin is therefore the cost floor and the trajectory oracle at
// once, which is what makes both assertions here exact rather than
// tolerance-based:
//
//   1. same seed -> byte-identical trajectories.  Skipping a rescan whose
//      inputs did not move is exact, not an approximation, and a skip that
//      dropped a real update would show here first.
//   2. same seed -> comparable wall time.  This is the defect.  The ratio
//      is ~1.0 with the value gate in place and ~200 without it, so the
//      bound below is loose by a factor of 40 in both directions and does
//      not depend on how fast the machine running it is.
//
// The opposite direction — a bare observable that DOES move must still
// force the rescan — is pinned twice over: analytically by
// local_fcn_global_obs_test's mixed-scope arm (Boff 56.5 live against 27.1
// frozen), and here by the add_molecules arm below, which moves `Vol`
// from outside the event loop entirely.  That last one is the case the
// value cache is most exposed to: no rule sees V as a reactant, so
// add_molecules' targeted rescan does not touch RA, and the gate is the
// only thing that can notice the population changed under it.

#include "rulemonkey/simulator.hpp"

#include <chrono>
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

constexpr double kTEnd = 100.0;
constexpr int kSteps = 20;

// Wall time of run() alone — the constructor's XML parse is identical work
// for both arms and would only dilute the ratio.
double timed_run(const std::string& xml, std::uint64_t seed, rulemonkey::Result& out) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto const t0 = std::chrono::steady_clock::now();
  out = sim.run({0.0, kTEnd, kSteps}, seed);
  auto const t1 = std::chrono::steady_clock::now();
  return std::chrono::duration<double>(t1 - t0).count();
}

const std::vector<double>& series(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i];
  throw std::runtime_error("observable not found: " + name);
}

// ---------------------------------------------------------------------------
// 1. The skip is exact: the gated rule's trajectory equals its cost floor's.
// ---------------------------------------------------------------------------
void test_trajectories_identical(const std::string& bare_xml, const std::string& folded_xml) {
  for (std::uint64_t const seed : {std::uint64_t{1}, std::uint64_t{4021}}) {
    rulemonkey::Result rb;
    rulemonkey::Result rf;
    timed_run(bare_xml, seed, rb);
    timed_run(folded_xml, seed, rf);

    check(rb.observable_names == rf.observable_names, "both arms report the same observables");
    if (rb.observable_names != rf.observable_names)
      continue;

    int mismatches = 0;
    for (const auto& name : rb.observable_names) {
      const auto& a = series(rb, name);
      const auto& b = series(rf, name);
      if (a.size() != b.size()) {
        ++mismatches;
        continue;
      }
      for (size_t k = 0; k < a.size(); ++k) {
        if (a[k] != b[k]) {
          if (mismatches == 0)
            std::fprintf(stderr, "  first divergence: seed %llu, %s[%zu] %g vs %g\n",
                         static_cast<unsigned long long>(seed), name.c_str(), k, a[k], b[k]);
          ++mismatches;
        }
      }
    }
    check(mismatches == 0, "seed " + std::to_string(seed) +
                               ": the gated rule's trajectory is identical to its folded twin's");

    // Same run, stated as the model's own invariants: `Vol` is what the
    // gate reads and it never moves, and the inert half — Aphos(x) = 0 —
    // must never react, which is what keeps the TAGGED half of the
    // classification pinned while the bare half is under test.
    check(series(rb, "Vol").back() == 50.0, "Vol conserved at 50 — the gate's input never moved");
    check(series(rb, "Aphos").back() == 20000.0, "Aphos conserved at 20000");
    check(series(rb, "Aoff").back() >= 2000.0, "the 2000 inert A(p~0) never reacted");
  }
}

// ---------------------------------------------------------------------------
// 2. And it is free: reading a constant observable costs what folding it in
//    costs.  This is the regression the issue is about.
// ---------------------------------------------------------------------------
void test_cost_matches_folded_twin(const std::string& bare_xml, const std::string& folded_xml) {
  constexpr int kTimingReps = 3;
  // Measured at ~1.0.  Unconditional rescanning puts it near 200, so this
  // bound separates the two by a factor of 40 either way — wide enough to
  // survive a loaded shared runner, narrow enough that the defect cannot
  // hide under it.
  constexpr double kMaxRatio = 5.0;

  double bare = 1e300;
  double folded = 1e300;
  for (int rep = 0; rep < kTimingReps; ++rep) {
    rulemonkey::Result r;
    bare = std::fmin(bare, timed_run(bare_xml, 7, r));
    folded = std::fmin(folded, timed_run(folded_xml, 7, r));
  }

  double const ratio = folded > 0.0 ? bare / folded : 0.0;
  std::fprintf(stderr,
               "cost: bare-observable %.3fs vs folded-constant %.3fs -> %.2fx (min of %d)\n", bare,
               folded, ratio, kTimingReps);
  check(bare <= kMaxRatio * folded,
        "a local rate over a constant bare observable costs no more than the "
        "same rate with that observable folded in");
}

// ---------------------------------------------------------------------------
// 3. The gate still fires when the global moves — including when it moves
//    from outside the event loop, which is where the value cache is the
//    only thing watching.
// ---------------------------------------------------------------------------
void test_add_molecules_moves_the_global(const std::string& bare_xml) {
  // Phase 1 runs at kf*Vol = 2e-4*50 = 0.01 per reactive molecule; adding
  // 50 more V doubles it to 0.02 for phase 2.  No rule has V as a
  // reactant, so add_molecules' own targeted rescan does not reach RA:
  // only the gate noticing `Vol` moved can pick this up.
  constexpr double kSplit = 50.0;
  double const reactive_at_split = 20000.0 * std::exp(-0.01 * kSplit);
  double const live = 2000.0 + (reactive_at_split * std::exp(-0.02 * kSplit));
  double const frozen = 2000.0 + (reactive_at_split * std::exp(-0.01 * kSplit));

  rulemonkey::RuleMonkeySimulator sim(bare_xml);
  sim.initialize(/*seed=*/4022);
  sim.step_to(kSplit);
  sim.add_molecules("V", 50);
  sim.step_to(2.0 * kSplit);

  auto const names = sim.observable_names();
  auto const vals = sim.get_observable_values();
  double vol = -1.0;
  double aoff = -1.0;
  for (size_t i = 0; i < names.size(); ++i) {
    if (names[i] == "Vol")
      vol = vals[i];
    if (names[i] == "Aoff")
      aoff = vals[i];
  }

  std::fprintf(stderr, "add_molecules arm: Vol %.0f, Aoff %.0f (live %.0f, frozen %.0f)\n", vol,
               aoff, live, frozen);
  check(vol == 100.0, "add_molecules doubled the bare observable the local rate reads");
  // The two predictions are ~2900 apart on a pool whose sampling spread is
  // under 100, so a midpoint split is decisive without pinning a variance.
  check(aoff < 0.5 * (live + frozen),
        "the doubled rate reached every instance — the rescan fired on the "
        "event after the observable moved");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: local_fcn_rescan_gate_test <ss_local_fcn_const_global.xml> "
                         "<local_fcn_const_global_folded_model.xml>\n");
    return 2;
  }
  try {
    test_trajectories_identical(argv[1], argv[2]);
    test_cost_matches_folded_twin(argv[1], argv[2]);
    test_add_molecules_moves_the_global(argv[1]);
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
