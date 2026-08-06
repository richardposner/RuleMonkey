// Rate-parity and coverage test for n-ary reactant rules (issue #24).
//
// A rule with three or more reactant patterns used to be silently inert.
// The engine carries two reactant slots, and every consumer of
// `reactant_pattern_starts` treats slot B as the whole tail
// `[starts[1], molecules.size())`, so `A + A + A -> P` collapsed its second
// and third patterns into one bond-free slot-B pattern.
// `count_multi_mol_fast_generic` can satisfy the unassigned molecules of
// such a pattern only from inside the seed's own complex, so three free
// monomers scored zero, `b_total` stayed zero, and the propensity was
// identically zero.  Mass was still conserved and every other rule behaved,
// so nothing in the trajectory looked wrong — the reporter found it only by
// running the same XML through NFsim.
//
// Rules whose every reactant pattern is a single molecule now take a
// separate n-ary path (see the NaryState comment in engine.cpp), whose
// propensity is the distinct-tuple sum expanded over the partition lattice.
// This test pins down three things:
//
//   1. The reproducer fires at all, and conserves mass exactly.
//   2. The realized *rate* matches the closed-form propensity — the check
//      that a wrong symmetry factor or a missed inclusion-exclusion term
//      would fail.  A depletion model cannot do this (it saturates: nearly
//      all A is consumed for any rate fast enough to matter), so the rate
//      fixture holds its population constant by putting the reactants back
//      on the product side.  Its propensity is then constant for the whole
//      run at
//
//          a = k · sf · N(N−1)(N−2) = 1e-6 · (1/6) · 400·399·398 = 10.5868
//
//      and the terminal P count is Poisson(a·t_end).
//   3. An n-ary shape the path does *not* implement — here a multi-molecule
//      reactant pattern — is still refused at load, rather than going back
//      to silently never firing.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

double final_value(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i].back();
  throw std::runtime_error("observable not found: " + name);
}

// (1) The issue's reproducer: A(s) + A(s) + A(s) -> P() with DeleteMolecules,
//     A seeded at 400, k = 1e-6, t = 5000.  Before the fix this held A at
//     exactly 400 and P at exactly 0 forever.
void test_reproducer_fires(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);

  // The rule is implemented now, so nothing about it may be reported.
  for (const auto& f : sim.unsupported_features()) {
    check(f.element != "ListOfReactantPatterns",
          "a single-molecule-pattern trimolecular rule is implemented and must "
          "not be reported as unsupported");
  }

  auto r = sim.run({0.0, 5000.0, 2}, /*seed=*/1);
  double const a_free = final_value(r, "A_free");
  double const p_made = final_value(r, "P_made");

  check(p_made > 0, "trimolecular rule must fire (it produced exactly 0 before the fix)");

  // Exact bookkeeping: each firing consumes three A and makes one P.
  check(std::abs((a_free + (3.0 * p_made)) - 400.0) < 1e-9,
        "mass must balance exactly: A_free + 3·P_made == 400, got A_free=" +
            std::to_string(a_free) + " P_made=" + std::to_string(p_made));

  // Deterministic depletion for dN/dt = −k·N³/2 gives N(5000) = 14.1, so
  // P ≈ 128.6; NFsim reports 126-129 on this model, and 20 RM seeds span
  // 127-130.  The band is wide enough not to be flaky and still catches a
  // rule that has stopped firing.
  check(p_made > 100 && p_made < 160,
        "trimolecular fire count should land near the deterministic 128.6, got " +
            std::to_string(p_made));
}

// (2) Rate parity against the closed-form propensity, on a fixture whose
//     population does not change (reactants reappear as products).
void test_rate_matches_closed_form(const std::string& xml) {
  constexpr double kN = 400.0;
  constexpr double kK = 1.0e-6;
  constexpr double kSym = 1.0 / 6.0; // BNG2's symmetry_factor for A+A+A
  constexpr double kTEnd = 100.0;
  constexpr int kReps = 40;

  // a = k · sf · N(N−1)(N−2), constant for the whole run.
  double const a = kK * kSym * kN * (kN - 1.0) * (kN - 2.0);
  double const expected = a * kTEnd;

  double sum = 0.0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r =
        sim.run({0.0, kTEnd, 2}, /*seed=*/std::uint64_t{500} + static_cast<std::uint64_t>(rep));

    // The population must be genuinely constant, or the closed form does
    // not apply and the comparison below would be meaningless.
    check(std::abs(final_value(r, "A_free") - kN) < 1e-9,
          "rate fixture must hold A at " + std::to_string(kN) + " (it is catalytic)");
    sum += final_value(r, "P_made");
  }
  double const mean = sum / kReps;

  // Counts are Poisson(expected) per replicate, so the mean's standard
  // error is √(expected/kReps) ≈ 5.1.  A 5σ band is ~2.4% wide: loose
  // enough never to flake, tight enough that any symmetry-factor or
  // inclusion-exclusion error (which move the rate by integer factors —
  // 6× for a dropped 1/6) fails it outright.
  double const sem = std::sqrt(expected / kReps);
  double const tol = 5.0 * sem;
  check(std::abs(mean - expected) < tol,
        "trimolecular rate must match k·sf·N(N-1)(N-2): expected " + std::to_string(expected) +
            " +- " + std::to_string(tol) + ", observed " + std::to_string(mean));
}

// (3) An n-ary shape outside the implemented path must still be refused,
//     not silently dropped.  `expected_reason` is the distinguishing phrase
//     the diagnostic has to carry, so a model author can tell which of the
//     three exclusions they hit without reading the engine source.
void test_unsupported_shape_still_refused(const std::string& xml,
                                          const std::string& expected_reason) {
  const rulemonkey::RuleMonkeySimulator sim(xml);

  const rulemonkey::UnsupportedFeature* hit = nullptr;
  for (const auto& f : sim.unsupported_features()) {
    if (f.element == "ListOfReactantPatterns") {
      hit = &f;
      break;
    }
  }

  check(hit != nullptr, "an unsupported n-ary shape must be reported (" + expected_reason + ")");
  if (hit != nullptr) {
    check(hit->severity == rulemonkey::Severity::Error,
          "unsupported n-ary shapes must be Error severity — a Warn would let the "
          "model run to a silently wrong trajectory");
    check(hit->feature.find("'r'") != std::string::npos,
          "diagnostic must name the offending rule so it is findable in a large model");
    check(hit->feature.find(expected_reason) != std::string::npos,
          "diagnostic must say which shape it refused (expected '" + expected_reason +
              "'), got: " + hit->feature);
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::fprintf(stderr, "Usage: trimolecular_test <reproducer_xml> <rate_xml> "
                         "<multimol_xml> <ratelaw_xml>\n");
    return 2;
  }

  try {
    test_reproducer_fires(argv[1]);
    test_rate_matches_closed_form(argv[2]);
    // The two exclusions reachable from BNG2-written XML.  A third — more
    // than 6 reactant patterns — shares the same code path.  (MM is not
    // reachable: BNG2 refuses to write XML for MM with 3 reactants.)
    test_unsupported_shape_still_refused(argv[3], "multi-molecule");
    test_unsupported_shape_still_refused(argv[4], "rate law");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: exception: %s\n", e.what());
    ++g_failures;
  }

  if (g_failures == 0)
    std::fprintf(stdout, "trimolecular_test: OK\n");
  return g_failures == 0 ? 0 : 1;
}
