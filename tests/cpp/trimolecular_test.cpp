// Rate-parity and coverage test for n-ary reactant rules (issues #24, #26).
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
// Such rules now take a separate n-ary path (see the NaryState comment in
// engine.cpp), whose propensity is the distinct-tuple sum expanded over the
// partition lattice.  A reactant pattern may be a single molecule or a
// `.`-joined complex, as it may on a bimolecular rule (#26).  This test pins
// down six things:
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
//   3. The multi-molecule reactant of #26 —
//      `A(s,d!1).D(d!1) + A(s) + A(s) -> P()` — fires and balances mass.
//   4. Overlapping patterns (a slot whose complex holds a molecule another
//      slot can also match) neither double-consume a molecule nor skew the
//      rate: the propensity counts those draws and the sampler rejects them
//      as null events, so the realized rate is the injective count exactly.
//      Checked with `-bscb` off, where that rejection is the only thing
//      standing between the two — under `-bscb` the same-complex check
//      already covers it.
//   5. A multi-molecule slot draws among the seed embeddings that actually
//      extend to a match, not among all of them: a seed molecule can offer
//      an embedding that dead-ends, and treating that as a null event would
//      halve the rate.
//   6. An n-ary shape the path does *not* implement is still refused at
//      load, rather than going back to silently never firing.

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

// (3) Issue #26's reproducer: a `.`-joined reactant complex alongside two
//     single-molecule patterns,
//
//         r: A(s,d!1).D(d!1) + A(s) + A(s) -> P()  k  DeleteMolecules
//
//     seeded with 200 A.D and 200 free A.  `A(s)` matches all 400 A — the
//     D-bound ones too, since `d` is unconstrained — so slots 1 and 2 draw
//     from a pool that includes the A inside slot 0's own complex.  RM used
//     to refuse this model outright.  Each firing takes three A and one D
//     and leaves one P, which is what makes mass balance a real check on the
//     overlap: consuming one molecule twice would break it.
//
//     NFsim on the same XML reports 110-124 P over 30 seeds (mean 117.0);
//     RM spans 109-122 (mean 117.4).
void test_multimol_reactant_fires(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);

  // The shape is implemented now, so nothing about it may be reported.
  for (const auto& f : sim.unsupported_features()) {
    check(f.element != "ListOfReactantPatterns",
          "a multi-molecule reactant pattern on an n-ary rule is implemented and "
          "must not be reported as unsupported");
  }

  auto r = sim.run({0.0, 5000.0, 2}, /*seed=*/1);
  double const a_tot = final_value(r, "A_tot");
  double const d_tot = final_value(r, "D_tot");
  double const p_made = final_value(r, "P_made");

  check(p_made > 0, "an n-ary rule with a multi-molecule reactant must fire "
                    "(the model was refused at load before issue #26)");

  check(std::abs((a_tot + (3.0 * p_made)) - 400.0) < 1e-9,
        "mass must balance exactly: A_tot + 3·P_made == 400, got A_tot=" + std::to_string(a_tot) +
            " P_made=" + std::to_string(p_made));
  check(std::abs((d_tot + p_made) - 200.0) < 1e-9,
        "mass must balance exactly: D_tot + P_made == 200, got D_tot=" + std::to_string(d_tot) +
            " P_made=" + std::to_string(p_made));

  check(p_made > 85 && p_made < 150,
        "fire count should land in NFsim's 110-124 band for this model, got " +
            std::to_string(p_made));
}

// (4) Rate parity for overlapping patterns, on a fixture whose population
//     does not change.  The rule is
//
//         r: D(d!1).A(s,d!1) + A(s) + A(s) -> (same) + P()  k
//
//     seeded with 5 A.D complexes and 5 free A, so slot 0 is seeded on a D
//     and drags in the A bonded to it — an A that slots 1 and 2 can draw
//     as well.  With c ≡ 1 everywhere:
//
//       seed-distinct sum (the propensity)  D     = 5 · 10 · 9 = 450
//       fully injective sum (what may fire) D_inj = 5 ·  9 · 8 = 360
//
//     so the realized rate is k·sf·D_inj = 0.05 · 0.5 · 360 = 9 per unit
//     time, and dropping the injectivity rejection would run 25% hot at
//     11.25.  Run with `-bscb` off: on, the same-complex check rejects the
//     same draws first and the injectivity check never gets to prove
//     itself.  (Both settings are exercised; they coincide here because no
//     complex in this model holds more than one A.)
void test_overlap_rate_matches_closed_form(const std::string& xml, bool bscb) {
  constexpr double kTEnd = 100.0;
  constexpr int kReps = 40;
  constexpr double kExpectedRate = 0.05 * 0.5 * 360.0; // k · sf · D_inj
  double const expected = kExpectedRate * kTEnd;

  double sum = 0.0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_block_same_complex_binding(bscb);
    auto r =
        sim.run({0.0, kTEnd, 2}, /*seed=*/std::uint64_t{900} + static_cast<std::uint64_t>(rep));

    check(std::abs(final_value(r, "A_tot") - 10.0) < 1e-9,
          "overlap fixture must hold A at 10 (it is catalytic)");
    sum += final_value(r, "P_made");
  }
  double const mean = sum / kReps;

  // Poisson(expected) per replicate → the mean's standard error is
  // √(expected/kReps) ≈ 4.7, so a 5σ band is ~2.6% wide.  The failure it
  // exists to catch — no injectivity rejection — is 25% off.
  double const tol = 5.0 * std::sqrt(expected / kReps);
  check(std::abs(mean - expected) < tol,
        std::string("overlapping n-ary rate must match k·sf·D_inj with bscb ") +
            (bscb ? "on" : "off") + ": expected " + std::to_string(expected) + " +- " +
            std::to_string(tol) + ", observed " + std::to_string(mean));
}

// (5) A multi-molecule slot whose seed offers embeddings that go nowhere.
//     `A(d,d)` bonded to one D and one E has two embeddings of the pattern's
//     seed molecule `A(d!1)` — one per bonded `d` — but only the one that
//     reaches the D extends to a match of `A(d!1).D(d!1)`.  `c_0(m)` counts
//     the extending one, so the sampler has to draw among those; drawing
//     over all seed embeddings and calling the dead end a null event would
//     halve the rate here.  The rule is
//
//         r: A(d!1).D(d!1) + X() + X() -> (same) + P()  k
//
//     with 100 A.D.E complexes and 100 X, catalytic, so
//
//         a = k · sf · T_0 · N_X(N_X−1) = 1e-6 · 0.5 · 100 · 100 · 99 = 0.495
//
//     NFsim makes 51.1 P over 10 seeds against this 49.5.
void test_deadend_seed_embeddings(const std::string& xml) {
  constexpr double kTEnd = 100.0;
  constexpr int kReps = 40;
  double const expected = 1.0e-6 * 0.5 * 100.0 * 100.0 * 99.0 * kTEnd;

  double sum = 0.0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r =
        sim.run({0.0, kTEnd, 2}, /*seed=*/std::uint64_t{300} + static_cast<std::uint64_t>(rep));
    check(std::abs(final_value(r, "AD") - 100.0) < 1e-9,
          "dead-end fixture must hold 100 A.D pairs (it is catalytic)");
    sum += final_value(r, "P_made");
  }
  double const mean = sum / kReps;

  // 5σ is ~11% wide here; the failure it exists to catch is 2×.
  double const tol = 5.0 * std::sqrt(expected / kReps);
  check(std::abs(mean - expected) < tol,
        "rate must count only the seed embeddings that extend: expected " +
            std::to_string(expected) + " +- " + std::to_string(tol) + ", observed " +
            std::to_string(mean));
}

// (6) The same overlap under DeleteMolecules: every firing deletes the D,
//     the A bonded to it, and the two A drawn for slots 1 and 2.  If a draw
//     that put the bonded A in slot 1 or 2 were allowed through, the rule
//     would delete one molecule twice and the books would not balance.
//     Run with `-bscb` off — with it on, the same-complex check would mask
//     the overlap before the injectivity check saw it.
void test_overlap_delete_conserves_mass(const std::string& xml) {
  for (int seed = 1; seed <= 8; ++seed) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_block_same_complex_binding(false);
    auto r = sim.run({0.0, 5000.0, 2}, static_cast<std::uint64_t>(seed));

    double const a_tot = final_value(r, "A_tot");
    double const d_tot = final_value(r, "D_tot");
    double const p_made = final_value(r, "P_made");
    std::string const at = " (seed " + std::to_string(seed) + ", A_tot=" + std::to_string(a_tot) +
                           " D_tot=" + std::to_string(d_tot) + " P_made=" + std::to_string(p_made) +
                           ")";

    check(p_made > 0, "overlapping DeleteMolecules rule must fire" + at);
    check(std::abs((a_tot + (3.0 * p_made)) - 400.0) < 1e-9, "A_tot + 3·P_made == 400" + at);
    check(std::abs((d_tot + p_made) - 200.0) < 1e-9, "D_tot + P_made == 200" + at);
  }
}

// (7) An n-ary shape outside the implemented path must still be refused,
//     not silently dropped.  `expected_reason` is the distinguishing phrase
//     the diagnostic has to carry, so a model author can tell which of the
//     exclusions they hit without reading the engine source.
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
  if (argc < 10) {
    std::fprintf(stderr, "Usage: trimolecular_test <reproducer_xml> <rate_xml> "
                         "<multimol_xml> <overlap_xml> <deadend_xml> "
                         "<overlap_delete_xml> <disjoint_xml> <toomany_xml> "
                         "<ratelaw_xml>\n");
    return 2;
  }

  try {
    test_reproducer_fires(argv[1]);
    test_rate_matches_closed_form(argv[2]);
    test_multimol_reactant_fires(argv[3]);
    test_overlap_rate_matches_closed_form(argv[4], /*bscb=*/false);
    test_overlap_rate_matches_closed_form(argv[4], /*bscb=*/true);
    test_deadend_seed_embeddings(argv[5]);
    test_overlap_delete_conserves_mass(argv[6]);
    // The three exclusions reachable from BNG2-written XML.  (MM is not
    // reachable: BNG2 refuses to write XML for MM with 3 reactants.)
    test_unsupported_shape_still_refused(argv[7], "disconnected complex");
    test_unsupported_shape_still_refused(argv[8], "n-ary limit");
    test_unsupported_shape_still_refused(argv[9], "rate law");
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: exception: %s\n", e.what());
    ++g_failures;
  }

  if (g_failures == 0)
    std::fprintf(stdout, "trimolecular_test: OK\n");
  return g_failures == 0 ? 0 : 1;
}
