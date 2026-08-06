// Rate parity for a bimolecular rule whose reactant is a multi-molecule
// complex, when the seed molecule offers embeddings that go nowhere.
//
// The sampler draws the seed molecule weighted by count_a, which
// count_multi_mol_fast defines as the number of seed embeddings that extend
// to a whole match — call it c.  It then had to pick one of the seed
// embeddings themselves, of which there are S >= c, and a pick that did not
// extend became a null event.  The realized rate therefore carried an extra
// factor c/S: the rule fired slow by exactly the fraction of its seed
// embeddings that dead-end, silently, with nothing in the trajectory looking
// wrong.
//
// The fixture makes that factor 1/2.  `A(d,d)` is bonded to a D on one `d`
// and an E on the other, so the pattern's seed molecule `A(d!1)` has two
// embeddings — one per bonded `d` — and only the one that lands on the D
// extends to a match of `A(d!1).D(d!1)`.  The rule is catalytic, so every
// population is constant and the propensity holds for the whole run at
//
//     a = k · sf · N_AD · N_X = 1e-4 · 1 · 100 · 100 = 1 per unit time
//
// giving 100 firings over t = 100, Poisson.  Before the fix RM produced 51.7
// on this model over 10 seeds where NFsim produced 98.8.
//
// The counting side is checked too: the AD observable is what count_a is
// derived from, and it has to read 100 (one per complex), not 200 (one per
// seed embedding).  A "fix" that inflated the count to S instead of drawing
// among the extending embeddings would land the rate correctly while making
// every Molecules observable on a multi-molecule pattern wrong, so pinning
// both is what separates the two.
//
// The second case here is the other way a multi-molecule pattern can go
// wrong: `mol_a != mol_b` separates the two *seeds*, but a molecule the
// pattern pulled in through its bonds can be the other slot's match as well.
// `A(s,d!1).A(s,d!1) + A(s) -> P()` with DeleteMolecules can draw the
// dimer's second A for slot B, delete it twice, and leave the books short.
// `-bscb` hides it — the two seeds are in one complex, so the same-complex
// check rejects the draw first — so the fixture runs with it off.

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

void test_deadend_seed_embeddings(const std::string& xml) {
  constexpr double kTEnd = 100.0;
  constexpr int kReps = 40;
  constexpr double kExpectedRate = 1.0e-4 * 1.0 * 100.0 * 100.0; // k · sf · N_AD · N_X
  double const expected = kExpectedRate * kTEnd;

  double sum = 0.0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    auto r =
        sim.run({0.0, kTEnd, 2}, /*seed=*/std::uint64_t{700} + static_cast<std::uint64_t>(rep));

    // Constant populations, or the closed form above does not apply.
    check(std::abs(final_value(r, "X_total") - 100.0) < 1e-9,
          "fixture must hold X at 100 (it is catalytic)");
    // One embedding per A.D complex — the pattern's second embedding of the
    // seed molecule reaches the E, not the D, and must not be counted.
    check(std::abs(final_value(r, "AD") - 100.0) < 1e-9,
          "A(d!1).D(d!1) must count 100, one per complex, not one per seed "
          "embedding; got " +
              std::to_string(final_value(r, "AD")));
    sum += final_value(r, "P_made");
  }
  double const mean = sum / kReps;

  // Poisson(expected) per replicate → the mean's standard error is
  // √(expected/kReps) = 1.58, so 5σ is ~8% wide.  The failure it exists to
  // catch is 2×.
  double const tol = 5.0 * std::sqrt(expected / kReps);
  check(std::abs(mean - expected) < tol,
        "rate must count only the seed embeddings that extend to a whole "
        "match: expected " +
            std::to_string(expected) + " +- " + std::to_string(tol) + ", observed " +
            std::to_string(mean));
}

// A molecule reached through a pattern's bonds must not also be the other
// slot's match.  Five dimers and nothing else, so `A(s)` for slot B draws
// from the same ten molecules the dimer pattern spans and the collision is
// frequent (~1 in 10 draws) rather than a tail event; every firing takes
// three A, so the books have to close at A_tot + 3·P_made == 10 exactly.
// Run with `-bscb` off — with it on, the same-complex check rejects the
// collision before injectivity is ever consulted.
void test_overlap_conserves_mass(const std::string& xml) {
  int fired = 0;
  for (int seed = 1; seed <= 40; ++seed) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_block_same_complex_binding(false);
    auto r = sim.run({0.0, 2000.0, 2}, static_cast<std::uint64_t>(seed));

    double const a_tot = final_value(r, "A_tot");
    double const p_made = final_value(r, "P_made");
    fired += static_cast<int>(p_made);
    check(std::abs((a_tot + (3.0 * p_made)) - 10.0) < 1e-9,
          "a molecule inside the reactant complex must not be consumed twice: "
          "A_tot + 3·P_made == 10 (seed " +
              std::to_string(seed) + "), got A_tot=" + std::to_string(a_tot) +
              " P_made=" + std::to_string(p_made));
  }
  // Guard against the fixture going inert and passing vacuously.
  check(fired > 40, "overlap fixture must actually fire (got " + std::to_string(fired) +
                        " firings over 40 seeds)");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: multimol_rate_test <deadend_xml> <overlap_xml>\n");
    return 2;
  }

  try {
    test_deadend_seed_embeddings(argv[1]);
    test_overlap_conserves_mass(argv[2]);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: exception: %s\n", e.what());
    ++g_failures;
  }

  if (g_failures == 0)
    std::fprintf(stdout, "multimol_rate_test: OK\n");
  return g_failures == 0 ? 0 : 1;
}
