// Pattern-keyed species-method integration test (issue #9 §1 step 3).
//
// Drives RuleMonkeySimulator's four §1 session methods —
// get_species_count / add_species / remove_species / set_species_count
// — through the public API against a live session:
//   - get_species_count agrees with enumerate_species() row counts;
//   - add_species / remove_species round-trip, single- and
//     multi-molecule (the A(a!1).A(a!1) homodimer exercises the
//     multi-molecule parse + connectivity + symmetric-complex path);
//   - set_species_count drives the live count to an exact target both
//     up and down;
//   - the SSA loop still runs after a structural mutation (the resync
//     left rule propensities valid);
//   - the documented error surface (no session, malformed /
//     under-specified pattern, non-positive counts, over-removal)
//     throws std::runtime_error.
//
// Model: A_plus_A.xml — molecule type A with one stateless component a,
// a reversible A+A binding rule.  The test is deterministic: dimers are
// created with add_species rather than relied upon to form by chance.

#include "rulemonkey/simulator.hpp"

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

// Assert that calling `fn` throws std::runtime_error.
template <typename Fn> void expect_throw(Fn&& fn, const std::string& what) {
  try {
    fn();
  } catch (const std::runtime_error&) {
    return;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "FAIL: %s threw a non-runtime_error: %s\n", what.c_str(), e.what());
    ++g_failures;
    return;
  }
  std::fprintf(stderr, "FAIL: %s did not throw\n", what.c_str());
  ++g_failures;
}

const char* kMonomer = "A(a)";
const char* kDimer = "A(a!1).A(a!1)";

// get_species_count must agree with the enumerate_species() row counts,
// and their sum must equal total_complex_count().
void test_agrees_with_enumeration(rulemonkey::RuleMonkeySimulator& sim) {
  auto rows = sim.enumerate_species();
  long sum = 0;
  for (const auto& row : rows) {
    const int viaPattern = sim.get_species_count(row.species);
    check(viaPattern == static_cast<int>(row.count),
          "get_species_count('" + row.species + "') should equal its enumerate_species row count");
    sum += row.count;
  }
  check(sum == sim.total_complex_count(),
        "enumerate_species row counts should sum to total_complex_count()");
}

// add_species / remove_species / set_species_count round-trip on the
// single-molecule species A(a).
void test_single_molecule_round_trip(rulemonkey::RuleMonkeySimulator& sim) {
  const int m0 = sim.get_species_count(kMonomer);
  const long total0 = sim.total_complex_count();

  sim.add_species(kMonomer, 7);
  check(sim.get_species_count(kMonomer) == m0 + 7, "add_species(A(a),7) raises the count by 7");
  check(sim.total_complex_count() == total0 + 7,
        "add_species of 7 monomers adds 7 complexes to the pool");

  sim.remove_species(kMonomer, 4);
  check(sim.get_species_count(kMonomer) == m0 + 3, "remove_species(A(a),4) lowers the count by 4");

  sim.set_species_count(kMonomer, m0);
  check(sim.get_species_count(kMonomer) == m0, "set_species_count(A(a),m0) drives the count down");

  sim.set_species_count(kMonomer, m0 + 20);
  check(sim.get_species_count(kMonomer) == m0 + 20, "set_species_count(A(a),m0+20) drives it up");

  sim.set_species_count(kMonomer, m0); // restore
  check(sim.get_species_count(kMonomer) == m0, "set_species_count restores the original count");
}

// The multi-molecule homodimer A(a!1).A(a!1): add/remove must round-trip
// and the count must agree with enumerate_species afterwards.
void test_multi_molecule_species(rulemonkey::RuleMonkeySimulator& sim) {
  check(sim.get_species_count(kDimer) == 0, "no A(a!1).A(a!1) dimers exist before add_species");

  const long total0 = sim.total_complex_count();
  sim.add_species(kDimer, 5);
  check(sim.get_species_count(kDimer) == 5, "add_species(A(a!1).A(a!1),5) creates 5 dimers");
  check(sim.total_complex_count() == total0 + 5,
        "5 dimers add 5 complexes (each dimer is one complex of two molecules)");
  test_agrees_with_enumeration(sim);

  sim.remove_species(kDimer, 2);
  check(sim.get_species_count(kDimer) == 3, "remove_species(A(a!1).A(a!1),2) leaves 3 dimers");
}

// After the structural mutations above, the SSA loop must still run:
// the propensity / observable resync left the engine in a valid state.
void test_ssa_runs_after_mutation(rulemonkey::RuleMonkeySimulator& sim) {
  bool threw = false;
  try {
    sim.simulate(sim.current_time(), sim.current_time() + 1.0, 3);
  } catch (const std::exception& e) {
    threw = true;
    std::fprintf(stderr, "FAIL: simulate() after a species mutation threw: %s\n", e.what());
    ++g_failures;
  }
  check(!threw, "the SSA loop runs after add_species / remove_species");
  // The session is still queryable after the segment.
  check(sim.get_species_count(kMonomer) >= 0, "get_species_count is callable post-simulate");
}

// Issue #13: get/match must be component-order-insensitive, consistent
// with the add path.  The species_order model declares X with a
// stateful `p~0~1` and a stateless `y`, seeded with 5000 X(p~0,y).  An
// exact pattern written with the components in either order denotes the
// same species, so get_species_count must return 5000 for both — and
// set_species_count must reach an exact target regardless of the order
// the caller wrote the pattern in (the silent wrong-baseline bug).
void test_component_order_insensitive(const std::string& order_xml) {
  rulemonkey::RuleMonkeySimulator sim(order_xml);
  sim.initialize(/*seed=*/11);

  const int n = sim.get_species_count("X(p~0,y)");
  check(n == 5000, "X(p~0,y) seeds 5000 copies");
  check(sim.get_species_count("X(y,p~0)") == n,
        "get_species_count is component-order-insensitive: X(y,p~0) == X(p~0,y)");

  // Mutate the *other* species (X(p~1,y)) through both orders and
  // confirm the match path tracks it either way.
  sim.add_species("X(y,p~1)", 25); // non-canonical order on the add path
  check(sim.get_species_count("X(p~1,y)") == 25,
        "add_species with swapped components creates the species the canonical order matches");
  check(sim.get_species_count("X(y,p~1)") == 25, "...and the swapped-order query agrees");

  // set_species_count with a NON-canonical order must reach exactly N,
  // not diff against a wrong-order baseline of 0 and overshoot.
  sim.set_species_count("X(p~1,y)", 100); // canonical order
  check(sim.get_species_count("X(y,p~1)") == 100, "set to 100 via canonical order lands on 100");
  sim.set_species_count("X(y,p~1)", 40); // non-canonical order, drives down
  check(sim.get_species_count("X(p~1,y)") == 40,
        "set_species_count with swapped components reaches the exact target");

  // remove_species via the swapped order targets the right species.
  sim.remove_species("X(y,p~1)", 40);
  check(sim.get_species_count("X(p~1,y)") == 0,
        "remove_species with swapped components removes the matched species");

  sim.destroy_session();
}

void test_error_surface(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);

  // No session yet.
  expect_throw([&] { sim.get_species_count(kMonomer); }, "get_species_count with no session");
  expect_throw([&] { sim.add_species(kMonomer, 1); }, "add_species with no session");

  sim.initialize(/*seed=*/5);

  // Malformed / unresolved / under-specified patterns.
  expect_throw([&] { sim.get_species_count("Zzz(q)"); }, "get_species_count of unknown type");
  expect_throw([&] { sim.get_species_count("A("); }, "get_species_count of malformed string");
  expect_throw([&] { sim.get_species_count("A()"); }, "get_species_count of under-specified A()");

  // Non-positive / negative counts.
  expect_throw([&] { sim.add_species(kMonomer, 0); }, "add_species with count 0");
  expect_throw([&] { sim.remove_species(kMonomer, 0); }, "remove_species with count 0");
  expect_throw([&] { sim.set_species_count(kMonomer, -1); }, "set_species_count with count -1");

  // Removing more copies than exist.
  expect_throw([&] { sim.remove_species(kMonomer, 1000000000); }, "remove_species over-removal");

  sim.destroy_session();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: species_methods_test <A_plus_A.xml> <species_order_model.xml>\n");
    return 2;
  }
  const std::string xml = argv[1];
  const std::string order_xml = argv[2];

  try {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.initialize(/*seed=*/7);

    test_agrees_with_enumeration(sim);
    test_single_molecule_round_trip(sim);
    test_multi_molecule_species(sim);
    test_ssa_runs_after_mutation(sim);
    sim.destroy_session();

    test_error_surface(xml);
    test_component_order_insensitive(order_xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: pattern-keyed species-method assertions all passed\n");
  return 0;
}
