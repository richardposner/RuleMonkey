// Post-load seed-species amount API (issue #23): initial_species,
// get_initial_amount, set_initial_amount, clear_initial_amount_overrides.
//
// The motivating case is a driver walking a dose-response curve on one
// loaded model.  Where the dose is parameter-driven, set_param is the
// right entry point and already carries the derivation (see
// set_param_test.cpp).  These methods cover what set_param cannot reach —
// a seed amount written as a bare literal, with no parameter behind it —
// and give a caller a way to read back what the next run will actually
// seed without having to run it first.
//
// Assertions are deterministic throughout: seeded counts are integers and
// every check is made at t=0 or against parameter values, so no
// trajectory statistics are involved.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdio>
#include <functional>
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

bool throws(const std::function<void()>& fn) {
  try {
    fn();
  } catch (const std::runtime_error&) {
    return true;
  } catch (...) { // NOLINT(bugprone-empty-catch)
  }
  return false;
}

double initial_value(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return r.observable_data[i].front();
  throw std::runtime_error("observable not found: " + name);
}

const rulemonkey::InitialSpeciesRow&
row_named(const std::vector<rulemonkey::InitialSpeciesRow>& rows, const std::string& name) {
  for (const auto& row : rows)
    if (row.name == name)
      return row;
  throw std::runtime_error("seed species not found: " + name);
}

// A_plus_A.xml: one seed species `A(a)` at concentration="A_tot", A_tot=1000.
void test_introspection(const std::string& xml) {
  const rulemonkey::RuleMonkeySimulator sim(xml);
  auto rows = sim.initial_species();
  check(rows.size() == 1, "A_plus_A should report exactly one seed species (got " +
                              std::to_string(rows.size()) + ")");
  const auto& a = rows.at(0);
  check(a.name == "A(a)", "seed species name should be the BNGL pattern (got '" + a.name + "')");
  check(!a.id.empty(), "seed species should carry its XML id");
  check(a.concentration_expr == "A_tot",
        "seed species should report its concentration attribute verbatim (got '" +
            a.concentration_expr + "')");
  check(a.amount == 1000.0,
        "seed amount should resolve to A_tot=1000 (got " + std::to_string(a.amount) + ")");
  check(!a.overridden, "a freshly loaded seed species should not be marked overridden");

  // Both keys address the same species.
  check(sim.get_initial_amount("A(a)") == 1000.0, "get_initial_amount by BNGL name");
  check(sim.get_initial_amount(a.id) == 1000.0, "get_initial_amount by XML id");
}

// The amount reported must track parameter overrides WITHOUT an intervening
// run — the same between-runs coherence get_parameter() offers.
void test_amount_tracks_param_override(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_param("A_tot", 250.0);
  check(sim.get_initial_amount("A(a)") == 250.0,
        "get_initial_amount should reflect a set_param override immediately (got " +
            std::to_string(sim.get_initial_amount("A(a)")) + ")");
  check(!sim.initial_species().at(0).overridden,
        "a parameter-driven amount is not a direct override");

  sim.clear_param_overrides();
  check(sim.get_initial_amount("A(a)") == 1000.0,
        "get_initial_amount should snap back after clear_param_overrides");
}

void test_set_initial_amount_reaches_engine(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_block_same_complex_binding(true);
  sim.set_initial_amount("A(a)", 37.0);
  check(sim.get_initial_amount("A(a)") == 37.0, "set_initial_amount should be readable back");
  check(sim.initial_species().at(0).overridden, "a pinned species should report overridden=true");

  auto r = sim.run({0.0, 1.0, 2}, /*seed=*/1);
  check(initial_value(r, "A_1") == 37.0,
        "set_initial_amount should seed exactly 37 monomers (got " +
            std::to_string(initial_value(r, "A_1")) + ")");

  // Sessions honor it too, and repeated runs stay pinned.
  sim.initialize(/*seed=*/1);
  check(sim.get_molecule_count("A") == 37,
        "initialize() should seed the pinned amount as well (got " +
            std::to_string(sim.get_molecule_count("A")) + ")");
  sim.destroy_session();

  // Real-valued amounts truncate toward zero, matching the engine's
  // NFsim-parity handling of a fractional <Species concentration=>.
  sim.set_initial_amount("A(a)", 12.9);
  auto r_frac = sim.run({0.0, 1.0, 2}, /*seed=*/1);
  check(initial_value(r_frac, "A_1") == 12.0,
        "a fractional pinned amount should truncate toward zero (got " +
            std::to_string(initial_value(r_frac, "A_1")) + ")");
}

// A pin is stated in molecules and must therefore outrank whatever the
// concentration expression derives — including a parameter override applied
// afterwards.  Clearing the pin hands control back to the expression.
void test_pin_outranks_param_and_clears(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_block_same_complex_binding(true);

  sim.set_initial_amount("A(a)", 60.0);
  sim.set_param("A_tot", 900.0);
  check(sim.get_initial_amount("A(a)") == 60.0,
        "a direct pin should outrank a later set_param on the driving parameter");
  check(sim.get_parameter("A_tot") == 900.0, "the parameter itself still reports its override");
  auto r_pinned = sim.run({0.0, 1.0, 2}, /*seed=*/1);
  check(initial_value(r_pinned, "A_1") == 60.0, "the pinned amount should reach the engine");

  // Clearing param overrides must not disturb the pin.
  sim.clear_param_overrides();
  check(sim.get_initial_amount("A(a)") == 60.0,
        "clear_param_overrides should leave the direct pin in force");

  // Releasing the pin restores the expression-derived amount, including in
  // the engine — the parsed model must be un-baked, not just the report.
  sim.clear_initial_amount_overrides();
  check(sim.get_initial_amount("A(a)") == 1000.0,
        "clear_initial_amount_overrides should restore the derived amount");
  check(!sim.initial_species().at(0).overridden, "the overridden flag should clear too");
  auto r_cleared = sim.run({0.0, 1.0, 2}, /*seed=*/1);
  check(initial_value(r_cleared, "A_1") == 1000.0,
        "clearing the pin after a run should restore the seeded count to 1000, not leave the "
        "baked pin at 60");
}

// A dose scan on one loaded model: the whole point of the API.
void test_repeated_scan_points(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_block_same_complex_binding(true);
  for (const double dose : {10.0, 250.0, 0.0, 777.0}) {
    sim.set_initial_amount("A(a)", dose);
    auto r = sim.run({0.0, 1.0, 2}, /*seed=*/1);
    check(initial_value(r, "A_1") == dose, "scan point A_tot=" + std::to_string(dose) +
                                               " should seed its own dose (got " +
                                               std::to_string(initial_value(r, "A_1")) + ")");
  }
}

void test_rejections(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);

  check(throws([&] { sim.set_initial_amount("Z(nope)", 5.0); }),
        "set_initial_amount on an unknown species should throw");
  check(throws([&] { (void)sim.get_initial_amount("Z(nope)"); }),
        "get_initial_amount on an unknown species should throw");
  check(throws([&] { sim.set_initial_amount("A(a)", -1.0); }),
        "set_initial_amount with a negative amount should throw");
  check(throws([&] { sim.set_initial_amount("A(a)", std::nan("")); }),
        "set_initial_amount with a non-finite amount should throw");

  // A rejected call must leave no trace.
  check(sim.get_initial_amount("A(a)") == 1000.0,
        "a rejected set_initial_amount should not perturb the seed amount");
  check(!sim.initial_species().at(0).overridden,
        "a rejected set_initial_amount should not mark the species overridden");

  // Same no-mid-session-mutation contract as the other configuration setters.
  sim.initialize(/*seed=*/1);
  check(throws([&] { sim.set_initial_amount("A(a)", 5.0); }),
        "set_initial_amount should throw during an active session");
  check(throws([&] { sim.clear_initial_amount_overrides(); }),
        "clear_initial_amount_overrides should throw during an active session");
  // Reads stay legal mid-session — they describe the initial state, which
  // does not change while a session runs.
  check(sim.get_initial_amount("A(a)") == 1000.0, "get_initial_amount should work mid-session");
  check(sim.initial_species().size() == 1, "initial_species should work mid-session");
  sim.destroy_session();
  check(!throws([&] { sim.set_initial_amount("A(a)", 5.0); }),
        "set_initial_amount should succeed again after destroy_session");
}

// ft_clamped_species exercises the two shapes set_param cannot serve.
//
//   <Species id="S3" concentration="0"      name="P()">        — a bare
//        literal with no parameter behind it, so there is nothing to
//        set_param.  This is the case the API exists for.
//   <Species id="S1" concentration="E_tot"  name="$E(s)" Fixed="1"> — a
//        clamped species: the engine holds its population at the seeded
//        value all run, so the pin must move the clamp target too rather
//        than let the run replenish toward the XML-time count.
void test_literal_and_clamped_species(const std::string& xml) {
  const rulemonkey::RuleMonkeySimulator probe(xml);
  check(probe.get_initial_amount("P()") == 0.0, "fixture sanity: P() seeds 0");
  check(probe.initial_species().size() == 3, "fixture sanity: three seed species");
  check(row_named(probe.initial_species(), "P()").concentration_expr == "0",
        "P() should report its literal concentration attribute, with no parameter behind it");

  // The literal case: no parameter names this amount, so set_param has no
  // handle on it and only set_initial_amount can move it.
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_initial_amount("P()", 25.0);
  auto r = sim.run({0.0, 0.0, 1}, /*seed=*/1);
  check(initial_value(r, "P_total") == 25.0,
        "set_initial_amount should seed a literal-concentration species (got " +
            std::to_string(initial_value(r, "P_total")) + ")");

  // The clamped case.  E is consumed into ES and released by catalysis, so
  // free E fluctuates; what the clamp pins is the total E population, and
  // at t=0 that is entirely free.  A clamp target still reading the
  // XML-time 10 would show up here as 10.
  check(probe.get_initial_amount("$E(s)") == 10.0, "fixture sanity: $E(s) seeds E_tot=10");
  rulemonkey::RuleMonkeySimulator clamp_sim(xml);
  clamp_sim.set_initial_amount("$E(s)", 17.0);
  clamp_sim.initialize(/*seed=*/1);
  check(clamp_sim.get_observable_values()[0] == 17.0,
        "a pinned clamped species should seed 17 free E at t=0 (got " +
            std::to_string(clamp_sim.get_observable_values()[0]) + ")");

  // Run it forward: the clamp must hold the TOTAL E population at the pin,
  // not at the XML-time 10.
  (void)clamp_sim.simulate(0.0, 5.0, 6);
  check(clamp_sim.get_molecule_count("E") == 17,
        "the Fixed clamp target should follow the pin, holding E at 17 (got " +
            std::to_string(clamp_sim.get_molecule_count("E")) + ")");
  clamp_sim.destroy_session();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <A_plus_A.xml> <ft_clamped_species.xml>\n", argv[0]);
    return 2;
  }
  const std::string a_plus_a = argv[1];
  const std::string clamped = argv[2];

  try {
    test_introspection(a_plus_a);
    test_amount_tracks_param_override(a_plus_a);
    test_set_initial_amount_reaches_engine(a_plus_a);
    test_pin_outranks_param_and_clears(a_plus_a);
    test_repeated_scan_points(a_plus_a);
    test_rejections(a_plus_a);
    test_literal_and_clamped_species(clamped);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: seed amounts are introspectable, pinnable, and releasable; pins "
                       "outrank parameter derivations, reach the engine, and un-bake on clear\n");
  return 0;
}
