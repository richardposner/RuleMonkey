// Regression test: set_param overrides reach the engine for the rate
// laws and initial-species concentrations that read parameter-derived
// fields baked at XML parse time.
//
// Before the apply_overrides re-resolve fix, set_param wrote into the
// override map but the engine read parse-time-baked rate_value /
// mm_kcat / mm_Km / SpeciesInit::concentration directly, so overrides
// were silent no-ops for the most common parameter shapes.
//
// Three models exercised:
//   A_plus_A.xml             — Ele rate (kp), initial conc (A_tot)
//   ft_mm_ratelaw.xml        — MM rate (kcat, Km)
//   derived_param_model.xml  — derived-parameter cascade
//                              (A_tot=A_base*A_factor, kp=kp_base*kp_mult)
//                              and get_parameter() coherence with set_param.
//
// All assertions are deterministic: zero-rate cases produce zero
// reaction events, and integer initial counts are exact.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
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
  int i = idx_of(r, name);
  if (i < 0)
    throw std::runtime_error("observable not found: " + name);
  return r.observable_data[i].back();
}

double initial_value(const rulemonkey::Result& r, const std::string& name) {
  int i = idx_of(r, name);
  if (i < 0)
    throw std::runtime_error("observable not found: " + name);
  return r.observable_data[i].front();
}

rulemonkey::Result run_with(const std::string& xml,
                            const std::vector<std::pair<std::string, double>>& overrides,
                            double t_end, int n_points, std::uint64_t seed) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_block_same_complex_binding(true);
  for (auto& [k, v] : overrides)
    sim.set_param(k, v);
  return sim.run({0.0, t_end, n_points}, seed);
}

// A_plus_A defaults: kp=0.001, km=0.5, A_tot=1000.
void test_ele_rate(const std::string& xml) {
  // kp=0 disables the forward reaction; AA_1 must remain at its
  // seed value (0) for the entire trajectory.
  auto r_zero = run_with(xml, {{"kp", 0.0}}, 10.0, 11, /*seed=*/1);
  check(final_value(r_zero, "AA_1") == 0, "Ele rate: set_param(kp, 0) should produce zero dimers");

  // kp = 1 is 1000x the default; with the same seed the trajectory
  // must reach more dimers than the unmodified default run.  If
  // set_param were a no-op the two runs would be bit-identical.
  auto r_high = run_with(xml, {{"kp", 1.0}}, 10.0, 11, /*seed=*/1);
  auto r_default = run_with(xml, {}, 10.0, 11, /*seed=*/1);
  check(final_value(r_high, "AA_1") > final_value(r_default, "AA_1"),
        "Ele rate: set_param(kp, 1.0) trajectory should exceed default kp=0.001");
}

void test_initial_concentration(const std::string& xml) {
  // SpeciesInit::concentration is truncated to int by the engine, so
  // for integer-valued overrides A_1(t=0) must match exactly.
  for (double a_tot : {0.0, 200.0, 5000.0}) {
    auto r = run_with(xml, {{"A_tot", a_tot}}, 10.0, 11, /*seed=*/1);
    check(initial_value(r, "A_1") == a_tot,
          "Init conc: set_param(A_tot, " + std::to_string((int)a_tot) +
              ") should give A_1(t=0)=" + std::to_string((int)a_tot));
  }
}

// ft_mm_ratelaw defaults: E0=5, S0=100, kcat=1.0, Km=10.0.
void test_mm_rate(const std::string& xml) {
  // kcat=0 zeros the MM rule's propensity; no product is ever made.
  auto r_zero = run_with(xml, {{"kcat", 0.0}}, 30.0, 31, /*seed=*/1);
  check(final_value(r_zero, "P_n") == 0, "MM rate: set_param(kcat, 0) should produce zero product");

  // Default kcat=1.0 produces substantial conversion by t=30.
  auto r_default = run_with(xml, {}, 30.0, 31, /*seed=*/1);
  check(final_value(r_default, "P_n") > 0,
        "MM rate: default kcat=1.0 should produce nonzero product");
}

void test_clear_param_overrides(const std::string& xml) {
  // Override A_tot, clear, run: should match the no-override default
  // (A_1(t=0) == 1000).
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_block_same_complex_binding(true);
  sim.set_param("A_tot", 0.0);
  sim.clear_param_overrides();
  auto r = sim.run({0.0, 10.0, 11}, /*seed=*/1);
  check(initial_value(r, "A_1") == 1000, "clear_param_overrides should restore A_tot=1000 default");
}

// Verifies that all four configuration mutators reject calls made
// while a session is active.  This locks in the symmetric "no
// mid-session mutation" contract — historically
// set_block_same_complex_binding silently mutated the underlying
// model whereas the others threw.
// Mid-session add_molecules: the BNGsim integration pattern is
// initialize → step_to (equilibrate) → add_molecules → simulate.  This
// test exercises that flow end-to-end and verifies mass conservation
// across the add boundary so a future regression in the targeted-rescan
// optimization (only affected rules rescanned on add) can't silently
// drift propensities.
void test_add_molecules_session(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/1);

  // Equilibrate at the default A_tot=1000 for a beat.
  auto seg1 = sim.simulate(0.0, 5.0, 6);
  int a1_idx_seg1 = idx_of(seg1, "A_1");
  int aa_idx_seg1 = idx_of(seg1, "AA_1");
  double a1_pre = seg1.observable_data[a1_idx_seg1].back();
  double aa_pre = seg1.observable_data[aa_idx_seg1].back();
  double total_pre = a1_pre + 2.0 * aa_pre;
  check(std::abs(total_pre - 1000.0) < 1e-9,
        "pre-add mass conservation: A_1 + 2*AA_1 should equal 1000 (got " +
            std::to_string(total_pre) + ")");

  // Mid-session bolus: add 500 free A monomers.
  sim.add_molecules("A", 500);
  check(sim.get_molecule_count("A") > 0,
        "add_molecules should produce a positive molecule count for A");

  // Continue.  Total mass should jump by exactly 500 and remain conserved.
  auto seg2 = sim.simulate(5.0, 10.0, 6);
  int a1_idx_seg2 = idx_of(seg2, "A_1");
  int aa_idx_seg2 = idx_of(seg2, "AA_1");
  double a1_post0 = seg2.observable_data[a1_idx_seg2].front();
  double aa_post0 = seg2.observable_data[aa_idx_seg2].front();
  double total_post0 = a1_post0 + 2.0 * aa_post0;
  check(std::abs(total_post0 - 1500.0) < 1e-9,
        "post-add mass conservation at segment start: should be 1500 (got " +
            std::to_string(total_post0) + ")");

  double a1_end = seg2.observable_data[a1_idx_seg2].back();
  double aa_end = seg2.observable_data[aa_idx_seg2].back();
  double total_end = a1_end + 2.0 * aa_end;
  check(std::abs(total_end - 1500.0) < 1e-9,
        "post-add mass conservation at segment end: should be 1500 (got " +
            std::to_string(total_end) + ")");

  sim.destroy_session();
}

void test_session_active_throws(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/1);

  auto throws = [](auto&& fn) {
    try {
      fn();
    } catch (const std::runtime_error&) {
      return true;
    } catch (...) {
    }
    return false;
  };

  check(throws([&] { sim.set_param("kp", 0.5); }), "set_param should throw during active session");
  check(throws([&] { sim.clear_param_overrides(); }),
        "clear_param_overrides should throw during active session");
  check(throws([&] { sim.set_molecule_limit(10000); }),
        "set_molecule_limit should throw during active session");
  check(throws([&] { sim.set_block_same_complex_binding(false); }),
        "set_block_same_complex_binding should throw during active session");

  // After destroying the session the mutators must be callable again.
  sim.destroy_session();
  check(!throws([&] { sim.set_param("kp", 0.5); }),
        "set_param should succeed after destroy_session");
  check(!throws([&] { sim.set_block_same_complex_binding(false); }),
        "set_block_same_complex_binding should succeed after destroy_session");
}

void test_setparam_to_default_is_noop(const std::string& xml) {
  // Overriding a parameter to its parsed default value must produce a
  // bit-identical trajectory to the no-override run.  Catches accidental
  // float drift or order-dependent re-resolution in apply_overrides.
  auto r_none = run_with(xml, {}, 10.0, 11, /*seed=*/1);
  auto r_atot = run_with(xml, {{"A_tot", 1000.0}}, 10.0, 11, /*seed=*/1);
  check(r_none.event_count == r_atot.event_count,
        "set_param to default should be bit-identical (event_count differs: " +
            std::to_string(r_none.event_count) + " vs " + std::to_string(r_atot.event_count) + ")");
  check(final_value(r_none, "AA_1") == final_value(r_atot, "AA_1"),
        "set_param to default should be bit-identical (AA_1 differs)");
}

// get_parameter must reflect set_param overrides immediately, BETWEEN the
// set_param call and the next run/initialize.  Pre-fix the simulator only
// synced model.parameters at apply_overrides time (driven by run/init/
// load_state) so get_parameter returned the parsed-at-load value.
void test_get_parameter_reflects_overrides(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  // Default parsed value first.
  check(sim.get_parameter("kp") == 0.001,
        "get_parameter should report parsed default before any override");

  sim.set_param("kp", 0.7);
  check(sim.get_parameter("kp") == 0.7,
        "get_parameter should reflect the most recent set_param override");

  sim.set_param("kp", 0.123);
  check(sim.get_parameter("kp") == 0.123,
        "get_parameter should reflect the latest of multiple set_param calls");

  sim.clear_param_overrides();
  check(sim.get_parameter("kp") == 0.001,
        "get_parameter should snap back to the parsed default after clear_param_overrides");
}

// set_param must reject names not declared in the loaded XML.  Pre-fix the
// override silently leaked into model.parameters as a brand-new entry; nothing
// referenced it, so a typo'd name was a silent no-op.
void test_setparam_rejects_unknown_name(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  bool threw = false;
  try {
    sim.set_param("kp_typo_no_such_param", 0.5);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "set_param on an unknown parameter name should throw");

  // After the rejected call the override map must be unchanged: a subsequent
  // run with no overrides must produce the default trajectory.
  auto r_default = run_with(xml, {}, 10.0, 11, /*seed=*/1);
  auto r_after_failed_set = sim.run({0.0, 10.0, 11}, /*seed=*/1);
  check(r_default.event_count == r_after_failed_set.event_count,
        "rejected set_param should not perturb the run that follows");
}

// Derived parameters declared in BNGL as expressions (B = 2*A, etc.) must
// re-resolve when their inputs are overridden.  Pre-fix the parameter map
// was rebuilt as base + override splat with no expression cascade, so a
// derived parameter kept its parsed-at-load numeric value regardless of
// what its inputs were set to.  Two layers exercised:
//   - A_tot = A_base * A_factor (derived, drives initial species count)
//   - kp    = kp_base * kp_mult  (derived, drives forward propensity)
void test_derived_parameter_cascade(const std::string& xml) {
  // Sanity: the parsed defaults match the fixture's hand-written values.
  rulemonkey::RuleMonkeySimulator probe(xml);
  check(probe.get_parameter("A_tot") == 500.0,
        "fixture sanity: parsed A_tot should be A_base*A_factor = 500");
  check(probe.get_parameter("kp") == 0.001,
        "fixture sanity: parsed kp should be kp_base*kp_mult = 0.001");

  // Overriding a base parameter must recompute the derived parameter that
  // references it AND propagate to the engine on the next run.
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_block_same_complex_binding(true);
  sim.set_param("A_base", 200.0); // A_tot should become 200*5 = 1000
  check(sim.get_parameter("A_tot") == 1000.0,
        "set_param(A_base, 200) should cascade to A_tot=1000 in get_parameter");

  auto r = sim.run({0.0, 10.0, 11}, /*seed=*/1);
  check(initial_value(r, "A_1") == 1000.0,
        "set_param(A_base, 200) should cascade through A_tot to A_1(t=0)=1000 in the engine");

  // Zeroing one factor must zero the derived product even when the other
  // factor is non-zero — verifies the override wins over the parsed value
  // and the cascade doesn't re-evaluate the base parameter from its expr
  // before applying the override.
  rulemonkey::RuleMonkeySimulator sim2(xml);
  sim2.set_param("A_factor", 0.0);
  check(sim2.get_parameter("A_tot") == 0.0, "set_param(A_factor, 0) should cascade to A_tot=0");
  auto r0 = sim2.run({0.0, 10.0, 11}, /*seed=*/1);
  check(initial_value(r0, "A_1") == 0.0,
        "set_param(A_factor, 0) should leave A_1(t=0)=0 (no monomers seeded)");

  // Cascade through the rate-constant chain: set_param(kp_base, 0) zeros
  // the forward propensity and the dimer count must remain 0 forever.
  rulemonkey::RuleMonkeySimulator sim3(xml);
  sim3.set_block_same_complex_binding(true);
  sim3.set_param("kp_base", 0.0);
  check(sim3.get_parameter("kp") == 0.0, "set_param(kp_base, 0) should cascade to kp=0");
  auto r_zero_kp = sim3.run({0.0, 10.0, 11}, /*seed=*/1);
  check(final_value(r_zero_kp, "AA_1") == 0,
        "set_param(kp_base, 0) should produce zero dimers (cascaded to zero kp)");

  // Overriding a derived parameter directly wins over the cascade — and
  // base parameters keep their parsed defaults.
  rulemonkey::RuleMonkeySimulator sim4(xml);
  sim4.set_param("A_tot", 42.0); // direct override of the derived param
  check(sim4.get_parameter("A_tot") == 42.0,
        "direct override of a derived parameter should win over its expression");
  check(sim4.get_parameter("A_base") == 100.0,
        "direct override of a derived parameter should not touch its inputs");
  auto r4 = sim4.run({0.0, 10.0, 11}, /*seed=*/1);
  check(initial_value(r4, "A_1") == 42.0,
        "direct override of A_tot should reach the engine via initial concentration");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::fprintf(stderr, "usage: %s <A_plus_A.xml> <ft_mm_ratelaw.xml> <derived_param_model.xml>\n",
                 argv[0]);
    return 2;
  }
  const std::string a_plus_a = argv[1];
  const std::string mm_model = argv[2];
  const std::string derived = argv[3];

  try {
    test_ele_rate(a_plus_a);
    test_initial_concentration(a_plus_a);
    test_mm_rate(mm_model);
    test_clear_param_overrides(a_plus_a);
    test_add_molecules_session(a_plus_a);
    test_session_active_throws(a_plus_a);
    test_setparam_to_default_is_noop(a_plus_a);
    test_get_parameter_reflects_overrides(a_plus_a);
    test_setparam_rejects_unknown_name(a_plus_a);
    test_derived_parameter_cascade(derived);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: set_param reaches Ele/MM/initial-conc, get_parameter is coherent, "
                       "unknown names throw, derived parameters cascade\n");
  return 0;
}
