// Regression test: set_param overrides reach the engine for the rate
// laws and initial-species concentrations that read parameter-derived
// fields baked at XML parse time.
//
// Before the apply_overrides re-resolve fix, set_param wrote into the
// override map but the engine read parse-time-baked rate_value /
// mm_kcat / mm_Km / SpeciesInit::concentration directly, so overrides
// were silent no-ops for the most common parameter shapes.
//
// Two models exercised:
//   A_plus_A.xml             — Ele rate (kp), initial conc (A_tot)
//   ft_mm_ratelaw.xml        — MM rate (kcat, Km)
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

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <A_plus_A.xml> <ft_mm_ratelaw.xml>\n", argv[0]);
    return 2;
  }
  const std::string a_plus_a = argv[1];
  const std::string mm_model = argv[2];

  try {
    test_ele_rate(a_plus_a);
    test_initial_concentration(a_plus_a);
    test_mm_rate(mm_model);
    test_clear_param_overrides(a_plus_a);
    test_add_molecules_session(a_plus_a);
    test_session_active_throws(a_plus_a);
    test_setparam_to_default_is_noop(a_plus_a);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ERROR: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: set_param reaches Ele rate, MM rate, and initial concentration\n");
  return 0;
}
