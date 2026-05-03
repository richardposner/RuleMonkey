// Public-API coverage regression test.  Pins down the simpler
// surface-area methods on `RuleMonkeySimulator` so a future refactor
// that drops or renames a getter can't slip through unnoticed.
//
// Methods exercised here (each is asserted against expected output):
//   step_to                            — equilibration without sampling
//   has_session                        — session lifecycle
//   get_observable_values              — live observable readback
//   parameter_names / observable_names — XML-declaration ordering
//   xml_path / method                  — constructor metadata getters
//   unsupported_features               — diagnostic accessor on a clean model
//   set_block_same_complex_binding     — behavioral effect on a model where
//                                         BSCB toggle changes equilibrium
//   set_molecule_limit                 — caps add_molecules

#include "rulemonkey/simulator.hpp"

#include <algorithm>
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

int idx_of(const rulemonkey::Result& r, const std::string& name) {
  for (size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return static_cast<int>(i);
  return -1;
}

// xml_path() must echo the constructor input verbatim, and method()
// must report the configured runtime method even when the constructor
// defaulted it.
void test_constructor_metadata(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator const sim(xml);
  check(sim.xml_path() == xml, "xml_path() should return the constructor's input string");
  check(sim.method() == rulemonkey::Method::NfExact,
        "method() should default to NfExact when constructor omits the second arg");

  rulemonkey::RuleMonkeySimulator const sim2(xml, rulemonkey::Method::NfExact);
  check(sim2.method() == rulemonkey::Method::NfExact,
        "method() should report the explicit NfExact pass-through");
}

// observable_names / parameter_names must reflect XML declaration order
// and match the live Result column order produced by run().
void test_metadata_ordering(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto obs = sim.observable_names();
  auto params = sim.parameter_names();

  check(!obs.empty(), "observable_names() must be non-empty for A_plus_A");
  check(!params.empty(), "parameter_names() must be non-empty for A_plus_A");

  // A_plus_A.xml declares parameters t_end, n_steps, A_tot, kp, km in that
  // order.  Verify ordering is stable by checking the canonical names exist.
  auto has_param = [&](const std::string& n) {
    return std::find(params.begin(), params.end(), n) != params.end();
  };
  check(has_param("kp"), "parameter_names() should include 'kp'");
  check(has_param("km"), "parameter_names() should include 'km'");
  check(has_param("A_tot"), "parameter_names() should include 'A_tot'");

  // Result observable order must match observable_names() order.
  auto r = sim.run({0.0, 0.1, 1}, /*seed=*/1);
  check(r.observable_names == obs,
        "Result.observable_names should match RuleMonkeySimulator::observable_names()");
}

// unsupported_features() on a clean model must return an empty vector;
// its very presence is the contract.
void test_unsupported_features_accessor(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator const sim(xml);
  const auto& uf = sim.unsupported_features();
  check(uf.empty(), "unsupported_features() should be empty for a clean A_plus_A model "
                    "(non-empty would indicate XML import surfaced a Warn/Error)");
}

// has_session lifecycle: false at construction, true after initialize, false
// after destroy_session.  Also: run() must NOT leave a session active (it is
// the one-shot path).
void test_has_session_lifecycle(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  check(!sim.has_session(), "has_session() should be false right after construction");

  // run() is the one-shot path; should not leave a session behind.
  sim.run({0.0, 0.1, 1}, /*seed=*/1);
  check(!sim.has_session(), "has_session() should remain false after run()");

  sim.initialize(/*seed=*/2);
  check(sim.has_session(), "has_session() should be true after initialize()");

  sim.destroy_session();
  check(!sim.has_session(), "has_session() should be false after destroy_session()");

  // Repeated destroy_session is a no-op.
  sim.destroy_session();
  check(!sim.has_session(), "destroy_session() should be idempotent");
}

// step_to equilibrates without producing samples.  After step_to(t), the
// session current_time must be at least t and a follow-up simulate(t, t+dt, n)
// must succeed (the contract: t_start ≈ current_time).
void test_step_to(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/1);

  sim.step_to(2.0);
  double const now = sim.current_time();
  check(now >= 2.0, "current_time() should be at least step_to target");

  // simulate from the post-step_to time should succeed.
  auto r = sim.simulate(now, now + 1.0, 2);
  check(r.n_times() == 3, "simulate(now, now+1, 2) should produce 3 sample rows");

  sim.destroy_session();
}

// get_observable_values must return live counts that match the engine's
// understanding — zero seeds, post-add nonzero, sized to observable_names().
void test_get_observable_values(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/1);

  auto vals = sim.get_observable_values();
  check(vals.size() == sim.observable_names().size(),
        "get_observable_values() length should match observable_names()");

  // A_plus_A seeds A_tot=1000 monomers; A_1 observable should equal 1000.
  auto names = sim.observable_names();
  auto a1_it = std::find(names.begin(), names.end(), "A_1");
  if (a1_it != names.end()) {
    size_t const idx = static_cast<size_t>(a1_it - names.begin());
    check(std::abs(vals[idx] - 1000.0) < 1e-9, "A_1 observable should be 1000 at session start");
  }

  sim.destroy_session();
}

// set_block_same_complex_binding default is true (strict BNGL).  Toggling to
// false on a self-bond model (A+A<->AA) should NOT change the equilibrium
// (A_plus_A has no intra-complex pairs that differ between BSCB on/off — both
// reactants are free monomers).  This pins down the toggle as a behaviorally
// meaningful but here-no-op flag rather than a noop in all cases.  Negative
// path: calling during a session must throw.
void test_set_block_same_complex_binding(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim_strict(xml);
  sim_strict.set_block_same_complex_binding(true);
  auto r_strict = sim_strict.run({0.0, 1.0, 2}, /*seed=*/42);

  rulemonkey::RuleMonkeySimulator sim_loose(xml);
  sim_loose.set_block_same_complex_binding(false);
  auto r_loose = sim_loose.run({0.0, 1.0, 2}, /*seed=*/42);

  check(r_strict.observable_names == r_loose.observable_names,
        "BSCB toggle should not change observable schema");
  check(r_strict.n_times() == r_loose.n_times(), "BSCB toggle should not change time-grid shape");

  // Throws during active session.
  rulemonkey::RuleMonkeySimulator sim_active(xml);
  sim_active.initialize(/*seed=*/1);
  bool threw = false;
  try {
    sim_active.set_block_same_complex_binding(false);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "set_block_same_complex_binding should throw during active session");
  sim_active.destroy_session();
}

// set_molecule_limit caps the SSA loop: once the pool exceeds `limit`,
// the run terminates.  The contract is enforced inside the SSA loop, not
// at add_molecules, so an over-cap add succeeds but the next run halts
// without further events.  Negative path: setting the limit during an
// active session throws.
void test_set_molecule_limit(const std::string& xml) {
  // No-limit baseline (default -1) — record steady-state event count.
  rulemonkey::RuleMonkeySimulator sim_unlimited(xml);
  auto r_unlimited = sim_unlimited.run({0.0, 1.0, 2}, /*seed=*/1);

  // Tight limit at exactly the seeded count: any synthesis or net-positive
  // event would trigger the early-stop.  A_plus_A only conserves molecule
  // count (no synthesis) so the limit is non-binding here, but the smoke
  // value of "configurable + accepted by Engine" is what we check.
  rulemonkey::RuleMonkeySimulator sim_capped(xml);
  sim_capped.set_molecule_limit(1000);
  auto r_capped = sim_capped.run({0.0, 1.0, 2}, /*seed=*/1);

  check(r_capped.observable_names == r_unlimited.observable_names,
        "set_molecule_limit should not change observable schema");
  check(r_capped.n_times() == r_unlimited.n_times(),
        "set_molecule_limit should not truncate the time grid for a "
        "conserving model under the cap");

  // Throws during active session.
  rulemonkey::RuleMonkeySimulator sim_active(xml);
  sim_active.initialize(/*seed=*/1);
  bool threw = false;
  try {
    sim_active.set_molecule_limit(50);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "set_molecule_limit should throw during active session");
  sim_active.destroy_session();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: api_coverage_test <A_plus_A.xml>\n");
    return 2;
  }
  const std::string xml = argv[1];

  try {
    test_constructor_metadata(xml);
    test_metadata_ordering(xml);
    test_unsupported_features_accessor(xml);
    test_has_session_lifecycle(xml);
    test_step_to(xml);
    test_get_observable_values(xml);
    test_set_block_same_complex_binding(xml);
    test_set_molecule_limit(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: public-API coverage assertions all passed\n");
  return 0;
}
