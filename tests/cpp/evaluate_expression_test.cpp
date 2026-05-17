// evaluate_expression session-API regression test (issue #9 §1).
//
// RuleMonkeySimulator::evaluate_expression(expr, extra) compiles an
// arbitrary BNGL expression and evaluates it against the active
// session's current state.  This test pins down its symbol scope and
// error surface:
//   - declared parameters resolve to get_parameter() values
//   - observables resolve to get_observable_values() values
//   - global functions resolve to get_function_values() values
//   - the `time()` builtin and bare `t` resolve to current_time()
//   - the `extra` map layers name=value bindings, shadowing model
//     symbols on a clash, without mutating session state
//   - a compile failure and a no-session call both throw
//
// Model: ft_nested_functions.xml declares parameter k0, observables
// S_U / S_P, and global functions frac_P() / modulator() / phos_rate().
// Every expected value is cross-checked against the session's own
// accessors at the same instant, so no hard-coded numbers are needed.

#include "rulemonkey/simulator.hpp"

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

int idx_of(const std::vector<std::string>& names, const std::string& name) {
  for (size_t i = 0; i < names.size(); ++i)
    if (names[i] == name)
      return static_cast<int>(i);
  return -1;
}

// No active session: evaluate_expression must throw rather than return
// a stale or zero value.
void test_no_session_throws(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  bool threw = false;
  try {
    sim.evaluate_expression("k0");
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "evaluate_expression() should throw when no session is active");
}

// Parameters, observables, and global functions must all resolve to the
// same values the dedicated session accessors report at the same instant.
void test_symbol_scope(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  const double k0 = sim.get_parameter("k0");

  sim.initialize(/*seed=*/3);
  sim.step_to(1.0); // advance so S_P is plausibly non-zero

  // Parameter resolves to get_parameter().
  check(std::fabs(sim.evaluate_expression("k0") - k0) < 1e-12,
        R"(evaluate_expression("k0") should equal get_parameter("k0"))");

  // Observables resolve to get_observable_values().
  auto ovals = sim.get_observable_values();
  const int s_u = idx_of(sim.observable_names(), "S_U");
  const int s_p = idx_of(sim.observable_names(), "S_P");
  check(s_u >= 0 && s_p >= 0, "S_U / S_P observables must be present");
  if (s_u >= 0 && s_p >= 0) {
    check(std::fabs(sim.evaluate_expression("S_U") - ovals[s_u]) < 1e-9,
          R"(evaluate_expression("S_U") should equal get_observable_values()[S_U])");
    // A compound expression mixing observables and a parameter.
    const double expr = sim.evaluate_expression("k0 * (S_U + S_P)");
    check(std::fabs(expr - (k0 * (ovals[s_u] + ovals[s_p]))) < 1e-9,
          "evaluate_expression should evaluate compound parameter/observable arithmetic");
  }

  // Global functions resolve to get_function_values().
  auto fvals = sim.get_function_values();
  const int phos = idx_of(sim.function_names(), "phos_rate");
  check(phos >= 0, "phos_rate function must be present");
  if (phos >= 0) {
    check(std::fabs(sim.evaluate_expression("phos_rate") - fvals[phos]) < 1e-9,
          R"(evaluate_expression("phos_rate") should equal get_function_values()[phos_rate])");
  }

  sim.destroy_session();
}

// The `time()` builtin and bare `t` must both resolve to the session's
// current logical time.
void test_time_symbol(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/5);
  sim.step_to(2.5);
  const double now = sim.current_time();

  check(std::fabs(sim.evaluate_expression("time()") - now) < 1e-9,
        R"msg(evaluate_expression("time()") should equal current_time())msg");
  check(std::fabs(sim.evaluate_expression("t") - now) < 1e-9,
        R"(evaluate_expression("t") should equal current_time())");
  sim.destroy_session();
}

// The `extra` map supplies extra bindings and shadows model symbols on a
// name clash; it must not mutate session state.
void test_extra_overrides(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  const double k0 = sim.get_parameter("k0");
  sim.initialize(/*seed=*/9);
  sim.step_to(1.0);

  // A free name unknown to the model resolves from `extra`.
  check(std::fabs(sim.evaluate_expression("2 * q", {{"q", 21.0}}) - 42.0) < 1e-12,
        "evaluate_expression should resolve a free symbol from the extra map");

  // `extra` shadows a model parameter for this call only.
  const double shadowed = sim.evaluate_expression("k0", {{"k0", k0 + 100.0}});
  check(std::fabs(shadowed - (k0 + 100.0)) < 1e-9,
        "an extra binding should shadow a model parameter of the same name");
  check(std::fabs(sim.get_parameter("k0") - k0) < 1e-12,
        "evaluate_expression with extra overrides must not mutate session parameters");
  check(std::fabs(sim.evaluate_expression("k0") - k0) < 1e-12,
        "a subsequent evaluate_expression must see the unmodified parameter");
  sim.destroy_session();
}

// A syntactically invalid or unresolved expression must throw.
void test_compile_error_throws(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/1);

  bool threw_syntax = false;
  try {
    sim.evaluate_expression("k0 +");
  } catch (const std::runtime_error&) {
    threw_syntax = true;
  }
  check(threw_syntax, "evaluate_expression should throw on a syntax error");

  bool threw_unknown = false;
  try {
    sim.evaluate_expression("no_such_identifier_xyz");
  } catch (const std::runtime_error&) {
    threw_unknown = true;
  }
  check(threw_unknown, "evaluate_expression should throw on an unresolved identifier");
  sim.destroy_session();
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: evaluate_expression_test <ft_nested_functions.xml>\n");
    return 2;
  }
  const std::string xml = argv[1];

  try {
    test_no_session_throws(xml);
    test_symbol_scope(xml);
    test_time_symbol(xml);
    test_extra_overrides(xml);
    test_compile_error_throws(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: evaluate_expression assertions all passed\n");
  return 0;
}
