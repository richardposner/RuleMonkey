// Global-function exposure regression test (issue #7).
//
// RuleMonkey evaluates BNGL global functions internally for rate laws.
// This test pins down the public surface that exposes their values to
// embedders:
//   Result::function_names / Result::function_data  — per-sample values
//   RuleMonkeySimulator::function_names()           — declaration order
//   RuleMonkeySimulator::get_function_values()      — live session readback
//
// Model: ft_nested_functions.xml declares three *global* functions
//   frac_P()    = S_P / ((S_U + S_P) + 1)
//   modulator() = 1.0 + 2.0 * frac_P()
//   phos_rate() = k0 * modulator()
// over observables S_U / S_P and parameter k0.  Because every exposed
// value is an exact algebraic function of observables already in the
// Result, the test cross-checks function_data against observable_data
// at every sample row rather than trusting a hard-coded number.
//
// A_plus_A.xml (no `begin functions` block) is the negative model: its
// function surface must be empty, not absent.

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

// function_names() must list exactly the three global functions in XML
// declaration order; the model has no local functions to filter out.
void test_function_names_ordering(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator const sim(xml);
  auto fns = sim.function_names();
  check(fns.size() == 3, "function_names() should list 3 global functions for ft_nested_functions");
  if (fns.size() == 3) {
    check(fns[0] == "frac_P" && fns[1] == "modulator" && fns[2] == "phos_rate",
          "function_names() should preserve XML declaration order frac_P,modulator,phos_rate");
  }
}

// Result::function_data must be column-major, parallel to observable_data:
// n_functions() rows, each n_times() long.  Every value is cross-checked
// against the observables in the same Result, so a regression in the
// settle order or the per-sample recording surfaces as a miscompare.
void test_result_function_data(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  const double k0 = sim.get_parameter("k0");

  auto r = sim.run({0.0, 5.0, 10}, /*seed=*/7);

  check(r.function_names == sim.function_names(),
        "Result.function_names should match RuleMonkeySimulator::function_names()");
  check(r.n_functions() == 3, "Result should carry 3 global functions");
  check(r.function_data.size() == r.n_functions(),
        "function_data outer length should equal n_functions()");
  for (const auto& col : r.function_data)
    check(col.size() == r.n_times(),
          "each function_data column should have one entry per sample time");

  const int s_u = idx_of(r.observable_names, "S_U");
  const int s_p = idx_of(r.observable_names, "S_P");
  const int frac = idx_of(r.function_names, "frac_P");
  const int modu = idx_of(r.function_names, "modulator");
  const int phos = idx_of(r.function_names, "phos_rate");
  check(s_u >= 0 && s_p >= 0, "S_U / S_P observables must be present");
  check(frac >= 0 && modu >= 0 && phos >= 0, "all three functions must be present");
  if (s_u < 0 || s_p < 0 || frac < 0 || modu < 0 || phos < 0)
    return;

  for (size_t t = 0; t < r.n_times(); ++t) {
    const double su = r.observable_data[s_u][t];
    const double sp = r.observable_data[s_p][t];
    const double exp_frac = sp / ((su + sp) + 1.0);
    const double exp_modu = 1.0 + (2.0 * exp_frac);
    const double exp_phos = k0 * exp_modu;
    check(std::fabs(r.function_data[frac][t] - exp_frac) < 1e-9,
          "frac_P column should equal S_P/((S_U+S_P)+1) at every sample");
    check(std::fabs(r.function_data[modu][t] - exp_modu) < 1e-9,
          "modulator column should equal 1+2*frac_P (nested-function settle)");
    check(std::fabs(r.function_data[phos][t] - exp_phos) < 1e-9,
          "phos_rate column should equal k0*modulator");
  }
}

// get_function_values() on a live session must return values consistent
// with get_observable_values() read at the same instant, and must be
// ordered like function_names().  Negative path: no active session throws.
void test_live_function_values(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  const double k0 = sim.get_parameter("k0");

  // No session yet — must throw, not return a stale/empty vector.
  bool threw = false;
  try {
    sim.get_function_values();
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "get_function_values() should throw when no session is active");

  sim.initialize(/*seed=*/3);
  sim.step_to(1.0); // advance so S_P is plausibly non-zero

  auto fvals = sim.get_function_values();
  auto ovals = sim.get_observable_values();
  check(fvals.size() == sim.function_names().size(),
        "get_function_values() length should match function_names()");

  const int s_u = idx_of(sim.observable_names(), "S_U");
  const int s_p = idx_of(sim.observable_names(), "S_P");
  const int frac = idx_of(sim.function_names(), "frac_P");
  const int phos = idx_of(sim.function_names(), "phos_rate");
  if (s_u >= 0 && s_p >= 0 && frac >= 0 && phos >= 0) {
    const double exp_frac = ovals[s_p] / ((ovals[s_u] + ovals[s_p]) + 1.0);
    const double exp_phos = k0 * (1.0 + (2.0 * exp_frac));
    check(std::fabs(fvals[frac] - exp_frac) < 1e-9,
          "live frac_P should equal S_P/((S_U+S_P)+1) against get_observable_values()");
    check(std::fabs(fvals[phos] - exp_phos) < 1e-9, "live phos_rate should equal k0*(1+2*frac_P)");
  }

  sim.destroy_session();
}

// A model with no `begin functions` block must present an empty — not
// absent — function surface: function_names() empty and Result carries
// zero function rows alongside its observables.
void test_no_functions_model(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  check(sim.function_names().empty(),
        "function_names() should be empty for a model with no functions");

  auto r = sim.run({0.0, 0.1, 2}, /*seed=*/1);
  check(r.n_functions() == 0, "Result.n_functions() should be 0 for a no-function model");
  check(r.function_data.empty(), "function_data should be empty for a no-function model");
  check(r.n_observables() > 0, "observable surface should still be populated");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: function_values_test <ft_nested_functions.xml> <A_plus_A.xml>\n");
    return 2;
  }
  const std::string fn_xml = argv[1];
  const std::string nofn_xml = argv[2];

  try {
    test_function_names_ordering(fn_xml);
    test_result_function_data(fn_xml);
    test_live_function_values(fn_xml);
    test_no_functions_model(nofn_xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: global-function exposure assertions all passed\n");
  return 0;
}
