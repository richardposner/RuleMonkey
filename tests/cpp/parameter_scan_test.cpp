// parameter_scan / bifurcate regression test (issue #8).
//
// Pins down the parameter-sweep API on `RuleMonkeySimulator`:
//   parameter_scan(ScanSpec, seed) -> ScanResult
//   bifurcate(ScanSpec, seed)      -> BifurcateResult
//
// What is asserted:
//   * value generation — explicit list, linear range, geometric (log)
//     range, single-point and zero-width-range degenerate cases;
//   * endpoint shape — one row per swept value, columns parallel to
//     observable_names / function_names;
//   * determinism — a fixed seed reproduces the sweep exactly, and (with
//     reset_conc) every point shares one random stream, so a flat sweep
//     (par_min == par_max) yields identical rows;
//   * monotone physical response — on A+A<->AA, raising the forward rate
//     kp must not decrease the AA endpoint;
//   * global-function columns — populated for a model with functions;
//   * bifurcate — forward and backward branches aligned to the same
//     ascending parameter axis, both with n_points rows;
//   * spec validation — unknown parameter, bad n_points, non-positive
//     log bound, and an active session each throw;
//   * no side effects — a sweep restores the swept parameter's override
//     state (a prior set_param survives; absent it the declared value is
//     left intact).
//
// argv[1] = A_plus_A.xml          (kp, km, A_tot; observables A_1, AA_1)
// argv[2] = ft_nested_functions.xml (k0; three global functions)

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

bool close(double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; }

int idx_of(const std::vector<std::string>& names, const std::string& name) {
  for (size_t i = 0; i < names.size(); ++i)
    if (names[i] == name)
      return static_cast<int>(i);
  return -1;
}

rulemonkey::ScanSpec range_spec(const std::string& param, double lo, double hi, int n, double t_end,
                                bool log_scale = false) {
  rulemonkey::ScanSpec spec;
  spec.parameter = param;
  spec.par_min = lo;
  spec.par_max = hi;
  spec.n_points = n;
  spec.log_scale = log_scale;
  spec.per_point = rulemonkey::TimeSpec{0.0, t_end, 4};
  return spec;
}

// A throwing call: `body` must raise std::runtime_error.
template <typename F> void check_throws(F&& body, const std::string& msg) {
  bool threw = false;
  bool wrong_type = false;
  try {
    body();
  } catch (const std::runtime_error&) {
    threw = true;
  } catch (...) {
    wrong_type = true; // threw, but not the documented std::runtime_error
  }
  check(threw && !wrong_type, msg);
}

// --- value generation -------------------------------------------------------

void test_linear_range(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto r = sim.parameter_scan(range_spec("kp", 0.001, 0.005, 5, 10.0));
  check(r.n_points() == 5, "linear range: 5 swept values");
  if (r.n_points() == 5) {
    const double expect[5] = {0.001, 0.002, 0.003, 0.004, 0.005};
    for (int k = 0; k < 5; ++k)
      check(close(r.param_values[k], expect[k]), "linear range: even spacing");
  }
  check(r.parameter == "kp", "ScanResult.parameter echoes the swept name");
  check(r.observable_endpoints.size() == r.n_points(), "one endpoint row per value");
  for (const auto& row : r.observable_endpoints)
    check(row.size() == r.observable_names.size(), "endpoint row width == observable count");
}

void test_log_range(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  // A_tot is an integer seed-count parameter; geometric spacing of it is
  // still a valid value sweep regardless of how the engine consumes it.
  auto r = sim.parameter_scan(range_spec("A_tot", 1.0, 100.0, 3, 5.0, /*log_scale=*/true));
  check(r.n_points() == 3, "log range: 3 swept values");
  if (r.n_points() == 3) {
    check(close(r.param_values[0], 1.0) && close(r.param_values[1], 10.0) &&
              close(r.param_values[2], 100.0),
          "log range: geometric 1,10,100 spacing");
  }
}

void test_explicit_values(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  rulemonkey::ScanSpec spec;
  spec.parameter = "kp";
  spec.values = {0.005, 0.001, 0.003}; // out of order on purpose
  spec.per_point = rulemonkey::TimeSpec{0.0, 8.0, 2};
  auto r = sim.parameter_scan(spec);
  check(r.n_points() == 3, "explicit values: 3 points");
  if (r.n_points() == 3)
    check(close(r.param_values[0], 0.005) && close(r.param_values[1], 0.001) &&
              close(r.param_values[2], 0.003),
          "explicit values preserved verbatim, including order");
}

void test_degenerate_range(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  // Single point: n_points == 1 with min == max.
  auto one = sim.parameter_scan(range_spec("kp", 0.002, 0.002, 1, 6.0));
  check(one.n_points() == 1 && close(one.param_values[0], 0.002), "single-point range");
  // Zero-width range, multiple points: par_min repeated.  With an
  // identical parameter and a shared random stream (reset_conc default),
  // every endpoint row must be bit-identical.
  auto flat = sim.parameter_scan(range_spec("kp", 0.002, 0.002, 4, 6.0));
  check(flat.n_points() == 4, "zero-width range: n_points rows");
  bool rows_identical = flat.n_points() == 4;
  for (size_t k = 1; k < flat.n_points(); ++k)
    for (size_t o = 0; o < flat.observable_names.size(); ++o)
      rows_identical = rows_identical &&
                       close(flat.observable_endpoints[k][o], flat.observable_endpoints[0][o], 0.0);
  check(rows_identical, "flat sweep: identical parameter + shared seed ⇒ identical endpoint rows");
}

// --- determinism & physics --------------------------------------------------

void test_determinism(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim_a(xml);
  rulemonkey::RuleMonkeySimulator sim_b(xml);
  auto a = sim_a.parameter_scan(range_spec("kp", 0.001, 0.005, 4, 12.0), /*seed=*/7);
  auto b = sim_b.parameter_scan(range_spec("kp", 0.001, 0.005, 4, 12.0), /*seed=*/7);
  bool same = a.n_points() == b.n_points();
  for (size_t k = 0; same && k < a.n_points(); ++k)
    for (size_t o = 0; o < a.observable_names.size(); ++o)
      same = same && close(a.observable_endpoints[k][o], b.observable_endpoints[k][o], 0.0);
  check(same, "parameter_scan is reproducible for a fixed seed");
}

void test_monotone_response(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  // A+A<->AA: raising the forward rate kp drives the equilibrium toward
  // the dimer, so the AA endpoint must be non-decreasing across the sweep.
  auto r = sim.parameter_scan(range_spec("kp", 0.0005, 0.01, 6, 40.0), /*seed=*/1);
  const int aa = idx_of(r.observable_names, "AA_1");
  check(aa >= 0, "A_plus_A exposes an AA_1 observable");
  if (aa >= 0) {
    bool monotone = true;
    for (size_t k = 1; k < r.n_points(); ++k)
      monotone =
          monotone && (r.observable_endpoints[k][aa] >= r.observable_endpoints[k - 1][aa] - 1.0);
    check(monotone, "AA endpoint is non-decreasing as kp rises");
  }
}

// --- global functions -------------------------------------------------------

void test_function_columns(const std::string& xml_with_functions) {
  rulemonkey::RuleMonkeySimulator sim(xml_with_functions);
  auto r = sim.parameter_scan(range_spec("k0", 0.01, 0.2, 4, 20.0));
  check(!r.function_names.empty(), "ft_nested_functions exposes global functions");
  check(r.function_endpoints.size() == r.n_points(), "one function-endpoint row per value");
  for (const auto& row : r.function_endpoints)
    check(row.size() == r.function_names.size(), "function row width == function count");
  // phos_rate() = k0 * modulator(); modulator() = 1 + 2*frac_P().  With
  // frac_P in [0,1], phos_rate must lie in [k0, 3*k0] at every point.
  const int pr = idx_of(r.function_names, "phos_rate");
  if (pr >= 0) {
    bool bounded = true;
    for (size_t k = 0; k < r.n_points(); ++k) {
      const double v = r.function_endpoints[k][pr];
      bounded = bounded && v >= r.param_values[k] - 1e-9 && v <= (3.0 * r.param_values[k]) + 1e-9;
    }
    check(bounded, "phos_rate endpoint stays within [k0, 3*k0]");
  }
}

// --- bifurcate --------------------------------------------------------------

void test_bifurcate_shape(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto b = sim.bifurcate(range_spec("kp", 0.001, 0.006, 5, 25.0));
  check(b.forward.n_points() == 5 && b.backward.n_points() == 5,
        "bifurcate: both branches have n_points rows");
  // Both branches must be reported on the same ascending parameter axis.
  bool aligned = b.forward.param_values.size() == b.backward.param_values.size();
  for (size_t k = 0; aligned && k < b.forward.param_values.size(); ++k)
    aligned = aligned && close(b.forward.param_values[k], b.backward.param_values[k]);
  check(aligned, "bifurcate: forward and backward share one ascending axis");
  bool ascending = true;
  for (size_t k = 1; k < b.forward.param_values.size(); ++k)
    ascending = ascending && b.forward.param_values[k] > b.forward.param_values[k - 1];
  check(ascending, "bifurcate: parameter axis is strictly ascending");
  for (const auto& row : b.backward.observable_endpoints)
    check(row.size() == b.backward.observable_names.size(),
          "bifurcate backward: endpoint row width == observable count");
}

// --- validation -------------------------------------------------------------

void test_validation(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);

  check_throws([&] { sim.parameter_scan(range_spec("not_a_param", 0.0, 1.0, 3, 5.0)); },
               "unknown parameter rejected");

  check_throws([&] { sim.parameter_scan(range_spec("kp", 0.001, 0.005, 1, 5.0)); },
               "n_points == 1 with distinct min/max rejected");

  check_throws([&] { sim.parameter_scan(range_spec("kp", 0.001, 0.005, 0, 5.0)); },
               "n_points == 0 rejected");

  check_throws([&] { sim.parameter_scan(range_spec("kp", -1.0, 5.0, 3, 5.0, /*log_scale=*/true)); },
               "log_scale with non-positive bound rejected");

  {
    rulemonkey::ScanSpec spec = range_spec("kp", 0.001, 0.005, 3, 5.0);
    spec.per_point.n_points = 0;
    check_throws([&] { sim.parameter_scan(spec); }, "per_point.n_points == 0 rejected");
  }

  check_throws(
      [&] {
        rulemonkey::RuleMonkeySimulator s(xml);
        s.initialize(1);
        s.parameter_scan(range_spec("kp", 0.001, 0.005, 3, 5.0));
      },
      "parameter_scan during an active session rejected");
}

// --- no side effects --------------------------------------------------------

void test_no_side_effects(const std::string& xml) {
  // A prior set_param override must survive a sweep over the same parameter.
  {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_param("kp", 0.1234);
    sim.parameter_scan(range_spec("kp", 0.001, 0.005, 3, 5.0));
    check(close(sim.get_parameter("kp"), 0.1234),
          "sweep restores a pre-existing set_param override");
  }
  // With no prior override, the declared value is left intact.
  {
    rulemonkey::RuleMonkeySimulator sim(xml);
    const double declared = sim.get_parameter("kp");
    sim.parameter_scan(range_spec("kp", 0.001, 0.005, 3, 5.0));
    check(close(sim.get_parameter("kp"), declared),
          "sweep leaves the declared parameter value intact when un-overridden");
    // The simulator is reusable afterwards — no lingering session.
    check(!sim.has_session(), "sweep leaves no active session behind");
    auto post = sim.run({0.0, 5.0, 3}, 3);
    check(post.n_times() == 4, "simulator still usable after a sweep");
  }
}

void test_reset_conc_carryover(const std::string& xml) {
  // reset_conc == false must run without error and still produce a full
  // sweep; the carry-over path is the same machinery bifurcate relies on.
  rulemonkey::RuleMonkeySimulator sim(xml);
  rulemonkey::ScanSpec spec = range_spec("kp", 0.001, 0.005, 4, 10.0);
  spec.reset_conc = false;
  auto r = sim.parameter_scan(spec, /*seed=*/5);
  check(r.n_points() == 4, "reset_conc=false sweep yields all points");
  check(!sim.has_session(), "reset_conc=false sweep leaves no session behind");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: parameter_scan_test <A_plus_A.xml> <ft_nested_functions.xml>\n");
    return 2;
  }
  const std::string xml = argv[1];
  const std::string xml_fn = argv[2];

  try {
    test_linear_range(xml);
    test_log_range(xml);
    test_explicit_values(xml);
    test_degenerate_range(xml);
    test_determinism(xml);
    test_monotone_response(xml);
    test_function_columns(xml_fn);
    test_bifurcate_shape(xml);
    test_validation(xml);
    test_no_side_effects(xml);
    test_reset_conc_carryover(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures > 0) {
    std::fprintf(stderr, "\n%d assertion(s) failed\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "OK: parameter_scan / bifurcate assertions all passed\n");
  return 0;
}
