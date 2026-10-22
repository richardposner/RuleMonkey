// The seed state sits at t_start (issue #98).
//
// run(ts) with ts.t_start > 0 used to start the clock at 0, simulate
// [0, t_start], and record its first row at t_start from whatever state that
// stretch left.  A run asked for [5, 7] therefore ran 7 time units, and its
// first row was not the seed state.  A session could not start anywhere but
// 0: initialize() put the clock there, and simulate() refuses a t_start that
// disagrees with the session.  Once a rate reads time() (issue #86), where
// the clock starts decides the answer.
//
// The seed state now sits at t_start on the model's own clock, as it does
// for BNG's network simulate (run_network -i t_start), so time() reads
// absolute model time from the first row on.  Sessions take the same origin
// from initialize(seed, t_start), and parameter_scan / bifurcate start every
// point, or the chain's first, there too.
//
// The fixture is one zero-order rule with nothing else to refresh it, so
// B(t) is Poisson with mean Lambda(t) = integral of the rate from t_start to
// t, known exactly at every sample.  Two shapes: a flat rate of 100, and a
// ramp 100*time(), whose Lambda depends on where the clock starts.

#include "rulemonkey/simulator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
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

constexpr int kReps = 100;
constexpr double kTStart = 5.0;
constexpr double kTEnd = 7.0;
constexpr int kPoints = 2; // samples at 5, 6, 7

std::uint64_t seed_of(int rep) { return std::uint64_t{2000} + static_cast<std::uint64_t>(rep); }

std::size_t index_of(const std::vector<std::string>& names, const std::string& name) {
  for (std::size_t i = 0; i < names.size(); ++i)
    if (names[i] == name)
      return i;
  throw std::runtime_error("name not found: " + name);
}

using Params = std::vector<std::pair<std::string, double>>;

const Params kFlat = {{"a", 0.0}, {"b", 100.0}};
const Params kRamp = {{"a", 100.0}};

double lambda_flat(double t) { return 100.0 * (t - kTStart); }
double lambda_ramp(double t) { return 50.0 * ((t * t) - (kTStart * kTStart)); }

void configure(rulemonkey::RuleMonkeySimulator& sim, const Params& params) {
  for (const auto& [name, value] : params)
    sim.set_param(name, value);
}

// A Poisson ensemble mean against Lambda, within four standard errors.  A
// sample where Lambda is 0 must read 0 in every seed.
void check_mean(const std::string& what, const std::vector<double>& values, double expect) {
  double sum = 0.0;
  for (double const v : values)
    sum += v;
  double const mean = sum / static_cast<double>(values.size());
  double const tol = 4.0 * std::sqrt(expect / static_cast<double>(values.size()));
  std::fprintf(stderr, "  %-34s mean=%8.3f  Lambda=%8.3f  tol=%.3f\n", what.c_str(), mean, expect,
               tol);
  check(std::fabs(mean - expect) <= tol,
        what + ": mean " + std::to_string(mean) + ", expected Lambda = " + std::to_string(expect));
}

// Btot at each sample of a one-shot run over [5, 7], one row per seed.
std::vector<std::vector<double>> run_window(const std::string& xml, const Params& params) {
  std::vector<std::vector<double>> out;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    configure(sim, params);
    auto const r = sim.run({kTStart, kTEnd, kPoints}, seed_of(rep));
    out.push_back(r.observable_data[index_of(r.observable_names, "Btot")]);
  }
  return out;
}

void check_window(const std::string& name, const std::vector<std::vector<double>>& runs,
                  const std::function<double(double)>& lambda) {
  for (int k = 0; k <= kPoints; ++k) {
    double const t = kTStart + ((kTEnd - kTStart) * k / kPoints);
    std::vector<double> col;
    col.reserve(runs.size());
    for (const auto& run : runs)
      col.push_back(run[k]);
    check_mean(name + " run B(" + std::to_string(t) + ")", col, lambda(t));
  }
}

// The flat rate over [5, 7] and over [0, 2] spends the same draws on the same
// two time units, so the two runs end on the same count in every seed.  Under
// the old pre-history the [5, 7] run spent 7 units.  The count is exact: the
// events go on the time elapsed since t_start (issue #110), which is the same
// sum of the same draws from either start.  It used to allow 1 for a partial
// sum that rounded to the other side of t_end on the shifted clock.
void check_flat_matches_shifted_run(const std::string& xml,
                                    const std::vector<std::vector<double>>& late) {
  int worst = 0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    configure(sim, kFlat);
    auto const r = sim.run({0.0, kTEnd - kTStart, kPoints}, seed_of(rep));
    double const early = r.observable_data[index_of(r.observable_names, "Btot")].back();
    worst = std::max(worst, static_cast<int>(std::fabs(late[rep].back() - early)));
  }
  check(worst == 0, "flat over [5, 7] and over [0, 2] differ by up to " + std::to_string(worst) +
                        " in one seed");
}

// clock() reads each sample's own time on the model's clock.
void check_clock_column(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  configure(sim, kRamp);
  auto const r = sim.run({kTStart, kTEnd, kPoints}, 3);
  std::size_t const ci = index_of(r.function_names, "clock");
  for (std::size_t k = 0; k < r.time.size(); ++k)
    check(r.function_data[ci][k] == r.time[k], "clock() at sample t=" + std::to_string(r.time[k]) +
                                                   " read " +
                                                   std::to_string(r.function_data[ci][k]));
}

// TimeSpec's contract: a sample time below t_start is recorded at the
// initial state.  With the pre-history it recorded the state the run had
// reached by then.
void check_sample_before_t_start(const std::string& xml) {
  int nonzero = 0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    configure(sim, kFlat);
    rulemonkey::TimeSpec ts;
    ts.t_start = kTStart;
    ts.t_end = kTEnd;
    ts.sample_times = {1.0, kTStart, kTEnd};
    auto const r = sim.run(ts, seed_of(rep));
    auto const& b = r.observable_data[index_of(r.observable_names, "Btot")];
    if (b[0] != 0.0 || b[1] != 0.0)
      ++nonzero;
  }
  check(nonzero == 0, "a sample at or before t_start read other than the seed state in " +
                          std::to_string(nonzero) + " of " + std::to_string(kReps) + " seeds");
}

// A reversed window ran nothing and labelled the state with every sample
// time; run() now refuses it, as simulate() and parameter_scan do.
void check_reversed_window(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  bool threw = false;
  try {
    sim.run({kTEnd, kTStart, kPoints}, 1);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "run() over [7, 5] should throw");
}

// initialize(seed, t_start) puts a session's seed state at t_start, and a
// session started there is the one-shot run, draw for draw.
void check_session(const std::string& xml) {
  {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.initialize(1);
    check(sim.current_time() == 0.0, "initialize(seed) should start the session at 0, not " +
                                         std::to_string(sim.current_time()));
  }

  std::vector<std::vector<double>> session_runs;
  int differ = 0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    configure(sim, kRamp);
    sim.initialize(seed_of(rep), kTStart);
    check(sim.current_time() == kTStart,
          "initialize(seed, 5) left the session at " + std::to_string(sim.current_time()));
    auto const r = sim.simulate(kTStart, kTEnd, kPoints);
    session_runs.push_back(r.observable_data[index_of(r.observable_names, "Btot")]);

    rulemonkey::RuleMonkeySimulator one_shot(xml);
    configure(one_shot, kRamp);
    auto const o = one_shot.run({kTStart, kTEnd, kPoints}, seed_of(rep));
    if (o.observable_data != r.observable_data || o.event_count != r.event_count)
      ++differ;
  }
  check(differ == 0, "a session started at 5 and the one-shot run over [5, 7] differ in " +
                         std::to_string(differ) + " of " + std::to_string(kReps) + " seeds");
  check_window("ramp session", session_runs, lambda_ramp);

  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(1, kTStart);
  bool threw = false;
  try {
    sim.simulate(0.0, kTEnd, kPoints);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  check(threw, "simulate(0, ...) on a session started at 5 should throw");

  // A failed initialize leaves no session behind.  One left uninitialized
  // would initialize itself at the next simulate()'s t_start instead.
  rulemonkey::RuleMonkeySimulator bad(xml);
  bad.initialize(1);
  bool threw_nan = false;
  try {
    bad.initialize(1, std::nan(""));
  } catch (const std::runtime_error&) {
    threw_nan = true;
  }
  check(threw_nan, "initialize(seed, NaN) should throw");
  check(!bad.has_session(), "a failed initialize left a session active");
}

// parameter_scan starts each fresh point at per_point.t_start, and a
// carry-over chain starts its first point there and continues on the same
// clock.  Under the pre-history a fresh point ran [0, 7] and the chain's
// first point ran [0, 2]: the two disagreed about the same window.
void check_scans(const std::string& xml) {
  rulemonkey::ScanSpec spec;
  spec.parameter = "a";
  spec.per_point = {kTStart, kTEnd, kPoints};

  std::vector<double> fresh;
  std::vector<double> chain0;
  std::vector<double> chain1;
  int differ = 0;
  int bad_clock = 0;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    spec.values = {100.0};
    spec.reset_conc = true;
    auto const f = sim.parameter_scan(spec, seed_of(rep));
    std::size_t const bi = index_of(f.observable_names, "Btot");
    std::size_t const ci = index_of(f.function_names, "clock");
    fresh.push_back(f.observable_endpoints[0][bi]);

    spec.values = {100.0, 100.0};
    spec.reset_conc = false;
    auto const c = sim.parameter_scan(spec, seed_of(rep));
    chain0.push_back(c.observable_endpoints[0][bi]);
    chain1.push_back(c.observable_endpoints[1][bi]);
    if (c.observable_endpoints[0][bi] != f.observable_endpoints[0][bi])
      ++differ;
    if (f.function_endpoints[0][ci] != kTEnd || c.function_endpoints[0][ci] != kTEnd ||
        c.function_endpoints[1][ci] != kTEnd + (kTEnd - kTStart))
      ++bad_clock;
  }
  check_mean("parameter_scan fresh B(7)", fresh, lambda_ramp(kTEnd));
  check_mean("parameter_scan chain point 0 B(7)", chain0, lambda_ramp(kTEnd));
  check_mean("parameter_scan chain point 1 B(9)", chain1, lambda_ramp(kTEnd + (kTEnd - kTStart)));
  check(differ == 0, "a fresh point and the chain's first point differ in " +
                         std::to_string(differ) + " of " + std::to_string(kReps) + " seeds");
  check(bad_clock == 0, "clock() endpoints were not 7, 7 and 9 in " + std::to_string(bad_clock) +
                            " of " + std::to_string(kReps) + " seeds");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: clock_origin_test <time_dependent_rate_model.xml>\n");
    return 2;
  }
  const std::string xml = argv[1];

  try {
    auto const flat = run_window(xml, kFlat);
    check_window("flat", flat, lambda_flat);
    int nonzero_first = 0;
    for (const auto& run : flat)
      if (run[0] != 0.0)
        ++nonzero_first;
    check(nonzero_first == 0, "the first row is not the seed state in " +
                                  std::to_string(nonzero_first) + " of " + std::to_string(kReps) +
                                  " seeds");
    check_flat_matches_shifted_run(xml, flat);

    check_window("ramp", run_window(xml, kRamp), lambda_ramp);
    check_clock_column(xml);
    check_sample_before_t_start(xml);
    check_reversed_window(xml);
    check_session(xml);
    check_scans(xml);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "EXCEPTION: %s\n", e.what());
    return 2;
  }

  if (g_failures != 0) {
    std::fprintf(stderr, "%d failure(s)\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "PASS\n");
  return 0;
}
