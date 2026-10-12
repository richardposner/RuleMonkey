// Rate laws that read the clock (issue #86).
//
// Between two events the pool is fixed, but a rate that reads `time` still
// moves.  RM used to hold every propensity at its value from the last event,
// which is exact only for a rate that does not move between events:
//
//   - `k*time()` is 0 at t = 0, so the run took the start as an absorbing
//     state and recorded B = 0 to the end, in every seed;
//   - `k*time() + 0.1` is merely small at t = 0, so the first waiting time
//     was drawn at rate 0.1 and outlasted the run in most seeds;
//   - `k*(time()-1)^2` passes through 0 at t = 1, and a run whose event
//     landed near there held the near-zero rate past the end.
//
// A pulse 0.01 wide at t = 1 pins the other half of the fix: no integration
// panel is longer than a thousandth of the run.  A rate that reads 0 at every
// node of a panel has no error to measure, so a step left to grow on the
// quiet stretch before the pulse would step over it whole.
//
// RM now integrates the clock-driven propensity between events and fires
// where the integral reaches the Exp(1) draw.  Every shape here is one
// zero-order rule with nothing else in the model to refresh it, so B(t) is
// Poisson with mean Lambda(t) = integral of the rate from 0 to t, and the
// expected trajectory is known exactly at every sample.
//
// Each shape puts Lambda(2) at 200.  The count at t = 2 is then the number
// of partial sums of the Exp(1) draws that stay below 200, a function of the
// draws alone, so every shape that consumes the RNG the same way must report
// the same B(2) in the same seed -- the time-change theorem, and a check on
// the integration in each realization rather than on average.

#include "rulemonkey/simulator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
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

constexpr int kReps = 100;
constexpr double kTEnd = 2.0;
constexpr int kPoints = 4; // samples at 0, 0.5, 1, 1.5, 2

std::size_t obs_index(const rulemonkey::Result& r, const std::string& name) {
  for (std::size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == name)
      return i;
  throw std::runtime_error("observable not found: " + name);
}

struct Shape {
  std::string name;
  std::string xml;
  std::string obs;
  std::vector<std::pair<std::string, double>> params; // set_param overrides
  std::function<double(double)> lambda;               // integrated rate
};

// One run per seed; returns the observable's trajectory for each.
std::vector<std::vector<double>> run_ensemble(const Shape& s) {
  std::vector<std::vector<double>> out;
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(s.xml);
    for (const auto& [name, value] : s.params)
      sim.set_param(name, value);
    auto const r =
        sim.run({0.0, kTEnd, kPoints}, std::uint64_t{1000} + static_cast<std::uint64_t>(rep));
    out.push_back(r.observable_data[obs_index(r, s.obs)]);
  }
  return out;
}

// The ensemble mean at each sample against Lambda(t), within four standard
// errors of a Poisson mean.  A sample where Lambda is 0 must read 0 in every
// seed.
void check_mean_trajectory(const Shape& s, const std::vector<std::vector<double>>& runs) {
  for (int k = 0; k <= kPoints; ++k) {
    double const t = kTEnd * k / kPoints;
    double sum = 0.0;
    for (const auto& run : runs)
      sum += run[k];
    double const mean = sum / static_cast<double>(runs.size());
    double const expect = s.lambda(t);
    double const tol = 4.0 * std::sqrt(expect / static_cast<double>(runs.size()));
    std::fprintf(stderr, "  %-6s t=%.1f  mean=%8.3f  Lambda=%8.3f  tol=%.3f\n", s.name.c_str(), t,
                 mean, expect, tol);
    check(std::fabs(mean - expect) <= tol, s.name + ": mean B(" + std::to_string(t) +
                                               ") = " + std::to_string(mean) +
                                               ", expected Lambda = " + std::to_string(expect));
  }
}

// A function of `time` reads the sample's own time, not the time of the last
// event before it.  Under the old recorder `clock() = time()` read the event
// time, and 0 at every sample of a run that had not fired yet.
void check_clock_output(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto const r = sim.run({0.0, kTEnd, kPoints}, 7);
  std::size_t ci = r.function_names.size();
  for (std::size_t i = 0; i < r.function_names.size(); ++i)
    if (r.function_names[i] == "clock")
      ci = i;
  check(ci < r.function_names.size(), "clock(): output function not found");
  if (ci == r.function_names.size())
    return;
  for (std::size_t k = 0; k < r.time.size(); ++k)
    check(r.function_data[ci][k] == r.time[k], "clock() at sample t=" + std::to_string(r.time[k]) +
                                                   " read " +
                                                   std::to_string(r.function_data[ci][k]));
}

// Continuing a session in process and continuing it from a saved state are
// the same continuation.  In process, the old loop resumed from whatever
// rate the last event had left, an instant before the save; load_state
// re-evaluates every rate at the saved time.  For a rate that reads the
// clock those differ, and so did the two runs.  This pins that the two
// agree, not that either reproduces an unsplit run: the draw that outlasts
// the first call is spent, and both continuations take a fresh one from
// t = 1 (issue #89).  split_run_test checks that a split run's count has the
// right distribution, and that the two agree across many boundaries.
void check_continuation(const std::string& xml) {
  const std::uint64_t seed = 11;
  auto configure = [](rulemonkey::RuleMonkeySimulator& sim) {
    sim.set_param("a", 100.0);
    sim.set_param("b", 0.1);
  };

  rulemonkey::RuleMonkeySimulator in_process(xml);
  configure(in_process);
  in_process.initialize(seed);
  auto const head = in_process.simulate(0.0, 1.0, 2);
  auto const tail_a = in_process.simulate(1.0, kTEnd, 2);

  // Unique per run, so two builds testing at once do not share the file.
  auto const tmp =
      std::filesystem::temp_directory_path() /
      ("rm_time_dependent_rate_test." +
       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".state");
  {
    rulemonkey::RuleMonkeySimulator first(xml);
    configure(first);
    first.initialize(seed);
    first.simulate(0.0, 1.0, 2);
    first.save_state(tmp.string());
  }
  rulemonkey::RuleMonkeySimulator resumed(xml);
  configure(resumed);
  resumed.load_state(tmp.string());
  auto const tail_b = resumed.simulate(1.0, kTEnd, 2);
  std::filesystem::remove(tmp);

  // Result::event_count counts the whole session, and the state file
  // carries it, so both tails report the count at t = 2.
  check(tail_a.event_count == tail_b.event_count,
        "continuation: in-process session ended on " + std::to_string(tail_a.event_count) +
            " events, the save/load one on " + std::to_string(tail_b.event_count));
  check(tail_a.observable_data == tail_b.observable_data,
        "continuation: in-process and save/load tails differ");
  check(tail_a.event_count > head.event_count, "continuation: the tail fired no events");
}

// A short call late in a run.  A thousandth of a 1e-8 call at t = 1e6 is
// under an ulp of t, so a step that long rounded back to t: the clock never
// moved, and the search gave up after a million empty steps -- or, with a
// rate large enough to fire, kept firing at a frozen time.
void check_late_short_call(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_param("a", 0.0);
  sim.set_param("d", 1e-3); // 1e-3 from t = 1 on
  sim.initialize(5);
  try {
    auto const head = sim.simulate(0.0, 1e6, 1);
    double const t0 = sim.current_time();
    auto const tail = sim.simulate(t0, t0 + 1e-8, 1);
    // Result::event_count counts the whole session, not the call.
    int64_t const fired = tail.event_count - head.event_count;
    check(sim.current_time() == t0 + 1e-8, "late short call: the clock did not reach t_end");
    check(fired <= 1,
          "late short call: fired " + std::to_string(fired) + " events in 1e-8 at rate 1e-3");
  } catch (const std::exception& e) {
    check(false, std::string("late short call threw: ") + e.what());
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::fprintf(stderr, "Usage: time_dependent_rate_test <time_dependent_rate_model.xml> "
                         "<tfun_time_rate_model.xml> <local_time_rate_model.xml>\n");
    return 2;
  }
  const std::string shapes_xml = argv[1];
  const std::string tfun_xml = argv[2];
  const std::string local_xml = argv[3];

  auto ramp = [](double t) { return 50.0 * t * t; };
  std::vector<Shape> const shapes = {
      {"ramp", shapes_xml, "Btot", {{"a", 100.0}}, ramp},
      {"small",
       shapes_xml,
       "Btot",
       {{"a", 100.0}, {"b", 0.1}},
       [](double t) { return (50.0 * t * t) + (0.1 * t); }},
      {"dip",
       shapes_xml,
       "Btot",
       {{"a", 0.0}, {"c", 300.0}},
       [](double t) { return 100.0 * (std::pow(t - 1.0, 3) + 1.0); }},
      {"step",
       shapes_xml,
       "Btot",
       {{"a", 0.0}, {"d", 200.0}},
       [](double t) { return t > 1.0 ? 200.0 * (t - 1.0) : 0.0; }},
      {"pulse",
       shapes_xml,
       "Btot",
       {{"a", 0.0}, {"e", 20000.0}},
       [](double t) { return 20000.0 * std::clamp(t - 1.0, 0.0, 0.01); }},
      {"tfun", tfun_xml, "Btot", {}, ramp},
      {"local", local_xml, "Ltot", {}, ramp},
  };

  try {
    std::vector<std::vector<std::vector<double>>> runs;
    for (const auto& s : shapes) {
      runs.push_back(run_ensemble(s));
      check_mean_trajectory(s, runs.back());
    }

    // ramp, dip, step, pulse and tfun: Lambda(2) = 200 and one rule drawing
    // the same numbers, so the same B(2) in every seed.  Allow 1 for a
    // partial sum that lands within the integration tolerance of 200.
    for (std::size_t const si : {std::size_t{2}, std::size_t{3}, std::size_t{4}, std::size_t{5}}) {
      int worst = 0;
      for (int rep = 0; rep < kReps; ++rep) {
        int const d = static_cast<int>(std::fabs(runs[si][rep].back() - runs[0][rep].back()));
        worst = std::max(worst, d);
      }
      check(worst <= 1, shapes[si].name + ": B(2) differs from ramp's by up to " +
                            std::to_string(worst) + " in one seed");
    }

    check_clock_output(shapes_xml);
    check_continuation(shapes_xml);
    check_late_short_call(shapes_xml);
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
