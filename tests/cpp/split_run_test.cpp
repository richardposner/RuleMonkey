// A session advanced in calls shorter than its waiting times (issue #89).
//
// Every call ends on a waiting-time draw that lands past its t_end.  RM used
// to roll the RNG back to before that draw, so the next call drew the same
// Exp(1) value again and spent it from its own start.  A value that has
// outlasted one call is already known to exceed the hazard over it, and
// spending it again is not a fresh draw: a session cut into calls short
// against its waiting times redrew the same long wait at every call and all
// but stopped firing.  `0 -> B() 1` run to t = 100 averaged 1.82 events in
// calls of 1.0 and 0.16 in calls of 0.1, against 100, and a `100*time()`
// ramp to t = 2 averaged 182 in calls of 0.1 and 6 in calls of 0.01,
// against 200.
//
// The draw is now spent and the next call takes a fresh one, which is exact:
// the pool and the clock at t_end are the whole state, and the waiting time
// has no memory.  Each model here is one zero-order rule, so the count is
// Poisson with a known mean however the run is cut:
//
//   - rate 1 to t = 100, in simulate() calls of 1.0 and of 0.1 and step_to()
//     calls of 0.1: the mean count is 100.  The counts in disjoint calls are
//     independent Poisson(dt), so the fraction of calls that fire nothing is
//     exp(-dt) as well, which a scheme that is right on average but not
//     memoryless call by call would miss.
//   - rate 100*time() to t = 2, in simulate() calls of 0.1 and step_to() calls
//     of 0.05: the mean count is Lambda(2) = 200.
//
// Result::event_count counts the whole session, not the call, so the events
// in a call are the difference of consecutive results.  Every event here is
// a birth, which makes the count after each call equal to event_count; that
// is checked too, since the empty-call fraction depends on it.
//
// A continuation from a saved state must still take the same numbers as one
// in process, across every boundary of a chain of short calls: the RNG is
// left past the spent draw, which is the state save_state writes.

#include "rulemonkey/simulator.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
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

constexpr int kSeeds = 100;

std::uint64_t seed_of(int rep) { return std::uint64_t{1000} + static_cast<std::uint64_t>(rep); }

using Params = std::vector<std::pair<std::string, double>>;

struct Case {
  std::string name;
  std::string xml;
  Params params; // set_param overrides
  double t_end;
  double mean; // expected count at t_end
};

std::size_t obs_index(const std::vector<std::string>& names) {
  for (std::size_t i = 0; i < names.size(); ++i)
    if (names[i] == "Btot")
      return i;
  throw std::runtime_error("observable not found: Btot");
}

rulemonkey::RuleMonkeySimulator make_sim(const Case& c) {
  rulemonkey::RuleMonkeySimulator sim(c.xml);
  for (const auto& [name, value] : c.params)
    sim.set_param(name, value);
  return sim;
}

// The ensemble mean count at t_end against the Poisson mean, within four
// standard errors.
void check_mean(const std::string& label, const std::vector<double>& counts, double expect) {
  double sum = 0.0;
  for (double const c : counts)
    sum += c;
  double const mean = sum / static_cast<double>(counts.size());
  double const tol = 4.0 * std::sqrt(expect / static_cast<double>(counts.size()));
  std::fprintf(stderr, "  %-44s mean=%8.3f  expect=%8.3f  tol=%.3f\n", label.c_str(), mean, expect,
               tol);
  check(std::fabs(mean - expect) <= tol,
        label + ": mean count " + std::to_string(mean) + ", expected " + std::to_string(expect));
}

// Advance one session per seed to t_end in n_calls equal simulate() calls.
// Checks the mean count at t_end, and, when `p_empty` is given, the fraction
// of calls that fired nothing against it.
void check_simulate_split(const Case& c, int n_calls, double p_empty = -1.0) {
  std::string const label = c.name + ", " + std::to_string(n_calls) + " simulate() calls";
  std::vector<double> counts;
  long empty = 0;
  long cumulative_mismatches = 0;
  for (int rep = 0; rep < kSeeds; ++rep) {
    auto sim = make_sim(c);
    sim.initialize(seed_of(rep));
    std::size_t const bi = obs_index(sim.observable_names());
    std::int64_t events = 0;
    double count = 0.0;
    for (int i = 1; i <= n_calls; ++i) {
      double const t1 = c.t_end * i / n_calls;
      auto const r = sim.simulate(sim.current_time(), t1, 1);
      count = r.observable_data[bi].back();
      if (static_cast<double>(r.event_count) != count)
        ++cumulative_mismatches;
      if (r.event_count == events)
        ++empty;
      events = r.event_count;
    }
    counts.push_back(count);
  }
  check(cumulative_mismatches == 0,
        label + ": Result::event_count differed from the session's count after " +
            std::to_string(cumulative_mismatches) + " calls; it should count the whole session");
  check_mean(label, counts, c.mean);
  if (p_empty >= 0.0) {
    double const n = static_cast<double>(kSeeds) * n_calls;
    double const frac = static_cast<double>(empty) / n;
    double const tol = 4.0 * std::sqrt(p_empty * (1.0 - p_empty) / n);
    std::fprintf(stderr, "  %-44s empty=%8.4f  expect=%8.4f  tol=%.4f\n", label.c_str(), frac,
                 p_empty, tol);
    check(std::fabs(frac - p_empty) <= tol, label + ": " + std::to_string(frac) +
                                                " of calls fired nothing, expected " +
                                                std::to_string(p_empty));
  }
}

// The same, advanced by step_to(), which returns no Result; the count is
// read from the session at t_end.
void check_step_to_split(const Case& c, int n_calls) {
  std::vector<double> counts;
  for (int rep = 0; rep < kSeeds; ++rep) {
    auto sim = make_sim(c);
    sim.initialize(seed_of(rep));
    for (int i = 1; i <= n_calls; ++i)
      sim.step_to(c.t_end * i / n_calls);
    check(sim.current_time() == c.t_end, c.name + ": step_to did not reach t_end");
    counts.push_back(sim.get_observable_values()[obs_index(sim.observable_names())]);
  }
  check_mean(c.name + ", " + std::to_string(n_calls) + " step_to() calls", counts, c.mean);
}

// One session advanced in n_calls in process, and the same session saved and
// loaded into a fresh simulator at every boundary: each call must return the
// same samples and the same event count in both.
void check_chain_agreement(const Case& c, int n_calls) {
  std::uint64_t const seed = 11;
  auto in_process = make_sim(c);
  in_process.initialize(seed);

  // Unique per run, so two builds testing at once do not share the file.
  auto const tmp =
      std::filesystem::temp_directory_path() /
      ("rm_split_run_test." +
       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".state");
  int mismatched_calls = 0;
  for (int i = 1; i <= n_calls; ++i) {
    double const t0 = c.t_end * (i - 1) / n_calls;
    double const t1 = c.t_end * i / n_calls;
    auto const a = in_process.simulate(in_process.current_time(), t1, 2);

    auto resumed = make_sim(c);
    if (i == 1)
      resumed.initialize(seed);
    else
      resumed.load_state(tmp.string());
    auto const b = resumed.simulate(t0, t1, 2);
    resumed.save_state(tmp.string());

    if (a.time != b.time || a.observable_data != b.observable_data ||
        a.event_count != b.event_count)
      ++mismatched_calls;
  }
  std::filesystem::remove(tmp);
  check(mismatched_calls == 0, c.name + ": " + std::to_string(mismatched_calls) + " of " +
                                   std::to_string(n_calls) +
                                   " calls differ between the in-process session and the "
                                   "save/load chain");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: split_run_test <zero_order_birth_model.xml> "
                         "<time_dependent_rate_model.xml>\n");
    return 2;
  }
  Case const constant{"rate 1, t=100", argv[1], {}, 100.0, 100.0};
  Case const ramp{"rate 100*time(), t=2", argv[2], {{"a", 100.0}}, 2.0, 200.0};

  try {
    check_simulate_split(constant, 100, std::exp(-1.0));
    check_simulate_split(constant, 1000, std::exp(-0.1));
    check_step_to_split(constant, 1000);
    check_simulate_split(ramp, 20);
    check_step_to_split(ramp, 40);

    check_chain_agreement(constant, 100);
    check_chain_agreement(ramp, 20);
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
