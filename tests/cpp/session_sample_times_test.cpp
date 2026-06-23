// Explicit output-time (sample_times) regression for the STATEFUL session
// `simulate(const TimeSpec&)` overload (issue #16, session counterpart).
//
// The stateless `run(TimeSpec)` sample_times contract is covered by
// sample_times_test.cpp. This pins down the same contract on the *session*
// API — the entry point in-process hosts (e.g. BNGsim) drive for multi-segment
// network-free protocols:
//
//   1. Non-invasive + bit-identical on the session.  A session sampled at an
//      explicit *subset* of a uniform grid (built with the exact double
//      arithmetic the engine uses) is bit-identical to a uniform session run
//      at the shared instants, same seed — and event_count is unchanged.
//
//   2. Mid-protocol, arbitrary non-uniform instants.  After advancing the live
//      session with a prior segment, an explicit non-uniform list starting at
//      the current session time produces exactly one row per requested time,
//      at exactly those times — i.e. it continues the trajectory rather than
//      restarting it.
//
//   3. Start-alignment contract.  A TimeSpec whose segment start disagrees with
//      the current session time throws std::runtime_error rather than silently
//      producing a degenerate trajectory.
//
// Fixture A+A<->AA (reference/nfsim/xml/A_plus_A.xml): observables are integer
// molecule counts, exactly representable as doubles, so the bit-identical
// comparison in (1) is a true equality test.

#include "rulemonkey/simulator.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stdout, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stdout, "Usage: session_sample_times_test <xml_path>\n");
    return 2;
  }
  const std::string xml_path = argv[1];
  constexpr std::uint64_t kSeed = 20260622;

  try {
    const double t_start = 0.0;
    const double t_end = 10.0;
    const int n_points = 100;
    const double dt = (t_end - t_start) / n_points;

    std::vector<double> grid;
    grid.reserve(n_points + 1);
    for (int i = 0; i <= n_points; ++i)
      grid.push_back(t_start + (i * dt));

    // ---- Contract 1: session explicit-subset == session uniform ----------
    rulemonkey::RuleMonkeySimulator uni_sim(xml_path);
    uni_sim.initialize(kSeed);
    const auto uniform = uni_sim.simulate(t_start, t_end, n_points);
    check(uniform.n_times() == static_cast<std::size_t>(n_points) + 1,
          "uniform session run should produce n_points+1 rows, got " +
              std::to_string(uniform.n_times()));

    constexpr int kStride = 10;
    std::vector<double> subset;
    for (int i = 0; i <= n_points; i += kStride)
      subset.push_back(grid[i]);

    rulemonkey::RuleMonkeySimulator exp_sim(xml_path);
    exp_sim.initialize(kSeed);
    rulemonkey::TimeSpec explicit_ts;
    explicit_ts.t_start = t_start;
    explicit_ts.t_end = t_end;
    explicit_ts.n_points = 999; // must be IGNORED when sample_times is set
    explicit_ts.sample_times = subset;
    const auto sampled = exp_sim.simulate(explicit_ts);

    check(sampled.n_times() == subset.size(),
          "explicit session run row count should equal sample_times.size() (" +
              std::to_string(subset.size()) + "), got " + std::to_string(sampled.n_times()));
    check(sampled.event_count == uniform.event_count,
          "event_count must be identical regardless of sampling schedule (uniform=" +
              std::to_string(uniform.event_count) +
              ", sampled=" + std::to_string(sampled.event_count) + ")");

    if (sampled.n_times() == subset.size() && sampled.n_observables() == uniform.n_observables()) {
      bool all_equal = true;
      for (std::size_t k = 0; k < subset.size(); ++k) {
        const std::size_t u = k * kStride;
        if (sampled.time[k] != uniform.time[u])
          all_equal = false;
        for (std::size_t o = 0; o < sampled.n_observables(); ++o)
          if (sampled.observable_data[o][k] != uniform.observable_data[o][u])
            all_equal = false;
      }
      check(all_equal, "explicit-subset session run must be bit-identical to the uniform "
                       "session run at the shared instants");
    }

    // ---- Contract 2: mid-protocol non-uniform instants -------------------
    rulemonkey::RuleMonkeySimulator mp_sim(xml_path);
    mp_sim.initialize(kSeed);
    (void)mp_sim.simulate(0.0, 5.0, 5); // advance the live session to t=5
    const std::vector<double> tail = {5.0, 6.5, 8.0, 10.0};
    rulemonkey::TimeSpec tail_ts;
    tail_ts.sample_times = tail;
    const auto mp = mp_sim.simulate(tail_ts);
    check(mp.n_times() == tail.size(),
          "mid-protocol explicit run should produce one row per requested time (" +
              std::to_string(tail.size()) + "), got " + std::to_string(mp.n_times()));
    if (mp.n_times() == tail.size()) {
      bool times_match = true;
      for (std::size_t k = 0; k < tail.size(); ++k)
        if (mp.time[k] != tail[k])
          times_match = false;
      check(times_match, "mid-protocol Result::time must echo the requested sample_times verbatim");
    }

    // ---- Contract 3: start-alignment is validated ------------------------
    bool threw = false;
    try {
      rulemonkey::RuleMonkeySimulator bad_sim(xml_path);
      bad_sim.initialize(kSeed);
      bad_sim.simulate(0.0, 3.0, 3); // session now at t=3
      rulemonkey::TimeSpec bad_ts;
      bad_ts.sample_times = {0.0, 1.0, 2.0}; // starts at 0, not the current t=3
      (void)bad_sim.simulate(bad_ts);
    } catch (const std::exception&) {
      threw = true;
    }
    check(threw, "a TimeSpec whose segment start disagrees with current_time must throw");

  } catch (const std::exception& e) {
    std::fprintf(stdout, "FAIL: unexpected exception: %s\n", e.what());
    return 1;
  }

  if (g_failures == 0) {
    std::fprintf(stdout, "PASS: session simulate(TimeSpec) honors explicit sample_times exactly, "
                         "non-invasively, mid-protocol, and validates start alignment\n");
    return 0;
  }
  std::fprintf(stdout, "FAILED: %d assertion(s)\n", g_failures);
  return 1;
}
