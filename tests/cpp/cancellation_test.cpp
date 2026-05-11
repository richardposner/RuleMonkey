// Cooperative-cancellation regression test for the embedder-facing
// CancelCallback hook on run() / simulate() / step_to().
//
// Pins down the four behavioral contracts that BNGsim's `timeout` kwarg
// (and any other embedder building wall-clock budgets, signal handling,
// or GUI cancel buttons on top of this hook) depends on:
//
//   1. A callback that returns `false` on its very first poll raises
//      `rulemonkey::Cancelled` out of run().  This is the "already
//      cancelled by the time we got here" boundary case — without the
//      first-iteration check, an embedder racing to cancel before
//      run() entered the loop would still pay a full simulation.
//   2. `Cancelled` inherits from `std::runtime_error`, so a generic
//      catch site still sees it.  Embedders that catch the dedicated
//      type can distinguish a cooperative stop from a genuine error.
//   3. Mid-segment cancellation of a session-based simulate() leaves
//      the session live, with current_time() somewhere inside the
//      requested window and not past it.  Embedders rely on this to
//      decide whether to resume or destroy_session() after a timeout.
//   4. A `should_continue` that always returns true (or is left empty)
//      produces the same numeric trajectory as the no-callback path —
//      installing a cancellation hook must not perturb the RNG, the
//      event schedule, or the recorded samples.
//
// The cancellation stride is intentionally not asserted here; the only
// public contract is "eventually" with sub-second granularity, and
// pinning the stride in a test would convert an internal tuning knob
// into an API guarantee.

#include "rulemonkey/simulator.hpp"

#include <atomic>
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

// (1) Pre-cancelled callback: should_continue returns false on first
// poll; run() must raise Cancelled before completing the segment.
// Verifies the engine honors a callback that's already cancelled at
// entry — otherwise a tight cancel/run race window could ignore the
// caller's request and run a full simulation.
void test_run_pre_cancelled(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  bool threw_cancelled = false;
  int poll_count = 0;
  try {
    sim.run({0.0, 100.0, 100}, /*seed=*/1, [&poll_count]() {
      ++poll_count;
      return false; // immediate cancel
    });
  } catch (const rulemonkey::Cancelled&) {
    threw_cancelled = true;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "  unexpected exception type: %s\n", e.what());
  }
  check(threw_cancelled, "run() with should_continue=false must throw rulemonkey::Cancelled");
  check(poll_count >= 1,
        "should_continue must be polled at least once before run() returns or throws");
}

// (2) Cancelled is-a runtime_error: a generic embedder that wraps RM
// calls in `catch (const std::runtime_error&)` (e.g. for logging)
// still catches the cooperative-cancellation case.  Without this
// guarantee, embedders would need a separate catch site purely for
// the new exception type — a silent breaking change.
void test_cancelled_is_runtime_error(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  bool caught_as_runtime_error = false;
  try {
    sim.run({0.0, 100.0, 100}, /*seed=*/1, []() { return false; });
  } catch (const std::runtime_error&) {
    caught_as_runtime_error = true;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "  caught wrong base type for Cancelled: %s\n", e.what());
  }
  check(caught_as_runtime_error,
        "rulemonkey::Cancelled must inherit from std::runtime_error so generic "
        "embedder catch handlers still see it");
}

// (3) Mid-session cancellation: simulate() throws Cancelled and leaves
// current_time() inside the requested window.  The exact stop time
// depends on RNG draws and the cancellation stride, so we only assert
// the boundary inequalities — not a specific value — and verify
// destroy_session()/re-initialize() works afterward.
void test_session_simulate_cancelled_midrun(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(/*seed=*/7);

  const double t_start = sim.current_time();
  const double t_end = t_start + 100.0;

  // Cancel after the callback has been polled twice; the first poll
  // happens at event_count=0 (loop entry) and lets the run proceed,
  // the second poll triggers cancellation deeper into the trajectory.
  // This exercises the "real work has happened, now stop" path rather
  // than the "stop immediately" path covered above.
  std::atomic<int> polls{0};
  bool threw_cancelled = false;
  try {
    sim.simulate(t_start, t_end, /*n_points=*/100, [&polls]() { return polls.fetch_add(1) < 1; });
  } catch (const rulemonkey::Cancelled&) {
    threw_cancelled = true;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "  unexpected exception type: %s\n", e.what());
  }
  check(threw_cancelled, "simulate() must throw Cancelled when should_continue returns false");
  check(sim.has_session(),
        "session must remain live after Cancelled so embedders can inspect or resume");
  if (sim.has_session()) {
    const double stop_t = sim.current_time();
    check(stop_t >= t_start, "current_time() after Cancelled must not regress below segment start");
    check(stop_t < t_end,
          "current_time() after Cancelled must be strictly less than the requested t_end");
  }

  // Recovery path: destroy_session, re-initialize, and run a fresh
  // segment to confirm no engine state is wedged by the prior
  // cancellation.
  sim.destroy_session();
  sim.initialize(/*seed=*/8);
  auto result = sim.simulate(0.0, 1.0, /*n_points=*/2);
  check(result.n_times() == 3,
        "post-cancellation re-initialize + simulate must produce a normal result");
}

// (4) Identity under always-true callback: a should_continue that
// returns true must not perturb the trajectory.  Tests that the
// polling logic has zero observable side effects on RNG state or
// sample timing — otherwise an embedder that opts into the hook
// "just in case" would silently get a different trajectory than one
// who didn't.
void test_always_continue_is_identity(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim_a(xml);
  rulemonkey::RuleMonkeySimulator sim_b(xml);

  auto r_a = sim_a.run({0.0, 5.0, 50}, /*seed=*/42);
  auto r_b = sim_b.run({0.0, 5.0, 50}, /*seed=*/42, []() { return true; });

  check(r_a.event_count == r_b.event_count,
        "event_count must match between no-callback and always-true-callback paths");
  check(r_a.n_times() == r_b.n_times(),
        "n_times() must match between no-callback and always-true-callback paths");

  const bool times_match = (r_a.time == r_b.time);
  check(times_match, "sample times must match bit-exact between the two paths");

  const bool obs_match = (r_a.observable_data == r_b.observable_data);
  check(obs_match, "observable trajectories must match bit-exact between the two paths");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: cancellation_test <xml_path>\n");
    return 1;
  }
  const std::string xml = argv[1];

  test_run_pre_cancelled(xml);
  test_cancelled_is_runtime_error(xml);
  test_session_simulate_cancelled_midrun(xml);
  test_always_continue_is_identity(xml);

  if (g_failures > 0) {
    std::fprintf(stderr, "cancellation_test: %d failure(s)\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "cancellation_test: all checks passed\n");
  return 0;
}
