// Replicate-scheduler regression (issue #113): pins down the ordering and
// failure contracts of detail::run_replicates, the scheduler behind
// run_ensemble, with stand-in replicates whose timing and failures the test
// controls.  Through the public API every replicate runs the same engine,
// so nothing there can say which replicate an exception came from.
//
//   * results are folded in replicate order whatever order they finish in,
//     for every thread count;
//   * the exception rethrown is the lowest-numbered replicate's, even when
//     a higher-numbered one fails first, and every replicate below it was
//     folded;
//   * a Cancelled a replicate throws itself is a failure, never dropped as
//     a replicate the scheduler stopped;
//   * a replicate's failure outranks a later user cancel.
//
// Timing only decides whether a test can tell a wrong scheduler from a
// right one; a right scheduler passes whatever the timing.

#include "ensemble.hpp"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stderr, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

// A stand-in replicate: polls `polls` times, `ms` apart, then throws
// runtime_error(err) if `err` is set, throws Cancelled itself if `cancel`
// is set, and otherwise returns a Result holding its own index.
struct Rep {
  int polls = 2;
  int ms = 1;
  const char* err = nullptr;
  bool cancel = false;
  std::atomic<bool>* failed = nullptr; // set just before `err` is thrown
};

rulemonkey::Result run_rep(int i, const rulemonkey::CancelCallback& poll, const Rep& rep) {
  for (int k = 0; k < rep.polls; ++k) {
    if (!poll())
      throw rulemonkey::Cancelled();
    std::this_thread::sleep_for(std::chrono::milliseconds(rep.ms));
  }
  if (rep.err != nullptr) {
    if (rep.failed != nullptr)
      rep.failed->store(true);
    throw std::runtime_error(rep.err);
  }
  if (rep.cancel)
    throw rulemonkey::Cancelled();
  rulemonkey::Result r;
  r.time = {0.0};
  r.observable_data = {{static_cast<double>(i)}};
  return r;
}

// Runs `reps` through the scheduler; returns the replicate indices folded,
// in fold order, and the exception's message ("" if none, "Cancelled" for
// a Cancelled).
struct Outcome {
  std::vector<int> folded;
  std::string error;
};

Outcome schedule(const std::vector<Rep>& reps, int n_threads,
                 const rulemonkey::CancelCallback& should_continue = {}) {
  Outcome out;
  try {
    rulemonkey::detail::run_replicates(
        static_cast<int>(reps.size()), n_threads, should_continue,
        [&](int i, const rulemonkey::CancelCallback& poll) {
          return run_rep(i, poll, reps[static_cast<std::size_t>(i)]);
        },
        [&](const rulemonkey::Result& r) {
          out.folded.push_back(static_cast<int>(r.observable_data[0][0]));
        });
  } catch (const rulemonkey::Cancelled&) {
    out.error = "Cancelled";
  } catch (const std::exception& e) {
    out.error = e.what();
  }
  return out;
}

std::vector<int> iota(int n) {
  std::vector<int> v;
  v.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i)
    v.push_back(i);
  return v;
}

std::string label(const std::string& what, int n_threads) {
  return what + " (n_threads=" + std::to_string(n_threads) + ")";
}

// Later replicates finish first; folding must still go 0, 1, 2, ...
void test_fold_order() {
  std::vector<Rep> reps(12);
  for (std::size_t i = 0; i < reps.size(); ++i)
    reps[i].polls = static_cast<int>(reps.size() - i);
  for (int const n_threads : {1, 3, 12, 0}) {
    auto const o = schedule(reps, n_threads);
    check(o.error.empty() && o.folded == iota(12),
          label("results must fold in replicate order", n_threads));
  }
}

// Replicate 5 fails at once, replicate 2 only after ~40 ms, so with several
// threads 5's failure comes first in time; 2's must be the one rethrown.
void test_lowest_failure_wins() {
  std::vector<Rep> reps(8);
  reps[5].polls = 0;
  reps[5].err = "rep 5";
  reps[2].polls = 40;
  reps[2].err = "rep 2";
  for (int const n_threads : {1, 2, 4, 8}) {
    auto const o = schedule(reps, n_threads);
    check(o.error == "rep 2",
          label("the lowest-numbered failure must be rethrown, got '" + o.error + "'", n_threads));
    check(o.folded == iota(2),
          label("every replicate below the lowest failure must be folded", n_threads));
  }
}

// Replicate 3 throws Cancelled itself, with nothing having asked it to stop.
void test_replicate_cancelled_is_a_failure() {
  std::vector<Rep> reps(8);
  reps[3].cancel = true;
  for (int const n_threads : {1, 4}) {
    auto const o = schedule(reps, n_threads);
    check(o.error == "Cancelled", label("a replicate's own Cancelled must be rethrown", n_threads));
    check(o.folded == iota(3), label("replicates 0-2 must be folded", n_threads));
  }
}

// The user cancels only once replicate 3 has failed; its error, not
// Cancelled, must come back.
void test_failure_outranks_cancel() {
  for (int const n_threads : {1, 4}) {
    std::atomic<bool> failed{false};
    std::vector<Rep> reps(8);
    for (auto& r : reps)
      r.polls = 20;
    reps[3].polls = 0;
    reps[3].err = "rep 3";
    reps[3].failed = &failed;
    auto const o = schedule(reps, n_threads, [&] { return !failed.load(); });
    check(o.error == "rep 3",
          label("a replicate's failure must outrank a later cancel, got '" + o.error + "'",
                n_threads));
  }
}

} // namespace

int main() {
  test_fold_order();
  test_lowest_failure_wins();
  test_replicate_cancelled_is_a_failure();
  test_failure_outranks_cancel();

  if (g_failures > 0) {
    std::fprintf(stderr, "ensemble_scheduler_test: %d failure(s)\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "ensemble_scheduler_test: all checks passed\n");
  return 0;
}
