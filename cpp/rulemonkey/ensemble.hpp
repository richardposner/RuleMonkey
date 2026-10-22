#pragma once

// Replicate scheduler behind RuleMonkeySimulator::run_ensemble (issue #113).
// It knows nothing of the engine, so its ordering and failure contracts are
// unit-tested with a stand-in replicate (ensemble_scheduler_test).

#include "rulemonkey/types.hpp"

#include <algorithm>
#include <atomic>
#include <exception>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace rulemonkey::detail {

// Runs replicates 0 .. n_reps-1 on up to `n_threads` workers (0 = one per
// hardware thread, 1 = the calling thread only, never more than n_reps).
// `run_one(i, poll)` runs replicate i and returns its Result; it must pass
// `poll` on as its cancel callback, which throws Cancelled out of it.
// `fold(Result&&)` is called once per replicate, in replicate order whatever
// order the workers finish in, never by two threads at once.
//
// Replicates numbered above the lowest failing one are stopped, every
// replicate below it runs to completion, and the lowest-numbered failure is
// rethrown, so the error is the same whatever the thread count or timing.
// `should_continue` is polled one call at a time under a lock; a false from
// it stops every replicate and throws Cancelled, unless a replicate had
// already failed, whose exception is thrown instead.  Returns normally only
// once every replicate has been folded.
template <class RunOne, class Fold>
void run_replicates(int n_reps, int n_threads, const CancelCallback& should_continue,
                    const RunOne& run_one, const Fold& fold) {
  if (n_threads == 0)
    n_threads = static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
  n_threads = std::min(n_threads, n_reps);

  // Replicates are claimed in index order from `next_rep`.  A replicate that
  // finishes ahead of its turn waits in `pending` until every lower-numbered
  // one has been folded in.  Everything below `mu` is guarded by it.
  std::atomic<int> next_rep{0};
  std::mutex mu;
  std::vector<std::unique_ptr<Result>> pending(static_cast<std::size_t>(n_reps));
  int next_fold = 0;

  // Replicates numbered at or above `stop_from` stop at their next poll.  A
  // failure at replicate i lowers it to i, so everything below i still runs
  // and the lowest failing replicate is found deterministically.  A user
  // cancel lowers it to 0.
  std::atomic<int> stop_from{n_reps};
  std::atomic<bool> user_cancelled{false};
  std::mutex cancel_mu;
  int first_error = n_reps;
  std::exception_ptr first_error_ptr;

  auto fold_ready = [&]() {
    while (next_fold < n_reps && pending[next_fold]) {
      std::unique_ptr<Result> const r = std::move(pending[next_fold]);
      ++next_fold;
      fold(std::move(*r));
    }
  };

  // Records replicate i's exception and stops every replicate above it.
  auto record_failure = [&](int i, const std::exception_ptr& e) {
    std::scoped_lock const lock(mu);
    if (i < first_error) {
      first_error = i;
      first_error_ptr = e;
    }
    int cur = stop_from.load();
    while (i < cur && !stop_from.compare_exchange_weak(cur, i)) {
    }
  };

  auto worker = [&]() {
    for (;;) {
      int const i = next_rep.fetch_add(1);
      if (i >= n_reps || i >= stop_from.load())
        return;
      // Everything that can throw stays inside the try: an exception
      // escaping a std::thread is std::terminate.
      try {
        CancelCallback const poll = [&, i]() {
          if (i >= stop_from.load())
            return false;
          if (!should_continue)
            return true;
          std::scoped_lock const lock(cancel_mu);
          if (user_cancelled.load())
            return false;
          if (should_continue())
            return true;
          user_cancelled.store(true);
          stop_from.store(0);
          return false;
        };
        auto r = std::make_unique<Result>(run_one(i, poll));
        std::scoped_lock const lock(mu);
        pending[i] = std::move(r);
        fold_ready();
      } catch (const Cancelled&) {
        // A replicate at or above `stop_from` was stopped here, by a user
        // cancel or a lower-numbered failure, and the ensemble is bound to
        // throw for that reason, so it can be dropped.  Any other Cancelled
        // (should_continue may throw one itself, as run() honors) is this
        // replicate's own failure: dropping it would leave the replicate
        // out of statistics that are then returned as complete.
        if (i < stop_from.load())
          record_failure(i, std::current_exception());
      } catch (...) {
        record_failure(i, std::current_exception());
      }
    }
  };

  if (n_threads == 1) {
    worker();
  } else {
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(n_threads));
    try {
      for (int t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    } catch (...) {
      // Could not start a thread: stop the ones already running.
      stop_from.store(0);
      for (auto& th : threads)
        th.join();
      throw;
    }
    for (auto& th : threads)
      th.join();
  }

  if (first_error_ptr)
    std::rethrow_exception(first_error_ptr);
  if (user_cancelled.load())
    throw Cancelled();
  // Statistics over fewer replicates than n_reps must never be returned.
  if (next_fold != n_reps)
    throw std::logic_error("run_ensemble: only " + std::to_string(next_fold) + " of " +
                           std::to_string(n_reps) + " replicates completed");
}

} // namespace rulemonkey::detail
