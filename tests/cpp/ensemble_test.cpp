// Replicate-ensemble regression (issue #113): pins down the run_ensemble()
// contracts an embedder computing error bars depends on.
//
//   * replicate i is bit-identical to run(ts, base_seed + i), parameter
//     overrides included — any one replicate can be reproduced alone;
//   * the mean / SD / event_count are bit-identical for every n_threads,
//     because statistics are folded in replicate order;
//   * the mean and sample SD (n - 1) agree with a two-pass computation over
//     the kept trajectories, for observables and global functions;
//   * n_reps == 1 gives the plain run() as the mean and an all-zero SD;
//   * invalid specs throw;
//   * an active session neither blocks run_ensemble nor changes its result,
//     and is left exactly where it was, as with run();
//   * a should_continue returning false throws Cancelled, the callback is
//     never entered by two workers at once, and the simulator stays usable;
//   * an exception escaping a replicate is rethrown as itself — not as
//     Cancelled, not std::terminate — with one thread or several;
//   * a Cancelled thrown by should_continue itself reaches the caller and is
//     never dropped as if run_ensemble had stopped that replicate.
//
// argv[1] = A_plus_A.xml (kp; exact integer counts, so equality is exact)
// argv[2] = ft_nested_functions.xml (three global functions)

#include "rulemonkey/simulator.hpp"

#include <atomic>
#include <cmath>
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

const rulemonkey::TimeSpec kTs{0.0, 5.0, 20};

rulemonkey::EnsembleSpec spec(int n_reps, int n_threads, bool keep = false) {
  rulemonkey::EnsembleSpec s;
  s.n_reps = n_reps;
  s.n_threads = n_threads;
  s.keep_trajectories = keep;
  return s;
}

// Replicate i == run(ts, base_seed + i), with an override in force so the
// test also covers overrides reaching every worker's engine.
void test_replicates_match_run(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_param("kp", 0.002);
  auto ens = sim.run_ensemble(kTs, spec(6, 3, /*keep=*/true), /*base_seed=*/100);

  check(ens.trajectories.size() == 6, "keep_trajectories must return one Result per replicate");
  check(ens.n_reps == 6 && ens.base_seed == 100, "n_reps / base_seed must be echoed");
  int64_t events = 0;
  for (std::size_t i = 0; i < ens.trajectories.size(); ++i) {
    auto r = sim.run(kTs, 100 + i);
    const auto& t = ens.trajectories[i];
    check(t.time == r.time && t.observable_data == r.observable_data &&
              t.function_data == r.function_data && t.event_count == r.event_count,
          "replicate " + std::to_string(i) + " must be bit-identical to run(ts, base_seed + i)");
    events += r.event_count;
  }
  check(ens.event_count == events, "event_count must be the sum over replicates");
  check(ens.time == ens.trajectories[0].time, "time must be the replicates' time grid");
}

void test_thread_count_invariance(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto base = sim.run_ensemble(kTs, spec(12, 1));
  for (int const n_threads : {2, 5, 0}) {
    auto e = sim.run_ensemble(kTs, spec(12, n_threads));
    check(e.observable_mean == base.observable_mean && e.observable_sd == base.observable_sd &&
              e.function_mean == base.function_mean && e.function_sd == base.function_sd &&
              e.event_count == base.event_count,
          "statistics must be bit-identical for n_threads=" + std::to_string(n_threads) +
              " and n_threads=1");
  }
}

bool close(double a, double b) { return std::fabs(a - b) <= 1e-9 * std::max(1.0, std::fabs(b)); }

// Two-pass mean / sample SD over `col(r)[c][t]` for every kept replicate.
template <class Col>
void check_stats(const rulemonkey::EnsembleResult& e, Col col,
                 const std::vector<std::vector<double>>& mean,
                 const std::vector<std::vector<double>>& sd, const std::string& what) {
  const auto n = static_cast<double>(e.trajectories.size());
  bool ok = mean.size() == col(e.trajectories[0]).size() && sd.size() == mean.size();
  for (std::size_t c = 0; ok && c < mean.size(); ++c) {
    ok = mean[c].size() == e.n_times() && sd[c].size() == e.n_times();
    for (std::size_t t = 0; ok && t < e.n_times(); ++t) {
      double m = 0;
      for (const auto& r : e.trajectories)
        m += col(r)[c][t];
      m /= n;
      double ss = 0;
      for (const auto& r : e.trajectories)
        ss += (col(r)[c][t] - m) * (col(r)[c][t] - m);
      ok = close(mean[c][t], m) && close(sd[c][t], std::sqrt(ss / (n - 1)));
    }
  }
  check(ok, what + " mean / sample SD must match a two-pass computation over the replicates");
}

void test_statistics(const std::string& xml, const std::string& label) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto e = sim.run_ensemble(kTs, spec(9, 4, /*keep=*/true));
  check(e.observable_names == e.trajectories[0].observable_names &&
            e.function_names == e.trajectories[0].function_names,
        label + ": column names must be the replicates' names");
  check_stats(
      e, [](const rulemonkey::Result& r) -> const auto& { return r.observable_data; },
      e.observable_mean, e.observable_sd, label + " observable");
  check_stats(
      e, [](const rulemonkey::Result& r) -> const auto& { return r.function_data; },
      e.function_mean, e.function_sd, label + " function");
}

void test_single_replicate(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  auto e = sim.run_ensemble(kTs, spec(1, 0), /*base_seed=*/7);
  auto r = sim.run(kTs, 7);
  check(e.observable_mean == r.observable_data, "n_reps=1 mean must be the plain run");
  bool zero = true;
  for (const auto& col : e.observable_sd)
    for (double const v : col)
      zero = zero && v == 0.0;
  check(zero, "n_reps=1 SD must be all zero");
  check(e.trajectories.empty(), "trajectories must be empty without keep_trajectories");
}

template <class F> bool throws_runtime_error(F f) {
  try {
    f();
  } catch (const std::runtime_error&) {
    return true;
  }
  return false;
}

void test_validation(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  check(throws_runtime_error([&] { sim.run_ensemble(kTs, spec(0, 1)); }), "n_reps=0 must throw");
  check(throws_runtime_error([&] { sim.run_ensemble(kTs, spec(2, -1)); }),
        "n_threads<0 must throw");
  check(throws_runtime_error([&] {
          sim.run_ensemble({5.0, 1.0, 4}, spec(2, 1));
        }),
        "t_end < t_start must throw");
}

void test_during_session(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator fresh(xml);
  auto const want = fresh.run_ensemble(kTs, spec(4, 2));

  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.initialize(1);
  sim.step_to(2.0);
  double const t_session = sim.current_time();
  auto const obs_session = sim.get_observable_values();
  auto const got = sim.run_ensemble(kTs, spec(4, 2));
  check(got.observable_mean == want.observable_mean && got.observable_sd == want.observable_sd &&
            got.event_count == want.event_count,
        "an active session must not change run_ensemble's result");
  check(sim.has_session() && sim.current_time() == t_session &&
            sim.get_observable_values() == obs_session,
        "run_ensemble must leave an active session where it was");
  sim.step_to(4.0);
  check(sim.current_time() == 4.0, "the session must still advance after run_ensemble");
}

void test_cancel(const std::string& xml) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  // Plain ints on purpose: the callback runs under run_ensemble's lock, so
  // two workers never enter it at once (a TSan build would flag them if
  // they did); `inside` catches an overlap in any build.
  int polls = 0;
  int inside = 0;
  bool overlapped = false;
  bool threw_cancelled = false;
  try {
    sim.run_ensemble({0.0, 1000.0, 10}, spec(16, 4), 1, [&] {
      overlapped = overlapped || inside != 0;
      ++inside;
      bool const go = ++polls < 6;
      --inside;
      return go;
    });
  } catch (const rulemonkey::Cancelled&) {
    threw_cancelled = true;
  }
  check(threw_cancelled, "a false should_continue must throw Cancelled");
  check(!overlapped, "should_continue must never be entered by two workers at once");
  auto e = sim.run_ensemble(kTs, spec(4, 2));
  check(e.n_times() == 21, "the simulator must stay usable after a cancelled ensemble");
}

void test_replicate_exception(const std::string& xml) {
  for (int const n_threads : {1, 4}) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    std::atomic<int> polls{0};
    bool threw_logic = false;
    try {
      sim.run_ensemble({0.0, 1000.0, 10}, spec(8, n_threads), 1, [&]() -> bool {
        if (polls.fetch_add(1) >= 3)
          throw std::logic_error("boom");
        return true;
      });
    } catch (const std::logic_error&) {
      threw_logic = true;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "  unexpected exception: %s\n", e.what());
    }
    check(threw_logic, "a replicate's own exception must be rethrown unchanged (n_threads=" +
                           std::to_string(n_threads) + ")");
    auto e = sim.run_ensemble(kTs, spec(2, n_threads));
    check(e.n_times() == 21, "the simulator must stay usable after a failed ensemble");
  }
}

// A Cancelled thrown by should_continue itself (run() passes it through) is
// that replicate's failure, not a replicate run_ensemble stopped.  Taking it
// for one dropped the replicate and returned statistics over the rest as
// if all n_reps had run.
void test_callback_throws_cancelled(const std::string& xml) {
  for (int const n_threads : {1, 4}) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    std::atomic<int> polls{0};
    bool threw_cancelled = false;
    try {
      sim.run_ensemble({0.0, 1000.0, 10}, spec(8, n_threads), 1, [&]() -> bool {
        if (polls.fetch_add(1) == 2)
          throw rulemonkey::Cancelled();
        return true;
      });
    } catch (const rulemonkey::Cancelled&) {
      threw_cancelled = true;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "  unexpected exception: %s\n", e.what());
    }
    check(threw_cancelled,
          "a Cancelled thrown by should_continue must reach the caller (n_threads=" +
              std::to_string(n_threads) + ")");
    auto e = sim.run_ensemble(kTs, spec(2, n_threads));
    check(e.n_times() == 21, "the simulator must stay usable after a callback throws Cancelled");
  }
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: ensemble_test <A_plus_A.xml> <ft_nested_functions.xml>\n");
    return 1;
  }
  const std::string aa = argv[1];
  const std::string fn = argv[2];

  test_replicates_match_run(aa);
  test_thread_count_invariance(aa);
  test_thread_count_invariance(fn);
  test_statistics(aa, "A_plus_A");
  test_statistics(fn, "ft_nested_functions");
  test_single_replicate(aa);
  test_validation(aa);
  test_during_session(aa);
  test_cancel(aa);
  test_replicate_exception(aa);
  test_callback_throws_cancelled(aa);

  if (g_failures > 0) {
    std::fprintf(stderr, "ensemble_test: %d failure(s)\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "ensemble_test: all checks passed\n");
  return 0;
}
