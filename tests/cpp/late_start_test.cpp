// A run from a late t_start places its events on the time elapsed since the
// call started (issue #110).
//
// Every waiting time used to be added to the model's clock in doubles.  From
// a large t_start the doubles near t are too far apart to hold it: at 1e12
// they are 1.2e-4 apart against waiting times of about 2e-5, and at 1e13
// they are 2e-3 apart.  A waiting time was rounded away or up to a whole
// spacing, and the run fired too slowly or not at all, with no error.  The
// tolerance for "this sample is at the start" was 1e-12 * |t_end| as well,
// 10 time units at 1e13, so a sample at t_start + 2 recorded the seed state.
// A rate that reads time() was integrated on panels placed on the same
// clock.
//
// The clock is now an origin, where the call starts, plus elapsed time, and
// the SSA advances the elapsed part from 0.  A model whose rates do not read
// the clock then runs the same draws at any t_start: every check below that
// compares a late start with a start at 0 asks for the same trajectory, seed
// for seed, not the same statistics.  Only a rate that reads time() sees the
// rounding, in the clock value it reads.
//
// The same change lets a session's clock reach t_end when nothing is left to
// fire.  It used to stay at the last event, and the next simulate() from
// t_end was refused for starting after the session's time.

#include "rulemonkey/simulator.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <limits>
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

constexpr int kReps = 20;
constexpr double kN0 = 1e4;
constexpr double kK = 0.5;
constexpr double kSpan = 2.0;

// Late starts whose t_start + 0.5*i labels are all exact doubles, and one
// (1e16, doubles 2 apart) where they are not.
const std::vector<double> kLate = {1e9, 1e12, 1e13, 1e15};
constexpr double kUnlabelled = 1e16;

std::uint64_t seed_of(int rep) { return std::uint64_t{7000} + static_cast<std::uint64_t>(rep); }

std::string fmt(double t) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%g", t);
  return buf;
}

std::string exact(double t) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.17g", t);
  return buf;
}

std::size_t index_of(const std::vector<std::string>& names, const std::string& name) {
  for (std::size_t i = 0; i < names.size(); ++i)
    if (names[i] == name)
      return i;
  throw std::runtime_error("name not found: " + name);
}

using Params = std::vector<std::pair<std::string, double>>;

rulemonkey::Result run(const std::string& xml, const Params& params, const rulemonkey::TimeSpec& ts,
                       std::uint64_t seed) {
  rulemonkey::RuleMonkeySimulator sim(xml);
  for (const auto& [name, value] : params)
    sim.set_param(name, value);
  return sim.run(ts, seed);
}

bool same_trajectory(const rulemonkey::Result& a, const rulemonkey::Result& b) {
  return a.observable_data == b.observable_data && a.event_count == b.event_count;
}

rulemonkey::TimeSpec grid(double t0) { return {t0, t0 + kSpan, 4, {}}; }

rulemonkey::TimeSpec listed(double t0) {
  rulemonkey::TimeSpec ts;
  ts.t_start = t0;
  ts.t_end = t0 + kSpan;
  ts.sample_times = {t0, t0 + 1.0, t0 + kSpan};
  return ts;
}

// A binomial ensemble mean of A0 against N0*p, within four standard errors.
void check_mean(const std::string& what, const std::vector<double>& values, double p) {
  double sum = 0.0;
  for (double const v : values)
    sum += v;
  double const mean = sum / static_cast<double>(values.size());
  double const expect = kN0 * p;
  double const tol = 4.0 * std::sqrt(kN0 * p * (1.0 - p) / static_cast<double>(values.size()));
  std::fprintf(stderr, "  %-44s mean=%9.2f  expect=%9.2f  tol=%.2f\n", what.c_str(), mean, expect,
               tol);
  check(std::fabs(mean - expect) <= tol,
        what + ": mean " + std::to_string(mean) + ", expected " + std::to_string(expect));
}

// The plain path: a run over [t0, t0 + 2] is the run over [0, 2], seed for
// seed, on a grid and on listed sample times, and labels its rows t0 + offset.
void check_plain(const std::string& xml) {
  int differ = 0;
  int bad_labels = 0;
  std::vector<double> at_end;
  for (int rep = 0; rep < kReps; ++rep) {
    auto const base_grid = run(xml, {}, grid(0.0), seed_of(rep));
    auto const base_listed = run(xml, {}, listed(0.0), seed_of(rep));
    for (double const t0 : kLate) {
      auto const g = run(xml, {}, grid(t0), seed_of(rep));
      auto const l = run(xml, {}, listed(t0), seed_of(rep));
      if (!same_trajectory(g, base_grid) || !same_trajectory(l, base_listed)) {
        ++differ;
        std::fprintf(stderr, "  seed %d, t0=%s: A0(t0+2) grid %g listed %g, at 0 %g\n", rep,
                     fmt(t0).c_str(), g.observable_data[0].back(), l.observable_data[0].back(),
                     base_grid.observable_data[0].back());
      }
      for (std::size_t i = 0; i < g.time.size(); ++i)
        if (g.time[i] != t0 + (static_cast<double>(i) * 0.5))
          ++bad_labels;
      if (t0 == 1e13)
        at_end.push_back(g.observable_data[0].back());
    }
  }
  check(differ == 0, "a late start and a start at 0 differ in " + std::to_string(differ) + " of " +
                         std::to_string(kReps * kLate.size()) + " (seed, t0) pairs");
  check(bad_labels == 0, std::to_string(bad_labels) + " grid labels were not t0 + 0.5*i");
  check_mean("plain A0(1e13 + 2)", at_end, std::exp(-kK * kSpan));
}

// From 1e16 the grid's labels cannot tell t0 + 0.5 from t0, but its points
// are placed at 0.5*i past t0, so the states are still the run from 0's.
void check_unlabelled_grid(const std::string& xml) {
  int differ = 0;
  for (int rep = 0; rep < kReps; ++rep)
    if (!same_trajectory(run(xml, {}, grid(kUnlabelled), seed_of(rep)),
                         run(xml, {}, grid(0.0), seed_of(rep))))
      ++differ;
  check(differ == 0, "a grid from 1e16 and one from 0 differ in " + std::to_string(differ) +
                         " of " + std::to_string(kReps) + " seeds");
}

// clock() reads each sample's own label, late or not.
void check_clock_column(const std::string& xml) {
  auto const r = run(xml, {}, grid(1e13), 1);
  std::size_t const ci = index_of(r.function_names, "clock");
  for (std::size_t i = 0; i < r.time.size(); ++i)
    check(r.function_data[ci][i] == r.time[i],
          "clock() at sample " + fmt(r.time[i]) + " read " + fmt(r.function_data[ci][i]));
}

// A session started late and advanced in two calls, or by step_to and a
// call, is the session started at 0, seed for seed; each call starts its
// elapsed clock again where the last one stopped.
void check_sessions(const std::string& xml) {
  auto segments = [&xml](double t0, std::uint64_t seed, bool step_first) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.initialize(seed, t0);
    std::vector<std::vector<std::vector<double>>> out;
    if (step_first)
      sim.step_to(t0 + 1.0);
    else
      out.push_back(sim.simulate(t0, t0 + 1.0, 2).observable_data);
    auto const r = sim.simulate(t0 + 1.0, t0 + kSpan, 2);
    out.push_back(r.observable_data);
    check(sim.current_time() == t0 + kSpan,
          "a session from " + fmt(t0) + " ended at " + fmt(sim.current_time()) + ", not t0 + 2");
    return std::make_pair(out, r.event_count);
  };
  int differ = 0;
  for (int rep = 0; rep < kReps; ++rep)
    for (bool const step_first : {false, true}) {
      auto const base = segments(0.0, seed_of(rep), step_first);
      for (double const t0 : kLate)
        if (segments(t0, seed_of(rep), step_first) != base)
          ++differ;
    }
  check(differ == 0, "a late session and one from 0 differ in " + std::to_string(differ) + " of " +
                         std::to_string(2 * kLate.size() * kReps) + " runs");
}

// parameter_scan's fresh points and its carry-over chain from a late
// per_point.t_start end where the scan from 0 does.
void check_scans(const std::string& xml) {
  int differ = 0;
  for (int rep = 0; rep < kReps; ++rep)
    for (bool const reset : {true, false}) {
      auto scan_from = [&](double t0) {
        rulemonkey::RuleMonkeySimulator sim(xml);
        rulemonkey::ScanSpec spec;
        spec.parameter = "k";
        spec.values = {kK, 2.0 * kK};
        spec.per_point = {t0, t0 + kSpan, 2};
        spec.reset_conc = reset;
        return sim.parameter_scan(spec, seed_of(rep)).observable_endpoints;
      };
      if (scan_from(1e13) != scan_from(0.0))
        ++differ;
    }
  check(differ == 0, "a scan from 1e13 and one from 0 differ in " + std::to_string(differ) +
                         " of " + std::to_string(2 * kReps) + " scans");
}

// The clocked path.  With a = b = 0 the rate is k exactly, and the late run
// is the run from 0, seed for seed.  Otherwise the rate reads time() - T0,
// which from 1e13 moves in steps of 2e-3, and the counts still match their
// closed forms: a ramp; a step at T0 + 1 that the panels cross with a step
// a few ulps of the elapsed clock long; and a pulse on (T0 + 1, T0 + 1.1)
// with the rate 0 elsewhere.  The pulse is the only hazard, so one search
// from T0 has to find it, and it does because no panel is longer than a
// thousandth of the call's span.  A thousandth of |t_end| from 1e13 is 1e10,
// one panel over the call, and no Lobatto node falls in the pulse.
void check_clocked(const std::string& clock_xml) {
  int differ = 0;
  for (int rep = 0; rep < kReps; ++rep) {
    auto const base = run(clock_xml, {{"T0", 0.0}}, grid(0.0), seed_of(rep));
    for (double const t0 : kLate)
      if (!same_trajectory(run(clock_xml, {{"T0", t0}}, grid(t0), seed_of(rep)), base))
        ++differ;
  }
  check(differ == 0, "a flat clocked rate from late and from 0 differ in " +
                         std::to_string(differ) + " of " + std::to_string(kReps * kLate.size()) +
                         " (seed, t0) pairs");

  constexpr double kA = 0.25;
  double const p = std::exp(-((kK * kSpan) + (0.5 * kA * kSpan * kSpan)));
  for (double const t0 : {0.0, 1e13}) {
    std::vector<double> at_end;
    std::vector<double> at_mid;
    for (int rep = 0; rep < kReps; ++rep) {
      auto const r = run(clock_xml, {{"a", kA}, {"T0", t0}}, grid(t0), seed_of(rep));
      at_mid.push_back(r.observable_data[0][2]);
      at_end.push_back(r.observable_data[0].back());
    }
    check_mean("ramp A0(" + fmt(t0) + " + 1)", at_mid, std::exp(-(kK + (0.5 * kA))));
    check_mean("ramp A0(" + fmt(t0) + " + 2)", at_end, p);
  }

  constexpr double kB = 0.3;
  for (double const t0 : {0.0, 1e13}) {
    std::vector<double> at_end;
    at_end.reserve(kReps);
    for (int rep = 0; rep < kReps; ++rep)
      at_end.push_back(run(clock_xml, {{"b", kB}, {"T0", t0}}, grid(t0), seed_of(rep))
                           .observable_data[0]
                           .back());
    check_mean("step A0(" + fmt(t0) + " + 2)", at_end, std::exp(-((kK * kSpan) + kB)));
  }

  // A pulse of hazard b*w = 3 on one molecule: it is gone by T0 + 2 with
  // probability 1 - e^-3 = 0.95.
  constexpr int kPulseReps = 40;
  constexpr double kPulseB = 30.0;
  constexpr double kPulseW = 0.1;
  double const p_gone = 1.0 - std::exp(-kPulseB * kPulseW);
  for (double const t0 : {0.0, 1e13}) {
    int gone = 0;
    for (int rep = 0; rep < kPulseReps; ++rep) {
      auto const r =
          run(clock_xml, {{"k", 0.0}, {"b", kPulseB}, {"w", kPulseW}, {"T0", t0}, {"N0", 1.0}},
              grid(t0), seed_of(rep));
      if (r.observable_data[0].back() == 0.0)
        ++gone;
    }
    double const frac = static_cast<double>(gone) / kPulseReps;
    double const tol = 4.0 * std::sqrt(p_gone * (1.0 - p_gone) / kPulseReps);
    std::fprintf(stderr, "  %-44s gone=%5.3f  expect=%9.3f  tol=%.3f\n",
                 ("pulse from " + fmt(t0)).c_str(), frac, p_gone, tol);
    check(std::fabs(frac - p_gone) <= tol, "pulse from " + fmt(t0) + ": the molecule was gone in " +
                                               std::to_string(gone) + " of " +
                                               std::to_string(kPulseReps) + " seeds");
  }
}

// A call leaves the clock at its t_end itself, not at origin + (t_end -
// origin), which for 0.7 and 2.9 is 2.9000000000000004: with events left to
// fire, on either path, and with none.
void check_ends_at_t_end(const std::string& xml, const std::string& clock_xml) {
  for (const auto& [path, n0] :
       {std::make_pair(xml, kN0), std::make_pair(clock_xml, kN0), std::make_pair(xml, 0.0)}) {
    rulemonkey::RuleMonkeySimulator sim(path);
    sim.set_param("N0", n0);
    sim.initialize(5, 0.7);
    sim.simulate(0.7, 2.9, 2);
    check(sim.current_time() == 2.9, "simulate(0.7, 2.9) with N0 = " + fmt(n0) +
                                         " left the clock at " + exact(sim.current_time()));
  }
}

// Each call starts its elapsed clock again at the session's time, so a call
// late in a session resolves its events as one at the start does.  A
// session parked at 1e13, which an empty pool reaches without drawing, then
// given molecules runs as one given them at 0; an origin fixed where the
// session began would place these events 1e13 past it.
void check_rebase_per_call(const std::string& xml) {
  int differ = 0;
  for (int rep = 0; rep < kReps; ++rep) {
    auto from = [&](double t0) {
      rulemonkey::RuleMonkeySimulator sim(xml);
      sim.set_param("N0", 0.0);
      sim.initialize(seed_of(rep), 0.0);
      if (t0 > 0.0)
        sim.step_to(t0);
      sim.add_molecules("A", static_cast<int>(kN0));
      auto const r = sim.simulate(t0, t0 + kSpan, 4);
      return std::make_pair(r.observable_data, r.event_count);
    };
    if (from(1e13) != from(0.0))
      ++differ;
  }
  check(differ == 0, "a session parked at 1e13 and one at 0 differ in " + std::to_string(differ) +
                         " of " + std::to_string(kReps) + " seeds");
}

// A TFUN whose counter is the clock reads the model's clock, not the call's
// elapsed one: in the second call of a session, at each sample's own time
// in its column, and in the rate law, which runs at f = 100*t.
void check_tfun_counter(const std::string& tfun_xml) {
  int bad_column = 0;
  std::vector<double> at_end;
  at_end.reserve(kReps);
  for (int rep = 0; rep < kReps; ++rep) {
    rulemonkey::RuleMonkeySimulator sim(tfun_xml);
    sim.initialize(seed_of(rep), 0.0);
    sim.simulate(0.0, 1.0, 2);
    auto const r = sim.simulate(1.0, kSpan, 2);
    std::size_t const fi = index_of(r.function_names, "f");
    for (std::size_t i = 0; i < r.time.size(); ++i)
      if (r.function_data[fi][i] != 100.0 * r.time[i])
        ++bad_column;
    at_end.push_back(r.observable_data[index_of(r.observable_names, "Btot")].back());
  }
  check(bad_column == 0, std::to_string(bad_column) + " samples of the TFUN column in a " +
                             "session's second call did not read 100*t");
  // B(2) is Poisson with mean Lambda(2) = 50*2^2 = 200.
  double sum = 0.0;
  for (double const v : at_end)
    sum += v;
  double const mean = sum / static_cast<double>(at_end.size());
  double const tol = 4.0 * std::sqrt(200.0 / static_cast<double>(at_end.size()));
  std::fprintf(stderr, "  %-44s mean=%9.2f  expect=%9.2f  tol=%.2f\n", "TFUN session B(2)", mean,
               200.0, tol);
  check(std::fabs(mean - 200.0) <= tol, "TFUN session B(2): mean " + std::to_string(mean));
}

// t_end must be finite: an infinite one has no grid to sample and, once
// nothing can fire, no time to leave the clock at.  step_to(inf) on a pool
// that decays out would leave current_time() at inf, a state load_state
// cannot read back; run({0, inf, 2}) labelled its rows NaN, inf, inf.  A
// refused call leaves the session where it was.  simulate(TimeSpec) refuses
// a t_end before the session's time, which bounds the run whatever the
// samples are.  A NaN segment start, which simulate()'s tolerance check let
// through, is refused too, and so is a NaN anywhere in sample_times.  A zero-length window whose
// t_start drifted an ulp past t_end is still accepted, as that tolerance says it is.
void check_refusals(const std::string& xml) {
  auto throws = [](const std::function<void()>& f) {
    try {
      f();
    } catch (const std::runtime_error&) {
      return true;
    }
    return false;
  };
  double const inf = std::numeric_limits<double>::infinity();
  rulemonkey::RuleMonkeySimulator one_shot(xml);
  check(throws([&] { one_shot.run({0.0, inf, 2}, 1); }), "run({0, inf, 2}) should throw");
  check(throws([&] { one_shot.run({0.0, std::nan(""), 2}, 1); }), "run({0, NaN, 2}) should throw");
  check(throws([&] {
          one_shot.run_ensemble({0.0, inf, 2}, {2, 1, false}, 1);
        }),
        "run_ensemble({0, inf, 2}) should throw");
  check(throws([&] {
          rulemonkey::ScanSpec spec;
          spec.parameter = "k";
          spec.values = {kK};
          spec.per_point = {0.0, inf, 2};
          one_shot.parameter_scan(spec, 1);
        }),
        "parameter_scan over [0, inf] should throw");

  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_param("N0", 10.0);
  sim.initialize(3, 0.0);
  check(throws([&] { sim.step_to(inf); }), "step_to(inf) should throw");
  check(throws([&] { sim.simulate(0.0, inf, 2); }), "simulate(0, inf, 2) should throw");
  check(sim.current_time() == 0.0,
        "a refused call moved the session to " + exact(sim.current_time()));
  check(!throws([&] { sim.simulate(0.0, 1.0, 2); }), "the session should go on after a refusal");

  rulemonkey::RuleMonkeySimulator back(xml);
  back.initialize(3, 10.0);
  rulemonkey::TimeSpec before;
  before.t_start = 10.0;
  before.t_end = 5.0;
  before.sample_times = {10.0, 11.0, 12.0};
  check(throws([&] { back.simulate(before); }),
        "simulate(TimeSpec) with t_end before the session's time should throw");

  check(throws([&] { sim.simulate(std::nan(""), 2.0, 2); }), "simulate(NaN, 2, 2) should throw");
  rulemonkey::TimeSpec nan_start;
  nan_start.t_end = 2.0;
  nan_start.sample_times = {std::nan(""), 2.0};
  check(throws([&] { sim.simulate(nan_start); }),
        "simulate(TimeSpec) whose first sample is NaN should throw");

  // A NaN later in sample_times passed the ordering check and stalled the
  // sample cursor, so every sample after it recorded the final state.
  rulemonkey::TimeSpec nan_mid;
  nan_mid.t_end = 4.0;
  nan_mid.sample_times = {0.0, std::nan(""), 1.0, 2.0, 4.0};
  check(throws([&] { one_shot.run(nan_mid, 1); }), "run() with a NaN sample time should throw");
  rulemonkey::TimeSpec nan_mid_session = nan_mid;
  nan_mid_session.sample_times = {1.0, std::nan(""), 2.0, 4.0};
  check(throws([&] { sim.simulate(nan_mid_session); }),
        "simulate(TimeSpec) with a NaN sample time should throw");
  check(sim.current_time() == 1.0,
        "a refused sample list moved the session to " + exact(sim.current_time()));

  rulemonkey::RuleMonkeySimulator drift(xml);
  drift.initialize(3, 0.3);
  double const drifted = 0.1 * 3.0; // 0.30000000000000004
  check(!throws([&] { drift.simulate(drifted, 0.3, 1); }),
        "a zero-length window whose t_start drifted past t_end should run");
}

// A session saved at a late time and loaded continues as the one left in
// process does, draw for draw: both start the next call's elapsed clock at
// the same current_time.
void check_late_save_load(const std::string& xml) {
  // Unique per run, so two builds testing at once do not share the file.
  auto const path =
      std::filesystem::temp_directory_path() /
      ("rm_late_start_test." +
       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".state");
  int differ = 0;
  try {
    for (int rep = 0; rep < kReps; ++rep) {
      rulemonkey::RuleMonkeySimulator in_process(xml);
      in_process.initialize(seed_of(rep), 1e13);
      in_process.simulate(1e13, 1e13 + 1.0, 2);
      in_process.save_state(path.string());
      auto const a = in_process.simulate(1e13 + 1.0, 1e13 + kSpan, 2);

      rulemonkey::RuleMonkeySimulator loaded(xml);
      loaded.load_state(path.string());
      check(loaded.current_time() == 1e13 + 1.0,
            "the loaded session is at " + exact(loaded.current_time()) + ", not 1e13 + 1");
      auto const b = loaded.simulate(1e13 + 1.0, 1e13 + kSpan, 2);
      if (!same_trajectory(a, b))
        ++differ;
    }
  } catch (...) {
    std::filesystem::remove(path);
    throw;
  }
  std::filesystem::remove(path);
  check(differ == 0, "a late session continued in process and from its saved state differ in " +
                         std::to_string(differ) + " of " + std::to_string(kReps) + " seeds");
}

// Once nothing can fire, the session's clock still goes to t_end, and the
// next call starts there.  Two shapes: a pool that decays out partway
// through the call, and one that is empty from the start.
void check_absorbing(const std::string& xml) {
  for (double const n0 : {10.0, 0.0}) {
    rulemonkey::RuleMonkeySimulator sim(xml);
    sim.set_param("N0", n0);
    sim.initialize(3, 0.0);
    std::string const what = "N0 = " + fmt(n0) + ": ";
    try {
      auto const r = sim.simulate(0.0, 100.0, 2);
      check(r.observable_data[0].back() == 0.0, what + "A0 did not decay out by t = 100");
      check(sim.current_time() == 100.0,
            what + "simulate(0, 100) left the clock at " + fmt(sim.current_time()));
      sim.simulate(100.0, 200.0, 2);
      check(sim.current_time() == 200.0,
            what + "simulate(100, 200) left the clock at " + fmt(sim.current_time()));
      sim.step_to(300.0);
      check(sim.current_time() == 300.0,
            what + "step_to(300) left the clock at " + fmt(sim.current_time()));
    } catch (const std::runtime_error& e) {
      check(false, what + "continuing past the absorbing state threw: " + e.what());
    }
  }

  // A carry-over chain's point that decays out hands the next point a clock
  // at its t_end, so point p ends at 100*(p + 1).  It used to hand on the
  // time of the last event, about 3.7, so every later point ran from there
  // and ended at about 103.7.
  rulemonkey::RuleMonkeySimulator sim(xml);
  sim.set_param("N0", 10.0);
  rulemonkey::ScanSpec spec;
  spec.parameter = "k";
  spec.values = {kK, kK, kK};
  spec.per_point = {0.0, 100.0, 1};
  spec.reset_conc = false;
  auto const scan = sim.parameter_scan(spec, 3);
  std::size_t const ci = index_of(scan.function_names, "clock");
  for (std::size_t p = 0; p < spec.values.size(); ++p)
    check(scan.function_endpoints[p][ci] == 100.0 * static_cast<double>(p + 1),
          "carry-over point " + std::to_string(p) +
              " ended at clock() = " + exact(scan.function_endpoints[p][ci]));
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::fprintf(stderr, "Usage: late_start_test <late_start_model.xml> "
                         "<late_start_clock_model.xml> <tfun_time_rate_model.xml>\n");
    return 2;
  }
  const std::string xml = argv[1];
  const std::string clock_xml = argv[2];
  const std::string tfun_xml = argv[3];

  // Each check runs on its own, so one that throws does not hide the rest.
  auto run_check = [](const char* name, const std::function<void()>& f) {
    try {
      f();
    } catch (const std::exception& e) {
      check(false, std::string(name) + " threw: " + e.what());
    }
  };
  run_check("check_plain", [&] { check_plain(xml); });
  run_check("check_unlabelled_grid", [&] { check_unlabelled_grid(xml); });
  run_check("check_clock_column", [&] { check_clock_column(xml); });
  run_check("check_sessions", [&] { check_sessions(xml); });
  run_check("check_scans", [&] { check_scans(xml); });
  run_check("check_clocked", [&] { check_clocked(clock_xml); });
  run_check("check_ends_at_t_end", [&] { check_ends_at_t_end(xml, clock_xml); });
  run_check("check_rebase_per_call", [&] { check_rebase_per_call(xml); });
  run_check("check_tfun_counter", [&] { check_tfun_counter(tfun_xml); });
  run_check("check_refusals", [&] { check_refusals(xml); });
  run_check("check_late_save_load", [&] { check_late_save_load(xml); });
  run_check("check_absorbing", [&] { check_absorbing(xml); });

  if (g_failures != 0) {
    std::fprintf(stderr, "%d failure(s)\n", g_failures);
    return 1;
  }
  std::fprintf(stderr, "PASS\n");
  return 0;
}
