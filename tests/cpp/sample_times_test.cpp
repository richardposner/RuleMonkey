// Explicit output-time (sample_times) regression for TimeSpec (issue #16).
//
// Pins down three contracts of the `TimeSpec::sample_times` override that
// `Engine::run_ssa` honors:
//
//   1. Non-invasive + bit-identical.  Sampling at an explicit *subset* of a
//      uniform grid (built with the exact same double arithmetic the engine
//      uses) yields output bit-identical to the corresponding rows of the
//      full uniform run, with the same seed — output times never draw from
//      the RNG or perturb reaction selection, so the realization (and
//      `Result::event_count`) is unchanged.  This is the headline claim from
//      the issue: explicit times honored exactly, single-pass, no drift.
//
//   2. Arbitrary non-uniform instants.  A hand-picked, non-uniform list
//      produces exactly one output row per requested time, at exactly those
//      times (the recorder pushes the requested sample value verbatim).
//
//   3. Ordering contract.  An out-of-order list throws std::runtime_error
//      rather than silently mis-sampling.
//
// The fixture is A+A<->AA (reference/nfsim/xml/A_plus_A.xml): its observables
// are molecule counts — integers exactly representable as doubles — so the
// bit-identical comparison in (1) is a true equality test, not an artifact of
// floating-point coincidence.

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
    std::fprintf(stdout, "Usage: sample_times_test <xml_path>\n");
    return 2;
  }
  const std::string xml_path = argv[1];
  constexpr std::uint64_t kSeed = 20260619;

  try {
    // ---- Contract 1: explicit subset == uniform grid, bit for bit --------
    //
    // Build the uniform grid here with the SAME expression the engine uses
    // (`t_start + i*dt`, i an int promoted to double) so the explicit list
    // carries doubles bit-identical to the engine's internal grid.  If they
    // differed even in the last bit, an event landing between the two
    // representations could be recorded on different sides and the
    // "bit-identical" comparison would be testing nothing.
    const double t_start = 0.0;
    const double t_end = 10.0;
    const int n_points = 100;
    const double dt = (t_end - t_start) / n_points;

    std::vector<double> grid;
    grid.reserve(n_points + 1);
    for (int i = 0; i <= n_points; ++i)
      grid.push_back(t_start + (i * dt));

    rulemonkey::RuleMonkeySimulator sim(xml_path);

    rulemonkey::TimeSpec uniform_ts;
    uniform_ts.t_start = t_start;
    uniform_ts.t_end = t_end;
    uniform_ts.n_points = n_points;
    const auto uniform = sim.run(uniform_ts, kSeed);

    check(uniform.n_times() == static_cast<std::size_t>(n_points) + 1,
          "uniform run should produce n_points+1 rows, got " + std::to_string(uniform.n_times()));

    // Explicit list = every 10th grid point (indices 0,10,...,100).
    constexpr int kStride = 10;
    std::vector<double> subset;
    for (int i = 0; i <= n_points; i += kStride)
      subset.push_back(grid[i]);

    rulemonkey::TimeSpec explicit_ts;
    explicit_ts.t_start = t_start;
    explicit_ts.t_end = t_end;
    explicit_ts.n_points = 999; // must be IGNORED when sample_times is set
    explicit_ts.sample_times = subset;
    const auto sampled = sim.run(explicit_ts, kSeed);

    check(sampled.n_times() == subset.size(),
          "explicit run row count should equal sample_times.size() (" +
              std::to_string(subset.size()) + "), got " + std::to_string(sampled.n_times()));
    check(sampled.n_observables() == uniform.n_observables(),
          "observable count must match between the two runs");
    // event_count is a property of the realization, which the output schedule
    // must not touch.
    check(sampled.event_count == uniform.event_count,
          "event_count must be identical regardless of sampling schedule (uniform=" +
              std::to_string(uniform.event_count) +
              ", sampled=" + std::to_string(sampled.event_count) + ")");

    if (sampled.n_times() == subset.size() && sampled.n_observables() == uniform.n_observables()) {
      bool all_equal = true;
      for (std::size_t k = 0; k < subset.size(); ++k) {
        const std::size_t u = k * kStride; // matching uniform-grid row
        if (sampled.time[k] != uniform.time[u]) {
          all_equal = false;
          std::fprintf(stdout, "  time mismatch at k=%zu: sampled=%.17g uniform=%.17g\n", k,
                       sampled.time[k], uniform.time[u]);
        }
        for (std::size_t o = 0; o < sampled.n_observables(); ++o) {
          if (sampled.observable_data[o][k] != uniform.observable_data[o][u]) {
            all_equal = false;
            std::fprintf(stdout, "  obs[%zu] mismatch at k=%zu: sampled=%.17g uniform=%.17g\n", o,
                         k, sampled.observable_data[o][k], uniform.observable_data[o][u]);
          }
        }
      }
      check(all_equal, "explicit-subset run must be bit-identical to the uniform run at "
                       "the shared instants");
    }

    // ---- Contract 2: arbitrary non-uniform instants ----------------------
    //
    // Distinct, non-uniform, not-on-any-grid times.  The recorder pushes the
    // requested sample value verbatim, so Result::time must echo the list,
    // and the row count must equal the list length even for a time beyond the
    // last reaction (flushed at the final state).
    const std::vector<double> nonuniform = {0.5, 1.7, 3.3, 7.0, 9.99};
    rulemonkey::TimeSpec nu_ts;
    nu_ts.t_start = 0.0;
    nu_ts.t_end = 10.0;
    nu_ts.sample_times = nonuniform;
    const auto nu = sim.run(nu_ts, kSeed);

    check(nu.n_times() == nonuniform.size(),
          "non-uniform run should produce one row per requested time (" +
              std::to_string(nonuniform.size()) + "), got " + std::to_string(nu.n_times()));
    if (nu.n_times() == nonuniform.size()) {
      bool times_match = true;
      for (std::size_t k = 0; k < nonuniform.size(); ++k)
        if (nu.time[k] != nonuniform[k])
          times_match = false;
      check(times_match, "Result::time must echo the requested sample_times verbatim");
    }

    // ---- Contract 3: ordering is validated -------------------------------
    bool threw = false;
    std::string what;
    try {
      rulemonkey::TimeSpec bad_ts;
      bad_ts.t_start = 0.0;
      bad_ts.t_end = 10.0;
      bad_ts.sample_times = {1.0, 0.5, 2.0}; // descending step at index 1
      (void)sim.run(bad_ts, kSeed);
    } catch (const std::exception& e) {
      threw = true;
      what = e.what();
    }
    check(threw, "an out-of-order sample_times must throw");
    check(what.find("ascending") != std::string::npos,
          "the throw message should mention the ascending-order contract; got: " + what);

  } catch (const std::exception& e) {
    std::fprintf(stdout, "FAIL: unexpected exception: %s\n", e.what());
    return 1;
  }

  if (g_failures == 0) {
    std::fprintf(stdout, "PASS: explicit sample_times honored exactly, non-invasive, "
                         "and order-validated\n");
    return 0;
  }
  std::fprintf(stdout, "FAILED: %d assertion(s)\n", g_failures);
  return 1;
}
