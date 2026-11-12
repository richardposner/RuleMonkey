// The deprecated loose mode, block_same_complex_binding = false (CLI
// -no-bscb), pinned until 4.0 removes it, and its deprecation warning
// (issue #139).
//
// The fixture (`no_bscb_ring_dimer_model.bngl/.xml`) is one
// A(x,y!1).B(x,y!1) dimer under `Ring: A(x) + B(x) <-> A(x!1).B(x!1)`.
// Its free x sites can only join inside the complex.  BNG2's network has
// no reaction that does that: its reactions join dimers into chains
// (generate_network never closes; max_iter=>1 gives only dimer + dimer ->
// tetramer), and with one dimer none of them can fire.  So the x ring
// never closes.  With the flag off RM tests neither half of the `+`: the
// forward rule closes the ring at kf, and the reverse opens it at kr
// although its products stay in one complex.  The ring is then a two-state
// Markov chain, closed a fraction kf/(kf+kr) of the time, as in NFsim run
// with neither -cb nor -bscb (0.46 and 0.48 at kf = kr = 1, two seeds,
// issue #139).
//
// Pinned here:
//   * with the flag on, by default or set explicitly, `Closed` is 0 at
//     every sample, and no deprecation warning is printed;
//   * with it off, at kf = kr = 1 and at kf = 3, kr = 1, the fraction of
//     samples with the ring closed matches kf/(kf+kr), which pins the ratio
//     of the rates, and the number of closings seen between samples matches
//     its expectation, which pins their sum;
//   * the first set_block_same_complex_binding(false) in the process prints
//     the deprecation warning, on one line naming -no-bscb, 4.0 and #139,
//     and later calls on that simulator or another print it no more.
//
// The stderr checks look only for the warning, so a debug print such as
// RM_PRINT_TIMING does not upset them.

#include "rulemonkey/simulator.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int g_failures = 0;

// Diagnostics go to stdout: stderr is redirected to the capture file.
void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stdout, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

std::size_t count_occurrences(const std::string& haystack, const std::string& needle) {
  std::size_t count = 0;
  std::size_t pos = 0;
  while ((pos = haystack.find(needle, pos)) != std::string::npos) {
    ++count;
    pos += needle.size();
  }
  return count;
}

std::string read_file(const std::string& path) {
  std::fflush(stderr);
  std::ifstream const f(path);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

constexpr double kTEnd = 2000.0;
constexpr int kPoints = 20000; // samples dt = 0.1 apart

// `Closed` at every sample after t = 0, which is the seed state (open).
std::vector<double> closed_series(rulemonkey::RuleMonkeySimulator& sim, std::uint64_t seed) {
  rulemonkey::TimeSpec ts;
  ts.t_start = 0.0;
  ts.t_end = kTEnd;
  ts.n_points = kPoints;
  auto const r = sim.run(ts, seed);
  for (std::size_t i = 0; i < r.observable_names.size(); ++i)
    if (r.observable_names[i] == "Closed")
      return {r.observable_data[i].begin() + 1, r.observable_data[i].end()};
  throw std::runtime_error("observable not found: Closed");
}

// With the flag off the ring is a two-state chain, closing at kf and
// opening at kr.  Over n samples dt apart:
//
//   * it is closed for a fraction p = kf/(kf+kr) of them.  Samples are
//     correlated by rho = exp(-(kf+kr) dt), so the fraction's standard
//     error is sqrt(p(1-p)/n * (1+rho)/(1-rho));
//   * an open sample is followed by a closed one with probability
//     p(1-rho), so n p(1-p)(1-rho) closings are seen.  That depends on
//     kf + kr through rho: doubling or halving both rates moves it by over
//     20 standard deviations.  Its spread is below that of the true number
//     of closings, a renewal count with variance T (1/kf^2 + 1/kr^2) /
//     (1/kf + 1/kr)^3 (19.7 against 22.4 at kf = kr = 1, 300 seeds of an
//     exact simulation of the chain), which is used as the scale.
//
// Each bound is 5 of its standard deviations.
void check_loose(rulemonkey::RuleMonkeySimulator& sim, double kf, double kr, std::uint64_t seed) {
  std::string const label =
      "-no-bscb, kf = " + std::to_string(kf) + ", kr = " + std::to_string(kr) + ": ";
  auto const closed = closed_series(sim, seed);
  auto const n = static_cast<double>(closed.size());
  double sum = 0.0;
  int closings = 0;
  int openings = 0;
  double prev = 0.0;
  for (double const c : closed) {
    check(c == 0.0 || c == 1.0, label + "Closed must be 0 or 1, got " + std::to_string(c));
    sum += c;
    closings += (prev == 0.0 && c == 1.0) ? 1 : 0;
    openings += (prev == 1.0 && c == 0.0) ? 1 : 0;
    prev = c;
  }
  double const frac = sum / n;
  double const p = kf / (kf + kr);
  double const rho = std::exp(-(kf + kr) * (kTEnd / kPoints));
  double const se = std::sqrt(p * (1.0 - p) / n * (1.0 + rho) / (1.0 - rho));
  double const mean_cycle = (1.0 / kf) + (1.0 / kr);
  double const expected_closings = n * p * (1.0 - p) * (1.0 - rho);
  double const closings_sd =
      std::sqrt(kTEnd * ((1.0 / (kf * kf)) + (1.0 / (kr * kr))) / std::pow(mean_cycle, 3));
  std::fprintf(stdout,
               "%sclosed fraction %.4f (expected %.4f +- %.4f), %d closings (expected %.0f +- "
               "%.0f), %d openings\n",
               label.c_str(), frac, p, 5.0 * se, closings, expected_closings, 5.0 * closings_sd,
               openings);
  check(std::fabs(frac - p) < 5.0 * se, label + "closed fraction " + std::to_string(frac) +
                                            ", expected " + std::to_string(p) + " +- " +
                                            std::to_string(5.0 * se));
  check(std::fabs(closings - expected_closings) < 5.0 * closings_sd,
        label + std::to_string(closings) + " closings, expected " +
            std::to_string(expected_closings) + " +- " + std::to_string(5.0 * closings_sd));
  // Each opening follows a closing.
  check(openings == closings || openings == closings - 1,
        label + std::to_string(closings) + " closings but " + std::to_string(openings) +
            " openings");
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stdout, "Usage: no_bscb_deprecation_test <no_bscb_ring_dimer_model.xml> "
                         "<stderr_capture_path>\n");
    return 2;
  }
  std::string const xml = argv[1];
  std::string const err_path = argv[2];

  // Not restored: the test ends right after its checks.
  if (std::freopen(err_path.c_str(), "w", stderr) == nullptr) {
    std::fprintf(stdout, "FAIL: could not redirect stderr to %s\n", err_path.c_str());
    return 1;
  }

  std::string const warning = "set_block_same_complex_binding(false) (-no-bscb) is deprecated";

  try {
    // The default, then the same set explicitly: BNGL's `+`, as BNG2.
    {
      rulemonkey::RuleMonkeySimulator sim(xml);
      for (bool const set_explicitly : {false, true}) {
        if (set_explicitly)
          sim.set_block_same_complex_binding(true);
        std::string const how = set_explicitly ? "set to true" : "by default";
        for (std::uint64_t const seed : {std::uint64_t{1}, std::uint64_t{2}}) {
          auto const closed = closed_series(sim, seed);
          double max_closed = 0.0;
          for (double const c : closed)
            max_closed = std::fmax(max_closed, c);
          check(max_closed == 0.0, "the x ring closed with block_same_complex_binding " + how +
                                       ", seed " + std::to_string(seed));
        }
      }
    }
    std::string const strict_err = read_file(err_path);
    check(strict_err.find(warning) == std::string::npos,
          "block_same_complex_binding = true printed the deprecation warning");

    // The loose mode.  The first `false` in the process warns, once.
    {
      rulemonkey::RuleMonkeySimulator sim(xml);
      sim.set_block_same_complex_binding(false);
      std::string const first = read_file(err_path);
      std::fprintf(stdout, "--- captured stderr ---\n%s--- end ---\n", first.c_str());
      check(count_occurrences(first, warning) == 1,
            "the first set_block_same_complex_binding(false) must print the deprecation "
            "warning once");
      // The warning's own line, from its start to the newline.
      std::size_t const at = first.find(warning);
      std::size_t const eol = first.find('\n', at == std::string::npos ? 0 : at);
      std::string const line = at == std::string::npos || eol == std::string::npos
                                   ? std::string()
                                   : first.substr(at, eol - at);
      check(line.find("RuleMonkey 4.0") != std::string::npos,
            "the deprecation warning must say, on its one line, that the setting goes in 4.0");
      check(line.find("#139") != std::string::npos,
            "the deprecation warning must point, on its one line, at issue #139");

      check_loose(sim, 1.0, 1.0, 1);

      sim.set_block_same_complex_binding(false);
      rulemonkey::RuleMonkeySimulator other(xml);
      other.set_block_same_complex_binding(false);
      other.set_param("kf", 3.0);
      check_loose(other, 3.0, 1.0, 1);

      std::string const all = read_file(err_path);
      check(count_occurrences(all, warning) == 1,
            "later set_block_same_complex_binding(false) calls must not print the warning again");
    }
  } catch (const std::exception& e) {
    std::fprintf(stdout, "FAIL: threw: %s\n", e.what());
    return 1;
  }

  if (g_failures == 0) {
    std::fprintf(stdout, "PASS: -no-bscb ring dimer pinned; deprecation warning printed once\n");
    return 0;
  }
  std::fprintf(stdout, "FAILED: %d assertion(s)\n", g_failures);
  return 1;
}
