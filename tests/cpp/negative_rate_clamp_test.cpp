// Regression test for the negative-propensity clamp diagnostic.
//
// Pins down two contracts of set_rule_propensity in engine.cpp:
//
//   1. When a function rate law evaluates to a negative number, the
//      affected rule's propensity is clamped to zero — the simulation
//      keeps running rather than throwing (the strict-validation
//      alternative NFsim picks).
//   2. The first clamp per rule emits one stderr WARN line so the
//      author sees the issue.  Further clamps on the same rule are
//      silent so an oscillator around the zero crossing doesn't spam.
//
// The fixture model (`negative_rate_clamp_model.bngl/.xml`) drives
// `Obs_A` from 10 down toward zero via a first-order decay, and a
// production rule with rate function `prod_rate() = Obs_A - setpoint`
// (setpoint=5).  Once Obs_A drops below 5, prod_rate goes negative on
// every subsequent propensity update for that rule.

#include "rulemonkey/simulator.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

namespace {

int g_failures = 0;

void check(bool ok, const std::string& msg) {
  if (!ok) {
    std::fprintf(stdout, "FAIL: %s\n", msg.c_str());
    ++g_failures;
  }
}

// Count how many times `needle` appears in `haystack`.
std::size_t count_occurrences(const std::string& haystack, const std::string& needle) {
  std::size_t count = 0;
  std::size_t pos = 0;
  while ((pos = haystack.find(needle, pos)) != std::string::npos) {
    ++count;
    pos += needle.size();
  }
  return count;
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::fprintf(stdout, "Usage: negative_rate_clamp_test <xml_path> <stderr_capture_path>\n");
    return 2;
  }
  const std::string xml_path = argv[1];
  const std::string err_path = argv[2];

  // Redirect stderr to a temp file for the duration of the run.  We do
  // not restore it: the test exits immediately after the assertions
  // below, and our own diagnostic messages go to stdout so they're
  // still visible in ctest output.
  if (std::freopen(err_path.c_str(), "w", stderr) == nullptr) {
    std::fprintf(stdout, "FAIL: could not redirect stderr to %s\n", err_path.c_str());
    return 1;
  }

  try {
    rulemonkey::RuleMonkeySimulator sim(xml_path);
    rulemonkey::TimeSpec ts;
    ts.t_start = 0.0;
    ts.t_end = 20.0;
    ts.n_points = 100;
    auto const result = sim.run(ts, /*seed=*/1);
    (void)result;
  } catch (const std::exception& e) {
    std::fflush(stderr);
    std::fprintf(stdout, "FAIL: sim.run threw: %s\n", e.what());
    return 1;
  }
  std::fflush(stderr);

  // Read back what the engine wrote to stderr.
  std::ifstream const ef(err_path);
  std::stringstream ss;
  ss << ef.rdbuf();
  const std::string err_blob = ss.str();

  // Surface the captured stderr to stdout so a failing run is debuggable
  // from the ctest log.
  std::fprintf(stdout, "--- captured stderr ---\n%s--- end ---\n", err_blob.c_str());

  // Contract 1: the clamp produced exactly one WARN line.  Two would
  // mean the per-rule latch failed (every clamp prints); zero would
  // mean either the trajectory never drove Obs_A below the setpoint
  // (fixture broken) or the warning was lost in a refactor.
  std::size_t const n_warn = count_occurrences(err_blob, "WARN: rule '");
  std::size_t const n_clamped = count_occurrences(err_blob, "clamped to 0");
  check(n_warn == 1, "exactly one WARN line expected, got " + std::to_string(n_warn));
  check(n_clamped == 1,
        "exactly one 'clamped to 0' marker expected, got " + std::to_string(n_clamped));

  // Contract 2: the message names the offending rate function so the
  // author can find it.  The fixture's function is `prod_rate`.
  check(err_blob.find("prod_rate") != std::string::npos,
        "WARN line should name the rate function 'prod_rate'");

  // Contract 3: sim.run() didn't throw (covered by reaching here).
  // No assertion needed beyond having gotten this far.

  if (g_failures == 0) {
    std::fprintf(stdout, "PASS: negative-rate clamp warning fires once and names the function\n");
    return 0;
  }
  std::fprintf(stdout, "FAILED: %d assertion(s)\n", g_failures);
  return 1;
}
