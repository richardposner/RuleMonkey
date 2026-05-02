#pragma once

// Interface-only contract snapshot for the RuleMonkey <-> BNGsim boundary.
// This header is intentionally free of engine code.

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace rulemonkey {

struct TimeSpec {
  double t_start = 0.0;
  double t_end = 0.0;
  int n_points = 0;
};

// Plain data carrier for one simulation run's output.
//
// `observable_data` is column-major: `observable_data[obs_idx][t_idx]`
// gives the value of the obs_idx-th observable at the t_idx-th sample
// time.  Both inner and outer vectors have the same length on every
// successful run (n_observables and n_times respectively).
//
// Construction and population is the engine's responsibility — host
// code reads but does not write to these fields.
struct Result {
  std::vector<double> time;
  std::vector<std::string> observable_names;
  std::vector<std::vector<double>> observable_data;
  int64_t event_count = 0; // SSA events fired during this run.

  std::size_t n_times() const noexcept { return time.size(); }
  std::size_t n_observables() const noexcept { return observable_names.size(); }
};

enum class Severity { Warn, Error };

struct UnsupportedFeature {
  Severity severity;
  std::string element; // XML element or attribute name
  std::string feature; // human-readable BNGL feature description
};

} // namespace rulemonkey
