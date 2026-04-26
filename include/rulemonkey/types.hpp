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

struct Result {
  std::vector<double> time;
  std::vector<std::string> observable_names;
  std::vector<std::vector<double>> observable_data;
  int64_t event_count = -1; // SSA events fired (-1 = not reported)

  void set_observable_names(const std::vector<std::string>& names) {
    observable_names = names;
    observable_data.clear();
    observable_data.reserve(observable_names.size());
  }

  void record_time_point(const double t, const std::vector<double>& values) {
    if (!observable_names.empty() && values.size() != observable_names.size()) {
      throw std::runtime_error("observable value count does not match observable name count");
    }
    if (observable_names.empty() && !values.empty()) {
      throw std::runtime_error("cannot record observable values without observable names");
    }
    if (observable_data.empty()) {
      observable_data.resize(values.size());
    }
    time.push_back(t);
    for (std::size_t i = 0; i < values.size(); ++i) {
      observable_data[i].push_back(values[i]);
    }
  }

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
