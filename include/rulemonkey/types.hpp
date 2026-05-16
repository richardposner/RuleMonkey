#pragma once

// Interface-only contract snapshot for the RuleMonkey <-> BNGsim boundary.
// This header is intentionally free of engine code.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace rulemonkey {

// Cooperative-cancellation hook for long-running run() / simulate() / step_to()
// calls.  The engine polls this callback periodically from the SSA event loop
// (granularity is roughly every 1024 SSA events, fine enough for wall-clock
// timeouts but coarse enough to keep per-event overhead negligible).  Return
// `true` to let the simulation keep running, `false` to cancel; cancellation
// causes the engine to throw a `Cancelled` exception out of the current call.
// A default-constructed (empty) std::function disables the check entirely —
// callers who don't need cancellation pay no per-event overhead.
using CancelCallback = std::function<bool()>;

// Thrown out of run() / simulate() / step_to() when the supplied
// `CancelCallback` returns false.  Inherits from std::runtime_error so generic
// catch handlers still see it, but a dedicated type lets embedders distinguish
// a cooperative cancellation (a clean stop requested from outside) from a
// genuine engine error.
class Cancelled : public std::runtime_error {
public:
  Cancelled() : std::runtime_error("RuleMonkey simulation cancelled") {}
  explicit Cancelled(const std::string& msg) : std::runtime_error(msg) {}
};

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
// `function_data` mirrors `observable_data` for BNGL global functions —
// the entries of a `begin functions ... end functions` block that have
// no local (per-molecule/per-species) arguments.  These are the derived
// quantities models routinely use as their measured/fitted outputs
// (e.g. `Clusters() = monomer + dimer + ...`).  The engine already
// evaluates them for rate laws; this is an exposure of values it
// already computes.  `function_data[fn_idx][t_idx]` gives the value of
// the fn_idx-th global function at the t_idx-th sample time, sampled at
// the same time points as the observables.  Local functions are
// excluded — they evaluate per-molecule and have no single global
// value — so `function_names` may be shorter than the model's full
// `begin functions` block.
//
// Construction and population is the engine's responsibility — host
// code reads but does not write to these fields.
struct Result {
  std::vector<double> time;
  std::vector<std::string> observable_names;
  std::vector<std::vector<double>> observable_data;
  std::vector<std::string> function_names;
  std::vector<std::vector<double>> function_data;
  int64_t event_count = 0; // SSA events fired during this run.

  std::size_t n_times() const noexcept { return time.size(); }
  std::size_t n_observables() const noexcept { return observable_names.size(); }
  std::size_t n_functions() const noexcept { return function_names.size(); }
};

enum class Severity { Warn, Error };

struct UnsupportedFeature {
  Severity severity;
  std::string element; // XML element or attribute name
  std::string feature; // human-readable BNGL feature description
};

} // namespace rulemonkey
