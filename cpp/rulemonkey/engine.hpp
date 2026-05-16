#pragma once

#include "model.hpp"
#include "rulemonkey/types.hpp"

#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace rulemonkey {

class Engine {
public:
  Engine(const Model& model, uint64_t seed, int molecule_limit = -1);
  ~Engine();

  // Engine owns a unique_ptr<Impl> with PIMPL semantics — copying it
  // would alias the per-engine RNG and pool state, which is never what
  // a caller wants.  Move is also disabled because callers always take
  // an Engine by reference (the live session is held inside the
  // surrounding Simulator).
  Engine(const Engine&) = delete;
  Engine& operator=(const Engine&) = delete;
  Engine(Engine&&) = delete;
  Engine& operator=(Engine&&) = delete;

  // Initialize agent pool from seed species; compute initial propensities.
  void initialize();

  // Run SSA from current time, recording observable values at sample points.
  // `should_continue`, if non-empty, is polled inside the SSA loop; returning
  // false throws `rulemonkey::Cancelled` out of run() at a safe between-event
  // point (no partial event state is left exposed).
  Result run(const TimeSpec& ts, const CancelCallback& should_continue = {});

  double current_time() const;
  std::vector<double> get_observable_values();

  // Global-function values for the current pool state, in the same order
  // as `function_names()`.  Non-const for the same reason as
  // get_observable_values(): the values are re-evaluated on demand.
  std::vector<double> get_function_values();

  // Names of the model's global (non-local) functions, in XML
  // declaration order.  Parallel to a Result's `function_names`.
  std::vector<std::string> function_names() const;

  int get_molecule_count(const std::string& type_name) const;
  void add_molecules(const std::string& type_name, int count);

  // Save full simulation state to file; load restores pool and derived state.
  void save_state(const std::string& path) const;
  void load_state(const std::string& path);

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace rulemonkey
