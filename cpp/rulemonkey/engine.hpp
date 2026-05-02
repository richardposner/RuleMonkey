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
  Result run(const TimeSpec& ts);

  double current_time() const;
  std::vector<double> get_observable_values();
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
