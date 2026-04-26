#pragma once

// Interface-only contract snapshot for the RuleMonkey <-> BNGsim boundary.
// This header is intentionally free of engine code.

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "rulemonkey/types.hpp"

namespace rulemonkey {

enum class Method {
  NfExact,
};

// Public embedding surface for in-process callers.
//
// Construction parses and validates the XML immediately. Each simulator owns
// its parsed model, parameter overrides, and molecule-limit configuration; no
// prior run state is reused across `run()` calls.
class RuleMonkeySimulator {
 public:
  // Parses and validates `xml_path` immediately.
  //
  // Throws std::runtime_error if the path is empty or missing, the method is
  // unsupported, or the XML cannot be imported by the active runtime.
  explicit RuleMonkeySimulator(const std::string& xml_path,
                               Method method = Method::NfExact);
  ~RuleMonkeySimulator();

  RuleMonkeySimulator(const RuleMonkeySimulator&) = delete;
  RuleMonkeySimulator& operator=(const RuleMonkeySimulator&) = delete;
  RuleMonkeySimulator(RuleMonkeySimulator&&) = delete;
  RuleMonkeySimulator& operator=(RuleMonkeySimulator&&) = delete;

  // Executes a fresh simulation from the parsed model plus the current
  // instance-local parameter overrides and molecule limit.
  //
  // Each call seeds a new RNG from `seed`; repeated calls with the same XML,
  // method, TimeSpec, overrides, molecule limit, and seed are reproducible
  // within the same runtime build. Prior successful or failed runs do not
  // mutate reusable simulator state.
  Result run(const TimeSpec& times, std::uint64_t seed = 42);

  // Creates or resets a live session from the parsed model plus the current
  // stored parameter overrides and molecule limit.
  void initialize(std::uint64_t seed = 42);

  // Advances the active session to absolute logical time `time` without
  // returning sampled output.
  void step_to(double time);

  // Samples a segment from the active session starting at the current session
  // time and ending at `t_end`.
  Result simulate(double t_start, double t_end, int n_points);

  // Adds `count` default, unbound molecules of the named imported
  // `MoleculeType` to the active session.
  void add_molecules(const std::string& molecule_type_name, int count);

  // Returns the current active-session molecule count for the named imported
  // `MoleculeType`.
  int get_molecule_count(const std::string& molecule_type_name) const;

  // Returns the current active-session observable values in the same order as
  // `observable_names()`.
  std::vector<double> get_observable_values() const;

  // Returns the numeric parameter value that the simulator would currently use
  // for the named declared parameter.
  double get_parameter(const std::string& name) const;

  // Saves the active session state to a file. The session must be active.
  void save_state(const std::string& path) const;

  // Creates a new session by loading state from a file. Replaces any existing
  // session. The model XML must match the one used to save the state.
  void load_state(const std::string& path);

  // Reports whether a live runtime/session state is currently active.
  bool has_session() const;

  // Destroys any live runtime/session state while preserving parsed-model
  // metadata and stored configuration. Safe to call repeatedly.
  void destroy_session();

  // Sets an instance-local override for a declared parameter. The override is
  // applied to subsequent `run()` calls and future `initialize()` calls.
  void set_param(const std::string& name, double value);

  // Clears all instance-local parameter overrides.
  void clear_param_overrides();

  // Sets an instance-local global molecule limit for subsequent runs and
  // future `initialize()` calls.
  void set_molecule_limit(int limit);

  // When true, bimolecular rules only fire between molecules in different
  // complexes (equivalent to NFsim's -bscb flag).  Default: false.
  void set_block_same_complex_binding(bool value);

  // Returns observable names in XML declaration order captured at
  // construction. The returned vector is a copy.
  std::vector<std::string> observable_names() const;

  // Returns parameter names in XML declaration order captured at
  // construction. The returned vector is a copy.
  std::vector<std::string> parameter_names() const;

  // Returns the original constructor XML path string.
  const std::string& xml_path() const;

  // Returns the validated runtime method for this simulator instance.
  Method method() const;

  // Returns warnings about unsupported BNGL features detected in the XML.
  const std::vector<UnsupportedFeature>& unsupported_warnings() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace rulemonkey
