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
  explicit RuleMonkeySimulator(const std::string& xml_path, Method method = Method::NfExact);
  ~RuleMonkeySimulator();

  RuleMonkeySimulator(const RuleMonkeySimulator&) = delete;
  RuleMonkeySimulator& operator=(const RuleMonkeySimulator&) = delete;

  // Move is canonical pImpl: only the unique_ptr<Impl> is rebound, no
  // engine state is copied.  Defined out-of-line in simulator.cpp because
  // Impl is incomplete in this header.
  RuleMonkeySimulator(RuleMonkeySimulator&&) noexcept;
  RuleMonkeySimulator& operator=(RuleMonkeySimulator&&) noexcept;

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

  // Returns the current active-session observable values in the same order
  // as `observable_names()`.  Non-const because the engine recomputes lazily
  // — calling this between SSA events forces a fresh evaluation of every
  // observable against current pool state.
  std::vector<double> get_observable_values();

  // Returns the numeric parameter value that the simulator would currently use
  // for the named declared parameter, accounting for any active overrides.
  // Derived parameter expressions (e.g., `B = 2*A` declared in the BNGL)
  // re-resolve when their inputs are overridden, so this reflects the
  // post-cascade value, not the parsed-at-load-time value.
  double get_parameter(const std::string& name) const;

  // Saves the active session state to a file. The session must be active.
  // The file embeds a 64-bit fingerprint of the model schema (molecule
  // types, components, allowed states) which `load_state` verifies on
  // read; parameter values, rate constants, and seed species do not
  // participate, so callers may legitimately mutate those between save
  // and load.
  void save_state(const std::string& path) const;

  // Creates a new session by loading state from a file. Replaces any
  // existing session.  Throws std::runtime_error if the file's schema
  // fingerprint does not match the model XML this simulator was
  // constructed from (the pool serialization is keyed by molecule-type
  // and component indices, so a mismatched XML would silently produce
  // corrupt trajectories).  The schema must match exactly; non-schema
  // fields (parameter values, rate constants, seed concentrations) may
  // differ between save and load.
  void load_state(const std::string& path);

  // Reports whether a live runtime/session state is currently active.
  bool has_session() const;

  // Destroys any live runtime/session state while preserving parsed-model
  // metadata and stored configuration. Safe to call repeatedly.
  void destroy_session();

  // Sets an instance-local override for a declared parameter. The override is
  // applied to subsequent `run()` calls and future `initialize()` calls, and
  // is reflected immediately in `get_parameter()` (including any derived
  // parameters that reference `name`).
  // Throws std::runtime_error if `name` is not a parameter declared in the
  // loaded XML, or if a session is currently active; call `destroy_session()`
  // first to mutate overrides, then re-`initialize()`.
  void set_param(const std::string& name, double value);

  // Clears all instance-local parameter overrides.
  // Throws std::runtime_error if a session is currently active.
  void clear_param_overrides();

  // Sets an instance-local global molecule limit for subsequent runs and
  // future `initialize()` calls.
  // Throws std::runtime_error if a session is currently active.
  void set_molecule_limit(int limit);

  // When true, bimolecular rules only fire between molecules in different
  // complexes (equivalent to NFsim's -bscb flag).  Default: true (strict
  // BNGL semantics).
  // Throws std::runtime_error if a session is currently active.
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

  // Returns the list of unsupported BNGL features detected in the XML at
  // load time.  Each entry has a `Severity` of `Warn` (best-effort run, may
  // still produce useful output) or `Error` (RM cannot honor BNGL semantics
  // for this construct — the rm_driver CLI refuses to run such models
  // unless --ignore-unsupported is passed; embedders should inspect the
  // severities themselves before deciding to call run()).
  const std::vector<UnsupportedFeature>& unsupported_features() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace rulemonkey
