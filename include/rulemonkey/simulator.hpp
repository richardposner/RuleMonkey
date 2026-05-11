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
  //
  // `should_continue`, if non-empty, is polled periodically from the SSA loop
  // (see CancelCallback in types.hpp).  A `false` return throws
  // `rulemonkey::Cancelled` out of this call; the simulator instance is left
  // in a usable state and a subsequent run() may be issued.
  Result run(const TimeSpec& ts, std::uint64_t seed = 42,
             const CancelCallback& should_continue = {});

  // Creates or resets a live session from the parsed model plus the current
  // stored parameter overrides and molecule limit.
  void initialize(std::uint64_t seed = 42);

  // Advances the active session to absolute logical time `time` and
  // discards the sampled trajectory.  Internally this still records
  // observable values at the segment endpoints (current_time and
  // `time`) — the in-process engine has no zero-sample fast path —
  // but the caller sees no `Result`.  Suitable for "equilibrate then
  // perturb then sample" flows where the equilibration trajectory
  // is uninteresting.
  //
  // `should_continue`, if non-empty, enables cooperative cancellation; a
  // `false` return throws `rulemonkey::Cancelled` mid-advance.  The session's
  // `current_time()` then reflects the last fully-applied event, so the
  // caller may inspect partial state, resume by calling step_to / simulate
  // again, or call `destroy_session()` to discard the run.
  void step_to(double time, const CancelCallback& should_continue = {});

  // Samples a segment from the active session starting at the current session
  // time and ending at `t_end`.
  //
  // `should_continue` has the same cooperative-cancellation semantics as on
  // `run()`; if it returns false mid-segment, `Cancelled` is thrown and the
  // session is left at the time of the last completed event.  No partial
  // Result is returned in that case.
  Result simulate(double t_start, double t_end, int n_points,
                  const CancelCallback& should_continue = {});

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
  //
  // *Portability caveat*: the RNG state is serialized via the C++ stdlib's
  // `operator<<` for `std::mt19937_64`, which is required by the standard
  // to round-trip within the same implementation but is NOT specified to
  // be byte-identical across different stdlib implementations (libc++ vs
  // libstdc++ vs MSVC STL).  Save/load is reliably reproducible only
  // between RuleMonkey binaries built against the same stdlib.  Crossing
  // stdlibs (e.g. saving with libstdc++ and loading with libc++) may
  // throw on read or, worse, succeed and continue from a divergent RNG
  // state.  If you need cross-stdlib state portability, generate the
  // checkpoint and resume on the same toolchain — RuleMonkey does not
  // currently hand-serialize the Mersenne Twister state.
  void save_state(const std::string& path) const;

  // Creates a new session by loading state from a file. Replaces any
  // existing session.  Throws std::runtime_error if the file's schema
  // fingerprint does not match the model XML this simulator was
  // constructed from (the pool serialization is keyed by molecule-type
  // and component indices, so a mismatched XML would silently produce
  // corrupt trajectories).  The schema must match exactly; non-schema
  // fields (parameter values, rate constants, seed concentrations) may
  // differ between save and load.
  //
  // *Important caveat*: the fingerprint covers molecule-type schema only
  // — molecule type names, component names, and allowed states.  It does
  // NOT cover ReactionRule patterns/operators or Observable patterns.
  // Loading a session into a simulator constructed from a different XML
  // that adds, removes, or restructures rules / observables (while keeping
  // the molecule-type schema unchanged) will succeed silently.  The pool
  // state is valid in absolute terms but the trajectory continuation
  // will reflect the *new* simulator's rules / observables, which may
  // not be what the caller intended.  Do not rely on save/load to pin
  // down a frozen rule set; instead, treat the source XML as the
  // canonical contract and only resume sessions across simulators built
  // from the same XML (or a strict superset that intentionally extends
  // the rule list).
  void load_state(const std::string& path);

  // Reports whether a live runtime/session state is currently active.
  bool has_session() const;

  // Returns the active session's current logical time.  Useful when
  // resuming after `load_state` to feed `simulate(t_start, …)` with
  // the matching segment-start time.  Throws if no session is active.
  double current_time() const;

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
