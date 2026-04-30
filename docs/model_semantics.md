# RuleMonkey BNGL Coverage

A reference for "will my BNGL model run on RuleMonkey?" RM consumes the
XML emitted by `BNG2.pl writeXML`; this document is keyed off that XML.

The full programmatic inventory lives in
`scan_unsupported()` in `cpp/rulemonkey/simulator.cpp` — every refusal and
warning below is emitted from that function and surfaced to embedders via
`RuleMonkeySimulator::unsupported_features()` (and to the CLI via
`rm_driver`'s startup banner).

The runtime severity model is two-level:

- **Error** — RM cannot honor BNGL semantics for the construct. The
  `rm_driver` CLI refuses to run such models with exit code 2 unless
  `--ignore-unsupported` is passed; embedders should inspect each
  feature's `Severity` themselves before deciding to call `run()`.
- **Warn** — best-effort run. The construct is parsed but partially
  honored (or ignored where ignoring is benign). Trajectories may differ
  from a fully-conformant simulator on edge cases.

## Supported

### Molecule and rule structure

- Multi-component molecule types with named components and enumerated
  per-component states (`A(s~U~P, b)`).
- Pattern matching with all four bond constraints: `Free` (`x`),
  `Bound` (`x!+` or `x!?`), `BoundTo` (`x!1` shared label),
  `Wildcard` (component omitted).
- Multi-molecule reactant patterns (`A(b!1).B(a!1)`), including
  multi-bond rings.
- Symmetric components (`P(s,s)`) with the correct combinatorial
  weighting (the engine carries a "same components" flag and a
  symmetry factor).
- Reactant constraints: `ListOfExcludeReactants`, `ListOfIncludeReactants`,
  `ListOfExcludeProducts`, `ListOfIncludeProducts`. Multiple patterns
  inside a single list are OR'd; the engine evaluates all of them on each
  candidate match.

### Reaction operations

- `AddBond` (including ring-closing bonds and bonds to molecules added
  in the same rule).
- `DeleteBond`, with `ensureConnected="1"` honored as a
  product-molecularity check at fire time.
- `StateChange`, including state-increment / state-decrement (the BNGL
  `s~U~P~PP+` shorthand).
- `AddMolecule` (synthesis from nothing or from existing reactants).
- `DeleteMolecule`, both `DeleteMolecules="1"` (single-molecule) and
  the default whole-species delete.

### Rate laws

- **`Ele`** — elementary mass-action. Rate constants resolve from the
  parameter map at run time (`set_param` overrides take effect, including
  through derived-parameter chains — see "Parameter overrides" below).
- **`Function`** — arbitrary expression evaluated against a live
  variable map. Supports parameters, observables, the special variables
  `time` / `t`, and references to other functions. Local functions
  (per-molecule or per-pattern arguments) are supported in both
  per-molecule and complex-wide scopes.
- **`MM(kcat, Km)`** — Michaelis-Menten via NFsim's QSS formula:
  `sFree = 0.5·((S − Km − E) + √((S − Km − E)² + 4·Km·S))`,
  `a = kcat·sFree·E/(Km + sFree)`. Mirrors `MMRxnClass::update_a` in
  the NFsim source.
- **`TFUN`** — table-function rate laws backed by external `.tfun`
  files. RM searches both the XML directory and one level up to handle
  both author-side and harness-side layouts. Counter sources may be
  `time`, a parameter, an observable, or another function.

### Observables

- `Molecules` (counts every matching molecule across all complexes).
- `Species` (counts complexes that match, treating each complex as one
  unit). Stoichiometric quantifiers `==N`, `>=N`, `<=N` on Species
  observables are supported (RM has a fast incremental path for these
  and substantially outperforms NFsim on count-relation Species
  observables — see `docs/timing_comparison.md`).
- Multi-pattern observables (BNGL `Molecules X A(),B()` with multiple
  patterns under one observable) — the engine sums embeddings across
  all listed patterns.

### Strict BNGL semantics

- `block_same_complex_binding` is **on by default** (matches NFsim's
  `-bscb` flag). Bimolecular rules only fire across distinct complexes;
  intramolecular ring-closure must be written as an explicit
  `A(x).B(y) -> A(x!1).B(y!1)` rule. Disable via
  `set_block_same_complex_binding(false)` or the CLI `-no-bscb` flag for
  parity with NFsim runs that omitted `-bscb`.

### Parameter overrides

- `set_param(name, value)` rejects names not declared in the loaded XML
  (typos throw rather than silently no-op).
- Overrides cascade through derived parameter expressions. If the BNGL
  declares `B = 2*A` (i.e. the XML emits
  `<Parameter id="B" value="2*A"/>`), `set_param("A", x)` recomputes
  `B` to `2*x` for the next run AND for `get_parameter("B")` queries
  in between runs. Overriding `B` directly wins over the cascade — the
  expression for `B` is skipped, and parameters that derive from `B`
  see the override.
- Cascade order is the parameter declaration order in the XML; BNG2
  emits parameters in dependency order, so a forward reference inside
  a derivation that the XML happens to expose without a dependency-
  ordered emit will still resolve as long as the chain is one level
  deep (the parser does a single forward-reference retry pass at load
  time; `apply_overrides` does not iterate to fixed point).
- `get_parameter(name)` reflects the current overrides + cascade
  immediately, without requiring a `run()` or `initialize()` call.

### Initial state and live mutation

- Seed species declared via `ListOfSpecies` with `concentration="N"` or
  `concentration="param_name"` — parameter-backed concentrations
  re-resolve through `set_param` overrides.
- Single-molecule `Fixed` species (BNGL `$` prefix) — the engine
  replenishes them after each event so their count is held at the
  initial value (matching BNG2's ODE semantics).
- `add_molecules(type_name, count)` for live perturbation between
  segments of a stateful session.

### Functional surface

- `set_param`, `clear_param_overrides`, `set_molecule_limit`,
  `set_block_same_complex_binding` — applied at the next `run()` /
  `initialize()` (throw if a session is currently active).
- `save_state(path)` / `load_state(path)` — full pool, RNG, and
  bookkeeping snapshot. The XML used at `load_state` time must match
  the one used at `save_state` time.

## Tier-0 refusals (Error severity)

These are emitted as `Severity::Error` entries on
`unsupported_features()`. The CLI exits 2 by default; embedders should
treat at least one of these as a signal to refuse the model.

| Trigger | Why refused |
|---|---|
| `<ListOfCompartments>` non-empty | RM does not implement compartment volume scaling — bimolecular rate constants would be silently incorrect. |
| Any rule with `RateLaw type="Arrhenius"` | eBNGL energy-pattern rate derivation is not implemented; rate constants would be silently wrong. (A bare `<ListOfEnergyPatterns>` with `Function`-type rate laws that inline the Boltzmann factors is fine — only `Arrhenius` is the trigger.) |
| Any rule with `RateLaw type="Sat"` | Deprecated; rewrite as `MM(kcat, Km)`. |
| Any rule with `RateLaw type="Hill"` | Network-only; use `generate_network()` + ODE/SSA instead of network-free. |
| Any rule with `RateLaw type="FunctionProduct"` | Not implemented; rewrite as a single multi-factor `Function`. |
| Any `<MoleculeType population="1">` | Hybrid particle-population SSA not implemented; would be silently treated as ordinary particles with diverging trajectories. |
| Multi-molecule or bonded `<Species Fixed="1">` | RM v1 only supports single-molecule, unbonded Fixed species. |
| Two or more `<Species Fixed="1">` of the same `MoleculeType` | RM v1 allows at most one Fixed species per molecule type to avoid matching-overlap ambiguity. |

The CLI's `--ignore-unsupported` flag downgrades these to runs-anyway
mode. Each error message includes the specific behavior change that
results — embedders inspecting the feature list can replicate that
decision per-feature.

## Best-effort warnings (Warn severity)

These are emitted as `Severity::Warn` and the run proceeds.

| Trigger | Effect |
|---|---|
| Any rule with `MoveConnected` operation | Requires compartments; emitted as a warning because RM ignores the operation entirely. |
| Any rule with a `priority` attribute | Honored as ordinary rule firing; the priority modifier is ignored. |

## Embedder integration pattern

```cpp
rulemonkey::RuleMonkeySimulator sim("model.xml");

bool any_errors = false;
for (const auto& f : sim.unsupported_features()) {
  if (f.severity == rulemonkey::Severity::Error) {
    log_error("RM cannot honor BNGL feature {} ({}): {}",
              f.element, f.feature, /* your model id */);
    any_errors = true;
  } else {
    log_warn("RM warning: {} — {}", f.element, f.feature);
  }
}
if (any_errors) {
  // Decide per-feature whether to refuse or fall back to a
  // different backend; do NOT call sim.run() if you'd misinterpret
  // an Error-level construct.
  return;
}

auto result = sim.run({0.0, t_end, n_points}, seed);
```

## Open work tracked elsewhere

- Compartment volume scaling — open work, no scheduled implementation.
- Arrhenius / energy-pattern rate derivation — same.
- Hybrid particle-population SSA — same.
- Multi-molecule Fixed species — would require pattern-based
  re-instantiation; not currently implemented (refused at Tier 0).
- Pattern canonical labeling (nauty integration) — flagged in
  `CHANGELOG.md` as the leading candidate for the next major
  performance work; no impact on coverage, only on speed for models
  with many isomorphic complexes.
