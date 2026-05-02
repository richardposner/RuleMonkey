# Feature Coverage Matrix

## RM-Supported Features (tested and passing)

| Feature | XML Element | Model(s) | Status |
|---------|-------------|----------|--------|
| Component states | Component/@state | all models | PASS |
| Bond operations (add/delete) | ListOfOperations | all binding models | PASS |
| Bond wildcards (!+, !?) | Component/@numberOfBonds | ft_bond_wildcards | PASS |
| State wildcards (~?) | (omitted state) | ft_state_wildcards | PASS |
| Symmetric components | multiple same-name components | ft_multi_site_binding, ft_blbr | PASS |
| Molecule creation (synthesis) | AddMolecule operation | ft_synthesis_degradation | PASS |
| Molecule deletion | DeleteMolecule operation | ft_delete_molecules | PASS |
| DeleteMolecules keyword | @DeleteMolecules="1" | ft_delete_molecules | PASS |
| Reversible rules | forward + reverse RateLaw | all <-> models | PASS |
| Elementary rate laws | RateLaw/@type="Ele" | most models | PASS |
| Functional rate laws | RateLaw/@type="Function" | ft_functional_rate | PASS |
| Local functions (tags) | Function with Arguments | ft_local_functions | PASS |
| Conditional if() | Expression with if() | ft_conditional_rate | PASS |
| Energy patterns | compiled to local functions | ft_energy_patterns | PASS |
| Molecules observable | Observable/@type="Molecules" | all models | PASS |
| Species observable | Observable/@type="Species" | ft_species_vs_molecules | PASS |
| MatchOnce | Observable/@MatchOnce | ft_match_once | PASS |
| Fixed/clamped species ($) | (constant concentration) | ft_clamped_species | PASS |
| Block-same-complex binding | -bscb flag | combo_symmetric_rings | PASS |
| Multi-mol unimolecular | multi-molecule reactant pattern | combo_multimol_unimol | PASS |
| Push-pull enzyme kinetics | Michaelis-Menten cycle pattern | ft_push_pull | PASS |
| Ring closure | intramolecular bonds | ft_ring_closure | PASS |
| Signaling cascade | multi-step recruitment | ft_signaling_cascade | PASS |
| State continuation (save/load) | --save-state/--load-state | ft_continue | PASS |
| exclude_reactants | ListOfExcludeReactants | ft_exclude_reactants, combo_exclude_with_complex | PASS (ODE verdict; NFsim ignores) |
| include_reactants | ListOfIncludeReactants | ft_include_reactants | PASS (ODE verdict; NFsim ignores) |
| exclude_products | ListOfExcludeProducts | ft_exclude_products | PASS (ODE verdict; NFsim ignores) |
| include_products | ListOfIncludeProducts | (tested via reverse rule in exclude_products) | PASS (ODE verdict; NFsim ignores) |

## RM-Supported but requires ODE verdict

NFsim ignores exclude/include reactant/product constraints in network-free mode.
For these features, RM is tested against BNG2 ODE (via `generate_network`) as the verdict reference.
The models are designed with strong constraint effects so RM-vs-NFsim diverges dramatically (confirming NFsim ignores the feature)
while RM-vs-ODE matches (confirming RM correctly implements the feature).

## RM-Supported but FAILING

| Feature | XML Element | Model | Status | Notes |
|---------|-------------|-------|--------|-------|
| TotalRate modifier | RateLaw/@totalrate="1" | ft_total_rate | **FAIL** | RM parses it but behavior differs from NFsim |

## RM Does NOT Support (silently ignored)

**These XML elements are present in some models but RM does not parse them.
Tests that appear to pass for these features are FALSE POSITIVES.**

| Feature | XML Element | Model(s) | Notes |
|---------|-------------|----------|-------|
| Compartments | ListOfCompartments | (none) | RM is well-mixed; no compartment support |
| MoveConnected | MoveConnected attribute | (none) | Requires compartments |
| Population types | ListOfPopulationTypes | (none) | Hybrid particle-population only |
| Rule priority | @priority | (none) | Not parsed |

## Not Applicable to Network-Free Simulation

| Feature | Notes |
|---------|-------|
| Sat(), Hill(), Arrhenius() rate laws | generate_network only |
| Pattern quantifiers (==, >=, <) | Not in NFsim XML format |
| Table functions (tfun/TFUN) | BNG2.pl side; not in standard XML |
| time() in functions | BNG2.pl can't generate XML with it |
| continue=>1 | NFsim doesn't support; RM has --save-state/--load-state |

## Recommended RM Improvements

1. **Emit warnings for unsupported XML elements**: RM should warn on stderr when it
   encounters ListOfCompartments, or other elements it doesn't handle.
   Silent incorrect results are worse than loud failures.

2. **Investigate TotalRate**: RM parses totalrate="1" but the behavior diverges from
   NFsim. The propensity calculation may need correction.
