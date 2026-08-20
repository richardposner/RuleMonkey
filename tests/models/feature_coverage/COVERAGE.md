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
| Reactant counts in a rate law (`reactant_N()`) | Function with an empty Expression | ft_reactant_count_rate | PASS (whole-rule rate functions; a local-function rate law, or an N the rule has no reactant pattern for, is refused at load) |
| Local functions (tags) | Function with Arguments | ft_local_functions | PASS |
| Local function on a bimolecular rule (DOR1) | Function with Arguments, 2 ReactantPatterns | ft_local_fcn_bimol | PASS |
| symmetry_factor on a local-function rate law | @symmetry_factor + Function with Arguments | ft_local_fcn_bimol_sym | PASS (ODE verdict; NFsim 2.9.3 drops the factor) |
| Michaelis-Menten rate law | RateLaw/@type="MM" | ft_mm_ratelaw | PASS |
| symmetry_factor on an MM rate law | @symmetry_factor + RateLaw/@type="MM" | ft_mm_ratelaw_sym | PASS (ODE verdict; NFsim 2.9.3 drops the factor) |
| Bare (global) observable inside a local function | Function with Arguments, untagged Observable reference | ft_local_fcn_global_obs | PASS |
| One observable at local scope in one function and global in another; bare observable moving mid-run | Function with Arguments, same Observable tagged and bare | ft_local_fcn_mixed_scope | PASS (ODE verdict; NFsim 2.9.3 never refreshes the bare observable) |
| Bare observable inside a local function that never moves, at scale | Function with Arguments, untagged Observable reference, 22000 instances | ss_local_fcn_const_global | PASS (cost guard: the per-event rescan must be skipped) |
| FunctionProduct (DOR2) | RateLaw/@type="FunctionProduct" | nf_function_product | PASS |
| Conditional if() | Expression with if() | ft_conditional_rate | PASS |
| Energy patterns | compiled to local functions | ft_energy_patterns | PASS |
| eBNGL Arrhenius (binding) | RateLaw/@type="Arrhenius" + energy patterns | ft_energy_arrhenius | PASS |
| eBNGL Arrhenius (context) | cooperative energy-pattern context conditions | ft_energy_arrhenius_coop | PASS |
| Molecules observable | Observable/@type="Molecules" | all models | PASS |
| Species observable | Observable/@type="Species" | ft_species_vs_molecules | PASS |
| MatchOnce | Observable/@MatchOnce | ft_match_once | PASS |
| Fixed/clamped species ($) | (constant concentration) | ft_clamped_species | PASS |
| Block-same-complex binding | -bscb flag | combo_symmetric_rings | PASS |
| Multi-mol unimolecular | multi-molecule reactant pattern | combo_multimol_unimol | PASS |
| N-ary rule with a complex reactant | 3 ReactantPatterns, one multi-molecule | ft_nary_complex_reactant | PASS |
| Multi-mol reactant with dead-end seed embeddings | multi-molecule ReactantPattern, seed with repeated components | ft_multimol_deadend_seed | PASS |
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

## Refused at load (Tier-0)

| Feature | XML Element | Model | Notes |
|---------|-------------|-------|-------|
| TotalRate modifier | RateLaw/@totalrate="1" | ft_total_rate | Refused because the construct has no agreed meaning — see below. `--ignore-unsupported` runs RM's reading, which is what the model in this suite exercises. |

### Why TotalRate is refused

BioNetGen does not implement TotalRate for network simulations. Its own source
says so, in `RxnRule.pm`:

    # TODO: implement TotalRate feature for Network simulations
    # (currently implemented only for XML network-free output)
    #   --Justin  2dec2010

Confirmed: `generate_network` writes the rate law into the `.net` as an ordinary
rate constant, so the ODE integrates plain mass action and the observable
crashes to zero where NFsim holds it flat. There is no BNG2 result to check a
TotalRate model against, which also means no ODE oracle for this suite.

That leaves NFsim as the only implementation, and it disagrees with RM. NFsim
expands a rule whose reactant pattern has interchangeable components into one
independent reaction class per permutation (`_R1_sym1`, `_R1_sym2`, ...). For an
ordinary rate law that is correct, since the matches partition across the classes
and sum back. Under TotalRate each class returns the whole total rate, so the
rule runs at

    rate x #{permutations whose reactant lists are all non-empty}

Measured on `C(s) + D(t)` at rate `kf`: 1.00x with one free site, 2.02x with two,
2.98x with three, and 1.00x again when every C has the same slot pre-bound, since
only one permutation is populated then. That factor counts NFsim's internal
reaction classes rather than anything in the model: it is capped by the
permutation count however many molecules exist, and it steps down as classes
empty. BNG2's network expansion writes the statistical factor per species and
live (`3*_rateLaw1`, `2*_rateLaw1`, `_rateLaw1`), which implies a third number
again.

So all three disagree, and RM refuses rather than pick one silently.
`--ignore-unsupported` runs the reading BioNetGen documents in `RateLaw.pm`
("If true, this ratelaw specifies the Total reaction rate") — the propensity is
the rate law's value. That agrees with NFsim on every rule whose reactant
patterns have no interchangeable components, which is what `ft_total_rate`
covers.

TotalRate cannot combine with local functions or with MM: BNG2 refuses both
("TotalRate keyword is not compatible with local functions" / "with MM type
RateLaw"), so those shapes never reach RM.

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
| Sat(), Hill() rate laws | generate_network only (Arrhenius energy rules are now supported for 2-reactant binding — see the eBNGL Arrhenius rows above) |
| Pattern quantifiers (==, >=, <) | Not in NFsim XML format |
| Table functions (tfun/TFUN) | BNG2.pl side; not in standard XML |
| time() in functions | BNG2.pl can't generate XML with it |
| continue=>1 | NFsim doesn't support; RM has --save-state/--load-state |

## Recommended RM Improvements

1. **Emit warnings for unsupported XML elements**: RM should warn on stderr when it
   encounters ListOfCompartments, or other elements it doesn't handle.
   Silent incorrect results are worse than loud failures.
