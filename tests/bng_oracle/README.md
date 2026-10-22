# BNG2-oracle regression suite

Models whose expected behaviour comes from **BioNetGen itself** — `BNG2.pl`'s
`generate_network()` expansion, integrated by the BioNetGen solver — rather than
from an NFsim ensemble.

## Why this is separate from `tests/reference/`

The rest of the corpus scores RuleMonkey against NFsim ensembles, which works
because RuleMonkey and NFsim agree almost everywhere. These models are the
exception: they pin behaviour where **NFsim is wrong**, so an NFsim reference
would be permanently red for the wrong reason and would pull RuleMonkey back
onto the defect it just fixed. The issues these models come from establish
ground truth by reading BNG2's generated reaction network directly, and this
suite keeps that same oracle.

A model belongs here when its reference has to be BioNetGen's own answer. A
model that merely happens to agree with NFsim belongs in the regular corpus.

## Reference kinds

Recorded per model in `generate_references.py`.

**`ode`** — `generate_network()` + `simulate(method=>"ode")`. The default, and
what the issue tables are read off. Valid when populations are large enough
that the mass-action ODE *is* the mean of the SSA; for a catalytic rule with a
conserved catalyst the substrate decay is pseudo-first-order, so it is exact.
`context_dor_price` qualifies for the same reason once its catalysts have
settled, which they do inside the first 1% of the run.

**`ssa`** — `generate_network()` + an ensemble of `simulate(method=>"ssa")`
reps, with mean and SEM committed. Needed for a model that deliberately runs at
a handful of complexes, where the ODE prices an `X + X` reaction as `[X]^2`
while the SSA prices it as `N(N-1)` — a factor of 2 at two copies. That is the
ordinary deterministic-limit gap and says nothing about the model; it just means
the deterministic solver is answering a different question than a 4-complex
stochastic model asks. Same tool, same expanded network, stochastic integrator.

`context_nary` takes it for the other reason: its rate is a product of two
substrate pools that deplete together and a catalyst-complex population that
fluctuates as dimers form and break, so the ODE of the means is not the mean of
the SSA. A stochastic reference sidesteps the question entirely — it is the same
process, not an approximation to it.

## Models

| model | reference | pins |
|---|---|---|
| `context_symmetry` | ode | An untransformed (context) reactant pattern is counted once per matching **complex**, not once per matching molecule — including the `Bind_sym` guard in the other direction, where the rule *does* transform the pattern and the 2x is correct (GH #33). |
| `context_sampler` | ssa | The same slot is **drawn** the way it is **counted**. Counting per complex while sampling per molecule leaves the same-complex rejection mis-weighted, which is a systematic bias wherever complexes hold unequal numbers of context matches (GH #33). |
| `context_nary` | ssa | The same per-complex count on a rule with **three or more** reactant patterns, which takes the n-ary distinct-tuple machinery rather than the two-slot propensity (GH #42). Pins it on a static homodimer (2x), a static homotrimer ring against a single-molecule pattern (3x), and a catalyst dimer that forms and breaks for the whole run — plus `Guard_trans`, an n-ary rule that *does* transform its multi-subunit slot and must keep its 2x. |
| `context_dor_price` | ode | **Which** molecule of the collapsed complex prices the surviving instance, when the context pattern also carries a per-molecule local function tag: the species' canonically-first match, as BNG2 expands it, and not the lowest live molecule id (GH #52). Three columns of the same chemistry that BNG gives the identical rate — seeded, assembled in-run, and one whose canonically-first molecule carries the *larger* value — so any disagreement between them is the defect. |
| `species_delete` | ode | **Deleting a species.** A whole-species delete (`A() -> 0`, no `DeleteMolecules`) fires once per matching species, because BNG2 makes the species its reaction center — on homodimers, trimers named through either subunit, the same trimer named `A().B()` and `B().A()`, and dimers that assemble during the run. NFsim counts one per match of the pattern's first molecule, which is 2x on most of these and order-dependent. A partial delete without the keyword (`C(x!1).D(x!1) -> D(x)`) drops only the named molecule and fires only if the remainder forms exactly the product patterns; NFsim refuses it at load. |
| `product_split` | ode | **The `+` between product patterns.** A bond the rule deletes between molecules it maps into different product patterns must leave them in different complexes once all of the rule's operations apply, whatever the number of reactant patterns (issue #102). Two-reactant splits and a split that also makes a bond, on ring pairs (never fire) and single pairs (fire). A rule deleting both bonds of a ring (fires). A deleted bond whose ends stay in one product pattern, held by a ring (fires) or not (never fires). Controls: the one-reactant split, and an energy (Arrhenius) unbinding rule, which RuleMonkey expands at load. NFsim `-bscb` gets six of the eleven columns wrong. |
| `product_split_full` | ode | **The `+` and `.` where a check bond by bond cannot decide them** (issue #104). A molecule deleted out of a complex under `DeleteMolecules`: a `+` between its neighbours on a chain (fires) and a ring (never fires), the same chain with an extra molecule on the deleted one (fires, the extra piece allowed), and a `.` on a ring (fires) and a chain (never fires). A delete without the keyword that leaves as many species as product patterns, but the wrong ones (never fires). A pattern cut into three pieces, under a `+` and under a `.`, on chains and rings. A molecule joined to its pattern only by `.`, on either side of a deleted bond, and held by a bond the rule keeps. |
| `products_never_form` | ode | **A rule whose right-hand side can never be true** (issue #105). Two molecules in one product pattern that start in different reactant patterns, one of them named by no component, or that the rule creates, with no bond added between them; and two molecules of one complex written into different product patterns with nothing deleted between them. BNG2 generates no reaction for any of them. Controls: a created molecule the rule bonds, and a pair the rule does part. |
| `arrhenius_formula` | ode | **The Arrhenius rate law as BNG2's network builds it** (issue #107): `exp(-(Ea0 + phi*dG))` with `dG` the reaction the rule writes, `1-phi` for the reverse half of a reversible rule, and no `RT`. phi = 0.3, and no parameter named `RT`. A reversible rule written as binding and one written as unbinding, unidirectional rules that make and break a bond, and a cooperative context pattern. And rules the modeler labelled `_reverse_Gone`, with no rule `Gone`, and `_reverse_Flip`, beside an unrelated rule `Flip`: rules as written, which take phi like any other, though BNG2 names the reverse half of a reversible rule that way (issue #109). And a cooperative pattern written with `!+`, which the load-time expansion prices exactly. |
| `arrhenius_rt` | ode | **A parameter named `RT` is ordinary** (issue #107): energies in kJ/mol, divided by the model's own `RT` where they are used, as BNG2's `Models2/energy_example1.bngl` does. The same rates as `arrhenius_formula`'s reversible rule and bond-breaking rule. |

## Running

`ctest -R bng_oracle` runs it, or directly:

```bash
python3 tests/bng_oracle/check.py --driver build/release/rm_driver
```

Scoring is a z on the difference of two means, combining RuleMonkey's ensemble
SEM with the reference's own in quadrature, so a stochastic reference is not
treated as exact. Each model fails by tens of sigma against the engine before
its fix, which is the property that makes them worth committing.

## Regenerating

Needs BNG2 (not vendored):

```bash
BNG2=~/path/to/BNG2.pl python3 tests/bng_oracle/generate_references.py
```

This rewrites both `xml/` and `reference/`. The SSA model runs its reps in
parallel and takes a few minutes.
