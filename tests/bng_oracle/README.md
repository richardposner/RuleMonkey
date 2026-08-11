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

**`ssa`** — `generate_network()` + an ensemble of `simulate(method=>"ssa")`
reps, with mean and SEM committed. Needed for a model that deliberately runs at
a handful of complexes, where the ODE prices an `X + X` reaction as `[X]^2`
while the SSA prices it as `N(N-1)` — a factor of 2 at two copies. That is the
ordinary deterministic-limit gap and says nothing about the model; it just means
the deterministic solver is answering a different question than a 4-complex
stochastic model asks. Same tool, same expanded network, stochastic integrator.

## Models

| model | reference | pins |
|---|---|---|
| `context_symmetry` | ode | An untransformed (context) reactant pattern is counted once per matching **complex**, not once per matching molecule — including the `Bind_sym` guard in the other direction, where the rule *does* transform the pattern and the 2x is correct (GH #33). |
| `context_sampler` | ssa | The same slot is **drawn** the way it is **counted**. Counting per complex while sampling per molecule leaves the same-complex rejection mis-weighted, which is a systematic bias wherever complexes hold unequal numbers of context matches (GH #33). |

## Running

`ctest -R bng_oracle` runs it, or directly:

```bash
python3 tests/bng_oracle/check.py --driver build/release/rm_driver
```

Scoring is a z on the difference of two means, combining RuleMonkey's ensemble
SEM with the reference's own in quadrature, so a stochastic reference is not
treated as exact. Both models fail by tens of sigma against the pre-fix engine,
which is the property that makes them worth committing.

## Regenerating

Needs BNG2 (not vendored):

```bash
BNG2=~/path/to/BNG2.pl python3 tests/bng_oracle/generate_references.py
```

This rewrites both `xml/` and `reference/`. The SSA model runs its reps in
parallel and takes a few minutes.
