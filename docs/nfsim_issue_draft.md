# Draft GitHub issue for `RuleWorld/NFsim`

Paste the body below into a new issue at
<https://github.com/RuleWorld/NFsim/issues/new>.

The repro models (`mlnr.bngl`, `pltr.bngl`, `testcase2a.bngl`,
`rm_tlbr.bngl`) all live in BioNetGen's own validation set, so
no extra files need to be attached.

---

## Title

> Performance: Species observables with count quantifiers (`X()=N`) re-walk every live complex per output sample — O(complexes × observables × members)

## Body

### Summary

`SpeciesObservable::isObservable(Complex*)` (`src/NFcore/observable.cpp:344-393`) is invoked once per `(complex, observable)` pair at every output sample. For models that declare a histogram of size-binned Species observables (e.g. `Species Size_N R()=N` for N=1..300), the per-sample cost grows as **O(C · O · M)** — *C* complexes, *O* observables, *M* avg complex size. With moderate corpus parameters this dominates total wall time.

### Reproducer

`mlnr.bngl` from BioNetGen's own model collection (Faeder/Hlavacek's 5-valent-ligand / 3-valent-receptor benchmark). Declares 300 Species observables of the form `R()=N` for N = 1..300:

```
begin observables
  ...
  Species Size_1   R()=1
  Species Size_2   R()=2
  ...
  Species Size_300 R()=300
end observables
```

```
NFsim -xml mlnr.xml -sim 1000 -oSteps 1000 -seed 1 -bscb -o /tmp/mlnr.gdat
```

Wall time on a 2024 Apple-silicon laptop, NFsim v1.14.3:

| Model       | NFsim wall (median of 3) | Same model in RuleMonkey 3.1 (cleanroom) | Ratio |
|---           | ---                      | ---                                       | ---   |
| `mlnr`       | 79.8 s                   | 3.4 s                                     | 23.7× |
| `pltr`       | 81.3 s                   | 3.3 s                                     | 24.8× |
| `testcase2a` | 30.4 s                   | 2.6 s                                     | 11.8× |
| `rm_tlbr`    | 32.2 s                   | 2.7 s                                     | 12.0× |

All four models share the same shape: 300 `R()=N` observables. The same chemistry without the histogram observables (e.g. `A_plus_B_rings`, `bench_tlbr_yang2008`, `bench_blbr_*`) finishes in ≤200 ms with NFsim — same engine, same machine, same kind of binding chemistry, but only ~half a dozen Species observables instead of 300. The asymmetry is in the observable-evaluation path, not the chemistry.

### Mechanism

`src/NFcore/system.cpp:693-711` (sample-time loop):

```cpp
for(obsIter = speciesObservables.begin(); obsIter != speciesObservables.end(); obsIter++)
    (*obsIter)->clear();

Complex * complex;
allComplexes.resetComplexIter();
while( (complex = allComplexes.nextComplex()) )
{
    if( complex->isAlive() )
    {
        for(obsIter = speciesObservables.begin(); obsIter != speciesObservables.end(); obsIter++)
        {
            match = (*obsIter)->isObservable( complex );
            for (int k=0; k<match; k++) (*obsIter)->straightAdd();
        }
    }
}
```

For each alive complex × each Species observable, `isObservable(complex)` walks the complex's member list to count template matches and apply the count relation. The result is recomputed from scratch every sample, even though no rule fired against most complexes between samples.

### Suggested fix

The count of any rigid template inside a complex changes only when the complex itself mutates (a rule fires AddBond / DeleteBond / AddMolecule / DeleteMolecule on it). The natural shape of the optimization:

1. Add a `dirty` flag on `Complex` set whenever a rule's transformation touches that complex.
2. Per-complex cache the latest computed `isObservable` value for each Species observable (or a more compact per-template count cache).
3. At sample time: only call `isObservable(complex)` for complexes whose `dirty` flag is set; subtract their stale contribution and add the fresh one. Clear the dirty bit.
4. Steady-state cost drops from O(C · O · M) to O(`dirty_complexes` · O), which is typically O(1) per sample step.

Equivalent caching is already used elsewhere in NFsim — `MoleculesObservable` updates incrementally on add/remove events via `MoleculesObservable::add(Molecule*)` / `straightAdd` — so the integration points exist. The change is structural but not invasive: the same per-event hooks that drive `MoleculesObservable` updates can flag `dirty` on the affected complex.

### Context

I'm working on a cleanroom RuleMonkey 3.0 / 3.1 rewrite (MIT-licensed, C++17) that lands an incremental Species-observable tracker along these lines and matches NFsim within stochastic noise on 144 of 173 corpus / feature-coverage / basicmodels models. Happy to share specifics of the implementation if helpful — RM's per-complex Species observable update is at `cpp/rulemonkey/engine.cpp` (the `dirty_cx` / `kSpeciesIncrObs` machinery, ~50 lines). The cleanroom is at <https://github.com/wshlavacek/RuleMonkey> if you want to look at the structure.

Reproducible measurements, scripts, and per-model data are at <https://github.com/wshlavacek/RuleMonkey/blob/main/docs/timing_comparison.md>.

### Notes

- Models without count-relation Species observables show NFsim and RuleMonkey at parity (RM-vs-NFsim median 0.83× across 173 models), so this isn't an across-the-board wall-time complaint — it's specifically the histogram-of-size-bins workload that gets slow.
- Aggregate-size histograms like the ones above are common in receptor-aggregation modeling, gel-formation studies, and any workflow that wants a cluster-size distribution as output. Worth fixing for that audience.

— Bill Hlavacek (`wshlavacek`)
