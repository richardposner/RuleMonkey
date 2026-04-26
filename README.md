# RuleMonkey

A network-free stochastic simulator for rule-based biochemical models written
in [BNGL](https://bionetgen.org/). RuleMonkey reads a BioNetGen-generated XML
model and runs an exact discrete-event simulation, producing observable
trajectories in `.gdat` format.

This is **RuleMonkey 3.0**, a cleanroom C++17 rewrite targeting integration
with [BNGsim](https://github.com/...). The legacy C implementation
(RuleMonkey 2.0.25), introduced in Colvin J, Monine MI, Gutenkunst RN,
Hlavacek WS, Von Hoff DD, Posner RG. *RuleMonkey: software for stochastic
simulation of rule-based models.* BMC Bioinformatics 11:404 (2010), is
preserved in the parent fork's git history.

## Building

Requires CMake ≥ 3.25, a C++17 compiler, and Ninja.

```bash
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The `rm_driver` executable is built at `build/release/rm_driver`.

## Usage

```bash
build/release/rm_driver <model.xml> <t_end> <n_steps> [seed] [-no-bscb] \
    [--ignore-unsupported] [--save-state <path>] [--load-state <path>] \
    [--t-start <time>]
```

Output is tab-separated `.gdat` on stdout (header line beginning with `#`,
then time + observable columns). CPU time is printed to stderr.

RuleMonkey defaults to strict BNGL semantics (`block_same_complex_binding`,
equivalent to NFsim's `-bscb` flag). Pass `-no-bscb` for compatibility with
NFsim runs that disabled this check.

## Layout

```
cpp/rulemonkey/      Core engine sources
cpp/cli/             rm_driver CLI
include/rulemonkey/  Public C++ API headers
tests/
  cpp/               Compiled smoke tests
  models/            BNGL test corpora (feature_coverage, real-world, basicModels)
  reference/         Gold-standard reference trajectories from NFsim and earlier RM
harness/             Python benchmarking + validation drivers
```




Reference trajectories under `tests/reference/nfsim/` are 100-replicate ensemble
means and standard deviations regenerated 2026-04-10 from a known NFsim build,
**with two exceptions** (`toy_jim`, `rm_tlbr_rings`) where NFsim was found to
produce incorrect output and the reference was regenerated from a hand-rolled
Gillespie SSA. See [`tests/reference/nfsim/PROVENANCE.md`](tests/reference/nfsim/PROVENANCE.md)
for details and the regeneration scripts under `harness/ssa/`.

## Public API

```cpp
#include <rulemonkey/simulator.hpp>

rulemonkey::RuleMonkeySimulator sim("model.xml");
sim.set_param("k1", 0.5);
sim.initialize(/*seed=*/42);
auto result = sim.simulate(0.0, 10.0, 100);
```

See `include/rulemonkey/simulator.hpp` for the complete contract.

## License

MIT. See [LICENSE](LICENSE).
