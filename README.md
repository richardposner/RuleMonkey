# RuleMonkey

A network-free stochastic simulator for rule-based biochemical models written
in [BNGL](https://bionetgen.org/). RuleMonkey reads a BioNetGen-generated XML
model and runs an exact discrete-event simulation, producing observable
trajectories in `.gdat` format.

This repository contains a cleanroom C++ implementation of RuleMonkey targeting
integration with [BNGsim](https://github.com/...). The original C
implementation by Yang, Monine, Faeder, and Hlavacek (2008) lives in repository
history prior to this release.

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

## Test corpora

| Suite | Models | Purpose |
|-------|--------|---------|
| `tests/models/feature_coverage/` | 51+ | BNGL feature coverage (invariants + golden values) |
| `tests/models/corpus/` | 71 | Real-world models (efficiency + correctness) |
| `tests/models/nfsim_basicmodels/` | 72 | NFsim-parity regression suite (Sneddon et al.) |

Reference trajectories under `tests/reference/nfsim/` are 100-replicate ensemble
means and standard deviations regenerated 2026-04-10 from a known NFsim build;
see `tests/reference/nfsim/PROVENANCE.md`.

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
