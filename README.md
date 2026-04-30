# RuleMonkey

A network-free stochastic simulator for rule-based biochemical models written
in [BNGL](https://bionetgen.org/). RuleMonkey reads a BioNetGen-generated XML
model and runs an exact discrete-event simulation, producing observable
trajectories in `.gdat` format.

This is the **RuleMonkey 3.x** cleanroom C++17 rewrite, designed to be
embedded as an in-process simulation kernel.  The engine exposes a
small C++ API (`include/rulemonkey/simulator.hpp`) that host
applications — most notably the BNGsim simulation engine — can drive
without spawning subprocesses or going through the file system.

The legacy C implementation (RuleMonkey 2.0.25), introduced in
Colvin J, Monine MI, Gutenkunst RN, Hlavacek WS, Von Hoff DD, Posner RG.
*RuleMonkey: software for stochastic simulation of rule-based models.*
BMC Bioinformatics 11:404 (2010), is preserved in the parent fork's
git history.

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

## Embedding (C++ API)

RuleMonkey is designed to be embedded directly into a host application.
The public surface is two headers — `include/rulemonkey/simulator.hpp`
and `include/rulemonkey/types.hpp` — and one library target,
`rulemonkey` (built by this project's CMake).

Stateless one-shot simulation:

```cpp
#include <rulemonkey/simulator.hpp>

rulemonkey::RuleMonkeySimulator sim("model.xml");
sim.set_block_same_complex_binding(true);
sim.set_param("kf", 0.5);

rulemonkey::TimeSpec ts{ /*t_start=*/0.0, /*t_end=*/10.0, /*n_points=*/101 };
auto result = sim.run(ts, /*seed=*/42);

for (size_t t = 0; t < result.n_times(); ++t) {
  // result.time[t], result.observable_data[obs_idx][t]
}
```

Stateful session for multi-segment workflows (initialize → simulate →
mutate → simulate → save state):

```cpp
sim.initialize(/*seed=*/42);
auto seg1 = sim.simulate(0.0, 100.0, 101);
sim.add_molecules("L", 500);                     // mid-run perturbation
auto seg2 = sim.simulate(100.0, 200.0, 101);
sim.save_state("checkpoint.bin");
sim.destroy_session();
```

See `include/rulemonkey/simulator.hpp` for the complete contract,
`include/rulemonkey/types.hpp` for the result/time-spec/feature-warning
structs, and `examples/embed.cpp` for a minimal compilable example
(plus a `find_package` / `add_subdirectory` CMake snippet in
`examples/CMakeLists.txt`).

## Model coverage

[`docs/model_semantics.md`](docs/model_semantics.md) is the reference for
"will my BNGL model run on RM?" — it lists every supported BNGL
construct, every Tier-0 refusal (compartments, Arrhenius / Sat / Hill /
FunctionProduct rate laws, population types, multi-molecule fixed
species), and the best-effort warnings the engine emits at load time.

[`docs/timing_comparison.md`](docs/timing_comparison.md) is the
companion performance reference: a 173-model RM-vs-NFsim wall-time
comparison with a "where does each engine win?" breakdown.

## License

MIT. See [LICENSE](LICENSE).
