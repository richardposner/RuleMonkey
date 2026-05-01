# `.gdat` output format

`rm_driver` writes its trajectory output to `stdout` in `.gdat`
format — the same shape NFsim produces.  This document is the
authoritative reference for that format.  CPU time goes to `stderr`.

## File shape

One header line followed by N data lines.  Every line is tab-separated.

```
# time<TAB>obs_1<TAB>obs_2<TAB>...<TAB>obs_K
t_0<TAB>v_1,0<TAB>v_2,0<TAB>...<TAB>v_K,0
t_1<TAB>v_1,1<TAB>v_2,1<TAB>...<TAB>v_K,1
...
t_N<TAB>v_1,N<TAB>v_2,N<TAB>...<TAB>v_K,N
```

- The header line begins with a literal `#` followed by a single
  space, then `time`, then a tab-separated list of observable names
  in the order they were declared in the BNGL `begin observables`
  block.
- Observable names are emitted verbatim from the XML's
  `<Observable name="...">` attribute.  RM never renames or
  re-encodes them.
- Every data line has exactly `1 + K` tab-separated fields, where
  `K` is the number of observables.
- For a `rm_driver model.xml T N` invocation, the file has
  **`N + 1` data lines** — samples at `t_start`, `t_start + dt`, …,
  `t_end`, where `dt = (t_end − t_start) / N`.

## Numeric precision

Time and observable values are emitted via the C++ `<<` operator
with `std::cout`'s default precision: **6 significant figures**.
This matches NFsim's default and is sufficient for the parity
benchmarks in this repo.  External consumers who need more digits
should:

- Round their reference data to six sig figs at compare time, or
- Patch `cpp/cli/rm_driver.cpp` to set `std::cout << std::setprecision(N)`
  for the desired precision (the `setprecision(17)` round-trip
  setting is appropriate for bit-exact regeneration).

The in-process C++ API (`rulemonkey::Result`) carries full
double-precision values; the precision loss is purely a property of
the CLI driver's text output.

## Standard error

```
ERROR: <message>          (only if loading or running failed)
WARN:  <unsupported feature> (during XML parse, if features at
                              Severity::Warn are present)
ERROR: <unsupported feature> (during XML parse, if features at
                              Severity::Error are present)
<cpu_seconds>             (last line, on success)
```

A successful run terminates with exit code 0 and prints a single
floating-point seconds value as the last `stderr` line.  If
unsupported BNGL features were detected, those are reported on
preceding `stderr` lines (Warn-severity always; Error-severity if
present, exit 2 unless `--ignore-unsupported`).

## Parsing recipe

```python
# Python — handles header normalisation and float conversion
def parse_gdat(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    headers = lines[0].lstrip("#").strip().split("\t")  # ['time', 'A', ...]
    rows = [list(map(float, l.split("\t"))) for l in lines[1:]]
    return headers, rows

headers, rows = parse_gdat("traj.gdat")
times = [r[0] for r in rows]
A = [r[headers.index("A")] for r in rows]
```

Two header-parsing notes:

- The `#` may be followed by a single space (RM emits `# time...`),
  or by no space (some NFsim builds emit `#time...`); `lstrip("#").strip()`
  handles both.
- The first column is always called `time` (lowercase) regardless of
  whether the BNGL declared an observable called `time`.

## Edge cases

- **Zero observables**: a model with no `begin observables` block
  produces a header of just `# time\n` and rows with a single time
  field per line.  Parsing the columns by tab still works.
- **Absorbing state**: when the system reaches zero total propensity
  (no rule can fire), RM fills the remaining sample points by
  repeating the last computed state, so the row count is unchanged.
- **`n_steps == 0`**: not supported on the CLI; pass `n_steps >= 1`.
  The in-process API treats `n_points == 0` as "step to t_end without
  recording" via `step_to`.

## See also

- [`quickstart.md`](quickstart.md) — full command-line walkthrough.
- [`include/rulemonkey/types.hpp`](../include/rulemonkey/types.hpp)
  — the in-process `Result` struct that holds full-precision values.
