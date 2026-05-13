# RM vs NFsim — time-horizon scaling (issue #5)

Tests whether RuleMonkey 3.x exhibits horizon-dependent slowdown (the artifact RGP observed in original RM).  Each row is the median of 3 reps; `n_steps=400` held constant across the ladder so observable-recording cost stays fixed and the comparison isolates per-event engine cost.

**Reading the table:** if per-event time stays roughly constant down each model's column, the engine scales linearly with horizon (no degradation).  A *rising* per-event time at longer horizons would be the original-RM signature.

Generated: `2026-05-13 00:21:03 MDT` (RuleMonkey build: `/Users/wish/Code/RuleMonkey/build/release/rm_driver`)

## `Lipniacki2006`

| t_end | RM wall | RM events | RM ns/event | NFsim wall | NFsim events | NFsim ns/event | RM speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 60.2 ms | 40,100 | 1.5 µs | 32.8 ms | 39,984 | 821 ns | 0.55× |
| 200 | 538 ms | 401,361 | 1.3 µs | 208 ms | 401,676 | 518 ns | 0.39× |
| 2000 | 5430 ms | 4,014,294 | 1.4 µs | 1968 ms | 4,014,759 | 490 ns | 0.36× |
| 20000 | 54.7 s | 40,155,223 | 1.4 µs | 19.6 s | 40,141,736 | 489 ns | 0.36× |

## `Samoilov2005_FutileCycle`

| t_end | RM wall | RM events | RM ns/event | NFsim wall | NFsim events | NFsim ns/event | RM speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.018 | 78.0 ms | 21,647 | 3.6 µs | 37.5 ms | 21,503 | 1.7 µs | 0.48× |
| 0.18 | 744 ms | 218,559 | 3.4 µs | 218 ms | 217,779 | 999 ns | 0.29× |
| 1.8 | 7614 ms | 2,186,005 | 3.5 µs | 1961 ms | 2,184,772 | 898 ns | 0.26× |
| 18 | 76.1 s | 21,853,591 | 3.5 µs | 18.9 s | 21,859,045 | 865 ns | 0.25× |

## `HBF1998_brusselator`

| t_end | RM wall | RM events | RM ns/event | NFsim wall | NFsim events | NFsim ns/event | RM speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 65.5 ms | 41,126 | 1.6 µs | 40.6 ms | 41,243 | 985 ns | 0.62× |
| 100 | 492 ms | 357,527 | 1.4 µs | 217 ms | 358,141 | 605 ns | 0.44× |
| 1000 | 4719 ms | 3,532,640 | 1.3 µs | 1959 ms | 3,534,818 | 554 ns | 0.42× |
| 10000 | 46.2 s | 35,238,009 | 1.3 µs | 18.9 s | 35,256,567 | 537 ns | 0.41× |

## `tlbr_macken1982`

| t_end | RM wall | RM events | RM ns/event | NFsim wall | NFsim events | NFsim ns/event | RM speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 60.3 ms | 361 | 167.0 µs | 42.5 ms | 343 | 123.9 µs | 0.71× |
| 300 | 114 ms | 2,711 | 42.0 µs | 82.3 ms | 2,725 | 30.2 µs | 0.72× |
| 3000 | 708 ms | 25,774 | 27.5 µs | 474 ms | 25,837 | 18.3 µs | 0.67× |
| 30000 | 6291 ms | 256,869 | 24.5 µs | 4304 ms | 256,944 | 16.8 µs | 0.68× |

## `blbr_dembo1978`

| t_end | RM wall | RM events | RM ns/event | NFsim wall | NFsim events | NFsim ns/event | RM speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 30 | 163 ms | 2,663 | 61.4 µs | 144 ms | 2,679 | 53.9 µs | 0.88× |
| 300 | 440 ms | 27,501 | 16.0 µs | 324 ms | 27,442 | 11.8 µs | 0.74× |
| 3000 | 3120 ms | 272,698 | 11.4 µs | 2247 ms | 271,989 | 8.3 µs | 0.72× |
| 30000 | 29.4 s | 2,719,927 | 10.8 µs | 21.7 s | 2,720,959 | 8.0 µs | 0.74× |

## Per-event time: smallest vs largest horizon

Drift > 1.5× (slower at longer horizon) would suggest engine cost growing with simulation state or time, not just with event count.

| Model | Engine | ns/evt @ short | ns/evt @ long | Drift (long/short) |
|---|---|---:|---:|---:|
| `Lipniacki2006` | rm | 1502 ns | 1363 ns | 0.91× |
| `Lipniacki2006` | nfsim | 821 ns | 489 ns | 0.60× |
| `Samoilov2005_FutileCycle` | rm | 3604 ns | 3481 ns | 0.97× |
| `Samoilov2005_FutileCycle` | nfsim | 1745 ns | 865 ns | 0.50× |
| `HBF1998_brusselator` | rm | 1592 ns | 1310 ns | 0.82× |
| `HBF1998_brusselator` | nfsim | 985 ns | 537 ns | 0.55× |
| `tlbr_macken1982` | rm | 167016 ns | 24490 ns | 0.15× |
| `tlbr_macken1982` | nfsim | 123943 ns | 16751 ns | 0.14× |
| `blbr_dembo1978` | rm | 61351 ns | 10804 ns | 0.18× |
| `blbr_dembo1978` | nfsim | 53901 ns | 7990 ns | 0.15× |

## Notes on methodology

- Wall time is end-to-end subprocess wall (process spawn + XML load + simulate + observable record + exit).  Same for both engines, so ratios are fair; absolute numbers shift per machine.
- RM event count comes from `RM_PRINT_TIMING=1` stderr (`events=N null=N`); NFsim event count from stdout (`You just simulated N reactions`).  Both count accepted SSA events.
- `-bscb` (Bind-Single-Copy-Bond) enabled on both engines.
- `n_steps` held constant; per-event cost should track engine work, not observable I/O.
- Each rep uses a distinct seed (1..N), so event totals differ between reps; medians are reported.
