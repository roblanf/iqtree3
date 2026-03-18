# sCFL simulation validation (Issue #137)

Simulation study to validate the fix for `--scfl` gap detection with
rate-heterogeneous models. See the issue for full results and figures.

## Prerequisites

- Built `iqtree3` binary at `build/iqtree3`
- Python 3 with `pandas`, `matplotlib`, `numpy`
- [AMAS](https://github.com/marekborowiec/AMAS) (`AMAS.py` on PATH)
- GNU parallel (optional, for faster analysis runs)

## Quick start

```bash
# 1. Simulate alignments (11 tree1% levels x 5 gap levels = 55 alignments)
bash tests/scfl_simulation/simulate.sh

# 2. Run all sCF analyses (1,100 runs, ~2-3 hours single-threaded)
bash tests/scfl_simulation/run_analyses.sh

# 3. Generate figures and validation checks
python3 tests/scfl_simulation/analyze.py
```

## What it does

Simulates 10,000-site DNA alignments on an 8-taxon tree under JC,
mixing known proportions of sites from the true tree and an NNI
alternative, then injects gaps at 0-20%. Runs sCF with 4 models
(GTR, GTR+I, GTR+I+G8, GTR+I+R4) and 4 methods (scf, scfl, scflg,
scfl --safe) on each alignment.

## Expected results

- **sN** depends on gap%, not tree1% (flat lines in Fig 2)
- **sCF** depends on tree1%, not gap% (monotonic increase in Fig 4)
- **sCFLG** deviates: gaps bias sCF downward (expected, by design)
- **sCFL = sCFL --safe**: max |diff| < 0.1 (validates SAFE_NUMERIC handling)

## Comparing old vs new binary

To test against the buggy binary, build before applying the fix:

```bash
cp build/iqtree3 build/iqtree3-old
```

`run_analyses.sh` will use `build/iqtree3-old` for comparison runs
if it exists. `analyze.py` reports the old vs new sN ratio for
multi-category models (expected: 10-100x improvement).
