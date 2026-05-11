# Benchmark: qsp-codegen vs MATLAB SimBiology

This directory contains a reproducible benchmark used to generate the
speed-up table in the JOSS paper.

## What it measures

For the same SBML model, two timing regimes are reported:

1. **Integration only.** MATLAB `tic`/`toc` around `sbiosimulate` inside
   one persistent MATLAB session (after a warm-up call). For C++,
   wall-clock per `qsp_sim` invocation minus a calibrated startup baseline.
   This isolates ODE-solver work from process-launch overhead.
2. **Wall-clock per invocation.** One full `matlab -batch` launch
   (load model + `sbiosimulate` + exit) versus one `qsp_sim` invocation.
   This is what an SBI workflow pays per call when not using a
   persistent MATLAB worker pool — MATLAB's ~6 s startup cost dominates.

Both sides use SUNDIALS as the integrator with matched tolerances
(`reltol=1e-6`, `abstol=1e-9`).

## Prerequisites

- MATLAB with SimBiology (tested on R2023b+).
- A working build of the `qsp_sim` binary (see top-level `README.md`).
- `qsp-codegen` installed in the active Python environment
  (`uv pip install -e .` from the repo root).
- `numpy` for the orchestrator script.

## Reproducing the paper table

From the repo root:

```bash
# 1. Build the toy SimBiology model and export it to SBML.
matlab -batch "run('paper/benchmark/build_model.m')"

# 2. Build the C++ runtime + qsp_sim binary (one-time setup).
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release
cmake --build cpp/build --target qsp_sim -j

# 3. Run the benchmark. ~30-60s plus ~6s per cold MATLAB launch.
python paper/benchmark/run_benchmark.py \
    --qsp-sim cpp/build/qsp_sim \
    --reps 20 \
    --cold-reps 3
```

Outputs:

- `results.csv` — raw per-rep timings (one row per rep, four columns).
- `results_table.md` — formatted Markdown table for pasting into
  `paper.md`.

## Using a different model

The default benchmark model is a small two-compartment PK with linear
elimination plus a Michaelis-Menten clearance term, defined in
[`build_model.m`](build_model.m). To benchmark against a production
QSP model, point both sides at a different SBML:

```bash
python paper/benchmark/run_benchmark.py \
    --sbml /path/to/PDAC_model.sbml \
    --qsp-sim cpp/build/qsp_sim \
    --t-end 365 --reps 20
```

The script does not require dosing or `evolve_to_diagnosis`; it
benchmarks the bare integration of a default-IC scenario, which is
the dominant cost in inference workloads.

## A note on time units

`qsp_sim` defaults to integrating in SI seconds (multiplies the
external `--t-end-days` by 86400 internally) on the assumption that
the model's rate constants are in `1/s`. Production QSP models that
ship explicit `<listOfUnitDefinitions>` in their SBML satisfy this
because `qsp-codegen`'s `_BARE_SI` table converts `1/day`, `1/hr`,
etc. to `1/s` at codegen time.

This benchmark's hand-built model is exported by SimBiology's
`sbmlexport` with no unit definitions, so its rate constants come
through unit-less in `1/day`. `run_benchmark.py` therefore passes
`--time-unit days` to `qsp_sim`, which sets the runtime time factor
to 1.0 and integrates in model-native days. Use `--time-unit
seconds` (or omit, since it's the default) for unit-annotated SBML.
The compile-time `MODEL_UNITS` define still exists as a way to flip
the *default* in builds shipping a particular convention, but the
runtime flag overrides it either way.

## Verifying correctness

The script runs one untimed `qsp_sim --csv-out` invocation per regime
after the timing loops finish, then diffs the resulting trajectory
against MATLAB's via `qsp_codegen.parity.compare`. The script exits
non-zero if any regime fails, so a passing exit code is part of the
benchmark's contract.

## Notes on fairness

- Both sides solve the same SBML file. The MATLAB side imports it via
  `sbmlimport`, the C++ side via `qsp-codegen` + the bundled
  `qsp_sim_core` runtime.
- The first MATLAB `sbiosimulate` call is treated as a warm-up and
  excluded from the median; SimBiology's first-call setup is amortised
  across the rest of the loop, just as it would be in any long-running
  inference job.
- The C++ baseline subtracts a median-of-five `qsp_sim` no-arg launch
  cost. This is the closest analogue to MATLAB's in-process timing
  that does not require modifying `qsp_sim` to accept a `--reps` flag.
- Codegen and CMake build time are explicitly excluded — they are
  one-time costs, not per-call costs.
