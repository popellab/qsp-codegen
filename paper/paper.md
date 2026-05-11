---
title: 'qsp-codegen: SBML to C++ code generation and a CVODE driver for quantitative systems pharmacology models'
tags:
  - Python
  - C++
  - quantitative systems pharmacology
  - SBML
  - code generation
  - CVODE
  - ODE solvers
authors:
  - name: Chase Christenson
    orcid: 0000-0000-0000-0000  # TODO: fill in before submission
    equal-contrib: true
    affiliation: 1
  - name: Joel Eliason
    orcid: 0000-0003-2227-8727
    equal-contrib: true
    affiliation: 1
  - name: Atul Deshpande
    orcid: 0000-0000-0000-0000  # TODO: fill in before submission
    affiliation: 2
  - name: Aleksander S. Popel
    orcid: 0000-0000-0000-0000  # TODO: fill in before submission
    affiliation: 1
affiliations:
  - name: Department of Biomedical Engineering, Johns Hopkins University, Baltimore, MD, USA
    index: 1
  - name: Department of Oncology, Johns Hopkins University School of Medicine, Baltimore, MD, USA
    index: 2
date: 6 May 2026
bibliography: paper.bib
---

# Summary

Quantitative systems pharmacology (QSP) models are mechanistic, ordinary differential equation (ODE)-based representations of disease biology and drug action used throughout pharmaceutical development [@Gadkar2016; @Allen2016; @Craig2023]. They are most often constructed in MATLAB's SimBiology, which provides a graphical and SBML-compatible authoring environment but a stiff, MATLAB-bound integrator that becomes a bottleneck when models are exercised in modern inference pipelines such as simulation-based inference (SBI) [@Cranmer2020; @Goncalves2020].

`qsp-codegen` is a Python package that takes a SimBiology-exported SBML Level 2 Version 4 file and emits a complete set of C++ sources for a CVODE-backed standalone simulator [@Hindmarsh2005]. The emitted code, together with a small bundled C++ runtime (`qsp_sim_core`) that the wheel ships and exposes through a CMake package, compiles into a `qsp_sim` executable that preserves the model's species, parameters, reactions, and initial state and substantially reduces per-call simulation cost relative to the original SimBiology workflow (see the benchmark below).

# Statement of need

SBML offers a portable specification for systems biology ODE models, and several mature tools generate executable code from SBML, including AMICI [@Frohlich2021], libRoadRunner [@Somogyi2015], SBMLtoODEpy [@Ruggiero2019], and COPASI [@Hoops2006]. These tools are excellent general-purpose simulators, but QSP applications place a few specific demands that motivate a smaller, focused code generator:

1. **SimBiology-exported SBML quirks.** Production QSP models frequently rely on patterns that are idiomatic in SimBiology but awkward for general SBML toolchains: dotted `Compartment.Species` identifiers, MathML where `max`/`min` are emitted as `<ci>max</ci>` identifier nodes rather than the standard `<max/>` operator, and `initialAssignment` rules that override XML `initialAmount`/`initialConcentration` values and must be evaluated in concentration space before being converted back to amounts for the integrator. `qsp-codegen` is written against the SBML dialect that SimBiology actually emits and normalizes these constructs explicitly.
2. **Optional dependency footprint.** General SBML simulators bring substantial dependency stacks. `qsp-codegen` requires only `sympy` and `numpy` at code-generation time; the runtime needs only CVODE.
3. **Reproducibility and distribution of QSP simulators.** SimBiology is the standard authoring environment for QSP models but is not always the most convenient deployment target — open-source CI, reviewers attempting to reproduce results, and downstream researchers using a different toolchain all benefit from a self-contained build that can be exercised independently of the authoring environment. `qsp-codegen`'s output is plain C++ depending only on CVODE, so a model authored interactively in SimBiology can be distributed as an executable that any collaborator can build and run.
4. **Embeddability inside larger C++ codebases.** The generated `ODE_system` is a plain C++ class with a small header surface, and `qsp_sim_core` links statically against CVODE alone — no Python interpreter and no LLVM JIT in the dependency closure. The codegen and runtime were extracted from a larger GPU agent-based cancer model [@spqsp_pdac], where a C++ host integrates the QSP submodel each step alongside a cell-level agent simulation; the emitted class is shaped by that embedded use case. libRoadRunner and AMICI both ship heavier runtimes (an LLVM JIT and a Python-facing API respectively) that are awkward to drop into a C++-native consumer.

`qsp-codegen` is not a competitor to AMICI or libRoadRunner for general systems-biology workflows. It occupies a narrower niche: it is the SBML-to-simulator step inside a QSP-specific stack that includes the `qsp-hpc-tools` HPC orchestration layer [@qsp_hpc_tools] and downstream Bayesian inference.

The reproducible benchmark in `paper/benchmark/` exports a 25-compartment, 73-reaction model from SimBiology and runs it under both backends with matched SUNDIALS tolerances (`reltol=1e-6`, `abstol=1e-9`) over a 365-day horizon. Both engines emit at their CVODE-internal step grid (MATLAB with `OutputTimes = []`, qsp-codegen with `CV_ONE_STEP` and a 2.4 h cadence floor). For the agreement check, MATLAB is then re-run with `OutputTimes = cpp_times` so trajectories are diffed row-by-row at identical timepoints with no interpolation; both regimes pass at `rtol=5e-3`, `atol=1e-6` before timings are reported. Three regimes are shown (medians of 30 repetitions; cold-start medians of 3, p25–p75 in parentheses):

| Mode | MATLAB SimBiology (s) | qsp-codegen / CVODE (s) | Speedup |
|---|---|---|---|
| Integration only, no dosing | 0.010 (0.009–0.011) | 0.0001 (0.0000–0.0002) | ≈130× |
| Integration only, 6-bolus schedule | 0.016 (0.015–0.017) | 0.0018 (0.0017–0.0019) | 8.8× |
| Wall-clock per invocation, no dosing | 8.467 (8.426–8.513) | 0.0091 (0.0085–0.0097) | ≈930× |

The integration-only regime isolates ODE-solver work. Both engines run CVODE with matched tolerances, so the one-to-two-order-of-magnitude advantage at this model size reflects per-call overhead we did not decompose, and is expected to grow with model size and stiffness. The no-dose C++ figure sits near the resolution of subprocess wall-clock timing minus a calibrated `qsp_sim` startup baseline; the ratio is reported to a leading digit rather than to three. The wall-clock regime is what an SBI workflow pays per call when not using a persistent MATLAB worker pool — MATLAB's process startup dominates and the gap stretches to three orders of magnitude end-to-end.

# Design and key features

## Code generation

`qsp-codegen` parses the input SBML file with `libsbml`, normalizes SimBiology-specific constructs (assignment rules, repeated assignments, and unit-factor conversions), and emits:

- `QSP_enum.h`: a strongly-typed enum of state variables and parameters.
- `ODE_system.h`, `ODE_system.cpp`: the ODE right-hand side and a CVODE-compatible Jacobian.
- `QSPParam.h`, `QSPParam.cpp`: a typed parameter container with default values from the SBML model.
- `qsp_params_xml_snippet.xml`: a parameter-XML fragment that the companion `qsp-refresh-param-xml` tool merges into a consumer's `param_all.xml` so that downstream parameter overrides remain centrally managed.

When `sympy` is available, the Jacobian is derived analytically and run through common-subexpression elimination to keep the emitted C++ compact; otherwise the runtime falls back to a finite-difference Jacobian.

## Bundled C++ runtime (`qsp_sim_core`)

The wheel ships a model-agnostic C++ runtime that consumer projects pull in via CMake by calling `python -m qsp_codegen.cmake --prefix` and then `find_package(qsp_sim_core CONFIG REQUIRED)`. The runtime provides:

- `CVODEBase` and `MolecularModelCVode` integrator wrappers around SUNDIALS/CVODE.
- `ParamBase`: a base class for the generated parameter container.
- A `qsp_sim` driver (`qsp_sim_main.cpp`) that parses YAML scenario and drug-metadata files, applies segmented dose schedules, writes CSV or binary trajectories, and reads/writes a packed evolve-cache used to share burn-in state across scenarios. The driver runs CVODE in `CV_ONE_STEP` mode and emits a row whenever the elapsed solver time since the last dump exceeds a configurable cadence floor (`--min-cadence-hours`, default 4 h), plus at every dose boundary and at `t_end`. This preserves the integrator's adaptive stepping end-to-end rather than forcing a fixed-grid output.
- A `model_hooks.h` callback named `evolve_to_diagnosis` that a consumer can override to integrate from a known starting state until a model-defined event (e.g. growing a tumor to a target diameter) before the main scenario begins. Models that don't need this phase use the default no-op.

## Parity harness

`qsp-codegen` ships a parity harness in `qsp_codegen.parity` that runs a SimBiology trajectory export side-by-side against the C++ simulator and reports per-species relative error. This makes regressions in the code generator easy to detect when SimBiology models are revised.

# Typical usage

```bash
# Generate C++ from an SBML export and refresh the parameter XML:
qsp-codegen --sbml model/PDAC_model.sbml --out-dir cpp/qsp/ode/

qsp-refresh-param-xml \
    --snippet cpp/qsp/ode/qsp_params_xml_snippet.xml \
    --xml cpp/sim/resource/param_all.xml
```

The emitted sources, plus the `qsp_sim_core` runtime located via CMake, build a `qsp_sim` binary that the `qsp-hpc-tools` `CppSimulator` invokes on an HPC cluster.

# Acknowledgements

Chase Christenson and Joel Eliason contributed equally to this work and are listed as co-first authors. This work was supported by the National Institutes of Health and the Lustgarten Foundation. The authors thank the Maryland Advanced Research Computing Center (MARCC) for HPC resources used to validate the generated simulator at scale.

# References
