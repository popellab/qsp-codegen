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
    affiliation: 1
  - name: Joel Eliason
    orcid: 0000-0003-2227-8727
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
  - name: Johns Hopkins University, Baltimore, MD, USA  # TODO: confirm department before submission
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
3. **Burn-in state reuse across parameter sweeps.** A common QSP workflow is to evolve a baseline (e.g., healthy tissue) until a diagnostic event, then run many treatment scenarios from that state for the same parameter vector. The bundled runtime serializes the post-burn-in CVODE integrator state to a small parameter-hashed binary cache, so the burn-in is paid once per parameter vector rather than once per scenario. General SBML simulators expose pre-equilibration or steady-state routines but not a portable checkpoint that can be reused across processes and HPC jobs.
4. **Reproducibility and distribution of QSP simulators.** SimBiology is the standard authoring environment for QSP models but is not always the most convenient deployment target — open-source CI, reviewers attempting to reproduce results, and downstream researchers using a different toolchain all benefit from a self-contained build that can be exercised independently of the authoring environment. `qsp-codegen`'s output is plain C++ depending only on CVODE, so a model authored interactively in SimBiology can be distributed as an executable that any collaborator can build and run.
5. **Embeddability inside larger C++ codebases.** The generated `ODE_system` is a plain C++ class with a small header surface, and `qsp_sim_core` links statically against CVODE alone — no Python interpreter and no LLVM JIT in the dependency closure. The codegen and runtime in this package were extracted from the SPQSP_PDAC GPU agent-based model [@spqsp_pdac], where a host-side `MolecularModelCVode<ODE_system>` integrates a lymph-central QSP submodel coupled each step to a FLAME GPU 2 cell-level ABM; the emitted class shape is set by that embedded use case. libRoadRunner and AMICI both bring substantially heavier runtimes (LLVM JIT and a Python-oriented Solver/Model/ReturnData stack respectively) that are awkward to drop into a C++-native consumer.

`qsp-codegen` is not a competitor to AMICI or libRoadRunner for general systems-biology workflows. It occupies a narrower niche: it is the SBML-to-simulator step inside a QSP-specific stack that includes the `qsp-hpc-tools` HPC orchestration layer [@qsp_hpc_tools] and downstream Bayesian inference.

The reproducible benchmark in `paper/benchmark/` exports a 25-compartment, 73-reaction model from SimBiology and runs it under both backends with matched SUNDIALS tolerances (`reltol=1e-6`, `abstol=1e-9`) over a 365-day horizon. Trajectories are checked for agreement (worst per-species relative error 1%) before timings are reported. Three regimes are shown (medians of 30 repetitions; cold-start medians of 3, p25–p75 in parentheses):

| Mode | MATLAB SimBiology (s) | qsp-codegen / CVODE (s) | Speedup |
|---|---|---|---|
| Integration only, no dosing | 0.011 (0.011–0.012) | 0.0058 (0.0056–0.0059) | 2.0× |
| Integration only, 6-bolus schedule | 0.017 (0.016–0.018) | 0.0069 (0.0067–0.0074) | 2.5× |
| Wall-clock per invocation, no dosing | 8.271 (8.251–8.358) | 0.0084 (0.0083–0.0100) | ≈980× |

The integration-only regime isolates ODE-solver work; both engines run CVODE with matched tolerances, so the ~2× advantage at this size reflects per-call effects we did not decompose (compiled vs. accelerator-JIT'd RHS dispatch, analytical vs. finite-difference Jacobian, state-vector marshalling) and is expected to grow with model size and stiffness. The wall-clock regime is what an SBI workflow pays per call when not using a persistent MATLAB worker pool — MATLAB's process startup dominates and the gap stretches to three orders of magnitude end-to-end.

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
- A `qsp_sim` driver (`qsp_sim_main.cpp`) that parses YAML scenario and drug-metadata files, applies segmented dose schedules, writes CSV or binary trajectories, and reads/writes a packed evolve-cache used to share burn-in state across scenarios.
- A `model_hooks.h` interface declaring the `evolve_to_diagnosis` extension point, with a default no-op implementation in the static library so models without a burn-in phase need do nothing. Models that do require one (for example, evolving healthy tissue until a tumor reaches a target diameter) provide a strong definition that wins at link time, leaving the generated ABI untouched.

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

This work was supported by the National Institutes of Health and the Lustgarten Foundation. The authors thank the Maryland Advanced Research Computing Center (MARCC) for HPC resources used to validate the generated simulator at scale.

# References
