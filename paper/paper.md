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
    affiliation: 1
  - name: Joel Eliason
    orcid: 0000-0003-2227-8727
    affiliation: 1
  - name: Aleksander S. Popel
    affiliation: 1
affiliations:
  - name: Department of Biomedical Engineering, Johns Hopkins University, Baltimore, MD, USA
    index: 1
date: 6 May 2026
bibliography: paper.bib
---

# Summary

Quantitative systems pharmacology (QSP) models are mechanistic, ordinary differential equation (ODE)-based representations of disease biology and drug action used throughout pharmaceutical development [@Gadkar2016; @Allen2016; @Craig2023]. They are most often constructed in MATLAB's SimBiology, which provides a graphical and SBML-compatible authoring environment but a stiff, MATLAB-bound integrator that becomes a bottleneck when models are exercised in modern inference pipelines such as simulation-based inference (SBI) [@Cranmer2020; @Goncalves2020].

`qsp-codegen` is a Python package that takes a SimBiology-exported SBML Level 2 Version 4 file and emits a complete set of C++ sources for a CVODE-backed standalone simulator [@Hindmarsh2005]. The emitted code, together with a small bundled C++ runtime (`qsp_sim_core`) that the wheel ships and exposes through a CMake package, compiles into a `qsp_sim` executable that integrates the model orders of magnitude faster than the original SimBiology workflow while preserving the model's species, parameters, reactions, and initial state.

# Statement of need

SBML offers a portable specification for systems biology ODE models, and several mature tools generate executable code from SBML, including AMICI [@Frohlich2021], libRoadRunner [@Somogyi2015], SBMLtoODEpy [@Ruggiero2020], and COPASI [@Hoops2006]. These tools are excellent general-purpose simulators, but QSP applications place a few specific demands that motivate a smaller, focused code generator:

1. **SimBiology-exported SBML quirks.** Production QSP models frequently use SimBiology-specific patterns (e.g., `repeatedAssignment` rules, particular unit annotations) that not every general SBML toolchain handles cleanly. `qsp-codegen` is written against the SBML dialect that SimBiology actually emits.
2. **Tight coupling to a hand-written C++ driver.** QSP simulations typically include model-specific machinery that lives outside the ODE itself: scenario-aware dosing schedules, a pre-treatment burn-in phase that evolves the patient until a diagnostic event, and consumer-specific initial-condition setup. These are awkward to express through a general-purpose simulator's API but natural in a thin C++ driver. `qsp-codegen` is designed to plug into such a driver via a stable hook interface (`evolve_to_diagnosis`, declared in `qsp_sim_core/model_hooks.h`) rather than around it.
3. **ABI stability for iterative inference workflows.** When QSP models are coupled to neural posterior estimation or other SBI methods, the simulator is rebuilt and called millions of times across model iterations. Keeping the generated ABI minimal and predictable (one `ODE_system` class plus a parameter container) lets consumer build systems cache aggressively and avoid spurious recompilation.
4. **Optional dependency footprint.** General SBML simulators bring substantial dependency stacks. `qsp-codegen` requires only `sympy` and `numpy` at code-generation time; the runtime needs only CVODE.

`qsp-codegen` is not a competitor to AMICI or libRoadRunner for general systems-biology workflows. It occupies a narrower niche: it is the SBML-to-simulator step inside a QSP-specific stack that includes the `qsp-hpc-tools` HPC orchestration layer [@qsp_hpc_tools] and downstream Bayesian inference. It is currently used in production for a pancreatic cancer QSP model, where the generated C++ simulator yields per-call speedups of roughly 25-87 times over the equivalent SimBiology integration on representative parameter sets, closing the gap that previously made full SBI workflows impractical in this setting.

# Design and key features

## Code generation

`qsp-codegen` parses the input SBML file with `libsbml`, normalizes SimBiology-specific constructs (assignment rules, repeated assignments, unit factors, conserved moiety patterns), and emits:

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

`qsp-codegen` ships a parity harness in `qsp_codegen.parity` that runs a SimBiology trajectory export side-by-side against the C++ simulator over a user-supplied parameter sweep and reports per-species relative error. This makes regressions in the code generator easy to detect when SimBiology models are revised.

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

This work was supported by the National Institutes of Health. The authors thank the Maryland Advanced Research Computing Center (MARCC) for HPC resources used to validate the generated simulator at scale.

# References
