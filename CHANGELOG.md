# Changelog

All notable changes to this project will be documented in this file. The
format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `qsp-codegen verify` subcommand: end-to-end self-test that generates C++,
  scaffolds and builds a minimal `qsp_sim`, and compares its trajectories
  against a MATLAB SimBiology reference (PASS/FAIL). The CLI now uses
  `generate`/`verify` subcommands; the legacy `qsp-codegen --sbml … --out-dir …`
  form still routes to `generate`.
- `generate` now also emits a complete, ready-to-run `param_all.xml` (ICs +
  parameters from the SBML), so the generated `qsp_sim` runs with no manual
  `<Param>` wrapping or `qsp-refresh-param-xml` merge step.
- SBML `<functionDefinition>` (user lambda) inlining at call sites.
- Codegen-time validation of generated C++: catches `AUX_VAR_*`
  use-before-definition (a dependency-ordering bug class) with a located,
  human-readable error instead of an opaque downstream C++ compiler error.
- Tests for the MathML→C++ converter, generated-code validation,
  `functionDefinition` inlining, and SBML event parsing.

### Fixed

- MathML→C++ converter now handles SimBiology's `<ci>nthroot</ci>` (and a
  broader set of `<ci>`-exported operators: `log2`/`log10`, `ceiling`, and the
  hyperbolic / inverse-trig families). Previously emitted a `/* unknown op */`
  comment that crashed the Jacobian `sympify`. Unrecognized operators now
  raise a clear, named error.
- Dependency-ordering bug for concentration assignment-rule species in a
  dynamic-volume compartment (e.g. a drug concentration in a growing tumor
  compartment): the compartment-volume temporary is now seeded into the
  init / `update_y_other` emission closures, fixing an undefined `AUX_VAR_*`
  reference in the generated C++.

- `qsp_sim` runtime flags `--time-unit days|seconds` and `--time-factor <N>`
  to override the compile-time time-scaling default per invocation. Lets a
  single binary integrate both unit-annotated SBML (default, runs in SI
  seconds) and SimBiology's vanilla `sbmlexport` output (`--time-unit days`).
- `paper/benchmark/`: reproducible benchmark harness comparing
  `qsp-codegen` against MATLAB SimBiology over a 25-compartment SBML model.
  Reports integration-only and wall-clock-per-invocation timings, and gates
  speed numbers on a parity check between the two backends.
- `CONTRIBUTING.md` and `CODE_OF_CONDUCT.md`.
- `pytest` configuration to scope test collection to `tests/` (avoids
  picking up third-party tests in benchmark scratch directories).

## [0.1.0] - 2026-05-06

Initial release. SBML → C++ code generation, CVODE-backed `qsp_sim` runtime,
parity harness, and `qsp-refresh-param-xml` companion tool.

DOI: [10.5281/zenodo.20059614](https://doi.org/10.5281/zenodo.20059614)
