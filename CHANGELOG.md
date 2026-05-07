# Changelog

All notable changes to this project will be documented in this file. The
format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

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
