# Contributing to qsp-codegen

Thanks for your interest. This is a focused tool — it generates a C++
CVODE simulator from a SimBiology-exported SBML model — so the bar for
adding new behavior is "does this make the generated code more correct
or more useful for QSP-on-HPC workflows?" rather than "does this add a
feature."

## Reporting issues

Open an issue at <https://github.com/popellab/qsp-codegen/issues> and
include:

- The SBML file (or a minimal reduction of it) that triggers the issue.
- The exact command line you ran and the full output.
- Your `qsp-codegen --version`, Python version, and OS.
- For numerical disagreement against MATLAB: a parity report from
  `qsp_codegen.parity.compare`, ideally for a small parameter set so we
  can reproduce.

## Asking questions / getting help

Preference order: GitHub issue (so the answer is searchable) → email
the authors listed in `CITATION.cff`.

## Submitting changes

1. Fork and create a feature branch off `main`.
2. Run the test suite locally before pushing:

   ```bash
   uv pip install -e ".[dev]"
   pytest tests/
   ```

3. If your change touches code generation, also re-run the parity
   harness against a representative model. The reproducible benchmark
   in `paper/benchmark/` is a good smoke test if you do not have a
   private model to point at.
4. Open a pull request describing what changed and why. Reference any
   issue it closes.

## Coding conventions

- Python: standard `pep8`, type hints encouraged on public functions.
- C++: C++17, follow the surrounding style in `cpp/`.
- Keep comments focused on *why* something is done a particular way
  rather than *what* the code does.

## Scope

Changes that fall *outside* the project's scope and are better filed in
a consumer repository:

- Model-specific RHS overrides or hooks (those belong in the consumer's
  C++ driver, plugged in via `model_hooks.h`).
- Simulator-orchestration logic (parameter sweeps, scheduling, caching
  policy) — those live in `qsp-hpc-tools`.
- Per-model parameter XML (those live in the consumer's resources
  directory and are kept in sync via `qsp-refresh-param-xml`).

If unsure, open an issue first to discuss the right home for a change.
