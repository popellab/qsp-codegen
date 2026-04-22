# qsp-codegen

SBML → C++ CVODE ODE code generator for QSP models.

Given a SimBiology-exported SBML Level 2 v4 file, emits the C++ sources
consumed by a CVODEBase-backed QSP simulator:

- `QSP_enum.h`, `ODE_system.h`, `ODE_system.cpp`
- `QSPParam.h`, `QSPParam.cpp`
- `qsp_params_xml_snippet.xml` (merged into consumer `param_all.xml` via
  `qsp-refresh-param-xml`)

## Install

```bash
pip install ~/Projects/qsp-codegen
```

## Usage

```bash
qsp-codegen --sbml path/to/PDAC_model.sbml --out-dir path/to/ode/

qsp-refresh-param-xml \
    --snippet path/to/ode/qsp_params_xml_snippet.xml \
    --xml path/to/param_all.xml \
    --xml path/to/param_all_test.xml
```

Run whenever the QSP model structure changes (species/parameters/reactions).
Not needed for parameter-value tweaks.

## Scope

- Parses SBML Level 2 v4 (SimBiology export dialect).
- Derives analytical Jacobian via sympy + CSE when sympy is available;
  falls back to numerical Jacobian otherwise.
- Consumer-side invariants (sync checks, ABM-specific param codegen) are
  *not* in scope — they live in the consumer repo.
