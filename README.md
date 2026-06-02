# qsp-codegen

SBML → C++ CVODE ODE code generator for QSP models.

Given a SimBiology-exported SBML Level 2 v4 file, emits the C++ sources
consumed by a CVODEBase-backed QSP simulator:

- `QSP_enum.h`, `ODE_system.h`, `ODE_system.cpp`
- `QSPParam.h`, `QSPParam.cpp`
- `param_all.xml` — a complete, ready-to-run parameter file (initial
  conditions + parameters taken straight from the SBML). Run the generated
  `qsp_sim` against it directly — no hand-editing, no merge step.
- `qsp_params_xml_snippet.xml` — the bare `<QSP>` block, for consumers that
  merge it into a maintained `param_all.xml` via `qsp-refresh-param-xml`.

## Prerequisites

- **Python ≥ 3.10** with `sympy` (for the analytical Jacobian) and `numpy`.
- **A C++17 compiler + CMake ≥ 3.18** to build the generated simulator. The
  first configure fetches and compiles SUNDIALS and yaml-cpp via CMake
  `FetchContent` (one-time, a few minutes).
- **MATLAB + SimBiology** only if you want to (a) export SBML from a
  SimBiology model or (b) run `qsp-codegen verify` against a SimBiology
  reference. Pure SBML-in → C++-out needs neither.

## Install

```bash
pip install ~/Projects/qsp-codegen        # or: uv pip install -e ~/Projects/qsp-codegen
```

## Quick start (SBML in → trajectories out)

```bash
# 1. Generate C++ + a ready-to-run param file from a SimBiology SBML export.
qsp-codegen generate --sbml MyModel.sbml --out-dir build/qsp/ode

# 2. Build a minimal qsp_sim against the generated ODE (driver + ODE, default
#    no-op init hook). See "Bundled C++ driver" below for the CMakeLists; or
#    let `qsp-codegen verify` scaffold and build it for you.
cmake -S build/sim -B build/sim/build -DCMAKE_BUILD_TYPE=Release
cmake --build build/sim/build --target qsp_sim -j4

# 3. Run it directly against the emitted param_all.xml — no hand-editing.
build/sim/build/qsp_sim build/qsp/ode/param_all.xml out.csv 365 4.0
```

The legacy form `qsp-codegen --sbml … --out-dir …` (no `generate` word) still
works. Re-run codegen whenever the model *structure* changes
(species/parameters/reactions); not needed for parameter-value tweaks.

## Self-test: does the C++ match SimBiology?

`qsp-codegen verify` runs the whole pipeline end-to-end — codegen, scaffold +
build a minimal `qsp_sim`, then compare its trajectories against a MATLAB
SimBiology reference over the same window — and reports PASS/FAIL:

```bash
qsp-codegen verify \
    --sbml MyModel.sbml \
    --matlab-dir /path/to/model/repo \      # has startup.m + the model script
    --matlab-script build_my_model \        # builds `model`; must NOT `clear`
    --stop-time 365
```

## Supported SBML features

Parses **SBML Level 2 v4** (the SimBiology export dialect):

- **Rate-law math**: the full SBML L2 operator set, including SimBiology's
  `<ci>`-exported named operators (`max`, `min`, `nthroot`, …), `<root>` with
  `<degree>`, and the trig / hyperbolic / log families.
- **User functions**: SBML `<functionDefinition>` lambdas are inlined at their
  call sites.
- **Rules**: assignment (`repeatedAssignment`) and initial assignments,
  including dynamic compartment volumes (e.g. a growing tumor compartment).
- **Events**: single-comparison triggers (`lt`/`leq`/`gt`/`geq`) with event
  assignments, mapped to CVODE root functions. *Not yet supported* (fail
  loudly with a clear message): event `<delay>`s and compound (`and`/`or`)
  triggers.
- **Jacobian**: analytical via sympy + CSE when sympy is present; numerical
  fallback otherwise.
- **Units**: converts `<listOfUnitDefinitions>` to SI for integration (see
  *Time units* below).

Anything outside this surface raises a **located, human-readable error** at
codegen time (naming the operator/function/event) rather than emitting broken
C++ — and the generator self-checks its output for undefined temporaries
before writing. Consumer-side invariants (sync checks, ABM-specific param
codegen) are out of scope — they live in the consumer repo.

## Bundled C++ driver (`qsp_sim_core`)

The wheel also ships the model-agnostic pieces of the CVODE-backed
standalone simulator under `qsp_codegen/cpp/`:

- `include/qsp_sim_core/CVODEBase.h`, `ParamBase.h`, `MolecularModelCVode.h`
- `include/qsp_sim_core/model_hooks.h` — declares the `evolve_to_diagnosis`
  override point used by the driver to set up a model-specific initial
  state before the scenario sim.
- `src/CVODEBase.cpp`, `ParamBase.cpp`, `default_hooks.cpp`
- `src/qsp_sim_main.cpp` — the `qsp_sim` driver (CVODE stepping, YAML
  scenario + drug-metadata parsing, segmented dose scheduling, CSV and
  binary trajectory output, QSTH evolve-cache I/O).
- `CMakeLists.txt`, `cmake/qsp_sim_coreConfig.cmake`, `cmake/QspSimCoreDeps.cmake`.

Consumer `CMakeLists.txt` pulls this in via Python:

```cmake
execute_process(COMMAND python -m qsp_codegen.cmake --prefix
                OUTPUT_VARIABLE QSP_SIM_CORE_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)
list(APPEND CMAKE_PREFIX_PATH "${QSP_SIM_CORE_PREFIX}")
find_package(qsp_sim_core CONFIG REQUIRED)

add_executable(qsp_sim
    ${QSP_SIM_CORE_DRIVER_SOURCE}
    qsp/ode/ODE_system.cpp                 # emitted by qsp-codegen
    qsp/ode/QSPParam.cpp                   # emitted by qsp-codegen
    sim/evolve_to_diagnosis.cpp            # consumer override (optional)
    sim/set_healthy_populations.cpp        # consumer-specific init
)
target_include_directories(qsp_sim PRIVATE qsp/ode sim)
target_link_libraries(qsp_sim PRIVATE qsp_sim_core::qsp_sim_core)
```

### Time units

`qsp_sim` integrates in SI seconds by default (multiplies external
`--t-end-days` and dose times by 86400 internally) on the assumption
that the SBML's `<listOfUnitDefinitions>` declares rates in `1/day`,
`1/hr`, etc. and the codegen has unit-converted them to `1/s`. This
matches QSP models authored with proper unit annotations.

Models exported by SimBiology's vanilla `sbmlexport` without unit
definitions come through with rate constants left in their authoring
unit (typically `1/day`). For these, pass `--time-unit days` at the
command line so the runtime integrates in model-native days; the same
binary handles both conventions:

```bash
# Unit-annotated SBML (default; production QSP models):
qsp_sim --param param_all.xml --csv-out out.csv --t-end-days 365

# SimBiology-vanilla export with unitless 1/day rates:
qsp_sim --param param_all.xml --csv-out out.csv --t-end-days 365 \
        --time-unit days
```

A consumer build can flip the *default* by compiling `qsp_sim` with
`-DMODEL_UNITS` (sets the default time factor to 1.0); the runtime
flag still wins on a per-invocation basis.

### Model-init hook (`evolve_to_diagnosis`)

Declared in `qsp_sim_core/model_hooks.h`:

```cpp
namespace CancerVCT {
struct EvolveOpts   { std::string yaml_path; double time_factor; bool verbose; /* ... */ };
struct EvolveResult { bool success; double t_diagnosis_days; double diameter_cm; std::string reject_reason; };
EvolveResult evolve_to_diagnosis(ODE_system& ode, const EvolveOpts& opts);
}
```

The driver calls this when `--evolve-to-diagnosis <yaml>` is passed on
the command line. A default no-op implementation ships inside the
`qsp_sim_core` static library and returns `{success=true, t_diag=0}`,
so models without a pre-scenario evolve phase don't need to do
anything. Models that do (e.g. evolve from healthy tissue until a
tumor diameter is reached) define their own `evolve_to_diagnosis`
in an object file compiled directly into the executable — the strong
definition wins over the archive member at link time.

The hook is a free function rather than a virtual on `ODE_system` so
that qsp-codegen's emitted ABI stays untouched as the hook interface
evolves.
