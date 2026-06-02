"""``qsp-codegen verify`` — end-to-end self-test for a generated model.

Answers the first question any new user has: *does the generated C++ actually
reproduce my SimBiology model?* In one command it:

  1. generates C++ from the SBML (reusing :func:`qsp_codegen.codegen.generate`),
  2. scaffolds a minimal ``qsp_sim`` CMake project (driver + generated ODE +
     the library's default no-op model-init hook) and builds it,
  3. assembles a ready-to-run ``param_all.xml`` from the codegen snippet,
  4. runs the compiled binary and a MATLAB SimBiology reference over the same
     window, and compares every matched species.

The build step needs a C++ toolchain + CMake (the first configure fetches and
compiles SUNDIALS/yaml-cpp); the reference step needs MATLAB on the path (or
``--matlab``). The MATLAB model script must build a ``model`` variable and must
NOT start with ``clear`` (it is ``run`` inside the parity harness's workspace).
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional

from .codegen import generate

# Minimal, model-agnostic consumer CMake project. Mirrors the thin
# pdac-build/cpp/sim consumer but with NO model-specific init hook: the
# library's default no-op evolve_to_diagnosis (default_hooks.cpp) is linked,
# which is exactly right for a default-IC parity. ``${ODE_DIR}`` is filled in.
_CMAKE_TEMPLATE = """cmake_minimum_required(VERSION 3.18)
project(qsp_sim CXX)

if(NOT DEFINED QSP_SIM_CORE_PREFIX)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -m qsp_codegen.cmake --prefix
        OUTPUT_VARIABLE QSP_SIM_CORE_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _qsp_rc)
    if(NOT _qsp_rc EQUAL 0)
        message(FATAL_ERROR "Could not locate qsp_sim_core (qsp-codegen wheel).")
    endif()
endif()
list(APPEND CMAKE_PREFIX_PATH "${QSP_SIM_CORE_PREFIX}")
find_package(qsp_sim_core CONFIG REQUIRED)

add_executable(qsp_sim
    ${QSP_SIM_CORE_DRIVER_SOURCE}
    @ODE_DIR@/ODE_system.cpp
    @ODE_DIR@/QSPParam.cpp
)
target_include_directories(qsp_sim PRIVATE @ODE_DIR@)
target_link_libraries(qsp_sim PRIVATE qsp_sim_core::qsp_sim_core)
"""


def assemble_param_xml(snippet_path: Path, out_path: Path) -> Path:
    """Wrap a codegen ``<QSP>`` snippet into a complete ``param_all.xml``.

    ``generate()`` now emits ``param_all.xml`` directly; this remains for
    callers that only have the snippet on hand.
    """
    from .codegen import wrap_param_xml

    out_path.write_text(wrap_param_xml(snippet_path.read_text()))
    return out_path


def _build_qsp_sim(ode_dir: Path, work_dir: Path, python_exe: str) -> Path:
    """Scaffold + build a minimal qsp_sim against the generated ODE."""
    sim_dir = work_dir / "sim"
    sim_dir.mkdir(parents=True, exist_ok=True)
    cmake_txt = _CMAKE_TEMPLATE.replace("@ODE_DIR@", str(ode_dir.resolve()))
    (sim_dir / "CMakeLists.txt").write_text(cmake_txt)
    build_dir = sim_dir / "build"

    print(f"  configuring + building qsp_sim in {build_dir} ...")
    subprocess.run(
        ["cmake", "-S", str(sim_dir), "-B", str(build_dir),
         "-DCMAKE_BUILD_TYPE=Release", f"-DPython3_EXECUTABLE={python_exe}"],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        ["cmake", "--build", str(build_dir), "--target", "qsp_sim", "-j", "4"],
        check=True, capture_output=True, text=True,
    )
    binary = build_dir / "qsp_sim"
    if not binary.exists():
        raise RuntimeError(f"build reported success but no binary at {binary}")
    return binary


def run_verify(
    sbml: Path,
    matlab_dir: Path,
    matlab_script: str,
    work_dir: Path,
    stop_time: float = 365.0,
    rtol: float = 0.05,
    atol: float = 1e-6,
    matlab_binary: str = "matlab",
    python_exe: Optional[str] = None,
) -> bool:
    """Run the full codegen → build → parity self-test. Returns pass/fail."""
    # Imported lazily so `generate`/validation work without numpy/matlab present.
    from .parity import compare, run_cpp_trajectories, run_matlab_trajectories

    python_exe = python_exe or sys.executable
    work_dir.mkdir(parents=True, exist_ok=True)
    ode_dir = work_dir / "qsp" / "ode"

    print(f"[1/4] codegen {sbml} -> {ode_dir}")
    generate(str(sbml), str(ode_dir))

    print("[2/4] build qsp_sim")
    qsp_sim = _build_qsp_sim(ode_dir, work_dir, python_exe)

    param_xml = assemble_param_xml(
        ode_dir / "qsp_params_xml_snippet.xml", work_dir / "param_all.xml"
    )

    cpp_csv = work_dir / "cpp.csv"
    matlab_csv = work_dir / "matlab.csv"
    print(f"[3/4] run C++ + MATLAB ({stop_time}d, default ICs, grid-pinned)")
    run_cpp_trajectories(
        qsp_sim=qsp_sim, param_xml=param_xml, out_csv=cpp_csv,
        t_end_days=stop_time, min_cadence_hours=4.0,
    )
    # Pin MATLAB to the C++ grid so the compare is row-aligned (no interp
    # artifacts on stiff early transients).
    import numpy as np
    cpp_times = np.loadtxt(cpp_csv, delimiter=",", skiprows=1, usecols=0)
    times_csv = work_dir / "cpp_times.csv"
    np.savetxt(times_csv, cpp_times)
    run_matlab_trajectories(
        matlab_model_dir=matlab_dir, matlab_model_script=matlab_script,
        sbml_path=sbml, param_xml=param_xml, out_csv=matlab_csv,
        stop_time=stop_time, output_times_csv=times_csv, matlab_binary=matlab_binary,
    )

    print(f"[4/4] compare (rtol={rtol}, atol={atol})")
    passed, report = compare(str(matlab_csv), str(cpp_csv), rtol=rtol, atol=atol)
    print(report)
    print("\nVERIFY:", "PASS ✅" if passed else "FAIL ❌")
    return passed


def add_subparser(subparsers) -> None:
    """Register the ``verify`` subcommand on a codegen argparse subparsers."""
    p = subparsers.add_parser(
        "verify",
        help="Codegen + build + C++↔MATLAB parity self-test for an SBML model.",
    )
    p.add_argument("--sbml", required=True, type=Path, help="SBML model file.")
    p.add_argument("--matlab-dir", required=True, type=Path,
                   help="Consumer repo root (has startup.m + the model script).")
    p.add_argument("--matlab-script", required=True,
                   help="Bare script name that builds `model` (no `.m`, no `clear`).")
    p.add_argument("--work-dir", type=Path, default=Path("qsp_verify_out"),
                   help="Scratch dir for generated code, build, and CSVs.")
    p.add_argument("--stop-time", type=float, default=365.0, help="Sim days.")
    p.add_argument("--rtol", type=float, default=0.05)
    p.add_argument("--atol", type=float, default=1e-6)
    p.add_argument("--matlab", default="matlab", help="MATLAB binary path.")
    p.set_defaults(_handler=_handle)


def _handle(args) -> int:
    ok = run_verify(
        sbml=args.sbml, matlab_dir=args.matlab_dir, matlab_script=args.matlab_script,
        work_dir=args.work_dir, stop_time=args.stop_time, rtol=args.rtol,
        atol=args.atol, matlab_binary=args.matlab,
    )
    return 0 if ok else 1
