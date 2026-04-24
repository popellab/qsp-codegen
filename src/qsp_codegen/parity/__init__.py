"""MATLAB ↔ C++ trajectory parity harness for QSP models.

Consumer repos drive their own ``qsp_sim`` binary and MATLAB model script
through the same harness so fixtures stay in-tree and reproducible. The
typical shape:

.. code-block:: python

    from qsp_codegen.parity import run_cpp_trajectories, run_matlab_trajectories, compare

    # Default-IC 365-day parity (no dosing, no natural-history evolution):
    cpp_csv = run_cpp_trajectories(
        qsp_sim=Path("cpp/sim/build/qsp_sim"),
        param_xml=Path("resources/cpp/param_all.xml"),
        out_csv=tmp_path / "cpp.csv",
        t_end_days=365,
        dt_days=0.1,
    )

    # Scenario parity with dosing + evolve_to_diagnosis on both sides:
    cpp_csv = run_cpp_trajectories(
        qsp_sim=Path("cpp/sim/build/qsp_sim"),
        param_xml=Path("resources/cpp/param_all.xml"),
        out_csv=tmp_path / "cpp.csv",
        t_end_days=30,
        dt_days=0.1,
        scenario_yaml=Path("scenarios/gvax_nivo_neoadjuvant.yaml"),
        drug_metadata_yaml=Path("resources/cpp/drug_metadata.yaml"),
        evolve_to_diagnosis_yaml=Path("resources/cpp/healthy_state.yaml"),
    )
    matlab_csv = run_matlab_trajectories(
        matlab_model_dir=Path("."),
        matlab_model_script="immune_oncology_model_PDAC",
        sbml_path=Path("PDAC_model.sbml"),
        param_xml=Path("resources/cpp/param_all.xml"),
        out_csv=tmp_path / "matlab.csv",
        stop_time=30,
        scenario_yaml=Path("scenarios/gvax_nivo_neoadjuvant.yaml"),
        evolve_function_name="evolve_to_diagnosis",
    )
    passed, report = compare(matlab_csv, cpp_csv, rtol=0.05)

``matlab_asset_dir()`` returns the package-internal path to the MATLAB
scripts (``export_matlab_trajectories.m``, ``yaml_read.m``) so consumers
don't have to vendor copies.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

from .compare_trajectories import compare, load_csv  # re-export

__all__ = [
    "compare",
    "load_csv",
    "matlab_asset_dir",
    "run_cpp_trajectories",
    "run_matlab_trajectories",
]


def matlab_asset_dir() -> Path:
    """Absolute path to the MATLAB harness scripts shipped with the package."""
    return Path(__file__).resolve().parent


def run_cpp_trajectories(
    qsp_sim: Path,
    param_xml: Path,
    out_csv: Path,
    t_end_days: float,
    dt_days: float,
    scenario_yaml: Optional[Path] = None,
    drug_metadata_yaml: Optional[Path] = None,
    evolve_to_diagnosis_yaml: Optional[Path] = None,
    timeout: float = 600.0,
) -> Path:
    """Run the compiled ``qsp_sim`` binary and return the output CSV path.

    Positional args always sent: ``<param_xml> <csv_out> <t_end> <dt>``.
    The three optional kwargs each map to a long flag on ``qsp_sim``:

    - ``scenario_yaml`` → ``--scenario``: apply the dosing schedule /
      config described by the YAML (typically built from a pdac-build-
      style ``scenarios/*.yaml``).
    - ``drug_metadata_yaml`` → ``--drug-metadata``: molecular-weight and
      dose-basis table the scenario references. Required by the binary
      whenever ``--scenario`` is used.
    - ``evolve_to_diagnosis_yaml`` → ``--evolve-to-diagnosis``: healthy-
      state YAML to evolve ICs from before the scenario window starts,
      matching the MATLAB side's ``evolve_to_diagnosis.m`` flow.

    Raises ``RuntimeError`` on nonzero exit or if the binary silently
    defaulted params to zero (stderr "QSP param not found"). Returns
    ``out_csv`` for convenience.
    """
    if not qsp_sim.exists():
        raise FileNotFoundError(f"qsp_sim binary not found: {qsp_sim}")
    if not param_xml.exists():
        raise FileNotFoundError(f"param XML not found: {param_xml}")
    if scenario_yaml is not None and drug_metadata_yaml is None:
        raise ValueError(
            "drug_metadata_yaml is required when scenario_yaml is provided — "
            "qsp_sim needs the drug MW/dose-basis table to interpret doses."
        )
    for kwarg_name, p in (
        ("scenario_yaml", scenario_yaml),
        ("drug_metadata_yaml", drug_metadata_yaml),
        ("evolve_to_diagnosis_yaml", evolve_to_diagnosis_yaml),
    ):
        if p is not None and not p.exists():
            raise FileNotFoundError(f"{kwarg_name} not found: {p}")

    argv = [str(qsp_sim), str(param_xml), str(out_csv),
            str(t_end_days), str(dt_days)]
    if scenario_yaml is not None:
        argv.extend(["--scenario", str(scenario_yaml)])
    if drug_metadata_yaml is not None:
        argv.extend(["--drug-metadata", str(drug_metadata_yaml)])
    if evolve_to_diagnosis_yaml is not None:
        argv.extend(["--evolve-to-diagnosis", str(evolve_to_diagnosis_yaml)])

    result = subprocess.run(
        argv, capture_output=True, text=True, timeout=timeout,
    )
    if result.returncode != 0:
        raise RuntimeError(f"qsp_sim failed (exit {result.returncode}):\n{result.stderr}")
    missing = [line for line in result.stderr.splitlines()
               if "QSP param not found" in line]
    if missing:
        raise RuntimeError(
            f"{len(missing)} QSP params silently defaulted to 0 — refresh "
            f"param_all.xml from the generated snippet. First few:\n  "
            + "\n  ".join(missing[:5])
        )
    if not out_csv.exists():
        raise RuntimeError(f"qsp_sim returned 0 but wrote no CSV: {out_csv}")
    return out_csv


def run_matlab_trajectories(
    matlab_model_dir: Path,
    matlab_model_script: str,
    sbml_path: Path,
    param_xml: Path,
    out_csv: Path,
    stop_time: float = 365.0,
    scenario_yaml: Optional[Path] = None,
    evolve_function_name: Optional[str] = None,
    matlab_binary: str = "matlab",
    timeout: float = 600.0,
) -> Path:
    """Run ``export_matlab_trajectories.m`` under MATLAB with the given inputs.

    Parameters mirror the script's expected workspace variables:

    - ``matlab_model_dir``: consumer repo root (contains ``startup.m``).
    - ``matlab_model_script``: bare script name (no ``.m``) that builds
      ``model`` in the MATLAB workspace (e.g., the live ``*_PDAC.m`` file).
    - ``sbml_path``: the SBML exported from that script (used for parity
      provenance only — MATLAB loads the live .m, not the SBML).
    - ``param_xml``: the same XML template the C++ side consumes.
    - ``evolve_function_name``: if provided, invoked on ``model`` before
      running the scenario (e.g., "evolve_to_diagnosis"). Must be on the
      MATLAB path.
    """
    asset_dir = matlab_asset_dir()
    cmd_parts = [
        f"output_csv='{out_csv}';",
        f"matlab_model_dir='{matlab_model_dir}';",
        f"matlab_model_script='{matlab_model_script}';",
        f"param_xml='{param_xml}';",
        f"sbml_path='{sbml_path}';",
        f"stop_time={stop_time};",
        f"addpath('{asset_dir}');",
    ]
    if scenario_yaml is not None:
        cmd_parts.append(f"scenario_yaml='{scenario_yaml}';")
    if evolve_function_name is not None:
        # MATLAB scripts use the `matlab_evolve_function` handle variable.
        cmd_parts.append(f"matlab_evolve_function=@{evolve_function_name};")
    cmd_parts.append(f"run('{asset_dir / 'export_matlab_trajectories.m'}')")
    matlab_cmd = " ".join(cmd_parts)
    result = subprocess.run(
        [matlab_binary, "-batch", matlab_cmd],
        capture_output=True, text=True, timeout=timeout,
        cwd=str(matlab_model_dir),
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"MATLAB export failed (exit {result.returncode}):\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    if not out_csv.exists():
        raise RuntimeError(
            f"MATLAB returned 0 but wrote no CSV: {out_csv}\n"
            f"STDOUT:\n{result.stdout}"
        )
    return out_csv
