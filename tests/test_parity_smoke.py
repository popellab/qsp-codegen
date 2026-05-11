"""Smoke tests for the parity harness module.

Doesn't invoke MATLAB or a compiled binary — those belong in consumer-repo
integration tests. Here we just check that the module imports, the MATLAB
assets ship with the package, and ``compare`` handles trivial cases
correctly.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from qsp_codegen import parity
from qsp_codegen import sync_checks


def test_matlab_asset_dir_contains_scripts():
    assets = parity.matlab_asset_dir()
    assert (assets / "export_matlab_trajectories.m").exists()
    assert (assets / "yaml_read.m").exists()


def test_compare_identical_trajectories(tmp_path: Path):
    csv = tmp_path / "traj.csv"
    csv.write_text("Time,V_T.CD8,V_T.Treg\n0,1,2\n1,3,4\n2,5,6\n")
    passed, report = parity.compare(str(csv), str(csv), rtol=1e-6, atol=1e-9)
    assert passed, report


def test_compare_detects_divergence(tmp_path: Path):
    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    a.write_text("Time,V_T.CD8\n0,1\n1,2\n")
    b.write_text("Time,V_T.CD8\n0,1\n1,10\n")
    passed, report = parity.compare(str(a), str(b), rtol=0.05, atol=1e-9)
    assert not passed
    assert "V_T.CD8" in report


def test_compare_irregular_cpp_grid_interpolates(tmp_path: Path):
    """Under qsp-codegen v3 (CV_ONE_STEP), the C++ side emits a non-uniform
    time grid that doesn't line up with MATLAB's fixed dt. compare() should
    linear-interpolate C++ onto MATLAB's grid and pass when the underlying
    species function agrees."""
    matlab_csv = tmp_path / "matlab.csv"
    cpp_csv = tmp_path / "cpp.csv"
    # MATLAB on a uniform 0.1 d grid; species value = 2*t exactly.
    matlab_csv.write_text(
        "Time,V_T.CD8\n"
        + "\n".join(f"{t:.4f},{2 * t:.6e}" for t in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        + "\n"
    )
    # C++ on an irregular grid that hits none of MATLAB's tick times; the
    # exact linear shape means np.interp recovers the species value at any
    # interior MATLAB tick to floating-point roundoff.
    cpp_csv.write_text(
        "Time,V_T.CD8\n"
        + "\n".join(f"{t:.4f},{2 * t:.6e}" for t in [0.0, 0.07, 0.23, 0.41, 0.5])
        + "\n"
    )
    passed, report = parity.compare(
        str(matlab_csv), str(cpp_csv), rtol=1e-6, atol=1e-12
    )
    assert passed, report
    assert "linear-interp" in report

    # Sanity: divergence on the same irregular C++ grid still trips the
    # comparator (scale species column by 2 → 100% rdiff at every point).
    cpp_csv.write_text(
        "Time,V_T.CD8\n"
        + "\n".join(f"{t:.4f},{4 * t:.6e}" for t in [0.0, 0.07, 0.23, 0.41, 0.5])
        + "\n"
    )
    passed, report = parity.compare(
        str(matlab_csv), str(cpp_csv), rtol=0.05, atol=1e-9
    )
    assert not passed


def test_run_cpp_trajectories_rejects_missing_binary(tmp_path: Path):
    with pytest.raises(FileNotFoundError, match="qsp_sim binary not found"):
        parity.run_cpp_trajectories(
            qsp_sim=tmp_path / "nope",
            param_xml=tmp_path / "nope.xml",
            out_csv=tmp_path / "out.csv",
            t_end_days=1, min_cadence_hours=4.0,
        )


def test_run_cpp_trajectories_requires_drug_metadata_with_scenario(tmp_path: Path):
    """Passing scenario_yaml without drug_metadata_yaml should be a hard fail
    (the binary needs the MW/dose-basis table to interpret doses)."""
    # fake the binary + xml so we reach the kwargs validation
    fake_bin = tmp_path / "qsp_sim"
    fake_bin.write_text("")
    fake_bin.chmod(0o755)
    fake_xml = tmp_path / "param.xml"
    fake_xml.write_text("<root/>")
    scenario = tmp_path / "s.yaml"
    scenario.write_text("stop_time: 1\n")
    with pytest.raises(ValueError, match="drug_metadata_yaml is required"):
        parity.run_cpp_trajectories(
            qsp_sim=fake_bin,
            param_xml=fake_xml,
            out_csv=tmp_path / "out.csv",
            t_end_days=1, min_cadence_hours=4.0,
            scenario_yaml=scenario,
        )


def test_run_cpp_trajectories_validates_kwarg_paths(tmp_path: Path):
    """Missing scenario / drug-metadata / evolve YAML surfaces a clear error."""
    fake_bin = tmp_path / "qsp_sim"
    fake_bin.write_text("")
    fake_bin.chmod(0o755)
    fake_xml = tmp_path / "param.xml"
    fake_xml.write_text("<root/>")
    with pytest.raises(FileNotFoundError, match="evolve_to_diagnosis_yaml not found"):
        parity.run_cpp_trajectories(
            qsp_sim=fake_bin,
            param_xml=fake_xml,
            out_csv=tmp_path / "out.csv",
            t_end_days=1, min_cadence_hours=4.0,
            evolve_to_diagnosis_yaml=tmp_path / "missing_healthy.yaml",
        )


def test_sync_check_sbml_newer_than_matlab_handles_missing(tmp_path: Path):
    ok, msg = sync_checks.check_sbml_newer_than_matlab(
        sbml=tmp_path / "does_not_exist.sbml",
        matlab_script=tmp_path / "also_missing.m",
    )
    # matlab script missing → skip (not a failure)
    assert ok and "skip" in msg


def test_sync_check_priors_csv_missing_names(tmp_path: Path):
    priors = tmp_path / "priors.csv"
    xml = tmp_path / "params.xml"
    priors.write_text("name,value\nfoo,1.0\nbar,2.0\nbaz,3.0\n")
    xml.write_text("<root><foo>1</foo><bar>2</bar></root>")
    ok, msg = sync_checks.check_priors_csv_names_in_param_xml(
        priors_csv=priors, param_xml=xml,
    )
    assert not ok
    assert "baz" in msg
