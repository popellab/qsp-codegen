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
