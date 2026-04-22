"""Smoke test: codegen against the live PDAC SBML (if present) should emit
all expected output files."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from qsp_codegen.codegen import generate

PDAC_SBML = Path.home() / "Projects" / "pdac-build" / "PDAC_model.sbml"

EXPECTED = {
    "QSP_enum.h",
    "ODE_system.h",
    "ODE_system.cpp",
    "QSPParam.h",
    "QSPParam.cpp",
    "qsp_params_xml_snippet.xml",
}


@pytest.mark.skipif(not PDAC_SBML.exists(), reason="PDAC SBML not found locally")
def test_generate_pdac(tmp_path):
    files = generate(str(PDAC_SBML), str(tmp_path))
    assert set(files.keys()) == EXPECTED
    for name in EXPECTED:
        p = tmp_path / name
        assert p.exists() and p.stat().st_size > 0


def test_cli_requires_sbml_and_out():
    from qsp_codegen.codegen import main

    with pytest.raises(SystemExit):
        main([])
