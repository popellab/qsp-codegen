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
    "param_all.xml",  # complete, ready-to-run param file (wraps the snippet)
}


@pytest.mark.skipif(not PDAC_SBML.exists(), reason="PDAC SBML not found locally")
def test_generate_pdac(tmp_path):
    files = generate(str(PDAC_SBML), str(tmp_path))
    assert set(files.keys()) == EXPECTED
    for name in EXPECTED:
        p = tmp_path / name
        assert p.exists() and p.stat().st_size > 0


def test_cli_generate_requires_sbml_and_out():
    from qsp_codegen.codegen import main

    # The `generate` subcommand still enforces its required args.
    with pytest.raises(SystemExit):
        main(["generate"])


def test_cli_no_subcommand_prints_help_and_returns_nonzero():
    from qsp_codegen.codegen import main

    assert main([]) == 1


def test_cli_legacy_bare_options_still_route_to_generate(tmp_path):
    # Back-compat: `qsp-codegen --sbml ... --out-dir ...` (no subcommand word)
    # must still reach `generate`. With a non-existent SBML it should fail-fast,
    # not silently no-op — assert it raises rather than returning 0.
    from qsp_codegen.codegen import main

    with pytest.raises(Exception):
        main(["--sbml", str(tmp_path / "nope.sbml"), "--out-dir", str(tmp_path)])
