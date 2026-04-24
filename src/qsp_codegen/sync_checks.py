"""Freshness + consistency checks across the SBML → codegen → binary → XML chain.

Consumers (pdac-build, other QSP model repos) pass their own paths in. Every
check returns ``(ok, message)``; messages are actionable so failures tell the
user exactly which command to run.

The checks:

* ``check_sbml_newer_than_matlab`` — SBML must be re-exported after edits to
  the live MATLAB model script.
* ``check_codegen_newer_than_sbml`` — ``ODE_system.cpp`` must be regenerated
  after the SBML changes.
* ``check_binary_newer_than_codegen`` — the built binary must be newer than
  its generated source.
* ``check_param_xml_contains_snippet`` — every leaf tag in the generated
  ``qsp_params_xml_snippet.xml`` must appear in the consumer's ``param_all.xml``.
* ``check_priors_csv_names_in_param_xml`` — every name in the consumer's
  priors CSV must appear in ``param_all.xml`` (so the sampled XML template
  render doesn't abort with ParamNotFoundError at simulation time).
"""
from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Tuple

Result = Tuple[bool, str]


def check_sbml_newer_than_matlab(
    sbml: Path,
    matlab_script: Path,
    regen_cmd: str = "matlab -batch \"run('scripts/export_sbml.m')\"",
) -> Result:
    """SBML must be re-exported after edits to the live MATLAB model script."""
    if not matlab_script.exists():
        return True, f"skip: MATLAB model script not found ({matlab_script})"
    if not sbml.exists():
        return False, f"SBML missing: {sbml}\n  Re-export: {regen_cmd}"
    drift = matlab_script.stat().st_mtime - sbml.stat().st_mtime
    if drift > 0:
        return False, (
            f"{sbml.name} is {drift:.0f}s older than {matlab_script.name}.\n"
            f"  Re-export: {regen_cmd}"
        )
    return True, "SBML is up to date with MATLAB model script"


def check_codegen_newer_than_sbml(
    ode_cpp: Path,
    sbml: Path,
    regen_cmd: str = "qsp-codegen --sbml <path> --out-dir <path>",
) -> Result:
    """``ODE_system.cpp`` must be regenerated after the SBML changes."""
    if not sbml.exists():
        return True, f"skip: SBML not found ({sbml})"
    if not ode_cpp.exists():
        return False, f"Generated ODE missing: {ode_cpp}\n  Run: {regen_cmd}"
    drift = sbml.stat().st_mtime - ode_cpp.stat().st_mtime
    if drift > 0:
        return False, (
            f"{ode_cpp.name} is {drift:.0f}s older than the SBML.\n"
            f"  Regenerate: {regen_cmd}"
        )
    return True, "Generated ODE is up to date with SBML"


def check_binary_newer_than_codegen(
    dump_bin: Path,
    ode_cpp: Path,
) -> Result:
    """The built binary must be newer than its generated source."""
    if not dump_bin.exists():
        return True, f"skip: binary not built ({dump_bin})"
    if not ode_cpp.exists():
        return True, "skip: no codegen output to compare against"
    drift = ode_cpp.stat().st_mtime - dump_bin.stat().st_mtime
    if drift > 0:
        return False, (
            f"{dump_bin.name} is {drift:.0f}s older than {ode_cpp.name}.\n"
            f"  Rebuild: cmake --build {dump_bin.parent} --target {dump_bin.name}"
        )
    return True, f"{dump_bin.name} is up to date"


def check_param_xml_contains_snippet(
    snippet: Path,
    param_xml: Path,
    regen_cmd: str = "qsp-refresh-param-xml --snippet <snippet> --xml <param_xml>",
) -> Result:
    """Every leaf tag in the generated snippet must appear in ``param_all.xml``."""
    if not snippet.exists():
        return True, f"skip: snippet not found ({snippet})"
    if not param_xml.exists():
        return False, f"param_all.xml missing: {param_xml}"
    snippet_names = set(re.findall(r"<([A-Za-z_][\w]*)>", snippet.read_text()))
    xml_names = set(re.findall(r"<([A-Za-z_][\w]*)>", param_xml.read_text()))
    missing = sorted(snippet_names - xml_names)
    if missing:
        return False, (
            f"{len(missing)} name(s) in {snippet.name} missing from "
            f"{param_xml.name}:\n  {missing[:20]}\n  Refresh: {regen_cmd}"
        )
    return True, f"{param_xml.name} contains all snippet names"


def check_priors_csv_names_in_param_xml(
    priors_csv: Path,
    param_xml: Path,
    name_column: str = "name",
) -> Result:
    """Every parameter name in the priors CSV must appear in ``param_all.xml``.

    The C++ worker substitutes each prior value into the XML template; a row
    whose name is absent from the template aborts the sim with
    ``ParamNotFoundError`` before the ODE solver runs. This check catches
    orphan prior rows for model components that were removed or renamed.
    """
    if not priors_csv.exists():
        return True, f"skip: priors CSV not found ({priors_csv})"
    if not param_xml.exists():
        return False, f"param_all.xml missing: {param_xml}"
    with open(priors_csv) as f:
        prior_names = {row[name_column] for row in csv.DictReader(f)}
    xml_names = set(re.findall(r"<([A-Za-z_][\w]*)>", param_xml.read_text()))
    missing = sorted(prior_names - xml_names)
    if missing:
        return False, (
            f"{len(missing)} prior(s) in {priors_csv.name} missing from "
            f"{param_xml.name}:\n  {missing[:20]}"
            f"{'...' if len(missing) > 20 else ''}\n"
            f"  Fix: drop the orphan rows, or re-export SBML + regenerate."
        )
    return True, f"all {len(prior_names)} prior names present in {param_xml.name}"


def run_all(
    checks: Iterable[Callable[[], Result]],
) -> List[Tuple[str, bool, str]]:
    """Run every zero-arg check callable and return ``(name, ok, message)``.

    Callers typically pass a list of ``functools.partial(check_fn, path=...)``
    bound to their repo's paths. Returns one triple per check so a CLI or test
    wrapper can decide how to format / abort.
    """
    return [(fn.__name__ if hasattr(fn, "__name__") else "check",
             *fn()) for fn in checks]
