"""Tests for SBML <functionDefinition> inlining in the MathML->C++ converter."""
import xml.etree.ElementTree as ET

import pytest

from qsp_codegen.codegen import SBMLModel

MATH = "http://www.w3.org/1998/Math/MathML"


def _model_with_hill():
    """Bare SBMLModel with a hand-built hill(x, k) = x / (x + k) lambda."""
    m = object.__new__(SBMLModel)
    m.id_to_name = {}
    body = ET.fromstring(
        f'<apply xmlns="{MATH}"><divide/><ci>x</ci>'
        f"<apply><plus/><ci>x</ci><ci>k</ci></apply></apply>"
    )
    m.function_defs = {"hill": {"bvars": ["x", "k"], "body": body}}
    return m


def _call(model, mathml):
    return model._mathml_to_infix(ET.fromstring(mathml))


def test_function_def_is_inlined():
    m = _model_with_hill()
    out = _call(
        m,
        f'<math xmlns="{MATH}"><apply><ci>hill</ci>'
        f'<ci>S</ci><cn type="integer">2</cn></apply></math>',
    )
    # Arguments substituted in; bound variable names must not leak.
    assert "S" in out and "2.0" in out and "/" in out
    assert "x" not in out and "k" not in out


def test_function_def_arg_count_mismatch_raises():
    m = _model_with_hill()
    with pytest.raises(ValueError, match="expects 2"):
        _call(m, f'<math xmlns="{MATH}"><apply><ci>hill</ci><ci>S</ci></apply></math>')


def test_nested_function_call_inlines():
    # hill(hill(S, 1), 2) — exercises reentrant bvar save/restore.
    m = _model_with_hill()
    out = _call(
        m,
        f'<math xmlns="{MATH}"><apply><ci>hill</ci>'
        f'<apply><ci>hill</ci><ci>S</ci><cn type="integer">1</cn></apply>'
        f'<cn type="integer">2</cn></apply></math>',
    )
    assert "S" in out and "1.0" in out and "2.0" in out
    assert "x" not in out and "k" not in out
