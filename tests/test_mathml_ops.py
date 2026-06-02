"""Unit tests for the MathML -> C++ infix converter (SBMLModel._mathml_to_infix).

Guards the operator coverage that SimBiology exercises via SBML export. The
motivating regression: SimBiology exports `nthroot` (and max/min) as
``<ci>name</ci>`` rather than a native MathML tag; an unsupported name used to
fall through to a ``/* unknown op */`` comment that only blew up much later in
the Jacobian ``sympify``. The converter now dispatches the full L2 math set and
fails loudly on anything it genuinely cannot translate.
"""
import xml.etree.ElementTree as ET

import pytest

from qsp_codegen.codegen import SBMLModel

MATH = "http://www.w3.org/1998/Math/MathML"


def _convert(mathml: str) -> str:
    # _mathml_to_infix only touches self.id_to_name; build a bare instance.
    model = object.__new__(SBMLModel)
    model.id_to_name = {}
    model.function_defs = {}
    return model._mathml_to_infix(ET.fromstring(mathml))


def _apply(op_xml: str, *operand_xml: str) -> str:
    operands = "".join(operand_xml)
    return _convert(f'<math xmlns="{MATH}"><apply>{op_xml}{operands}</apply></math>')


X = "<ci>x</ci>"
N3 = '<cn type="integer">3</cn>'
N2 = '<cn type="integer">2</cn>'


def test_nthroot_ci_becomes_pow_reciprocal():
    # nthroot(x, 3) == x^(1/3)  — the TNBC vasculature K^(2/3) term.
    out = _apply("<ci>nthroot</ci>", X, N3)
    assert "std::pow(x, 1.0 / (3.0))" in out


def test_power_tag():
    assert _apply("<power/>", X, N2) == "std::pow(x, 2.0)"


def test_root_with_degree():
    out = _apply("<root/>", f"<degree>{N3}</degree>", X)
    assert "std::pow(x, 1.0 / 3.0)" in out


def test_ci_exported_max():
    # SimBiology exports max as <ci>max</ci>.
    assert "std::max" in _apply("<ci>max</ci>", X, '<cn type="integer">0</cn>')


def test_hyperbolic_ci():
    assert "std::tanh(x)" == _apply("<ci>tanh</ci>", X)


def test_unknown_ci_function_raises_with_context():
    with pytest.raises(ValueError, match="mysteryFunc"):
        _apply("<ci>mysteryFunc</ci>", X)
