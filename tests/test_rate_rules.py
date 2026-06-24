"""Tests for SBML <rateRule> on a constant=false parameter.

SimBiology models an abstract non-physical state (e.g. a dimensionless niche
fraction) as a parameter with ConstantValue=false plus a rate rule, rather than
as a species. The generator has no native rate-rule path, so it promotes the
parameter to an amount-tracked pseudo-species and injects a single synthetic
production reaction whose rate law IS the rate-rule expression. The existing
stoichiometry -> ydot -> Jacobian -> RHS machinery then evolves it like any
other state.

Regression: before this support, such a rate rule was silently dropped — the
parameter stayed frozen at its initial value, with no derivative emitted.
"""
import xml.etree.ElementTree as ET

from qsp_codegen.codegen import SBMLModel, SBML_NS, MATH_NS

SNS = SBML_NS.strip("{}")
MNS = MATH_NS.strip("{}")


def _build(params_xml: str, rules_xml: str) -> SBMLModel:
    doc = (
        f'<?xml version="1.0"?>'
        f'<sbml xmlns="{SNS}" level="2" version="4"><model id="m">'
        f"<listOfParameters>{params_xml}</listOfParameters>"
        f"<listOfRules>{rules_xml}</listOfRules>"
        f"</model></sbml>"
    )
    root = ET.fromstring(doc)
    m = object.__new__(SBMLModel)
    m.model = root.find(f"{SBML_NS}model")
    m.id_to_name = {}
    m.name_to_id = {}
    m.function_defs = {}
    m.species = []
    m.parameters = []
    m.reactions = []
    m._parse_parameters()
    m._parse_rate_rules()
    return m


_RATE_RULE = (
    '<rateRule variable="TLA">'
    f'<math xmlns="{MNS}"><apply><minus/>'
    "<apply><times/><ci>kf</ci>"
    "<apply><minus/><cn>1</cn><ci>TLA</ci></apply></apply>"
    "<apply><times/><ci>kd</ci><ci>TLA</ci></apply>"
    "</apply></math></rateRule>"
)
_PARAMS = (
    '<parameter id="TLA" name="TLA" value="0" constant="false"/>'
    '<parameter id="kf" name="kf" value="0.5" constant="true"/>'
    '<parameter id="kd" name="kd" value="0.05" constant="true"/>'
)


def test_rate_ruled_parameter_promoted_to_state():
    m = _build(_PARAMS, _RATE_RULE)

    # The rate-ruled parameter is removed from parameters (it's now a state);
    # the rate constants remain parameters.
    assert {p["name"] for p in m.parameters} == {"kf", "kd"}

    # It is added as an amount-tracked dimensionless pseudo-species.
    tla = [s for s in m.species if s["name"] == "TLA"]
    assert len(tla) == 1
    assert tla[0]["is_rate_ruled"] is True
    assert tla[0]["has_only_substance_units"] is True  # no volume scaling
    assert tla[0]["initial_value"] == 0.0

    # A single synthetic production reaction (∅ -> TLA) carries the rate-rule
    # expression as its rate law, so ydot[TLA] = +rate = the ODE.
    rxn = [r for r in m.reactions if r["product_names"] == ["TLA"]]
    assert len(rxn) == 1
    assert rxn[0]["reactant_names"] == []
    expr = rxn[0]["rate_law"].replace(" ", "")
    assert "kf" in expr and "kd" in expr and "TLA" in expr


def test_rate_rule_on_non_parameter_is_left_alone():
    # A rate rule whose target is not a parameter (e.g. a species, already a
    # state) is not promoted — the normal reaction path handles it.
    rr = (
        '<rateRule variable="some_species_id">'
        f'<math xmlns="{MNS}"><ci>kf</ci></math></rateRule>'
    )
    m = _build('<parameter id="kf" name="kf" value="0.5" constant="true"/>', rr)
    assert m.species == []
    assert m.reactions == []
    assert {p["name"] for p in m.parameters} == {"kf"}
