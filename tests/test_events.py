"""Tests for SBML <event> parsing: supported single-comparison triggers +
assignments, and the loud failures for delays / compound triggers.

End-to-end event behavior is additionally covered by the TNBC parity run
(3 events, including a `V_T.C < 0.5*cell -> V_T.C = 0.01*cell` reset), which
matches the SimBiology reference to within tolerance.
"""
import xml.etree.ElementTree as ET

import pytest

from qsp_codegen.codegen import SBMLModel, SBML_NS, MATH_NS

SNS = SBML_NS.strip("{}")
MNS = MATH_NS.strip("{}")


def _parse_events(event_xml: str):
    doc = (
        f'<?xml version="1.0"?>'
        f'<sbml xmlns="{SNS}" level="2" version="4"><model id="m">'
        f"<listOfEvents>{event_xml}</listOfEvents>"
        f"</model></sbml>"
    )
    root = ET.fromstring(doc)
    m = object.__new__(SBMLModel)
    m.model = root.find(f"{SBML_NS}model")
    m.id_to_name = {}
    m.function_defs = {}
    m.events = []
    m._parse_events()
    return m.events


def _trigger(op="lt", delay=""):
    return (
        f'<event id="e1"><trigger>'
        f'<math xmlns="{MNS}"><apply><{op}/><ci>C</ci><cn>0.5</cn></apply></math>'
        f"</trigger>{delay}"
        f'<listOfEventAssignments><eventAssignment variable="C">'
        f'<math xmlns="{MNS}"><cn>0.01</cn></math></eventAssignment>'
        f"</listOfEventAssignments></event>"
    )


def test_single_comparison_event_parsed():
    events = _parse_events(_trigger("lt"))
    assert len(events) == 1
    ev = events[0]
    assert ev["trigger_op"] == "lt"
    assert ev["trigger_left"] == "C"
    assert ev["trigger_right"] == "0.5"
    assert ev["assignments"][0]["variable_id"] == "C"
    assert "0.01" in ev["assignments"][0]["expression"]


@pytest.mark.parametrize("op", ["lt", "leq", "gt", "geq"])
def test_supported_comparison_ops(op):
    assert _parse_events(_trigger(op))[0]["trigger_op"] == op


def test_delay_raises():
    delay = f'<delay><math xmlns="{MNS}"><cn>1</cn></math></delay>'
    with pytest.raises(NotImplementedError, match="delay"):
        _parse_events(_trigger("lt", delay=delay))


def test_compound_trigger_raises():
    ev = (
        f'<event id="e1"><trigger><math xmlns="{MNS}"><apply><and/>'
        f"<apply><lt/><ci>C</ci><cn>0.5</cn></apply>"
        f"<apply><gt/><ci>C</ci><cn>0.1</cn></apply>"
        f'</apply></math></trigger>'
        f'<listOfEventAssignments><eventAssignment variable="C">'
        f'<math xmlns="{MNS}"><cn>0.01</cn></math></eventAssignment>'
        f"</listOfEventAssignments></event>"
    )
    with pytest.raises(NotImplementedError, match="unsupported trigger op"):
        _parse_events(ev)
