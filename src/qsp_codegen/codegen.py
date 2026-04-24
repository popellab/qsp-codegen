"""qsp_codegen.codegen — SBML → C++ CVODE ODE code generator.

Given a SimBiology-exported SBML Level 2 v4 file, emit the C++ sources
consumed by the CVODEBase-backed QSP simulator:

  QSP_enum.h     — species enum, parameter enum, QSP file param enum
  ODE_system.h   — class declaration (CVODEBase interface)
  ODE_system.cpp — RHS function, setup_class_parameters, init assignments
  QSPParam.h     — parameter reader declaration
  QSPParam.cpp   — parameter reader implementation
  qsp_params_xml_snippet.xml — XML snippet merged into consumer param files

Run only when QSP model structure changes, not for parameter value tweaks.
"""

import argparse
import os
import re
import sys
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Set, Tuple

SBML_NS = "{http://www.sbml.org/sbml/level2/version4}"
MATH_NS = "{http://www.w3.org/1998/Math/MathML}"


# =========================================================================
# SBML Parser
# =========================================================================

class SBMLModel:
    """Parse an SBML Level 2 v4 file into code-generation-ready structures."""

    def __init__(self, sbml_path: str):
        tree = ET.parse(sbml_path)
        root = tree.getroot()
        self.model = root.find(f"{SBML_NS}model")
        if self.model is None:
            raise ValueError("No <model> element found in SBML")

        # ID → name lookup (covers compartments, species, parameters, rules)
        self.id_to_name: Dict[str, str] = {}
        self.id_to_comp: Dict[str, str] = {}   # species id → compartment name
        self.name_to_id: Dict[str, str] = {}

        self.compartments: List[dict] = []
        self.species: List[dict] = []
        self.parameters: List[dict] = []
        self.reactions: List[dict] = []
        self.assignment_rules: List[dict] = []    # repeatedAssignment
        self.initial_assignments: List[dict] = []
        self.events: List[dict] = []
        self.unit_defs: Dict[str, float] = {}     # unit_id → SI factor

        self._parse_units()
        self._parse_compartments()
        self._parse_parameters()  # before species (species may ref param units)
        self._parse_species()
        self._parse_reactions()
        self._parse_rules()
        self._parse_initial_assignments()
        self._parse_events()

    # --- Units -----------------------------------------------------------

    def _parse_units(self):
        """Parse <listOfUnitDefinitions> → SI conversion factor per unit ID.

        Converts all quantities to SI base units (mole, metre, second, kg).
        This matches MATLAB SimBiology's internal unit resolution.
        The SBML model defines cell = 1.66e-24 moles (1/NA multiplier),
        so cells are stored as moles internally — same as MATLAB.
        """
        for ud in self.model.findall(f".//{SBML_NS}unitDefinition"):
            uid = ud.get("id")
            factor = 1.0
            for u in ud.findall(f".//{SBML_NS}unit"):
                kind = u.get("kind", "dimensionless")
                exp = float(u.get("exponent", "1"))
                scale = int(u.get("scale", "0"))
                mult = float(u.get("multiplier", "1"))
                base_si = {"metre": 1.0, "second": 1.0, "mole": 1.0,
                           "kilogram": 1.0, "gram": 1e-3,
                           "dimensionless": 1.0, "item": 1.0 / 6.02214076e23,
                           "litre": 1e-3, "liter": 1e-3}
                base = base_si.get(kind, 1.0)
                component = (mult * (10.0 ** scale) * base) ** exp
                factor *= component
            self.unit_defs[uid] = factor

        # Built-in units always available
        self.unit_defs.setdefault("dimensionless", 1.0)
        self.unit_defs.setdefault("substance", 1.0)
        self.unit_defs.setdefault("volume", 1.0)
        self.unit_defs.setdefault("time", 1.0)

    # Fallback SI factors for bare unit names not in listOfUnitDefinitions
    _BARE_SI = {"litre": 1e-3, "liter": 1e-3, "metre": 1.0, "meter": 1.0,
                "second": 1.0, "kilogram": 1.0, "gram": 1e-3,
                "mole": 1.0, "item": 1.0 / 6.02214076e23,
                "dimensionless": 1.0}

    def get_si_factor(self, unit_id: str) -> float:
        """Get SI conversion factor for a unit ID. Returns 1.0 if unknown."""
        if not unit_id:
            return 1.0
        factor = self.unit_defs.get(unit_id)
        if factor is not None:
            return factor
        factor = self._BARE_SI.get(unit_id)
        if factor is not None:
            return factor
        print(f"  WARNING: unknown unit '{unit_id}', assuming factor=1",
              file=sys.stderr)
        return 1.0

    # --- Compartments ----------------------------------------------------

    def _parse_compartments(self):
        for c in self.model.findall(f".//{SBML_NS}compartment"):
            cid = c.get("id")
            name = c.get("name")
            self.id_to_name[cid] = name
            self.name_to_id[name] = cid
            self.compartments.append({
                "id": cid,
                "name": name,
                "size": float(c.get("size", "1")),
                "spatial_dim": int(c.get("spatialDimensions", "3")),
                "units": c.get("units", ""),
            })

    # --- Parameters ------------------------------------------------------

    def _parse_parameters(self):
        for p in self.model.findall(
                f"{SBML_NS}listOfParameters/{SBML_NS}parameter"):
            pid = p.get("id")
            name = p.get("name")
            self.id_to_name[pid] = name
            self.name_to_id[name] = pid
            self.parameters.append({
                "id": pid,
                "name": name,
                "value": float(p.get("value", "0")),
                "units": p.get("units", ""),
                "constant": p.get("constant", "true") == "true",
            })

    # --- Species ---------------------------------------------------------

    def _parse_species(self):
        comp_id_to_name = {c["id"]: c["name"] for c in self.compartments}
        for sp in self.model.findall(f".//{SBML_NS}species"):
            sid = sp.get("id")
            name = sp.get("name")
            comp_id = sp.get("compartment")
            comp_name = comp_id_to_name.get(comp_id, comp_id)
            full_name = f"{comp_name}.{name}"

            self.id_to_name[sid] = full_name
            self.id_to_comp[sid] = comp_name
            self.name_to_id[full_name] = sid

            init_amount = sp.get("initialAmount")
            init_conc = sp.get("initialConcentration")
            if init_amount is not None:
                init = float(init_amount)
                is_concentration = False
            elif init_conc is not None:
                init = float(init_conc)
                is_concentration = True
            else:
                init = 0.0
                is_concentration = False
            self.species.append({
                "id": sid,
                "name": full_name,
                "base_name": name,
                "compartment": comp_name,
                "initial_value": init,
                "is_initial_concentration": is_concentration,
                "units": sp.get("substanceUnits", ""),
                "has_only_substance_units":
                    sp.get("hasOnlySubstanceUnits", "false") == "true",
            })

    # --- Reactions -------------------------------------------------------

    def _parse_reactions(self):
        for rxn in self.model.findall(f".//{SBML_NS}reaction"):
            rl = rxn.find(f"{SBML_NS}listOfReactants")
            pl = rxn.find(f"{SBML_NS}listOfProducts")

            reactant_ids = [sr.get("species")
                            for sr in (rl if rl is not None else [])]
            product_ids = [sr.get("species")
                           for sr in (pl if pl is not None else [])]

            # Parse MathML rate law
            kl = rxn.find(f"{SBML_NS}kineticLaw")
            rate_expr = ""
            if kl is not None:
                math = kl.find(f"{MATH_NS}math")
                if math is not None:
                    rate_expr = self._mathml_to_infix(math)

            self.reactions.append({
                "name": rxn.get("name", f"rxn_{len(self.reactions)}"),
                "reactant_ids": reactant_ids,
                "product_ids": product_ids,
                "reactant_names": [self.id_to_name.get(r, r) for r in reactant_ids],
                "product_names": [self.id_to_name.get(p, p) for p in product_ids],
                "rate_law": rate_expr,
            })

    # --- Assignment Rules -----------------------------------------------

    def _parse_rules(self):
        for ar in self.model.findall(f".//{SBML_NS}assignmentRule"):
            var_id = ar.get("variable")
            var_name = self.id_to_name.get(var_id, var_id)
            math = ar.find(f"{MATH_NS}math")
            expr = self._mathml_to_infix(math) if math is not None else "0"
            self.assignment_rules.append({
                "variable_id": var_id,
                "variable_name": var_name,
                "expression": expr,
            })

    # Supported comparison ops → (C++ op, root convention).
    # For `lt(a, b)` (a < b), the root function `b - a` crosses 0 upward
    # when a drops past b → rootsFound=+1 (matches "trigger went false→true").
    _CMP_OPS = {
        "lt":  ("<",  ("right", "left")),   # gout = right - left
        "leq": ("<=", ("right", "left")),
        "gt":  (">",  ("left", "right")),
        "geq": (">=", ("left", "right")),
    }

    def _parse_events(self):
        """Parse <listOfEvents>.

        Only supports single-comparison triggers (lt/leq/gt/geq) with no
        <delay>. Complex triggers (and/or of multiple conditions) would
        require multiple root functions per event; unsupported for now.
        """
        for ev in self.model.findall(f".//{SBML_NS}event"):
            name = ev.get("name") or ev.get("id")
            trig_el = ev.find(f"{SBML_NS}trigger")
            if trig_el is None:
                continue
            if ev.find(f"{SBML_NS}delay") is not None:
                raise NotImplementedError(
                    f"Event '{name}' has a <delay>; not supported by codegen."
                )
            math = trig_el.find(f"{MATH_NS}math")
            if math is None or len(list(math)) == 0:
                continue
            apply_node = list(math)[0]
            if apply_node.tag.replace(MATH_NS, "") != "apply":
                raise NotImplementedError(
                    f"Event '{name}' trigger must be a single comparison."
                )
            kids = list(apply_node)
            op = kids[0].tag.replace(MATH_NS, "")
            if op not in self._CMP_OPS:
                raise NotImplementedError(
                    f"Event '{name}' uses unsupported trigger op '{op}'. "
                    f"Supported: {sorted(self._CMP_OPS)}."
                )
            if len(kids) != 3:
                raise NotImplementedError(
                    f"Event '{name}' trigger must have exactly two arguments."
                )
            left = self._mathml_to_infix(kids[1])
            right = self._mathml_to_infix(kids[2])

            # SBML L2V4 default: initialValue="true" means the trigger is
            # considered already satisfied at t=0⁻. Events only fire on
            # subsequent false→true transitions. If initialValue="false"
            # and the trigger evaluates true at t=0⁺, it fires immediately.
            initial_value = trig_el.get("initialValue", "true") == "true"

            assignments = []
            for ea in ev.findall(f".//{SBML_NS}eventAssignment"):
                vid = ea.get("variable")
                vname = self.id_to_name.get(vid, vid)
                mr = ea.find(f"{MATH_NS}math")
                expr = self._mathml_to_infix(mr) if mr is not None else "0"
                assignments.append({
                    "variable_id": vid,
                    "variable_name": vname,
                    "expression": expr,
                })

            self.events.append({
                "name": name,
                "trigger_op": op,
                "trigger_left": left,
                "trigger_right": right,
                "trigger_initial_value": initial_value,
                "assignments": assignments,
            })

    def _parse_initial_assignments(self):
        for ia in self.model.findall(f".//{SBML_NS}initialAssignment"):
            sym_id = ia.get("symbol")
            sym_name = self.id_to_name.get(sym_id, sym_id)
            math = ia.find(f"{MATH_NS}math")
            expr = self._mathml_to_infix(math) if math is not None else "0"
            self.initial_assignments.append({
                "variable_id": sym_id,
                "variable_name": sym_name,
                "expression": expr,
            })

    # --- MathML → Infix Converter ----------------------------------------

    def _mathml_to_infix(self, node) -> str:
        """Convert a MathML <math> or <apply> tree to infix C++ string.

        Resolves SBML IDs to human-readable names.
        """
        tag = node.tag.replace(MATH_NS, "")

        if tag == "math":
            children = list(node)
            if len(children) == 1:
                return self._mathml_to_infix(children[0])
            # Multiple children at top level — shouldn't happen, join with ;
            return "; ".join(self._mathml_to_infix(c) for c in children)

        if tag == "cn":
            # Numeric constant
            typ = node.get("type", "real")
            text = (node.text or "").strip()
            if typ == "e-notation":
                # <cn type="e-notation"> mantissa <sep/> exponent </cn>
                sep = node.find(f"{MATH_NS}sep")
                if sep is not None:
                    mantissa = (node.text or "").strip()
                    exponent = (sep.tail or "").strip()
                    return f"{mantissa}e{exponent}"
            if typ == "integer":
                return f"{text}.0"  # Avoid C++ integer division
            return text

        if tag == "ci":
            # Identifier — resolve SBML ID to name
            sbml_id = (node.text or "").strip()
            return self.id_to_name.get(sbml_id, sbml_id)

        if tag == "apply":
            children = list(node)
            if not children:
                return "0"
            op_node = children[0]
            op_tag = op_node.tag.replace(MATH_NS, "")

            # SimBiology exports max/min as <ci>max</ci> instead of <max/>
            if op_tag == "ci":
                func_name = (op_node.text or "").strip()
                _CI_FUNCTIONS = {"max", "min", "abs", "floor", "ceil",
                                 "exp", "log", "ln", "sqrt", "power"}
                if func_name in _CI_FUNCTIONS:
                    op_tag = func_name

            # Handle <root> specially: extract <degree> child
            if op_tag == "root":
                degree_val = "2"  # default: square root
                value_nodes = []
                for c in children[1:]:
                    ct = c.tag.replace(MATH_NS, "")
                    if ct == "degree":
                        # <degree><cn>N</cn></degree>
                        inner = list(c)
                        if inner:
                            degree_val = self._mathml_to_infix(inner[0])
                    else:
                        value_nodes.append(c)
                if value_nodes:
                    val = self._mathml_to_infix(value_nodes[0])
                else:
                    val = "0"
                if degree_val == "2":
                    return f"std::sqrt({val})"
                return f"std::pow({val}, 1.0 / {degree_val})"

            args = [self._mathml_to_infix(c) for c in children[1:]]
            return self._apply_op(op_tag, args)

        if tag == "piecewise":
            return self._convert_piecewise(node)

        if tag in ("sep", "degree"):
            return ""

        # Fallback
        return f"/* unknown MathML tag: {tag} */"

    def _apply_op(self, op: str, args: List[str]) -> str:
        """Convert a MathML <apply> operation to infix C++."""
        # Binary arithmetic
        if op == "plus":
            return " + ".join(f"({a})" if " - " in a or " + " in a else a
                              for a in args) if len(args) > 1 else args[0]
        if op == "minus":
            if len(args) == 1:
                return f"(-{args[0]})"
            return f"({args[0]} - {args[1]})"
        if op == "times":
            return " * ".join(f"({a})" if (" + " in a or " - " in a) else a
                              for a in args)
        if op == "divide":
            num = f"({args[0]})" if (" + " in args[0] or " - " in args[0]) else args[0]
            den = f"({args[1]})" if (" + " in args[1] or " - " in args[1]
                                     or " * " in args[1] or " / " in args[1]) else args[1]
            return f"{num} / {den}"
        if op == "power":
            return f"std::pow({args[0]}, {args[1]})"

        # Comparison
        if op == "lt":
            return f"({args[0]} < {args[1]})"
        if op == "leq":
            return f"({args[0]} <= {args[1]})"
        if op == "gt":
            return f"({args[0]} > {args[1]})"
        if op == "geq":
            return f"({args[0]} >= {args[1]})"
        if op == "eq":
            return f"({args[0]} == {args[1]})"
        if op == "neq":
            return f"({args[0]} != {args[1]})"

        # Logic
        if op == "and":
            return " && ".join(f"({a})" for a in args)
        if op == "or":
            return " || ".join(f"({a})" for a in args)
        if op == "not":
            return f"(!({args[0]}))"

        # Functions
        func_map = {
            "ln": "std::log", "log": "std::log", "log2": "std::log2",
            "log10": "std::log10", "exp": "std::exp", "sqrt": "std::sqrt",
            "abs": "std::abs", "floor": "std::floor", "ceiling": "std::ceil",
            "sin": "std::sin", "cos": "std::cos", "tan": "std::tan",
        }
        if op in func_map:
            return f"{func_map[op]}({', '.join(args)})"

        if op == "max":
            if len(args) == 2:
                return f"std::max({args[0]}, {args[1]})"
            # Nested max for >2 args
            result = args[0]
            for a in args[1:]:
                result = f"std::max({result}, {a})"
            return result
        if op == "min":
            if len(args) == 2:
                return f"std::min({args[0]}, {args[1]})"
            result = args[0]
            for a in args[1:]:
                result = f"std::min({result}, {a})"
            return result

        if op == "root":
            # <apply><root/><degree><cn>2</cn></degree><ci>x</ci></apply>
            # args[0] is degree, args[1] is the value (if degree present)
            if len(args) == 2:
                return f"std::pow({args[1]}, 1.0 / {args[0]})"
            return f"std::sqrt({args[0]})"

        return f"/* unknown op: {op} */({', '.join(args)})"

    def _convert_piecewise(self, node) -> str:
        """Convert MathML <piecewise> to C++ ternary."""
        pieces = node.findall(f"{MATH_NS}piece")
        otherwise = node.find(f"{MATH_NS}otherwise")

        parts = []
        for piece in pieces:
            children = list(piece)
            if len(children) >= 2:
                val = self._mathml_to_infix(children[0])
                cond = self._mathml_to_infix(children[1])
                parts.append((cond, val))

        default = "0.0"
        if otherwise is not None:
            children = list(otherwise)
            if children:
                default = self._mathml_to_infix(children[0])

        # Build nested ternary
        result = default
        for cond, val in reversed(parts):
            result = f"({cond} ? {val} : {result})"
        return result


# =========================================================================
# Identifier classification for C++ macro wrapping
# =========================================================================

def _sanitize(name: str) -> str:
    """V_T.C1 → V_T_C1"""
    return name.replace(".", "_")


def _collect_rule_deps(
    exprs: List[str],
    rules_by_name: Dict[str, dict],
    sbml: "SBMLModel" = None,
) -> List[dict]:
    """Return the assignment rules transitively referenced by `exprs`, in
    topological order (dependencies first). Rules not referenced are omitted.
    """
    needed: Set[str] = set()

    def walk(vn: str):
        if vn in needed:
            return
        needed.add(vn)
        if vn in rules_by_name:
            rhs = rules_by_name[vn]["expression"]
            for other in rules_by_name:
                if other == vn:
                    continue
                pat = re.escape(other)
                if re.search(r"(?<![a-zA-Z0-9_.])" + pat + r"(?![a-zA-Z0-9_.])", rhs):
                    walk(other)

    for expr in exprs:
        for vn in rules_by_name:
            pat = re.escape(vn)
            if re.search(r"(?<![a-zA-Z0-9_.])" + pat + r"(?![a-zA-Z0-9_.])", expr):
                walk(vn)

    needed_list = [r for r in rules_by_name.values() if r["variable_name"] in needed]
    return order_rules(needed_list, needed, sbml=sbml)


def _emit_event_code(lines: List[str], sbml: "SBMLModel") -> None:
    """Emit g(), triggerComponentEvaluate(), eventEvaluate(), eventExecution()."""
    if not sbml.events:
        lines.append('bool ODE_system::triggerComponentEvaluate(int, realtype, bool){ return false; }')
        lines.append('bool ODE_system::eventEvaluate(int){ return false; }')
        lines.append('bool ODE_system::eventExecution(int, bool, realtype&){ return false; }')
        lines.append('')
        return

    # Two mappings: (1) RHS mapping for g() where `y` is in scope → SPVAR;
    # (2) `_y`-reading mapping for member funcs (no y param) where species
    # references become NV_DATA_S(_y)[SP_xxx] (same pattern as eval_init_assignment).
    rhs_mapping = classify_identifiers(sbml)
    rule_vars = {r["variable_name"] for r in sbml.assignment_rules}
    mem_mapping = {}
    for sp in sbml.species:
        san = _sanitize(sp["name"])
        if sp.get("has_only_substance_units", False):
            mem_mapping[sp["name"]] = f'NV_DATA_S(_y)[SP_{san}]'
        else:
            comp_name = sp["compartment"]
            comp_san = _sanitize(comp_name)
            vol_expr = (f'AUX_VAR_{comp_san}' if comp_name in rule_vars
                        else f'PARAM(P_{comp_san})')
            mem_mapping[sp["name"]] = f'(NV_DATA_S(_y)[SP_{san}] / {vol_expr})'
    for r in sbml.assignment_rules:
        mem_mapping[r["variable_name"]] = f'AUX_VAR_{_sanitize(r["variable_name"])}'
    for p in sbml.parameters:
        if p["name"] not in rule_vars:
            mem_mapping[p["name"]] = f'PARAM(P_{_sanitize(p["name"])})'
    for c in sbml.compartments:
        if c["name"] in rule_vars:
            pass
        elif c["name"] not in mem_mapping:
            mem_mapping[c["name"]] = f'PARAM(P_{_sanitize(c["name"])})'

    rules_by_name = {r["variable_name"]: r for r in sbml.assignment_rules}

    def emit_aux_vars(exprs: List[str], mapping: dict) -> List[str]:
        """Emit the AUX_VAR computations needed to evaluate `exprs`."""
        ordered = _collect_rule_deps(exprs, rules_by_name, sbml=sbml)
        out = []
        for r in ordered:
            vn = r["variable_name"]
            san = _sanitize(vn)
            cpp_expr = wrap_expression(r["expression"], mapping)
            out.append(f'    realtype AUX_VAR_{san} = {cpp_expr};')
        return out

    # --- g() ----------------------------------------------------------
    lines.append('int ODE_system::g(realtype t, N_Vector y, realtype* gout, void* user_data){')
    trigger_exprs = []
    for ev in sbml.events:
        trigger_exprs.extend([ev["trigger_left"], ev["trigger_right"]])
    lines.extend(emit_aux_vars(trigger_exprs, rhs_mapping))
    for i, ev in enumerate(sbml.events):
        _, (minus_side, plus_side) = SBMLModel._CMP_OPS[ev["trigger_op"]]
        # gout = <minus_side> - <plus_side> so it increases when trigger becomes true
        expr_minus = wrap_expression(ev[f"trigger_{minus_side}"], rhs_mapping)
        expr_plus = wrap_expression(ev[f"trigger_{plus_side}"], rhs_mapping)
        lines.append(f'    // Event {i}: {ev["name"]} — trigger {ev["trigger_left"]} '
                     f'{ev["trigger_op"]} {ev["trigger_right"]}')
        lines.append(f'    gout[{i}] = ({expr_minus}) - ({expr_plus});')
    lines.append('    return 0;')
    lines.append('}')
    lines.append('')

    # --- triggerComponentEvaluate() ----------------------------------
    lines.append('bool ODE_system::triggerComponentEvaluate(int i, realtype, bool){')
    lines.extend(emit_aux_vars(trigger_exprs, mem_mapping))
    lines.append('    switch (i) {')
    for i, ev in enumerate(sbml.events):
        cpp_op, _ = SBMLModel._CMP_OPS[ev["trigger_op"]]
        left = wrap_expression(ev["trigger_left"], mem_mapping)
        right = wrap_expression(ev["trigger_right"], mem_mapping)
        lines.append(f'    case {i}: return ({left}) {cpp_op} ({right});')
    lines.append('    }')
    lines.append('    return false;')
    lines.append('}')
    lines.append('')

    # --- eventEvaluate() ---------------------------------------------
    # Single-component-per-event: 1:1 with _trigger_element_satisfied.
    lines.append('bool ODE_system::eventEvaluate(int i){')
    lines.append('    return (i >= 0 && i < _nroot) ? _trigger_element_satisfied[i] : false;')
    lines.append('}')
    lines.append('')

    # --- eventExecution() --------------------------------------------
    # Apply assignments. For amount species (hasOnlySubstanceUnits=true),
    # the RHS evaluates directly to SI amount (because PARAMs are stored
    # in SI). For concentration species, multiply by compartment to get
    # amount. Writes both _species_var and _y so post-event RHS reads see
    # the new value.
    lines.append('bool ODE_system::eventExecution(int i, bool, realtype&){')
    # Collect all assignment RHS expressions for dep-emission
    assign_exprs = [a["expression"] for ev in sbml.events for a in ev["assignments"]]
    lines.extend(emit_aux_vars(assign_exprs, mem_mapping))
    lines.append('    switch (i) {')
    sp_by_name = {sp["name"]: sp for sp in sbml.species}
    for i, ev in enumerate(sbml.events):
        lines.append(f'    case {i}:  // {ev["name"]}')
        for a in ev["assignments"]:
            vn = a["variable_name"]
            sp = sp_by_name.get(vn)
            if sp is None:
                raise NotImplementedError(
                    f"Event '{ev['name']}' assigns to non-species '{vn}'; not supported."
                )
            san = _sanitize(vn)
            cpp_expr = wrap_expression(a["expression"], mem_mapping)
            if sp.get("has_only_substance_units", False):
                rhs = cpp_expr
            else:
                comp_name = sp["compartment"]
                comp_san = _sanitize(comp_name)
                vol = (f'AUX_VAR_{comp_san}' if comp_name in rule_vars
                       else f'PARAM(P_{comp_san})')
                rhs = f'({cpp_expr}) * {vol}'
            lines.append(f'        _species_var[SP_{san}] = {rhs};')
            lines.append(f'        NV_DATA_S(_y)[SP_{san}] = _species_var[SP_{san}];')
        lines.append('        break;')
    lines.append('    }')
    lines.append('    return false;  // no delay')
    lines.append('}')
    lines.append('')


def classify_identifiers(sbml: SBMLModel) -> dict:
    """Build identifier → C++ macro mapping.

    Returns dict: human_name → "SPVAR(SP_xxx)" | "PARAM(P_xxx)" | "AUX_VAR_xxx"

    SBML convention: species with hasOnlySubstanceUnits=false appear as
    CONCENTRATION in kinetic law expressions (not amount). CVODE tracks
    amounts. So concentration-based species are mapped to SPVAR/V_comp,
    giving concentration when substituted into rate law formulas.
    """
    mapping = {}

    rule_vars = {r["variable_name"] for r in sbml.assignment_rules}

    for sp in sbml.species:
        san = _sanitize(sp["name"])
        if sp.get("has_only_substance_units", False):
            # Amount-based species: SPVAR stores amount → use directly
            mapping[sp["name"]] = f'SPVAR(SP_{san})'
        else:
            # Concentration-based species: SPVAR stores amount, but kinetic
            # laws expect concentration → divide by compartment volume/area.
            comp_name = sp["compartment"]
            comp_san = _sanitize(comp_name)
            if comp_name in rule_vars:
                vol_expr = f'AUX_VAR_{comp_san}'
            else:
                vol_expr = f'PARAM(P_{comp_san})'
            mapping[sp["name"]] = f'(SPVAR(SP_{san}) / {vol_expr})'

    for r in sbml.assignment_rules:
        vn = r["variable_name"]
        mapping[vn] = f'AUX_VAR_{_sanitize(vn)}'

    # Parameters that are NOT overridden by rules
    for p in sbml.parameters:
        if p["name"] not in rule_vars:
            mapping[p["name"]] = f'PARAM(P_{_sanitize(p["name"])})'

    # Compartment names that appear as rule variables → AUX_VAR
    # Compartment names that are just fixed params → PARAM
    for c in sbml.compartments:
        if c["name"] in rule_vars:
            pass  # already handled above
        elif c["name"] not in mapping:
            mapping[c["name"]] = f'PARAM(P_{_sanitize(c["name"])})'

    return mapping


def wrap_expression(expr: str, mapping: dict) -> str:
    """Replace human-readable names in an expression with C++ macros."""
    # Sort by length descending for longest-match-first
    sorted_names = sorted(mapping.keys(), key=len, reverse=True)

    result = expr
    replacements = []
    for name in sorted_names:
        pattern = re.escape(name)
        full_pattern = r'(?<![a-zA-Z0-9_\.])' + pattern + r'(?![a-zA-Z0-9_\.])'
        placeholder = f"__PH{len(replacements)}__"
        if re.search(full_pattern, result):
            result = re.sub(full_pattern, placeholder, result)
            replacements.append((placeholder, mapping[name]))

    for ph, repl in replacements:
        result = result.replace(ph, repl)

    return result


# =========================================================================
# Rule dependency ordering
# =========================================================================

def order_rules(rules: List[dict], all_rule_names: Set[str],
                sbml: "SBMLModel" = None) -> List[dict]:
    """Topologically sort assignment rules so dependencies come first.

    Also accounts for implicit dependencies introduced by concentration-based
    species: when a species X in compartment C is used in a rule expression,
    and X is concentration-based (hasOnlySubstanceUnits=False), the codegen
    emits (SPVAR(SP_X) / AUX_VAR_C). So the rule implicitly depends on the
    compartment rule C. We only add this dependency when C is itself a rule
    AND the species with that compartment is concentration-based.
    """
    name_to_rule = {r["variable_name"]: r for r in rules}

    # Build mapping: compartment_name → set of concentration-based species names
    # that are prefixed "compartment_name.species_short" in SBML expressions.
    # We only add a compartment dependency if at least one concentration-based
    # species in that compartment appears in the rule expression.
    conc_species_by_comp: dict = {}  # comp_name → set of full species names
    if sbml is not None:
        for sp in sbml.species:
            if not sp.get("has_only_substance_units", False):
                comp = sp["compartment"]
                conc_species_by_comp.setdefault(comp, set()).add(sp["name"])

    deps = {}
    for r in rules:
        vn = r["variable_name"]
        deps[vn] = set()
        expr = r["expression"]
        for other_vn in all_rule_names:
            if other_vn == vn:
                continue
            pattern = re.escape(other_vn)
            # Standard dependency: other_vn appears as a standalone token
            if re.search(r'(?<![a-zA-Z0-9_\.])' + pattern + r'(?![a-zA-Z0-9_\.])',
                         expr):
                deps[vn].add(other_vn)
            # Compartment-prefix dependency: other_vn is a compartment that is
            # also a rule, and at least one concentration-based species in that
            # compartment appears in this expression (as "other_vn.speciesName").
            # The codegen will emit (SPVAR / AUX_VAR_other_vn) for those species,
            # so other_vn must be defined before this rule.
            elif other_vn in name_to_rule and other_vn in conc_species_by_comp:
                for sp_name in conc_species_by_comp[other_vn]:
                    sp_pat = re.escape(sp_name)
                    if re.search(r'(?<![a-zA-Z0-9_\.])' + sp_pat + r'(?![a-zA-Z0-9_\.])',
                                 expr):
                        deps[vn].add(other_vn)
                        break

    # Kahn's algorithm
    in_degree = {vn: len(deps[vn]) for vn in name_to_rule}
    queue = sorted([vn for vn in name_to_rule if in_degree[vn] == 0])
    ordered = []
    while queue:
        node = queue.pop(0)
        ordered.append(name_to_rule[node])
        for other in name_to_rule:
            if node in deps.get(other, set()):
                deps[other].discard(node)
                in_degree[other] -= 1
                if in_degree[other] == 0:
                    queue.append(other)
        queue.sort()

    if len(ordered) != len(rules):
        remaining = set(name_to_rule.keys()) - {r["variable_name"] for r in ordered}
        print(f"WARNING: circular rule dependencies: {remaining}", file=sys.stderr)
        for vn in sorted(remaining):
            ordered.append(name_to_rule[vn])

    return ordered


# =========================================================================
# Volume scaling for ydot
# =========================================================================

def needs_volume_scaling(species: dict, sbml: SBMLModel) -> Optional[str]:
    """Check if a species' ydot needs 1/V_compartment scaling.

    Since we track AMOUNTS for all species (concentration * V for conc species),
    and rate laws substitute concentration = SPVAR/V, the resulting fluxes are
    in substance/time. ydot = sum(fluxes) with no volume scaling needed.

    Returns None (no scaling) for all species.
    """
    return None


# =========================================================================
# Stoichiometry
# =========================================================================

def build_stoichiometry(reactions: list, species_names: set) -> Dict[str, List[Tuple[int, int]]]:
    """species_name → [(flux_1based, +1/-1), ...]"""
    stoich: Dict[str, List[Tuple[int, int]]] = {sp: [] for sp in species_names}
    for i, rxn in enumerate(reactions):
        flux_idx = i + 1
        for sp in rxn["reactant_names"]:
            if sp in stoich:
                stoich[sp].append((flux_idx, -1))
        for sp in rxn["product_names"]:
            if sp in stoich:
                stoich[sp].append((flux_idx, +1))
    return stoich


# =========================================================================
# Symbolic Jacobian (for KLU sparse linear solver)
# =========================================================================
# Builds df/dy from the SBML reaction rates + assignment rules using sympy,
# applies CSE to compress output, and emits a CVODE-compatible jac() function
# plus a static CSC sparsity pattern. Measured sparsity on PDAC_model.sbml is
# ~1296 nnz out of 164² = ~4.8% density, so sparse LU beats dense by ~20x
# in flops.

def _build_sympy_context(sbml: SBMLModel):
    """Build sympy symbol registry + a C++-infix-to-sympy parser.

    The sympy path differs from the main codegen's string-substitution path
    in one critical way: SBML concentration species appear in rate laws as
    CONCENTRATION (= amount / compartment_volume), but CVODE integrates
    AMOUNT. For the analytical Jacobian ∂f_i/∂y_j to be correct under y=amount,
    we substitute concentration references with (amount_symbol / comp_symbol)
    in sympy, then let sympy's expansion/differentiation handle the chain rule
    correctly when the compartment rule gets substituted.

    Returns a dict with:
        - orig_to_san: map SBML name ("V_T.C1") → sanitized ("V_T_C1")
        - symbols:     map sanitized name → sp.Symbol
        - species_syms / species_set: the species symbols in enum order
        - to_sympy(expr_str) → sp.Expr
        - sympy_wrap_mapping: orig-name → C++ macro for wrapping ccode
          output in the jac() function (uses raw SPVAR, no
          compartment-division, since concentration was done in sympy)
    """
    import sympy as sp

    all_names: Set[str] = set()
    for s_ in sbml.species:     all_names.add(s_["name"])
    for p  in sbml.parameters:  all_names.add(p["name"])
    for c  in sbml.compartments: all_names.add(c["name"])
    for r  in sbml.assignment_rules: all_names.add(r["variable_name"])

    orig_to_san = {n: _sanitize(n) for n in all_names}
    symbols: Dict[str, "sp.Symbol"] = {_sanitize(n): sp.Symbol(_sanitize(n)) for n in all_names}

    # Per-identifier substitution strings for the sympy parser. Dotted
    # concentration species expand to (amount / compartment); amount species
    # and plain identifiers just sanitize their dots.
    id_infix: Dict[str, str] = {}
    rule_var_set = {r["variable_name"] for r in sbml.assignment_rules}
    for sp_ in sbml.species:
        nm = sp_["name"]
        san = _sanitize(nm)
        if sp_.get("has_only_substance_units", False):
            id_infix[nm] = san
        else:
            comp = sp_["compartment"]
            comp_san = _sanitize(comp)
            id_infix[nm] = f"({san} / {comp_san})"
    # All non-species identifiers just sanitize (dotted params/rule vars are rare
    # but possible; stay robust).
    for nm in all_names:
        if nm not in id_infix:
            id_infix[nm] = _sanitize(nm)

    # Longest first so V_T.C1 wins over V_T.
    sorted_orig = sorted(id_infix.keys(), key=len, reverse=True)

    _fn_rewrites = [
        (re.compile(r"\bstd::"), ""),
        (re.compile(r"\bmax\b"), "Max"),
        (re.compile(r"\bmin\b"), "Min"),
        (re.compile(r"\babs\b"), "Abs"),
        (re.compile(r"\bln\b"), "log"),
        (re.compile(r"\bceil\b"), "ceiling"),
    ]

    def to_sympy(expr: str) -> "sp.Expr":
        s = expr
        for pat, repl in _fn_rewrites:
            s = pat.sub(repl, s)
        # Two-phase substitution: first insert placeholders for each known
        # identifier, then realize placeholders to their sympy infix. This
        # avoids 'V_T' eating part of 'V_T_C1' that we just produced.
        slots: List[Tuple[str, str]] = []
        for orig in sorted_orig:
            pat = r'(?<![A-Za-z0-9_\.])' + re.escape(orig) + r'(?![A-Za-z0-9_\.])'
            if re.search(pat, s):
                ph = f"__SYMIDX{len(slots)}__"
                s = re.sub(pat, ph, s)
                slots.append((ph, id_infix[orig]))
        for ph, infix in slots:
            s = s.replace(ph, infix)
        return sp.sympify(s, locals=symbols)

    species_syms = [symbols[orig_to_san[sp_["name"]]] for sp_ in sbml.species]

    # Mapping for wrapping the sympy-emitted C++ code back to SPVAR/PARAM
    # macros. Differs from classify_identifiers': concentration species are
    # just their raw amount (SPVAR(SP_X)) because the compartment division
    # was already done in sympy-space and has been expanded symbolically.
    sympy_wrap_mapping: Dict[str, str] = {}
    for sp_ in sbml.species:
        sympy_wrap_mapping[sp_["name"]] = f'SPVAR(SP_{_sanitize(sp_["name"])})'
    for p in sbml.parameters:
        if p["name"] not in rule_var_set:
            sympy_wrap_mapping[p["name"]] = f'PARAM(P_{_sanitize(p["name"])})'
    for c in sbml.compartments:
        if c["name"] not in rule_var_set and c["name"] not in sympy_wrap_mapping:
            sympy_wrap_mapping[c["name"]] = f'PARAM(P_{_sanitize(c["name"])})'
    # Assignment-rule variables shouldn't appear post-expand, but map them as
    # a defensive fallback so sympy output referencing them is visibly broken
    # rather than silently emitting a bare identifier.
    for r in sbml.assignment_rules:
        sympy_wrap_mapping.setdefault(r["variable_name"],
                                      f'/*UNEXPANDED_RULE_{_sanitize(r["variable_name"])}*/')

    return {
        "orig_to_san": orig_to_san,
        "symbols": symbols,
        "species_syms": species_syms,
        "species_set": set(species_syms),
        "to_sympy": to_sympy,
        "sympy_wrap_mapping": sympy_wrap_mapping,
    }


def _postprocess_ccode(cpp: str) -> str:
    """Translate sympy's default ccode output to match the rest of the
    codegen's conventions (std:: prefix on math functions)."""
    # Order matters: log2/log10 before log.
    subs = [
        (r"\blog10\(", "std::log10("),
        (r"\blog2\(",  "std::log2("),
        (r"\blog\(",   "std::log("),
        (r"\bexp\(",   "std::exp("),
        (r"\bsqrt\(",  "std::sqrt("),
        (r"\bpow\(",   "std::pow("),
        (r"\bfmax\(",  "std::fmax("),
        (r"\bfmin\(",  "std::fmin("),
        (r"\bfabs\(",  "std::fabs("),
        (r"\bsin\(",   "std::sin("),
        (r"\bcos\(",   "std::cos("),
        (r"\btan\(",   "std::tan("),
    ]
    for pat, repl in subs:
        cpp = re.sub(pat, repl, cpp)
    return cpp


def _sympy_wrap(cpp_code: str, mapping: Dict[str, str], orig_to_san: Dict[str, str]) -> str:
    """Replace sanitized identifiers in sympy's ccode output with the
    C++ macros (SPVAR/PARAM/AUX_VAR) from the codegen's mapping."""
    # mapping is keyed by original names; we need sanitized-key lookups
    # because sympy emits sanitized identifier names.
    san_to_macro: Dict[str, str] = {}
    for orig, macro in mapping.items():
        san_to_macro[orig_to_san.get(orig, orig.replace(".", "_"))] = macro

    # Two-pass placeholder substitution (same trick wrap_expression uses).
    sorted_sans = sorted(san_to_macro, key=len, reverse=True)
    placeholders: List[Tuple[str, str]] = []
    for san in sorted_sans:
        pat = r'(?<![A-Za-z0-9_])' + re.escape(san) + r'(?![A-Za-z0-9_])'
        if re.search(pat, cpp_code):
            ph = f"__JPH{len(placeholders)}__"
            cpp_code = re.sub(pat, ph, cpp_code)
            placeholders.append((ph, san_to_macro[san]))
    for ph, macro in placeholders:
        cpp_code = cpp_code.replace(ph, macro)
    return cpp_code


def compute_jacobian(sbml: SBMLModel, mapping: Dict[str, str]) -> dict:
    """Derive df/dy symbolically, run CSE, prepare C++ strings.

    Returns:
        {
          "nnz":         int,                         # nonzero count
          "col_ptrs":    List[int] length neq+1,      # CSC column pointers
          "row_indices": List[int] length nnz,        # CSC row indices
          "cse_items":   List[Tuple[str, str]],       # [(aux_var, cpp_expr)]
          "entry_cpp":   List[str] length nnz,        # J values in CSC order
        }
    """
    import sympy as sp

    ctx = _build_sympy_context(sbml)
    symbols = ctx["symbols"]
    orig_to_san = ctx["orig_to_san"]
    species_syms = ctx["species_syms"]
    species_set = ctx["species_set"]
    to_sympy = ctx["to_sympy"]

    rate_sym = [to_sympy(r["rate_law"]) for r in sbml.reactions]
    rule_exprs = {symbols[orig_to_san[r["variable_name"]]]: to_sympy(r["expression"])
                  for r in sbml.assignment_rules}

    def expand(expr: "sp.Expr") -> "sp.Expr":
        # Fixpoint-substitute until no rule symbol remains. Assignment-rule
        # DAG is acyclic by SBML contract; 50-pass cap is a safety net.
        prev = None
        e = expr
        for _ in range(50):
            e = e.subs(rule_exprs)
            if e == prev:
                return e
            prev = e
        return e

    rate_expanded = [expand(r) for r in rate_sym]

    sp_names = {sp_["name"] for sp_ in sbml.species}
    stoich = build_stoichiometry(sbml.reactions, sp_names)

    ydot: Dict["sp.Symbol", "sp.Expr"] = {}
    for sp_ in sbml.species:
        sym = symbols[orig_to_san[sp_["name"]]]
        entries = stoich.get(sp_["name"], [])
        if not entries:
            ydot[sym] = sp.Integer(0)
            continue
        flux_sum = sum((sign * rate_expanded[fi-1] for fi, sign in entries), sp.Integer(0))
        comp_name = needs_volume_scaling(sp_, sbml)
        if comp_name is not None:
            v_sym = symbols[orig_to_san[comp_name]]
            v_expr = expand(v_sym) if v_sym in rule_exprs else v_sym
            ydot[sym] = flux_sum / v_expr
        else:
            ydot[sym] = flux_sum

    # Walk columns (j) so col_ptrs drops out naturally.
    entries_by_col: Dict[int, List[Tuple[int, "sp.Expr"]]] = {
        j: [] for j in range(len(sbml.species))
    }
    for i, sp_i in enumerate(sbml.species):
        sym_i = symbols[orig_to_san[sp_i["name"]]]
        expr_i = ydot[sym_i]
        if expr_i == 0:
            continue
        free = expr_i.free_symbols & species_set
        for sym_j in free:
            j = species_syms.index(sym_j)
            d = sp.diff(expr_i, sym_j)
            if d != 0:
                entries_by_col[j].append((i, d))

    # KLU + CVODE work on M = I - gamma*J for BDF-Newton iteration. SUNDIALS
    # stores M in-place inside J's allocated sparsity pattern, which means
    # we must include every diagonal (j, j) entry — even ones where the
    # analytical Jacobian is exactly 0. Otherwise SUNMatScaleAddI has no
    # slot to place the `1` and the linear-solver setup fails with
    # "At t=0, the setup routine failed in an unrecoverable manner."
    for j in range(len(sbml.species)):
        if not any(i == j for (i, _) in entries_by_col[j]):
            entries_by_col[j].append((j, sp.Integer(0)))

    col_ptrs: List[int] = [0]
    row_indices: List[int] = []
    entry_exprs: List["sp.Expr"] = []
    for j in range(len(sbml.species)):
        for (i, d) in sorted(entries_by_col[j], key=lambda p: p[0]):
            row_indices.append(i)
            entry_exprs.append(d)
        col_ptrs.append(len(row_indices))

    # CSE compresses ~1.7 MB of naive ccode down to ~40 KB by extracting
    # shared subexpressions. optimizations=None (rather than "basic")
    # because "basic" aggressively factors out 1/x terms across multiple
    # entries — e.g. it turns (n-1)*pow(x, n-1)*A into aux=(n-1)/x,
    # pow(x, n)*A*aux. That aux diverges at x=0 even though the original
    # entry was finite, which breaks KLU's numeric factorization whenever
    # a species touches zero. optimizations=None keeps each entry
    # arithmetically isolated; the extra ~250 aux vars are negligible,
    # and this CSE pass is also ~70x faster (~0.5s vs ~38s).
    replacements, reduced = sp.cse(entry_exprs, optimizations=None)

    wrap_map = ctx["sympy_wrap_mapping"]

    cse_items: List[Tuple[str, str]] = []
    for aux_sym, aux_rhs in replacements:
        cpp = _postprocess_ccode(sp.ccode(aux_rhs))
        cpp = _sympy_wrap(cpp, wrap_map, orig_to_san)
        cse_items.append((str(aux_sym), cpp))

    entry_cpp: List[str] = []
    for r in reduced:
        cpp = _postprocess_ccode(sp.ccode(r))
        cpp = _sympy_wrap(cpp, wrap_map, orig_to_san)
        entry_cpp.append(cpp)

    return {
        "nnz": len(row_indices),
        "col_ptrs": col_ptrs,
        "row_indices": row_indices,
        "cse_items": cse_items,
        "entry_cpp": entry_cpp,
    }


def gen_jacobian_cpp(sbml: SBMLModel, info: dict) -> str:
    """Emit the static CSC arrays + the ODE_system::jac function body."""
    neq = len(sbml.species)
    nnz = info["nnz"]

    def _emit_int_array(name: str, values: List[int], per_line: int = 16) -> List[str]:
        out = [f"const sunindextype ODE_system::{name}[] = {{"]
        for i in range(0, len(values), per_line):
            chunk = values[i:i+per_line]
            out.append("    " + ", ".join(str(x) for x in chunk) + ",")
        out.append("};")
        return out

    lines: List[str] = []
    lines.append("// --- Analytical Jacobian (generated) ---")
    lines.extend(_emit_int_array("_jac_col_ptrs", info["col_ptrs"]))
    lines.append("")
    lines.extend(_emit_int_array("_jac_row_indices", info["row_indices"]))
    lines.append("")
    lines.append("int ODE_system::jac(realtype t, N_Vector y, N_Vector fy,")
    lines.append("                     SUNMatrix J, void *user_data,")
    lines.append("                     N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){")
    lines.append("    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);")
    lines.append("    sunindextype *colptrs = SUNSparseMatrix_IndexPointers(J);")
    lines.append("    sunindextype *rowvals = SUNSparseMatrix_IndexValues(J);")
    lines.append("    realtype *data = SUNSparseMatrix_Data(J);")
    lines.append("")
    lines.append("    // Restamp the sparsity pattern every call: cheap, and")
    lines.append("    // robust against CVODE reallocating between solves.")
    lines.append(f"    for (sunindextype k = 0; k <= {neq}; ++k) colptrs[k] = _jac_col_ptrs[k];")
    lines.append(f"    for (sunindextype k = 0; k < {nnz}; ++k) rowvals[k] = _jac_row_indices[k];")
    lines.append("")
    lines.append("    // Common subexpressions extracted by sympy.cse() from the")
    lines.append("    // ~1.7 MB naive Jacobian — compresses to ~40 KB of output.")
    for aux_name, cpp in info["cse_items"]:
        lines.append(f"    realtype {aux_name} = {cpp};")
    lines.append("")
    lines.append("    // Jacobian values (CSC order, matches _jac_row_indices / _jac_col_ptrs):")
    for k, cpp in enumerate(info["entry_cpp"]):
        lines.append(f"    data[{k}] = {cpp};")
    lines.append("")
    lines.append("    // Boundary-state safety clamp. Analytical Jacobian entries can produce")
    lines.append("    // NaN/Inf when a species is exactly zero and its rate law contains a")
    lines.append("    // pow(x, n<1) or 1/x factor that the formal derivative inherits — e.g.")
    lines.append("    // d/dx[x^(n-1)/x] at x=0. CVODE's FD Jacobian avoids this naturally,")
    lines.append("    // but we traded that robustness for speed. Clamping to 0 is safe: the")
    lines.append("    // step controller rejects the step if Newton fails to converge, and")
    lines.append("    // once species move off zero on subsequent substeps the entries become")
    lines.append("    // finite. This is the same remedy AMICI and libroadrunner use.")
    lines.append(f"    for (sunindextype k = 0; k < {nnz}; ++k) {{")
    lines.append("        if (!std::isfinite(data[k])) data[k] = 0.0;")
    lines.append("    }")
    lines.append("    return 0;")
    lines.append("}")
    return "\n".join(lines) + "\n"


# =========================================================================
# Code generators
# =========================================================================

def _get_constant_compartments(sbml: SBMLModel) -> list:
    """Compartments whose volumes are constant (not set by assignment rules)."""
    rule_vars = {r["variable_name"] for r in sbml.assignment_rules}
    param_names = {p["name"] for p in sbml.parameters}
    return [c for c in sbml.compartments
            if c["name"] not in rule_vars and c["name"] not in param_names]


def gen_enum_h(sbml: SBMLModel) -> str:
    const_comps = _get_constant_compartments(sbml)

    lines = [
        "#pragma once",
        "",
        "// Auto-generated by qsp_codegen.py from SBML — do not edit manually",
        "",
        "#include <utility>",
        "",
        "namespace CancerVCT{",
        "",
        "// QSP Species Enum (ODE state vector indices)",
        "enum QSPSpeciesEnum",
        "{",
    ]
    for sp in sbml.species:
        lines.append(f'SP_{_sanitize(sp["name"])},')
    lines.append(f"QSP_SPECIES_COUNT  // = {len(sbml.species)}")
    lines.append("};")
    lines.append("")
    lines.append("enum QSPNonSpeciesEnum { QSP_NON_SPECIES_COUNT = 0 };")
    lines.append("")
    lines.append("// QSP Parameter Enum (class parameters, SBML-native units)")
    lines.append("enum QSPParamEnum")
    lines.append("{")
    for p in sbml.parameters:
        lines.append(f'P_{_sanitize(p["name"])},')
    for c in const_comps:
        lines.append(f'P_{_sanitize(c["name"])},')
    lines.append(f"QSP_PARAM_COUNT")
    lines.append("};")
    lines.append("")

    # Unified file param enum: compartment ICs + species ICs + model parameters
    # This matches XML hierarchy: init_value.Compartment, init_value.Species,
    # init_value.Parameter
    lines.append("// QSP File Param Enum (XML file indices)")
    lines.append("// Order: Compartment ICs, Species ICs, Model Parameters")
    lines.append("enum QSPFileParamEnum")
    lines.append("{")
    lines.append("// --- Compartment initial values ---")
    for c in sbml.compartments:
        lines.append(f'QSP_{_sanitize(c["name"])},')
    lines.append("// --- Species initial values ---")
    for sp in sbml.species:
        lines.append(f'QSP_{_sanitize(sp["name"])},')
    lines.append("// --- Model parameters ---")
    for p in sbml.parameters:
        lines.append(f'QSP_{_sanitize(p["name"])},')
    lines.append(f"QSP_FILE_PARAM_COUNT")
    lines.append("};")
    lines.append("")
    lines.append("}  // namespace CancerVCT")
    return "\n".join(lines) + "\n"


def gen_ode_h(jac_nnz: int = 0) -> str:
    jac_decl = ''
    if jac_nnz > 0:
        jac_decl = (
            '\n    // --- Analytical Jacobian for KLU sparse linear solver ---\n'
            '    // Emitted by qsp_codegen.py via sympy symbolic differentiation of f.\n'
            '    // Sparsity pattern is static (derived from SBML reaction/rule structure)\n'
            '    // so it is safe to allocate the SUNSparseMatrix with this nnz once per\n'
            '    // solver instance. CSC format: col_ptrs has length neq+1, row_indices\n'
            '    // has length nnz.\n'
            f'    static constexpr sunindextype _jac_nnz = {jac_nnz};\n'
            '    static const sunindextype _jac_col_ptrs[];\n'
            '    static const sunindextype _jac_row_indices[];\n'
            '    static int jac(realtype t, N_Vector y, N_Vector fy,\n'
            '                   SUNMatrix J, void *user_data,\n'
            '                   N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);\n'
            '    // Hooks CVODEBase::setupCVODE consults to wire up the sparse\n'
            '    // linsol + jac callback (only active when built with USE_KLU).\n'
            '    sunindextype getJacobianNnz() const override { return _jac_nnz; }\n'
            '    CVLsJacFn getJacobianFn() const override { return &ODE_system::jac; }\n'
        )

    return ('''#ifndef __CancerVCT_ODE__
#define __CancerVCT_ODE__

// Auto-generated by qsp_codegen.py from SBML — do not edit manually

#include "../cvode/CVODEBase.h"
#include "QSPParam.h"
#include <sunmatrix/sunmatrix_sparse.h>
#include <string>
#include <vector>

namespace CancerVCT{

class ODE_system :
    public CVODEBase
{
public:
    static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
    static int g(realtype t, N_Vector y, realtype *gout, void *user_data);''' + jac_decl + '''
    static std::string getHeader();
    static void setup_class_parameters(QSPParam& param);
    static double get_class_param(unsigned int i);
    static void set_class_param(unsigned int i, double v);

    template<class Archive>
    static void classSerialize(Archive & ar, const unsigned int version);

    static bool use_steady_state;
    static bool use_resection;
    static double _QSP_weight;
private:
    static state_type _class_parameter;
public:
    ODE_system();
    ODE_system(const ODE_system& c);
    ~ODE_system();
    void setup_instance_tolerance(QSPParam& param);
    void setup_instance_variables(QSPParam& param);
    void eval_init_assignment(void);
    unsigned int get_num_variables(void)const { return _species_var.size(); };
    unsigned int get_num_params(void)const { return _class_parameter.size(); };
    // Evaluate a dynamic (assignment-rule) compartment's volume at the
    // current state, returned in its SBML-native unit (e.g. mL for V_T).
    // Returns the static compartment size for constant compartments.
    // Throws std::out_of_range for an unknown name.
    realtype get_compartment_volume(const std::string& name) const;
    // Evaluate a (non-compartment) SBML assignment rule's current value
    // in its native unit. Use this for derived quantities like
    // phi_collagen, C_total, etc. that downstream calibration targets
    // expect to read alongside species. Throws std::out_of_range for an
    // unknown name.
    realtype get_assignment_rule_value(const std::string& name) const;
    // Static enumeration of compartment names (for output column
    // ordering). Matches the order returned by get_compartment_volume's
    // case branches.
    static std::vector<std::string> getCompartmentNames();
    // Static enumeration of assignment-rule names whose values are
    // exposed via get_assignment_rule_value (excludes rules whose
    // target is also a compartment — those go through
    // get_compartment_volume).
    static std::vector<std::string> getAssignmentRuleNames();

protected:
    void setupVariables(void);
    void setupEvents(void);
    void initSolver(realtype t0);
    void update_y_other(void);
    void adjust_hybrid_variables(void);
    bool triggerComponentEvaluate(int i, realtype t, bool curr);
    bool eventEvaluate(int i);
    bool eventExecution(int i, bool delay, realtype& dt);
    realtype get_unit_conversion_species(int i) const;
    realtype get_unit_conversion_nspvar(int i) const;
    double getVarOriginalUnit(int i) const override;
    bool allow_negative(int i)const{ return false; };
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/);
};

template<class Archive>
inline void ODE_system::serialize(Archive & ar, const unsigned int){
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(CVODEBase);
}
template<class Archive>
void ODE_system::classSerialize(Archive & ar, const unsigned int){
    ar & BOOST_SERIALIZATION_NVP(_class_parameter);
}
inline double ODE_system::get_class_param(unsigned int i){
    return _class_parameter[i];
}
inline void ODE_system::set_class_param(unsigned int i, double v){
    _class_parameter[i] = v;
}

}  // namespace CancerVCT
#endif
''')


def gen_ode_cpp(sbml: SBMLModel, mapping: dict) -> str:
    lines = []
    lines.append('// Auto-generated by qsp_codegen.py from SBML — do not edit manually')
    lines.append('#include "ODE_system.h"')
    lines.append('#include "QSP_enum.h"')
    lines.append('')
    lines.append('#define SPVAR(x) NV_DATA_S(y)[x]')
    lines.append('#define NSPVAR(x) ptrOde->_nonspecies_var[x]')
    lines.append('#define PARAM(x) _class_parameter[x]')
    lines.append('#define PFILE(x) param.getVal(x)')
    lines.append('')
    lines.append('namespace CancerVCT{')
    lines.append('#define QSP_W ODE_system::_QSP_weight')
    lines.append('')
    lines.append('bool ODE_system::use_steady_state = false;')
    lines.append('bool ODE_system::use_resection = false;')
    lines.append('double ODE_system::_QSP_weight = 1.0;')
    lines.append('')

    # Constructor/destructor
    lines.append('ODE_system::ODE_system() :CVODEBase() {')
    lines.append('    setupVariables(); setupEvents(); setupCVODE(); update_y_other();')
    lines.append('}')
    lines.append('ODE_system::ODE_system(const ODE_system& c) { setupCVODE(); }')
    lines.append('ODE_system::~ODE_system() {}')
    lines.append('')
    lines.append('void ODE_system::initSolver(realtype t){')
    lines.append('    restore_y(); int flag;')
    lines.append('    flag = CVodeInit(_cvode_mem, f, t, _y);')
    lines.append('    check_flag(&flag, "CVodeInit", 1);')
    lines.append('    flag = CVodeRootInit(_cvode_mem, _nroot, g);')
    lines.append('    check_flag(&flag, "CVodeRootInit", 1);')
    lines.append('}')
    lines.append('')

    # Count total params (parameters + compartment-as-params)
    rule_vars = {r["variable_name"] for r in sbml.assignment_rules}
    param_names_set = {p["name"] for p in sbml.parameters}
    extra_comp_params = [c for c in sbml.compartments
                         if c["name"] not in rule_vars and c["name"] not in param_names_set]
    total_params = len(sbml.parameters) + len(extra_comp_params)

    lines.append(f'state_type ODE_system::_class_parameter = state_type({total_params}, 0);')
    lines.append('')

    # setup_class_parameters — convert to SI for dimensional consistency
    # SBML parameters use mixed time units (1/day, 1/s, 1/min, 1/hr)
    # SI factors ensure all parameters are in a single consistent system
    lines.append('void ODE_system::setup_class_parameters(QSPParam& param){')
    for i, p in enumerate(sbml.parameters):
        factor = sbml.get_si_factor(p["units"])
        san = _sanitize(p["name"])
        lines.append(f'    //{p["name"]}, index: {i}, units: {p["units"]}')
        lines.append(f'    _class_parameter[P_{san}] = PFILE(QSP_{san}) * {factor:.15g};')
    for j, c in enumerate(extra_comp_params):
        factor = sbml.get_si_factor(c["units"])
        san = _sanitize(c["name"])
        idx = len(sbml.parameters) + j
        lines.append(f'    //{c["name"]} (compartment), index: {idx}, units: {c["units"]}')
        lines.append(f'    _class_parameter[P_{san}] = PFILE(QSP_{san}) * {factor:.15g};')
    lines.append('}')
    lines.append('')

    # setupVariables
    lines.append('void ODE_system::setupVariables(void){')
    lines.append(f'    _species_var = std::vector<realtype>({len(sbml.species)}, 0);')
    lines.append(f'    _nonspecies_var = std::vector<realtype>(0, 0);')
    lines.append(f'    _species_other = std::vector<realtype>(0, 0);')
    lines.append('}')
    lines.append('')

    # Stubs (non-event). update_y_other is emitted later, after
    # eval_init_assignment, so it can share the assignment-rule writeback
    # block (see `emit_ar_species_writeback` in the init-assignment emit).
    lines.append('void ODE_system::adjust_hybrid_variables(void){ }')

    # --- Events --------------------------------------------------------
    # One root per event (single-comparison triggers only — see _parse_events).
    # Index i in _trigger_element_satisfied corresponds to event i.
    n_events = len(sbml.events)
    lines.append('void ODE_system::setupEvents(void){')
    lines.append(f'    _nevent = {n_events};')
    lines.append(f'    _nroot = {n_events};')
    lines.append(f'    _trigger_element_satisfied = std::vector<bool>({n_events}, false);')
    # TRIGGER_NON_INSTANT = standard root-finding event (triggers on sign
    # change of g()). Transient EQ/NEQ types are not currently emitted.
    lines.append(f'    _trigger_element_type = std::vector<EVENT_TRIGGER_ELEM_TYPE>('
                 f'{n_events}, TRIGGER_NON_INSTANT);')
    # Per-event initial trigger state from SBML trigger.initialValue (default
    # true). initialValue=true means the trigger is considered already
    # satisfied at t=0⁻, so events only fire on subsequent false→true
    # transitions. initialValue=false gives "fire at t=0 if condition holds".
    iv_list = ", ".join("true" if ev["trigger_initial_value"] else "false"
                        for ev in sbml.events)
    lines.append(f'    _event_triggered = std::vector<bool>{{{iv_list}}};')
    lines.append('}')
    lines.append('')

    _emit_event_code(lines, sbml)
    # get_unit_conversion_species — return SI factor per species
    # Used for: (1) tolerance scaling, (2) CSV output in SBML native units
    lines.append('realtype ODE_system::get_unit_conversion_species(int i) const {')
    lines.append('    static const realtype factors[] = {')
    for sp in sbml.species:
        factor = sbml.get_si_factor(sp["units"])
        lines.append(f'        {factor:.15g}, // {sp["name"]}')
    lines.append('    };')
    lines.append('    return factors[i];')
    lines.append('}')
    lines.append('realtype ODE_system::get_unit_conversion_nspvar(int) const { return 1.0; }')
    lines.append('')

    # getVarOriginalUnit — override to output SBML-native units
    # Amount species: _species_var / substance_factor → amount in SBML substance units
    # Concentration species: additionally divide by compartment volume to get concentration
    # in SBML concentration units (= substance / volume), matching SimBiology output.
    rule_vars_for_output = {r["variable_name"] for r in sbml.assignment_rules}
    comp_by_name_out = {c["name"]: c for c in sbml.compartments}

    # Identify dynamic compartments (those defined by assignment rules)
    dynamic_comps = set()
    for c in sbml.compartments:
        if c["name"] in rule_vars_for_output:
            dynamic_comps.add(c["name"])

    # Check if any concentration species are in dynamic compartments
    need_dynamic_vol = any(
        not sp.get("has_only_substance_units", False) and sp["compartment"] in dynamic_comps
        for sp in sbml.species
    )

    lines.append('double ODE_system::getVarOriginalUnit(int i) const {')
    lines.append('    // Base: amount in SI moles → SBML substance units')
    lines.append('    realtype v = (i < _neq) ? _species_var[i] : _species_other[i - _neq];')
    lines.append('    realtype sub_factor = get_unit_conversion_species(i);')
    lines.append('')

    # Identify assignment-rule species (their _species_var is stale; compute on-the-fly)
    sp_names_set = {sp["name"] for sp in sbml.species}
    ar_species_names = {r["variable_name"] for r in sbml.assignment_rules
                        if r["variable_name"] in sp_names_set}

    # Collect all assignment rules needed for output:
    # 1) Dynamic compartment volume rules (needed for concentration species)
    # 2) Assignment rules targeting species (needed to compute their values)
    # 3) All transitive dependencies of the above
    rules_by_name = {r["variable_name"]: r for r in sbml.assignment_rules}
    needed_rules = set()

    def collect_deps(vn):
        if vn in needed_rules:
            return
        needed_rules.add(vn)
        if vn in rules_by_name:
            expr = rules_by_name[vn]["expression"]
            for other_vn in rules_by_name:
                if other_vn == vn:
                    continue
                pattern = re.escape(other_vn)
                if re.search(r'(?<![a-zA-Z0-9_\.])' + pattern + r'(?![a-zA-Z0-9_\.])', expr):
                    collect_deps(other_vn)

    if need_dynamic_vol:
        for dc in dynamic_comps:
            collect_deps(dc)
    for ar_sp in ar_species_names:
        collect_deps(ar_sp)

    if needed_rules:
        # Build a mapping that uses _species_var[] instead of SPVAR()/NV_DATA_S(y)[]
        # since getVarOriginalUnit is a const member function (no y parameter)
        output_mapping = {}
        for sp in sbml.species:
            san = _sanitize(sp["name"])
            if sp.get("has_only_substance_units", False):
                output_mapping[sp["name"]] = f'_species_var[SP_{san}]'
            else:
                comp_name = sp["compartment"]
                comp_san = _sanitize(comp_name)
                if comp_name in rule_vars_for_output:
                    vol_expr = f'AUX_VAR_{comp_san}'
                else:
                    vol_expr = f'_class_parameter[P_{comp_san}]'
                output_mapping[sp["name"]] = f'(_species_var[SP_{san}] / {vol_expr})'
        for r in sbml.assignment_rules:
            vn = r["variable_name"]
            output_mapping[vn] = f'AUX_VAR_{_sanitize(vn)}'
        for p in sbml.parameters:
            if p["name"] not in rule_vars_for_output:
                output_mapping[p["name"]] = f'_class_parameter[P_{_sanitize(p["name"])}]'
        for c in sbml.compartments:
            if c["name"] in rule_vars_for_output:
                pass  # handled by assignment rules above
            elif c["name"] not in output_mapping:
                output_mapping[c["name"]] = f'_class_parameter[P_{_sanitize(c["name"])}]'

        # Order needed rules topologically
        needed_rule_list = [r for r in sbml.assignment_rules if r["variable_name"] in needed_rules]
        ordered_needed = order_rules(needed_rule_list,
                                     {r["variable_name"] for r in needed_rule_list},
                                     sbml=sbml)

        lines.append('    // Compute assignment rules from current state')
        for r in ordered_needed:
            vn = r["variable_name"]
            san = _sanitize(vn)
            cpp_expr = wrap_expression(r["expression"], output_mapping)
            lines.append(f'    realtype AUX_VAR_{san} = {cpp_expr};')
        lines.append('')

    # Generate switch for concentration species
    has_conc_species = any(not sp.get("has_only_substance_units", False) for sp in sbml.species)

    if has_conc_species:
        lines.append('    // For concentration species: divide by compartment volume')
        lines.append('    // Assignment-rule species: return computed AUX_VAR directly')
        lines.append('    switch (i) {')
        for idx, sp in enumerate(sbml.species):
            if sp.get("has_only_substance_units", False):
                # Check if this amount species is an assignment-rule target
                if sp["name"] in ar_species_names:
                    san = _sanitize(sp["name"])
                    # AUX_VAR is in SI moles, convert to SBML substance units
                    lines.append(f'    case {idx}: // {sp["name"]} (amount, assignment rule)')
                    lines.append(f'        return AUX_VAR_{san} / sub_factor;')
                continue
            comp_name = sp["compartment"]
            comp = comp_by_name_out[comp_name]
            comp_vol_factor = sbml.get_si_factor(comp["units"])
            comp_san = _sanitize(comp_name)

            if sp["name"] in ar_species_names:
                # Assignment-rule species: AUX_VAR is in SI concentration (mol/m² or mol/L)
                # Convert to SBML units: multiply by comp_vol_factor / sub_factor
                lines.append(f'    case {idx}: // {sp["name"]} (conc, assignment rule)')
                lines.append(f'        return AUX_VAR_{_sanitize(sp["name"])} * {comp_vol_factor:.15g} / sub_factor;')
            elif comp_name in dynamic_comps:
                # Dynamic compartment: volume is AUX_VAR computed above (in SI)
                lines.append(f'    case {idx}: // {sp["name"]} (conc, dynamic {comp_name})')
                lines.append(f'        return v * {comp_vol_factor:.15g} / (sub_factor * AUX_VAR_{comp_san});')
            else:
                # Fixed compartment: V_comp_SI = PARAM(P_comp)
                lines.append(f'    case {idx}: // {sp["name"]} (conc, fixed {comp_name})')
                lines.append(f'        return v * {comp_vol_factor:.15g} / (sub_factor * PARAM(P_{comp_san}));')
        lines.append('    default:')
        lines.append('        break;')
        lines.append('    }')
        lines.append('')

    lines.append('    // Amount species: just divide by substance factor')
    lines.append('    return v / sub_factor;')
    lines.append('}')
    lines.append('')

    lines.append('void ODE_system::setup_instance_tolerance(QSPParam&){')
    lines.append('    // Per-species absolute tolerance: abstol_base × SI_factor')
    lines.append('    // This ensures tolerance is meaningful in each species\' native scale.')
    lines.append('    // E.g., 1e-12 cells × (1/NA) = 1.66e-36 mol → controls to ~1e-12 cells.')
    lines.append('    realtype reltol = 1e-6;')
    lines.append('    realtype abstol_base = 1e-12;')
    lines.append('    N_Vector abstol = N_VNew_Serial(_neq, _sunctx);')
    lines.append('    for (int i = 0; i < _neq; i++) {')
    lines.append('        NV_DATA_S(abstol)[i] = abstol_base * get_unit_conversion_species(i);')
    lines.append('    }')
    lines.append('    int flag = CVodeSVtolerances(_cvode_mem, reltol, abstol);')
    lines.append('    check_flag(&flag, "CVodeSVtolerances", 1);')
    lines.append('    N_VDestroy(abstol);')
    lines.append('}')
    lines.append('')
    # Build compartment lookup: name → compartment dict
    comp_by_name = {c["name"]: c for c in sbml.compartments}

    # ---- get_compartment_volume ----------------------------------------
    # Returns the current volume of a named compartment in its SBML-native
    # unit. Dynamic (rule-driven) compartments recompute via the same AUX_VAR
    # chain used in f(). Constant compartments return _class_parameter[P_…]
    # converted back to native by dividing by its SI factor.
    rules_by_name_cv = {r["variable_name"]: r for r in sbml.assignment_rules}
    rule_vars_cv = set(rules_by_name_cv)
    mem_mapping_cv = {}
    for sp in sbml.species:
        san = _sanitize(sp["name"])
        if sp.get("has_only_substance_units", False):
            mem_mapping_cv[sp["name"]] = f'_species_var[SP_{san}]'
        else:
            cn = sp["compartment"]
            cs = _sanitize(cn)
            vol = f'AUX_VAR_{cs}' if cn in rule_vars_cv else f'PARAM(P_{cs})'
            mem_mapping_cv[sp["name"]] = f'(_species_var[SP_{san}] / {vol})'
    for r in sbml.assignment_rules:
        mem_mapping_cv[r["variable_name"]] = f'AUX_VAR_{_sanitize(r["variable_name"])}'
    for p in sbml.parameters:
        if p["name"] not in rule_vars_cv:
            mem_mapping_cv[p["name"]] = f'PARAM(P_{_sanitize(p["name"])})'
    for c in sbml.compartments:
        if c["name"] in rule_vars_cv:
            pass
        elif c["name"] not in mem_mapping_cv:
            mem_mapping_cv[c["name"]] = f'PARAM(P_{_sanitize(c["name"])})'

    lines.append('realtype ODE_system::get_compartment_volume(const std::string& name) const {')
    # Dynamic compartments: emit a case block per rule-driven compartment.
    dyn_comps = [c for c in sbml.compartments if c["name"] in rule_vars_cv]
    for c in dyn_comps:
        cname = c["name"]
        csan = _sanitize(cname)
        native_factor = sbml.get_si_factor(c["units"])
        # Collect rule chain needed for this compartment's rule.
        needed = _collect_rule_deps([cname], rules_by_name_cv, sbml=sbml)
        lines.append(f'    if (name == "{cname}") {{')
        for r in needed:
            rn = r["variable_name"]
            rs = _sanitize(rn)
            expr = wrap_expression(r["expression"], mem_mapping_cv)
            lines.append(f'        realtype AUX_VAR_{rs} = {expr};')
        # AUX_VAR value is in SI (m³/m² for compartments). Convert to native.
        lines.append(f'        return AUX_VAR_{csan} / {native_factor:.15g};')
        lines.append('    }')
    # Constant compartments: PARAM(P_name) stored in SI → divide by SI factor.
    for c in sbml.compartments:
        if c["name"] in rule_vars_cv:
            continue
        cname = c["name"]
        csan = _sanitize(cname)
        native_factor = sbml.get_si_factor(c["units"])
        lines.append(f'    if (name == "{cname}") return PARAM(P_{csan}) / {native_factor:.15g};')
    lines.append('    throw std::out_of_range("unknown compartment: " + name);')
    lines.append('}')
    lines.append('')

    # ---- get_assignment_rule_value -------------------------------------
    # Mirror of get_compartment_volume but for non-compartment assignment
    # rules — derived quantities like phi_collagen, C_total, K_T_Treg
    # that calibration-target functions expect to read alongside species
    # values. Compartment-targeted rules are excluded here (they are
    # served by get_compartment_volume).
    #
    # Dep collection runs in two passes. _collect_rule_deps walks the
    # SBML expression for *rule names*; that misses rules that get
    # introduced by mem_mapping substitution (e.g. a concentration
    # species "V_T.NO" expands to "(_species_var[SP_V_T_NO] / AUX_VAR_V_T)",
    # implicitly requiring V_T's rule even though "V_T" isn't a token in
    # the SBML source). A second pass scans each substituted expression
    # for AUX_VAR_X references and pulls in the matching rule until the
    # set is closed.
    comp_names_set = {c["name"] for c in sbml.compartments}
    param_by_name = {p["name"]: p for p in sbml.parameters}
    # Exclude rules whose LHS is a compartment (handled by get_compartment_volume)
    # OR a species (already emitted in the species dump via getVarOriginalUnit's
    # switch; emitting again here produces duplicate CSV columns and a stale
    # second copy because the standalone get_assignment_rule_value path doesn't
    # go through the same scaling fix-ups).
    sp_names_set_ar = {sp["name"] for sp in sbml.species}
    non_comp_rules = [
        r for r in sbml.assignment_rules
        if r["variable_name"] not in comp_names_set
        and r["variable_name"] not in sp_names_set_ar
    ]
    # Pre-build sanitized → variable_name lookup so we can resolve
    # AUX_VAR_<sanitized> tokens back to the original rule name.
    san_to_rule_name = {_sanitize(r["variable_name"]): r["variable_name"]
                        for r in sbml.assignment_rules}
    aux_var_re = re.compile(r'AUX_VAR_([a-zA-Z0-9_]+)')

    def _close_under_substitution(seed_names):
        """Return the rule dict list (topo-ordered) that includes every
        rule needed to evaluate the seeds, accounting for AUX_VAR
        references introduced by mem_mapping substitution."""
        needed_names = set()
        queue = list(seed_names)
        while queue:
            vn = queue.pop()
            if vn in needed_names or vn not in rules_by_name_cv:
                continue
            needed_names.add(vn)
            substituted = wrap_expression(
                rules_by_name_cv[vn]["expression"], mem_mapping_cv
            )
            for match in aux_var_re.finditer(substituted):
                ref_name = san_to_rule_name.get(match.group(1))
                if ref_name and ref_name != vn and ref_name not in needed_names:
                    queue.append(ref_name)
        rules_list = [rules_by_name_cv[n] for n in needed_names]
        return order_rules(rules_list, needed_names, sbml=sbml)

    lines.append(
        'realtype ODE_system::get_assignment_rule_value(const std::string& name) const {'
    )
    for r in non_comp_rules:
        rname = r["variable_name"]
        rsan = _sanitize(rname)
        # Resolve native-unit factor: rules typically write to a
        # parameter, in which case its declared unit is the target unit.
        # Fall back to dimensionless (factor 1) if we can't find one.
        target_param = param_by_name.get(rname)
        if target_param is not None and target_param.get("units"):
            native_factor = sbml.get_si_factor(target_param["units"])
        else:
            native_factor = 1.0
        needed = _close_under_substitution([rname])
        lines.append(f'    if (name == "{rname}") {{')
        for dep in needed:
            dn = dep["variable_name"]
            ds = _sanitize(dn)
            expr = wrap_expression(dep["expression"], mem_mapping_cv)
            lines.append(f'        realtype AUX_VAR_{ds} = {expr};')
        lines.append(f'        return AUX_VAR_{rsan} / {native_factor:.15g};')
        lines.append('    }')
    lines.append(
        '    throw std::out_of_range("unknown assignment rule: " + name);'
    )
    lines.append('}')
    lines.append('')

    # ---- name lists for output column ordering -----------------------
    lines.append('std::vector<std::string> ODE_system::getCompartmentNames() {')
    lines.append('    return {')
    for c in sbml.compartments:
        lines.append(f'        "{c["name"]}",')
    lines.append('    };')
    lines.append('}')
    lines.append('')

    lines.append('std::vector<std::string> ODE_system::getAssignmentRuleNames() {')
    lines.append('    return {')
    for r in non_comp_rules:
        lines.append(f'        "{r["variable_name"]}",')
    lines.append('    };')
    lines.append('}')
    lines.append('')

    # setup_instance_variables — convert to SI amounts
    # For AMOUNT species (hasOnlySubstanceUnits=true):
    #   PFILE = initialAmount in SBML units → × factor → SI moles
    # For CONCENTRATION species (hasOnlySubstanceUnits=false, initialConcentration):
    #   PFILE = initialConcentration → need amount = conc × comp_size
    #   → PFILE × PFILE(comp) × factor → SI moles
    # Note: species with initialAssignment rules get overridden by eval_init_assignment().
    lines.append('void ODE_system::setup_instance_variables(QSPParam& param){')
    for i, sp in enumerate(sbml.species):
        factor = sbml.get_si_factor(sp["units"])
        san = _sanitize(sp["name"])
        is_conc = sp.get("is_initial_concentration", False)
        hosu = sp.get("has_only_substance_units", False)
        lines.append(f'    //{sp["name"]}, index: {i}, units: {sp["units"]}')
        if is_conc and not hosu:
            # Concentration species: amount = concentration × compartment_size
            # PFILE(species) is initialConcentration in native units (e.g. nM)
            # PFILE(comp) is compartment size in native units (e.g. mm³)
            # Their product = amount in substance units (e.g. nM·mm³ = nmol)
            # factor converts substance units → SI (e.g. nmol → mol)
            comp_san = _sanitize(sp["compartment"])
            lines.append(f'    _species_var[SP_{san}] = PFILE(QSP_{san}) * PFILE(QSP_{comp_san}) * {factor:.15g};')
        else:
            lines.append(f'    _species_var[SP_{san}] = PFILE(QSP_{san}) * {factor:.15g};')
    lines.append('}')
    lines.append('')

    # getHeader
    header = ",".join(sp["name"] for sp in sbml.species)
    lines.append('std::string ODE_system::getHeader(){')
    lines.append(f'    return "{header}";')
    lines.append('}')
    lines.append('')

    # g() is emitted by _emit_event_code above; stub it only if there
    # are no events (keeps setupCVODE's CVodeRootInit(_nroot=0) happy).
    if not sbml.events:
        lines.append('int ODE_system::g(realtype, N_Vector, realtype*, void*){ return 0; }')
        lines.append('')

    # Build species lookup: name → species dict (for eval_init_assignment)
    sp_by_name = {sp["name"]: sp for sp in sbml.species}

    # eval_init_assignment — write to _species_var (not _y) so values persist
    # The caller should call restore_y() after this to sync to CVODE's _y.
    # Expression references to species use SPVAR() which reads NV_DATA_S(y),
    # so we first restore_y() to sync _species_var → _y, then write results
    # back to _species_var.
    lines.append('void ODE_system::eval_init_assignment(void){')
    lines.append('    // Sync _species_var → _y so SPVAR() reads see current values')
    lines.append('    restore_y();')
    lines.append('')
    for ia in sbml.initial_assignments:
        vn = ia["variable_name"]
        san = _sanitize(vn)
        cpp_expr = wrap_expression(ia["expression"], mapping)
        # For init, concentration-based species appearing in the expression
        # are already mapped to (SPVAR/V_comp) by classify_identifiers.
        # The expression evaluates to CONCENTRATION in SI (mol/m² or mol/m³).
        # To store AMOUNT: multiply by compartment volume/area.
        sp_info = sp_by_name.get(vn)
        if sp_info is not None and not sp_info.get("has_only_substance_units", False):
            comp = comp_by_name.get(sp_info["compartment"], {})
            comp_san = _sanitize(sp_info["compartment"])
            comp_name_str = sp_info["compartment"]
            rule_names_local = {r["variable_name"] for r in sbml.assignment_rules}
            if comp_name_str in rule_names_local:
                vol_str = f'AUX_VAR_{comp_san}'
            else:
                vol_str = f'PARAM(P_{comp_san})'
            lines.append(f'    // {vn} = ({ia["expression"][:60]}) * {comp_name_str}')
            lines.append(f'    _species_var[SP_{san}] = ({cpp_expr}) * {vol_str};')
            lines.append(f'    NV_DATA_S(_y)[SP_{san}] = _species_var[SP_{san}];')
        else:
            lines.append(f'    // {vn} = {ia["expression"][:80]}')
            lines.append(f'    _species_var[SP_{san}] = {cpp_expr};')
            lines.append(f'    NV_DATA_S(_y)[SP_{san}] = _species_var[SP_{san}];')
    # Evaluate assignment rules that target species, so _species_var is correct
    # before any CSV output at step 0.
    ar_species = [(r, sp_by_name[r["variable_name"]])
                  for r in sbml.assignment_rules
                  if r["variable_name"] in sp_by_name]

    def _emit_ar_species_writeback(out_lines: List[str], header_comment: str):
        """Emit rule-species evaluation + writeback block.

        Reads SBML species/parameter values via `_species_var[]` and
        `_class_parameter[]` (no `y` in scope — fits eval_init_assignment
        and update_y_other). For each assignment rule that targets a
        species, writes the rule-computed value to both `_species_var[]`
        and `NV_DATA_S(_y)[]` so post-step CSV dumps and downstream
        consumers see the correct value.

        Called from two sites:
          - eval_init_assignment: initialize rule species at t=0.
          - update_y_other: after each CVode step, re-evaluate rule
            species since save_y() copies _y -> _species_var and
            overwrites the rule value f() wrote (ydot=0 for rule species
            means _y never changes, so _species_var is pinned at the
            initial XML value — 0 for new rule species).
        """
        if not ar_species:
            return
        # Build mapping using _species_var[] instead of SPVAR() (no y in scope)
        ar_mapping = {}
        rule_names_ar = {r["variable_name"] for r in sbml.assignment_rules}
        for sp in sbml.species:
            san = _sanitize(sp["name"])
            if sp.get("has_only_substance_units", False):
                ar_mapping[sp["name"]] = f'_species_var[SP_{san}]'
            else:
                comp_name = sp["compartment"]
                comp_san = _sanitize(comp_name)
                if comp_name in rule_names_ar:
                    vol_expr = f'AUX_VAR_{comp_san}'
                else:
                    vol_expr = f'_class_parameter[P_{comp_san}]'
                ar_mapping[sp["name"]] = f'(_species_var[SP_{san}] / {vol_expr})'
        for r in sbml.assignment_rules:
            vn = r["variable_name"]
            ar_mapping[vn] = f'AUX_VAR_{_sanitize(vn)}'
        for p in sbml.parameters:
            if p["name"] not in rule_names_ar:
                ar_mapping[p["name"]] = f'_class_parameter[P_{_sanitize(p["name"])}]'
        for c in sbml.compartments:
            if c["name"] in rule_names_ar:
                pass
            elif c["name"] not in ar_mapping:
                ar_mapping[c["name"]] = f'_class_parameter[P_{_sanitize(c["name"])}]'

        # Topologically order all rules needed by the species assignment rules
        needed_rules = set()
        rules_by_name = {r["variable_name"]: r for r in sbml.assignment_rules}
        def collect_ar_deps(vn):
            if vn in needed_rules:
                return
            needed_rules.add(vn)
            if vn in rules_by_name:
                expr = rules_by_name[vn]["expression"]
                for other_vn in rules_by_name:
                    if other_vn == vn:
                        continue
                    pattern = re.escape(other_vn)
                    if re.search(r'(?<![a-zA-Z0-9_\.])' + pattern + r'(?![a-zA-Z0-9_\.])', expr):
                        collect_ar_deps(other_vn)
        for r, sp in ar_species:
            collect_ar_deps(r["variable_name"])
        needed_list = [r for r in sbml.assignment_rules if r["variable_name"] in needed_rules]
        ordered_needed = order_rules(needed_list, needed_rules, sbml=sbml)

        out_lines.append('')
        out_lines.append(f'    // {header_comment}')
        for r in ordered_needed:
            vn = r["variable_name"]
            san = _sanitize(vn)
            cpp_expr = wrap_expression(r["expression"], ar_mapping)
            out_lines.append(f'    realtype AUX_VAR_{san} = {cpp_expr};')
            # Write back to _species_var if this rule targets a species
            if vn in sp_by_name:
                sp = sp_by_name[vn]
                if sp.get("has_only_substance_units", False):
                    out_lines.append(f'    _species_var[SP_{san}] = AUX_VAR_{san};')
                else:
                    comp_name = sp["compartment"]
                    comp_san = _sanitize(comp_name)
                    if comp_name in rule_names_ar:
                        vol_expr = f'AUX_VAR_{comp_san}'
                    else:
                        vol_expr = f'_class_parameter[P_{comp_san}]'
                    out_lines.append(f'    _species_var[SP_{san}] = AUX_VAR_{san} * {vol_expr};')
                out_lines.append(f'    NV_DATA_S(_y)[SP_{san}] = _species_var[SP_{san}];')

    _emit_ar_species_writeback(
        lines,
        'Evaluate assignment-rule species so _species_var is correct at step 0',
    )

    lines.append('}')
    lines.append('')

    # ---------------------------------------------------------------
    # update_y_other: called by CVODEBase after every integration step.
    # save_y() copies _y -> _species_var, which wipes the rule-species
    # values that f() wrote during the step (ydot=0 for rule species, so
    # _y stays at its initial XML value — 0 for newly-introduced rule
    # species). Re-evaluate rules here so operator<< / getVarOriginalUnit
    # emit correct values in the CSV dump.
    # ---------------------------------------------------------------
    lines.append('void ODE_system::update_y_other(void){')
    _emit_ar_species_writeback(
        lines,
        'Re-evaluate assignment-rule species after save_y() overwrote them',
    )
    lines.append('}')
    lines.append('')

    # ---------------------------------------------------------------
    # RHS function f()
    # ---------------------------------------------------------------
    lines.append('int ODE_system::f(realtype t, N_Vector y, N_Vector ydot, void *user_data){')
    lines.append('')
    lines.append('    ODE_system* ptrOde = static_cast<ODE_system*>(user_data);')
    lines.append('')

    # 1) Assignment rules (topologically ordered)
    lines.append('    //Assignment rules:')
    lines.append('')
    rule_names = {r["variable_name"] for r in sbml.assignment_rules}
    ordered = order_rules(sbml.assignment_rules, rule_names, sbml=sbml)
    # Identify which assignment rules target species (need writeback to _species_var)
    sp_names_set = {sp["name"] for sp in sbml.species}
    sp_by_name = {sp["name"]: sp for sp in sbml.species}
    for r in ordered:
        vn = r["variable_name"]
        san = _sanitize(vn)
        cpp_expr = wrap_expression(r["expression"], mapping)
        lines.append(f'    realtype AUX_VAR_{san} = {cpp_expr};')
        # Write assignment-rule value back to _species_var so CSV output is correct
        if vn in sp_names_set:
            sp = sp_by_name[vn]
            if sp.get("has_only_substance_units", False):
                # Amount species: AUX_VAR is already in SI moles
                lines.append(f'    ptrOde->_species_var[SP_{san}] = AUX_VAR_{san};')
            else:
                # Concentration species: AUX_VAR is SI_moles/comp_vol,
                # multiply by compartment volume to get SI moles
                comp_name = sp["compartment"]
                comp_san = _sanitize(comp_name)
                if comp_name in rule_names:
                    vol_expr = f'AUX_VAR_{comp_san}'
                else:
                    vol_expr = f'PARAM(P_{comp_san})'
                lines.append(f'    ptrOde->_species_var[SP_{san}] = AUX_VAR_{san} * {vol_expr};')
        lines.append('')

    # 2) Reaction fluxes
    lines.append('    //Reaction fluxes:')
    lines.append('')
    for i, rxn in enumerate(sbml.reactions):
        flux_idx = i + 1
        cpp_rate = wrap_expression(rxn["rate_law"], mapping)
        lines.append(f'    realtype ReactionFlux{flux_idx} = {cpp_rate};')
        lines.append('')

    # 3) ydot assembly
    lines.append('    //ODE right-hand side:')
    lines.append('')
    sp_names = {sp["name"] for sp in sbml.species}
    stoich = build_stoichiometry(sbml.reactions, sp_names)

    for sp in sbml.species:
        spn = sp["name"]
        san = _sanitize(spn)
        entries = stoich.get(spn, [])

        if not entries:
            lines.append(f'    NV_DATA_S(ydot)[SP_{san}] = 0.0;')
            continue

        flux_terms = []
        for fidx, sign in entries:
            flux_terms.append(f"{'+ ' if sign > 0 else '- '}ReactionFlux{fidx}")
        flux_expr = " ".join(flux_terms)
        if flux_expr.startswith("+ "):
            flux_expr = flux_expr[2:]

        comp_name = needs_volume_scaling(sp, sbml)
        if comp_name is not None:
            comp_san = _sanitize(comp_name)
            if comp_name in rule_names:
                vol = f"AUX_VAR_{comp_san}"
            else:
                vol = f"PARAM(P_{comp_san})"
            lines.append(f'    NV_DATA_S(ydot)[SP_{san}] = 1/{vol}*({flux_expr});')
        else:
            lines.append(f'    NV_DATA_S(ydot)[SP_{san}] = {flux_expr};')

    lines.append('')
    lines.append('    return 0;')
    lines.append('}')
    lines.append('')

    # Analytical Jacobian (KLU-ready sparse CSC format). Generated via sympy
    # symbolic differentiation of the same rate/rule expressions that
    # produced f() above. The CVODEBase setupCVODE path will attach this
    # via CVodeSetJacFn when built with KLU support; dense builds ignore it.
    if gen_ode_cpp._jacobian_info is not None:
        lines.append(gen_jacobian_cpp(sbml, gen_ode_cpp._jacobian_info))
    lines.append('}  // namespace CancerVCT')
    return "\n".join(lines) + "\n"


# Set by main() before gen_ode_cpp runs, so gen_ode_h (called first) can
# embed the nnz and gen_ode_cpp can reuse the sympy work.
gen_ode_cpp._jacobian_info = None


def gen_qsp_param_h() -> str:
    return '''#ifndef __CancerVCT_QSPParam__
#define __CancerVCT_QSPParam__

// Auto-generated by qsp_codegen.py from SBML — do not edit manually

#include "../../core/ParamBase.h"

namespace CancerVCT{

class QSPParam : public SP_QSP_IO::ParamBase
{
public:
    QSPParam();
    ~QSPParam(){};
    double getVal(int n) const;
    void printParam(void) const;
    bool readParamsFromXml(std::string inFileName) override {
        _readParameters(inFileName);
        return true;
    }
    void setupParam() override {}
    void processInternalParams() override {}
private:
    std::vector<double> _param;
    void _readParameters(const std::string& filename);
    static const char* _xml_paths[];
};

}  // namespace CancerVCT
#endif
'''


def gen_qsp_param_cpp(sbml: SBMLModel) -> str:
    total = len(sbml.compartments) + len(sbml.species) + len(sbml.parameters)

    lines = []
    lines.append('// Auto-generated by qsp_codegen.py from SBML — do not edit manually')
    lines.append('#include "QSPParam.h"')
    lines.append('#include "QSP_enum.h"')
    lines.append('#include <boost/property_tree/xml_parser.hpp>')
    lines.append('#include <iostream>')
    lines.append('')
    lines.append('namespace CancerVCT{')
    lines.append('')

    # XML paths: Compartment ICs, Species ICs, Model Parameters
    # Must match QSPFileParamEnum order exactly
    lines.append('const char* QSPParam::_xml_paths[] = {')
    lines.append('    // Compartment initial values')
    for c in sbml.compartments:
        san = _sanitize(c["name"])
        lines.append(f'    "Param.QSP.init_value.Compartment.{san}",')
    lines.append('    // Species initial values')
    for sp in sbml.species:
        san = _sanitize(sp["name"])
        lines.append(f'    "Param.QSP.init_value.Species.{san}",')
    lines.append('    // Model parameters')
    for p in sbml.parameters:
        san = _sanitize(p["name"])
        lines.append(f'    "Param.QSP.init_value.Parameter.{san}",')
    lines.append('};')
    lines.append('')
    lines.append(f'QSPParam::QSPParam() :_param({total}, 0) {{}}')
    lines.append('')
    lines.append('double QSPParam::getVal(int n) const { return _param[n]; }')
    lines.append('')
    lines.append('void QSPParam::_readParameters(const std::string& filename){')
    lines.append('    namespace pt = boost::property_tree;')
    lines.append('    pt::ptree tree;')
    lines.append('    pt::read_xml(filename, tree, pt::xml_parser::trim_whitespace);')
    lines.append(f'    for (int i = 0; i < {total}; i++){{')
    lines.append('        try { _param[i] = tree.get<double>(_xml_paths[i]); }')
    lines.append('        catch (const pt::ptree_bad_path&) {')
    lines.append('            std::cerr << "WARNING: QSP param not found: "')
    lines.append('                      << _xml_paths[i] << std::endl;')
    lines.append('        }')
    lines.append('    }')
    lines.append('}')
    lines.append('')
    lines.append('void QSPParam::printParam(void) const {')
    lines.append(f'    for (int i = 0; i < {total}; i++)')
    lines.append('        std::cout << _xml_paths[i] << " = " << _param[i] << std::endl;')
    lines.append('}')
    lines.append('')
    lines.append('}  // namespace CancerVCT')
    return "\n".join(lines) + "\n"


def gen_xml_snippet(sbml: SBMLModel) -> str:
    """Generate XML <QSP> section for param_all_test.xml."""
    lines = ['<QSP>']
    lines.append('  <simulation>')
    lines.append('    <weight_qsp>0.8</weight_qsp>')
    lines.append('    <t_steadystate>5000</t_steadystate>')
    lines.append('    <use_resection>0</use_resection>')
    lines.append('    <t_resection>1000</t_resection>')
    lines.append('    <presimulation_diam_frac>1.00</presimulation_diam_frac>')
    lines.append('    <start>0</start>')
    lines.append('    <step>1</step>')
    lines.append('    <n_step>360</n_step>')
    lines.append('    <tol_rel>1e-09</tol_rel>')
    lines.append('    <tol_abs>1e-12</tol_abs>')
    lines.append('  </simulation>')
    lines.append('  <init_value>')

    # Compartments
    lines.append('    <Compartment>')
    for c in sbml.compartments:
        san = _sanitize(c["name"])
        lines.append(f'      <{san}>{c["size"]}</{san}>')
    lines.append('    </Compartment>')

    # Species
    lines.append('    <Species>')
    for sp in sbml.species:
        san = _sanitize(sp["name"])
        lines.append(f'      <{san}>{sp["initial_value"]}</{san}>')
    lines.append('    </Species>')

    # Parameters
    lines.append('    <Parameter>')
    for p in sbml.parameters:
        san = _sanitize(p["name"])
        lines.append(f'      <{san}>{p["value"]}</{san}>')
    lines.append('    </Parameter>')

    lines.append('  </init_value>')
    lines.append('</QSP>')
    return "\n".join(lines)


# =========================================================================
# Main
# =========================================================================

def generate(sbml_path: str, out_dir: str) -> Dict[str, str]:
    """Run codegen for the given SBML, writing all outputs under ``out_dir``.

    Returns a dict of ``{filename: content}`` so callers (e.g. tests) can
    inspect generated code without re-reading from disk.
    """
    print(f"Parsing SBML: {sbml_path}")
    sbml = SBMLModel(sbml_path)
    print(f"  Species:          {len(sbml.species)}")
    print(f"  Compartments:     {len(sbml.compartments)}")
    print(f"  Reactions:        {len(sbml.reactions)}")
    print(f"  Parameters:       {len(sbml.parameters)}")
    print(f"  Assignment rules: {len(sbml.assignment_rules)}")
    print(f"  Initial assigns:  {len(sbml.initial_assignments)}")
    print(f"  Unit definitions: {len(sbml.unit_defs)}")

    mapping = classify_identifiers(sbml)
    print(f"  Identifier mappings: {len(mapping)}")

    try:
        import sympy  # noqa: F401
        print("\nDeriving analytical Jacobian (sympy + CSE)...")
        jac_info = compute_jacobian(sbml, mapping)
        print(f"  nnz = {jac_info['nnz']} "
              f"({100.0 * jac_info['nnz'] / (len(sbml.species) ** 2):.2f}% density)")
        print(f"  CSE subexpressions: {len(jac_info['cse_items'])}")
        gen_ode_cpp._jacobian_info = jac_info
        jac_nnz = jac_info["nnz"]
    except ImportError:
        print("\nsympy not installed — skipping analytical Jacobian emission. "
              "Install with `uv pip install sympy` to enable KLU sparse solver.")
        gen_ode_cpp._jacobian_info = None
        jac_nnz = 0

    print("\nGenerating C++ files...")
    files = {
        "QSP_enum.h": gen_enum_h(sbml),
        "ODE_system.h": gen_ode_h(jac_nnz=jac_nnz),
        "ODE_system.cpp": gen_ode_cpp(sbml, mapping),
        "QSPParam.h": gen_qsp_param_h(),
        "QSPParam.cpp": gen_qsp_param_cpp(sbml),
        "qsp_params_xml_snippet.xml": gen_xml_snippet(sbml),
    }

    os.makedirs(out_dir, exist_ok=True)
    for fname, content in files.items():
        path = os.path.join(out_dir, fname)
        with open(path, "w") as f:
            f.write(content)
        print(f"  {fname}: {content.count(chr(10))} lines")

    print(f"\nDone! Output: {out_dir}")
    return files


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="SBML → C++ CVODE ODE code generator.",
    )
    parser.add_argument(
        "--sbml",
        required=True,
        help="Path to SBML Level 2 v4 file (e.g. PDAC_model.sbml).",
    )
    parser.add_argument(
        "--out-dir",
        required=True,
        help="Directory where generated C++ sources are written.",
    )
    args = parser.parse_args(argv)
    generate(args.sbml, args.out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
