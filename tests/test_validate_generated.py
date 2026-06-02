"""Tests for validate_generated_cpp — the codegen-time guard against
use-before-definition of AUX_VAR_* temporaries in generated C++.

Regression context: a concentration rule-species in a dynamic-volume
compartment emitted `AUX_VAR_V_T` undefined in the init/update_y_other blocks,
surfacing only as an opaque C++ compiler error on machine-generated lines. This
guard turns that class of dependency-ordering bug into a clear codegen-time error.
"""
import pytest

from qsp_codegen.codegen import validate_generated_cpp as validate


def test_clean_passes():
    cpp = (
        "namespace X {\n"
        "void f(){\n"
        "    realtype AUX_VAR_V_T = 1.0;\n"
        "    realtype c = AUX_VAR_V_T * 2.0;\n"
        "}\n"
        "}\n"
    )
    validate(cpp)  # no raise


def test_redeclaration_across_functions_ok():
    # Same temp legitimately redeclared in separate function scopes.
    cpp = (
        "namespace X {\n"
        "void f(){ realtype AUX_VAR_V_T = 1.0; realtype z = AUX_VAR_V_T; }\n"
        "void g(){ realtype AUX_VAR_V_T = 2.0; realtype w = AUX_VAR_V_T; }\n"
        "}\n"
    )
    validate(cpp)  # no raise


def test_use_before_decl_multiline_raises():
    cpp = (
        "namespace X {\n"
        "void f(){\n"
        "    realtype y = AUX_VAR_V_T * 2.0;\n"   # used here
        "    realtype AUX_VAR_V_T = 3.0;\n"       # declared after
        "}\n"
        "}\n"
    )
    with pytest.raises(ValueError, match="AUX_VAR_V_T"):
        validate(cpp)


def test_use_in_function_without_decl_raises():
    # The exact shape of the real bug: declared in one function, used in
    # another that never declares it.
    cpp = (
        "namespace X {\n"
        "void f(){ realtype AUX_VAR_V_T = 1.0; realtype a = AUX_VAR_V_T; }\n"
        "void g(){ realtype b = AUX_VAR_V_T * 2.0; }\n"  # no decl in g
        "}\n"
    )
    with pytest.raises(ValueError, match="line"):
        validate(cpp)
