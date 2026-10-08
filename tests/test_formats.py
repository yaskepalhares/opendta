import math

import pytest

from opendta.core import missing as M
from opendta.core.formats import format_value, number_to_macro


# Valores observados nos logs do Stata 14 (compat/expected/0001 e 0002).
@pytest.mark.parametrize("value, expected", [
    (1 / 3, ".33333333"),
    (-1 / 3, "-.33333333"),
    (math.pi, "3.1415927"),
    (-math.pi, "-3.1415927"),
    (123456789, "1.235e+08"),
    (-123456789, "-1.235e+08"),
    (1234567890, "1.235e+09"),
    (0.00000123, "1.230e-06"),
    (123.456, "123.456"),
    (0.5, ".5"),
    (-0.5, "-.5"),
    (1e10, "1.000e+10"),
    (1e-10, "1.000e-10"),
    (12345.678901, "12345.679"),
    (99999.99999, "100000"),
    (100000, "100000"),
    (1000000000, "1.000e+09"),
    (2, "2"),
])
def test_display_default(value, expected):
    from opendta.core.formats import DISPLAY_NUMERIC
    assert format_value(value, DISPLAY_NUMERIC, pad=False) == expected


@pytest.mark.parametrize("value, fmt, expected", [
    (1 / 3, "%10.0g", " .33333333"),
    (1 / 3, "%5.0g", " .333"),
    (123456789, "%8.0g", " 1.2e+08"),
])
def test_general_widths(value, fmt, expected):
    assert format_value(value, fmt) == expected


def test_decimal_comma():
    from opendta.core.formats import set_decimal_comma
    set_decimal_comma(True)
    try:
        assert format_value(1 / 3, "%10.0g", pad=False) == ",33333333"
        assert format_value(1234567.891, "%12.2fc") == "1.234.567,89"
        assert format_value(1 / 3, "%9.0g", pad=False, dp=False) == ".3333333"
    finally:
        set_decimal_comma(False)


@pytest.mark.parametrize("value, fmt, expected", [
    (math.pi, "%9.2f", "     3.14"),
    (math.pi, "%-9.2f", "3.14     "),
    (math.pi, "%09.2f", "000003.14"),
    (1234567.891, "%12.2fc", "1,234,567.89"),
    (1234.5, "%10.3e", " 1.235e+03"),      # meio para cima, como no Stata
    ("ab", "%-5s", "ab   "),
    ("ab", "%5s", "   ab"),
    (0, "%td", "01jan1960"),
])
def test_explicit(value, fmt, expected):
    assert format_value(value, fmt) == expected


def test_missing_codes():
    assert M.SYSMISS == 2.0 ** 1023
    assert M.EXTENDED["a"] > M.SYSMISS
    assert M.EXTENDED["z"] > M.EXTENDED["a"]
    assert format_value(M.EXTENDED["b"], "%9.2f") == "       .b"
    assert M.MAXDOUBLE < M.SYSMISS


def test_number_to_macro():
    assert number_to_macro(5.0) == "5"
    assert number_to_macro(1 / 3) == ".3333333333333333"
    assert number_to_macro(1 / 7) == ".1428571428571428"
    assert number_to_macro(1e20) == "1.00000000000e+20"
    assert number_to_macro(-0.25) == "-.25"
    assert number_to_macro(M.SYSMISS) == "."
