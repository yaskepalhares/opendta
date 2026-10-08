import math

import pytest

from opendta.core import missing as M
from opendta.core.formats import format_value, number_to_macro


@pytest.mark.parametrize("value, expected", [
    (1 / 3, ".33333333"),
    (math.pi, "3.1415927"),
    (123456789, "123456789"),
    (1234567890, "1.235e+09"),
    (0.5, ".5"),
    (123.456, "123.456"),
    (2, "2"),
])
def test_general_default(value, expected):
    assert format_value(value, pad=False, sign_outside_width=True) == expected


@pytest.mark.parametrize("value, fmt, expected", [
    (math.pi, "%9.2f", "     3.14"),
    (math.pi, "%-9.2f", "3.14     "),
    (math.pi, "%09.2f", "000003.14"),
    (1234567.891, "%12.2fc", "1,234,567.89"),
    (1234.5, "%10.3e", " 1.234e+03"),
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
    assert number_to_macro(-0.25) == "-.25"
    assert number_to_macro(M.SYSMISS) == "."
