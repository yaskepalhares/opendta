import math

import pytest

from opendta.core import missing as M
from opendta.core.errors import StataError


@pytest.mark.parametrize("expr, expected", [
    ("1+2*3", 7),
    ("(1+2)*3", 9),
    ("2^3^2", 64),        # ^ associativa à esquerda
    ("-2^2", -4),         # ^ antes do menos unário
    ("10/4", 2.5),
    ("!0", 1),
    ("!5", 0),
    ("3>2 & 2>1", 1),
    ("0 | 0", 0),
    ("1 == 1", 1),
    ("1 != 1", 0),
    ("1 ~= 2", 1),
    ("0^0", 1),
])
def test_arithmetic(run, expr, expected):
    assert run.eval(expr) == expected


def test_missing_rules(run):
    assert M.is_missing(run.eval("1/0"))
    assert M.is_missing(run.eval(". + 1"))
    assert M.is_missing(run.eval("sqrt(-1)"))
    assert M.is_missing(run.eval("(-8)^(1/3)"))
    assert run.eval(". > 1e300") == 1
    assert run.eval(".a > .") == 1
    assert run.eval(".z > .a") == 1
    assert run.eval("!.") == 0           # missing é verdadeiro
    assert run.eval(". == .") == 1
    assert run.eval(".a == .") == 0


def test_strings(run):
    assert run.eval('"a" + "b"') == "ab"
    assert run.eval('"ab" * 3') == "ababab"
    assert run.eval('"abc" < "abd"') == 1
    with pytest.raises(StataError) as e:
        run.eval('"a" + 1')
    assert e.value.rc == 109


def test_functions(run):
    assert run.eval("abs(-3)") == 3
    assert run.eval("round(2.5)") == 3
    assert run.eval("round(1.234, .01)") == pytest.approx(1.23)
    assert run.eval("mod(7, 3)") == 1
    assert run.eval("mod(-7, 3)") == 2
    assert run.eval("max(1, ., 3)") == 3
    assert run.eval("min(., 2)") == 2
    assert run.eval("cond(1, 10, 20)") == 10
    assert run.eval("inlist(2, 1, 2, 3)") == 1
    assert run.eval("inrange(5, 1, 10)") == 1
    assert run.eval("missing(.)") == 1
    assert run.eval('missing("")') == 1
    assert run.eval('substr("abcdef", 2, 3)') == "bcd"
    assert run.eval('substr("abcdef", -2, .)') == "ef"
    assert run.eval('strpos("abc", "c")') == 3
    assert run.eval('upper("abc")') == "ABC"
    assert run.eval('word("a b c", 2)') == "b"
    assert run.eval('wordcount("a b  c")') == 3
    assert run.eval('real("1.5")') == 1.5
    assert run.eval('string(1/4)') == ".25"
    assert run.eval('subinstr("aaa", "a", "b", 2)') == "bba"
    assert run.eval("normal(0)") == pytest.approx(0.5)
    assert run.eval("invnormal(0.975)") == pytest.approx(1.959964, abs=1e-6)
    assert run.eval("_pi") == math.pi
    assert run.eval("mdy(1, 1, 1960)") == 0
    assert run.eval("year(mdy(10, 8, 2026))") == 2026
    assert run.eval('date("08oct2026", "DMY")') == run.eval("mdy(10, 8, 2026)")
    assert run.eval('regexm("abc123", "[0-9]+")') == 1
    # comportamentos observados no Stata (compat/expected/0003)
    assert run.eval('proper("joão da silva")') == "JoãO Da Silva"
    assert run.eval('abbrev("variavel_muito_longa", 10)') == "variavel~a"
    assert run.eval("round(-2.5)") == -2 and run.eval("round(-1.5)") == -1


def test_unknown_function(run):
    with pytest.raises(StataError) as e:
        run.eval("nope(1)")
    assert e.value.rc == 133
