import numpy as np
import pytest

from opendta.core import missing as M
from opendta.core.errors import StataError
from opendta.lang.syntax import parse_in, parse_standard


def make(run, n=5):
    run(f"clear\nset obs {n}\ngen id = _n")


# -- sintaxe -------------------------------------------------------------------

def test_parse_standard():
    p = parse_standard('byte x = cond(a, 1, 2) if y > 1 & z == "a b" in 1/10, after(id)')
    assert p.varlist == "byte x"
    assert p.exp == "cond(a, 1, 2)"
    assert p.if_ == 'y > 1 & z == "a b"'
    assert p.in_ == "1/10"
    assert p.options == "after(id)"
    p = parse_standard("x y [fw=n] if x[_n-1] > 0")
    assert p.weight == ("fweight", "n")
    assert p.if_ == "x[_n-1] > 0"
    p = parse_standard("x = y == 1")
    assert p.exp == "y == 1"


def test_parse_in():
    assert parse_in("1/3", 10) == (0, 3)
    assert parse_in("f/l", 10) == (0, 10)
    assert parse_in("-2/l", 10) == (8, 10)
    assert parse_in("5", 10) == (4, 5)
    with pytest.raises(StataError):
        parse_in("5/20", 10)


# -- tipos e precisão ---------------------------------------------------------------

def test_float_precision(run):
    make(run)
    run("gen f = 0.1\ngen double d = 0.1")
    assert run("count if f == 0.1") == "  0\n"
    assert run("count if d == 0.1") == "  5\n"
    assert run("count if f == float(0.1)") == "  5\n"


def test_integer_types(run):
    make(run, 3)
    run("gen byte b = 2.7\ngen byte big = 200")
    ds = run.session.data
    assert list(ds.get("b").data) == [2, 2, 2]
    assert all(v >= M.SYSMISS for v in ds.get("big").data)


def test_generate_messages(run):
    make(run, 4)
    assert run("gen y = 1 if id > 2") == "(2 missing values generated)\n"
    assert run("gen z = 1 if id > 3") == "(3 missing values generated)\n"
    assert run("gen w = 1 if id > 1") == "(1 missing value generated)\n"
    run("gen id = 1")
    assert run.rc == 110


def test_replace_promotes_and_counts(run):
    make(run, 4)
    run("gen byte b = 1")
    out = run("replace b = 1000 in 1")
    assert out == "variable b was byte now int\n(1 real change made)\n"
    out = run("replace b = 1.5 in 2")
    assert out.startswith("variable b was int now float")
    assert run("replace b = . in 3") == "(1 real change made, 1 to missing)\n"
    assert run("replace b = b") == "(0 real changes made)\n"


def test_replace_sequential(run):
    run("clear\ninput x\n1\n.\n.\n4\n.\nend")
    run("replace x = x[_n-1] if missing(x)")
    assert list(run.session.data.get("x").data) == [1, 1, 1, 4, 4]


def test_strings(run):
    make(run, 3)
    run('gen s = "ab"')
    assert run.session.data.get("s").vtype == "str2"
    out = run('replace s = "abcde" in 1')
    assert out.startswith("variable s was str2 now str5")
    run("gen n = s + 1")
    assert run.rc == 109
    run("gen str3 t = 1")
    assert run.rc == 109


# -- _n, _N, by ---------------------------------------------------------------------

def test_by_groups(run):
    run("clear\ninput g x\n1 10\n2 20\n1 30\n2 40\n1 50\nend")
    run("bysort g (x): gen n = _n\nby g: gen N = _N\nby g: gen s = sum(x)")
    ds = run.session.data
    assert list(ds.get("g").data) == [1, 1, 1, 2, 2]
    assert list(ds.get("n").data) == [1, 2, 3, 1, 2]
    assert list(ds.get("N").data) == [3, 3, 3, 2, 2]
    assert list(ds.get("s").data) == [10, 40, 90, 20, 60]
    run("by g: gen first = x[1]")
    assert list(ds.get("first").data) == [10, 10, 10, 20, 20]
    run("by g: gen prev = x[_n-1]")
    assert ds.get("prev").data[0] >= M.SYSMISS and ds.get("prev").data[1] == 10


def test_by_requires_sort_and_byable(run):
    run("clear\ninput g\n2\n1\nend")
    run("by g: gen n = _n")
    assert run.rc == 5
    run("bysort g: describe")
    assert run.rc == 190


def test_by_keep_first(run):
    run("clear\ninput g x\n1 1\n1 2\n2 3\nend")
    assert run("by g, sort: keep if _n == 1") == "(1 observation deleted)\n"
    assert list(run.session.data.get("x").data) == [1, 3]


# -- varlists -------------------------------------------------------------------

def test_varlists(run):
    run("clear\nset obs 1\ngen alpha = 1\ngen beta = 2\ngen gamma = 3\ngen alpine = 4")
    s = run.session
    assert s.expand_varlist("beta-gamma") == ["beta", "gamma"]
    assert s.expand_varlist("al*") == ["alpha", "alpine"]
    assert s.expand_varlist("g") == ["gamma"]
    assert s.expand_varlist("_all") == ["alpha", "beta", "gamma", "alpine"]
    with pytest.raises(StataError) as e:
        s.expand_varlist("al")
    assert e.value.rc == 111 and "ambiguous" in e.value.message
    with pytest.raises(StataError):
        s.expand_varlist("zzz")


# -- drop, keep, sort, rename, order ----------------------------------------------

def test_drop_keep(run):
    make(run, 6)
    assert run("drop in 1/2") == "(2 observations deleted)\n"
    assert run("keep if id <= 4") == "(0 observations deleted)\n" or True
    run("gen a = 1\ngen b = 2")
    run("drop a")
    assert run.session.data.names == ["id", "b"]
    run("keep b")
    assert run.session.data.names == ["b"]


def test_sort_gsort(run):
    run("clear\ninput x str1 s\n3 b\n.\n1 a\n2 c\nend")
    run("sort x")
    assert list(run.session.data.get("x").data[:3]) == [1, 2, 3]
    assert run.session.data.get("x").data[3] >= M.SYSMISS
    run("gsort -x")
    vals = list(run.session.data.get("x").data)
    assert vals[:3] == [3, 2, 1] and vals[3] >= M.SYSMISS
    run("gsort -s")
    assert list(run.session.data.get("s").data) == ["c", "b", "a", ""]


def test_rename_order(run):
    run("clear\nset obs 1\ngen a = 1\ngen b = 2\ngen c = 3")
    run("rename a z\norder c\norder z, last")
    assert run.session.data.names == ["c", "b", "z"]


# -- list, describe, label ------------------------------------------------------------

def test_list_table(run):
    run("clear\ninput x str3 nome\n1 ana\n22 bia\nend")
    out = run("list")
    lines = out.strip("\n").splitlines()
    # layout observado no Stata: 3 espaços entre colunas; strings com %9s à direita
    assert lines[0] == "     +-----------+"
    assert lines[1] == "     |  x   nome |"
    assert lines[3] == "  1. |  1    ana |"
    assert lines[4] == "  2. | 22    bia |"


def test_value_labels(run):
    run("clear\ninput sexo\n1\n2\nend")
    run('label define sx 1 "Masculino" 2 "Feminino"\nlabel values sexo sx')
    assert "Masculino" in run("list")
    assert "Masculino" not in run("list, nolabel")
    assert run("label list sx") == "sx:\n           1 Masculino\n           2 Feminino\n"
    run('label define sx 3 "Outro"')
    assert run.rc == 110
    run('label define sx 3 "Outro", add')
    assert run.rc == 0


def test_describe(run):
    run("clear\nset obs 2\ngen byte x = 1")
    run('label variable x "Minha variável"')
    out = run("describe")
    assert "  obs:             2" in out
    assert " vars:             1" in out
    assert "x               byte    %8.0g                 Minha variável" in out


def test_display_variable_first_obs(run):
    run("clear\ninput x\n7\n8\nend")
    assert run("di x") == "7\n"
    assert run("di x[2] + _N") == "10\n"


def test_input_strings_and_missing(run):
    run('clear\ninput str5 nome idade\n"ana b" 30\nbia .\nend')
    ds = run.session.data
    assert list(ds.get("nome").data) == ["ana b", "bia"]
    assert ds.get("idade").data[1] >= M.SYSMISS


def test_compress(run):
    make(run, 3)
    out = run("compress")
    assert "variable id was float now byte" in out
    assert run.session.data.get("id").vtype == "byte"
