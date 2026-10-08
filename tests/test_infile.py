"""infile (livre e com dicionário) e infix."""

import pytest

from opendta.core import missing as M
from opendta.io.fixed import parse_infix_spec, tokenize_free


@pytest.fixture
def here(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_tokenize():
    assert tokenize_free('a "b c", 3\n4,,5 .a') == ["a", "b c", "3", "4", "", "5", ".a"]


def test_infile_free(run, here):
    (here / "f.raw").write_text('Ana 30 1.5\n"Bia Lima" . 2\nCaio, 41, abc\nDavi 7\n')
    out = run("infile str10 nome idade renda using f")
    assert "'abc' cannot be read as a number for renda[3]" in out
    assert "(eof not at end of obs)" in out and out.endswith("(4 observations read)\n")
    d = run.session.data
    assert [v.vtype for v in d.vars] == ["str10", "float", "float"]
    assert list(d.get("nome").data) == ["Ana", "Bia Lima", "Caio", "Davi"]
    assert M.is_missing(d.get("renda").data[2]) and M.is_missing(d.get("renda").data[3])
    run("infile str10 nome int(idade renda) using f.raw if idade < 40, clear automatic")
    d = run.session.data
    assert d.nobs == 2 and d.get("idade").vtype == "int"
    assert d.value_labels == {"renda": {1: "abc"}} and d.get("renda").value_label == "renda"
    out = run("infile a using f.raw")
    assert run.rc == 4


def test_infix(run, here):
    (here / "fx.raw").write_text("000123Ana       031150\n000124Bia Lima  .  2000\n")
    assert run("infix long id 1-6 str nome 7-16 idade 17-19 renda 20-23 using fx, clear") \
        == "(2 observations read)\n"
    d = run.session.data
    assert list(d.get("id").data) == [123, 124] and d.get("nome").vtype == "str10"
    assert list(d.get("nome").data) == ["Ana", "Bia Lima"] and M.is_missing(d.get("idade").data[1])
    (here / "ml.raw").write_text("1001Ana\n  25 1500\n1002Bia\n  31  200\n")
    run("infix 2 lines 1: id 1-4 str nome 5-7 2: idade 1-4 renda 5-9 using ml.raw, clear")
    d = run.session.data
    assert list(d.get("renda").data) == [1500, 200] and list(d.get("nome").data) == ["Ana", "Bia"]
    lay = parse_infix_spec("3 firstlineoffile 2 lines id 1-4 / x 1-2")
    assert lay.first == 3 and lay.lines == 2 and [f.line for f in lay.fields] == [1, 2]


def test_dictionaries(run, here):
    (here / "fx.raw").write_text("000123Ana       031150\n000124Bia Lima  .  2000\n")
    (here / "d.dct").write_text(
        'dictionary using fx.raw {\n* comentário\n_column(1) long id %6f "Identificador"\n'
        '  str10 nome :nm %10s "Nome"\n  int idade %3f\n  renda %4.1f\n}\n')
    run("infile using d, clear")
    d = run.session.data
    assert d.names == ["id", "nome", "idade", "renda"]
    assert d.get("id").label == "Identificador" and d.get("nome").value_label == "nm"
    assert list(d.get("renda").data) == [15, 200]            # 1 casa decimal implícita
    (here / "i.dct").write_text('infile dictionary {\n str5 uf\n int pop\n}\n"MG" 200\n"SP" 450\n')
    run("infile using i.dct, clear")
    assert list(run.session.data.get("uf").data) == ["MG", "SP"]
    (here / "x.dct").write_text("infix dictionary using fx.raw {\n long id 1-6\n str nome 7-16\n}\n")
    run("infix using x.dct, clear")
    assert list(run.session.data.get("id").data) == [123, 124]
