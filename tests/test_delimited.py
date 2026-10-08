"""import/export delimited, insheet, outsheet e type."""

import pytest

from opendta.io.delimited import ReadOptions, number_text, read_delimited


@pytest.fixture
def here(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


CSV = 'Nome,Idade,Renda Mensal,uf\nAna,30,1500.5,SP\n"Bia, a ""grande""",,2000,RJ\nCaio,41,.,"MG"\n'


def test_read_types_and_names(here):
    (here / "p.csv").write_text(CSV)
    ds = read_delimited(here / "p.csv")
    assert ds.names == ["nome", "idade", "renda_mensal", "uf"]
    assert [v.vtype for v in ds.vars] == ["str15", "byte", "float", "str2"]
    assert ds.get("nome").data[1] == 'Bia, a "grande"'
    assert ds.get("renda_mensal").label == "Renda Mensal"
    assert ds.get("idade").label == ""


def test_read_tab_and_options(here):
    (here / "t.txt").write_text("a\tb\n1\t70000\n2\t3.25\n")
    ds = read_delimited(here / "t.txt")
    assert [v.vtype for v in ds.vars] == ["byte", "float"]
    ds = read_delimited(here / "t.txt", ReadOptions(asdouble=True, stringcols={1}))
    assert [v.vtype for v in ds.vars] == ["str1", "double"]
    ds = read_delimited(here / "t.txt", ReadOptions(varnames=0))
    assert ds.names == ["v1", "v2"] and ds.nobs == 3


def test_number_text():
    assert number_text(0.5, "double") == ".5"
    assert number_text(-0.25, "float") == "-.25"
    assert number_text(1500.5, "float") == "1500.5"
    assert number_text(0.1, "float") == ".1"
    assert number_text(70000, "long") == "70000"
    assert number_text(1e20, "double") == "1e+20"


def test_import_export_commands(run, here):
    (here / "p.csv").write_text(CSV)
    assert run("import delimited p.csv, clear") == "(4 vars, 3 obs)\n"
    assert run("export delimited out, replace") == ""
    assert (here / "out.csv").read_text() == (
        'nome,idade,renda_mensal,uf\nAna,30,1500.5,SP\n"Bia, a ""grande""",,2000,RJ\nCaio,41,,MG\n')
    out = run("export delimited out")
    assert run.rc == 602 and "file out.csv already exists" in out
    run('label define s 30 "trinta"\nlabel values idade s')
    run("export delimited nome idade using q if idade < ., quote")
    assert (here / "q.csv").read_text() == 'nome,idade\n"Ana","trinta"\n"Caio",41\n'
    run("export delimited nome idade using q, replace nolabel novarnames")
    assert (here / "q.csv").read_text().splitlines()[0] == "Ana,30"
    run("gen x = 1")
    out = run("import delimited p.csv")
    assert run.rc == 4


def test_insheet_outsheet(run, here):
    (here / "p.csv").write_text(CSV)
    run("insheet using p.csv, clear")
    run("outsheet nome uf using o, replace")
    assert (here / "o.out").read_text().splitlines()[1] == '"Ana"\t"SP"'
    run("outsheet nome uf using o.csv, comma noquote replace")
    assert (here / "o.csv").read_text().splitlines()[1] == "Ana,SP"


def test_type(run, here):
    (here / "a.txt").write_text("linha 1\na\tb\n")
    assert run("type a.txt") == "linha 1\na       b\n"
    assert run("type a.txt, showtabs") == "linha 1\na<T>b\n"
    out = run("type nada.txt")
    assert run.rc == 601 and "file nada.txt not found" in out
