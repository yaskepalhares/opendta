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
    assert run("export delimited out, replace") == \
        "(note: file out.csv not found)\nfile out.csv saved\n"
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


# -- leitura em blocos (CHUNK_ROWS pequeno força vários blocos) ------------------------

@pytest.fixture
def small_chunks(monkeypatch):
    import opendta.io.delimited as D
    monkeypatch.setattr(D, "CHUNK_ROWS", 3)


def test_chunked_late_text_and_ragged(here, small_chunks):
    # 'cod' é numérico nos primeiros blocos e vira texto na 8ª linha;
    # a 3ª coluna só aparece na 5ª linha; linhas em branco no fim somem
    lines = ["id,cod"] + [f"{i},0{i}" for i in range(1, 8)] + ["8,X9"]
    lines[5] += ",extra"
    (here / "c.csv").write_text("\n".join(lines) + "\n,\n\n")
    ds = read_delimited(here / "c.csv")
    assert ds.nobs == 8 and ds.names == ["id", "cod", "v3"]
    assert ds.get("cod").vtype == "str2"
    assert list(ds.get("cod").data) == ["01", "02", "03", "04", "05", "06", "07", "X9"]
    assert ds.get("v3").data[4] == "extra" and ds.get("v3").data[0] == ""


def test_chunked_numbers_and_missing(here, small_chunks):
    (here / "n.csv").write_text("a,b,c\n1,.a,nan\n2,,x\n3,1d2,y\n4,.,z\n")
    ds = read_delimited(here / "n.csv")
    assert list(ds.get("a").data) == [1, 2, 3, 4]
    from opendta.core import missing as M
    b = ds.get("b").data
    assert M.missing_name(b[0]) == ".a" and M.is_missing(b[1]) and b[2] == 100
    assert ds.get("c").is_string                    # 'nan' não é número no Stata


def test_rowrange_colrange_latin1(here, small_chunks):
    (here / "l.csv").write_bytes("nome,x,y\nJosé,1,2\nAna,3,4\nBia,5,6\nCaio,7,8\n".encode("latin-1"))
    ds = read_delimited(here / "l.csv", ReadOptions(rowrange=(2, 3), colrange=(1, 2)))
    assert ds.names == ["nome", "x"] and list(ds.get("nome").data) == ["Ana", "Bia"]
    assert read_delimited(here / "l.csv").get("nome").data[0] == "José"


def test_export_import_round_trip(run, here, small_chunks):
    run("clear\nset obs 10\ngen id = _n\ngen double d = _n / 4 - 1\ngen float f = -_n / 3\n"
        'gen str5 s = "a,b" in 1/3\nreplace s = `"q"q"\' in 4\nreplace d = .b in 5')
    run("export delimited t, replace")
    text = (here / "t.csv").read_text().splitlines()
    assert text[1] == '1,-.75,-.33333334,"a,b"'
    assert text[4] == '4,0,-1.3333334,"q""q"'
    assert text[5].split(",")[1] == ".b"
    run("import delimited t, clear")
    d = run.session.data
    assert list(d.get("id").data) == list(range(1, 11))
    assert d.get("s").data[3] == 'q"q'
    assert d.get("d").data[0] == -0.75
