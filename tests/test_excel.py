"""import excel / export excel."""

import datetime as dt

import pytest

openpyxl = pytest.importorskip("openpyxl")

from opendta.core import missing as M  # noqa: E402
from opendta.io.excel import ExcelOptions, col_letters, col_number, parse_range, read_excel  # noqa: E402


@pytest.fixture
def here(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def make_book(path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Dados"
    ws.append(["Nome", "Idade", "Renda (R$)", "Nascimento", "Hora"])
    ws.append(["Ana", 30, 1500.5, dt.date(2020, 1, 1), dt.datetime(2020, 1, 1, 12, 30)])
    ws.append(["Bia", None, 2000, None, dt.datetime(1960, 1, 1, 0, 0, 1)])
    ws.append(["Caio", 41, "n/d", dt.date(1959, 12, 31), None])
    ws.append([None, None, None, None, None])
    wb.create_sheet("Vazia")
    wb.save(path)


def test_letters():
    assert [col_letters(n) for n in (1, 26, 27, 702)] == ["A", "Z", "AA", "ZZ"]
    assert col_number("AA") == 27
    assert parse_range("B2:D10") == ((2, 2), (10, 4))
    assert parse_range("C3") == ((3, 3), None)


def test_read_types(here):
    make_book(here / "b.xlsx")
    ds = read_excel(here / "b.xlsx", ExcelOptions(firstrow=True))
    assert ds.names == ["Nome", "Idade", "RendaR", "Nascimento", "Hora"]
    assert ds.nobs == 3                                 # linha vazia no fim some
    assert [v.vtype for v in ds.vars][:3] == ["str4", "byte", "str6"]
    assert ds.get("RendaR").label == "Renda (R$)"
    nasc = ds.get("Nascimento")
    assert nasc.fmt == "%tdnn/dd/CCYY" and list(nasc.data[[0, 2]]) == [21915, -1] and M.is_missing(nasc.data[1])
    hora = ds.get("Hora")
    assert hora.fmt == "%tcnn/dd/CCYY_hh:MM:SS" and hora.vtype == "double" and hora.data[1] == 1000
    plain = read_excel(here / "b.xlsx")
    assert plain.names == ["A", "B", "C", "D", "E"] and plain.nobs == 4


def test_cellrange_and_sheet(here):
    make_book(here / "b.xlsx")
    ds = read_excel(here / "b.xlsx", ExcelOptions(cellrange="B2:C3"))
    assert ds.names == ["B", "C"] and list(ds.get("C").data) == [1500.5, 2000]
    with pytest.raises(Exception):
        read_excel(here / "b.xlsx", ExcelOptions(sheet="Nada"))


def test_commands_round_trip(run, here):
    run('clear\ninput str5 nome idade double renda d\n"Ana" 30 1500.5 21915\n"Bia" . 2.25 .\nend\n'
        'format d %td\nlabel define id 30 "trinta"\nlabel values idade id\n'
        'label variable renda "Renda mensal"')
    assert run("export excel t.xlsx, firstrow(variables)") == "file t.xlsx saved\n"
    out = run("export excel t.xlsx")
    assert run.rc == 602 and "already exists" in out
    run('export excel nome using t.xlsx, sheet("Nomes") sheetmodify')
    out = run("import excel t.xlsx, describe")
    assert "Sheet1 | A1:D3" in out and "Nomes | A1:A2" in out
    assert run("import excel t.xlsx, clear firstrow") == ""      # o Stata não mostra nada
    d = run.session.data
    assert list(d.get("idade").data) == ["trinta", ""]       # rótulo exportado como texto
    assert list(d.get("renda").data) == [1500.5, 2.25]
    assert d.get("d").fmt == "%tdnn/dd/CCYY" and d.get("d").data[0] == 21915
    run("export excel t2.xlsx, firstrow(varlabels) nolabel")
    wb = openpyxl.load_workbook(here / "t2.xlsx")
    assert [c.value for c in wb.active[1]] == ["nome", "idade", "renda", "d"]
    assert wb.active["B2"].value == "trinta"
