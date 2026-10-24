"""Block "Data management: options": binary files, infile byvariable(), describe and list options,
putexcel to .xls and the SAS, SPSS and XPORT readers (compat case 1701)."""

import datetime as dt
from types import SimpleNamespace

import numpy as np
import pytest

from opendta.core import binfmt
from opendta.core import missing as M
from opendta.io import statpkg


# -- binary formats ----------------------------------------------------------------------------

@pytest.mark.parametrize("fmt,value,size", [
    ("%1b", 100, 1), ("%1b", -127, 1), ("%2b", 32740, 2), ("%4b", -2147483647, 4),
    ("%1bu", 255, 1), ("%2bs", -32767, 2), ("%4bu", 4294967295, 4),
    ("%8z", 3.14159, 8), ("%8z", -1e300, 8)])
@pytest.mark.parametrize("order", [binfmt.HILO, binfmt.LOHI])
def test_binfmt_round_trip(fmt, value, size, order):
    spec = binfmt.parse_spec(fmt)
    raw = binfmt.encode(value, spec, order)
    assert len(raw) == size
    assert binfmt.decode(raw, spec, order) == value


def test_binfmt_missing_and_ranges():
    b = binfmt.parse_spec("%1b")
    assert binfmt.encode(M.SYSMISS, b) == bytes([101])
    assert binfmt.encode(M.EXTENDED["a"], b) == bytes([102])
    assert binfmt.decode(bytes([127]), b) == M.EXTENDED["z"]
    assert binfmt.decode(binfmt.encode(2137, b), b) == M.SYSMISS          # out of range: missing
    assert binfmt.decode(binfmt.encode(2137, binfmt.parse_spec("%1bs")), binfmt.parse_spec("%1bs")) == 127
    assert binfmt.decode(binfmt.encode(2137, binfmt.parse_spec("%1bu")), binfmt.parse_spec("%1bu")) == 255
    assert binfmt.decode(binfmt.encode(124.75, binfmt.parse_spec("%2b")), binfmt.parse_spec("%2b")) == 124
    f = binfmt.parse_spec("%4z")
    assert binfmt.decode(binfmt.encode(M.SYSMISS, f), f) == M.SYSMISS
    assert binfmt.decode(binfmt.encode(M.EXTENDED["c"], f), f) == M.EXTENDED["c"]
    assert binfmt.decode(binfmt.encode(0.5, f), f) == 0.5
    d = binfmt.parse_spec("%8z")
    assert binfmt.decode(binfmt.encode(M.EXTENDED["q"], d), d) == M.EXTENDED["q"]


def test_binfmt_strings_and_bad_formats():
    s = binfmt.parse_spec("%5s")
    assert binfmt.encode("ab", s) == b"ab\0\0\0"
    assert binfmt.encode("abcdefg", s) == b"abcde"
    assert binfmt.decode(b"ab\0\0\0", s) == "ab"
    for bad in ("%3b", "%9z", "%0s", "%f", "x"):
        with pytest.raises(ValueError):
            binfmt.parse_spec(bad)


def test_file_binary(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run("file open h using b.bin, write binary\n"
        "file write h %14s \"mymatout 1.0.1\" %1b (byteorder())\n"
        "file write h %2b (3) %4b (-5) %8z (3.14159) %4z (.) %1b (.a) %1bu (300) %3s \"abcdef\"\n"
        "file seek h query")
    assert run.eval("r(loc)") == 14 + 1 + 2 + 4 + 8 + 4 + 1 + 1 + 3
    run("file close h")
    assert (tmp_path / "b.bin").stat().st_size == 38
    run("file open h using b.bin, read binary\n"
        "file read h %14s sig\nfile read h %1b bo\nlocal o = bo\nfile set h byteorder `o'")
    assert run.session.macros.get_local("sig") == "mymatout 1.0.1"
    run("file read h %2b v")
    assert run.eval("v") == 3
    run("file read h %4b v")
    assert run.eval("v") == -5
    run("file read h %8z v")
    assert run.eval("v") == pytest.approx(3.14159)
    run("file read h %4z v")
    assert run.eval("v") == M.SYSMISS
    run("file read h %1b v")
    assert run.eval("v") == M.EXTENDED["a"]
    run("file read h %1bu v")
    assert run.eval("v") == 255
    run("file read h %3s t")
    assert run.session.macros.get_local("t") == "abc"
    assert run.eval("r(eof)") == 0
    run("file read h %2b v")
    assert run.eval("r(eof)") == 1
    run("file seek h tof\nfile read h %2s t")
    assert run.session.macros.get_local("t") == "my"
    run("file seek h 3")
    assert run.eval("r(loc)") == 3
    run("file query")
    assert run.eval("r(N)") == 1
    run("file close h")
    run("file query")
    assert run.eval("r(N)") == 0


def test_file_binary_byteorder_and_read_write(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run("file open h using o.bin, write binary\nfile set h byteorder hilo\nfile write h %2b (258)\nfile close h")
    assert (tmp_path / "o.bin").read_bytes() == bytes([1, 2])
    run("file open h using o.bin, write binary replace\nfile set h byteorder 2\nfile write h %2b (258)\nfile close h")
    assert (tmp_path / "o.bin").read_bytes() == bytes([2, 1])
    run("file open h using rw.bin, read write binary\nfile write h %4b (7) %2b (9)\nfile seek h tof\nfile read h %4b v")
    assert run.eval("v") == 7
    run("file close all")
    run("file open h using o.bin, read binary")
    run("file write h %1b (1)")
    assert run.rc == 608
    run("file close h")
    run("file open h using o.bin, read binary\nfile write h %q (1)")
    assert run.rc == 608
    run("file close h")
    run("file open h using x.bin, write binary\nfile write h %3s (5)")
    assert run.rc == 109


# -- infile, byvariable() ------------------------------------------------------------------------

def test_infile_byvariable(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "time.raw").write_text("1980 1981 1982 1983\n14 17 25 30\n120 135 150\n180\n")
    out = run("infile year amount cost using time, byvariable(4) clear")
    assert "(4 observations read)" in out
    d = run.session.data
    assert list(d.get("year").data) == [1980, 1981, 1982, 1983]
    assert list(d.get("amount").data) == [14, 17, 25, 30]
    assert list(d.get("cost").data) == [120, 135, 150, 180]
    (tmp_path / "t2.raw").write_text("1980 1981 1982 1983;\n14 17 25 ;\n120 135 150 180;\n")
    out = run("infile year amount cost using t2, byvariable(100) clear")     # as in Stata 14 (compat 1701)
    assert "'1983;' cannot be read as a number for year[4]" in out and "(7 observations read)" in out
    d = run.session.data
    amount = d.get("amount").data
    assert list(d.get("year").data[4:]) == [14, 17, 25] and list(amount[:3]) == [120, 135, 150]
    assert amount[3] == M.SYSMISS and d.get("cost").data[0] == M.SYSMISS
    run("infile year amount cost using t2, byvariable(x) clear")
    assert run.rc == 198
    run("infile str5 s using t2, byvariable(2) clear")
    assert run.rc == 0 and run.session.data.nobs == 2


# -- describe and list ---------------------------------------------------------------------------

DATA = """clear
input str8 make price mpg
"AMC" 4099 22
"Buick" 5788 18
"Cad" 7827 .
end
label variable price "Price (USD)"
label define yn 1 "one"
label values mpg yn
sort price
"""


def test_describe_numbers_varlist(run):
    run(DATA)
    out = run("describe price mpg, numbers")
    assert "  2. price" in out and "  3. mpg" in out and "make" not in out
    run("describe, varlist")
    assert run.session.r["varlist"] == "make price mpg"
    assert run.session.r["sortlist"] == "price"
    assert run.eval("r(k)") == 3


def test_describe_replace(run):
    run(DATA)
    run("describe, replace")
    assert run.rc == 4                          # the data were not saved
    out = run("describe, replace clear")
    assert out == ""
    d = run.session.data
    assert d.names == ["position", "name", "type", "isnumeric", "format", "vallab", "varlab"]
    assert list(d.get("name").raw) == ["make", "price", "mpg"]
    assert list(d.get("isnumeric").data) == [0, 1, 1]
    assert list(d.get("vallab").raw) == ["", "", "yn"]
    assert list(d.get("varlab").raw) == ["", "Price (USD)", ""]
    run("describe, clear")
    assert run.rc == 198


def test_describe_using(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run(DATA + "save d, replace\nclear")
    out = run("describe using d")
    assert "Contains data" in out and "from d.dta" not in out and "Price (USD)" in out
    assert run.session.data.nvars == 0                                   # nothing was loaded
    assert run("describe using d, simple").split() == ["make", "price", "mpg"]
    out = run("describe make using d, varlist")
    assert "price" not in out and run.session.r["varlist"] == "make price mpg"
    assert "changed" not in run.session.r
    run("describe using nothing")
    assert run.rc == 601


def test_list_header_and_nodotz(run):
    run("clear\ninput x y\n1 2\n3 .a\n5 .\n7 8\nend")
    out = run("list, header(2)").splitlines()
    heads = [k for k, ln in enumerate(out) if "| x " in ln]
    row = {n: next(k for k, ln in enumerate(out) if ln.lstrip().startswith(f"{n}.")) for n in (2, 3)}
    assert len(heads) == 2 and row[2] < heads[1] < row[3]
    clean = run("list, clean header(2)")
    assert clean.count("x") == 2
    assert ".a" in run("list") and ".a" in run("list, nodotz")
    run("input z\n.z\n1\nend")
    assert ".z" in run("list z") and ".z" not in run("list z, nodotz")
    run("list, noheader header")
    assert run.rc == 198


def test_list_summary_rows(run):
    run("clear\ninput g c p\n1 47 5\n2 123 10\n3 22 2\nend")
    out = run("list g c p, sum(c p) labvar(g)")
    last = [ln for ln in out.splitlines() if "Sum" in ln][0]
    assert last.split() == ["|", "Sum", "192", "17", "|"]
    out = run("list, mean N").splitlines()
    assert any(ln.split()[:2] == ["Mean", "|"] and "5.666667" in ln for ln in out)
    assert any(ln.split()[:2] == ["N", "|"] for ln in out)
    run("list, nocompress")
    assert run.rc == 0
    out = run("list, labvar(zz) sum")
    assert run.rc == 198 and "variable zz not found" in out


# -- SAS, SPSS and XPORT ---------------------------------------------------------------------------

def _meta(**kw):
    base = dict(column_names=[], column_labels=None, variable_value_labels={}, file_label="")
    base.update(kw)
    return SimpleNamespace(**base)


def test_dataset_from_columns():
    cols = {"ID": [1.0, 2.0, 3.0], "Weight kg": [1.5, float("nan"), 3.0], "name": ["a", "bcd", None],
            "sex": [1.0, 2.0, float("nan")], "when": [dt.date(1960, 1, 2), dt.date(2000, 1, 1), None],
            "at": [dt.datetime(1960, 1, 1, 0, 0, 1)] * 3, "big": [100000.0, 2.0, 3.0],
            "empty": [float("nan")] * 3}
    meta = _meta(column_names=list(cols), column_labels=["Id", None, "Name", "Sex", "d", "t", "P", "e"],
                 variable_value_labels={"sex": {1.0: "male", 2.0: "female"}}, file_label="hello")
    ds = statpkg.dataset_from_columns(cols, meta, case="lower")
    assert ds.names == ["id", "v2", "name", "sex", "when", "at", "big", "empty"]
    assert [v.vtype for v in ds.vars] == ["byte", "double", "str3", "byte", "double", "double", "long", "byte"]
    assert ds.get("v2").data[1] == M.SYSMISS and ds.get("v2").data[2] == 3.0
    assert ds.get("name").raw[2] == ""
    assert ds.get("when").fmt == "%td" and list(ds.get("when").data[:2]) == [1.0, 14610.0]
    assert ds.get("at").fmt == "%tc" and ds.get("at").data[0] == 1000.0
    assert ds.get("sex").value_label == "labels1" and ds.value_labels["labels1"] == {1: "male", 2: "female"}
    assert ds.get("sex").data[2] == M.SYSMISS
    assert ds.chars["v2"]["original_name"] == "Weight kg"
    assert ds.get("id").label == "Id" and ds.label == "hello"
    assert statpkg.dataset_from_columns(cols, meta, case="upper").names[0] == "ID"
    assert statpkg.dataset_from_columns(cols, meta).names[0] == "ID"
    assert statpkg.select(statpkg.dataset_from_columns(cols, meta), ["SEX", "id"]).names == ["ID", "sex"]
    with pytest.raises(Exception):
        statpkg.select(statpkg.dataset_from_columns(cols, meta), ["nothing"])


def test_legal_names():
    names, original = statpkg.legal_names(["a", "1x", "a", "if", "x" * 40], "preserve")
    assert names == ["a", "v2", "v3", "v4", "v5"]
    assert original["v3"] == "a" and original["v5"] == "x" * 40


def test_import_spss(run, tmp_path, monkeypatch):
    pyreadstat = pytest.importorskip("pyreadstat")
    pd = pytest.importorskip("pandas")
    monkeypatch.chdir(tmp_path)
    df = pd.DataFrame({"Id": [1.0, 2.0, 3.0], "grp": [1.0, 2.0, 1.0], "name": ["a", "bb", "c"],
                       "w": [1.5, np.nan, 2.5]})
    pyreadstat.write_sav(df, str(tmp_path / "t.sav"), column_labels=["Identifier", "Group", "Name", "W"],
                         variable_value_labels={"grp": {1: "one", 2: "two"}})
    out = run("import spss using t")
    assert "(4 vars, 3 obs)" in out
    assert run.eval("r(N)") == 3 and run.eval("r(k)") == 4
    d = run.session.data
    assert d.names == ["Id", "grp", "name", "w"]
    assert d.get("Id").label == "Identifier" and d.get("grp").value_label == "labels1"
    assert d.value_labels["labels1"] == {1: "one", 2: "two"}
    assert d.get("w").data[1] == M.SYSMISS and d.get("name").vtype == "str2"
    run("import spss using t")
    assert run.rc == 4                         # data in memory would be lost (unsaved)
    run("import spss id w if w < 2 using t, clear case(lower)")
    assert run.session.data.names == ["id", "w"] and run.session.data.nobs == 1
    run("import spss using nothing, clear")
    assert run.rc == 601


def test_import_export_sasxport(run, tmp_path, monkeypatch):
    pytest.importorskip("pyreadstat")
    pytest.importorskip("pandas")
    monkeypatch.chdir(tmp_path)
    run("clear\ninput id str6 name weight\n1 \"ana\" 60.5\n2 \"bruno\" .\n3 \"carla\" 72\nend\n"
        "label variable weight \"Weight (kg)\"")
    assert "odta_x.xpt saved" in run("export sasxport odta_x") or (tmp_path / "odta_x.xpt").exists()
    run("export sasxport odta_x")
    assert run.rc == 602
    run("export sasxport odta_x, replace")
    assert run.rc == 0
    run("clear")
    out = run("import sasxport odta_x")
    assert out.strip() == ""                    # import sasxport is silent in Stata 14 (compat 1701)
    d = run.session.data
    assert d.names == ["id", "name", "weight"]
    assert [v.vtype for v in d.vars] == ["double", "str6", "double"]
    assert list(d.get("id").data) == [1, 2, 3] and d.get("weight").data[1] == M.SYSMISS
    assert d.get("weight").label == "Weight (kg)" and list(d.get("name").raw) == ["ana", "bruno", "carla"]
    run("clear\ninput averylongname\n1\nend")
    run("export sasxport odta_y")
    assert run.rc == 198
    run("clear\ninput x\n1\nend\nlabel define l 1 \"one\"\nlabel values x l")
    run("export sasxport odta_z")
    assert run.rc == 0 and (tmp_path / "formats.xpf").exists()


def test_sasxport_value_labels_members_and_rename(run, tmp_path, monkeypatch):
    pytest.importorskip("pyreadstat")
    pytest.importorskip("pandas")
    monkeypatch.chdir(tmp_path)
    run("clear\ninput id sex averylongname\n1 1 5\n2 2 .\nend")
    run('label define sexlab 1 "Female" 2 "Male"\nlabel values sex sexlab')
    run("export sasxport odta_y")
    assert run.rc == 198
    run("export sasxport odta_y, rename vallabfile(both)")
    assert run.rc == 0 and (tmp_path / "odta_y.sas").exists()
    out = run("import sasxport odta_y, describe")
    assert "member ODTA_Y" in out and "averylon" in out and "SEXLAB" in out
    run("import sasxport odta_y, clear")
    d = run.session.data
    assert d.names == ["id", "sex", "averylon"]
    assert d.get("sex").value_label == "sexlab" and d.value_labels["sexlab"] == {1: "Female", 2: "Male"}
    run("import sasxport odta_y, clear novallabels member(odta_y)")
    assert run.session.data.get("sex").value_label == ""
    run("import sasxport odta_y, clear member(nope)")
    assert run.rc == 198


def test_xport_reader_reads_every_member(tmp_path):
    pyreadstat = pytest.importorskip("pyreadstat")
    pd = pytest.importorskip("pandas")
    from opendta.io import xport
    pyreadstat.write_xport(pd.DataFrame({"A": [1.5, None]}), str(tmp_path / "a.xpt"), table_name="ONE",
                           file_format_version=5)
    pyreadstat.write_xport(pd.DataFrame({"B": ["x", "yz"]}), str(tmp_path / "b.xpt"), table_name="TWO",
                           file_format_version=5)
    both = (tmp_path / "a.xpt").read_bytes() + (tmp_path / "b.xpt").read_bytes()[3 * 80:]
    (tmp_path / "c.xpt").write_bytes(both)
    members = xport.read(tmp_path / "c.xpt")
    assert [m.name for m in members] == ["ONE", "TWO"]
    assert members[0].rows == [[1.5], [None]] and members[1].rows == [["x"], ["yz"]]


# -- putexcel to .xls ----------------------------------------------------------------------------------

def test_putexcel_xls(run, tmp_path, monkeypatch):
    xlrd = pytest.importorskip("xlrd")
    pytest.importorskip("xlwt")
    monkeypatch.chdir(tmp_path)
    run('putexcel A1=("id") B1=("score") A2=(1) B2=(2.5) using t.xls, replace')
    run("putexcel set t.xls, modify")
    run('putexcel C1=("n") C2=(3)')
    sh = xlrd.open_workbook(str(tmp_path / "t.xls")).sheet_by_index(0)
    assert [sh.cell_value(0, j) for j in range(3)] == ["id", "score", "n"]
    assert [sh.cell_value(1, j) for j in range(3)] == [1, 2.5, 3]
    run('putexcel D1=("s2") using t.xls, sheet("two")')
    book = xlrd.open_workbook(str(tmp_path / "t.xls"))
    assert book.sheet_names() == ["Sheet1", "two"] and book.sheet_by_name("two").cell_value(0, 3) == "s2"
    run("clear\ninput a b\n1 2\nend\nexport excel using t.xls, sheet(\"three\") sheetmodify")
    run("matrix m = (1, 2 \\ 3, 4)\nputexcel A5=matrix(m) using t.xls")
    book = xlrd.open_workbook(str(tmp_path / "t.xls"))
    assert "three" in book.sheet_names() and book.sheet_by_index(0).cell_value(5, 1) == 4


def test_list_string_and_display(run):
    run('clear\ninput str20 make price\n"AMC Concord" 4099\n"Buick LeSabre" 5788\nend')
    out = run("list, string(5)")
    assert "AMC C" in out and "Concord" not in out
    out = run("list in 1, display")
    assert out.splitlines()[2] == "  1. |        make | price |"
    assert out.splitlines()[3] == "     | AMC Concord |  4099 |"
    run("set linesize 20")
    out = run("list in 1")                     # wider than the line: display format
    assert "     |------" in out


def test_table_supercolumns_and_by(run):
    run("clear\nset obs 12\ngen a = 1 + mod(_n,2)\ngen b = 1 + mod(_n,3)\ngen c = 1 + mod(_n,2)")
    out = run("table a b c")
    lines = out.splitlines()
    assert lines[2].split("|")[1].strip() == "c and b"
    assert lines[3].split("|")[1].split() == ["------", "1", "-------", "------", "2", "-------"]
    out = run("table a, by(c) row")
    assert "c and a   |" in out and "1         |" in out and "Total |" in out
