"""Leitura e gravação de .dta (formatos 117, 118, 119).

Três camadas de verificação:
1. ida e volta pelo próprio leitor/gravador;
2. leitura cruzada pela ReadStat (pyreadstat), que é uma implementação
   independente da especificação — se ela lê igual, o arquivo está conforme;
3. leitura, pelo OpenDTA, de arquivos gravados pela ReadStat (versões 8–15).
O juiz final é o Stata (roteiro em compat/do/0104_arquivos.do).
"""

import numpy as np
import pytest

from opendta.core import missing as M
from opendta.io.dta import read_dta, write_dta

RELEASES = (117, 118, 119)


def build(run, release):
    accent = "Ana Lúcia" if release != 117 else "José"
    run(f"""clear
input byte b int i long l float f double d str12 s
1 300 70000 1.5 0.1 "{accent}"
2 -5 -70000 . 1e300 "Bia"
. 7 1 2.25 -3 ""
end
replace b = .a in 2
replace d = .z in 3
gen strL longo = "x" * 3000 in 1
replace longo = "curto" in 2
label data "Teste dta"
label variable s "Nome"
label define sim 1 "Sim" 2 "Não", replace
label values b sim
format d %12.4f
notes: nota do arquivo
notes i: nota da variável
sort l
""")
    return run.session.data


@pytest.fixture(params=RELEASES)
def release(request):
    return request.param


def _check_same(a, b):
    assert a.names == b.names
    assert a.nobs == b.nobs
    for va, vb in zip(a.vars, b.vars):
        assert va.vtype == vb.vtype, va.name
        assert va.fmt == vb.fmt, va.name
        assert va.label == vb.label, va.name
        assert va.value_label == vb.value_label, va.name
        if va.is_string:
            assert list(va.data) == list(vb.data), va.name
        else:
            assert np.array_equal(va.data.view(np.int64), vb.data.view(np.int64)), va.name


def test_round_trip(run, tmp_path, release):
    ds = build(run, release)
    p = tmp_path / f"t{release}.dta"
    ts = write_dta(ds, p, release=release)
    back = read_dta(p)
    _check_same(ds, back)
    assert back.label == "Teste dta"
    assert back.timestamp == ts
    assert back.value_labels["sim"] == ds.value_labels["sim"]
    assert back.sortlist == ["l"]
    assert back.chars["_dta"]["note1"] == "nota do arquivo"
    assert back.chars["i"]["note1"] == "nota da variável"
    assert M.missing_name(back.get("b").data[0]) == ".a"
    assert M.missing_name(back.get("d").data[1]) == ".z"
    assert back.get("longo").data[2] == "x" * 3000


def test_header_bytes(run, tmp_path, release):
    ds = build(run, release)
    p = tmp_path / "h.dta"
    write_dta(ds, p, release=release)
    raw = p.read_bytes()
    assert raw.startswith(f"<stata_dta><header><release>{release}</release>".encode())
    assert raw.endswith(b"</stata_dta>")


def test_latin1_in_117(run, tmp_path):
    ds = build(run, 117)
    p = tmp_path / "l.dta"
    write_dta(ds, p, release=117)
    assert "José".encode("latin-1") in p.read_bytes()
    assert read_dta(p).get("s").data[2] == "José"


# -- verificação independente pela ReadStat -------------------------------------------

def test_readstat_reads_ours(run, tmp_path, release):
    pyreadstat = pytest.importorskip("pyreadstat")
    ds = build(run, release)
    p = tmp_path / f"r{release}.dta"
    write_dta(ds, p, release=release)
    df, meta = pyreadstat.read_dta(str(p), user_missing=True, output_format="dict")
    assert list(df) == ds.names
    assert meta.file_label == "Teste dta"
    assert meta.column_names_to_labels["s"] == "Nome"
    assert meta.variable_value_labels["b"] == {1.0: "Sim", 2.0: "Não"}
    assert meta.readstat_variable_types == {
        "b": "int8", "i": "int16", "l": "int32", "f": "float", "d": "double",
        "s": "string", "longo": "string"}
    assert meta.original_variable_types["d"] == "%12.4f"
    assert list(df["l"]) == [-70000, 1, 70000]
    assert list(df["i"]) == [-5, 7, 300]
    assert df["f"][2] == 1.5
    assert df["b"][0] == "a"                 # .a vira "a" com user_missing
    assert df["d"][1] == "z"
    assert df["longo"][2] == "x" * 3000


@pytest.mark.parametrize("version", [8, 10, 11, 12, 13, 14, 15])
def test_we_read_readstat(tmp_path, version):
    pyreadstat = pytest.importorskip("pyreadstat")
    pd = pytest.importorskip("pandas")      # a gravação pela ReadStat pede um DataFrame
    df = pd.DataFrame({"n": [1.5, None, 3.0], "k": [1, 2, 3], "t": ["a", "b", "ç"]})
    p = tmp_path / f"v{version}.dta"
    pyreadstat.write_dta(df, str(p), version=version, file_label="rs",
                         column_labels=["num", "int", "txt"])
    ds = read_dta(p)
    assert ds.names == ["n", "k", "t"]
    assert ds.get("n").data[0] == 1.5
    assert M.is_missing(ds.get("n").data[1])
    assert list(ds.get("k").data) == [1, 2, 3]
    assert ds.get("t").data[2] == "ç"
    assert ds.get("t").label == "txt"
    assert ds.label == "rs"


# -- comandos ------------------------------------------------------------------------

@pytest.fixture
def here(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_save_messages(run, here):
    run("clear\nset obs 2\ngen x = _n")
    assert run("save a") == "file a.dta saved\n"
    out = run("save a")
    assert run.rc == 602 and "file a.dta already exists" in out
    assert run("save a, replace") == "file a.dta saved\n"
    assert run("save b, replace") == "(note: file b.dta not found)\nfile b.dta saved\n"
    out = run("save")
    assert run.rc == 602 and "file b.dta already exists" in out
    assert run("save, replace") == "file b.dta saved\n"


def test_use_messages(run, here):
    run('clear\nset obs 2\ngen x = _n\nlabel data "Rótulo"\nsave a')
    run("gen y = 1")
    out = run("use a")
    assert run.rc == 4 and "no; data in memory would be lost" in out
    assert run("use a, clear") == "(Rótulo)\n"
    assert run.session.data.names == ["x"]
    out = run("use nada, clear")
    assert run.rc == 601 and "file nada.dta not found" in out


def test_use_subset(run, here):
    run("clear\nset obs 4\ngen x = _n\ngen y = 10*_n\nsave a")
    run("use y using a if x > 2, clear")
    d = run.session.data
    assert d.names == ["y"]
    assert list(d.get("y").data) == [30, 40]
    run("use a in 2/3, clear")
    assert list(run.session.data.get("x").data) == [2, 3]


def test_saveold_writes_117(run, here):
    run("clear\nset obs 1\ngen x = 1")
    out = run("saveold v13")
    assert out.endswith("file v13.dta saved\n")
    assert (here / "v13.dta").read_bytes().startswith(b"<stata_dta><header><release>117")
    run("save v14")
    assert (here / "v14.dta").read_bytes().startswith(b"<stata_dta><header><release>118")


def test_notes_and_char(run, here):
    run("clear\nset obs 1\ngen x = 1\nnotes: primeira\nnotes: segunda\nnotes x: da x")
    out = run("notes")
    assert "_dta:" in out and "  1.  primeira" in out and "  2.  segunda" in out
    assert "x:" in out and "  1.  da x" in out
    run("notes drop _dta")
    assert "_dta:" not in run("notes")
    run('char x[fonte] "censo"')
    assert "x[fonte]:" in run("char list")
    run("save c\nuse c, clear")
    assert run.session.data.chars["x"]["fonte"] == "censo"
    assert "* indicated variables have notes" in run("describe")


def test_numeric_edges(run, tmp_path, release):
    """Negativos e extremos de cada tipo (um float negativo já virou missing)."""
    import numpy as np

    from opendta.core.dataset import Dataset, Variable
    ds = Dataset()
    f32 = lambda xs: np.array(xs, dtype=np.float32).astype(np.float64)
    cols = {
        "byte": np.array([-127, 100, -1, 0]),
        "int": np.array([-32767, 32740, -1, 0]),
        "long": np.array([-2147483647, 2147483620, -1, 0]),
        "float": f32([-1.7014117331926443e38, 1.7014117331926443e38, -0.1791305, -1e-30]),
        "double": np.array([-8.98846567431158e307, 8.98846567431158e307, -0.5, -1e-300]),
    }
    for t, vals in cols.items():
        ds.vars.append(Variable(f"x_{t}", t, np.asarray(vals, dtype=np.float64)))
    ds.nobs = 4
    p = tmp_path / "edges.dta"
    write_dta(ds, p, release=release)
    back = read_dta(p)
    for a, b in zip(ds.vars, back.vars):
        assert np.array_equal(a.data, b.data), a.name


def test_chunked_io_and_atomic_replace(run, tmp_path, release, monkeypatch):
    """Blocos minúsculos forçam vários blocos de leitura e gravação."""
    import opendta.io.dta as D
    monkeypatch.setattr(D, "CHUNK_BYTES", 64)
    ds = build(run, release)
    p = tmp_path / "c.dta"
    write_dta(ds, p, release=release)
    _check_same(ds, read_dta(p))
    # falha no meio da gravação: o arquivo antigo continua intacto
    before = p.read_bytes()
    monkeypatch.setattr(D, "_write_body", lambda *a, **k: (_ for _ in ()).throw(OSError("disco cheio")))
    with pytest.raises(OSError):
        write_dta(ds, p, release=release)
    assert p.read_bytes() == before
    assert not list(tmp_path.glob("*.opendta-tmp"))


def test_native_storage_sizes(run):
    run("clear\nset obs 1000\ngen byte b = 1\ngen int i = 1\ngen long l = 1\n"
        "gen float f = 1\ngen double d = 1")
    d = run.session.data
    assert [d.get(n).raw.nbytes for n in "bilfd"] == [1000, 2000, 4000, 4000, 8000]


@pytest.mark.parametrize("old", [114, 115])
def test_old_formats(run, tmp_path, old):
    """114/115: gravação própria, leitura própria e leitura pela ReadStat."""
    ds = build(run, 117)                     # acentos em Latin-1, como no Stata 12
    p = tmp_path / f"o{old}.dta"
    write_dta(ds, p, release=old)
    raw = p.read_bytes()
    assert raw[0] == old and raw[1] == 2
    back = read_dta(p)
    assert back.names == ds.names
    for a, b in zip(ds.vars, back.vars):
        if a.vtype == "strL":
            assert b.vtype == "str244" and b.data[2] == "x" * 244
            continue
        assert a.vtype == b.vtype and a.fmt == b.fmt and a.label == b.label
        if a.is_string:
            assert list(a.data) == list(b.data)
        else:
            assert np.array_equal(a.data.view(np.int64), b.data.view(np.int64)), a.name
    assert back.value_labels == ds.value_labels and back.sortlist == ["l"]
    assert back.chars["i"]["note1"] == "nota da variável"
    pyreadstat = pytest.importorskip("pyreadstat")
    df, meta = pyreadstat.read_dta(str(p), user_missing=True, output_format="dict")
    assert list(df) == ds.names and meta.file_label == "Teste dta"
    assert list(df["l"]) == [-70000, 1, 70000] and df["b"][0] == "a"
    assert meta.variable_value_labels["b"] == {1.0: "Sim", 2.0: "Não"}


def test_saveold_versions(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run('clear\nset obs 1\ngen strL t = "y" * 300\ngen x = 1')
    out = run("saveold v12, version(12)")
    assert "(saving in Stata 12 format)" in out and "variable t truncated to str244" in out
    assert (tmp_path / "v12.dta").read_bytes()[0] == 115
    run("saveold v11, version(11)")
    assert (tmp_path / "v11.dta").read_bytes()[0] == 114
    run("use v11, clear")
    assert run.session.data.get("t").vtype == "str244"
    out = run("saveold v9, version(9)")
    assert run.rc == 198
