"""Data-management commands: variables and observations, describing and comparing
data, files and folders, statsby and putmata, example datasets, ICD codes,
business calendars, unicode and odbc."""

import gzip

import numpy as np
import pytest

BASE = """clear
input id grp x y
1 1 1 10
2 1 2 .
3 1 4 30
4 2 1 5
5 2 3 .
6 2 5 25
end
label variable x "Dose given"
"""


def col(run, name):
    return run.session.data.get(name).data


# ---------------------------------------------------------------------------
# variables and observations
# ---------------------------------------------------------------------------

def test_clonevar_keeps_attributes(run):
    run(BASE + 'label define yn 0 "no" 1 "yes"\nlabel values grp yn\nformat x %5.1f')
    run("clonevar x2 = x if id < 4")
    v = run.session.data.get("x2")
    assert v.label == "Dose given" and v.fmt == "%5.1f"
    assert list(v.data[:3]) == [1, 2, 4] and v.data[3] >= 8.98e307
    run("clonevar g2 = grp")
    assert run.session.data.get("g2").value_label == "yn"


def test_ds_and_lookfor(run):
    run(BASE + "gen str3 s = \"a\"")
    run("ds, has(type string)")
    assert run.session.r["varlist"] == "s"
    run("lookfor dose")
    assert run.session.r["varlist"] == "x"


def test_insobs_and_range(run):
    run(BASE)
    run("insobs 2, before(1)")
    d = run.session.data
    assert d.nobs == 8 and col(run, "id")[0] >= 8.98e307 and col(run, "id")[2] == 1
    run("clear\nset obs 3\nrange t 0 1 5")
    assert run.session.data.nobs == 5
    assert list(col(run, "t")) == pytest.approx([0, .25, .5, .75, 1])


def test_separate(run):
    run(BASE)
    run("separate x, by(grp)")
    x1, x2 = col(run, "x1"), col(run, "x2")
    assert list(x1[:3]) == [1, 2, 4] and all(v >= 8.98e307 for v in x1[3:])
    assert list(x2[3:]) == [1, 3, 5]


def test_ipolate_by_group(run):
    run(BASE)
    run("by grp: ipolate y x, gen(yi)")
    yi = col(run, "yi")
    assert yi[1] == pytest.approx(10 + (2 - 1) * (30 - 10) / (4 - 1))
    assert yi[4] == pytest.approx(5 + (3 - 1) * (25 - 5) / (5 - 1))


def test_expandcl(run):
    run(BASE)
    run("expandcl 2, cluster(grp) generate(newcl)")
    assert run.session.data.nobs == 12
    assert sorted(set(col(run, "newcl"))) == [1, 2, 3, 4]


def test_stack_and_xpose(run):
    run("clear\nset obs 3\ngen a = _n\ngen b = 10*_n\nstack a b, into(z) clear")
    assert run.session.data.nobs == 6
    assert list(col(run, "_stack")) == [1, 1, 1, 2, 2, 2]
    assert list(col(run, "z")) == [1, 2, 3, 10, 20, 30]
    run("clear\nset obs 2\ngen a = _n\ngen b = 10*_n\nxpose, clear varname")
    d = run.session.data
    assert d.nobs == 2 and list(col(run, "v1")) == [1, 10] and list(col(run, "v2")) == [2, 20]


def test_corr2data_exact_moments(run):
    run("clear\nmatrix C = (1, .3 \\ .3, 1)\ncorr2data a b, n(200) corr(C) means(5 10) sds(2 3)")
    a, b = col(run, "a"), col(run, "b")
    assert a.mean() == pytest.approx(5, abs=1e-4) and b.std(ddof=1) == pytest.approx(3, rel=1e-5)
    assert np.corrcoef(a, b)[0, 1] == pytest.approx(.3, abs=1e-5)


def test_drawnorm_distribution(run):
    run("clear\nmatrix C = (1, .5 \\ .5, 1)\ndrawnorm a b, n(5000) corr(C) seed(7)")
    a, b = col(run, "a"), col(run, "b")
    assert abs(a.mean()) < .1 and np.corrcoef(a, b)[0, 1] == pytest.approx(.5, abs=.05)


# ---------------------------------------------------------------------------
# describing and comparing
# ---------------------------------------------------------------------------

def test_inspect_counts(run):
    run("clear\ninput x\n-1\n0\n2\n2.5\n.\nend")
    out = run("inspect x")
    assert "Negative" in out and "Missing" in out
    assert run.session.r["N_neg"] == 1 and run.session.r["N_unique"] == 4


def test_labelbook_lists_definitions(run):
    run(BASE + 'label define yn 0 "no" 1 "yes"\nlabel values grp yn')
    out = run("labelbook yn")
    assert "value label yn" in out and "yes" in out and "grp" in out


def test_cf_reports_mismatch(run, tmp_path):
    f = tmp_path / "a.dta"
    run(BASE + f'save "{f}"')
    run("replace x = 99 in 1")
    run(f'capture cf _all using "{f}"')
    assert run.session.eval("_rc") == 9
    run("replace x = 1 in 1")
    run(f'cf _all using "{f}"')
    assert run.rc == 0


def test_datasignature_confirm(run):
    run(BASE)
    run("datasignature set")
    run("datasignature confirm")
    assert run.rc == 0
    run("replace x = 7 in 1")
    run("capture datasignature confirm")
    assert run.session.eval("_rc") != 0


def test_snapshot_restore(run):
    run(BASE)
    run('snapshot save, label("before")')
    run("drop if grp == 2")
    assert run.session.data.nobs == 3
    run("snapshot restore 1")
    assert run.session.data.nobs == 6


def test_compare_counts(run):
    run(BASE + "gen x2 = x\nreplace x2 = x + 1 in 1")
    out = run("compare x x2")
    assert "x<x2" in out.replace(" ", "")


# ---------------------------------------------------------------------------
# files and folders
# ---------------------------------------------------------------------------

def test_files_and_folders(run, tmp_path):
    src = tmp_path / "a.txt"
    src.write_text("one;two\nthree;four\n", encoding="utf-8")
    d = tmp_path / "sub"
    run(f'mkdir "{d}"')
    assert d.is_dir()
    run(f'copy "{src}" "{d / "b.txt"}"')
    assert (d / "b.txt").read_text(encoding="utf-8") == src.read_text(encoding="utf-8")
    run(f'filefilter "{src}" "{tmp_path / "c.txt"}", from(";") to(",")')
    assert (tmp_path / "c.txt").read_text(encoding="utf-8") == "one,two\nthree,four\n"
    out = run(f'dir "{d}"')
    assert "b.txt" in out
    run(f'checksum "{src}"')
    assert run.session.r["filelen"] == src.stat().st_size
    (d / "b.txt").unlink()
    run(f'rmdir "{d}"')
    assert not d.exists()


def test_zip_roundtrip(run, tmp_path):
    import os
    (tmp_path / "f.txt").write_text("hello", encoding="utf-8")
    cwd = os.getcwd()
    run(f'cd "{tmp_path}"')
    try:
        run("zipfile f.txt, saving(z.zip)")
        (tmp_path / "f.txt").unlink()
        run("unzipfile z.zip")
        assert (tmp_path / "f.txt").read_text(encoding="utf-8") == "hello"
    finally:
        run(f'cd "{cwd}"')


def test_shell_runs_a_command(run):
    out = run("shell echo opendta_shell_ok")
    assert "opendta_shell_ok" in out


def test_outfile_infile_roundtrip(run, tmp_path):
    f = tmp_path / "d.raw"
    run(BASE + f'outfile id x y using "{f}"')
    run(f'clear\ninfile id x y using "{f}"')
    assert list(col(run, "x")) == [1, 2, 4, 1, 3, 5]
    assert col(run, "y")[1] >= 8.98e307


# ---------------------------------------------------------------------------
# statsby, putmata
# ---------------------------------------------------------------------------

def test_statsby_collects_results(run):
    run(BASE)
    run("statsby m=r(mean) n=r(N), by(grp) clear nodots: summarize x")
    assert list(col(run, "m")) == pytest.approx([7 / 3, 3])
    assert list(col(run, "n")) == [3, 3]


def test_putmata_matrix(run):
    run(BASE)
    run("putmata X=(id x)")
    out = run("mata: rows(X), cols(X)")
    assert "6" in out and "2" in out


# ---------------------------------------------------------------------------
# example datasets
# ---------------------------------------------------------------------------

def test_sysuse_own_datasets(run):
    out = run("sysuse dir")
    assert "diet.dta" in out
    run("sysuse diet, clear")
    assert run.session.data.nobs > 0
    run("capture sysuse auto, clear")
    assert run.session.eval("_rc") == 601


def test_webuse_default_address(run):
    out = run("webuse query")
    assert "stata-press.com" in out
    run("webuse set https://example.org/data")
    assert "https://example.org/data/" in run("webuse query")


# ---------------------------------------------------------------------------
# ICD codes (public-domain CMS tables)
# ---------------------------------------------------------------------------

def test_icd_tables_are_shipped():
    from opendta.core import icd
    assert len(icd.table("icd9")) > 14000 and len(icd.table("icd10")) > 90000
    assert icd.defined("icd9", "250.00") and icd.defined("icd9", "250")
    assert icd.defined("icd10", "E11.9") and not icd.defined("icd10", "E11.99")
    assert icd.describe("icd10", "E119") == "Type 2 diabetes mellitus without complications"


def test_icd_check_and_describe(run):
    run('clear\ninput str8 dx\n"E11.9"\n"E11.99"\n"I10"\n"1X"\nend')
    run("icd10 check dx, generate(prob)")
    prob = col(run, "prob")
    assert prob[0] == 0 and prob[2] == 0
    assert prob[1] == 9            # well formed but not defined
    assert prob[3] > 0
    run("icd10 generate d = dx, description")
    assert run.session.data.get("d").raw[0] == "Type 2 diabetes mellitus without complications"
    out = run("icd10 search diabetes kidney")
    assert "E11.22" in out


def test_icd9p_has_no_list(run):
    run("capture icd9p lookup 01.0")
    assert run.session.eval("_rc") == 198


# ---------------------------------------------------------------------------
# business calendars, unicode, odbc
# ---------------------------------------------------------------------------

def test_bcal_create_and_functions(run, tmp_path):
    import os
    cwd = os.getcwd()
    run(f'cd "{tmp_path}"')
    try:
        run("""clear
input str9 ds
"01nov2011"
"02nov2011"
"04nov2011"
"07nov2011"
end
gen d = date(ds, "DMY")
format d %td
bcal create mycal, from(d) generate(bd) replace""")
        assert (tmp_path / "mycal.stbcal").exists()
        bd = col(run, "bd")
        assert list(bd) == [0, 1, 2, 3]
        assert run.session.eval('dofb(2, "mycal")') == run.session.eval('mdy(11, 4, 2011)')
        assert run.session.eval('bofd("mycal", mdy(11, 7, 2011))') == 3
    finally:
        run(f'cd "{cwd}"')


def test_unicode_convertfile(run, tmp_path):
    src, dst = tmp_path / "a.txt", tmp_path / "b.txt"
    src.write_bytes("ação".encode("latin-1"))
    run(f'unicode convertfile "{src}" "{dst}", srcencoding(latin1) dstencoding(utf-8)')
    assert dst.read_text(encoding="utf-8") == "ação"


def test_odbc_explains_missing_pyodbc(run):
    pytest.importorskip("opendta")
    try:
        import pyodbc  # noqa: F401
        pytest.skip("pyodbc installed")
    except ImportError:
        pass
    out = run("capture noisily odbc list")
    assert "pyodbc" in out
