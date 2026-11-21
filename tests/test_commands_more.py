"""Commands that research do-files use (compat 2005 and 2008): earlier names, small data and
matrix utilities, tests, quantile regressions, pk and forecast."""

import pytest


@pytest.mark.parametrize("cmd, rc, msg", [
    ("hettest", 301, "last estimates not found"), ("ovtest", 301, "last estimates not found"),
    ("lstat", 301, "last estimates not found"), ("treatreg", 301, "last estimates not found"),
    ("prop", 301, "last estimates not found"), ("svylc", 198, "invalid syntax"),
    ("svydes", 111, "no variables defined"),
    ("vif", 301, "use option uncentered to get uncentered VIFs"),
])
def test_earlier_names_without_arguments(run, cmd, rc, msg):
    out = run(f"clear\ncapture noisily {cmd}\ndisplay _rc")
    lines = out.splitlines()
    assert lines[-1] == str(rc) and lines[-2] == msg


def _regdata(run):
    run("clear\nset obs 30\ngen x1 = mod(_n * 7, 11)\ngen x2 = x1 + mod(_n, 4)\n"
        "gen x3 = mod(_n * 3, 5)\ngen y = x1 + 2 * x2 - x3 + mod(_n * 13, 7)")


def test_estat_vif(run):
    import numpy as np
    _regdata(run)
    run("quietly regress y x1 x2 x3")
    out = run("estat vif")
    assert "    Variable |       VIF       1/VIF  " in out and "Mean VIF" in out
    X = np.array(run.session.data.get("x1").data), np.array(run.session.data.get("x2").data)
    x1, x2 = X
    x3 = np.array(run.session.data.get("x3").data)
    Z = np.column_stack([x2, x3, np.ones(30)])
    e = x1 - Z @ np.linalg.lstsq(Z, x1, rcond=None)[0]
    vif1 = ((x1 - x1.mean()) @ (x1 - x1.mean())) / (e @ e)
    names = [run.eval(f"r(name_{i})") for i in (1, 2, 3)]
    vals = {n: run.eval(f"r(vif_{i})") for i, n in enumerate(names, 1)}
    assert vals["x1"] == pytest.approx(vif1)
    assert run.eval("r(vif_1) >= r(vif_2) & r(vif_2) >= r(vif_3)") == 1
    assert run("vif") == out
    run("quietly regress y x1 x2 x3, noconstant")
    assert run("capture estat vif\ndisplay _rc").strip() == "301"
    assert "_cons" not in run("estat vif, uncentered")


def test_hettest_is_estat_hettest(run):
    _regdata(run)
    run("quietly regress y x1 x2 x3")
    assert run("hettest") == run("estat hettest")
    assert run("ovtest") == run("estat ovtest")


def test_prop_abbreviation_replays(run):
    run("clear\nset obs 10\ngen g = mod(_n, 2)\nquietly proportion g")
    assert run("prop") == run("proportion")


def _vars(run):
    return [v.name for v in run.session.data.vars]


def test_aorder_move_renpfix(run):
    assert run("clear\ncapture noisily aorder\ndisplay _rc").strip() == "0"
    run("clear\nset obs 1\ngen c = 1\ngen a = 2\ngen B = 3\ngen x1 = 4\ngen x2 = 5")
    run("aorder")
    assert _vars(run) == ["a", "B", "c", "x1", "x2"]
    run("move x2 a")
    assert _vars(run) == ["x2", "a", "B", "c", "x1"]
    assert run("capture noisily move x2\ndisplay _rc").splitlines()[-1] == "102"
    run("renpfix x z")
    assert _vars(run) == ["z2", "a", "B", "c", "z1"]
    assert run.eval("r(varlist)") == "z2 z1"
    out = run("clear\ncapture noisily renpfix\ndisplay _rc").splitlines()
    assert out[-2:] == ["variable * not found", "111"]


def test_numlabel(run):
    out = run("clear\ncapture noisily numlabel\ndisplay _rc").splitlines()
    assert out[-2:] == ["option add or remove should be specified", "198"]
    run('label define yn 0 "No" 1 "Yes"\nlabel define g 1 "One"\nnumlabel yn, add')
    vl = run.session.data.value_labels
    assert vl["yn"] == {0: "0. No", 1: "1. Yes"} and vl["g"] == {1: "One"}
    run("numlabel yn, add")
    assert vl["yn"][0] == "0. No"
    run('numlabel g, add mask("[#] ")')
    assert vl["g"] == {1: "[1] One"}
    run("numlabel yn, remove")
    assert vl["yn"] == {0: "No", 1: "Yes"}


def test_mkmat_svmat(run):
    out = run("clear\ncapture noisily mkmat\ndisplay _rc").splitlines()
    assert out[-2:] == ["varlist required", "100"]
    out = run("capture noisily svmat\ndisplay _rc").splitlines()
    assert out[-2:] == ["matrix  not found", "111"]
    run('clear\ninput a b str3 id\n1 2 "x"\n3 . "y"\n5 6 "z"\nend')
    run("mkmat a b, matrix(A) rownames(id)")
    A = run.session.matrices["A"]
    assert A.rownames == ["x", "y", "z"] and A.colnames == ["a", "b"]
    assert A.data[2].tolist() == [5, 6]
    run("mkmat a b, matrix(B) nomissing")
    assert run.session.matrices["B"].rows == 2
    run("mkmat a")
    assert run.session.matrices["a"].colnames == ["a"]
    run("matrix C = (1, 2 \\ 3, 4 \\ 5, 6 \\ 7, 8)\nmatrix colnames C = u v\nsvmat C")
    assert run.session.data.nobs == 4 and run.eval("C2[4]") == 8
    run("svmat C, names(col)")
    assert run.eval("u[2] + v[1]") == 5


@pytest.mark.parametrize("arg, result", [
    ('"1/4"', "1 2 3 4"), ('"1(2)7 10"', "1 3 5 7 10"), ('"3 1 2", sort', "1 2 3"),
    ('".5 1.5"', ".5 1.5"),
])
def test_numlist(run, arg, result):
    run(f"numlist {arg}")
    assert run.eval("r(numlist)") == result


@pytest.mark.parametrize("arg, rc", [
    ("", 122), ('"1/3", max(2)', 123), ('"3 1", ascending', 124), ('"1/3", range(>=2)', 125),
    ('"1.5", integer', 126), ('"1 ."', 127), ('"a"', 121),
])
def test_numlist_errors(run, arg, rc):
    assert run(f"capture numlist {arg}\ndisplay _rc").strip() == str(rc)


def test_rcof(run):
    out = run("capture noisily rcof\ndisplay _rc").splitlines()
    assert out[-3:] == ["rcof:  _rc  *NOT TRUE* from", "       _rc == 0", "9"]
    assert run('rcof "noisily display 1/" == 198\ndisplay "ok"').splitlines()[-1] == "ok"
    assert run('capture rcof "display 1" != 0\ndisplay _rc').strip() == "9"


def test_query_memory_and_set_defaults(run):
    run("set linesize 80")
    lines = run("query memory").splitlines()
    assert lines[0] == "-" * 80 and lines[1] == "    Memory settings"
    assert lines[2] == "      set maxvar           5000       2048-32767; max. vars allowed"
    assert lines[-1] == "      set segmentsize      32m        1m-32g"
    out = run("capture noisily set_defaults\ndisplay _rc").splitlines()
    assert out[-3:-1] == ["set category required", "    Nothing was found where memory, output, "
                          "interface, efficiency, network, update, trace, mata, unicode, other, "
                          "or _all expected."]
    run("set level 90\nset eolchar mac\nset_defaults output")
    assert run.eval("c(level)") == 95 and run.eval("c(eolchar)") == "unix"


def test_creturn_list_layout(run):
    run("set linesize 80")
    lines = run("creturn list").splitlines()
    assert lines[1] == "System values" and lines[3] == "    " + "-" * 76
    assert "           c(rmsg_time) = 0                          (seconds, from set rmsg)" in lines
    assert '            c(matalibs) = "lmatabase;lmatamc.."      (set matalibs)' in lines
    assert "Trace (program debugging) settings" in lines
    assert "         c(segmentsize) = 33554432                   (set segmentsize)" in lines


def test_fvset_base(run):
    assert run("clear\ncapture noisily fvset\ndisplay _rc").strip() == "0"
    run("clear\nset obs 30\ngen g = mod(_n, 3) + 1\ngen y = g * 2 + mod(_n * 7, 5)")
    run("fvset base 2 g\nquietly regress y i.g")
    assert run.eval("_b[2.g]") == 0 and run.eval("_b[1.g]") != 0
    assert "2" in run("fvset report")
    run("fvset base last g\nquietly regress y i.g")
    assert run.eval("_b[3.g]") == 0
    run("quietly regress y ib1.g")
    assert run.eval("_b[1.g]") == 0
    run("fvset clear g\nquietly regress y i.g")
    assert run.eval("_b[1.g]") == 0


def test_bsample(run):
    assert run("clear\ncapture noisily bsample\ndisplay _rc").strip() == "0"
    run("clear\nset obs 20\ngen id = _n\ngen h = id > 10\ngen c = ceil(id / 2)\nset seed 1")
    run("preserve\nbsample")
    assert run.session.data.nobs == 20
    vals = run.session.data.get("id").data
    assert set(vals) <= set(range(1, 21)) and len(set(vals)) < 20
    run("restore\npreserve\nbsample 5, strata(h)")
    h = run.session.data.get("h").data
    assert len(h) == 10 and h.sum() == 5
    run("restore\npreserve\nbsample, cluster(c) idcluster(nc)")
    ds = run.session.data
    assert ds.nobs == 20 and sorted(set(ds.get("nc").data)) == list(range(1, 11))
    run("restore\ngen w = .\nbsample, weight(w)")
    assert run.session.data.nobs == 20 and run.session.data.get("w").data.sum() == 20
    assert run("capture bsample 30\ndisplay _rc").strip() == "498"


def test_hexdump(run, tmp_path, monkeypatch):
    out = run("capture noisily hexdump\ndisplay _rc").splitlines()
    assert out[-2:] == ["invalid file specification", "198"]
    monkeypatch.chdir(tmp_path)
    (tmp_path / "a.txt").write_bytes(b"Ab 1\r\nx\ty\n\xe9")
    out = run("hexdump a.txt")
    assert "4162 2031 0d0a 7809 790a e9" in out and "| Ab 1..x.y.." in out
    run("hexdump a.txt, results")
    r = {k: run.eval(f"r({k})") for k in ("Windows", "Unix", "uc", "lc", "digit", "tab",
                                          "blank", "extended", "filesize")}
    assert r == {"Windows": 1, "Unix": 1, "uc": 1, "lc": 3, "digit": 1, "tab": 1, "blank": 1,
                 "extended": 1, "filesize": 11}
    assert "Line-end characters" in run("hexdump a.txt, analyze")


def test_sdtesti(run):
    out = run("capture noisily sdtesti\ndisplay _rc").splitlines()
    assert out[-2:] == ["either four arguments for a one-sample test or six arguments for a "
                        "two-sample test must be specified", "198"]
    run("sdtesti 20 5 2 1.5")
    assert run.eval("r(chi2)") == pytest.approx(19 * 4 / 2.25)
    run("sdtesti 20 5 2 15 4 3")
    assert run.eval("r(F)") == pytest.approx(4 / 9) and run.eval("r(df_2)") == 14
    run("clear\nset obs 20\ngen x = mod(_n * 7, 11)\ngen g = _n > 12")
    out = run("sdtest x, by(g)")
    assert "Variance ratio test" in out and run.eval("r(df_1)") == 11
