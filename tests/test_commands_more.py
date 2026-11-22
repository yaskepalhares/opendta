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
    # compat 2008: return list from the last result; the values are floats; not after logit
    run("quietly regress y x1 x2 x3\nquietly estat vif")
    words = run("return list").split()
    assert [w for w in words if w.startswith("r(vif")] == ["r(vif_3)", "r(vif_2)", "r(vif_1)"]
    assert run.eval("r(vif_1)") == pytest.approx(8.756916999816895, abs=1e-12)
    run("gen yb = y > 15\nquietly logit yb x1")
    out = run("capture noisily estat vif\ndisplay _rc").splitlines()
    assert out[-2:] == ["subcommand estat vif is unrecognized", "321"]


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
    run("clear\nset obs 1\ngen c = 1\ngen a = 2\ngen B = 3\ngen x10 = 4\ngen x2 = 5\n"
        "rename x10 x1")
    run("gen x10 = 6\naorder")
    # compat 2008: upper case first, digits compared as numbers
    assert _vars(run) == ["B", "a", "c", "x1", "x2", "x10"]
    run("drop x10\nmove x2 a")
    assert _vars(run) == ["B", "x2", "a", "c", "x1"]
    assert run("capture noisily move x2\ndisplay _rc").splitlines()[-1] == "102"
    run("renpfix x z")
    assert _vars(run) == ["B", "z2", "a", "c", "z1"]
    assert run.eval("r(varlist)") == "z2 z1"
    out = run("clear\ncapture noisily renpfix\ndisplay _rc").splitlines()
    assert out[-2:] == ["variable * not found", "111"]


def test_numlabel(run):
    out = run("clear\ncapture noisily numlabel\ndisplay _rc").splitlines()
    assert out[-2:] == ["option add or remove should be specified", "198"]
    run('label define yn 0 "No" 1 "Yes"\nlabel define g 1 "One"\nnumlabel yn, add')
    vl = run.session.data.value_labels
    assert vl["yn"] == {0: "0. No", 1: "1. Yes"} and vl["g"] == {1: "One"}
    for opt in ("add", "add force"):
        assert run(f"numlabel yn, {opt}").strip() == "(no value label to be modified)"
        assert vl["yn"][0] == "0. No"
    out = run("capture noisily numlabel yn, add remove\ndisplay _rc").splitlines()
    assert out[-2:] == ["options add and remove must not be both specified", "198"]
    run('numlabel g, add mask("[#] ")')
    assert vl["g"] == {1: "[1] One"}
    run("numlabel yn, remove")
    assert vl["yn"] == {0: "No", 1: "Yes"}


def test_mkmat_svmat_names(run):
    # compat 2008: default row names r1 r2..., rowprefix() replaces the r; eqcol joins directly
    run("clear\ninput a b\n1 2\n3 .\n5 6\nend")
    run("mkmat a b, matrix(B) obs\nmkmat a b, matrix(C) nomissing rowprefix(p)")
    assert "r1" in run("matrix list B") and run("matrix list C").splitlines()[-1].startswith("p2")
    run("matrix D = (1, 2 \\ 3, 4)\nmatrix colnames D = u v\nmatrix coleq D = e1 e2")
    run("svmat D, names(eqcol)\nsvmat D, names(matcol)")
    assert _vars(run)[-4:] == ["e1u", "e2v", "Du", "Dv"]


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
    out = run("capture display 1/\ncapture noisily rcof\ndisplay _rc").splitlines()
    assert out[-4:] == ["rcof:  _rc  *NOT TRUE* from", "       ", "       _rc == 0", "9"]
    out = run('capture noisily rcof "display 1/" == 0').splitlines()
    assert out[-3:] == ["rcof:  _rc  == 0 *NOT TRUE* from", "       display 1/",
                        "       _rc == 198"]
    run('numlist "1/2000"')
    assert len(run.eval("r(numlist)").split()) == 2000
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
    out = run("set level 90\nset eolchar mac\nset_defaults output").splitlines()
    assert out[:4] == ["-> set more on", "-> set rmsg off", "-> set dp period",
                       "-> (linesize left at current value)"]
    assert out[-2:] == ["-> set logtype smcl", "(preferences reset)"]
    # eolchar is not among the settings set_defaults output resets (compat 2008)
    assert run.eval("c(level)") == 95 and run.eval("c(eolchar)") == "mac"
    out = run("capture noisily set_defaults foo\ndisplay _rc").splitlines()
    assert out[-2].startswith('    foo found where" memory, output,') and out[-1] == "198"


def test_query_output_layout(run):
    run("set linesize 80")
    lines = run("query output").splitlines()
    assert lines[:4] == ["-" * 80, "    Output settings", "      set more             off",
                         "      set rmsg             off"]
    assert "      set dp               period     may be period or comma" in lines
    assert lines[7] == "      " + "-" * 74
    assert run("set httpproxyhost abc").strip() == "(set httpproxyhost preference recorded)"


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
    # compat 2008: columns as wide as their entries plus 4, between rules
    assert run("fvset report").splitlines()[1:] == [
        "  Variable    Base    Design", "  " + "-" * 26, "  g           2", "  " + "-" * 26]
    out = run("capture noisily fvset base bad g\ndisplay _rc").splitlines()
    assert out[-2:] == ["invalid base specification", "198"]
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
    # compat 2008: the layouts and r() of Stata 14
    (tmp_path / "b.txt").write_bytes(b"Ab 1\nx\ty\n")
    out = run("hexdump b.txt").splitlines()
    assert out[2] == "         address |  0 1  2 3  4 5  6 7  8 9  a b  c d  e f | 0123456789abcdef"
    assert out[4] == "               0 | 4162 2031 0a78 0979 0a                  | Ab 1.x.y."
    out = run("hexdump b.txt, analyze").splitlines()
    assert out[2] == ("    \\r\\n         (Windows)              0      minimum"
                      "                        3")
    assert out[-1] == "     \\t \\n blank 1 A b x y"
    assert run.eval("r(lnum)") == 2 and run.eval("r(eoleof)") == 1 and run.eval("r(l1)") == 4
    assert run.eval("r(format)") == "ASCII"


def test_sdtesti(run):
    out = run("capture noisily sdtesti\ndisplay _rc").splitlines()
    assert out[-2:] == ["either four arguments for a one-sample test or six arguments for a "
                        "two-sample test must be specified", "198"]
    run("sdtesti 20 5 2 1.5")
    assert run.eval("r(chi2)") == pytest.approx(19 * 4 / 2.25)
    out = run("return list").split("scalars:")[1].split()
    assert [w for w in out if w.startswith("r(")] == ["r(p)", "r(p_l)", "r(p_u)", "r(sd)",
                                                     "r(chi2)", "r(df)", "r(N)"]
    out = run("sdtesti 20 5 2 15 4 3").splitlines()
    # the immediate form has no heading over the names (compat 2008)
    assert any(line.startswith("         |     Obs") for line in out)
    assert run.eval("r(F)") == pytest.approx(4 / 9) and run.eval("r(df_2)") == 14
    run("clear\nset obs 20\ngen x = mod(_n * 7, 11)\ngen g = _n > 12")
    out = run("sdtest x, by(g)")
    assert "Variance ratio test" in out and run.eval("r(df_1)") == 11


def test_robvar(run):
    import numpy as np
    from scipy import stats
    out = run("clear\ncapture noisily robvar\ndisplay _rc").splitlines()
    assert out[-2:] == ["varlist required", "100"]
    run("clear\nset obs 40\ngen g = mod(_n, 3)\ngen x = mod(_n * 7, 13) * (1 + g)")
    out = run("robvar x, by(g)")
    assert "Summary of x" in out and "df(2, 37)" in out
    assert "W0  =  7.5583963   df(2, 37)     Pr > F = 0.00176866" in out   # compat 2008
    x = np.array(run.session.data.get("x").data)
    g = np.array(run.session.data.get("g").data)
    groups = [x[g == k] for k in range(3)]
    assert run.eval("r(w0)") == pytest.approx(stats.levene(*groups, center="mean").statistic)
    assert run.eval("r(w50)") == pytest.approx(stats.levene(*groups, center="median").statistic)


def test_sfrancia(run):
    import numpy as np
    from scipy import stats
    out = run("clear\ncapture noisily sfrancia\ndisplay _rc").splitlines()
    assert out[-2:] == ["varlist required", "100"]
    run("clear\nset obs 50\ngen x = invnormal((_n - .5) / 50) + (_n > 45) * 3")
    out = run("sfrancia x")
    assert "Shapiro-Francia W' test for normal data" in out
    x = np.sort(np.array(run.session.data.get("x").data))
    m = stats.norm.ppf((np.arange(1, 51) - .375) / 50.25)
    W = np.corrcoef(m, x)[0, 1] ** 2
    assert run.eval("r(W)") == pytest.approx(W)
    assert run.eval("r(z)") > 0 and run.eval("r(p)") < .05
    # with ties, Blom's scores of the average ranks (compat 2008: W' = 0.94865)
    run("clear\nset obs 40\ngen g = mod(_n, 3)\ngen x = mod(_n * 7, 13) * (1 + g)")
    assert "        40    0.94865      2.250     1.512    0.06526" in run("sfrancia x")
    out = run("capture noisily sfrancia x if _n < 5\ndisplay _rc").splitlines()
    assert "           x |         4          .          .         .          ." in out
    assert out[-1] == "0"


@pytest.mark.parametrize("cmd", ["qqplot", "qladder"])
def test_quantile_graphs_need_varlist(run, cmd):
    out = run(f"clear\ncapture noisily {cmd}\ndisplay _rc").splitlines()
    assert out[-2:] == ["varlist required", "100"]


def test_qqplot_qladder(run):
    import numpy as np
    from opendta.graphics.baseplots import matched_quantiles
    a, b = matched_quantiles(np.array([3.0, 1, 2]), np.array([10.0, 20, 30, 40, 50, 60]))
    assert a.tolist() == [1, 2, 3] and b.tolist() == [15, 35, 55]
    run("clear\nset obs 30\ngen x = _n\ngen y = sqrt(_n)")
    run("qqplot y x, name(q1)\nqladder x, name(q2)")
    from opendta.graphics.store import store
    assert store(run.session).get("q2").family == "combine"
    assert len(store(run.session).get("q2").children) == 9


@pytest.mark.parametrize("cmd", ["iqreg", "sqreg", "bsqreg"])
def test_quantile_regressions_replay(run, cmd):
    out = run(f"clear\ncapture noisily {cmd}\ndisplay _rc").splitlines()
    assert out[-2:] == ["last estimates not found", "301"]


def test_quantile_regressions(run):
    run("clear\nset obs 60\ngen x = mod(_n*7,13)\ngen z = mod(_n*5,7)\n"
        "gen y = 2 + x + .5*z + mod(_n*11,9)\nset seed 1")
    run("quietly qreg y x z")
    b50 = run.eval("_b[x]")
    run("quietly qreg y x z, q(.25)")
    b25 = run.eval("_b[x]")
    run("quietly qreg y x z, q(.75)")
    b75 = run.eval("_b[x]")
    out = run("bsqreg y x z, reps(10)")
    assert "Median regression, bootstrap(10) SEs" in out and run.eval("_b[x]") == pytest.approx(b50)
    assert run.eval("_se[x]") > 0 and run.eval("e(reps)") == 10
    run("quietly iqreg y x z, reps(10)")
    assert run.eval("_b[x]") == pytest.approx(b75 - b25)
    run("quietly sqreg y x z, q(.25 .75) reps(10)")
    assert run.eval("_b[q25:x]") == pytest.approx(b25) and run.eval("_b[q75:x]") == pytest.approx(b75)
    assert run.eval("e(n_q)") == 2 and run.eval('e(eqnames) == "q25 q75"') == 1
    # compat 2008: names, labels and e() of Stata 14
    assert "(fitting base model)" not in run("iqreg y x z, reps(10) nodots")
    assert run.eval("e(sumrdev0)") == pytest.approx(82.625) and run.eval('e(predict) == "qreg_p"')
    out = run("sqreg y x z, q(.1 .333) reps(10) nodots").splitlines()
    assert "q100         |" in out and "q333         |" in out
    assert any(".10 Pseudo R2 =" in line for line in out)
    assert run('ereturn list').count("e(vce)") == 0


def test_input_leading_decimal_point(run):
    run("clear\ninput x\n.5\n-.25\n.a\nend")
    assert list(run.session.data.get("x").data[:2]) == [.5, -.25]
    assert run.eval("x[3] == .a") == 1


_PK = ("clear\ninput id t c\n" + "\n".join(
    f"{i} {t} {c}" for i, rows in ((1, [(0, 0), (.5, 3.1), (1, 4.8), (2, 4.2), (4, 3), (8, 1.6), (12, .8)]),
                                   (2, [(0, 0), (.5, 2.5), (1, 4.1), (2, 4.5), (4, 3.4), (8, 1.9), (12, 1.1)]))
    for t, c in rows) + "\nend")


@pytest.mark.parametrize("cmd", ["pkexamine", "pksumm", "pkshape", "pkcross", "pkequiv",
                                 "pkcollapse"])
def test_pk_needs_varlist(run, cmd):
    out = run(f"clear\ncapture noisily {cmd}\ndisplay _rc").splitlines()
    assert out[-2:] == ["varlist required", "100"]


def test_pkexamine_measures(run):
    import math
    import numpy as np
    run(_PK)
    run("pkexamine t c if id == 1, trapezoid")
    t = np.array([0, .5, 1, 2, 4, 8, 12])
    c = np.array([0, 3.1, 4.8, 4.2, 3, 1.6, .8])
    assert run.eval("r(auc)") == pytest.approx(np.sum((c[1:] + c[:-1]) / 2 * np.diff(t)), rel=1e-6)
    b = np.polyfit(t[-3:], np.log(c[-3:]), 1)[0]
    assert run.eval("r(ke)") == pytest.approx(-b, rel=1e-6)
    assert run.eval("r(half)") == pytest.approx(math.log(2) / -b, rel=1e-6)
    assert run.eval("r(cmax)") == pytest.approx(4.8) and run.eval("r(tomc)") == 1
    assert run.eval("r(auc_ln) - r(auc)") == pytest.approx(.8 / -b, rel=1e-6)
    # compat 2008: auc kept as a float, the line triangle from the last concentration, layout
    run("pkexamine t c if id == 1, trapezoid fit(4)")
    assert run.eval("r(auc)") == pytest.approx(28.45000076293945, abs=1e-12)
    assert run.eval("r(auc_line)") == pytest.approx(29.13979678642299, rel=1e-9)
    out = run("pkexamine t c if id == 1 & t > 0").splitlines()
    assert "    Warning: the point (0, 0) is not in your data.  It will be added." in out
    assert " " * 42 + "Maximum concentration =       4.8" in out
    assert "auclog" in run("pksumm id t c")
    run("pkcollapse t c, id(id) stat(auc cmax)")
    assert [v.name for v in run.session.data.vars] == ["id", "auc_c", "cmax_c"]
    assert run.session.data.nobs == 2


def test_pkshape_pkcross_pkequiv(run):
    run("clear\ninput id seq p1 p2\n1 1 10 12\n2 1 11 13\n3 1 9 11\n4 2 12 10\n5 2 13 12\n"
        "6 2 11 10\nend\npkshape id seq p1 p2, order(ab ba)")
    ds = run.session.data
    assert ds.nobs == 12 and [v.name for v in ds.vars] == ["id", "sequence", "outcome", "treat",
                                                           "carry", "period"]
    # compat 2008: the first period of every subject, then the second
    assert list(ds.get("period").data[5:7]) == [1, 2] and list(ds.get("treat").data[6:10]) == [2, 2, 2, 1]
    assert list(ds.get("carry").data[9:]) == [2, 2, 2]
    out = run("pkcross outcome")
    assert "Treatment effect" in out and "Carryover effect" not in out
    assert "    Omnibus measure of separability of treatment and carryover =   29.2893%" in out
    out = run("pkequiv outcome treat period sequence id")
    assert "      Note: Reference treatment = 1" in out
    assert run.eval("r(lci)") == pytest.approx(1.311358867621937, rel=1e-8)
    assert run.eval("r(uci)") == pytest.approx(2.021974465711403, rel=1e-8)
    assert run.eval("r(stddev)") == pytest.approx(.2041241459923532, rel=1e-8)


def test_forecast(run):
    out = run("clear\ncapture noisily forecast\ndisplay _rc").splitlines()
    assert out[-3] == "missing forecast subcommand" and out[-1] == "198"
    assert out[-2].startswith("    forecast subcommands are create, estimates, identity")
    run("clear\nset obs 30\ngen t = _n\ntsset t\ngen x = mod(_n*7, 11)\ngen y = 5 if _n == 1\n"
        "replace y = 2 + .5*L.y + .3*x + mod(_n*3,5)/10 if _n > 1\nreplace y = . if _n > 25\n"
        "quietly regress y L.y x\nestimates store m1")
    b = {k: run.eval(f"_b[{k}]") for k in ("L.y", "x", "_cons")}
    run("forecast create m\nforecast estimates m1\nforecast identity z = y + x, generate\n"
        "forecast exogenous x\nforecast solve")
    d = run.session.data
    fy, x, y = d.get("f_y").data, d.get("x").data, d.get("y").data
    assert fy[25] == pytest.approx(b["_cons"] + b["L.y"] * y[24] + b["x"] * x[25], rel=1e-6)
    assert fy[26] == pytest.approx(b["_cons"] + b["L.y"] * fy[25] + b["x"] * x[26], rel=1e-6)
    assert d.get("f_z").data[26] == pytest.approx(fy[26] + x[26], rel=1e-6)
    assert fy[24] >= 8.9e307
    run("forecast adjust y = y + 1\nforecast solve, prefix(a_) log(off)")
    assert d.get("a_y").data[25] == pytest.approx(fy[25] + 1, rel=1e-6)
    run("forecast solve, prefix(s_) static begin(20) log(off)")
    sy = run.session.data.get("s_y").data
    assert sy[20] == pytest.approx(b["_cons"] + b["L.y"] * y[19] + b["x"] * x[20] + 1, rel=1e-6)
    assert run("forecast query").strip() == "Forecast model m has been started."
    run("forecast drop\nforecast clear")
    assert run("forecast query").strip() == "No forecast model has been started."
