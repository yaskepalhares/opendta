"""egen: funções de grupo e de linha, conferidas com numpy."""

import numpy as np
import pytest

MISS = 8.98846567431158e307


def col(run, name):
    return np.array(run.session.data.get(name).data)


def strcol(run, name):
    return list(run.session.data.get(name).raw)


@pytest.fixture
def base(run):
    run("""set obs 10
gen g = mod(_n, 3)
gen x = _n * 1.5
replace x = . in 4
gen y = 10 - _n
replace y = . in 7
""")
    return run


def _groups(run):
    g = col(run, "g")
    x = col(run, "x")
    return g, x, x < MISS


def test_group_stats_match_numpy(base):
    run = base
    run("egen m = mean(x), by(g)\negen s = sd(x), by(g)\negen t = total(x), by(g)\n"
        "egen c = count(x), by(g)\negen lo = min(x), by(g)\negen hi = max(x), by(g)")
    g, x, ok = _groups(run)
    for k in (0, 1, 2):
        sel = (g == k) & ok
        v = x[sel]
        assert np.allclose(col(run, "m")[g == k], v.mean(), rtol=1e-6)
        assert np.allclose(col(run, "s")[g == k], v.std(ddof=1), rtol=1e-6)
        assert np.allclose(col(run, "t")[g == k], v.sum(), rtol=1e-6)
        assert (col(run, "c")[g == k] == len(v)).all()
        assert (col(run, "lo")[g == k] == np.float32(v.min())).all()
        assert (col(run, "hi")[g == k] == np.float32(v.max())).all()
    assert run.session.data.get("c").vtype == "float"      # tipo padrão (Stata 14)


def test_by_prefix_equals_by_option(base):
    run = base
    run("egen a = mean(x), by(g)\nbysort g: egen b = mean(x)")
    assert np.array_equal(col(run, "a"), col(run, "b"))


def test_median_pctile_iqr(base):
    run = base
    run("egen md = median(x)\negen p25 = pctile(x), p(25)\negen q = iqr(x)")
    x = np.sort(col(run, "x")[col(run, "x") < MISS])   # 9 valores
    assert col(run, "md")[0] == pytest.approx(x[4], rel=1e-6)
    # P = 9·0,25 = 2,25 → 3º menor valor
    assert col(run, "p25")[0] == pytest.approx(x[2], rel=1e-6)
    assert col(run, "q")[0] == pytest.approx(x[6] - x[2], rel=1e-6)


def test_mode_ties(run):
    run("set obs 6\ngen v = 1 in 1/2\nreplace v = 2 in 3/4\nreplace v = 3 in 5")
    out = run("egen a = mode(v)\negen b = mode(v), minmode\negen c = mode(v), maxmode")
    assert "6 missing values generated" in out
    assert col(run, "b")[0] == 1 and col(run, "c")[0] == 2


def test_rank_variants(run):
    run("set obs 5\ngen v = 3\nreplace v = 1 in 2\nreplace v = 5 in 4\n"
        "egen r = rank(v)\negen f = rank(v), field\negen t = rank(v), track\negen u = rank(v), unique")
    # v = 3 1 3 5 3
    assert list(col(run, "r")) == [3, 1, 3, 5, 3]
    assert list(col(run, "f")) == [2, 5, 2, 1, 2]
    assert list(col(run, "t")) == [2, 1, 2, 5, 2]
    assert sorted(col(run, "u")) == [1, 2, 3, 4, 5]


def test_std_and_skew(base):
    run = base
    run("egen z = std(x)\negen z2 = std(x), mean(50) sd(10)\negen k = skew(x)")
    x = col(run, "x")
    ok = x < MISS
    v = x[ok]
    assert np.allclose(col(run, "z")[ok], (v - v.mean()) / v.std(ddof=1), atol=1e-6)
    assert np.allclose(col(run, "z2")[ok], 50 + 10 * (v - v.mean()) / v.std(ddof=1), atol=1e-4)
    from scipy.stats import skew
    assert col(run, "k")[0] == pytest.approx(skew(v), rel=1e-5)
    assert (col(run, "z")[~ok] >= MISS).all()


def test_group_tag_and_label(run):
    run("""set obs 6
gen str1 s = "b"
replace s = "a" in 2/3
replace s = "" in 6
gen n = _n > 3
egen gid = group(s n), label
egen tg = tag(s)
egen gm = group(s), missing
""")
    # (b,0)=2 (a,0)=1 (a,0)=1 (b,1)=3 (b,1)=3 (missing)=.
    assert list(col(run, "gid")[:5]) == [2, 1, 1, 3, 3]
    assert col(run, "gid")[5] >= MISS
    assert run.session.data.value_labels["gid"][3] == "b 1"
    assert list(col(run, "tg")) == [1, 1, 0, 0, 0, 0]
    assert list(col(run, "gm")) == [3, 2, 2, 3, 3, 1]


def test_seq_fill_cut(run):
    run("""set obs 7
gen g = _n > 4
egen a = seq(), by(g)
egen b = seq(), from(1) to(3)
egen c = seq(), block(2)
egen d = fill(5 10)
egen e = fill(1 1 2 2)
egen f = fill(7 8 9 7 8 9)
gen x = _n
egen ct = cut(x), at(1 3 6)
egen ic = cut(x), at(1 3 6) icodes
egen gr = cut(x), group(2)
""")
    assert list(col(run, "a")) == [1, 2, 3, 4, 1, 2, 3]
    assert list(col(run, "b")) == [1, 2, 3, 1, 2, 3, 1]
    assert list(col(run, "c")) == [1, 1, 2, 2, 3, 3, 4]
    assert list(col(run, "d")) == [5, 10, 15, 20, 25, 30, 35]
    assert list(col(run, "e")) == [1, 1, 2, 2, 3, 3, 4]
    assert list(col(run, "f")) == [7, 8, 9, 7, 8, 9, 7]
    ct = col(run, "ct")
    assert list(ct[:5]) == [1, 1, 3, 3, 3] and (ct[5:] >= MISS).all()
    assert list(col(run, "ic")[:5]) == [0, 0, 1, 1, 1]
    assert list(col(run, "gr")) == [0, 0, 0, 1, 1, 1, 1]


def test_row_functions(base):
    run = base
    run("""egen rm = rowmean(x y)
egen rt = rowtotal(x y)
egen rtm = rowtotal(x y), missing
egen rn = rownonmiss(x y)
egen rmi = rowmiss(x y)
egen rmin = rowmin(x y)
egen rmax = rowmax(x y)
egen rf = rowfirst(x y)
egen rl = rowlast(x y)
egen rs = rowsd(x y g)
""")
    x, y, g = col(run, "x"), col(run, "y"), col(run, "g")
    X = np.column_stack([x, y])
    Xn = np.where(X >= MISS, np.nan, X)
    assert np.allclose(col(run, "rm"), np.nanmean(Xn, 1), rtol=1e-6)
    assert np.allclose(col(run, "rt"), np.nansum(Xn, 1), rtol=1e-6)
    assert list(col(run, "rn")) == list((~np.isnan(Xn)).sum(1))
    assert list(col(run, "rmi")) == list(np.isnan(Xn).sum(1))
    assert np.allclose(col(run, "rmin"), np.nanmin(Xn, 1), rtol=1e-6)
    assert np.allclose(col(run, "rmax"), np.nanmax(Xn, 1), rtol=1e-6)
    assert col(run, "rf")[3] == pytest.approx(y[3]) and col(run, "rl")[6] == pytest.approx(x[6])
    Z = np.column_stack([x, y, g])
    Zn = np.where(Z >= MISS, np.nan, Z)
    assert np.allclose(col(run, "rs"), np.nanstd(Zn, 1, ddof=1), rtol=1e-6)
    assert run.session.data.get("rn").vtype == "float"


def test_rowtotal_missing_option(run):
    run("set obs 2\ngen a = 1 in 1\ngen b = 2 in 1\n"
        "egen t = rowtotal(a b)\negen tm = rowtotal(a b), missing")
    assert list(col(run, "t")) == [3, 0]
    assert col(run, "tm")[1] >= MISS


def test_any_functions(run):
    run("""set obs 4
gen a = _n
gen b = 5 - _n
egen n = anycount(a b), values(1 2)
egen m = anymatch(a b), values(4)
egen v = anyvalue(a), values(2 3)
""")
    assert list(col(run, "n")) == [1, 1, 1, 1]
    assert list(col(run, "m")) == [1, 0, 0, 1]
    v = col(run, "v")
    assert v[1] == 2 and v[2] == 3 and v[0] >= MISS and v[3] >= MISS


def test_diff(run):
    run("set obs 3\ngen a = 1\ngen b = _n\negen d = diff(a b)")
    assert list(col(run, "d")) == [0, 1, 1]


def test_concat_and_ends(run):
    run("""set obs 2
gen n = _n
gen str10 s = "ab cd ef"
label define yn 1 "yes" 2 "no"
label values n yn
egen c = concat(s n), punct(_)
egen cd = concat(n s), decode punct(" ")
egen h = ends(s)
egen l = ends(s), last
egen t = ends(s), tail
""")
    assert strcol(run, "c") == ["ab cd ef_1", "ab cd ef_2"]
    assert strcol(run, "cd") == ["yes ab cd ef", "no ab cd ef"]
    assert strcol(run, "h")[0] == "ab"
    assert strcol(run, "l")[0] == "ef"
    assert strcol(run, "t")[0] == "cd ef"


def test_if_in_and_errors(base):
    run = base
    out = run("egen m = mean(x) if g == 1")
    m = col(run, "m")
    g = col(run, "g")
    assert (m[g != 1] >= MISS).all()
    assert "missing values generated" in out
    run("egen m = mean(x)")
    assert run.rc == 110
    run("egen q = nada(x)")
    assert run.rc == 133
    run("egen q = mean(x")
    assert run.rc == 198
