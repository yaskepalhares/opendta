"""collapse, contract, expand, fillin, append, merge, joinby, cross e reshape."""

import numpy as np
import pytest

MISS = 8.98846567431158e307


def col(run, name):
    return list(np.array(run.session.data.get(name).data))


def names(run):
    return run.session.data.names


@pytest.fixture
def files(run, tmp_path):
    d = tmp_path.as_posix()
    run(f"""clear
input id str5 nome renda
1 "ana" 10
2 "bia" .
3 "caio" 30
5 "eva" 50
end
label define sim 1 "um"
label values id sim
save "{d}/m"
clear
input id renda byte idade
1 11 20
2 22 30
4 44 40
end
save "{d}/u"
""")
    assert run.rc == 0
    return d


# -- collapse ---------------------------------------------------------------

def test_collapse_stats_by_group(run):
    run("""set obs 9
gen g = mod(_n, 3)
gen x = _n * 2
replace x = . in 4
gen str3 s = "a" + string(_n)
collapse (mean) mx=x (sd) sx=x (count) n=x (median) md=x (first) s (max) x (sum) t=x, by(g)
""")
    assert run.rc == 0
    assert names(run) == ["g", "mx", "sx", "n", "md", "s", "x", "t"]
    assert col(run, "mx") == [12, 8, 10]
    assert col(run, "sx")[1] == pytest.approx(np.std([2, 14], ddof=1), rel=1e-6)
    assert col(run, "n") == [3, 2, 3]
    assert col(run, "t") == [36, 16, 30]
    assert list(run.session.data.get("s").raw) == ["a3", "a1", "a2"]
    ds = run.session.data
    assert ds.get("n").vtype == "long" and ds.get("mx").vtype == "float"
    assert ds.get("mx").label == "(mean) x"
    assert ds.sortlist == ["g"]
    out = run("list in 1")
    assert run.rc == 0 and out


def test_collapse_weights_and_percentiles(run):
    run("""set obs 6
gen g = _n > 3
gen x = _n
gen w = _n
collapse (mean) m=x (sum) t=x (count) n=x (p25) q=x [fw=w], by(g)
""")
    x = np.arange(1, 7)
    w = x
    assert col(run, "m")[0] == pytest.approx((x[:3] * w[:3]).sum() / w[:3].sum(), rel=1e-6)
    assert col(run, "t") == [(x[:3] * w[:3]).sum(), (x[3:] * w[3:]).sum()]
    assert col(run, "n") == [6, 15]
    # g=0: valores 1,2,2,3,3,3 → P = 6·0,25 = 1,5 → 2º valor = 2
    assert col(run, "q")[0] == 2


def test_collapse_aweight_sum_normalized(run):
    run("set obs 2\ngen x = _n\ngen w = cond(_n==1, 1, 3)\ncollapse (sum) x (mean) m=x [aw=w]")
    # pesos normalizados: 0,5 e 1,5 → soma = 0,5·1 + 1,5·2
    assert col(run, "x") == [3.5]
    assert col(run, "m")[0] == pytest.approx(1.75)


def test_collapse_errors(run):
    run("set obs 2\ngen x = _n\ngen str1 s = \"a\"")
    run("collapse (mean) s")
    assert run.rc == 109
    run("collapse (bogus) x")
    assert run.rc == 198
    run("collapse x if x > 5")
    assert run.rc == 2000


# -- contract / expand / fillin -----------------------------------------------

def test_contract(run):
    run("""set obs 6
gen a = _n > 2
gen b = mod(_n, 2)
replace b = 0 in 3
contract a b, percent(p) cfreq(cf) zero
""")
    # combinações (a,b): (0,1) (0,0) (1,0) (1,1) (1,0) (1,0)
    assert names(run) == ["a", "b", "_freq", "cf", "p"]
    assert col(run, "_freq") == [1, 1, 3, 1]
    assert col(run, "cf") == [1, 2, 5, 6]
    assert col(run, "p")[2] == pytest.approx(50)
    run("clear\nset obs 3\ngen a = 1\ngen b = _n\ncontract a b, zero freq(n)")
    assert col(run, "n") == [1, 1, 1]


def test_expand(run):
    run("set obs 3\ngen k = _n - 1")
    out = run("expand k, generate(novo)")
    assert out == "(1 observation created)\n"     # k < 1 mantém a obs. sem cópias
    assert col(run, "k") == [0, 1, 2, 2]
    assert col(run, "novo") == [0, 0, 0, 1]
    out = run("expand 3 if k == 1")
    assert "(2 observations created)" in out
    assert run.session.data.nobs == 6


def test_fillin(run):
    run("set obs 3\ngen a = _n\ngen b = _n > 1\ngen x = 10 * _n\nfillin a b")
    assert col(run, "a") == [1, 1, 2, 2, 3, 3]
    assert col(run, "b") == [0, 1, 0, 1, 0, 1]
    assert col(run, "_fillin") == [0, 1, 1, 0, 1, 0]
    assert col(run, "x")[1] >= MISS
    run("fillin a b")
    assert run.rc == 110


# -- append -------------------------------------------------------------------

def test_append(run, files):
    run(f'use "{files}/m", clear\nappend using "{files}/u", generate(fonte)')
    assert run.rc == 0
    assert names(run) == ["id", "nome", "renda", "idade", "fonte"]
    assert col(run, "fonte") == [0, 0, 0, 0, 1, 1, 1]
    assert col(run, "renda")[4:] == [11, 22, 44]
    assert list(run.session.data.get("nome").raw)[4:] == ["", "", ""]
    assert col(run, "idade")[0] >= MISS


def test_append_promotes_and_mismatch(run, tmp_path):
    d = tmp_path.as_posix()
    run(f"""clear
set obs 1
gen byte v = 1
gen str2 s = "ab"
save "{d}/a"
clear
set obs 1
gen v = 1.5
gen str5 s = "abcde"
save "{d}/b"
use "{d}/a", clear
""")
    out = run(f'append using "{d}/b"')
    assert "variable v was byte, now float" in out
    ds = run.session.data
    assert ds.get("v").vtype == "float" and ds.get("s").vtype == "str5"
    run(f'clear\nset obs 1\ngen str1 v = "x"\nappend using "{d}/b"')
    assert run.rc == 106
    run(f'append using "{d}/b", force')
    assert run.rc == 0


# -- merge --------------------------------------------------------------------

def test_merge_1to1_report(run, files):
    run(f'use "{files}/m", clear')
    out = run(f'merge 1:1 id using "{files}/u"')
    expected = (
        "\n    Result                           # of obs.\n"
        "    -----------------------------------------\n"
        "    not matched                             3\n"
        "        from master                         2  (_merge==1)\n"
        "        from using                          1  (_merge==2)\n"
        "\n"
        "    matched                                 2  (_merge==3)\n"
        "    -----------------------------------------\n")
    assert out == expected
    assert col(run, "id") == [1, 2, 3, 4, 5]
    assert col(run, "_merge") == [3, 3, 1, 2, 1]
    assert col(run, "renda")[:4] == [10, MISS, 30, 44]
    assert col(run, "idade")[:2] == [20, 30]
    ds = run.session.data
    assert ds.value_labels["_merge"][3] == "matched (3)"
    assert ds.sortlist == ["id"]


def test_merge_update_replace(run, files):
    run(f'use "{files}/m", clear\nmerge 1:1 id using "{files}/u", update')
    assert col(run, "_merge")[:2] == [5, 4]
    assert col(run, "renda")[:2] == [10, 22]
    run(f'use "{files}/m", clear\nmerge 1:1 id using "{files}/u", update replace')
    assert col(run, "renda")[:2] == [11, 22]


def test_merge_keep_keepusing_nogen(run, files):
    run(f'use "{files}/m", clear')
    out = run(f'merge 1:1 id using "{files}/u", keep(match) keepusing(idade) nogenerate')
    assert "matched                                 2  (_merge==3)" in out
    assert names(run) == ["id", "nome", "renda", "idade"]
    assert col(run, "renda") == [10, MISS]


def test_merge_m1_and_1m(run, tmp_path):
    d = tmp_path.as_posix()
    run(f"""clear
input g v
1 100
2 200
end
save "{d}/look"
clear
input g x
1 1
1 2
2 3
3 4
end
merge m:1 g using "{d}/look", noreport
""")
    assert col(run, "v")[:3] == [100, 100, 200]
    assert col(run, "_merge") == [3, 3, 3, 1]
    run(f'use "{d}/look", clear\ndrop in 2\nsave "{d}/one", replace\n'
        f'clear\ninput g x\n1 1\n1 2\nend\nsave "{d}/many"\nuse "{d}/one", clear\n'
        f'merge 1:m g using "{d}/many", noreport')
    assert col(run, "x") == [1, 2] and col(run, "v") == [100, 100]


def test_merge_errors(run, files, tmp_path):
    d = tmp_path.as_posix()
    run(f'clear\ninput id\n1\n1\nend\nmerge 1:1 id using "{files}/u"')
    assert run.rc == 459
    run(f'merge m:1 id using "{files}/u", assert(match)')
    assert run.rc == 9
    run(f'use "{files}/m", clear\nmerge 1:1 nome using "{files}/u"')
    assert run.rc == 111
    run(f'merge id using "{files}/u"')
    assert run.rc == 198
    run(f'merge 1:1 id using "{d}/nada"')
    assert run.rc == 601


def test_merge_by_n(run, files):
    run(f'use "{files}/m", clear\nmerge 1:1 _n using "{files}/u", noreport')
    # renda vem da master nas linhas casadas
    assert col(run, "renda")[:3] == [10, MISS, 30]
    assert col(run, "_merge") == [3, 3, 3, 1]


# -- joinby / cross -------------------------------------------------------------

def test_joinby_and_cross(run, tmp_path):
    d = tmp_path.as_posix()
    run(f"""clear
input k str3 a
1 "x"
1 "y"
2 "z"
end
save "{d}/j1"
clear
input k b
1 10
1 20
3 30
end
joinby k using "{d}/j1"
""")
    assert run.session.data.nobs == 4
    assert col(run, "b") == [10, 10, 20, 20]
    run(f'clear\ninput k b\n1 10\n3 30\nend\njoinby k using "{d}/j1", unmatched(both) _merge(m)')
    assert col(run, "m") == [3, 3, 2, 1]
    run(f'clear\nset obs 2\ngen c = _n\ncross using "{d}/j1"')
    assert run.session.data.nobs == 6
    assert col(run, "c") == [1, 1, 1, 2, 2, 2]


# -- reshape ------------------------------------------------------------------

_WIDE = """clear
input id sex inc80 inc81 inc82 ue80 ue81
1 0 5000 5500 6000 0 1
2 1 2000 2200 3300 1 0
3 0 3000 2000 1000 0 0
end
"""


def test_reshape_long_output(run):
    run(_WIDE)
    out = run("reshape long inc ue, i(id) j(year)")
    assert out.startswith("(note: j = 80 81 82)\n(note: ue82 not found)\n\n"
                          "Data                               wide   ->   long\n")
    assert "Number of obs.                        3   ->       9\n" in out
    assert "j variable (3 values)                     ->   year\n" in out
    assert "                      inc80 inc81 inc82   ->   inc\n" in out
    assert names(run) == ["id", "year", "sex", "inc", "ue"]
    assert col(run, "year") == [80, 81, 82] * 3
    assert col(run, "inc")[3:6] == [2000, 2200, 3300]
    assert col(run, "ue")[2] >= MISS


def test_reshape_roundtrip_and_string_j(run):
    run(_WIDE + "reshape long inc ue, i(id) j(year)")
    out = run("reshape wide")
    assert "j variable (3 values)              year   ->   (dropped)\n" in out
    assert names(run) == ["id", "inc80", "ue80", "inc81", "ue81", "inc82", "ue82", "sex"]
    assert col(run, "inc81") == [5500, 2200, 2000]
    run("reshape long inc@ ue, i(id) j(t) string")
    assert run.rc == 0
    assert list(run.session.data.get("t").raw[:3]) == ["80", "81", "82"]


def test_reshape_errors(run):
    run(_WIDE + "replace id = 1 in 2\nreshape long inc ue, i(id) j(year)")
    assert run.rc == 9
    run(_WIDE + "reshape long inc ue, i(id) j(year)\nreplace sex = 1 in 1\nreshape wide inc ue, i(id) j(year)")
    assert run.rc == 9
    run(_WIDE + "reshape long inc ue, i(id) j(year)\nreplace year = 80 in 2\nreshape wide inc ue, i(id) j(year)")
    assert run.rc == 9
