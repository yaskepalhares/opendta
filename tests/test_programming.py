"""Block 3, Programming: varlist utilities, mark, postfile, matlist, tabdisp, sortpreserve,
plugin, pause/sleep/window, viewsource, signestimationsample and makecns."""

DATA = """
set obs 20
gen rep78 = 1+mod(_n,5)
gen foreign = mod(_n,2)
gen mpg = _n
gen turn = 30+_n
gen t = _n
tsset t
"""


def test_unab_and_errors(run):
    run(DATA)
    run("unab x : mp tur")
    assert run('display "`x\'"').strip() == "mpg turn"
    assert "variable junk not found" in run("capture noisily unab x : junk")
    assert "too many variables" in run("capture noisily unab x : mpg turn, max(1)")
    assert "too few variables" in run("capture noisily unab x : mpg, min(3)")


def test_tsunab_fvunab(run):
    run(DATA)
    run("tsunab l : l(1/3).mpg")
    assert run('display "`l\'"').strip() == "L.mpg L2.mpg L3.mpg"
    run("fvunab f : i.rep78 mpg rep78#foreign")
    assert run('display "`f\'"').strip() == "i.rep78 mpg i.rep78#i.foreign"


def test_unabcmd(run):
    out = run("unabcmd gen\nreturn list")
    assert '"generate"' in out
    assert "not found" in run("capture noisily unabcmd ka")


def test_fvexpand_fvrevar(run):
    run(DATA)
    run("fvexpand i.rep78 mpg")
    assert run("display r(varlist)").strip() == "1b.rep78 2.rep78 3.rep78 4.rep78 5.rep78 mpg"
    run("fvrevar i(2,3)bn.rep78 mpg, list")
    assert run("display r(varlist)").strip() == "mpg rep78"
    run("tsrevar l.mpg turn")
    assert len(run("display r(varlist)").split()) == 2
    run("tsrevar l.mpg d.mpg turn, list")
    assert run("display r(varlist)").split()[-1] == "turn"


def test_rmcoll(run):
    run(DATA)
    run("gen tt = turn + mpg")
    out = run("_rmcoll turn mpg tt")
    assert "omitted because of collinearity" in out
    assert "o." in run("display r(varlist)")


def test_mark_markin_markout(run):
    run("clear\nset obs 10\ngen x = _n\nreplace x = . in 3\nreplace x = .a in 4")
    run("mark m1 if x > 2")
    assert run("count if m1").strip() == "8"
    run("markin if x > 5, name(rng)")
    assert run('display "`rng\'"').strip() == "in 3/10"
    run("markin in 2/6")
    assert run('display "`in\'"').strip() == "in 2/6"
    run("gen byte mk = 1\nmarkout mk x")
    assert run("count if mk").strip() == "8"
    run("gen byte mk2 = 1\nmarkout mk2 x, sysmissok")
    assert run("count if mk2").strip() == "9"
    assert "already defined" in run("capture noisily mark m1")
    run("gen mkf = 1\nmarkout mkf x")
    assert run.rc == 198


def test_mark_weights(run):
    run("clear\nset obs 3\ngen w = cond(_n==2, 0, 1)")
    run("mark m2 [fweight=w]")
    assert run("count if m2").strip() == "2"
    run("mark m3 [fweight=w], zeroweight")
    assert run("count if m3").strip() == "3"


def test_marksample_with_fv_ts(run):
    run(DATA)
    run("""
program define myp
    syntax [varlist(fv ts)] [if] [in]
    marksample touse
    count if `touse'
end
""")
    assert run("myp L.mpg").strip() == "19"
    assert run("myp mpg").strip() == "20"


def test_postfile(run, tmp_path):
    f = str(tmp_path / "res").replace("\\", "/")
    run(f"""
tempname m
postfile `m' a b str5 c double(d e) using "{f}", replace
post `m' (1) (2) ("abcdefgh") (3.5) (.)
postclose `m'
use "{f}", clear
""")
    out = run("describe")
    assert "str5" in out and "double" in out
    assert run("count").strip() == "1"
    assert run("display c[1]").strip() == "abcde"


def test_post_errors(run, tmp_path):
    f = str(tmp_path / "r2").replace("\\", "/")
    run(f'tempname m\npostfile `m\' a b using "{f}", replace')
    assert "2 expressions expected and only 1 found" in run("capture noisily post `m' (1)")
    # the error closes the postfile
    assert "not found" in run("capture noisily postclose `m'")
    assert "already exists" in run(f'capture noisily postfile x z using "{f}"')


def test_postfile_loop_in_program(run, tmp_path):
    f = str(tmp_path / "sim").replace("\\", "/")
    run(f"""
program lnsim
    tempname sim
    postfile `sim' mean using "{f}", replace
    forvalues i = 1/5 {{
        post `sim' (`i')
    }}
    postclose `sim'
end
lnsim
use "{f}", clear
""")
    assert run("count").strip() == "5"
    assert run("quietly summarize mean\ndisplay r(mean)").strip() == "3"


def test_matlist(run):
    run("matrix A = (1,2\\3,4)\nmatrix rownames A = r1 r2\nmatrix colnames A = c1 c2")
    out = run("matlist A")
    assert "c1" in out and "r2" in out and "4" in out
    out = run("matlist A, format(%5.1f) title(Hi)")
    assert "Hi" in out and "4.0" in out
    assert run("matlist A, noheader") != out


def test_tabdisp(run):
    run("clear\nset obs 6\ngen a = mod(_n,2)\ngen b = 1+mod(_n,3)\ngen v = _n")
    out = run("tabdisp a b, cell(v)")
    assert "b" in out and "|" in out
    assert [ln.split("|")[-1].split() for ln in out.splitlines() if ln.strip().startswith(("0", "1"))] \
        == [["6", "4", "2"], ["3", "1", "5"]]
    out2 = run("tabdisp a, cell(v)")
    assert "v" in out2


def test_sortpreserve(run):
    run("""
clear
set obs 5
gen x = 6 - _n
program define srt
    syntax
    sort x
end
""")
    run("program define srt2, sortpreserve\n sort x\nend\nsrt2")
    assert run("display x[1]").strip() == "5"
    run("program define srt3\n sort x\nend\nsrt3")
    assert run("display x[1]").strip() == "1"


def test_plugin_explains(run, tmp_path):
    out = run("program foo, plugin using(odta_zzz_missing.plugin)")
    assert run.rc == 601 and "not found" in out
    f = tmp_path / "odta_p.plugin"
    f.write_text("x")
    out = run(f'program foo, plugin using("{f.as_posix()}")')
    assert run.rc == 198 and "cannot be loaded" in out


def test_sleep_pause_window_batch(run):
    run("sleep 10")
    assert run.rc == 0
    run("pause off\npause")
    assert run.rc == 0
    run("window push dir")
    assert run.rc == 0


def test_viewsource(run):
    assert run("viewsource unab.ado").strip()
    run("viewsource nonexistent_zzz.ado")
    assert run.rc == 601


def test_signestimationsample_and_makecns(run):
    run("clear\nset obs 30\nset seed 1\ngen x = rnormal()\ngen z = rnormal()\ngen y = x + z + rnormal()")
    run("regress y x z")
    run("signestimationsample x z")
    assert run("display e(datasignature)").strip() != ""
    run("checkestimationsample")
    assert run.rc == 0
    run("replace x = x + 1 in 1")
    run("checkestimationsample")
    assert run.rc == 459
    run("matrix C = (1, -1, 0, 0)")
    run("makecns C")
    assert run.rc == 0


def test_tempnames_are_hex_and_freed_by_programs(run):
    run("tempname a b c d e f g h i j k")
    assert run("display \"`k'\"").strip() == "__00000A"
    run("program define pp\n tempvar t\n display \"`t'\"\nend")
    first = run("pp").strip()
    assert run("pp").strip() == first


def test_rmcoll_omits_the_last(run):
    run("clear\nset obs 20\ngen a = _n\ngen b = mod(_n,7)\ngen c = a + b")
    out = run("_rmcoll a b c")
    assert "c omitted" in out
    assert run("display r(varlist)").strip() == "a b o.c"


def test_list_duplicates_and_noobs_clean(run):
    run("clear\nset obs 2\ngen x = _n")
    out = run("list x x, noobs clean")
    assert out.split()[:2] == ["x", "x"]
    assert out.splitlines()[1].startswith("    ")


def test_makecns_returns(run):
    run("clear\nset obs 30\nset seed 2\ngen x1 = rnormal()\ngen x2 = rnormal()\ngen y = x1 + x2 + rnormal()")
    run("quietly regress y x1 x2\nmatrix C = (1, -1, 0, 0)\nmakecns C, r")
    assert run("display r(cns1)").strip() == "x1 - x2 = 0"
    assert run("display r(k)").strip() == "1"
    out = run("makecns Zzz")
    assert "constraint matrix 'Zzz' caused error r(111)" in out


def test_makecns_with_a_constraint_list(run):
    run("clear\nset obs 30\nset seed 3\ngen x1 = rnormal()\ngen x2 = rnormal()\ngen y = x1 + x2 + rnormal()")
    run("quietly regress y x1 x2")
    run("constraint 1 x1 = x2\nconstraint define 2 _cons = 1\nconstraint 3 2*x1 = 2*x2\nconstraint 4 zz = 0")
    out = run("makecns 1-4, displaycns r")
    assert "(note: constraint number 3 caused error r(412))" in out      # redundant with 1
    assert "(note: constraint number 4 caused error r(111))" in out      # zz is not a coefficient
    assert " ( 1)  x1 - x2 = 0\n ( 2)  _cons = 1\n" in out
    assert run('display "`r(clist)\'" " " r(k)').strip() == "1 2 2"
    assert run("display el(e(Cns), 2, 4)").strip() == "1"
    assert run("makecns 1 4, nocnsnotes") == ""
    assert run("matrix dispCns") == " ( 1)  x1 - x2 = 0\n"
    run("matrix dispCns, r")
    assert run("display r(k) \" \" r(cns1)").strip() == "1 x1 - x2 = 0"


def test_makecns_names_equations(run):
    run("matrix b = (1, 2, 3, 4)\nmatrix colnames b = y:x y:_cons z:x z:_cons")
    run("ereturn post b\nconstraint 1 [y]x = [z]x")
    run("makecns 1, r")
    assert run("display r(cns1)").strip() == "[y]x - [z]x = 0"


def test_matlist_equation_options(run):
    run("matrix A = (1,2,3,4 \\ 2,5,6,7 \\ 3,6,8,9 \\ 4,7,9,10)")
    run("matrix rownames A = y:x1 y:_cons z:x1 z:_cons\nmatrix colnames A = y:x1 y:_cons z:x1 z:_cons")
    out = run("matlist A, lines(eq)")
    lines = out.splitlines()
    assert lines[1] == "             |          y           |          z"
    assert lines[2] == "             |        x1      _cons |        x1      _cons"
    assert lines[4].startswith("y            |")
    assert lines[7].startswith("-------------+")
    out = run("matlist A, showcoleq(each) roweqonly")
    assert out.splitlines()[1].split("|")[1].split() == ["y", "y", "z", "z"]
    assert out.splitlines()[4].startswith("           y |")
    assert "_cons" not in run("matlist A, coleqonly").splitlines()[1]
    run("matlist A, showcoleq(sideways)")
    assert run.rc == 198
