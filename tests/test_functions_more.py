"""Functions of manual [FN] that research do-files use (compat 2005 and 2008): dates and times,
strings, mathematics, distributions, constants and files."""

import math

import pytest

from opendta.core import missing as M


@pytest.mark.parametrize("expr, value", [
    ("ym(2020,1)", 720), ("yq(2020,1)", 240), ("yw(2020,1)", 3120), ("yh(2020,1)", 120),
    ("mofd(21915)", 720), ("dofm(720)", 21915), ("qofd(21915)", 240), ("yofd(21915)", 2020),
    ("wofd(21915)", 3120), ("week(21915)", 1), ("halfyear(21915)", 1),
    ("dofq(241)", 22006), ("dofh(121)", 22097), ("dofw(3121)", 21922), ("dofy(2020)", 21915),
    ("hofd(22097)", 121), ("week(mdy(12,31,2020))", 52),
    ("daily(\"01jan2020\",\"DMY\")", 21915), ("monthly(\"2020m1\",\"YM\")", 720),
    ("quarterly(\"2020q3\",\"YQ\")", 242), ("halfyearly(\"2020h2\",\"YH\")", 121),
    ("weekly(\"2020w2\",\"YW\")", 3121), ("yearly(\"2020\",\"Y\")", 2020),
    ("date(\"20200105\",\"YMD\")", 21919), ("date(\"1/2/99\",\"MDY\",2020)", 14246),
    ("date(\"1/2/99\",\"MD19Y\")", 14246),
    ("hours(3600000)", 1), ("dhms(1,1,1,1)", 90061000), ("msofhours(1)", 3600000),
    ("cofd(21915)", 1893456000000), ("dofc(1893456000000)", 21915),
    ("mdyhms(1,1,2020,1,1,1)", 1893456000000 + 3661000),
    ("clock(\"01jan2020 10:00\",\"DMYhm\")", 1893456000000 + 36000000),
    ("clock(\"5mar2020 3:04:05.5 pm\",\"DMYhms\")", (21979 * 86400 + 15 * 3600 + 245.5) * 1000),
    ("hh(clock(\"5mar2020 3:04:05.5 pm\",\"DMYhms\"))", 15),
    ("ss(clock(\"5mar2020 3:04:05.5 pm\",\"DMYhms\"))", 5.5),
    # 27 leap seconds from 1972 to 2016
    ("Cofc(1893456000000)", 1893456026000), ("cofC(1893456027000)", 1893456001000),
    ("dofC(Cofd(21915))", 21915),
])
def test_date_functions(run, expr, value):
    assert run.eval(expr) == pytest.approx(value)


@pytest.mark.parametrize("expr", [
    "ym(2020,13)", "daily(\"31feb2020\",\"DMY\")", "daily(\"1/2/99\",\"MDY\")",
    "monthly(\"2020m13\",\"YM\")", "clock(\"01jan2020 25:00\",\"DMYhm\")", "mofd(.)",
])
def test_date_functions_missing(run, expr):
    assert M.is_missing(run.eval(expr))


def test_date_functions_in_generate(run):
    out = run('clear\nset obs 3\ngen m = ym(2020, _n)\ngen d = dofm(m)\nformat d %td\nlist, noobs')
    assert "01mar2020" in out


@pytest.mark.parametrize("expr, value", [
    ('strtoname("a b")', "a_b"), ('strtoname("1x")', "_1x"), ('strtoname("1x", 0)', "1x"), ('strtoname("")', ""),
    ('strtoname("' + "a" * 40 + '")', "a" * 32),
    ('plural(2,"cat")', "cats"), ('plural(1,"cat")', "cat"), ('plural(2,"box","+es")', "boxes"),
    ('plural(2,"glasses","-es")', "glass"), ('plural(2,"man","men")', "men"),
    ('soundex("robert")', "R163"), ('soundex("Rupert")', "R163"), ('soundex("Tymczak")', "T522"),
    ('soundex_nara("Ashcraft")', "A261"),
    ('strrpos("abcb","b")', 4), ('strrpos("abc","z")', 0), ('indexnot("a","b")', 1),
    ('indexnot("aab","a")', 3), ('indexnot("aa","a")', 0),
    ('ustrpos("abc","b")', 2), ('ustrpos("abcb","b",3)', 4), ('ustrrpos("abcb","b")', 4),
    ('usubstr("ação",2,2)', "çã"), ('ustrleft("abc",1)', "a"), ('ustrright("abc",2)', "bc"),
    ('ustrtrim(" a ")', "a"), ('ustrltrim(" a ")', "a "), ('ustrreverse("ação")', "oãça"),
    ('ustrword("a b c",-1)', "c"), ('ustrword("a b c",2)', "b"), ('ustrwordcount("a b c")', 3),
    ('ustrregexm("abc","b")', 1), ('ustrregexm("ABC","b",1)', 1), ('ustrregexm("ABC","b")', 0),
    ('ustrregexra("abab","b","x")', "axax"), ('ustrregexrf("abab","b","x")', "axab"),
    ('uchar(231)', "ç"),
])
def test_string_functions(run, expr, value):
    assert run.eval(expr) == value


def test_ustrregexs(run):
    run('display ustrregexm("ab12", "([a-z]+)([0-9]+)")')
    assert run.eval("ustrregexs(2)") == "12"


@pytest.mark.parametrize("expr, value", [
    ("binomialp(10,5,.5)", .24609375), ("poissonp(2,2)", .27067057), ("reldif(1,1.1)", .04761905),
    ("trigamma(1)", math.pi ** 2 / 6), ("logistic(0)", .5), ("tden(5,0)", .37960669),
    ("invibeta(1,1,.5)", .5), ("nchi2(1,1,1)", .47724987), ("invnchi2(1,1,.5)", 1.1036433),
    ("hypergeometricp(10,5,3,2)", .41666667), ("nbinomial(5,2,.5)", .2265625),
    ("binomial(10,4,invbinomial(10,4,.3))", .3), ("binomialtail(10,4,invbinomialtail(10,4,.3))", .3),
    ("poisson(invpoisson(3,.4),3)", .4), ("poissontail(invpoissontail(3,.4),3)", .4),
    ("nbinomial(5,2,invnbinomial(5,2,.3))", .3), ("nbinomialtail(5,2,invnbinomialtail(5,2,.3))", .3),
    ("nchi2(3,npnchi2(3,4,.2),4)", .2), ("gammap(2,invgammap(2,.3))", .3),
    ("ibetatail(2,3,.4) + ibeta(2,3,.4)", 1), ("cloglog(invcloglog(.7))", .7),
    ("binormal(0,0,0)", .25), ("hypergeometric(10,5,3,3)", 1),
    ("nttail(5,0,1)", .18160873), ("nF(2,3,0,1) - F(2,3,1) + 1", 1), ("lnnormalden(0)", -0.91893853),
])
def test_math_and_distribution_functions(run, expr, value):
    assert run.eval(expr) == pytest.approx(value, rel=1e-6)


@pytest.mark.parametrize("expr, value", [
    ("maxbyte()", 100), ("maxint()", 32740), ("maxlong()", 2147483620), ("minbyte()", -127),
    ("minint()", -32767), ("minlong()", -2147483647), ("float(maxfloat()) == maxfloat()", 1),
    ("recode(5,1,3,5)", 5), ("recode(2,1,3,5)", 3), ("recode(9,1,3,5)", 5),
    ("irecode(5,1,3,5)", 2), ("irecode(0,1,3,5)", 0), ("irecode(9,1,3,5)", 3),
    ("autocode(5,4,0,10)", 5), ("autocode(5.1,4,0,10)", 7.5), ("autocode(11,4,0,10)", 10),
    ('fileexists("nada.xyz")', 0),
])
def test_constants_recode_files(run, expr, value):
    assert run.eval(expr) == pytest.approx(value)


def test_file_functions(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert run.eval('filewrite("a.txt", "hello")') == 5
    assert run.eval('fileexists("a.txt") + filereadable("a.txt")') == 2
    assert run.eval('filewrite("a.txt", "x")') < 0
    assert run.eval('filewrite("a.txt", " world", 2)') == 6
    assert run.eval('fileread("a.txt")') == "hello world"
    assert run.eval('fileread("none.txt")').startswith("fileread() error")


def test_storage_type_functions_are_unknown(run):
    # double(), byte(), long() and strdup() are not functions in Stata 14 (compat 2005)
    out = run("capture noisily display double(5)\ndisplay _rc")
    assert out.splitlines()[-2:] == ["unknown function double()", "133"]
