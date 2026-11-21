"""System and dataset state (compat 2003, 2006 and 2007): Stata error codes for failed commands,
c() values, c(changed) and the sort marker."""

import pytest


@pytest.mark.parametrize("cmd", ["ksmirnov", "ladder", "prtest", "runtest", "sdtest"])
def test_bare_command_needs_varlist(run, cmd):
    out = run(f"clear\nset obs 5\ngen x = _n\ncapture noisily {cmd}\ndisplay _rc")
    assert out.splitlines()[-2:] == ["varlist required", "100"]


@pytest.mark.parametrize("cmd", [
    'copy "{src}" "{bad}/b.txt", replace',
    'file open fh using "{bad}/a.txt", write replace',
    'outfile x using "{bad}/a.txt", replace',
    'log using "{bad}/x.log", replace name(probe)',
    'save "{bad}/a.dta", replace',
    'export delimited "{bad}/a.csv", replace',
    'postfile pf x using "{bad}/p.dta", replace',
])
def test_missing_folder_is_r603(run, tmp_path, cmd):
    src = tmp_path / "src.txt"
    src.write_text("a\n")
    bad = (tmp_path / "missing").as_posix()
    run("clear\nset obs 5\ngen x = _n")
    out = run("capture noisily " + cmd.format(src=src.as_posix(), bad=bad) + "\ndisplay _rc")
    lines = out.splitlines()
    assert lines[-1] == "603"
    assert lines[-2].endswith(" could not be opened") and lines[-2].startswith("file ")
    assert lines[-3].startswith("(note: file ") and lines[-3].endswith(" not found)")


def test_python_errors_become_stata_errors(run, monkeypatch):
    from opendta.commands import registry

    def boom(s, args):
        raise KeyError("x")
    import dataclasses
    spec = registry.REGISTRY["describe"]
    monkeypatch.setitem(registry.REGISTRY, "describe", dataclasses.replace(spec, fn=boom))
    assert run("capture describe\ndisplay _rc").strip() == "198"


@pytest.mark.parametrize("name", ["sortlist", "eps", "fastscroll", "lastcmd", "adoplace", "iterlog",
                                  "nao_existe_xyz"])
def test_unknown_c_name_is_r133(run, name):
    out = run(f'capture noisily display "[" c({name}) "]"\ndisplay _rc')
    assert out.splitlines() == [f"[c({name}) undefined", "133"]


def test_c_system_values(run):
    for name, shown in [("maxbyte", "100"), ("minint", "-32767"), ("maxlong", "2.147e+09"),
                        ("epsfloat", "1.192e-07"), ("maxfloat", "1.701e+38"), ("namelen", "32"),
                        ("maxstrvarlen", "2045"), ("matsize", "400"), ("maxvar", "5000"),
                        ("SE", "1"), ("MP", "0"), ("byteorder", "lohi"), ("logtype", "smcl"),
                        ("emptycells", "keep"), ("tracedepth", "32000"), ("noisily", "1"),
                        ("mode", ""), ("showbaselevels", ""), ("filedate", "")]:
        assert run(f'display "[" c({name}) "]"').strip() == f"[{shown}]", name
    assert run('display c(Weekdays)').startswith("Sunday Monday")
    assert run("quietly display c(noisily)") == ""
    assert run("display c(tmpdir)").strip().endswith("/")


def test_s_globals_through_a_local(run):
    out = run('foreach m in S_level S_ADO S_StataSE S_OS {\ndisplay "`m\'=[${`m\'}]"\n}')
    assert out.splitlines()[:3] == ["S_level=[95]", "S_ADO=[BASE;SITE;.;PERSONAL;PLUS;OLDPLACE]",
                                    "S_StataSE=[SE]"]


DATA = """clear
set seed 5
set obs 50
gen id = _n
gen g = mod(_n,3)
gen x = round(invnormal(runiform()),.01)
gen str20 s = "abc"
"""


def test_changed_stays_0_after_sort_compress_and_keep(run, tmp_path):
    f = (tmp_path / "d.dta").as_posix()
    run(DATA + f'save "{f}"')
    for cmd in ("sort g x", "compress", "keep if id > 0"):
        run(cmd)
        assert run("display c(changed)").strip() == "0", cmd
    assert "has changed" not in run("describe, short")
    run("replace x = 1")
    assert "Note: Dataset has changed since last saved." in run("describe, short")


def test_sort_marker_keeps_repeated_names(run):
    run(DATA + "sort g g x")
    assert "Sorted by: g  g  x\n" in run("describe, short")


def test_merge_1_1_is_sorted_by_the_key(run, tmp_path):
    f = (tmp_path / "d.dta").as_posix()
    run(DATA + f'save "{f}"\nsort g x\nmerge 1:1 id using "{f}", nogenerate')
    assert "Sorted by: id\n" in run("describe, short")
    assert run("list id in 1/3, noobs clean").split()[1:] == ["1", "2", "3"]


def test_merge_with_using_only_rows_is_not_sorted(run, tmp_path):
    f = (tmp_path / "u.dta").as_posix()
    run(f'clear\ninput id u\n4 1\n1 2\nend\nsave "{f}"')
    run('clear\ninput id m\n3 1\n1 2\nend\nmerge 1:1 id using "' + f + '", nogenerate')
    assert run("list id, noobs clean").split()[1:] == ["1", "3", "4"]
    assert "Sorted by: \n" in run("describe, short")


def test_label_list_after_drop_merge(run, tmp_path):
    f = (tmp_path / "d.dta").as_posix()
    run(DATA + 'label define lab 0 "zero"\nlabel values g lab\n' + f'save "{f}"')
    run(f'merge 1:1 id using "{f}"')
    assert "_merge:" in run("label list")
    run("drop _merge")
    out = run("label list")
    assert "_merge" not in out and "lab:" in out


_LONG = "zero com texto bem comprido para testar quebra"


def _labeled(run):
    run(f'clear\nset obs 6\ngen g = mod(_n, 3)\ngen x = _n\ngen str3 s = "abc"\n'
        f'label define lg 0 "{_LONG}" 1 "um" 2 "dois"\nlabel values g lg')


def test_tab_long_labels(run):
    # compat 2007: one-way tables widen the column up to 39 characters, two-way to 21
    _labeled(run)
    lines = run("tab g").splitlines()
    assert lines[1] == " " * 38 + "g |      Freq.     Percent        Cum."
    assert lines[3].startswith(_LONG[:39] + " |")
    lines = run("tab g s").splitlines()
    assert lines[2] == " " * 20 + "g |       abc |     Total"
    assert lines[4].startswith(_LONG[:21] + " |")


def test_tabstat_one_variable_by(run):
    # compat 2007: with one variable the statistics are the columns; labels cut to 16
    _labeled(run)
    lines = run("tabstat x, by(g) stat(mean sd n)").splitlines()
    assert lines[1:3] == ["Summary for variables: x", "     by categories of: g"]
    assert lines[4] == "               g |      mean        sd         N"
    assert lines[6].startswith(_LONG[:16] + " |")
    assert lines[-2].startswith("           Total |")


def test_compress_keeps_numeric_format(run):
    run('clear\nset obs 3\ngen id = _n\ngen str20 s = "abc"\nquietly compress')
    d = run.session.data
    assert (d.get("id").vtype, d.get("id").fmt) == ("byte", "%9.0g")
    assert (d.get("s").vtype, d.get("s").fmt) == ("str3", "%9s")


def test_use_subset_keeps_sort_prefix_and_changed(run, tmp_path):
    # compat 2001: data sorted by a b c read as `use a c` is sorted by a, and use leaves
    # c(changed) at 0 even with if
    f = (tmp_path / "t.dta").as_posix()
    run(f'clear\nset obs 5\ngen a = _n\ngen b = _n\ngen c = _n\nsort a b c\nsave "{f}"')
    run(f'use a c using "{f}", clear')
    assert run.session.data.sortlist == ["a"] and run.eval("c(changed)") == 0
    run(f'use a c using "{f}" if a > 2, clear')
    assert run.eval("c(changed)") == 0
    run(f'use "{f}", clear\ndrop b')
    assert run.session.data.sortlist == ["a"]


def test_batch_log_is_the_unnamed_log(run):
    run.session.batch_log = True
    out = run('capture noisily log using "/no/such/x.log", text replace\ndisplay _rc').splitlines()
    assert out[-2:] == ["log file already open", "604"]
