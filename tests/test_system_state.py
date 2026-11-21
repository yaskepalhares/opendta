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
