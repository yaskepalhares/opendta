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
