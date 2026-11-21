"""Differences from Stata 14 seen when running a survey do-file: memory, S_DATE/S_TIME, c(flavor),
log using in a missing folder, use varlist and describe."""

import re

from opendta.io.dta import _valid_sortlist


def test_memory_layout(run):
    out = run("clear\nset obs 4\ngen x = _n\nmemory")
    assert run.rc == 0
    assert "  Memory usage" in out and "grand total" in out and "set maxvar usage" in out
    assert re.search(r"data  +16 +16\n", out)
    run("memory x")
    assert run.rc == 198


def test_system_globals_and_flavor(run):
    out = run('display "[$S_DATE] [$S_TIME] [`c(flavor)\']"')
    assert re.search(r"\[ ?\d{1,2} [A-Z][a-z]{2} \d{4}\] \[\d\d:\d\d:\d\d\] \[IC\]", out)


def test_log_using_missing_folder(run, tmp_path):
    run(f'log using "{tmp_path}/nao/existe/x.log", text')
    assert run.rc == 603
    run(f'capture log using "{tmp_path}/nao/existe/x.log", text\ndisplay _rc')
    assert run.rc == 0


def test_use_varlist_does_not_mark_changed(run, tmp_path):
    f = tmp_path / "d.dta"
    run(f'clear\nset obs 3\ngen a = _n\ngen b = _n\nsave "{f}"\nuse a using "{f}", clear')
    assert run.eval("c(changed)") == 0
    run(f'use a using "{f}" if a > 1, clear')
    assert run.eval("c(changed)") == 1


def test_sortlist_with_repeated_variable_is_unsorted():
    assert _valid_sortlist(["a", "b", "b"]) == []
    assert _valid_sortlist(["a", "b"]) == ["a", "b"]


def test_describe_long_value_label_name(run):
    run("clear\nset obs 2\ngen x = _n\nlabel define lbl_muito_longo 1 a\n"
        "label values x lbl_muito_longo\nlabel variable x rotulo")
    out = run("describe")
    lines = out.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith("x "))
    assert lines[i].rstrip().endswith("lbl_muito_longo")
    assert lines[i + 1] == " " * 46 + "rotulo"
