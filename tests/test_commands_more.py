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
