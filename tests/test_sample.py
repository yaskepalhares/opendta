"""sample: simple random sampling without replacement."""

import numpy as np
import pytest

from opendta.session import Session

BASE = "clear\nset seed 310\nset obs 40\ngen id = _n\ngen g = mod(_n, 3)\n"


def _run(extra: str):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(BASE + extra)
    return s, "".join(lines), rc


def test_percentage_matches_stata_draws_and_order():
    # Stata 14 (compat 0310): the same ten observations, left in the order of the draws
    s, out, rc = _run("sample 25")
    assert rc == 0 and "(30 observations deleted)" in out
    assert s.data.get("id").data.tolist() == [16, 40, 7, 30, 22, 6, 17, 34, 36, 12]
    assert s.data.sortlist == []
    s, out, _ = _run("keep in 1/37\nsample 10")      # 3.7 -> 4
    assert s.data.nobs == 4


def test_count_and_same_seed_same_sample():
    s1, _, _ = _run("sample 7, count")
    s2, _, _ = _run("sample 7, count")
    a, b = s1.data.get("id").data, s2.data.get("id").data
    assert len(a) == 7 and np.array_equal(a, b)
    s3, _, _ = _run("sample 100, count")             # more than N: all kept
    assert s3.data.nobs == 40


def test_if_keeps_the_other_observations():
    s, out, _ = _run("sample 50 if g == 0")
    g = s.data.get("g").data
    assert np.sum(g != 0) == 27 and np.sum(g == 0) == 7     # 13 * .5 = 6.5 -> 7


def test_by_option_and_by_prefix():
    s, _, _ = _run("sample 2, count by(g)")
    assert sorted(np.bincount(s.data.get("g").data.astype(int)).tolist()) == [2, 2, 2]
    s, _, _ = _run("bysort g: sample 50")
    assert np.bincount(s.data.get("g").data.astype(int)).tolist() == [7, 7, 7]


@pytest.mark.parametrize("cmd,msg,rc", [
    ("sample 120", "120:  # must be between 0 and 100", "198"),
    ("sample 2.5, count", "'2.5' found where integer expected", "126"),
    ("sample", "invalid syntax", "198")])
def test_errors(cmd, msg, rc):
    _s, out, _ = _run(f"capture noisily {cmd}\ndisplay _rc")
    assert msg in out and f"\n{rc}\n" in out


def test_negative_count_explains_and_drops_all():
    # Stata 14 carries on after the message (compat 0310)
    s, out, _ = _run("sample -1, count\ndisplay _rc")
    assert "-1:  # must be 0 or positive" in out and "(40 observations deleted)" in out
    assert s.data.nobs == 0


def test_sample_100_draws_nothing():
    s, out, _ = _run("sample 100\ndisplay runiform()")
    assert "deleted" not in out
    t, out2, _ = _run("display runiform()")
    assert out.strip().splitlines()[-1] == out2.strip().splitlines()[-1]
