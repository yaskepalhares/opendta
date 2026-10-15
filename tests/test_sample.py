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


def test_percentage_rounds_to_nearest_and_keeps_order():
    s, out, rc = _run("sample 25")
    assert rc == 0 and "(30 observations deleted)" in out
    ids = s.data.get("id").data
    assert len(ids) == 10 and np.all(np.diff(ids) > 0) and len(set(ids)) == 10
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


@pytest.mark.parametrize("cmd,msg", [("sample 120", "percentage must be between 0 and 100"),
                                     ("sample -1, count", "invalid observation count"),
                                     ("sample", "invalid syntax")])
def test_errors(cmd, msg):
    _s, out, rc = _run(f"capture noisily {cmd}\ndisplay _rc")
    assert msg in out and "198" in out
