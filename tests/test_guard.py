"""Confirmation before risky commands (guard.py)."""

import pytest

from opendta import guard
from opendta.session import Session


def _session(n: int = 20_000) -> Session:
    s = Session()
    s.run(f"quietly set obs {n}")
    s.run("quietly generate year = 2000 + mod(_n, 20)")
    s.run("quietly generate x = runiform()")
    return s


def test_keep_and_drop_need_many_observations_and_a_large_share():
    s = _session()
    risk = guard.check(s, "keep if year >= 2016")           # removes 16,000 of 20,000
    assert risk is not None and risk.kind == "drop_obs" and "16,000" in risk.message
    assert guard.check(s, "quietly keep if year >= 2016") is not None
    assert guard.check(s, "keep if year >= 2000") is None    # nothing removed
    assert guard.check(s, "drop if year == 2000") is None    # 1,000 observations: under the limit
    assert guard.check(_session(100), "keep if year >= 2016") is None   # small data set
    s.run("local n = 2016")
    assert guard.check(s, "keep if year >= `n'") is not None  # macros are expanded first
    assert guard.check(s, "keep if nosuchvar > 1") is None   # the command reports its own error
    assert guard.check(s, "summarize x") is None
    assert risk.command == "keep" and guard.check(s, "drop if year < 2016").command == "drop"   # 16,000 of 20,000
    assert guard.check(s, "keep if year >= 2016", skip=frozenset({"keep"})) is None    # switched off
    assert guard.check(s, "drop if year < 2016", skip=frozenset({"keep"})) is not None


def test_clear_asks_only_with_unsaved_changes():
    s = _session()
    assert guard.check(s, "clear").kind == "clear"
    assert guard.check(s, "clear all").kind == "clear"
    s.data.changed = False
    assert guard.check(s, "clear") is None
    assert guard.check(Session(), "clear") is None


def test_dropping_most_of_the_variables():
    s = Session()
    s.run("quietly set obs 5")
    for k in range(1, 13):
        s.run(f"quietly generate v{k} = {k}")
    assert guard.check(s, "drop v1-v12").kind == "drop_vars"
    assert guard.check(s, "keep v1").kind == "drop_vars"
    assert guard.check(s, "drop v1") is None


def test_save_replace_over_a_much_larger_file(tmp_path):
    s = _session(400_000)
    s.run("quietly generate double y = x * 2")
    path = tmp_path / "big.dta"
    s.run(f'quietly save "{path}"')
    assert path.stat().st_size > guard.MIN_FILE_BYTES
    assert guard.check(s, f'save "{path}", replace') is None          # same data
    assert guard.check(s, f'save "{path}"') is None                   # no replace: Stata itself refuses
    s.run("quietly keep in 1/1000")
    risk = guard.check(s, f'save "{path}", replace')
    assert risk is not None and risk.kind == "save" and "big.dta" in risk.message
    assert guard.check(s, f'save "{tmp_path / "new.dta"}", replace') is None   # no file to lose


@pytest.mark.parametrize("line", ["", "   ", "1 + 1", "by x: keep if x > 0"])
def test_ignores_what_is_not_risky(line):
    assert guard.check(_session(), line) is None
