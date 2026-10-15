"""tools/compare.py: limits of the known differences (ignore and tolerance)."""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))

import compare as C  # noqa: E402


def test_tolerance_accepts_small_numeric_differences_only():
    k = C.Known("why", tol=1e-5)
    a = ["       _cons |   10.12385   .6077607    16.66   0.000      8.93266    11.31504"]
    b = ["       _cons |   10.12385   .6077604    16.66   0.000     8.932661    11.31504"]
    assert k.accepts(a, b)
    assert not k.accepts(a, ["       _cons |   10.12385   .6077604    16.66   0.000     8.932661    11.40000"])
    assert not k.accepts(a, ["       _kons |   10.12385   .6077604    16.66   0.000     8.932661    11.31504"])


def test_ignore_drops_lines_from_both_logs():
    k = C.Known("why", ignore=re.compile(r"^Iteration \d+:"))
    a = ["Iteration 0:   ll = -1", "Iteration 1:   ll = -2", "Log likelihood = -2"]
    b = ["Iteration 0:   ll = -1.5", "Log likelihood = -2"]
    assert k.accepts(a, b)
    assert not k.accepts(a, ["Log likelihood = -3"])


def test_parse_limits(tmp_path, monkeypatch):
    f = tmp_path / "known_differences.txt"
    f.write_text("x  ignore(/^It/) tolerance(1e-3) reason here\ny plain reason\n", encoding="utf-8")
    monkeypatch.setattr(C, "COMPAT", tmp_path)
    known = C._known_differences()
    assert known["x"].why == "reason here" and known["x"].tol == 1e-3
    assert known["x"].ignore.pattern == "^It"
    assert not known["y"].limited()
