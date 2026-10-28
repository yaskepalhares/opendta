"""tools/compare.py: limits of the known differences (ignore, skip and tolerance)."""

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


def test_skip_drops_the_output_of_matching_commands_only():
    k = C.Known("why", skip=re.compile(r"xtunitroot llc"))
    a = [". xtunitroot llc rw", "Levin-Lin-Chu unit-root test", "-----", ". xtunitroot ht rw", "rho 0.5"]
    b = [". xtunitroot llc rw", "not available yet", ". xtunitroot ht rw", "rho 0.5"]
    assert k.accepts(a, b)
    assert not k.accepts(a, [". xtunitroot llc rw", "not available yet", ". xtunitroot ht rw", "rho 0.6"])


def test_parse_limits(tmp_path, monkeypatch):
    f = tmp_path / "known_differences.txt"
    f.write_text("x  ignore(/^It/) tolerance(1e-3) skip(/llc/) reason here\ny plain reason\n", encoding="utf-8")
    monkeypatch.setattr(C, "COMPAT", tmp_path)
    known = C._known_differences()
    assert known["x"].why == "reason here" and known["x"].tol == 1e-3
    assert known["x"].ignore.pattern == "^It" and known["x"].skip.pattern == "llc"
    assert not known["y"].limited()


def test_brief_diff_leaves_out_iterations_and_is_cut(tmp_path, capsys, monkeypatch):
    """--brief --show-diff: only the cases that differ, without iteration lines, at most
    --max-diff lines of diff (the CI annotation, printed only under GitHub Actions, is left out)."""
    monkeypatch.delenv("GITHUB_ACTIONS", raising=False)
    exp = tmp_path / "expected"
    exp.mkdir()
    out = C.run_opendta(C.COMPAT / "do" / "1910_programming.do", C.COMPAT / "out", "")
    lines = out.read_text(encoding="utf-8").splitlines()
    k = next(i for i, ln in enumerate(lines) if ln.startswith(". constraint"))
    lines[k + 1:k + 1] = ["Iteration 1:   log likelihood = -1"] + [f"extra {j}" for j in range(10)]
    (exp / "1910_programming.log").write_text("\n".join(lines) + "\n", encoding="utf-8")
    C.main(["1910", "--brief", "--show-diff", "--max-diff", "6", "--expected", str(exp)])
    text = capsys.readouterr().out
    assert "XX 1910_programming" in text and "Iteration 1:" not in text
    assert "more diff lines" in text
