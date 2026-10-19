"""Compares OpenDTA's output with reference logs produced in Stata 14.

Usage:
    python tools/compare.py              # runs all compat/do/*.do
    python tools/compare.py 0001         # only those starting with 0001
    python tools/compare.py --show-diff  # shows the differences line by line

For each compat/do/NAME.do:
  1. runs it in OpenDTA in batch mode, writing compat/out/NAME.log;
  2. if compat/expected/NAME.log exists (produced in Stata), compares the two.

Before comparing, the log header and footer, `. do ...`/`. log close` lines
and trailing spaces are dropped; blank lines are ignored.

Numbers written with 15 or more characters on either side (return list,
ereturn list) may differ in the last digit: the order of floating-point sums is not
the same as Stata's. These lines count as equal when the rest of the text is
identical and the relative difference is below 1e-11 (Wald statistics and
computations with numerical Hessians accumulate error in the last two digits).
"""

from __future__ import annotations

import argparse
import difflib
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

COMPAT = ROOT / "compat"
_HEADER = re.compile(r"^\s*(name|log|log type|opened on|closed on):")
_SKIP = re.compile(r"^\.\s+((capture\s+)?(noisily\s+)?(do|run)|log close|log using|set linesize|quietly set (dp|hints|linesize|superscript|numerics))\b")


# .dta date stamp (describe): changes on every run
_STAMP = re.compile(r"\b\d{1,2} [A-Z][a-z]{2} \d{4} \d{2}:\d{2}\b")
# absolute path of the case's temporary files (graph export to PDF shows the full path, which
# changes from machine to machine): absolute path of odta_* files (macOS/Linux "/..."; Windows
# "C:\\...")
# e(datasignature): Stata's checksum algorithm is undocumented; N and k are compared, the two
# checksums are not
_DATASIG = re.compile(r'(e\(datasignature\) : "\d+:\d+):\d+:\d+"')
_ABSPATH = re.compile(r"(?<=\(file )(?:[A-Za-z]:)?[\\/]\S*[\\/](odta_\S+)")


def normalize(text: str) -> list[str]:
    out: list[str] = []
    lines = text.splitlines()
    for k, raw in enumerate(lines):
        line = _ABSPATH.sub(r"<folder>/\1", _STAMP.sub("<date>", raw.rstrip()))
        line = _DATASIG.sub(r'\1:<checksum>"', line)
        if not line or _HEADER.match(line) or _SKIP.match(line):
            continue
        if set(line) <= {"-"}:
            # only the dashed lines of the log header and footer are dropped; those of the tables
            # are part of the output
            near = lines[max(0, k - 1):k] + lines[k + 1:k + 2]
            if any(_HEADER.match(n) for n in near) or k == 0 or k == len(lines) - 1:
                continue
        if line in (".",):
            continue
        out.append(line)
    # batch mode ends an interrupted do-file with "end of do-file" and r(#); in Stata the
    # references run under capture noisily, without the final r(#)
    if len(out) >= 2 and out[-2] == "end of do-file" and re.fullmatch(r"r\(\d+\);", out[-1]):
        out.pop()
    return out


_NUM = re.compile(r"-?(?:\d+\.?\d*|\.\d+)(?:e[-+]\d+)?")


def _digits(tok: str) -> int:
    mant = tok.lstrip("-").split("e")[0].replace(".", "").lstrip("0")
    return len(mant)


def _close_line(a: str, b: str) -> bool:
    """Same line, up to noise in the 16th digit of long numbers."""
    if _NUM.sub("#", a) != _NUM.sub("#", b):
        return False
    na, nb = _NUM.findall(a), _NUM.findall(b)
    for x, y in zip(na, nb):
        if x == y:
            continue
        # long numbers (%18.0g): 15+ characters, such as .0032959009263453, on either side
        # (112.8 against 112.8000000000001)
        if max(len(x.lstrip("-")), len(y.lstrip("-"))) < 15:
            return False
        fx, fy = float(x), float(y)
        if abs(fx - fy) > 1e-11 * max(abs(fx), abs(fy)):
            return False
    return True


def same_logs(a: list[str], b: list[str]) -> bool:
    return len(a) == len(b) and all(x == y or _close_line(x, y) for x, y in zip(a, b))


def run_opendta(dofile: Path, outdir: Path, setup: str = "") -> Path:
    from opendta.cli import run_batch

    outdir.mkdir(parents=True, exist_ok=True)
    # OpenDTA's error explanations and superscript exponent turned off: Stata does not have them.
    # `setup` (e.g. set dp comma) runs first, without appearing in the compared log
    target = outdir / dofile.name
    # gerar_esperados.do runs the cases with set linesize 255
    pre = ("quietly set hints off\nquietly set superscript off\nquietly set numerics stata\n"
           "quietly set linesize 255\n"
           + (f"quietly {setup}\n" if setup else ""))
    target.write_text(pre + dofile.read_text(encoding="utf-8"), encoding="utf-8")
    # leftovers of an interrupted case (odta_*) must not affect the next one; gerar_esperados.do
    # does the same cleanup in Stata
    for leftover in outdir.glob("odta_*"):
        if leftover.is_file():
            leftover.unlink()
    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        run_batch(str(target), [])
    finally:
        os.chdir(cwd)
    return outdir / (dofile.stem + ".log")


_LIMIT = re.compile(r"^(ignore\(/(.*?)/\)|tolerance\(([0-9.eE+-]+)\))\s*")


@dataclass
class Known:
    """A known difference: the reason and, optionally, its limits. ignore(/REGEX/) drops the
    lines matching REGEX from both logs; tolerance(REL) accepts numbers that differ by at most
    REL relative to the larger of |a|, |b| and 1. With limits, the rest of the log must still be
    equal, so layout regressions are still caught."""
    why: str
    ignore: re.Pattern | None = None
    tol: float | None = None

    def limited(self) -> bool:
        return self.ignore is not None or self.tol is not None

    def accepts(self, a: list[str], b: list[str]) -> bool:
        if not self.limited():
            return True
        if self.ignore is not None:
            a = [x for x in a if not self.ignore.search(x)]
            b = [x for x in b if not self.ignore.search(x)]
        if len(a) != len(b):
            return False
        return all(x == y or _close_line(x, y) or (self.tol is not None and _within(x, y, self.tol))
                   for x, y in zip(a, b))


def _within(a: str, b: str, tol: float) -> bool:
    # numbers with fewer digits shift the columns: spacing is not compared here
    if " ".join(_NUM.sub("#", a).split()) != " ".join(_NUM.sub("#", b).split()):
        return False
    for x, y in zip(_NUM.findall(a), _NUM.findall(b)):
        fx, fy = float(x), float(y)
        if abs(fx - fy) > tol * max(abs(fx), abs(fy), 1.0):
            return False
    return True


def _known_differences() -> dict[str, Known]:
    """compat/known_differences.txt: cases that may differ, with the reason. The reason may
    start with ignore(/REGEX/) and/or tolerance(REL) (see Known)."""
    path = COMPAT / "known_differences.txt"
    out: dict[str, Known] = {}
    if path.exists():
        for line in path.read_text(encoding="utf-8").splitlines():
            if line.strip() and not line.lstrip().startswith("#"):
                name, _, why = line.strip().partition(" ")
                k = Known(why.strip())
                while True:
                    m = _LIMIT.match(k.why)
                    if not m:
                        break
                    if m.group(2) is not None:
                        k.ignore = re.compile(m.group(2))
                    else:
                        k.tol = float(m.group(3))
                    k.why = k.why[m.end():]
                out[name] = k
    return out


def main(argv: list[str] | None = None) -> int:
    # the Windows console (cp1252) cannot print "≈", "·" etc. from the descriptions
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass
    p = argparse.ArgumentParser()
    p.add_argument("prefix", nargs="?", default="")
    p.add_argument("--show-diff", action="store_true")
    p.add_argument("--expected", default=str(COMPAT / "expected"),
                   help="folder with the reference logs (default: compat/expected)")
    p.add_argument("--setup", default="",
                   help="command run before each do-file, e.g. 'set dp comma'")
    ns = p.parse_args(argv)

    files = sorted(f for f in (COMPAT / "do").glob("*.do") if f.name.startswith(ns.prefix))
    known = _known_differences()
    passed = failed = missing = tolerated = 0
    for f in files:
        got = run_opendta(f, COMPAT / "out", ns.setup)
        exp = Path(ns.expected) / (f.stem + ".log")
        if not exp.exists():
            print(f"  ?  {f.stem:<32} no Stata reference log")
            missing += 1
            continue
        a = normalize(exp.read_text(encoding="utf-8", errors="replace"))
        b = normalize(got.read_text(encoding="utf-8", errors="replace"))
        if same_logs(a, b):
            print(f"  ok {f.stem}")
            passed += 1
        else:
            k = known.get(f.stem)
            if k is not None and k.accepts(a, b):
                print(f"  ~~ {f.stem:<32} known difference: {k.why}")
                tolerated += 1
            else:
                print(f"  XX {f.stem}")
                failed += 1
                if os.environ.get("GITHUB_ACTIONS"):
                    # GitHub annotation with the first differing lines, to read the difference
                    # without opening the whole CI log
                    diff = [ln for ln in difflib.unified_diff(a, b, "Stata 14", "OpenDTA", lineterm="", n=0)
                            if ln[:1] in "+-" and not ln.startswith(("---", "+++"))][:12]
                    body = "%0A".join(ln.replace("%", "%25") for ln in diff)
                    print(f"::error title=compat {f.stem}::{body}")
            if ns.show_diff:
                for line in difflib.unified_diff(a, b, "Stata 14", "OpenDTA", lineterm="", n=1):
                    print("     " + line)
    extra = f", {tolerated} with a known difference" if tolerated else ""
    print(f"\n{passed} equal, {failed} different{extra}, {missing} without a reference")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
