"""hexdump ([D] hexdump): the bytes of a file, in hexadecimal and as characters, or a summary
of the kinds of characters it has."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lang.words import strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def dump_lines(data: bytes, start: int = 0) -> list[str]:
    """16 bytes per line: address, the bytes in groups of two and the printable characters
    (others as a dot). VERIFY: Stata's header and spacing."""
    out = []
    for off in range(0, len(data), 16):
        chunk = data[off:off + 16]
        hexes = " ".join(chunk[i:i + 2].hex() for i in range(0, len(chunk), 2))
        text = "".join(chr(b) if 32 <= b < 127 else "." for b in chunk)
        out.append(f"{start + off:>11}   {hexes:<39}  | {text}")
    return out


def analyze(data: bytes) -> dict[str, float]:
    """Counts of the kinds of bytes: line ends (Windows \\r\\n, Mac \\r, Unix \\n), tabs,
    blanks, letters, digits, other ASCII, control and extended (128-255) bytes."""
    win = data.count(b"\r\n")
    mac = data.count(b"\r") - win
    unix = data.count(b"\n") - win
    r = {"Windows": win, "Mac": mac, "Unix": unix, "tab": data.count(b"\t"),
         "blank": data.count(b" "), "uc": 0, "lc": 0, "digit": 0, "special": 0, "ctl": 0,
         "null": data.count(b"\0"), "extended": 0, "filesize": len(data)}
    for b in data:
        if 65 <= b <= 90:
            r["uc"] += 1
        elif 97 <= b <= 122:
            r["lc"] += 1
        elif 48 <= b <= 57:
            r["digit"] += 1
        elif 33 <= b <= 126:
            r["special"] += 1
        elif b >= 128:
            r["extended"] += 1
        elif b not in (9, 10, 13, 32, 0):
            r["ctl"] += 1
    return {k: float(v) for k, v in r.items()}


@command("hexdump")
def cmd_hexdump(s: "Session", args: str) -> None:
    """hexdump filename [, analyze results from(#) to(#)]. VERIFY: the layout of the dump and of
    analyze, and the names of the r() results."""
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    name = strip_outer_quotes(body.strip())
    if not name:
        raise StataError(198, "invalid file specification")   # compat 2005
    o = match_options(opt, {"analyze": 1, "results": 1, "from": 4, "to": 2, "tabulate": 3,
                            "noextended": 5})
    if o.get("tabulate"):
        raise StataError(198, "option tabulate is not available yet in OpenDTA")
    if not os.path.isfile(name):
        raise StataError(601, f"file {name} not found")
    with open(name, "rb") as fh:
        data = fh.read()
    lo = int(o["from"]) if o.get("from") else 0
    hi = int(o["to"]) + 1 if o.get("to") else len(data)
    part = data[lo:hi]
    out = s.output
    if not (o.get("analyze") or o.get("results")):
        out.write("\n" + " " * 11 + "   " + f"{'hex representation':<39}  | character\n", "text")
        out.write("-" * 56 + "+" + "-" * 17 + "\n", "text")
        for line in dump_lines(part, lo):
            out.write(line + "\n", "result")
        return
    r = analyze(part)
    if o.get("analyze"):
        out.write("\n  Line-end characters\n", "text")
        for label, k in (("\\r\\n (Windows)", "Windows"), ("\\r by itself (Mac)", "Mac"),
                         ("\\n by itself (Unix)", "Unix")):
            out.write(f"    {label:<34}{int(r[k]):>12,}\n", "result")
        out.write("  Space/separator characters\n", "text")
        for label, k in (("[blank]", "blank"), ("[tab]", "tab")):
            out.write(f"    {label:<34}{int(r[k]):>12,}\n", "result")
        out.write("  Other characters\n", "text")
        for label, k in (("A-Z", "uc"), ("a-z", "lc"), ("0-9", "digit"),
                         ("Special (!@#$ etc.)", "special"), ("Control", "ctl"),
                         ("Binary 0", "null"), ("Extended (128-255)", "extended")):
            out.write(f"    {label:<34}{int(r[k]):>12,}\n", "result")
        out.write(f"  {'Total':<36}{int(r['filesize']):>12,}\n", "result")
    s.r = r
