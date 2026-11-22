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
    (others as a dot), as in compat 2008. VERIFY: whether the address is in hexadecimal."""
    out = []
    for off in range(0, len(data), 16):
        chunk = data[off:off + 16]
        hexes = " ".join(chunk[i:i + 2].hex() for i in range(0, len(chunk), 2))
        text = "".join(chr(b) if 32 <= b < 127 else "." for b in chunk)
        out.append(f"{start + off:>16x} | {hexes:<39} | {text}")
    return out


def _lines(data: bytes) -> list[int]:
    """Lengths of the lines (a tab counts 1) ended by \r\n, \r or \n; a last line with no
    line end counts too. VERIFY: that last line."""
    parts = data.replace(b"\r\n", b"\n").replace(b"\r", b"\n").split(b"\n")
    if parts and parts[-1] == b"":
        parts = parts[:-1]
    return [len(p) for p in parts]


def analyze(data: bytes) -> dict:
    """Counts of the kinds of bytes ([D] hexdump, analyze): line ends (Windows \r\n, Mac \r,
    Unix \n), separators (blank, tab, comma), control bytes (binary 0, other control
    characters but \r \n \t, DEL, 128-159 and 255), printable ASCII (A-Z, a-z, 0-9, other)
    and extended printable bytes (160-254); and the lengths of the lines."""
    win = data.count(b"\r\n")
    r = {"Windows": win, "Mac": data.count(b"\r") - win, "Unix": data.count(b"\n") - win,
         "blank": data.count(b" "), "tab": data.count(b"\t"), "comma": data.count(b","),
         "null": data.count(b"\0"), "ctl": 0, "del": data.count(b"\x7f"), "ext_ctl": 0,
         "uc": 0, "lc": 0, "digit": 0, "special": 0, "extended": 0, "filesize": len(data)}
    for b in data:
        if 65 <= b <= 90:
            r["uc"] += 1
        elif 97 <= b <= 122:
            r["lc"] += 1
        elif 48 <= b <= 57:
            r["digit"] += 1
        elif 33 <= b <= 126 and b != 44:
            r["special"] += 1
        elif 160 <= b <= 254:
            r["extended"] += 1
        elif 128 <= b <= 159 or b == 255:
            r["ext_ctl"] += 1
        elif 1 <= b <= 31 and b not in (9, 10, 13):
            r["ctl"] += 1
    lens = _lines(data)
    r.update({"lmin": min(lens) if lens else None, "lmax": max(lens) if lens else None,
              "lnum": len(lens), "eoleof": int(data.endswith((b"\n", b"\r"))),
              "first": lens[:5]})
    # VERIFY: the names of the formats other than ASCII
    if r["null"] or r["ctl"] or r["del"] or r["ext_ctl"]:
        r["format"] = "BINARY"
    elif r["extended"]:
        r["format"] = "EXTENDED ASCII"
    else:
        r["format"] = "ASCII"
    r["observed"] = sorted(set(data))
    return r


_NAMES = {9: "\\t", 10: "\\n", 13: "\\r", 32: "blank"}


def _report(r: dict) -> list[str]:
    """The two columns of hexdump, analyze (compat 2008): 41 characters on the left, 38 on the
    right."""
    def n(x):
        return "." if x is None else f"{x:,}"
    first = r["first"] + [None] * (5 - len(r["first"]))
    left = [("  Line-end characters", None), ("    \\r\\n         (Windows)", r["Windows"]),
            ("    \\r by itself (Mac)", r["Mac"]), ("    \\n by itself (Unix)", r["Unix"]),
            ("  Space/separator characters", None), ("    [blank]", r["blank"]),
            ("    [tab]", r["tab"]), ("    [comma] (,)", r["comma"]),
            ("  Control characters", None), ("    binary 0", r["null"]),
            ("    CTL excl. \\r, \\n, \\t", r["ctl"]), ("    DEL", r["del"]),
            ("    Extended (128-159,255)", r["ext_ctl"]), ("  ASCII printable", None),
            ("    A-Z", r["uc"]), ("    a-z", r["lc"]), ("    0-9", r["digit"]),
            ("    Special (!@#$ etc.)", r["special"]), ("    Extended (160-254)", r["extended"])]
    right = {0: ("    Line length (tab=1)", ""), 1: ("      minimum", n(r["lmin"])),
             2: ("      maximum", n(r["lmax"])), 4: ("    Number of lines", n(r["lnum"])),
             5: ("      EOL at EOF?", "yes" if r["eoleof"] else "no"),
             7: ("    Length of first 5 lines", ""), 15: ("    File format", r["format"])}
    for k in range(5):
        right[8 + k] = (f"      Line {k + 1}", n(first[k]))
    out = []
    for i, (label, value) in enumerate(left):
        line = label if value is None else f"{label:<26}{n(value):>15}"
        if i in right:
            rl, rv = right[i]
            line = f"{line:<41}" + (f"{rl:<23}{rv:>15}" if rv else rl)
        out.append(line)
    out += [" " * 26 + "-" * 15, f"{'  Total':<26}{n(r['filesize']):>15}", "",
            "  Observed were:",
            "     " + " ".join(_NAMES.get(b, chr(b) if 33 <= b < 127 else f"\\x{b:02x}")
                           for b in r["observed"])]   # VERIFY: bytes with no name; wrapping
    return out


@command("hexdump")
def cmd_hexdump(s: "Session", args: str) -> None:
    """hexdump filename [, analyze results from(#) to(#)]. Layout and r() of compat 2008."""
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
        out.write(" " * 17 + "|" + " " * 41 + "|    character\n"
                  + " " * 17 + "|           hex representation            |  representation\n"
                  + "         address |  0 1  2 3  4 5  6 7  8 9  a b  c d  e f | 0123456789abcdef\n"
                  + "-" * 17 + "+" + "-" * 41 + "+" + "-" * 17 + "\n", "text")
        for line in dump_lines(part, lo):
            out.write(line + "\n", "result")
        return
    r = analyze(part)
    out.write("\n" + "\n".join(_report(r)) + "\n", "result")
    # return list shows them in this order (compat 2008); missing line lengths are left out
    pairs = [("Windows", r["Windows"]), ("lmin", r["lmin"]), ("Mac", r["Mac"]),
             ("lmax", r["lmax"]), ("Unix", r["Unix"]), ("lnum", r["lnum"]),
             ("blank", r["blank"]), ("eoleof", r["eoleof"]), ("tab", r["tab"]),
             ("comma", r["comma"])]
    pairs += [(f"l{k}", v) for k, v in enumerate(r["first"], 1)]
    pairs += [(k, r[k]) for k in ("ctl", "uc", "lc", "digit", "special", "extended", "filesize")]
    res = {k: float(v) for k, v in reversed(pairs) if v is not None}
    res["format"] = r["format"]
    s.r = res
