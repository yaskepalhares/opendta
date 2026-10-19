"""file open/write/read/close ([P] file), the basics to read and write text.

Binary files (file open ..., binary; formats in core/binfmt.py) support file read, file write,
file seek and file set byteorder; file query counts the open files.

file write directives: "text", `"text"', (expression), %fmt (expression),
_newline[(#)] (_n), _tab[(#)], _skip[(#)], _char(#), _dup(#).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from ..core import binfmt
from ..core.errors import StataError
from ..core.formats import DISPLAY_NUMERIC, format_value, parse_format
from ..lang.syntax import match_options
from ..lang.words import strip_outer_quotes
from .display import _literal
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_DIR = re.compile(r"_(newline|n|tab|skip|char|dup)\b(\s*\(\s*([^)]*)\s*\))?")
_FMT = re.compile(r"%-?~?0?\d*(\.\d+)?[a-zA-Z]+c?")


class Handle:
    def __init__(self, path: Path, mode: str, fh, binary: bool = False):
        self.path = path
        self.mode = mode          # "r", "w" or "rw"
        self.fh = fh
        self.binary = binary
        self.byteorder = binfmt.NATIVE_ORDER


def _handles(s: "Session") -> dict[str, Handle]:
    if not hasattr(s, "file_handles"):
        s.file_handles = {}
    return s.file_handles


def _get(s: "Session", name: str) -> Handle:
    h = _handles(s).get(name)
    if h is None:
        raise StataError(603, f"file handle {name} not found")   # VERIFY
    return h


def _open(s: "Session", rest: str) -> None:
    head, comma, opts_text = rest.partition(",")
    m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s+using\s+(.+)$", head.strip())
    if not m:
        raise StataError(198, "invalid syntax")
    name, file = m.group(1), strip_outer_quotes(m.group(2).strip())
    o = match_options(opts_text, {"read": 1, "write": 1, "text": 1, "binary": 1,
                                  "replace": 3, "append": 3, "all": 3}) if comma else {}
    if name in _handles(s):
        raise StataError(110, f"file handle {name} already exists")   # VERIFY
    path = Path(file).expanduser()
    binary = bool(o.get("binary"))
    if o.get("read") and o.get("write"):
        # read write: opens an existing file or creates a new one, positioned at the top
        fh = open(path, "r+b" if path.exists() else "w+b") if binary else \
            open(path, "r+" if path.exists() else "w+", encoding="utf-8", newline="")
        mode = "rw"
    elif o.get("write"):
        if path.exists() and not (o.get("replace") or o.get("append")):
            raise StataError(602, f"file {file} already exists")
        if o.get("replace") and not path.exists():
            s.output.write(f"(note: file {file} not found)\n", "text")
        mode = "a" if o.get("append") else "w"
        fh = open(path, mode + "b") if binary else open(path, mode, encoding="utf-8", newline="")
        mode = "w"
    elif o.get("read"):
        if not path.exists():
            raise StataError(601, f"file {file} not found")
        fh = open(path, "rb") if binary else open(path, "r", encoding="utf-8", errors="replace", newline="")
        mode = "r"
    else:
        raise StataError(198, "option read or write required")
    _handles(s)[name] = Handle(path, mode, fh, binary)


def _write(s: "Session", rest: str) -> None:
    m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*(.*)$", rest, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    h = _get(s, m.group(1))
    if h.mode == "r":
        raise StataError(608, f"file {h.path} is read-only")   # VERIFY
    if h.binary:
        _write_binary(s, h, m.group(2))
        return
    text = m.group(2)
    out: list[str] = []
    fmt = None
    dup = 1
    i, n = 0, len(text)
    while i < n:
        while i < n and text[i] in " \t":
            i += 1
        if i >= n:
            break
        rest_ = text[i:]
        d = _DIR.match(rest_)
        if d:
            kind, arg = d.group(1), d.group(3)
            num = int(s.eval(arg)) if arg and arg.strip() else None
            i += d.end()
            if kind in ("newline", "n"):
                out.append("\n" * (num if num is not None else 1))
            elif kind == "tab":
                out.append("\t" * (num if num is not None else 1))
            elif kind == "skip":
                out.append(" " * (num if num is not None else 1))
            elif kind == "char":
                out.append(chr(num or 0) * dup)
                dup = 1
            elif kind == "dup":
                dup = max(num or 0, 0)
            continue
        f = _FMT.match(rest_)
        if f and rest_.startswith("%"):
            fmt = parse_format(f.group(0))
            i += f.end()
            continue
        lit = _literal(rest_)
        if lit is not None:
            value, used = lit
            i += used
            out.append((format_value(value, fmt) if fmt is not None and fmt.kind == "s" else value) * dup)
            fmt, dup = None, 1
            continue
        if rest_.startswith("("):
            depth, j = 0, 0
            while j < len(rest_):
                ch = rest_[j]
                if ch == '"':
                    k = rest_.find('"', j + 1)
                    j = len(rest_) - 1 if k == -1 else k
                elif ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            value = s.eval(rest_[1:j])
            i += j + 1
            if isinstance(value, str):
                piece = format_value(value, fmt) if fmt is not None and fmt.kind == "s" else value
            else:
                piece = format_value(value, fmt if fmt is not None else DISPLAY_NUMERIC, pad=fmt is not None)
                piece = piece if fmt is not None else piece.strip()
            out.append(piece * dup)
            fmt, dup = None, 1
            continue
        raise StataError(198, "invalid syntax")
    h.fh.write("".join(out))


def _write_binary(s: "Session", h: Handle, text: str) -> None:
    """%fmt (exp), %#s "text", %#s `"text"', %#s (exp)."""
    i, n = 0, len(text)
    while i < n:
        while i < n and text[i].isspace():
            i += 1
        if i >= n:
            break
        f = re.match(r"%\d*[a-z]+", text[i:])
        if not f:
            raise StataError(198, "invalid syntax")
        try:
            spec = binfmt.parse_spec(f.group(0))
        except ValueError:
            raise StataError(198, f"{f.group(0)} invalid format")   # VERIFY
        i += f.end()
        while i < n and text[i].isspace():
            i += 1
        rest_ = text[i:]
        lit = _literal(rest_)
        if lit is not None:
            value, used = lit
            i += used
        elif rest_.startswith("("):
            depth, j = 0, 0
            while j < len(rest_):
                ch = rest_[j]
                if ch == '"':
                    k = rest_.find('"', j + 1)
                    j = len(rest_) - 1 if k == -1 else k
                elif ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            value = s.eval(rest_[1:j])
            i += j + 1
        else:
            raise StataError(198, "invalid syntax")
        if (spec.kind == "s") != isinstance(value, str):
            raise StataError(109, "type mismatch")
        h.fh.write(binfmt.encode(value, spec, h.byteorder))


def _read_binary(s: "Session", h: Handle, rest: str) -> None:
    """file read h %fmt name: numbers go to a scalar, strings to a local macro."""
    specs = re.findall(r"(%\d*[a-z]+)\s+(\S+)", rest)
    if not specs:
        raise StataError(198, "invalid syntax")
    eof = False
    for fmt, name in specs:
        try:
            spec = binfmt.parse_spec(fmt)
        except ValueError:
            raise StataError(198, f"{fmt} invalid format")   # VERIFY
        raw = h.fh.read(spec.size)
        if len(raw) < spec.size:        # VERIFY: a partial value at the end counts as eof
            eof = True
            break
        value = binfmt.decode(raw, spec, h.byteorder)
        if spec.kind == "s":
            s.macros.set_local(name, value)
        else:
            s.scalars[name] = value
    s.r = {"eof": 1.0 if eof else 0.0}


def _seek(s: "Session", rest: str) -> None:
    parts = rest.split()
    if len(parts) != 2:
        raise StataError(198, "invalid syntax")
    h = _get(s, parts[0])
    where = parts[1]
    if where == "query":
        pass
    elif where == "tof":
        h.fh.seek(0)
    elif where == "eof":
        h.fh.seek(0, 2)
    else:
        try:
            pos = int(float(s.eval(where)))
        except (ValueError, TypeError):
            raise StataError(198, "invalid syntax")
        if pos < 0:
            raise StataError(198, "invalid file position")   # VERIFY
        h.fh.seek(pos)
    s.r = {"loc": float(h.fh.tell())}


def _set(s: "Session", rest: str) -> None:
    parts = rest.split()
    if len(parts) != 3 or parts[1] != "byteorder":
        raise StataError(198, "invalid syntax")
    h = _get(s, parts[0])
    order = {"hilo": binfmt.HILO, "1": binfmt.HILO, "lohi": binfmt.LOHI, "2": binfmt.LOHI}.get(parts[2])
    if order is None:
        raise StataError(198, "byteorder must be hilo, lohi, 1 or 2")
    h.byteorder = order


def _read(s: "Session", rest: str) -> None:
    parts = rest.split()
    h = _get(s, parts[0]) if parts else None
    if h is not None and h.binary:
        if h.mode == "w":
            raise StataError(608, f"file {h.path} is write-only")   # VERIFY
        _read_binary(s, h, rest.split(None, 1)[1] if len(parts) > 1 else "")
        return
    if len(parts) != 2:
        raise StataError(198, "invalid syntax")
    if h.mode == "w":
        raise StataError(608, f"file {h.path} is write-only")   # VERIFY
    line = h.fh.readline()
    eof = line == ""
    s.macros.set_local(parts[1], line.rstrip("\r\n"))
    s.r = {"eof": 1.0 if eof else 0.0}


def _close(s: "Session", rest: str) -> None:
    names = rest.split()
    handles = _handles(s)
    if names in (["_all"], ["all"]):
        names = list(handles)
    for name in names:
        h = _get(s, name)
        h.fh.close()
        del handles[name]


@command("file")
def cmd_file(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    if sub in ("open",):
        _open(s, rest)
    elif sub in ("write", "w"):
        from ..core.formats import plain_numbers
        with plain_numbers():          # the file stores the number as Stata does
            _write(s, rest)
    elif sub in ("read", "r"):
        _read(s, rest)
    elif sub == "close":
        _close(s, rest)
    elif sub == "seek":
        _seek(s, rest)
    elif sub == "set":
        _set(s, rest)
    elif sub == "query":
        s.r = {"N": float(len(_handles(s)))}
    else:
        raise StataError(198, f"file {sub}: subcommand not yet supported by OpenDTA")
