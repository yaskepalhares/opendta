"""Mata library ([M-5]): strings and output (sprintf, printf, display)."""

from __future__ import annotations

import re

import numpy as np

from ..core import missing as M
from ..core.formats import format_value
from .libbase import int_arg, lib, real_arg, str_arg
from .values import MV, SYS, MataError, real, string

# ---------------------------------------------------------------------------
# strings
# ---------------------------------------------------------------------------

def _plain(x: float, fmt: str) -> str:
    """Number as text that becomes data: no superscript exponent."""
    from ..core.formats import plain_numbers
    with plain_numbers():
        return format_value(x, fmt, pad=False).strip()


@lib("strofreal", 1, 2)
def _strofreal(eng, v, r):
    fmt = v[1].str_scalar() if len(v) > 1 else "%9.0g"
    a = real_arg(v[0])
    out = np.empty(a.shape, dtype=object)
    for idx in np.ndindex(a.shape):
        x = float(a[idx])
        out[idx] = M.missing_name(x) if x >= SYS else _plain(x, fmt)
    return MV(out, "string")


@lib("strtoreal", 1)
def _strtoreal(eng, v, r):
    a = str_arg(v[0])
    out = np.empty(a.shape)
    for idx in np.ndindex(a.shape):
        t = str(a[idx]).strip()
        try:
            out[idx] = M.missing_code(t) if t.startswith(".") and len(t) <= 2 and not t[1:].isdigit() \
                else float(t)
        except ValueError:
            out[idx] = SYS
    return real(out)


@lib("tokens", 1, 2)
def _tokens(eng, v, r):
    s = v[0].str_scalar()
    if len(v) > 1 and v[1] is not None:
        seps = v[1].str_scalar()
        out, cur = [], ""
        for ch in s:
            if ch == " ":
                if cur:
                    out.append(cur)
                cur = ""
            elif ch in seps:
                if cur:
                    out.append(cur)
                out.append(ch)
                cur = ""
            else:
                cur += ch
        if cur:
            out.append(cur)
    else:
        from ..lang.words import split_words
        out = split_words(s)
    arr = np.empty((1, len(out)), dtype=object)
    for k, w in enumerate(out):
        arr[0, k] = w
    return MV(arr, "string")


@lib("invtokens", 1, 2)
def _invtokens(eng, v, r):
    sep = v[1].str_scalar() if len(v) > 1 else " "
    return string(sep.join(str(x) for x in str_arg(v[0]).ravel()))


@lib("strdup", 2)
def _strdup(eng, v, r):
    return string(v[0].str_scalar() * max(0, int_arg(v[1])))


@lib("ascii", 1)
def _ascii(eng, v, r):
    return real(np.array([[float(b) for b in v[0].str_scalar().encode("latin-1", "replace")]]).reshape(1, -1))


@lib("char", 1)
def _char(eng, v, r):
    return string("".join(chr(int(x)) for x in real_arg(v[0]).ravel() if x < SYS))


@lib("strcat", 2)
def _strcat(eng, v, r):
    from .ops import binary
    return binary("+", v[0], v[1])


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------

_SMCL_TAGS = re.compile(r"\{(txt|text|res|result|err|error|inp|input|com|sf|it|bf|ul off|ul on)\}")


def _smcl_strip(text: str) -> str:
    text = _SMCL_TAGS.sub("", text)
    text = re.sub(r"\{hline (\d+)\}", lambda m: "-" * int(m.group(1)), text)
    text = text.replace("{c |}", "|").replace("{c -}", "-").replace("{c +}", "+")
    return text


_FMT = re.compile(r"%(-|~)?(0)?(\d*)(?:\.(\d+))?(f|g|e|s|gc|fc|x|t[a-zA-Z]*)")


def sprintf(fmt: str, args: list[MV]) -> str:
    out = []
    i = 0
    k = 0
    while i < len(fmt):
        ch = fmt[i]
        if ch == "\\" and i + 1 < len(fmt):
            nxt = fmt[i + 1]
            out.append({"n": "\n", "t": "\t", "\\": "\\", "r": "\r"}.get(nxt, "\\" + nxt))
            i += 2
            continue
        if ch == "%":
            if fmt.startswith("%%", i):
                out.append("%")
                i += 2
                continue
            m = _FMT.match(fmt, i)
            if not m:
                out.append(ch)
                i += 1
                continue
            if k >= len(args):
                raise MataError(3001, "too few arguments")
            arg = args[k]
            k += 1
            spec = m.group(0)
            just, zero, width, dec, kind = m.groups()
            x = arg.scalar() if arg.is_scalar else None
            if x is None:
                raise MataError(3204, f"{arg.orgtype()} found where scalar required")
            if kind == "s":
                if arg.t != "string":
                    raise MataError(3254, "nonstring found where string required")
                text = str(x)
                if width:
                    text = text.ljust(int(width)) if just == "-" else text.rjust(int(width))
                out.append(text)
            else:
                if arg.t != "real":
                    raise MataError(3253, "nonreal found where real required")
                use = spec
                if not width:
                    use = "%" + (just or "") + "9" + ("." + dec if dec else ".0") + kind if kind != "g" else \
                        "%9.0g"
                    text = format_value(float(x), use, pad=False).strip()
                else:
                    text = format_value(float(x), use, pad=True)
                out.append(text)
            i = m.end()
            continue
        out.append(ch)
        i += 1
    return "".join(out)


@lib("printf", 1, 99)
def _printf(eng, v, r):
    text = _smcl_strip(sprintf(v[0].str_scalar(), v[1:]))
    eng.s.output.write(text, "result")
    return None


@lib("sprintf", 1, 99)
def _sprintf(eng, v, r):
    from ..core.formats import plain_numbers
    with plain_numbers():              # becomes data: no superscript exponent
        return string(sprintf(v[0].str_scalar(), v[1:]))


@lib("errprintf", 1, 99)
def _errprintf(eng, v, r):
    eng.s.output.write(_smcl_strip(sprintf(v[0].str_scalar(), v[1:])), "error")
    return None


@lib("display", 1, 2)
def _display(eng, v, r):
    a = str_arg(v[0])
    for x in a.ravel():
        text = str(x) if len(v) > 1 and v[1] is not None and v[1].real_scalar() else _smcl_strip(str(x))
        eng.s.output.write(text + "\n", "result")
    return None
