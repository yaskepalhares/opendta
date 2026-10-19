"""Binary formats of `file read` and `file write` ([P] file): %1b, %2b, %4b (Stata integers
with missing values), %1bs/%1bu ... (plain signed and unsigned integers), %4z and %8z (floats)
and %#s (strings of fixed width).

Plain functions on numbers and bytes, with no Session. Numbers use OpenDTA's double encoding of
missing values (core/missing.py); the codes on disk are the ones of the .dta specification.
"""

from __future__ import annotations

import re
import struct
import sys
from dataclasses import dataclass

from . import missing as M

HILO, LOHI = 1, 2
NATIVE_ORDER = LOHI if sys.byteorder == "little" else HILO

_SPEC = re.compile(r"%(?:(\d+)s|([124])b([su]?)|([48])z)$")

# first missing code (.) and the largest valid value of the integer types of the .dta format
_INT_MISS = {1: (101, 100, -127), 2: (32741, 32740, -32767), 4: (2147483621, 2147483620, -2147483647)}
_FLOAT_DOT = 0x7F000000
_FLOAT_STEP = 0x00000800


@dataclass(frozen=True)
class Spec:
    kind: str        # "s" string, "int" integer, "float" float or double
    size: int        # bytes
    sign: str = ""   # "", "s" or "u" for the integers


def parse_spec(fmt: str) -> Spec:
    m = _SPEC.match(fmt.strip())
    if not m or (m.group(1) is not None and int(m.group(1)) < 1):
        raise ValueError(fmt)
    if m.group(1) is not None:
        return Spec("s", int(m.group(1)))
    if m.group(2) is not None:
        return Spec("int", int(m.group(2)), m.group(3))
    return Spec("float", int(m.group(4)))


def _end(order: int) -> str:
    return ">" if order == HILO else "<"


def _ext_index(x: float) -> int:
    """0 for ., 1 for .a ... 26 for .z."""
    name = M.missing_name(x)
    return 0 if name == "." else ord(name[1]) - 96


def encode(value, spec: Spec, order: int = NATIVE_ORDER) -> bytes:
    """Bytes of a value; strings are cut or padded with binary zeros to the width."""
    if spec.kind == "s":
        raw = str(value).encode("utf-8")[:spec.size]
        return raw.ljust(spec.size, b"\0")
    x = float(value)
    miss = M.is_missing(x)
    e = _end(order)
    if spec.kind == "float":
        if spec.size == 8:
            return struct.pack(e + "d", x)
        if miss:
            return struct.pack(e + "I", _FLOAT_DOT + _ext_index(x) * _FLOAT_STEP)
        if abs(x) > 1.7e38:                                  # out of range: written as missing
            return struct.pack(e + "I", _FLOAT_DOT)
        return struct.pack(e + "f", x)
    n = spec.size
    if spec.sign:
        lo, hi = ((-(2 ** (8 * n - 1)) + 1, 2 ** (8 * n - 1) - 1) if spec.sign == "s"
                  else (0, 2 ** (8 * n) - 1))
        v = hi if miss else min(max(int(x), lo), hi)         # out of range: clipped
        return v.to_bytes(n, "big" if order == HILO else "little", signed=spec.sign == "s")
    first, top, bottom = _INT_MISS[n]
    if miss:
        v = first + _ext_index(x)
    else:
        v = int(x)                                           # noninteger values are truncated
        if v > top or v < bottom:
            v = first
    return v.to_bytes(n, "big" if order == HILO else "little", signed=True)


def decode(raw: bytes, spec: Spec, order: int = NATIVE_ORDER):
    """Value of the bytes (a string or a number, with Stata's missing values)."""
    if spec.kind == "s":
        return raw.split(b"\0", 1)[0].decode("utf-8", errors="replace")
    e = _end(order)
    if spec.kind == "float":
        if spec.size == 8:
            x = struct.unpack(e + "d", raw)[0]
            return M.normalize(x)
        bits = struct.unpack(e + "I", raw)[0]
        if bits >= _FLOAT_DOT:
            k = (bits - _FLOAT_DOT) // _FLOAT_STEP
            return M.SYSMISS if k == 0 or k > 26 else M.EXTENDED[chr(96 + k)]
        return float(struct.unpack(e + "f", raw)[0])
    n = spec.size
    big = order == HILO
    if spec.sign:
        return float(int.from_bytes(raw, "big" if big else "little", signed=spec.sign == "s"))
    v = int.from_bytes(raw, "big" if big else "little", signed=True)
    first, top, _ = _INT_MISS[n]
    if v > top:
        k = v - first
        return M.SYSMISS if k <= 0 or k > 26 else M.EXTENDED[chr(96 + k)]
    return float(v)
