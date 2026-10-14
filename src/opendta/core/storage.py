"""Native storage of numeric variables.

Each variable is kept in memory in its own type, with 1 byte per value for byte,
2 for int, 4 for long and float, 8 for double, as in Stata. Missing values use the
same encoding as the .dta format, which makes reading and writing files a
direct copy:

    byte   int8     .  = 101,         .a..z = 102..127
    int    int16    .  = 32741,       .a..z = 32742..32767
    long   int32    .  = 2147483621,  .a..z = 2147483622..2147483647
    float  float32  .  = bits 0x7f000000, .a..z = + k * 0x800
    double float64  the same codes used in the rest of OpenDTA (missing.py)

Computations stay in float64: `decode` produces the float64 vector with
OpenDTA's missing codes and `encode` goes the other way.
"""

from __future__ import annotations

import numpy as np

from . import missing as M

DTYPE = {"byte": np.dtype(np.int8), "int": np.dtype(np.int16), "long": np.dtype(np.int32),
         "float": np.dtype(np.float32), "double": np.dtype(np.float64)}
INT_MISS = {"byte": 101, "int": 32741, "long": 2147483621}
INT_RANGES = {"byte": (-127, 100), "int": (-32767, 32740), "long": (-2147483647, 2147483620)}
FLOAT_MISS_BITS = 0x7F000000
FLOAT_MISS_STEP = 0x00000800
FLOAT_MAX = 1.7014117331926443e38   # largest non-missing float (bits 0x7effffff)

_DBL_BASE = np.uint64(0x7FE0000000000000)
_DBL_SHIFT = np.uint64(40)                 # step of 0x0000010000000000 between .a, .b, ...
_CODES = np.array([M.SYSMISS] + [M.EXTENDED[c] for c in "abcdefghijklmnopqrstuvwxyz"])


def _missing_index(v: np.ndarray, miss: np.ndarray) -> np.ndarray:
    """0 for '.', 1..26 for .a..z (unknown values become '.')."""
    bits = v[miss].view(np.uint64)
    k = ((bits - _DBL_BASE) >> _DBL_SHIFT).astype(np.int64)
    known = (k >= 0) & (k <= 26) & (bits == _DBL_BASE + (k.astype(np.uint64) << _DBL_SHIFT))
    return np.where(known, k, 0)


def encode(values, vtype: str) -> np.ndarray:
    """float64 (OpenDTA codes) -> vector in the native type.

    Values outside the type's range become '.', and non-integers are truncated
    toward zero, as in fit_numeric; values normally arrive already
    adjusted."""
    v = np.asarray(values, dtype=np.float64)
    if vtype == "double":
        return np.array(v, dtype=np.float64, copy=True)
    miss = (v >= M.SYSMISS) | (v != v)
    if vtype == "float":
        over = ~miss & (np.abs(v) > FLOAT_MAX)
        out = np.where(miss | over, 0.0, v).astype(np.float32)
        allmiss = miss | over
        if allmiss.any():
            idx = np.zeros(len(v), dtype=np.int64)
            if miss.any():
                idx[miss] = _missing_index(np.where(v != v, M.SYSMISS, v), miss)
            bits = out.view(np.uint32)
            bits[allmiss] = (FLOAT_MISS_BITS + idx[allmiss] * FLOAT_MISS_STEP).astype(np.uint32)
        return out
    lo, hi = INT_RANGES[vtype]
    over = ~miss & ((v < lo) | (v > hi))
    allmiss = miss | over
    out = np.where(allmiss, 0.0, v).astype(DTYPE[vtype])
    if allmiss.any():
        idx = np.zeros(len(v), dtype=np.int64)
        if miss.any():
            idx[miss] = _missing_index(np.where(v != v, M.SYSMISS, v), miss)
        out[allmiss] = (INT_MISS[vtype] + idx[allmiss]).astype(DTYPE[vtype])
    return out


def decode(raw: np.ndarray, vtype: str) -> np.ndarray:
    """Vector in the native type -> float64 with OpenDTA's missing codes."""
    if vtype == "double":
        return np.array(raw, dtype=np.float64, copy=True)
    v = raw.astype(np.float64)
    if vtype == "float":
        bits = raw.view(np.uint32)
        # only the positive values above the largest valid float; with the sign bit (bit 31) it is
        # a number
        miss = (bits >= FLOAT_MISS_BITS) & (bits < 0x80000000)
        if miss.any():
            k = ((bits[miss] - FLOAT_MISS_BITS) // FLOAT_MISS_STEP).astype(np.int64).clip(0, 26)
            v[miss] = _CODES[k]
        return v
    base = INT_MISS[vtype]
    miss = raw >= base
    if miss.any():
        v[miss] = _CODES[(raw[miss].astype(np.int64) - base).clip(0, 26)]
    return v


def decode_one(x, vtype: str) -> float:
    """A single value (list, Data Browser, sequential replace): without creating vectors."""
    if vtype == "double":
        return float(x)
    if vtype == "float":
        bits = int(np.float32(x).view(np.uint32))
        if FLOAT_MISS_BITS <= bits < 0x80000000:
            return float(_CODES[min((bits - FLOAT_MISS_BITS) // FLOAT_MISS_STEP, 26)])
        return float(x)
    k = int(x)
    base = INT_MISS[vtype]
    return float(k) if k < base else float(_CODES[min(k - base, 26)])


def encode_one(x: float, vtype: str):
    return encode(np.array([x], dtype=np.float64), vtype)[0]


def missing_raw(vtype: str):
    """Native value of '.' in the type."""
    return encode_one(M.SYSMISS, vtype)
