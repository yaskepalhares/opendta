"""What every module of the Mata library ([M-5]) shares: the table of functions and the
argument checks.

Each function receives (engine, values, references) and returns an MV or None (void
functions). `refs[k]` is a reference to the variable passed as the k-th argument, used by the
functions that write results into their arguments (symeigensystem, svd, _sort...).
"""

from __future__ import annotations

import numpy as np

from .values import MV, SYS, MataError

# name -> (function, minimum and maximum number of arguments, positions written through refs)
LIBRARY: dict[str, tuple] = {}


def lib(name: str, lo: int, hi: int | None = None, outs: tuple = ()):
    def deco(fn):
        LIBRARY[name] = (fn, lo, lo if hi is None else hi, set(outs))
        return fn
    return deco


def real_arg(v: MV) -> np.ndarray:
    if v is None:
        raise MataError(3499, "argument missing")
    if v.t != "real":
        raise MataError(3253, "nonreal found where real required")
    return v.a


def str_arg(v: MV) -> np.ndarray:
    if v.t != "string":
        raise MataError(3254, f"{v.t} found where string required")
    return v.a


def int_arg(v: MV) -> int:
    return v.int_scalar()


def nan_missing(a: np.ndarray) -> np.ndarray:
    """Missing values as NaN (for numpy functions that skip NaN)."""
    return np.where(a >= SYS, np.nan, a)


def zero_missing(a: np.ndarray) -> np.ndarray:
    return np.where(a >= SYS, 0.0, a)
