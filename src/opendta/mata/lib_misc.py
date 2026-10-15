"""Mata library ([M-5]): errors, exit and utilities."""

from __future__ import annotations

import math

import numpy as np

from ..core import missing as M
from .libbase import int_arg, lib, real_arg
from .values import SYS, MataError, MataExit, clean, real, string

# ---------------------------------------------------------------------------
# errors and utilities
# ---------------------------------------------------------------------------

@lib("exit", 0, 1)
def _exit(eng, v, r):
    rc = int_arg(v[0]) if v else 0
    raise MataExit(rc)


@lib("_error", 1, 2)
def _error_(eng, v, r):
    if v[0].t == "string":
        raise MataError(3498, v[0].str_scalar())
    rc = int_arg(v[0])
    msg = v[1].str_scalar() if len(v) > 1 else ""
    raise MataError(rc, msg)


@lib("error", 1)
def _error(eng, v, r):
    rc = int_arg(v[0])
    # VERIFY: Mata prints Stata's message for the code
    return real(float(rc))


@lib("assert", 1)
def _assert(eng, v, r):
    if real_arg(v[0]).size == 0 or not np.all(real_arg(v[0]) != 0):
        raise MataError(3498, "assertion is false")
    return None


@lib("args", 0)
def _args(eng, v, r):
    if not eng.argc:
        raise MataError(3000, "args() used outside of a function")   # VERIFY
    return real(float(eng.argc[-1]))


@lib("c", 1)
def _c(eng, v, r):
    val = eng.s.creturn(v[0].str_scalar())
    return string(val) if isinstance(val, str) else real(float(val))


@lib("swap", 2, outs=(0, 1))
def _swap(eng, v, r):
    if r[0] and r[1]:
        a, b = r[0].get(), r[1].get()
        r[0].set(b)
        r[1].set(a)
    return None


@lib("st_isname", 1)
def _isname(eng, v, r):
    from ..core.dataset import NAME_RE, RESERVED
    s = v[0].str_scalar()
    return real(1.0 if NAME_RE.match(s) and s not in RESERVED else 0.0)


@lib("floatround", 1)
def _floatround(eng, v, r):
    return real(clean(np.where(real_arg(v[0]) >= SYS, SYS, real_arg(v[0]).astype(np.float32).astype(np.float64))))


@lib("epsilon", 1)
def _epsilon(eng, v, r):
    return real(np.abs(real_arg(v[0])) * 2.0 ** -52)


@lib("maxdouble", 0)
def _maxdouble(eng, v, r):
    return real(M.MAXDOUBLE)


@lib("mindouble", 0)
def _mindouble(eng, v, r):
    return real(M.MINDOUBLE)


@lib("smallestdouble", 0)
def _smallestdouble(eng, v, r):
    return real(2.0 ** -1022)


@lib("pi", 0)
def _pi(eng, v, r):
    return real(math.pi)
