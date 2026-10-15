"""Mata library ([M-5]): the functions Mata shares with Stata's expressions (sqrt, ln, normal,
substr, strupper...), applied element by element to c-conformable arguments.

They are the same functions as in generate and replace, evaluated the same way
(lang/vexpr.py): the vectorized form when there is one and every argument is real, otherwise
lang/functions.py once per distinct combination of arguments.
"""

from __future__ import annotations

import math

import numpy as np

from ..core.errors import StataError
from ..lang import functions as F
from ..lang import vexpr
from .libbase import lib, real_arg
from .values import MV, SYS, clean, conformability, real, type_mismatch


def _broadcast(vals: list[MV]) -> tuple[list[np.ndarray], tuple]:
    shapes = [v.a.shape for v in vals]
    for x in shapes:
        for y in shapes:
            for dx, dy in zip(x, y):
                if dx != dy and dx != 1 and dy != 1:
                    raise conformability()
    shape = np.broadcast_shapes(*shapes)
    return [np.broadcast_to(v.a, shape) for v in vals], shape


def _elementwise(fname: str):
    fast = vexpr._FAST.get(fname)

    def fn(eng, vals, refs):
        arrays, shape = _broadcast(vals)
        if not np.prod(shape):
            return MV(np.zeros(shape), "real")
        flat = [a.ravel() for a in arrays]
        if fast is not None and all(v.t == "real" for v in vals):
            return real(np.asarray(fast(*[f.astype(np.float64) for f in flat])).reshape(shape))
        try:
            out = vexpr._elementwise(fname, [f.astype(object) if v.t == "string" else f.astype(np.float64)
                                             for f, v in zip(flat, vals)], len(flat[0]))
        except StataError as e:
            if e.rc == 109:
                raise type_mismatch()
            raise
        if isinstance(out[0], str):
            return MV(np.array(out, dtype=object).reshape(shape), "string")
        return real(np.array(out, dtype=np.float64).reshape(shape))
    return fn


for _name in ["abs", "sqrt", "exp", "ln", "log", "log10", "sin", "cos", "tan", "asin", "acos", "atan",
              "atan2", "sinh", "cosh", "tanh", "floor", "ceil", "trunc", "round", "sign", "mod", "comb",
              "lnfactorial", "lngamma", "digamma", "normal", "normalden", "invnormal", "lnnormal",
              "chi2", "chi2tail", "invchi2", "invchi2tail", "t", "ttail", "invt", "invttail", "F",
              "Ftail", "invF", "invFtail", "binomial", "binomialtail", "poisson", "poissontail",
              "ibeta", "betaden", "gammap", "logit", "invlogit", "strlen", "ustrlen", "substr",
              "strupper", "strlower", "strproper", "strtrim", "strltrim", "strrtrim", "stritrim",
              "strpos", "subinstr", "subinword", "strreverse", "abbrev", "regexm", "regexr",
              "regexs", "strmatch", "ustrupper", "ustrlower", "ustrtitle", "date", "mdy", "year",
              "month", "day", "dow", "doy", "quarter"]:
    _nargs = F.FUNCTIONS[_name][1:]
    lib(_name, _nargs[0], _nargs[1])(_elementwise(_name))


@lib("factorial", 1)
def _factorial(eng, v, r):
    a = real_arg(v[0])
    with np.errstate(all="ignore"):
        out = np.where((a >= 0) & (a == np.trunc(a)) & (a < SYS),
                       np.vectorize(lambda x: math.gamma(x + 1) if x < 171 else np.inf)(np.where(a < SYS, a, 0)),
                       np.nan)
    return real(clean(out))
