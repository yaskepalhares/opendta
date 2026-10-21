"""Shared helpers of the multivariate commands: reading the estimation sample, posting e() for
commands that have no e(b)/e(V), and the small formatting functions of their tables."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lazy import lazy
from ._util import touse
from .summarize import weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def numeric_names(s: "Session", varlist: str, minimum: int = 1) -> list[str]:
    names = s.expand_varlist(varlist)
    if len(names) < minimum:
        raise StataError(102, "too few variables specified")
    for n in names:
        if s.data.get(n).is_string:
            raise StataError(109, "type mismatch")
    return names


def sample(s: "Session", p, names: list[str], allowed=("fweight", "aweight"), extra=()):
    """(X, w, weight type, mask): the rows with the variables (and the `extra` ones) present."""
    ds = s.data
    mask = touse(s, p)
    for n in list(names) + list(extra):
        mask &= ds.get(n).data < SYS
    w, wtype, mask = weights(s, p, mask, allowed)
    if not mask.any():
        raise StataError(2000, "no observations")
    X = np.column_stack([ds.get(n).data[mask] for n in names]).astype(float)
    return X, (None if w is None else w[mask]), wtype, mask


def matrix(data, rows=None, cols=None):
    from .matrix import Matrix
    d = np.atleast_2d(np.asarray(data, dtype=float))
    return Matrix(d, list(rows or []), list(cols or []))


def post_e(s: "Session", cmd: str, scalars, macros, mats, mask, N: int, *, estat=None, predict=None,
           predict_multi=None, predict_opts=None, extra=None) -> None:
    """Stores e() for a command without e(b) and e(V). `scalars`, `macros` and `mats` are lists of
    (name, value) in the order `ereturn list` must show them; the dictionary is filled in the
    reverse order, as `estimation.results.post` does."""
    from ..estimation.results import Estimates
    e: dict = {}
    for k, v in reversed(macros):
        if v != "":
            e[k] = v
    for k, v in reversed(scalars):
        e[k] = float(v) if v is not None else SYS
    for k, m in reversed(mats):
        e[k] = m
    s.e = e
    s.data.esample = np.asarray(mask, dtype=bool).copy()
    est = Estimates(cmd=cmd, depvar="", names=[], b=np.zeros(0), V=np.zeros((0, 0)), N=int(N))
    est.estat = estat or {}
    est.predict = predict
    est.predict_multi = predict_multi or {}
    est.predict_opts = predict_opts or {}
    est.extra = extra or {}
    est.no_b = True
    s._estimates = est


def g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def f(x: float, width: int, dec: int) -> str:
    if x is None or not np.isfinite(x):
        return f"{'.':>{width}}"
    return f"{x:{width}.{dec}f}"


def chi2_p(x: float, df: float) -> float:
    if not np.isfinite(x) or df <= 0:
        return float("nan")
    return float(st.chi2.sf(x, df))
