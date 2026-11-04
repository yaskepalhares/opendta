"""The tools of ml evaluators ([R] ml, [P] mleval, [P] mlsum): mleval, mlsum, mlvecsum and
mlmatsum, which programs of method d0, d1 and d2 call to build the likelihood and its
derivatives."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from .mlmodel import Eq, Model, _full, _mats

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _tool_split(args: str, opts: dict):
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    return body, (match_options(opt, opts) if opt.strip() else {})


def _eq_index(m: Model, v) -> int:
    k = int(v) if isinstance(v, str) else 1
    if not 1 <= k <= m.K:
        raise StataError(198, f"equation {k} not found")
    return k - 1


def cmd_mleval(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"eq": 2, "scalar": 6})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mleval can only be used within an ml evaluator")
    mm = re.fullmatch(r"\s*(\S+?)\s*=\s*(\S+)\s*", body)
    if not mm:
        raise StataError(198, "invalid syntax")
    new, bname = mm.groups()
    e = m.eqs[_eq_index(m, o.get("eq"))]
    bvec = np.asarray(_mats(s)[bname].data, dtype=float).ravel()
    coef = bvec[e.start:e.start + e.k]
    if o.get("scalar"):
        s.scalars[new] = float(coef[-1] if e.cons else 0.0)
        return
    from ..core.dataset import Variable
    xb = e.X @ coef + e.off
    full = _full(m, xb, s.data.nobs)
    if s.data.has(new):
        s.data.get(new).data = full
    else:
        s.data.add(Variable(new, "double", full))


def _expr(s: "Session", m: Model, text: str):
    mm = re.fullmatch(r"\s*(.+?)\s*", text)
    body, _, ifp = mm.group(1).partition(" if ")
    from ..commands._util import eval_vector
    v = np.broadcast_to(np.asarray(eval_vector(s, body), dtype=float), (s.data.nobs,))
    sel = m.touse.copy()
    if ifp.strip():
        from ..lang.vexpr import truth
        t = np.broadcast_to(np.asarray(truth(eval_vector(s, ifp)), dtype=bool), (s.data.nobs,))
        sel &= t
    return v, sel


def cmd_mlsum(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"noweight": 3})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mlsum can only be used within an ml evaluator")
    name, eq_sign, expr = body.partition("=")
    if not eq_sign:
        raise StataError(198, "invalid syntax")
    v, sel = _expr(s, m, expr)
    x = v[sel]
    if np.any(x >= SYS):
        s.scalars[name.strip()] = SYS
        return
    wfull = np.ones(s.data.nobs)
    if m.w is not None and not o.get("noweight"):
        wfull = np.zeros(s.data.nobs)
        wfull[m.touse] = m.w
    s.scalars[name.strip()] = float((x * wfull[sel]).sum())


def _x_full(m: Model, e: Eq, nobs: int) -> np.ndarray:
    X = np.zeros((nobs, e.k))
    X[m.touse] = e.X
    return X


def cmd_mlvecsum(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"eq": 2})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mlvecsum can only be used within an ml evaluator")
    left, eq_sign, expr = body.partition("=")
    names = left.split()
    if len(names) != 2 or not eq_sign:
        raise StataError(198, "invalid syntax")
    lnf, vec = names
    e = m.eqs[_eq_index(m, o.get("eq"))]
    v, sel = _expr(s, m, expr)
    x = v[sel]
    if np.any(x >= SYS):
        s.scalars[lnf] = SYS
        return
    wfull = np.ones(s.data.nobs)
    if m.w is not None:
        wfull = np.zeros(s.data.nobs)
        wfull[m.touse] = m.w
    X = _x_full(m, e, s.data.nobs)[sel]
    from ..commands.matrix import Matrix
    row = (x * wfull[sel]) @ X
    _mats(s)[vec] = Matrix(row.reshape(1, -1), ["r1"],
                           [c for q, c in m.names if q == e.name])


def cmd_mlmatsum(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"eq": 2})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mlmatsum can only be used within an ml evaluator")
    left, eq_sign, expr = body.partition("=")
    names = left.split()
    if len(names) != 2 or not eq_sign:
        raise StataError(198, "invalid syntax")
    lnf, mat = names
    ev = o.get("eq")
    if isinstance(ev, str) and "," in ev:
        i1, i2 = (int(x) for x in ev.split(","))
    else:
        i1 = int(ev) if isinstance(ev, str) else 1
        i2 = i1
    e1, e2 = m.eqs[i1 - 1], m.eqs[i2 - 1]
    v, sel = _expr(s, m, expr)
    x = v[sel]
    if np.any(x >= SYS):
        s.scalars[lnf] = SYS
        return
    wfull = np.ones(s.data.nobs)
    if m.w is not None:
        wfull = np.zeros(s.data.nobs)
        wfull[m.touse] = m.w
    X1 = _x_full(m, e1, s.data.nobs)[sel]
    X2 = _x_full(m, e2, s.data.nobs)[sel]
    from ..commands.matrix import Matrix
    _mats(s)[mat] = Matrix(X1.T @ (X2 * (x * wfull[sel])[:, None]),
                           [c for q, c in m.names if q == e1.name], [c for q, c in m.names if q == e2.name])
