"""The model of ml ([R] ml): its equations and the evaluation of the likelihood for each method
(lf, lf0-lf2, d0-d2, gf0-gf2), with the numerical derivatives the methods leave to ml."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from .ml import MLResult


if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# evaluators called as "prog todo b lnfj ...": per-observation log likelihood
_LFAMILY = ("lf0", "lf1", "lf2", "gf0", "gf1", "gf2")


@dataclass
class Eq:
    name: str
    depvars: list[str]
    indep: list[str]
    cons: bool = True
    offset: str = ""
    exposure: str = ""
    scalar: bool = False                # /name: constant only, a parameter of the likelihood
    X: np.ndarray | None = None         # design matrix over the sample (constant last)
    off: np.ndarray | None = None
    start: int = 0

    @property
    def k(self) -> int:
        return len(self.indep) + (1 if self.cons else 0)


@dataclass
class Model:
    method: str
    prog: str
    eqs: list[Eq]
    touse: np.ndarray
    w: np.ndarray | None
    wtype: str
    opts: dict
    cmdline: str
    cluster: np.ndarray | None = None
    clustvar: str = ""
    b0: np.ndarray | None = None
    result: MLResult | None = None
    V: np.ndarray | None = None
    N: int = 0
    names: list[tuple[str, str]] = field(default_factory=list)    # (equation, column)
    ysets: list[str] = field(default_factory=list)
    sampvar: str = ""
    wvar: str = ""

    @property
    def K(self) -> int:
        return len(self.eqs)

    @property
    def nb(self) -> int:
        return sum(e.k for e in self.eqs)


def _mats(s: "Session") -> dict:
    from ..commands.matrix import _store
    return _store(s)


def _theta(m: Model, b: np.ndarray) -> np.ndarray:
    th = np.empty((m.N, m.K))
    for j, e in enumerate(m.eqs):
        th[:, j] = e.X @ b[e.start:e.start + e.k] + e.off
    return th


def _setup(s: "Session", m: Model) -> None:
    """Globals and sample variable read by the evaluator programs."""
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    ds = s.data
    if not m.sampvar or not ds.has(m.sampvar):
        m.sampvar = _temp_name(s)
        ds.add(Variable(m.sampvar, "double", np.zeros(ds.nobs)))
    ds.get(m.sampvar).data = m.touse.astype(float)
    g = s.macros.globals
    g["ML_samp"] = m.sampvar
    for i, y in enumerate(m.ysets, start=1):
        g[f"ML_y{i}"] = y
    if m.w is not None:
        if not m.wvar or not ds.has(m.wvar):
            m.wvar = _temp_name(s)
            ds.add(Variable(m.wvar, "double", np.zeros(ds.nobs)))
        full = np.zeros(ds.nobs)
        full[m.touse] = m.w
        ds.get(m.wvar).data = full
        g["ML_w"] = m.wvar
        m.wvar_posted = m.wvar          # e(wvar) of ml, svy (the variable is dropped at the end)
    else:
        g.pop("ML_w", None)
    s._ml_active = m


def _full(m: Model, vec: np.ndarray, nobs: int) -> np.ndarray:
    out = np.full(nobs, SYS)
    out[m.touse] = vec
    return out


def _call(s: "Session", text: str) -> None:
    with s.output.quietly():
        s.interp._execute_expanded(text)


def _lnf_obs(s: "Session", m: Model, th: np.ndarray) -> np.ndarray:
    """lf: the program fills lnf for the given linear predictors; returns lnf per observation."""
    from .mlexp import _mlexp_obs
    if m.method == "mlexp":
        return _mlexp_obs(s, m, th)
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    ds = s.data
    names = getattr(m, "_tmp", None)
    if names is None or any(not ds.has(n) for n in names):
        names = []
        for _ in range(m.K + 1):
            nm = _temp_name(s)
            ds.add(Variable(nm, "double", np.full(ds.nobs, SYS)))
            names.append(nm)
        m._tmp = names
    for j in range(m.K):
        ds.get(names[j + 1]).data = _full(m, th[:, j], ds.nobs)
    ds.get(names[0]).data = np.full(ds.nobs, SYS)
    _call(s, f"{m.prog} " + " ".join(names))
    v = ds.get(names[0]).data[m.touse]
    return np.where(v >= SYS, np.nan, v)


def _total(m: Model, lnf: np.ndarray) -> float:
    if np.any(np.isnan(lnf)):
        return float("nan")
    return float((lnf if m.w is None else lnf * m.w).sum())


def _eval_d(s: "Session", m: Model, b: np.ndarray, todo: int):
    mats = _mats(s)
    from ..commands.matrix import Matrix
    bm = Matrix(b.reshape(1, -1), ["y1"], [c for _, c in m.names])
    bm.coleq = [q for q, _ in m.names]
    mats["__mlb"] = bm
    s.scalars["__mllnf"] = SYS
    extra = ""
    if todo >= 1:
        mats.pop("__mlg", None)
        extra = " __mlg"
    if todo >= 2:
        mats.pop("__mlH", None)
        extra += " __mlH"
    _call(s, f"{m.prog} {todo} __mlb __mllnf{extra}")
    ll = float(s.scalars.get("__mllnf", SYS))
    ll = float("nan") if ll >= SYS else ll
    g = H = None
    if todo >= 1 and "__mlg" in mats:
        g = np.asarray(mats["__mlg"].data, dtype=float).ravel()
    if todo >= 2 and "__mlH" in mats:
        H = np.asarray(mats["__mlH"].data, dtype=float)
    return ll, g, H


def _num_fun(f0, nb: int):
    """fun(b, todo) with fourth-order numerical gradient and Hessian of the scalar f0(b)."""
    def num_grad(b):
        g = np.zeros(nb)
        for j in range(nb):
            hj = 5e-3 * max(1.0, abs(b[j]))
            e = np.zeros(nb)
            e[j] = hj
            g[j] = (-f0(b + 2 * e) + 8 * f0(b + e) - 8 * f0(b - e) + f0(b - 2 * e)) / (12 * hj)
        return g

    def num_hess(b, f00):
        H = np.zeros((nb, nb))
        hs = [5e-3 * max(1.0, abs(b[j])) for j in range(nb)]
        for j in range(nb):
            e = np.zeros(nb)
            e[j] = hs[j]
            H[j, j] = (-f0(b + 2 * e) + 16 * f0(b + e) - 30 * f00 + 16 * f0(b - e) - f0(b - 2 * e)) / (12 * hs[j] ** 2)
            for k in range(j):
                def D(a):
                    ej = np.zeros(nb)
                    ek = np.zeros(nb)
                    ej[j], ek[k] = a * hs[j], a * hs[k]
                    return (f0(b + ej + ek) - f0(b + ej - ek) - f0(b - ej + ek) + f0(b - ej - ek)) / (4 * a * a * hs[j] * hs[k])
                H[j, k] = H[k, j] = (4 * D(1) - D(2)) / 3
        return H

    def fun(b, todo):
        ll = f0(b)
        if todo == 0:
            return (ll if np.isfinite(ll) else -np.inf), None, None
        if not np.isfinite(ll):
            return -np.inf, np.zeros(nb), np.full((nb, nb), np.nan)
        return ll, num_grad(b), num_hess(b, ll)
    return fun


def _d_fun(s: "Session", m: Model):
    """fun(b, todo) -> (ll, g, H) for maximize()."""
    nb = m.nb
    use_g = m.method in ("d1", "d2")
    use_H = m.method == "d2"
    if not use_g:
        return _num_fun(lambda b: _eval_d(s, m, b, 0)[0], nb)

    def fun(b, todo):
        if todo == 0:
            ll = _eval_d(s, m, b, 0)[0]
            return (ll if np.isfinite(ll) else -np.inf), None, None
        ll, g, H = _eval_d(s, m, b, 2 if use_H else 1)
        if g is None:
            raise StataError(198, f"{m.prog} did not return the gradient")   # VERIFY
        if not np.isfinite(ll):
            return -np.inf, g, np.full((nb, nb), np.nan)
        if H is None:
            H = np.zeros((nb, nb))
            for j in range(nb):
                hj = 5e-3 * max(1.0, abs(b[j]))
                e = np.zeros(nb)
                e[j] = hj
                gp, gm = _eval_d(s, m, b + e, 1)[1], _eval_d(s, m, b - e, 1)[1]
                gp2, gm2 = _eval_d(s, m, b + 2 * e, 1)[1], _eval_d(s, m, b - 2 * e, 1)[1]
                H[:, j] = (-gp2 + 8 * gp - 8 * gm + gm2) / (12 * hj)
            H = (H + H.T) / 2
        return ll, g, H
    return fun


def _lf_derivs(s: "Session", m: Model, th: np.ndarray):
    """Per-observation first and second derivatives of lnf with respect to the linear
    predictors (fourth-order central differences)."""
    n, K = th.shape
    h = 5e-3 * np.maximum(1.0, np.abs(th))

    def f(t):
        return _lnf_obs(s, m, t)

    def shift(cols):
        t = th.copy()
        for j, a in cols:
            t[:, j] = t[:, j] + a * h[:, j]
        return t
    f0 = f(th)
    G = np.zeros((n, K))
    H = np.zeros((n, K, K))
    for j in range(K):
        fp1, fm1 = f(shift([(j, 1)])), f(shift([(j, -1)]))
        fp2, fm2 = f(shift([(j, 2)])), f(shift([(j, -2)]))
        G[:, j] = (-fp2 + 8 * fp1 - 8 * fm1 + fm2) / (12 * h[:, j])
        H[:, j, j] = (-fp2 + 16 * fp1 - 30 * f0 + 16 * fm1 - fm2) / (12 * h[:, j] ** 2)
        for k in range(j):
            def D(a):
                return (f(shift([(j, a), (k, a)])) - f(shift([(j, a), (k, -a)])) - f(shift([(j, -a), (k, a)]))
                        + f(shift([(j, -a), (k, -a)]))) / (4 * a * a * h[:, j] * h[:, k])
            H[:, j, k] = H[:, k, j] = (4 * D(1) - D(2)) / 3
    return f0, G, H


def _lf_fun(s: "Session", m: Model):
    w = m.w if m.w is not None else np.ones(m.N)

    def fun(b, todo):
        th = _theta(m, b)
        if todo == 0:
            ll = _total(m, _lnf_obs(s, m, th))
            return (ll if np.isfinite(ll) else -np.inf), None, None
        lnf, G, Hh = _lf_derivs(s, m, th)
        ll = _total(m, lnf)
        if not np.isfinite(ll) or not np.all(np.isfinite(G)):
            return -np.inf, np.zeros(m.nb), np.full((m.nb, m.nb), np.nan)
        g = np.zeros(m.nb)
        H = np.zeros((m.nb, m.nb))
        for j, ej in enumerate(m.eqs):
            sj = slice(ej.start, ej.start + ej.k)
            g[sj] = ej.X.T @ (w * G[:, j])
            for k, ek in enumerate(m.eqs):
                sk = slice(ek.start, ek.start + ek.k)
                H[sj, sk] = ej.X.T @ (ek.X * (w * Hh[:, j, k])[:, None])
        m._G = G
        return ll, g, H
    return fun


def _eval_l(s: "Session", m: Model, b: np.ndarray, todo: int):
    """lf0, lf1, lf2, gf0, gf1, gf2: "prog todo b lnfj [g1 ... gK | S] [H]". Returns lnfj over the
    sample, the per-observation scores with respect to the linear predictors (lf1, lf2; n x K) or
    to the coefficients (gf1, gf2; S, n x nb), and the Hessian (lf2, gf2)."""
    from ..commands.matrix import Matrix
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    ds = s.data
    mats = _mats(s)
    bm = Matrix(b.reshape(1, -1), ["y1"], [c for _, c in m.names])
    bm.coleq = [q for q, _ in m.names]
    mats["__mlb"] = bm
    ng = m.K if m.method in ("lf1", "lf2") else 0
    names = getattr(m, "_tmp", None)
    if names is None or len(names) != 1 + ng or any(not ds.has(n) for n in names):
        names = []
        for _ in range(1 + ng):
            nm = _temp_name(s)
            ds.add(Variable(nm, "double", np.full(ds.nobs, SYS)))
            names.append(nm)
        m._tmp = names
    for nm in names:
        ds.get(nm).data = np.full(ds.nobs, SYS)
    args = f"{m.prog} {todo} __mlb {names[0]}"
    if m.method in ("lf1", "lf2"):
        args += " " + " ".join(names[1:])
    elif m.method in ("gf1", "gf2"):
        mats.pop("__mlS", None)
        args += " __mlS"
    if m.method in ("lf2", "gf2"):
        mats.pop("__mlH", None)
        args += " __mlH"
    _call(s, args)
    lnf = ds.get(names[0]).data[m.touse]
    if m.method.startswith("gf"):
        lnf = np.where(lnf >= SYS, 0.0, lnf)         # gf: rows without a contribution (VERIFY)
    lnf = np.where(lnf >= SYS, np.nan, lnf)
    G = S = H = None
    if todo >= 1 and m.method in ("lf1", "lf2"):
        G = np.column_stack([ds.get(nm).data[m.touse] for nm in names[1:]])
        G = np.where(G >= SYS, np.nan, G)
    if todo >= 1 and m.method in ("gf1", "gf2") and "__mlS" in mats:
        S = np.asarray(mats["__mlS"].data, dtype=float)
    if todo >= 2 and m.method in ("lf2", "gf2") and "__mlH" in mats:
        H = np.asarray(mats["__mlH"].data, dtype=float)
    return lnf, G, S, H


def _l_fun(s: "Session", m: Model):
    """fun(b, todo) for the lf0 ... gf2 evaluators; sets m._S, the weighted per-observation
    scores with respect to the coefficients, for the robust variance."""
    w = m.w if m.w is not None else np.ones(m.N)
    nb = m.nb

    def obs(b):
        return _eval_l(s, m, b, 0)[0]

    def num_scores(b):
        S = np.zeros((m.N, nb))
        for j in range(nb):
            hj = 1e-4 * max(1.0, abs(b[j]))
            e = np.zeros(nb)
            e[j] = hj
            S[:, j] = (obs(b + e) - obs(b - e)) / (2 * hj)
        return S

    def scores(b):
        if m.method in ("lf0", "gf0"):
            return num_scores(b)
        lnf, G, S, _ = _eval_l(s, m, b, 1)
        if m.method in ("lf1", "lf2"):
            if G is None:
                raise StataError(198, f"{m.prog} did not return the scores")      # VERIFY
            return np.column_stack([e.X * G[:, [j]] for j, e in enumerate(m.eqs)])
        if S is None or S.shape != (m.N, nb):
            raise StataError(503, f"{m.prog} did not return a {m.N} x {nb} matrix of scores")   # VERIFY
        return S

    def fun(b, todo):
        if todo == 0:
            ll = _total(m, obs(b))
            return (ll if np.isfinite(ll) else -np.inf), None, None
        lnf = obs(b)
        ll = _total(m, lnf)
        if not np.isfinite(ll):
            return -np.inf, np.zeros(nb), np.full((nb, nb), np.nan)
        S = scores(b)
        g = (S * w[:, None]).sum(axis=0)
        H = None
        if m.method in ("lf2", "gf2"):
            H = _eval_l(s, m, b, 2)[3]
            if H is None:
                raise StataError(198, f"{m.prog} did not return the Hessian")    # VERIFY
        if H is None:
            # numerical Hessian from the gradient
            H = np.zeros((nb, nb))
            for j in range(nb):
                hj = 1e-4 * max(1.0, abs(b[j]))
                e = np.zeros(nb)
                e[j] = hj
                H[:, j] = ((scores(b + e) - scores(b - e)) * w[:, None]).sum(axis=0) / (2 * hj)
            H = (H + H.T) / 2
        m._S = S * w[:, None]
        return ll, g, H
    return fun


def _fun(s: "Session", m: Model):
    _setup(s, m)
    if m.method in _LFAMILY:
        return _l_fun(s, m)
    return _lf_fun(s, m) if m.method in ("lf", "mlexp") else _d_fun(s, m)
