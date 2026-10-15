"""Mata's moptimize() ([M-5] moptimize()).

Model: dependent variables, equations (independent variables, constant,
offset), weights and an evaluator of type lf, d0/d1/d2, lf0/lf1/lf2 or
gf0/gf1/gf2. The coefficient vector is the concatenation of the equations
(their variables, then _cons). Maximization uses the same Newton-Raphson
steps as ml (estimation/ml.py); missing derivatives are numerical.
Variance: oim, opg or robust (with optional clusters).

The iteration log and the result display follow ml's layout (VERIFY the
exact text that moptimize prints in Stata 14: compat/do/0408_mata_moptimize.do).
"""

from __future__ import annotations

import numpy as np

from ..core.errors import StataError
from ..core.formats import format_value
from .advanced import _box, _num_grad, _num_hess, _unbox
from .library import lib
from .values import MV, SYS, MataError, empty, real, string


class MoptState:
    def __init__(self):
        self.evaluator: str | None = None
        self.evaltype = "d0"
        self.depvars: dict[int, MV] = {}
        self.depnames: dict[int, list[str]] = {}
        self.eqs: dict[int, dict] = {}
        self.weight: np.ndarray | None = None
        self.wtype = ""
        self.which = "max"
        self.technique = "nr"
        self.trace = "value"
        self.maxiter = 16000
        self.ptol, self.vtol, self.nrtol = 1e-6, 1e-7, 1e-5
        self.vce = "oim"
        self.cluster: np.ndarray | None = None
        self.userinfo: dict[int, MV] = {}
        self.touse: np.ndarray | None = None
        self.valueid = "f(p)"
        self.search = "off"
        self.result: dict = {}
        self.b: np.ndarray | None = None
        self.N = 0
        # shifts added to the linear indexes while the derivatives of an lf evaluator are taken
        # with respect to them (see _lf_derivs)
        self.theta_shift: dict[int, float] = {}

    # layout ----------------------------------------------------------------
    def neq(self) -> int:
        return max(self.eqs) if self.eqs else 0

    def eq(self, i: int) -> dict:
        return self.eqs.setdefault(i, {"X": None, "xnames": [], "cons": True, "name": "",
                                       "offset": None, "coefs": None})

    def nobs(self) -> int:
        for e in self.eqs.values():
            if e["X"] is not None:
                return e["X"].shape[0]
        for d in self.depvars.values():
            return d.rows
        return 0

    def columns(self, i: int) -> int:
        e = self.eq(i)
        return (0 if e["X"] is None else e["X"].shape[1]) + (1 if e["cons"] else 0)

    def indices(self, i: int) -> tuple[int, int]:
        start = sum(self.columns(j) for j in range(1, i))
        return start, start + self.columns(i)

    def k(self) -> int:
        return sum(self.columns(j) for j in range(1, self.neq() + 1))

    def design(self, i: int) -> np.ndarray:
        e = self.eq(i)
        n = self.N
        parts = []
        if e["X"] is not None:
            parts.append(e["X"])
        if e["cons"]:
            parts.append(np.ones((n, 1)))
        return np.column_stack(parts) if parts else np.zeros((n, 0))

    def xb(self, i: int, b: np.ndarray) -> np.ndarray:
        a, z = self.indices(i)
        v = self.design(i) @ b[a:z] + self.theta_shift.get(i, 0.0)
        off = self.eq(i)["offset"]
        return v + off if off is not None else v

    def w(self) -> np.ndarray:
        return self.weight if self.weight is not None else np.ones(self.N)

    def names(self) -> tuple[list[str], list[str]]:
        eqn, coln = [], []
        for i in range(1, self.neq() + 1):
            e = self.eq(i)
            nm = e["name"] or f"eq{i}"
            p = 0 if e["X"] is None else e["X"].shape[1]
            xn = e["xnames"] if len(e["xnames"]) == p else [f"x{j + 1}" for j in range(p)]
            for c in xn + (["_cons"] if e["cons"] else []):
                eqn.append(nm)
                coln.append(c)
        return eqn, coln


def _st(v: MV) -> MoptState:
    return _unbox(v, MoptState)


def _data_matrix(eng, spec: MV) -> tuple[np.ndarray, list[str]]:
    """A real matrix, or a string of variable names read from the data (with touse)."""
    if spec.t == "string":
        names = eng.s.expand_varlist(spec.str_scalar())
        ds = eng.s.data
        cols = [ds.get(n).data.astype(np.float64) for n in names]
        X = np.column_stack(cols) if cols else np.zeros((ds.nobs, 0))
        X = np.where(X >= SYS, np.nan, X)
        return X, names
    return spec.a.astype(np.float64), []


@lib("moptimize_init", 0)
def _init(eng, v, r):
    return _box(MoptState())


@lib("moptimize_init_evaluator", 1, 2)
def _evaluator(eng, v, r):
    S = _st(v[0])
    if len(v) == 1:
        return string(S.evaluator or "")
    p = v[1]
    target = p.scalar().get() if p.t == "pointer" and p.scalar() is not None else None
    if not isinstance(target, str):
        raise MataError(3000, "function pointer required")   # VERIFY
    S.evaluator = target
    return None


def _setter(attr, conv):
    def f(eng, v, r):
        S = _st(v[0])
        if len(v) == 1:
            val = getattr(S, attr)
            return string(val) if isinstance(val, str) else real(float(val))
        setattr(S, attr, conv(v[1]))
        return None
    return f


lib("moptimize_init_evaluatortype", 1, 2)(_setter("evaltype", lambda x: x.str_scalar().strip()))
lib("moptimize_init_which", 1, 2)(_setter("which", lambda x: x.str_scalar().strip()))
lib("moptimize_init_technique", 1, 2)(_setter("technique", lambda x: x.str_scalar().split()[0]))
lib("moptimize_init_tracelevel", 1, 2)(_setter("trace", lambda x: x.str_scalar().strip()))
lib("moptimize_init_conv_maxiter", 1, 2)(_setter("maxiter", lambda x: x.int_scalar()))
lib("moptimize_init_conv_ptol", 1, 2)(_setter("ptol", lambda x: x.real_scalar()))
lib("moptimize_init_conv_vtol", 1, 2)(_setter("vtol", lambda x: x.real_scalar()))
lib("moptimize_init_conv_nrtol", 1, 2)(_setter("nrtol", lambda x: x.real_scalar()))
lib("moptimize_init_vcetype", 1, 2)(_setter("vce", lambda x: x.str_scalar().strip().lower()))
lib("moptimize_init_valueid", 1, 2)(_setter("valueid", lambda x: x.str_scalar()))
lib("moptimize_init_search", 1, 2)(_setter("search", lambda x: x.str_scalar()))
lib("moptimize_init_weighttype", 1, 2)(_setter("wtype", lambda x: x.str_scalar().strip()))


@lib("moptimize_init_depvar", 2, 3)
def _depvar(eng, v, r):
    S = _st(v[0])
    j = v[1].int_scalar()
    if len(v) == 2:
        return S.depvars.get(j, empty("real", 0, 0))
    X, names = _data_matrix(eng, v[2])
    S.depvars[j] = real(X)
    S.depnames[j] = names
    return None


@lib("moptimize_init_eq_indepvars", 2, 3)
def _indepvars(eng, v, r):
    S = _st(v[0])
    e = S.eq(v[1].int_scalar())
    if len(v) == 2:
        return real(e["X"]) if e["X"] is not None else empty("real", 0, 0)
    X, names = _data_matrix(eng, v[2])
    e["X"] = X if X.size or X.shape[1] else None
    e["xnames"] = names
    return None


@lib("moptimize_init_eq_cons", 2, 3)
def _cons(eng, v, r):
    S = _st(v[0])
    e = S.eq(v[1].int_scalar())
    if len(v) == 2:
        return string("on" if e["cons"] else "off")
    e["cons"] = v[2].str_scalar().strip().lower() != "off"
    return None


@lib("moptimize_init_eq_name", 2, 3)
def _eqname(eng, v, r):
    S = _st(v[0])
    e = S.eq(v[1].int_scalar())
    if len(v) == 2:
        return string(e["name"])
    e["name"] = v[2].str_scalar()
    return None


@lib("moptimize_init_eq_offset", 2, 3)
def _offset(eng, v, r):
    S = _st(v[0])
    e = S.eq(v[1].int_scalar())
    X, _ = _data_matrix(eng, v[2])
    e["offset"] = X.ravel()
    return None


@lib("moptimize_init_eq_exposure", 3)
def _exposure(eng, v, r):
    S = _st(v[0])
    e = S.eq(v[1].int_scalar())
    X, _ = _data_matrix(eng, v[2])
    e["offset"] = np.log(X.ravel())
    return None


@lib("moptimize_init_eq_coefs", 2, 3)
def _eqcoefs(eng, v, r):
    S = _st(v[0])
    e = S.eq(v[1].int_scalar())
    if len(v) == 2:
        return real(e["coefs"].reshape(1, -1)) if e["coefs"] is not None else empty("real", 0, 0)
    e["coefs"] = v[2].a.astype(np.float64).ravel()
    return None


@lib("moptimize_init_weight", 2)
def _weight(eng, v, r):
    S = _st(v[0])
    X, _ = _data_matrix(eng, v[1])
    S.weight = X.ravel()
    if not S.wtype:
        S.wtype = "iweight"
    return None


@lib("moptimize_init_cluster", 2)
def _cluster(eng, v, r):
    S = _st(v[0])
    X, _ = _data_matrix(eng, v[1])
    S.cluster = X.ravel()
    return None


@lib("moptimize_init_touse", 2)
def _touse(eng, v, r):
    S = _st(v[0])
    name = v[1].str_scalar()
    S.touse = eng.s.data.get(name).data != 0
    return None


@lib("moptimize_init_userinfo", 3)
def _userinfo(eng, v, r):
    S = _st(v[0])
    S.userinfo[v[1].int_scalar()] = v[2].copy()
    return None


@lib("moptimize_util_userinfo", 2)
def _util_userinfo(eng, v, r):
    S = _st(v[0])
    return S.userinfo.get(v[1].int_scalar(), empty("real", 0, 0))


@lib("moptimize_init_nuserinfo", 1, 2)
def _nuserinfo(eng, v, r):
    S = _st(v[0])
    return real(float(len(S.userinfo))) if len(v) == 1 else None


@lib("moptimize_init_ndepvars", 1, 2)
def _ndepvars(eng, v, r):
    S = _st(v[0])
    return real(float(len(S.depvars))) if len(v) == 1 else None


# ---------------------------------------------------------------------------
# utilities for evaluators
# ---------------------------------------------------------------------------

def _b(v: MV) -> np.ndarray:
    return v.a.astype(np.float64).ravel()


@lib("moptimize_util_xb", 3)
def _util_xb(eng, v, r):
    S = _st(v[0])
    i = v[2].int_scalar()
    e = S.eq(i)
    if (e["X"] is None or e["X"].shape[1] == 0) and e["offset"] is None:
        # a constant-only equation gives a scalar ([M-5] moptimize(); compat 0408)
        a, z = S.indices(i)
        b = _b(v[1])
        if z - a == 1:
            return real(float(b[a]) + S.theta_shift.get(i, 0.0))
    return real(S.xb(i, _b(v[1])).reshape(-1, 1))


@lib("moptimize_util_depvar", 2)
def _util_depvar(eng, v, r):
    S = _st(v[0])
    j = v[1].int_scalar()
    if j not in S.depvars:
        raise MataError(3498, f"dependent variable {j} not defined")   # VERIFY
    return S.depvars[j]


@lib("moptimize_util_sum", 2)
def _util_sum(eng, v, r):
    S = _st(v[0])
    x = v[1].a.astype(np.float64).ravel()
    w = S.w()
    if x.size == 1:
        return real(float(x[0] * w.sum()))
    if np.any(x >= SYS) or np.any(~np.isfinite(x)):
        return real(SYS)
    return real(float(np.sum(w * x)))


@lib("moptimize_util_vecsum", 4, 4, outs=(3,))
def _util_vecsum(eng, v, r):
    S = _st(v[0])
    i = v[1].int_scalar()
    s = v[2].a.astype(np.float64).ravel()
    X = S.design(i)
    if s.size == 1:
        s = np.full(S.N, s[0])
    if np.any(s >= SYS) and r[3] is not None:
        r[3].set(real(SYS))
    return real((X.T @ (S.w() * s)).reshape(1, -1))


@lib("moptimize_util_matsum", 5, 5, outs=(4,))
def _util_matsum(eng, v, r):
    S = _st(v[0])
    i, j = v[1].int_scalar(), v[2].int_scalar()
    d = v[3].a.astype(np.float64).ravel()
    if d.size == 1:
        d = np.full(S.N, d[0])
    if np.any(d >= SYS) and r[4] is not None:
        r[4].set(real(SYS))
    Xi, Xj = S.design(i), S.design(j)
    return real((Xi * (S.w() * d)[:, None]).T @ Xj)


@lib("moptimize_util_eq_indices", 2, 3)
def _util_eq_indices(eng, v, r):
    S = _st(v[0])
    a, z = S.indices(v[1].int_scalar())
    return real(np.array([[a + 1, z]], dtype=float))


# ---------------------------------------------------------------------------
# the optimization
# ---------------------------------------------------------------------------

def _call(eng, S: MoptState, M: MV, args: list, outs: int):
    f = eng.funcs.get(S.evaluator or "")
    if f is None:
        raise MataError(3499, f"{S.evaluator}() not found")
    vals = eng.call_with_values(f, [M] + args)
    return vals[-outs:]


def _evaluate(eng, S: MoptState, M: MV, b: np.ndarray, todo: int):
    """(ll, g, H, obs-level scores wrt parameters or None)."""
    et = S.evaltype.lower().replace("debug", "")
    brow = real(b.reshape(1, -1))
    w = S.w()
    if et == "lf":
        (lnfj,) = _call(eng, S, M, [brow, real(SYS)], 1)
        x = lnfj.a.astype(np.float64).ravel() if lnfj is not None else np.array([SYS])
        if x.size != S.N or np.any(x >= SYS) or np.any(~np.isfinite(x)):
            return -np.inf, None, None, None
        return float(np.sum(w * x)), None, None, None
    if et in ("d0", "d1", "d2"):
        lnf, g, H = _call(eng, S, M, [real(float(todo)), brow, real(SYS), real(SYS), real(SYS)], 3)
        val = float(lnf.a.ravel()[0]) if lnf is not None and lnf.a.size else SYS
        if val >= SYS or not np.isfinite(val):
            return -np.inf, None, None, None
        gg = g.a.astype(np.float64).ravel() if todo >= 1 and g is not None and g.a.size == b.size else None
        HH = H.a.astype(np.float64) if todo >= 2 and H is not None and H.a.shape == (b.size, b.size) else None
        return val, gg, HH, None
    # lf0/lf1/lf2 and gf0/gf1/gf2: f(M, todo, b, lnfj, S, H)
    lnfj, Sc, H = _call(eng, S, M, [real(float(todo)), brow, real(SYS), real(SYS), real(SYS)], 3)
    x = lnfj.a.astype(np.float64).ravel() if lnfj is not None else np.array([SYS])
    if np.any(x >= SYS) or np.any(~np.isfinite(x)):
        return -np.inf, None, None, None
    ll = float(np.sum(w * x)) if x.size == S.N else float(x.sum())
    scores = None
    gg = None
    if todo >= 1 and Sc is not None and Sc.a.size:
        sc = Sc.a.astype(np.float64)
        if et.startswith("lf"):
            # scores with respect to the linear indexes -> parameters
            cols = []
            for i in range(1, S.neq() + 1):
                cols.append(S.design(i) * sc[:, i - 1][:, None])
            scores = np.column_stack(cols)
        else:
            scores = sc
        gg = (w[:, None] * scores).sum(axis=0)
    HH = H.a.astype(np.float64) if todo >= 2 and H is not None and H.a.shape == (b.size, b.size) else None
    return ll, gg, HH, scores


def _lf_derivs(eng, S: MoptState, M: MV, b: np.ndarray):
    """Gradient and Hessian of an lf evaluator from derivatives of each observation's log
    likelihood with respect to the linear indexes, by central differences, as ml does for
    linear-form evaluators ([R] ml, Methods and formulas): g = sum_i X_i' dl/dtheta_i and
    H_ij = X_i' diag(d2l/dtheta_i dtheta_j) X_j. Much more precise than differences of the
    total log likelihood with respect to b (compat 0408)."""
    w = S.w()
    neq = S.neq()

    def lnfj(shift: dict) -> np.ndarray:
        S.theta_shift = shift
        try:
            (x,) = _call(eng, S, M, [real(b.reshape(1, -1)), real(SYS)], 1)
        finally:
            S.theta_shift = {}
        return x.a.astype(np.float64).ravel()

    l0 = lnfj({})
    hs = {}
    for i in range(1, neq + 1):
        th = S.xb(i, b)
        hs[i] = 1e-4 * (float(np.mean(np.abs(th))) + 1.0)
    d1, d2 = {}, {}
    plus, minus = {}, {}
    for i in range(1, neq + 1):
        h = hs[i]
        plus[i], minus[i] = lnfj({i: h}), lnfj({i: -h})
        d1[i] = (plus[i] - minus[i]) / (2 * h)
        d2[(i, i)] = (plus[i] - 2 * l0 + minus[i]) / h ** 2
    for i in range(1, neq + 1):
        for j in range(i + 1, neq + 1):
            hi, hj = hs[i], hs[j]
            v = (lnfj({i: hi, j: hj}) - lnfj({i: hi, j: -hj}) - lnfj({i: -hi, j: hj})
                 + lnfj({i: -hi, j: -hj})) / (4 * hi * hj)
            d2[(i, j)] = d2[(j, i)] = v
    k = b.size
    g = np.zeros(k)
    H = np.zeros((k, k))
    for i in range(1, neq + 1):
        ai, zi = S.indices(i)
        Xi = S.design(i)
        g[ai:zi] = Xi.T @ (w * d1[i])
        for j in range(1, neq + 1):
            aj, zj = S.indices(j)
            H[ai:zi, aj:zj] = (Xi * (w * d2[(i, j)])[:, None]).T @ S.design(j)
    scores = np.column_stack([S.design(i) * d1[i][:, None] for i in range(1, neq + 1)])
    return g, (H + H.T) / 2, scores


def _obs_scores(eng, S, M, b) -> np.ndarray | None:
    et = S.evaltype.lower().replace("debug", "")
    if et.startswith(("lf1", "lf2", "gf1", "gf2")):
        return _evaluate(eng, S, M, b, 1)[3]
    if et == "lf":
        return _lf_derivs(eng, S, M, b)[2]
    if et in ("lf0", "gf0"):
        # numerical observation-level scores
        k = b.size
        out = np.zeros((S.N, k))
        brow = b.copy()

        def lnfj(bb):
            if et == "lf":
                (x,) = _call(eng, S, M, [real(bb.reshape(1, -1)), real(SYS)], 1)
            else:
                x = _call(eng, S, M, [real(0.0), real(bb.reshape(1, -1)), real(SYS), real(SYS), real(SYS)], 3)[0]
            return x.a.astype(np.float64).ravel()
        for j in range(k):
            h = (abs(brow[j]) + 1e-3) * 1e-5
            bp, bm = brow.copy(), brow.copy()
            bp[j] += h
            bm[j] -= h
            out[:, j] = (lnfj(bp) - lnfj(bm)) / (2 * h)
        return out
    return None


@lib("moptimize", 1)
def _moptimize(eng, v, r):
    S = _st(v[0])
    M = v[0]
    _prepare(eng, S)
    if S.evaluator is None:
        raise MataError(3498, "evaluator not set")   # VERIFY
    sign = 1.0 if S.which == "max" else -1.0
    et = S.evaltype.lower().replace("debug", "")
    order = 2 if et in ("d2", "lf2", "gf2") else 1 if et in ("d1", "lf1", "gf1") else 0
    k = S.k()

    def fval(b):
        return sign * _evaluate(eng, S, M, b, 0)[0]

    def fun(b, todo):
        ll, g, H, _ = _evaluate(eng, S, M, b, min(todo, order))
        ll = sign * ll
        if todo == 0 or not np.isfinite(ll):
            return ll, None, None
        if et == "lf":
            g, H, _sc = _lf_derivs(eng, S, M, b)
            return ll, sign * g, sign * H
        g = sign * g if g is not None else _num_grad(fval, b, ll)
        if H is not None:
            H = sign * H
        else:
            def grad(bb):
                if order >= 1:
                    gg = _evaluate(eng, S, M, bb, 1)[1]
                    if gg is not None:
                        return sign * gg
                return _num_grad(fval, bb)
            H = _num_hess(fval, grad, b, analytic_grad=order >= 1)
        return ll, g, H

    b0 = np.zeros(k)
    for i in range(1, S.neq() + 1):
        c = S.eq(i)["coefs"]
        if c is not None:
            a, z = S.indices(i)
            if c.size != z - a:
                raise MataError(3200, "conformability error")
            b0[a:z] = c
    from ..estimation.ml import maximize
    try:
        res = maximize(fun, b0, maxiter=S.maxiter, tol=S.ptol, ltol=S.vtol, nrtol=S.nrtol)
    except StataError as e:
        raise MataError(e.rc, e.message)
    out = eng.s.output
    if S.trace != "none":
        for it, ll, note in res.log:
            out.write(f"{'Iteration ' + str(it) + ':':<14} {S.valueid} = "
                      f"{format_value(sign * ll, '%10.0g', pad=False).strip():>10}{note}\n", "text")
        if not res.converged:
            out.write("convergence not achieved\n", "error")
    b = res.b
    H = sign * res.H
    try:
        Vo = np.linalg.inv(-H)
        Vo = (Vo + Vo.T) / 2
    except np.linalg.LinAlgError:
        Vo = np.full((k, k), SYS)
    V = Vo
    vce = S.vce
    if vce in ("opg", "robust"):
        sc = _obs_scores(eng, S, M, b)
        if sc is None:
            raise MataError(3498, f"vce({vce}) requires an lf, lf# or gf# evaluator")   # VERIFY
        w = S.w()
        if vce == "opg":
            V = np.linalg.inv((sc * w[:, None]).T @ sc)
        else:
            ws = sc * w[:, None]
            if S.cluster is not None:
                uq, inv = np.unique(S.cluster, return_inverse=True)
                G = np.zeros((len(uq), k))
                np.add.at(G, inv, ws)
                meat = G.T @ G * len(uq) / (len(uq) - 1)
            else:
                n_eff = float(w.sum()) if S.wtype == "fweight" else float(S.N)
                meat = (ws.T @ (sc if S.wtype == "fweight" else ws)) * n_eff / (n_eff - 1)
            V = Vo @ meat @ Vo
        V = (V + V.T) / 2
    S.b = b
    S.result = {"coefs": b, "value": sign * res.ll, "gradient": sign * res.g, "H": H, "V": V,
                "iterations": res.iterations, "converged": res.converged, "N": S.N}
    return None


def _prepare(eng, S: MoptState) -> None:
    """Applies touse and drops observations with missing values in the model variables."""
    n = S.nobs()
    keep = np.ones(n, dtype=bool) if S.touse is None or len(S.touse) != n else S.touse.copy()
    for d in S.depvars.values():
        a = d.a.astype(np.float64)
        keep &= ~np.any((a >= SYS) | np.isnan(a), axis=1)
    for e in S.eqs.values():
        if e["X"] is not None:
            keep &= ~np.any(np.isnan(e["X"]) | (e["X"] >= SYS), axis=1)
    if not keep.all():
        for j, d in list(S.depvars.items()):
            S.depvars[j] = real(d.a[keep])
        for e in S.eqs.values():
            if e["X"] is not None:
                e["X"] = e["X"][keep]
            if e["offset"] is not None:
                e["offset"] = e["offset"][keep]
        if S.weight is not None:
            S.weight = S.weight[keep]
        if S.cluster is not None:
            S.cluster = S.cluster[keep]
        S.touse = None
    S.N = int(keep.sum())


def _res(key):
    def f(eng, v, r):
        S = _st(v[0])
        if key not in S.result:
            raise MataError(3498, "moptimize() has not been run")   # VERIFY
        x = S.result[key]
        if isinstance(x, bool):
            return real(1.0 if x else 0.0)
        if isinstance(x, (int, float)):
            return real(float(x))
        x = np.asarray(x, dtype=np.float64)
        return real(x.reshape(1, -1) if x.ndim == 1 else x)
    return f


lib("moptimize_result_coefs", 1)(_res("coefs"))
lib("moptimize_result_value", 1)(_res("value"))
lib("moptimize_result_gradient", 1)(_res("gradient"))
lib("moptimize_result_Hessian", 1)(_res("H"))
lib("moptimize_result_V", 1)(_res("V"))
lib("moptimize_result_iterations", 1)(_res("iterations"))
lib("moptimize_result_converged", 1)(_res("converged"))


@lib("moptimize_result_eq_coefs", 1, 2)
def _res_eq(eng, v, r):
    S = _st(v[0])
    i = v[1].int_scalar() if len(v) > 1 else 1
    a, z = S.indices(i)
    return real(S.result["coefs"][a:z].reshape(1, -1))


@lib("moptimize_result_colstripe", 1, 2)
def _colstripe(eng, v, r):
    S = _st(v[0])
    eqn, coln = S.names()
    return MV(np.array([[e, c] for e, c in zip(eqn, coln)], dtype=object).reshape(-1, 2), "string")


@lib("moptimize_result_post", 1, 2)
def _post(eng, v, r):
    """ereturn post b V: e(b), e(V), e(N), e(ll) and e(sample)-less results."""
    from ..commands.matrix import Matrix
    S = _st(v[0])
    eqn, coln = S.names()
    b = S.result["coefs"].reshape(1, -1)
    V = np.asarray(S.result["V"])
    mb = Matrix(b, ["y1"], list(coln))
    mV = Matrix(V, list(coln), list(coln))
    mb.roweq, mb.coleq = ["_"], list(eqn)
    mV.roweq, mV.coleq = list(eqn), list(eqn)
    eng.s.e = {"N": float(S.N), "ll": float(S.result["value"]), "k": float(len(coln)),
               "k_eq": float(S.neq()), "converged": 1.0 if S.result["converged"] else 0.0,
               "ic": float(S.result["iterations"]), "rank": float(np.linalg.matrix_rank(V)),
               "vce": S.vce, "opt": "moptimize", "b": mb, "V": mV}
    return None


@lib("moptimize_result_display", 1, 3)
def _display(eng, v, r):
    from ..estimation.results import CoefRow, coef_table
    S = _st(v[0])
    if "coefs" not in S.result:
        raise MataError(3498, "moptimize() has not been run")   # VERIFY
    b = S.result["coefs"]
    V = np.asarray(S.result["V"])
    eqn, coln = S.names()
    out = eng.s.output
    neq = S.neq()
    # Stata shows only the number of observations above the table (compat 0408)
    out.write("\n", "text")
    out.write(f"{'':<48}{'Number of obs':<18}= ", "text")
    out.write(f"{S.N:>10,}\n\n", "result")
    se = np.sqrt(np.maximum(np.diag(V), 0))
    trows = []
    for i in range(1, neq + 1):
        a, z = S.indices(i)
        if neq > 1:
            trows.append(CoefRow("eq", eqn[a]))
        for j in range(a, z):
            trows.append(CoefRow("coef", coln[j], float(b[j]), float(se[j]), name=coln[j]))
        if neq > 1 and i < neq:
            trows.append(CoefRow("sep"))
    level = float(eng.s.settings.get("level", "95"))
    vlab = {"robust": "Robust", "opg": "OPG"}.get(S.vce, "")
    # the first dependent variable heads the table, with one or several equations (compat 0408)
    dep = ""
    if S.depnames.get(1):
        dep = S.depnames[1][0]
    coef_table(eng.s, dep, trows, stat="z", df=None, level=level, vcetype=vlab)
    return None
