"""Arrays associativos (asarray) e otimização (optimize) do Mata.

optimize() segue [M-5] optimize(): avaliadores d0, d1 e d2 (função,
gradiente e hessiana), técnicas nr (Newton–Raphson) e bfgs, maximização ou
minimização, critérios de convergência ptol, vtol e nrtol. As derivadas
numéricas usam diferenças centrais; os passos do Stata não são
documentados, então os valores de cada iteração podem diferir nas últimas
casas (VERIFICAR o registro de iterações em compat/do/0406_mata_optimize.do).
"""

from __future__ import annotations

import numpy as np

from ..core.formats import format_value
from .library import lib
from .values import MV, MataError, SYS, empty, real, string


class AssocArray:
    def __init__(self, keytype: str, keydim: int):
        self.keytype = keytype
        self.keydim = keydim
        self.data: dict = {}
        self.notfound: MV = empty("real", 0, 0)


def _box(obj) -> MV:
    arr = np.empty((1, 1), dtype=object)
    arr[0, 0] = obj
    return MV(arr, "struct")


def _unbox(v: MV, cls):
    obj = v.a[0, 0] if v.a.size == 1 else None
    if not isinstance(obj, cls):
        raise MataError(3000, "invalid argument")   # VERIFICAR
    return obj


def _key(A: AssocArray, k: MV):
    vals = k.a.ravel().tolist()
    if len(vals) != A.keydim:
        raise MataError(3200, "conformability error")
    return tuple(str(x) if k.t == "string" else float(x) for x in vals)


@lib("asarray_create", 0, 3)
def _asarray_create(eng, v, r):
    keytype = v[0].str_scalar() if v and v[0] is not None else "string"
    keydim = v[1].int_scalar() if len(v) > 1 and v[1] is not None else 1
    return _box(AssocArray(keytype, keydim))


@lib("asarray", 2, 3)
def _asarray(eng, v, r):
    A = _unbox(v[0], AssocArray)
    key = _key(A, v[1])
    if len(v) == 3:
        A.data[key] = v[2].copy()
        return None
    return A.data.get(key, A.notfound).copy()


@lib("asarray_contains", 2)
def _asarray_contains(eng, v, r):
    A = _unbox(v[0], AssocArray)
    return real(1.0 if _key(A, v[1]) in A.data else 0.0)


@lib("asarray_remove", 2)
def _asarray_remove(eng, v, r):
    A = _unbox(v[0], AssocArray)
    A.data.pop(_key(A, v[1]), None)
    return None


@lib("asarray_elements", 1)
def _asarray_elements(eng, v, r):
    return real(float(len(_unbox(v[0], AssocArray).data)))


@lib("asarray_keys", 1)
def _asarray_keys(eng, v, r):
    A = _unbox(v[0], AssocArray)
    keys = sorted(A.data)          # VERIFICAR: o Stata devolve na ordem do hash
    if A.keytype == "string":
        out = np.empty((len(keys), A.keydim), dtype=object)
        for i, k in enumerate(keys):
            out[i, :] = list(k)
        return MV(out, "string") if keys else empty("string", 0, A.keydim)
    return real(np.array(keys, dtype=np.float64).reshape(len(keys), A.keydim))


@lib("asarray_notfound", 1, 2)
def _asarray_notfound(eng, v, r):
    A = _unbox(v[0], AssocArray)
    if len(v) == 1:
        return A.notfound.copy()
    A.notfound = v[1].copy()
    return None


# ---------------------------------------------------------------------------
# optimize()
# ---------------------------------------------------------------------------

class OptState:
    def __init__(self):
        self.evaluator: str | None = None
        self.evaltype = "d0"
        self.params: np.ndarray | None = None
        self.which = "max"
        self.technique = "nr"
        self.args: dict[int, MV] = {}
        self.trace = "value"
        self.maxiter = 300
        self.ptol, self.vtol, self.nrtol = 1e-6, 1e-7, 1e-5
        self.result: dict = {}


def _state(v: MV) -> OptState:
    return _unbox(v, OptState)


@lib("optimize_init", 0)
def _optimize_init(eng, v, r):
    return _box(OptState())


@lib("optimize_init_evaluator", 1, 2)
def _opt_evaluator(eng, v, r):
    S = _state(v[0])
    if len(v) == 1:
        return string(S.evaluator or "")
    p = v[1]
    if p.t != "pointer" or p.scalar() is None:
        raise MataError(3000, "function pointer required")   # VERIFICAR
    target = p.scalar().get()
    if not isinstance(target, str):
        raise MataError(3000, "function pointer required")
    S.evaluator = target
    return None


def _setter(attr, conv):
    def f(eng, v, r):
        S = _state(v[0])
        if len(v) == 1:
            val = getattr(S, attr)
            return string(val) if isinstance(val, str) else real(float(val))
        setattr(S, attr, conv(v[1]))
        return None
    return f


lib("optimize_init_evaluatortype", 1, 2)(_setter("evaltype", lambda x: x.str_scalar()))
lib("optimize_init_which", 1, 2)(_setter("which", lambda x: x.str_scalar()))
lib("optimize_init_technique", 1, 2)(_setter("technique", lambda x: x.str_scalar().split()[0]))
lib("optimize_init_tracelevel", 1, 2)(_setter("trace", lambda x: x.str_scalar()))
lib("optimize_init_conv_maxiter", 1, 2)(_setter("maxiter", lambda x: x.int_scalar()))
lib("optimize_init_conv_ptol", 1, 2)(_setter("ptol", lambda x: x.real_scalar()))
lib("optimize_init_conv_vtol", 1, 2)(_setter("vtol", lambda x: x.real_scalar()))
lib("optimize_init_conv_nrtol", 1, 2)(_setter("nrtol", lambda x: x.real_scalar()))


@lib("optimize_init_params", 1, 2)
def _opt_params(eng, v, r):
    S = _state(v[0])
    if len(v) == 1:
        return real(S.params.reshape(1, -1)) if S.params is not None else empty("real", 0, 0)
    S.params = v[1].a.astype(np.float64).ravel().copy()
    return None


@lib("optimize_init_argument", 3)
def _opt_argument(eng, v, r):
    S = _state(v[0])
    S.args[v[1].int_scalar()] = v[2].copy()
    return None


@lib("optimize_init_narguments", 1, 2)
def _opt_nargs(eng, v, r):
    S = _state(v[0])
    if len(v) == 1:
        return real(float(len(S.args)))
    return None


def _evaluate(eng, S: OptState, p: np.ndarray, todo: int):
    """Chama o avaliador: f(todo, p, args..., v, g, H) e devolve (v, g, H)."""
    f = eng.funcs.get(S.evaluator or "")
    if f is None:
        raise MataError(3499, f"{S.evaluator}() not found")
    extra = [S.args[k] for k in sorted(S.args)]
    values = [real(float(todo)), real(p.reshape(1, -1))] + extra + [real(SYS), real(SYS), real(SYS)]
    out = eng.call_with_values(f, values)
    n = len(values)
    vv, gg, HH = out[n - 3], out[n - 2], out[n - 1]
    val = float(vv.a.ravel()[0]) if vv is not None and vv.a.size else SYS
    g = gg.a.astype(np.float64).ravel() if todo >= 1 and gg is not None and gg.a.size == p.size else None
    H = HH.a.astype(np.float64) if todo >= 2 and HH is not None and HH.a.shape == (p.size, p.size) else None
    return val, g, H


def _num_grad(fun, p: np.ndarray) -> np.ndarray:
    g = np.zeros_like(p)
    for i in range(p.size):
        h = 1e-5 * (abs(p[i]) + 1)
        a, b = p.copy(), p.copy()
        a[i] += h
        b[i] -= h
        g[i] = (fun(a) - fun(b)) / (2 * h)
    return g


def _num_hess(fun, grad, p: np.ndarray) -> np.ndarray:
    """Hessiana por segundas diferenças centrais da própria função."""
    k = p.size
    H = np.zeros((k, k))
    hs = 1e-4 * (np.abs(p) + 1)
    f0 = fun(p)
    for i in range(k):
        for j in range(i, k):
            if i == j:
                a, b = p.copy(), p.copy()
                a[i] += hs[i]
                b[i] -= hs[i]
                H[i, i] = (fun(a) - 2 * f0 + fun(b)) / hs[i] ** 2
            else:
                pp, pm, mp, mm = (p.copy() for _ in range(4))
                pp[i] += hs[i]; pp[j] += hs[j]          # noqa: E702
                pm[i] += hs[i]; pm[j] -= hs[j]          # noqa: E702
                mp[i] -= hs[i]; mp[j] += hs[j]          # noqa: E702
                mm[i] -= hs[i]; mm[j] -= hs[j]          # noqa: E702
                H[i, j] = H[j, i] = (fun(pp) - fun(pm) - fun(mp) + fun(mm)) / (4 * hs[i] * hs[j])
    return H


@lib("optimize", 1)
def _optimize(eng, v, r):
    S = _state(v[0])
    if S.params is None:
        raise MataError(111, "initial values not set")   # VERIFICAR
    sign = 1.0 if S.which == "max" else -1.0
    order = {"d0": 0, "d1": 1, "d2": 2}.get(S.evaltype.lower().replace("debug", ""), 0)

    def fval(p):
        return sign * _evaluate(eng, S, p, 0)[0]

    def grad(p):
        if order >= 1:
            _, g, _ = _evaluate(eng, S, p, 1)
            if g is not None:
                return sign * g
        return _num_grad(fval, p)

    def hess(p):
        if order >= 2:
            _, _, H = _evaluate(eng, S, p, 2)
            if H is not None:
                return sign * H
        return _num_hess(fval, grad, p)

    out = eng.s.output
    p = S.params.astype(np.float64).copy()
    f = fval(p)
    if f >= SYS or not np.isfinite(f):
        raise MataError(3498, "could not calculate numerical derivatives -- missing values encountered")
    converged = False
    it = 0
    Binv = np.eye(p.size)
    g = grad(p)
    while True:
        concave = True
        if S.technique == "bfgs":
            step = Binv @ g
        else:
            H = hess(p)
            try:
                np.linalg.cholesky(-H)
            except np.linalg.LinAlgError:
                concave = False
            Hm = H.copy()
            if not concave:
                # força definida negativa somando um múltiplo da identidade
                lam = max(1e-8, float(np.max(np.linalg.eigvalsh(H))) * 1.1 + 1e-6)
                Hm = H - lam * np.eye(p.size)
            step = -np.linalg.solve(Hm, g)
        if S.trace != "none":
            note = "  (not concave)" if not concave else ""
            out.write(f"Iteration {it}:   f(p) = {format_value(sign * f, '%10.0g', pad=False).strip():>10}{note}\n",
                      "text")
        if it >= S.maxiter:
            break
        # passo com recuo (step halving)
        t = 1.0
        for _ in range(60):
            pn = p + t * step
            fn = fval(pn)
            if fn < SYS and np.isfinite(fn) and fn >= f - 1e-12 * max(1.0, abs(f)):
                break
            t /= 2
        else:
            break
        gn = grad(pn)
        if S.technique == "bfgs":
            s_ = pn - p
            y = -(gn - g)
            sy = float(s_ @ y)
            if sy > 1e-12:
                rho = 1 / sy
                Id = np.eye(p.size)
                Binv = (Id - rho * np.outer(s_, y)) @ Binv @ (Id - rho * np.outer(y, s_)) + rho * np.outer(s_, s_)
        dp = np.max(np.abs(pn - p) / (np.abs(p) + 1))
        dv = abs(fn - f) / (abs(f) + 1)
        p, f, g = pn, fn, gn
        it += 1
        try:
            Hc = hess(p) if S.technique != "bfgs" else -np.linalg.inv(Binv)
            nr = float(abs(g @ np.linalg.solve(Hc, g)))
        except np.linalg.LinAlgError:
            nr = np.inf
        if (dp < S.ptol or dv < S.vtol) and nr < S.nrtol:
            converged = True
            if S.trace != "none":
                out.write(f"Iteration {it}:   f(p) = {format_value(sign * f, '%10.0g', pad=False).strip():>10}\n",
                          "text")
            break
    H = hess(p)
    try:
        V = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        V = np.full(H.shape, SYS)
    S.params = p
    S.result = {"params": p, "value": sign * f, "gradient": sign * g, "H": sign * H, "V": V,
                "iterations": it, "converged": converged}
    return real(p.reshape(1, -1))


def _result(key):
    def f(eng, v, r):
        S = _state(v[0])
        if key not in S.result:
            raise MataError(3498, "optimize() has not been run")   # VERIFICAR
        x = S.result[key]
        if isinstance(x, bool):
            return real(1.0 if x else 0.0)
        if isinstance(x, (int, float)):
            return real(float(x))
        x = np.asarray(x, dtype=np.float64)
        return real(x.reshape(1, -1) if x.ndim == 1 else x)
    return f


lib("optimize_result_params", 1)(_result("params"))
lib("optimize_result_value", 1)(_result("value"))
lib("optimize_result_gradient", 1)(_result("gradient"))
lib("optimize_result_Hessian", 1)(_result("H"))
lib("optimize_result_V", 1)(_result("V"))
lib("optimize_result_iterations", 1)(_result("iterations"))
lib("optimize_result_converged", 1)(_result("converged"))
