"""Biblioteca de funções do Mata ([M-5]).

Cada função recebe (engine, valores, referências) e devolve um MV ou None
(funções void). `refs[k]` é uma referência à variável passada como k-ésimo
argumento, usada pelas funções que gravam resultados nos argumentos
(symeigensystem, svd, _sort...).

Funções elementares com o mesmo nome e comportamento que as do Stata
(sqrt, ln, normal, substr, strupper...) reaproveitam lang/functions.py,
aplicadas elemento a elemento.
"""

from __future__ import annotations

import math
import re

import numpy as np

from ..core import missing as M
from ..core.formats import format_value
from ..lang import functions as F
from .values import MV, MataError, SYS, clean, conformability, empty, real, string, type_mismatch

LIBRARY: dict[str, tuple] = {}


def lib(name: str, lo: int, hi: int | None = None, outs: tuple = ()):
    def deco(fn):
        LIBRARY[name] = (fn, lo, lo if hi is None else hi, set(outs))
        return fn
    return deco


# ---------------------------------------------------------------------------
# auxiliares
# ---------------------------------------------------------------------------

def _r(v: MV) -> np.ndarray:
    if v is None:
        raise MataError(3499, "argument missing")
    if v.t != "real":
        raise MataError(3253 if v.t != "string" else 3253, "nonreal found where real required")
    return v.a


def _s(v: MV) -> np.ndarray:
    if v.t != "string":
        raise MataError(3254, f"{v.t} found where string required")
    return v.a


def _int(v: MV) -> int:
    return v.int_scalar()


def _nan(a: np.ndarray) -> np.ndarray:
    return np.where(a >= SYS, np.nan, a)


def _miss(a: np.ndarray) -> np.ndarray:
    return clean(a)


def _broadcast(vals: list[MV]) -> tuple[list[np.ndarray], tuple]:
    arrays = [v.a for v in vals]
    for x in arrays:
        for y in arrays:
            for dx, dy in zip(x.shape, y.shape):
                if dx != dy and dx != 1 and dy != 1:
                    raise conformability()
    shape = np.broadcast_shapes(*[a.shape for a in arrays])
    return [np.broadcast_to(a, shape) for a in arrays], shape


def _stata_elementwise(fname: str):
    """Função do Stata aplicada elemento a elemento (argumentos c-conformes)."""
    def fn(eng, vals, refs):
        arrays, shape = _broadcast(vals)
        out = np.empty(shape, dtype=object)
        for idx in np.ndindex(shape):
            args = []
            for a, v in zip(arrays, vals):
                x = a[idx]
                args.append(str(x) if v.t == "string" else float(x))
            try:
                out[idx] = F.call(fname, args)
            except Exception as e:     # noqa: BLE001
                from ..core.errors import StataError
                if isinstance(e, StataError) and e.rc == 109:
                    raise type_mismatch()
                raise
        if out.size and isinstance(out.flat[0], str):
            return MV(out, "string")
        if not out.size:
            return MV(np.zeros(shape), "real")
        return real(out.astype(np.float64))
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
    lib(_name, _nargs[0], _nargs[1])(_stata_elementwise(_name))


@lib("factorial", 1)
def _factorial(eng, v, r):
    a = _r(v[0])
    with np.errstate(all="ignore"):
        out = np.where((a >= 0) & (a == np.trunc(a)) & (a < SYS),
                       np.vectorize(lambda x: math.gamma(x + 1) if x < 171 else np.inf)(np.where(a < SYS, a, 0)),
                       np.nan)
    return real(clean(out))


# ---------------------------------------------------------------------------
# tamanho e tipo
# ---------------------------------------------------------------------------

@lib("rows", 1)
def _rows(eng, v, r):
    return real(v[0].rows)


@lib("cols", 1)
def _cols(eng, v, r):
    return real(v[0].cols)


@lib("length", 1)
def _length(eng, v, r):
    return real(v[0].rows * v[0].cols)


@lib("eltype", 1)
def _eltype(eng, v, r):
    return string(v[0].t)


@lib("orgtype", 1)
def _orgtype(eng, v, r):
    return string(v[0].orgtype())


for _t in ("real", "string", "complex", "pointer"):
    lib("is" + _t, 1)(lambda eng, v, r, _t=_t: real(1.0 if v[0].t == _t else 0.0))


@lib("isrealvalues", 1)
def _isrealvalues(eng, v, r):
    if v[0].t == "real":
        return real(1.0)
    if v[0].t == "complex":
        return real(1.0 if np.all(v[0].a.imag == 0) else 0.0)
    return real(0.0)


@lib("isview", 1)
def _isview(eng, v, r):
    return real(0.0)


@lib("issymmetric", 1)
def _issym(eng, v, r):
    a = v[0].a
    return real(1.0 if a.shape[0] == a.shape[1] and np.array_equal(a, a.T) else 0.0)


@lib("isdiagonal", 1)
def _isdiag(eng, v, r):
    a = v[0].a
    return real(1.0 if a.shape[0] == a.shape[1] and np.array_equal(a, np.diag(np.diag(a))) else 0.0)


@lib("sizeof", 1)
def _sizeof(eng, v, r):
    return real(float(v[0].a.nbytes))


# ---------------------------------------------------------------------------
# construção
# ---------------------------------------------------------------------------

@lib("J", 3)
def _J(eng, v, r):
    nr, nc, val = _int(v[0]), _int(v[1]), v[2]
    if nr < 0 or nc < 0:
        raise MataError(3300, "argument out of range")
    if val.is_scalar:
        return MV(np.full((nr, nc), val.a[0, 0], dtype=val.a.dtype), val.t)
    return MV(np.tile(val.a, (nr, nc)), val.t)


@lib("I", 1, 2)
def _I(eng, v, r):
    n = _int(v[0])
    m = _int(v[1]) if len(v) > 1 else n
    return real(np.eye(n, m))


@lib("e", 2)
def _e(eng, v, r):
    i, n = _int(v[0]), _int(v[1])
    if not 1 <= i <= n:
        raise MataError(3300, "argument out of range")
    out = np.zeros((1, n))
    out[0, i - 1] = 1
    return real(out)


@lib("range", 3)
def _range(eng, v, r):
    a, b, d = v[0].real_scalar(), v[1].real_scalar(), abs(v[2].real_scalar())
    if d == 0:
        raise MataError(3300, "argument out of range")
    n = int(math.floor(abs(b - a) / d + 1e-12)) + 1
    step = d if b >= a else -d
    return real((a + step * np.arange(n)).reshape(-1, 1))


@lib("rangen", 3)
def _rangen(eng, v, r):
    a, b, n = v[0].real_scalar(), v[1].real_scalar(), _int(v[2])
    if n <= 0:
        return empty("real", 0, 1)
    if n == 1:
        return real(np.array([[a]]))
    return real(np.linspace(a, b, n).reshape(-1, 1))


@lib("runiform", 2, 4)
def _runiform(eng, v, r):
    from ..core import rng
    nr, nc = _int(v[0]), _int(v[1])
    vals = rng.draw("runiform", [x.real_scalar() for x in v[2:]], nr * nc)
    return real(vals.reshape(nr, nc))


@lib("rnormal", 4)
def _rnormal(eng, v, r):
    from ..core import rng
    nr, nc = _int(v[0]), _int(v[1])
    m, s = v[2], v[3]
    shape = (nr, nc)
    mm = np.broadcast_to(_r(m), shape) if m.is_scalar or m.a.shape == shape else None
    ss = np.broadcast_to(_r(s), shape) if s.is_scalar or s.a.shape == shape else None
    if mm is None or ss is None:
        raise conformability()
    vals = rng.draw("rnormal", [mm.ravel(), ss.ravel()], nr * nc)
    return real(vals.reshape(shape))


@lib("rseed", 0, 1)
def _rseed(eng, v, r):
    from ..core import rng
    if not v:
        return string(rng.RNG.state())
    x = v[0]
    if x.t == "string":
        rng.RNG.set_state(x.str_scalar())
    else:
        rng.RNG.seed(_int(x) % 2**64)
    return None


# ---------------------------------------------------------------------------
# somas, extremos, médias
# ---------------------------------------------------------------------------

def _zero_missing(a):
    return np.where(a >= SYS, 0.0, a)


@lib("sum", 1, 2)
def _sum(eng, v, r):
    a = _r(v[0])
    return real(float(_zero_missing(a).sum()))


@lib("quadsum", 1, 2)
def _quadsum(eng, v, r):
    return _sum(eng, v, r)


@lib("rowsum", 1, 2)
def _rowsum(eng, v, r):
    return real(_zero_missing(_r(v[0])).sum(axis=1, keepdims=True))


@lib("colsum", 1, 2)
def _colsum(eng, v, r):
    return real(_zero_missing(_r(v[0])).sum(axis=0, keepdims=True))


LIBRARY["quadrowsum"] = LIBRARY["rowsum"]
LIBRARY["quadcolsum"] = LIBRARY["colsum"]


def _ext(fn, axis):
    def f(eng, v, r):
        a = _r(v[0])
        if a.size == 0:
            return real(SYS) if axis is None else empty("real", 0 if axis == 1 else 1, 0 if axis == 0 else 1)
        x = _nan(a)
        import warnings
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            res = fn(x) if axis is None else fn(x, axis=axis, keepdims=True)
        return real(clean(res))
    return f


lib("max", 1)(_ext(np.nanmax, None))
lib("min", 1)(_ext(np.nanmin, None))
lib("rowmax", 1)(_ext(np.nanmax, 1))
lib("rowmin", 1)(_ext(np.nanmin, 1))
lib("colmax", 1)(_ext(np.nanmax, 0))
lib("colmin", 1)(_ext(np.nanmin, 0))


@lib("minmax", 1)
def _minmax(eng, v, r):
    lo = LIBRARY["min"][0](eng, v, r)
    hi = LIBRARY["max"][0](eng, v, r)
    return real(np.array([[lo.a[0, 0], hi.a[0, 0]]]))


@lib("colminmax", 1)
def _colminmax(eng, v, r):
    return real(np.vstack([LIBRARY["colmin"][0](eng, v, r).a, LIBRARY["colmax"][0](eng, v, r).a]))


@lib("rowminmax", 1)
def _rowminmax(eng, v, r):
    return real(np.hstack([LIBRARY["rowmin"][0](eng, v, r).a, LIBRARY["rowmax"][0](eng, v, r).a]))


def _complete_rows(*arrays) -> np.ndarray:
    ok = np.ones(arrays[0].shape[0], dtype=bool)
    for a in arrays:
        ok &= ~np.any(a >= SYS, axis=1)
    return ok


def _weights(v, n):
    if len(v) > 1 and v[1] is not None:
        w = _r(v[1]).reshape(-1, 1)
        if w.shape[0] == 1 and n != 1:
            w = np.full((n, 1), w[0, 0])
        return w
    return np.ones((n, 1))


@lib("mean", 1, 2)
def _mean(eng, v, r):
    X = _r(v[0])
    w = _weights(v, X.shape[0])
    ok = _complete_rows(X, w)
    Xo, wo = X[ok], w[ok]
    if Xo.shape[0] == 0 or wo.sum() == 0:
        return real(np.full((1, X.shape[1]), SYS))
    return real(((wo * Xo).sum(0) / wo.sum()).reshape(1, -1))


def _var(X, w, corr=False):
    ok = _complete_rows(X, w)
    Xo, wo = X[ok], w[ok]
    N = wo.sum()
    k = X.shape[1]
    if N <= 1:
        return np.full((k, k), SYS)
    m = (wo * Xo).sum(0) / N
    D = Xo - m
    V = (wo * D).T @ D / (N - 1)
    if corr:
        sd = np.sqrt(np.diag(V))
        with np.errstate(all="ignore"):
            V = V / np.outer(sd, sd)
    return clean(V)


@lib("variance", 1, 2)
def _variance(eng, v, r):
    X = _r(v[0])
    return real(_var(X, _weights(v, X.shape[0])))


@lib("quadvariance", 1, 2)
def _qvariance(eng, v, r):
    return _variance(eng, v, r)


@lib("correlation", 1, 2)
def _correlation(eng, v, r):
    X = _r(v[0])
    return real(_var(X, _weights(v, X.shape[0]), corr=True))


@lib("meanvariance", 1, 2)
def _meanvariance(eng, v, r):
    return real(np.vstack([_mean(eng, v, r).a, _variance(eng, v, r).a]))


@lib("cross", 2, 4)
def _cross(eng, v, r):
    """cross(X, Z), cross(X, w, Z): X'Z (ou X'diag(w)Z); linhas com missing
    em X ou Z ficam de fora."""
    if len(v) == 2:
        X, Z, w = _r(v[0]), _r(v[1]), None
    elif len(v) == 3:
        X, w, Z = _r(v[0]), _r(v[1]), _r(v[2])
    else:
        raise MataError(3001, "expected 2 or 3 arguments")
    if X.shape[0] != Z.shape[0]:
        raise conformability()
    ok = _complete_rows(X, Z)
    W = np.ones((X.shape[0], 1)) if w is None else (np.full((X.shape[0], 1), w[0, 0]) if w.size == 1
                                                     else w.reshape(-1, 1))
    ok &= ~np.any(W >= SYS, axis=1)
    return real((W[ok] * X[ok]).T @ Z[ok])


LIBRARY["quadcross"] = LIBRARY["cross"]


@lib("crossdev", 4, 5)
def _crossdev(eng, v, r):
    if len(v) == 4:
        X, x, Z, z = (_r(a) for a in v)
        return real((X - x).T @ (Z - z))
    X, x, w, Z, z = (_r(a) for a in v)
    return real(((X - x) * w.reshape(-1, 1)).T @ (Z - z))


@lib("runningsum", 1)
def _runningsum(eng, v, r):
    a = _r(v[0])
    flat = np.cumsum(_zero_missing(a.ravel()))
    return real(flat.reshape(a.shape))


LIBRARY["quadrunningsum"] = LIBRARY["runningsum"]


# ---------------------------------------------------------------------------
# missing
# ---------------------------------------------------------------------------

def _is_miss(v: MV) -> np.ndarray:
    if v.t == "string":
        return np.vectorize(lambda s: s == "", otypes=[bool])(v.a) if v.a.size else np.zeros(v.a.shape, bool)
    if v.t == "real":
        return v.a >= SYS
    return np.zeros(v.a.shape, dtype=bool)


@lib("missing", 1)
def _missing(eng, v, r):
    return real(float(_is_miss(v[0]).sum()))


@lib("nonmissing", 1)
def _nonmissing(eng, v, r):
    return real(float((~_is_miss(v[0])).sum()))


@lib("hasmissing", 1)
def _hasmissing(eng, v, r):
    return real(1.0 if _is_miss(v[0]).any() else 0.0)


@lib("rowmissing", 1)
def _rowmissing(eng, v, r):
    return real(_is_miss(v[0]).sum(1, keepdims=True).astype(float))


@lib("colmissing", 1)
def _colmissing(eng, v, r):
    return real(_is_miss(v[0]).sum(0, keepdims=True).astype(float))


@lib("rownonmissing", 1)
def _rownonmissing(eng, v, r):
    return real((~_is_miss(v[0])).sum(1, keepdims=True).astype(float))


@lib("colnonmissing", 1)
def _colnonmissing(eng, v, r):
    return real((~_is_miss(v[0])).sum(0, keepdims=True).astype(float))


@lib("editmissing", 2)
def _editmissing(eng, v, r):
    out = v[0].copy()
    m = _is_miss(out)
    out.a[m] = v[1].scalar()
    return out


@lib("_editmissing", 2, outs=(0,))
def _editmissing_(eng, v, r):
    res = _editmissing(eng, v, r)
    r[0].set(res)
    return None


@lib("editvalue", 3)
def _editvalue(eng, v, r):
    out = v[0].copy()
    out.a[out.a == v[1].scalar()] = v[2].scalar()
    return out


@lib("all", 1)
def _all(eng, v, r):
    return real(1.0 if np.all(_r(v[0]) != 0) else 0.0)


@lib("any", 1)
def _any(eng, v, r):
    return real(1.0 if np.any(_r(v[0]) != 0) else 0.0)


# ---------------------------------------------------------------------------
# forma e seleção
# ---------------------------------------------------------------------------

@lib("rowshape", 2)
def _rowshape(eng, v, r):
    a, n = v[0], _int(v[1])
    if n == 0 or a.a.size % n:
        raise conformability()
    return MV(a.a.reshape(n, -1).copy(), a.t)


@lib("colshape", 2)
def _colshape(eng, v, r):
    a, n = v[0], _int(v[1])
    if n == 0 or a.a.size % n:
        raise conformability()
    return MV(a.a.reshape(-1, n).copy(), a.t)


@lib("vec", 1)
def _vec(eng, v, r):
    return MV(v[0].a.T.reshape(-1, 1).copy(), v[0].t)


@lib("vech", 1)
def _vech(eng, v, r):
    a = v[0].a
    if a.shape[0] != a.shape[1]:
        raise conformability()
    out = [a[i, j] for j in range(a.shape[1]) for i in range(j, a.shape[0])]
    return MV(np.array(out, dtype=a.dtype).reshape(-1, 1), v[0].t)


@lib("invvech", 1)
def _invvech(eng, v, r):
    x = _r(v[0]).ravel()
    n = int((math.isqrt(8 * x.size + 1) - 1) // 2)
    if n * (n + 1) // 2 != x.size:
        raise conformability()
    out = np.zeros((n, n))
    k = 0
    for j in range(n):
        for i in range(j, n):
            out[i, j] = out[j, i] = x[k]
            k += 1
    return real(out)


@lib("select", 2)
def _select(eng, v, r):
    X, s = v[0], _r(v[1])
    if s.shape[1] == 1 and s.shape[0] == X.rows:
        keep = s[:, 0] != 0
        return MV(X.a[keep, :].copy(), X.t)
    if s.shape[0] == 1 and s.shape[1] == X.cols:
        keep = s[0, :] != 0
        return MV(X.a[:, keep].copy(), X.t)
    raise conformability()


@lib("selectindex", 1)
def _selectindex(eng, v, r):
    a = _r(v[0])
    idx = np.flatnonzero(a.ravel() != 0) + 1.0
    return real(idx.reshape(1, -1) if a.shape[0] == 1 else idx.reshape(-1, 1))


def _sort_order(X: MV, idx: np.ndarray) -> np.ndarray:
    keys = []
    for k in idx.ravel()[::-1]:
        k = int(k)
        col = X.a[:, abs(k) - 1]
        if X.t == "string":
            _, inv = np.unique(col.astype(str), return_inverse=True)
            key = inv.astype(np.float64)
        else:
            key = col.astype(np.float64)
        keys.append(-key if k < 0 else key)
    if not keys:
        return np.arange(X.rows)
    return np.lexsort(keys)


@lib("sort", 2)
def _sort(eng, v, r):
    X = v[0]
    order = _sort_order(X, _r(v[1]))
    return MV(X.a[order].copy(), X.t)


@lib("_sort", 2, outs=(0,))
def _sort_in_place(eng, v, r):
    res = _sort(eng, v, r)
    if r[0] is not None:
        r[0].set(res)
    return None


@lib("order", 2)
def _order(eng, v, r):
    return real((_sort_order(v[0], _r(v[1])) + 1.0).reshape(-1, 1))


@lib("revorder", 1)
def _revorder(eng, v, r):
    p = _r(v[0])
    return real(p[::-1].copy() if p.shape[1] == 1 else p[:, ::-1].copy())


@lib("uniqrows", 1)
def _uniqrows(eng, v, r):
    X = v[0]
    if X.rows == 0:
        return X.copy()
    if X.t == "string":
        rows = sorted({tuple(row) for row in X.a.tolist()})
        return MV(np.array(rows, dtype=object).reshape(len(rows), X.cols), "string")
    u = np.unique(X.a, axis=0)
    return real(u)


@lib("jumble", 1)
def _jumble(eng, v, r):
    from ..core import rng
    X = v[0]
    u = rng.RNG.uniform(X.rows)
    return MV(X.a[np.argsort(u, kind="stable")].copy(), X.t)


# ---------------------------------------------------------------------------
# álgebra linear
# ---------------------------------------------------------------------------

def _square(a: np.ndarray) -> None:
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise conformability()


@lib("diag", 1)
def _diag(eng, v, r):
    a = _r(v[0])
    if a.shape[0] == 1 or a.shape[1] == 1:
        return real(np.diag(a.ravel()))
    _square(a)
    return real(np.diag(np.diag(a)))


@lib("diagonal", 1)
def _diagonal(eng, v, r):
    return MV(np.diag(v[0].a).reshape(-1, 1).copy(), v[0].t)


@lib("trace", 1)
def _trace(eng, v, r):
    a = _r(v[0])
    _square(a)
    return real(float(_zero_missing(np.diag(a)).sum()) if a.size else 0.0)


@lib("det", 1)
def _det(eng, v, r):
    a = _r(v[0])
    _square(a)
    if np.any(a >= SYS):
        return real(SYS)
    return real(float(np.linalg.det(a)) if a.size else 1.0)


@lib("rank", 1, 2)
def _rank(eng, v, r):
    a = _r(v[0])
    if np.any(a >= SYS):
        return real(SYS)
    return real(float(np.linalg.matrix_rank(a)) if a.size else 0.0)


@lib("invsym", 1, 2)
def _invsym(eng, v, r):
    """Inversa generalizada de matriz simétrica: varre as colunas em ordem e
    zera as linhas/colunas que se mostram colineares (como o invsym() do
    Stata: em matriz singular, as colunas redundantes viram zero)."""
    a = _r(v[0]).astype(np.float64)
    _square(a)
    if np.any(a >= SYS):
        return real(np.full(a.shape, SYS))
    n = a.shape[0]
    A = (a + a.T) / 2
    keep: list[int] = []
    tol = 1e-13 * max(1.0, float(np.max(np.abs(np.diag(A)))) if n else 1.0)
    order = [int(x) - 1 for x in _r(v[1]).ravel()] if len(v) > 1 and v[1] is not None else []
    seq_ = order + [k for k in range(n) if k not in order]
    for k in seq_:
        trial = keep + [k]
        sub = A[np.ix_(trial, trial)]
        # pivô do sweep: variância residual da coluna k dado o que já entrou
        if keep:
            S = A[np.ix_(keep, keep)]
            b = A[np.ix_(keep, [k])]
            try:
                resid = float(A[k, k] - (b.T @ np.linalg.solve(S, b))[0, 0])
            except np.linalg.LinAlgError:
                resid = 0.0
        else:
            resid = float(A[k, k])
        if abs(resid) > tol * max(1.0, abs(A[k, k])):
            keep = trial
        del sub
    out = np.zeros((n, n))
    if keep:
        out[np.ix_(keep, keep)] = np.linalg.inv(A[np.ix_(keep, keep)])
    return real(out)


@lib("luinv", 1)
def _luinv(eng, v, r):
    a = _r(v[0])
    _square(a)
    try:
        if np.any(a >= SYS) or (a.size and abs(np.linalg.det(a)) < 1e-300):
            raise np.linalg.LinAlgError
        return real(np.linalg.inv(a) if a.size else a.copy())
    except np.linalg.LinAlgError:
        return real(np.full(a.shape, SYS))


LIBRARY["inv"] = LIBRARY["luinv"]


@lib("cholinv", 1)
def _cholinv(eng, v, r):
    a = _r(v[0])
    _square(a)
    try:
        np.linalg.cholesky(a)
        return real(np.linalg.inv(a))
    except np.linalg.LinAlgError:
        return real(np.full(a.shape, SYS))


@lib("pinv", 1)
def _pinv(eng, v, r):
    return real(np.linalg.pinv(_nan(_r(v[0]))))


@lib("cholesky", 1)
def _cholesky(eng, v, r):
    a = _r(v[0])
    _square(a)
    try:
        return real(np.linalg.cholesky(a))
    except np.linalg.LinAlgError:
        return real(np.full(a.shape, SYS))


def _solver(kind):
    def f(eng, v, r):
        A, B = _r(v[0]), _r(v[1])
        if A.shape[0] != B.shape[0]:
            raise conformability()
        try:
            if kind == "qr":
                x = np.linalg.lstsq(A, B, rcond=None)[0]
            else:
                if kind == "chol":
                    np.linalg.cholesky(A)
                x = np.linalg.solve(A, B)
            return real(clean(x))
        except np.linalg.LinAlgError:
            return real(np.full((A.shape[1], B.shape[1]), SYS))
    return f


lib("lusolve", 2, 3)(_solver("lu"))
lib("cholsolve", 2, 3)(_solver("chol"))
lib("qrsolve", 2, 3)(_solver("qr"))
lib("svsolve", 2, 3)(_solver("qr"))


@lib("symeigensystem", 3, 3, outs=(1, 2))
def _symeigensystem(eng, v, r):
    a = _r(v[0])
    w, X = np.linalg.eigh(a)
    order = np.argsort(-w, kind="stable")
    r[1].set(real(X[:, order])) if r[1] else None
    r[2].set(real(w[order].reshape(1, -1))) if r[2] else None
    return None


@lib("symeigenvalues", 1)
def _symeigenvalues(eng, v, r):
    w = np.linalg.eigvalsh(_r(v[0]))
    return real(np.sort(w)[::-1].reshape(1, -1))


@lib("eigensystem", 3, 5, outs=(1, 2))
def _eigensystem(eng, v, r):
    w, X = np.linalg.eig(_r(v[0]))
    order = np.argsort(-np.abs(w), kind="stable")
    if r[1]:
        r[1].set(MV(X[:, order], "complex"))
    if r[2]:
        r[2].set(MV(w[order].reshape(1, -1), "complex"))
    return None


@lib("eigenvalues", 1)
def _eigenvalues(eng, v, r):
    w = np.linalg.eigvals(_r(v[0]))
    return MV(w[np.argsort(-np.abs(w), kind="stable")].reshape(1, -1), "complex")


@lib("svd", 4, 4, outs=(1, 2, 3))
def _svd(eng, v, r):
    U, s, Vt = np.linalg.svd(_r(v[0]), full_matrices=False)
    for k, val in ((1, real(U)), (2, real(s.reshape(-1, 1))), (3, real(Vt))):
        if r[k]:
            r[k].set(val)
    return None


@lib("svdsv", 1)
def _svdsv(eng, v, r):
    return real(np.linalg.svd(_r(v[0]), compute_uv=False).reshape(-1, 1))


@lib("norm", 1, 2)
def _norm(eng, v, r):
    a = _r(v[0])
    p = v[1].real_scalar() if len(v) > 1 else 2
    if a.shape[0] == 1 or a.shape[1] == 1:
        x = a.ravel()
        return real(float(np.max(np.abs(x))) if p >= SYS else float(np.sum(np.abs(x) ** p) ** (1 / p)))
    if p == 2:
        return real(float(np.linalg.norm(a, 2)))
    return real(float(np.linalg.norm(a, 1 if p == 1 else np.inf if p >= SYS else "fro")))


@lib("lowertriangle", 1, 2)
def _lowertriangle(eng, v, r):
    a = _r(v[0]).copy()
    out = np.tril(a)
    if len(v) > 1 and v[1] is not None:
        np.fill_diagonal(out, v[1].real_scalar())
    return real(out)


@lib("uppertriangle", 1, 2)
def _uppertriangle(eng, v, r):
    out = np.triu(_r(v[0]).copy())
    if len(v) > 1 and v[1] is not None:
        np.fill_diagonal(out, v[1].real_scalar())
    return real(out)


@lib("makesymmetric", 1)
def _makesymmetric(eng, v, r):
    a = _r(v[0])
    _square(a)
    return real(np.tril(a) + np.tril(a, -1).T)


@lib("_makesymmetric", 1, outs=(0,))
def _makesymmetric_(eng, v, r):
    if r[0]:
        r[0].set(_makesymmetric(eng, v, r))
    return None


@lib("Re", 1)
def _Re(eng, v, r):
    return real(np.real(v[0].a).astype(np.float64))


@lib("Im", 1)
def _Im(eng, v, r):
    return real(np.imag(v[0].a).astype(np.float64) if v[0].t == "complex" else np.zeros(v[0].a.shape))


@lib("C", 1, 2)
def _C(eng, v, r):
    if len(v) == 2:
        return MV(_r(v[0]) + 1j * _r(v[1]), "complex")
    return MV(v[0].a.astype(np.complex128), "complex")


# ---------------------------------------------------------------------------
# strings
# ---------------------------------------------------------------------------

@lib("strofreal", 1, 2)
def _strofreal(eng, v, r):
    fmt = v[1].str_scalar() if len(v) > 1 else "%9.0g"
    a = _r(v[0])
    out = np.empty(a.shape, dtype=object)
    for idx in np.ndindex(a.shape):
        x = float(a[idx])
        out[idx] = M.missing_name(x) if x >= SYS else format_value(x, fmt, pad=False).strip()
    return MV(out, "string")


@lib("strtoreal", 1)
def _strtoreal(eng, v, r):
    a = _s(v[0])
    out = np.empty(a.shape)
    for idx in np.ndindex(a.shape):
        t = str(a[idx]).strip()
        try:
            out[idx] = M.missing_code(t) if t.startswith(".") and len(t) <= 2 and not t[1:].isdigit() \
                else float(t)
        except ValueError:
            out[idx] = SYS
    return real(out)


@lib("tokens", 1, 2)
def _tokens(eng, v, r):
    s = v[0].str_scalar()
    if len(v) > 1 and v[1] is not None:
        seps = v[1].str_scalar()
        out, cur = [], ""
        for ch in s:
            if ch == " ":
                if cur:
                    out.append(cur)
                cur = ""
            elif ch in seps:
                if cur:
                    out.append(cur)
                out.append(ch)
                cur = ""
            else:
                cur += ch
        if cur:
            out.append(cur)
    else:
        from ..lang.words import split_words
        out = split_words(s)
    arr = np.empty((1, len(out)), dtype=object)
    for k, w in enumerate(out):
        arr[0, k] = w
    return MV(arr, "string")


@lib("invtokens", 1, 2)
def _invtokens(eng, v, r):
    sep = v[1].str_scalar() if len(v) > 1 else " "
    return string(sep.join(str(x) for x in _s(v[0]).ravel()))


@lib("strdup", 2)
def _strdup(eng, v, r):
    return string(v[0].str_scalar() * max(0, _int(v[1])))


@lib("ascii", 1)
def _ascii(eng, v, r):
    return real(np.array([[float(b) for b in v[0].str_scalar().encode("latin-1", "replace")]]).reshape(1, -1))


@lib("char", 1)
def _char(eng, v, r):
    return string("".join(chr(int(x)) for x in _r(v[0]).ravel() if x < SYS))


@lib("strcat", 2)
def _strcat(eng, v, r):
    return LIBRARY_OP("+", v[0], v[1])


def LIBRARY_OP(op, a, b):
    from .ops import binary
    return binary(op, a, b)


# ---------------------------------------------------------------------------
# saída
# ---------------------------------------------------------------------------

_SMCL_TAGS = re.compile(r"\{(txt|text|res|result|err|error|inp|input|com|sf|it|bf|ul off|ul on)\}")


def _smcl_strip(text: str) -> str:
    text = _SMCL_TAGS.sub("", text)
    text = re.sub(r"\{hline (\d+)\}", lambda m: "-" * int(m.group(1)), text)
    text = text.replace("{c |}", "|").replace("{c -}", "-").replace("{c +}", "+")
    return text


_FMT = re.compile(r"%(-|~)?(0)?(\d*)(?:\.(\d+))?(f|g|e|s|gc|fc|x|t[a-zA-Z]*)")


def sprintf(fmt: str, args: list[MV]) -> str:
    out = []
    i = 0
    k = 0
    while i < len(fmt):
        ch = fmt[i]
        if ch == "\\" and i + 1 < len(fmt):
            nxt = fmt[i + 1]
            out.append({"n": "\n", "t": "\t", "\\": "\\", "r": "\r"}.get(nxt, "\\" + nxt))
            i += 2
            continue
        if ch == "%":
            if fmt.startswith("%%", i):
                out.append("%")
                i += 2
                continue
            m = _FMT.match(fmt, i)
            if not m:
                out.append(ch)
                i += 1
                continue
            if k >= len(args):
                raise MataError(3001, "too few arguments")
            arg = args[k]
            k += 1
            spec = m.group(0)
            just, zero, width, dec, kind = m.groups()
            x = arg.scalar() if arg.is_scalar else None
            if x is None:
                raise MataError(3204, f"{arg.orgtype()} found where scalar required")
            if kind == "s":
                if arg.t != "string":
                    raise MataError(3254, "nonstring found where string required")
                text = str(x)
                if width:
                    text = text.ljust(int(width)) if just == "-" else text.rjust(int(width))
                out.append(text)
            else:
                if arg.t != "real":
                    raise MataError(3253, "nonreal found where real required")
                use = spec
                if not width:
                    use = "%" + (just or "") + "9" + ("." + dec if dec else ".0") + kind if kind != "g" else \
                        "%9.0g"
                    text = format_value(float(x), use, pad=False).strip()
                else:
                    text = format_value(float(x), use, pad=True)
                out.append(text)
            i = m.end()
            continue
        out.append(ch)
        i += 1
    return "".join(out)


@lib("printf", 1, 99)
def _printf(eng, v, r):
    text = _smcl_strip(sprintf(v[0].str_scalar(), v[1:]))
    eng.s.output.write(text, "result")
    return None


@lib("sprintf", 1, 99)
def _sprintf(eng, v, r):
    return string(sprintf(v[0].str_scalar(), v[1:]))


@lib("errprintf", 1, 99)
def _errprintf(eng, v, r):
    eng.s.output.write(_smcl_strip(sprintf(v[0].str_scalar(), v[1:])), "error")
    return None


@lib("display", 1, 2)
def _display(eng, v, r):
    a = _s(v[0])
    for x in a.ravel():
        text = str(x) if len(v) > 1 and v[1] is not None and v[1].real_scalar() else _smcl_strip(str(x))
        eng.s.output.write(text + "\n", "result")
    return None


# ---------------------------------------------------------------------------
# erros e utilidades
# ---------------------------------------------------------------------------

class MataExit(Exception):
    def __init__(self, rc: int):
        self.rc = rc


@lib("exit", 0, 1)
def _exit(eng, v, r):
    rc = _int(v[0]) if v else 0
    raise MataExit(rc)


@lib("_error", 1, 2)
def _error_(eng, v, r):
    if v[0].t == "string":
        raise MataError(3498, v[0].str_scalar())
    rc = _int(v[0])
    msg = v[1].str_scalar() if len(v) > 1 else ""
    raise MataError(rc, msg)


@lib("error", 1)
def _error(eng, v, r):
    rc = _int(v[0])
    # VERIFICAR: o Mata imprime a mensagem do Stata para o código
    return real(float(rc))


@lib("assert", 1)
def _assert(eng, v, r):
    if _r(v[0]).size == 0 or not np.all(_r(v[0]) != 0):
        raise MataError(3498, "assertion is false")
    return None


@lib("args", 0)
def _args(eng, v, r):
    if not eng.argc:
        raise MataError(3000, "args() used outside of a function")   # VERIFICAR
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
    return real(clean(np.where(_r(v[0]) >= SYS, SYS, _r(v[0]).astype(np.float32).astype(np.float64))))


@lib("epsilon", 1)
def _epsilon(eng, v, r):
    return real(np.abs(_r(v[0])) * 2.0 ** -52)


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


from . import stata_api  # noqa: E402,F401  (st_* entram na mesma tabela)
