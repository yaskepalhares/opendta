"""Mata library ([M-5]): size and type, construction, sums and means, missing values,
shape and selection."""

from __future__ import annotations

import math

import numpy as np

from .libbase import LIBRARY, int_arg, lib, nan_missing, real_arg, zero_missing
from .values import MV, SYS, MataError, clean, conformability, empty, real, string

# ---------------------------------------------------------------------------
# size and type
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
# construction
# ---------------------------------------------------------------------------

@lib("J", 3)
def _J(eng, v, r):
    nr, nc, val = int_arg(v[0]), int_arg(v[1]), v[2]
    if nr < 0 or nc < 0:
        raise MataError(3300, "argument out of range")
    if val.is_scalar:
        return MV(np.full((nr, nc), val.a[0, 0], dtype=val.a.dtype), val.t)
    return MV(np.tile(val.a, (nr, nc)), val.t)


@lib("I", 1, 2)
def _I(eng, v, r):
    n = int_arg(v[0])
    m = int_arg(v[1]) if len(v) > 1 else n
    return real(np.eye(n, m))


@lib("e", 2)
def _e(eng, v, r):
    i, n = int_arg(v[0]), int_arg(v[1])
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
    a, b, n = v[0].real_scalar(), v[1].real_scalar(), int_arg(v[2])
    if n <= 0:
        return empty("real", 0, 1)
    if n == 1:
        return real(np.array([[a]]))
    return real(np.linspace(a, b, n).reshape(-1, 1))


@lib("runiform", 2, 4)
def _runiform(eng, v, r):
    from ..core import rng
    nr, nc = int_arg(v[0]), int_arg(v[1])
    vals = rng.draw("runiform", [x.real_scalar() for x in v[2:]], nr * nc)
    return real(vals.reshape(nr, nc))


@lib("rnormal", 4)
def _rnormal(eng, v, r):
    from ..core import rng
    nr, nc = int_arg(v[0]), int_arg(v[1])
    m, s = v[2], v[3]
    shape = (nr, nc)
    mm = np.broadcast_to(real_arg(m), shape) if m.is_scalar or m.a.shape == shape else None
    ss = np.broadcast_to(real_arg(s), shape) if s.is_scalar or s.a.shape == shape else None
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
        rng.RNG.seed(int_arg(x) % 2**64)
    return None


# ---------------------------------------------------------------------------
# sums, extremes, means
# ---------------------------------------------------------------------------


@lib("sum", 1, 2)
def _sum(eng, v, r):
    a = real_arg(v[0])
    return real(float(zero_missing(a).sum()))


@lib("quadsum", 1, 2)
def _quadsum(eng, v, r):
    return _sum(eng, v, r)


@lib("rowsum", 1, 2)
def _rowsum(eng, v, r):
    return real(zero_missing(real_arg(v[0])).sum(axis=1, keepdims=True))


@lib("colsum", 1, 2)
def _colsum(eng, v, r):
    return real(zero_missing(real_arg(v[0])).sum(axis=0, keepdims=True))


LIBRARY["quadrowsum"] = LIBRARY["rowsum"]
LIBRARY["quadcolsum"] = LIBRARY["colsum"]


def _ext(fn, axis):
    def f(eng, v, r):
        a = real_arg(v[0])
        if a.size == 0:
            return real(SYS) if axis is None else empty("real", 0 if axis == 1 else 1, 0 if axis == 0 else 1)
        x = nan_missing(a)
        import warnings
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            res = fn(x) if axis is None else fn(x, axis=axis, keepdims=True)
        return real(clean(res))
    return f


_max = lib("max", 1)(_ext(np.nanmax, None))
_min = lib("min", 1)(_ext(np.nanmin, None))
_rowmax = lib("rowmax", 1)(_ext(np.nanmax, 1))
_rowmin = lib("rowmin", 1)(_ext(np.nanmin, 1))
_colmax = lib("colmax", 1)(_ext(np.nanmax, 0))
_colmin = lib("colmin", 1)(_ext(np.nanmin, 0))


@lib("minmax", 1)
def _minmax(eng, v, r):
    lo, hi = _min(eng, v, r), _max(eng, v, r)
    return real(np.array([[lo.a[0, 0], hi.a[0, 0]]]))


@lib("colminmax", 1)
def _colminmax(eng, v, r):
    return real(np.vstack([_colmin(eng, v, r).a, _colmax(eng, v, r).a]))


@lib("rowminmax", 1)
def _rowminmax(eng, v, r):
    return real(np.hstack([_rowmin(eng, v, r).a, _rowmax(eng, v, r).a]))


def _complete_rows(*arrays) -> np.ndarray:
    ok = np.ones(arrays[0].shape[0], dtype=bool)
    for a in arrays:
        ok &= ~np.any(a >= SYS, axis=1)
    return ok


def _weights(v, n):
    if len(v) > 1 and v[1] is not None:
        w = real_arg(v[1]).reshape(-1, 1)
        if w.shape[0] == 1 and n != 1:
            w = np.full((n, 1), w[0, 0])
        return w
    return np.ones((n, 1))


@lib("mean", 1, 2)
def _mean(eng, v, r):
    X = real_arg(v[0])
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
    X = real_arg(v[0])
    return real(_var(X, _weights(v, X.shape[0])))


@lib("quadvariance", 1, 2)
def _qvariance(eng, v, r):
    return _variance(eng, v, r)


@lib("correlation", 1, 2)
def _correlation(eng, v, r):
    X = real_arg(v[0])
    return real(_var(X, _weights(v, X.shape[0]), corr=True))


@lib("meanvariance", 1, 2)
def _meanvariance(eng, v, r):
    return real(np.vstack([_mean(eng, v, r).a, _variance(eng, v, r).a]))


@lib("cross", 2, 4)
def _cross(eng, v, r):
    """cross(X, Z), cross(X, w, Z): X'Z (or X'diag(w)Z); rows with missing values
    in X or Z are left out."""
    if len(v) == 2:
        X, Z, w = real_arg(v[0]), real_arg(v[1]), None
    elif len(v) == 3:
        X, w, Z = real_arg(v[0]), real_arg(v[1]), real_arg(v[2])
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
        X, x, Z, z = (real_arg(a) for a in v)
        return real((X - x).T @ (Z - z))
    X, x, w, Z, z = (real_arg(a) for a in v)
    return real(((X - x) * w.reshape(-1, 1)).T @ (Z - z))


@lib("runningsum", 1)
def _runningsum(eng, v, r):
    a = real_arg(v[0])
    flat = np.cumsum(zero_missing(a.ravel()))
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
    return real(1.0 if np.all(real_arg(v[0]) != 0) else 0.0)


@lib("any", 1)
def _any(eng, v, r):
    return real(1.0 if np.any(real_arg(v[0]) != 0) else 0.0)


# ---------------------------------------------------------------------------
# shape and selection
# ---------------------------------------------------------------------------

@lib("rowshape", 2)
def _rowshape(eng, v, r):
    a, n = v[0], int_arg(v[1])
    if n == 0 or a.a.size % n:
        raise conformability()
    return MV(a.a.reshape(n, -1).copy(), a.t)


@lib("colshape", 2)
def _colshape(eng, v, r):
    a, n = v[0], int_arg(v[1])
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
    x = real_arg(v[0]).ravel()
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
    X, s = v[0], real_arg(v[1])
    if s.shape[1] == 1 and s.shape[0] == X.rows:
        keep = s[:, 0] != 0
        return MV(X.a[keep, :].copy(), X.t)
    if s.shape[0] == 1 and s.shape[1] == X.cols:
        keep = s[0, :] != 0
        return MV(X.a[:, keep].copy(), X.t)
    raise conformability()


@lib("selectindex", 1)
def _selectindex(eng, v, r):
    a = real_arg(v[0])
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
    order = _sort_order(X, real_arg(v[1]))
    return MV(X.a[order].copy(), X.t)


@lib("_sort", 2, outs=(0,))
def _sort_in_place(eng, v, r):
    res = _sort(eng, v, r)
    if r[0] is not None:
        r[0].set(res)
    return None


@lib("order", 2)
def _order(eng, v, r):
    return real((_sort_order(v[0], real_arg(v[1])) + 1.0).reshape(-1, 1))


@lib("revorder", 1)
def _revorder(eng, v, r):
    p = real_arg(v[0])
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
