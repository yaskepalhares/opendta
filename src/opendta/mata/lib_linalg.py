"""Mata library ([M-5]): linear algebra."""

from __future__ import annotations

import numpy as np

from .libbase import LIBRARY, lib, nan_missing, real_arg, zero_missing
from .values import MV, SYS, clean, conformability, real

# ---------------------------------------------------------------------------
# linear algebra
# ---------------------------------------------------------------------------

def _square(a: np.ndarray) -> None:
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise conformability()


@lib("diag", 1)
def _diag(eng, v, r):
    a = real_arg(v[0])
    if a.shape[0] == 1 or a.shape[1] == 1:
        return real(np.diag(a.ravel()))
    _square(a)
    return real(np.diag(np.diag(a)))


@lib("diagonal", 1)
def _diagonal(eng, v, r):
    return MV(np.diag(v[0].a).reshape(-1, 1).copy(), v[0].t)


@lib("trace", 1)
def _trace(eng, v, r):
    a = real_arg(v[0])
    _square(a)
    return real(float(zero_missing(np.diag(a)).sum()) if a.size else 0.0)


@lib("det", 1)
def _det(eng, v, r):
    a = real_arg(v[0])
    _square(a)
    if np.any(a >= SYS):
        return real(SYS)
    return real(float(np.linalg.det(a)) if a.size else 1.0)


@lib("rank", 1, 2)
def _rank(eng, v, r):
    a = real_arg(v[0])
    if np.any(a >= SYS):
        return real(SYS)
    return real(float(np.linalg.matrix_rank(a)) if a.size else 0.0)


@lib("invsym", 1, 2)
def _invsym(eng, v, r):
    """Generalized inverse of a symmetric matrix: sweeps the columns in order and
    zeroes the rows/columns that turn out collinear (like Stata's invsym():
    for a singular matrix, the redundant columns become zero)."""
    a = real_arg(v[0]).astype(np.float64)
    _square(a)
    if np.any(a >= SYS):
        return real(np.full(a.shape, SYS))
    n = a.shape[0]
    A = (a + a.T) / 2
    keep: list[int] = []
    tol = 1e-13 * max(1.0, float(np.max(np.abs(np.diag(A)))) if n else 1.0)
    order = [int(x) - 1 for x in real_arg(v[1]).ravel()] if len(v) > 1 and v[1] is not None else []
    seq_ = order + [k for k in range(n) if k not in order]
    for k in seq_:
        trial = keep + [k]
        sub = A[np.ix_(trial, trial)]
        # sweep pivot: residual variance of column k given what has already entered
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
    a = real_arg(v[0])
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
    a = real_arg(v[0])
    _square(a)
    try:
        np.linalg.cholesky(a)
        return real(np.linalg.inv(a))
    except np.linalg.LinAlgError:
        return real(np.full(a.shape, SYS))


@lib("pinv", 1)
def _pinv(eng, v, r):
    return real(np.linalg.pinv(nan_missing(real_arg(v[0]))))


@lib("cholesky", 1)
def _cholesky(eng, v, r):
    a = real_arg(v[0])
    _square(a)
    try:
        return real(np.linalg.cholesky(a))
    except np.linalg.LinAlgError:
        return real(np.full(a.shape, SYS))


def _solver(kind):
    def f(eng, v, r):
        A, B = real_arg(v[0]), real_arg(v[1])
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
    a = real_arg(v[0])
    w, X = np.linalg.eigh(a)
    order = np.argsort(-w, kind="stable")
    r[1].set(real(X[:, order])) if r[1] else None
    r[2].set(real(w[order].reshape(1, -1))) if r[2] else None
    return None


@lib("symeigenvalues", 1)
def _symeigenvalues(eng, v, r):
    w = np.linalg.eigvalsh(real_arg(v[0]))
    return real(np.sort(w)[::-1].reshape(1, -1))


@lib("eigensystem", 3, 5, outs=(1, 2))
def _eigensystem(eng, v, r):
    w, X = np.linalg.eig(real_arg(v[0]))
    order = np.argsort(-np.abs(w), kind="stable")
    if r[1]:
        r[1].set(MV(X[:, order], "complex"))
    if r[2]:
        r[2].set(MV(w[order].reshape(1, -1), "complex"))
    return None


@lib("eigenvalues", 1)
def _eigenvalues(eng, v, r):
    w = np.linalg.eigvals(real_arg(v[0]))
    return MV(w[np.argsort(-np.abs(w), kind="stable")].reshape(1, -1), "complex")


@lib("svd", 4, 4, outs=(1, 2, 3))
def _svd(eng, v, r):
    U, s, Vt = np.linalg.svd(real_arg(v[0]), full_matrices=False)
    for k, val in ((1, real(U)), (2, real(s.reshape(-1, 1))), (3, real(Vt))):
        if r[k]:
            r[k].set(val)
    return None


@lib("svdsv", 1)
def _svdsv(eng, v, r):
    return real(np.linalg.svd(real_arg(v[0]), compute_uv=False).reshape(-1, 1))


@lib("norm", 1, 2)
def _norm(eng, v, r):
    a = real_arg(v[0])
    p = v[1].real_scalar() if len(v) > 1 else 2
    if a.shape[0] == 1 or a.shape[1] == 1:
        x = a.ravel()
        return real(float(np.max(np.abs(x))) if p >= SYS else float(np.sum(np.abs(x) ** p) ** (1 / p)))
    if p == 2:
        return real(float(np.linalg.norm(a, 2)))
    return real(float(np.linalg.norm(a, 1 if p == 1 else np.inf if p >= SYS else "fro")))


@lib("lowertriangle", 1, 2)
def _lowertriangle(eng, v, r):
    a = real_arg(v[0]).copy()
    out = np.tril(a)
    if len(v) > 1 and v[1] is not None:
        np.fill_diagonal(out, v[1].real_scalar())
    return real(out)


@lib("uppertriangle", 1, 2)
def _uppertriangle(eng, v, r):
    out = np.triu(real_arg(v[0]).copy())
    if len(v) > 1 and v[1] is not None:
        np.fill_diagonal(out, v[1].real_scalar())
    return real(out)


@lib("makesymmetric", 1)
def _makesymmetric(eng, v, r):
    a = real_arg(v[0])
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
        return MV(real_arg(v[0]) + 1j * real_arg(v[1]), "complex")
    return MV(v[0].a.astype(np.complex128), "complex")
