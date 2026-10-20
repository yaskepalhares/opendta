"""Likelihood maximization in Stata's style ([R] maximize, [R] ml).

Newton-Raphson with the step observed in Mata's optimize() (compat 0406):
if the full step improves, it advances in increments of 1/8, 1/4, 1/2...
while it improves; if it gets worse, it halves ("backed up"). Convergence:
(mreldif(b) < tolerance or reldif(ll) < ltolerance) and g H⁻¹ g' < nrtolerance.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from ..core.errors import StataError
from ..core.formats import format_value


@dataclass
class MLResult:
    b: np.ndarray
    ll: float
    g: np.ndarray
    H: np.ndarray
    iterations: int
    converged: bool
    log: list[tuple[int, float, str]] = field(default_factory=list)


def _g10(x: float) -> str:
    return format_value(x, "%10.0g", pad=False).strip()


def maximize(fun: Callable[[np.ndarray, int], tuple], b0: np.ndarray, *, maxiter: int = 16000,
             tol: float = 1e-6, ltol: float = 1e-7, nrtol: float = 1e-5, start_iter: int = 0,
             first_ll: float | None = None, stop_at_zero: bool = False,
             nr_stop: float | None = None) -> MLResult:
    """fun(b, todo) -> (ll, g, H) with todo 0 (ll only) or 2 (ll, g and H).
    Returns the iterations for the log ("Iteration k:   log likelihood = ...")."""
    b = np.asarray(b0, dtype=np.float64).copy()
    ll, g, H = fun(b, 2)
    if not np.isfinite(ll):
        raise StataError(430, "initial values not feasible")   # VERIFY
    log: list[tuple[int, float, str]] = []
    it = start_iter
    converged = False
    k = b.size
    # already at the maximum (e.g. gaussian glm from IRLS): Stata stops at iteration 0 (compat
    # 0509)
    if stop_at_zero and k and np.all(np.abs(g) <= 1e-9 * (1 + abs(ll))):
        log.append((it, ll, ""))
        return MLResult(b, ll, g, H, it, True, log)
    while True:
        concave = True
        try:
            np.linalg.cholesky(-H)
        except np.linalg.LinAlgError:
            concave = False
        Hm = H
        if not concave:
            # non-concave Hessian: adds a multiple of the identity (VERIFY: Stata uses a
            # combination of gradient and modified Hessian)
            lam = max(1e-8, float(np.max(np.linalg.eigvalsh(H))) * 1.1 + 1e-6)
            Hm = H - lam * np.eye(k)
        log.append((it, ll, "  (not concave)" if not concave else ""))
        if it - start_iter >= maxiter:
            break
        step = -np.linalg.solve(Hm, g) if k else np.zeros(0)
        t = 1.0
        backed = False
        for _ in range(60):
            bn = b + t * step
            fn = fun(bn, 0)[0]
            if np.isfinite(fn) and fn >= ll - 1e-12 * max(1.0, abs(ll)):
                break
            t /= 2
            backed = True
        else:
            break
        if t == 1.0 and fn > ll:
            inc = 0.125
            for _ in range(60):
                bt = b + (t + inc) * step
                ft = fun(bt, 0)[0]
                if not (np.isfinite(ft) and ft > fn):
                    break
                t, bn, fn = t + inc, bt, ft
                inc *= 2
        if backed:
            log[-1] = (log[-1][0], log[-1][1], log[-1][2] + "  (backed up)")
        dp = float(np.max(np.abs(bn - b) / (np.abs(b) + 1))) if k else 0.0
        dv = abs(fn - ll) / (abs(ll) + 1)
        b = bn
        ll, g, H = fun(b, 2)
        it += 1
        try:
            nr = float(abs(g @ np.linalg.solve(H, g))) if k else 0.0
            np.linalg.cholesky(-H)
            concave_now = True
        except np.linalg.LinAlgError:
            nr, concave_now = np.inf, False
        # nr_stop: the iteration is final when g H^-1 g' falls below it, whatever the change since
        # the previous iteration (stcox; compat 0902)
        if ((dp < tol or dv < ltol) and nr < nrtol or nr_stop is not None and nr < nr_stop) \
                and concave_now:
            converged = True
            log.append((it, ll, ""))
            break
    res = MLResult(b, ll, g, H, it, converged, log)
    if converged and k:
        from .numerics import precise
        if precise():
            refine(fun, res)      # set numerics precise
    return res


def refine(fun, res: MLResult, steps: int = 5) -> None:
    """Newton steps outside the iteration log, while the log likelihood does not fall and the
    gradient keeps shrinking: they take b to the exact maximum (set numerics precise; stcox's
    "Refining estimates")."""
    for _ in range(steps):
        try:
            bn = res.b - np.linalg.solve(res.H, res.g)
        except np.linalg.LinAlgError:
            return
        lln, gn, Hn = fun(bn, 2)
        if not np.isfinite(lln) or lln < res.ll - 1e-12 * max(1.0, abs(res.ll)) or \
                np.linalg.norm(gn) >= np.linalg.norm(res.g):
            return
        res.b, res.ll, res.g, res.H = bn, lln, gn, Hn


def options(o: dict, maxiter: int = 16000) -> dict:
    """maximize() keywords from the iterate(), tolerance(), ltolerance() and nrtolerance()
    options ([R] maximize), with Stata's defaults."""
    return dict(maxiter=int(o["iterate"]) if o.get("iterate") else maxiter,
                tol=float(o["tolerance"]) if o.get("tolerance") else 1e-6,
                ltol=float(o["ltolerance"]) if o.get("ltolerance") else 1e-7,
                nrtol=float(o["nrtolerance"]) if o.get("nrtolerance") else 1e-5)


def maximize_switching(llf: Callable[[np.ndarray], float], scores: Callable[[np.ndarray], np.ndarray],
                       b0: np.ndarray, *, plan=(("bhhh", 5), ("bfgs", 10)), maxiter: int = 16000,
                       tol: float = 1e-6, ltol: float = 1e-7, nrtol: float = 1e-5) -> MLResult:
    """technique(bhhh 5 bfgs 10), the default of arima ([R] maximize): BHHH steps with the outer
    product of the scores, then BFGS updates, alternating until convergence. The log carries the
    "(setting optimization to ...)" and "(switching optimization to ...)" lines as notes with
    iteration number -1. VERIFY: Stata's step-size rules are undocumented; the maximum is the same
    but the iteration values may differ."""
    b = np.asarray(b0, dtype=np.float64).copy()
    k = b.size
    ll = llf(b)
    if not np.isfinite(ll):
        raise StataError(430, "initial values not feasible")   # VERIFY

    def grad(x):
        S_ = scores(x)
        return S_.sum(axis=0), S_

    g, S = grad(b)
    log: list[tuple[int, float, str]] = []
    it = 0
    converged = False
    stage = 0
    tech, left = plan[0]
    log.append((-1, 0.0, f"(setting optimization to {tech.upper()})"))
    W = None            # BFGS approximation of the inverse of -H
    stalled = 0
    while True:
        log.append((it, ll, ""))
        if it >= maxiter:
            break
        if tech == "bhhh" or W is None:
            B = S.T @ S
            try:
                W = np.linalg.inv(B)
            except np.linalg.LinAlgError:
                W = np.linalg.pinv(B)
        step = W @ g
        t = 1.0
        backed = False
        for _ in range(60):
            bn = b + t * step
            fn = llf(bn)
            if np.isfinite(fn) and fn >= ll - 1e-12 * max(1.0, abs(ll)):
                break
            t /= 2
            backed = True
        else:
            break
        if backed:
            log[-1] = (log[-1][0], log[-1][1], "  (backed up)")
        gn, Sn = grad(bn)
        dp = float(np.max(np.abs(bn - b) / (np.abs(b) + 1))) if k else 0.0
        dv = abs(fn - ll) / (abs(ll) + 1)
        if tech == "bfgs":
            sv, yv = bn - b, -(gn - g)          # y: change of the gradient of -ll
            sy = float(sv @ yv)
            if sy > 1e-12:
                rho = 1.0 / sy
                I = np.eye(k)
                W = (I - rho * np.outer(sv, yv)) @ W @ (I - rho * np.outer(yv, sv)) + rho * np.outer(sv, sv)
        b, ll, g, S = bn, fn, gn, Sn
        it += 1
        left -= 1
        stalled = stalled + 1 if dv < 1e-13 else 0
        if stalled >= 25:
            break               # no progress at all (an optimum on the boundary of the parameter space)
        if left == 0:
            # the switch is announced before the next iteration, even when that iteration is the
            # last one (compat 1101)
            stage = (stage + 1) % len(plan)
            tech, left = plan[stage]
            log.append((-1, 0.0, f"(switching optimization to {tech.upper()})"))
        if tech == "bhhh":
            B = S.T @ S
            try:
                nr = float(g @ np.linalg.solve(B, g))
            except np.linalg.LinAlgError:
                nr = np.inf
        else:
            nr = float(g @ W @ g)
        if (dp < tol or dv < ltol) and abs(nr) < nrtol:
            converged = True
            log.append((it, ll, ""))
            break
    return MLResult(b, ll, g, np.zeros((k, k)), it, converged, log)


def print_log(out, res: MLResult, label: str = "log likelihood") -> None:
    for it, ll, note in res.log:
        if it < 0:
            out.write(note + "\n", "text")
            continue
        lead = f"Iteration {it}:"
        out.write(f"{lead:<14} {label} = {_g10(ll):>10}{note}  \n", "text")
    if not res.converged and res.log:
        out.write("convergence not achieved\n", "error")


def ilog(res: MLResult) -> np.ndarray:
    """e(ilog): log likelihood of the first 20 iterations (zeros afterward)."""
    v = np.zeros((1, 20))
    vals = [ll for it, ll, _ in res.log if it >= 0][:20]
    v[0, :len(vals)] = vals
    return v
