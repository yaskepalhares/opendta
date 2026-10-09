"""Maximização da verossimilhança no estilo do Stata ([R] maximize, [R] ml).

Newton-Raphson com o passo observado no optimize() do Mata (compat 0406):
se o passo inteiro melhora, avança em incrementos de 1/8, 1/4, 1/2...
enquanto melhorar; se piora, divide por 2 ("backed up"). Convergência:
(mreldif(b) < tolerance ou reldif(ll) < ltolerance) e g H⁻¹ g' < nrtolerance.
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


def maximize(fun: Callable[[np.ndarray, int], tuple], b0: np.ndarray, *, maxiter: int = 300,
             tol: float = 1e-6, ltol: float = 1e-7, nrtol: float = 1e-5, start_iter: int = 0,
             first_ll: float | None = None) -> MLResult:
    """fun(b, todo) -> (ll, g, H) com todo 0 (só ll) ou 2 (ll, g e H).
    Devolve as iterações para o log ("Iteration k:   log likelihood = ...")."""
    b = np.asarray(b0, dtype=np.float64).copy()
    ll, g, H = fun(b, 2)
    if not np.isfinite(ll):
        raise StataError(430, "initial values not feasible")   # VERIFICAR
    log: list[tuple[int, float, str]] = []
    it = start_iter
    converged = False
    k = b.size
    while True:
        concave = True
        try:
            np.linalg.cholesky(-H)
        except np.linalg.LinAlgError:
            concave = False
        Hm = H
        if not concave:
            # Hessiana não côncava: soma um múltiplo da identidade (VERIFICAR:
            # o Stata usa uma combinação de gradiente e Hessiana modificada)
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
        if (dp < tol or dv < ltol) and nr < nrtol and concave_now:
            converged = True
            log.append((it, ll, ""))
            break
    if converged and k:
        from .numerics import precise
        if precise():
            # set numerics precise: passos de Newton extras (fora do log) até o
            # gradiente parar de diminuir, levando b ao máximo exato
            for _ in range(5):
                try:
                    bn = b - np.linalg.solve(H, g)
                except np.linalg.LinAlgError:
                    break
                lln, gn, Hn = fun(bn, 2)
                if not np.isfinite(lln) or lln < ll - 1e-12 * max(1.0, abs(ll)) or \
                        np.linalg.norm(gn) >= np.linalg.norm(g):
                    break
                b, ll, g, H = bn, lln, gn, Hn
    return MLResult(b, ll, g, H, it, converged, log)


def print_log(out, res: MLResult, label: str = "log likelihood") -> None:
    for it, ll, note in res.log:
        out.write(f"Iteration {it}:   {label} = {_g10(ll):>10}{note}  \n", "text")
