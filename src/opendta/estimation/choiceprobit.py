"""Multinomial probit models by simulated maximum likelihood: mprobit, asmprobit and asroprobit
([R] mprobit, [R] asmprobit, [R] asroprobit). Probabilities of the chosen alternative (or of the
ranking) are orthant probabilities of normal differences, simulated with the GHK algorithm on
Hammersley points (core/choice.py). Stata's draws and normalizations are not documented in full:
the estimates agree with Stata's up to simulation error, and the layout lines marked VERIFY wait
for the Stata 14 log of case 1908."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import choice as C
from ..core import missing as M
from ..core.errors import StataError
from ..core.varlist import resolve_name
from . import fvars
from .choice import (_Alt, _asc_design, _asc_finish, _cluster_by_group, _replay, _restrict, display_eq)
from .mleq import (Block, Engine, add_base_row, aux_row, block_from_text, finish_model, lab, names_in, peek, std_right)
from .models2 import _Prep, g
from .regress import comma

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _draws(n: int, dim: int) -> np.ndarray:
    """Hammersley points: (i + .5)/n in the first dimension, Halton in the others."""
    if dim <= 0:
        return np.zeros((n, 1)) + 0.5
    H = C.halton(n, dim)
    H[:, 0] = (np.arange(n) + 0.5) / n
    return H


def _omega(theta: np.ndarray, J: int, base: int, structure: str):
    """Covariance of the J errors for the parameters theta. The base alternative's error is the
    reference (zero); the others carry the covariance Omega_D of the differences to it."""
    m = J - 1
    if structure == "identity":
        # iid errors N(0, 1): covariance of the differences to the base is I + 11'
        Od = np.eye(m) + np.ones((m, m))
    else:
        L = np.zeros((m, m))
        t = 0
        for i in range(m - 1):
            for j in range(i):
                L[i, j] = theta[t]
                t += 1
            L[i, i] = np.exp(theta[t])
            t += 1
        f = np.r_[theta[t:t + m - 1], 1.0]
        L[m - 1, :] = np.sqrt(2.0) * f / np.linalg.norm(f)     # scale: Var of the last difference = 2
        Od = L @ L.T
    Om = np.zeros((J, J))
    nb = [j for j in range(J) if j != base]
    Om[np.ix_(nb, nb)] = Od
    return Om


def _theta_names(J: int) -> list[str]:
    m = J - 1
    out = []
    for i in range(m - 1):
        for j in range(i):
            out.append(f"l{i + 2}_{j + 2}")
        out.append(f"lnl{i + 2}_{i + 2}")
    for j in range(m - 1):
        out.append(f"l{m + 1}_{j + 2}")
    return out


def _ntheta(J: int, structure: str) -> int:
    return 0 if structure == "identity" else (J - 1) * J // 2 - 1


def _chol(S: np.ndarray) -> np.ndarray:
    """Cholesky factor; a covariance that the optimizer pushed to the edge of positive
    definiteness (platform-dependent rounding) is repaired by clipping its eigenvalues."""
    try:
        return np.linalg.cholesky(S)
    except np.linalg.LinAlgError:
        w, v = np.linalg.eigh((S + S.T) / 2)
        return np.linalg.cholesky((v * np.maximum(w, 1e-10)) @ v.T)


def _case_ll(eta: np.ndarray, kind: str, key: np.ndarray, Om: np.ndarray, U: np.ndarray) -> np.ndarray:
    """Simulated log probability per case. kind "choice": key = chosen alternative; kind "rank":
    key = (cases x J) alternatives ordered from best to worst."""
    n, J = eta.shape
    out = np.zeros(n)
    if kind == "choice":
        for k in range(J):
            rows = np.flatnonzero(key == k)
            if not len(rows):
                continue
            others = [j for j in range(J) if j != k]
            T = np.zeros((J - 1, J))
            for r, j in enumerate(others):
                T[r, j], T[r, k] = 1.0, -1.0
            S = T @ Om @ T.T
            L = _chol(S + 1e-12 * np.eye(J - 1))
            mu = eta[rows][:, [k]] - eta[rows][:, others]
            out[rows] = np.log(np.maximum(C.ghk_upper(mu, L, U), 1e-300))
    else:
        seen = {}
        for i in range(n):
            seen.setdefault(tuple(key[i]), []).append(i)
        for order, rows in seen.items():
            rows = np.array(rows)
            T = np.zeros((J - 1, J))
            for r in range(J - 1):
                T[r, order[r + 1]], T[r, order[r]] = 1.0, -1.0       # u_(r+1) - u_(r) < 0
            S = T @ Om @ T.T
            L = _chol(S + 1e-12 * np.eye(J - 1))
            mu = -(eta[rows][:, [order[r + 1] for r in range(J - 1)]]
                   - eta[rows][:, [order[r] for r in range(J - 1)]])
            out[rows] = np.log(np.maximum(C.ghk_upper(mu, L, U), 1e-300))
    return out


def _numeric_lik(fun, etas_shape_rows, n_theta):
    """Wraps `fun(eta_row_matrix, theta) -> per-case ll` into the (ll, d_eta, d_theta) triple the
    engine needs, by central differences (the simulated likelihood is smooth: draws are fixed)."""
    raise NotImplementedError


# ---------------------------------------------------------------------------
# asmprobit and asroprobit
# ---------------------------------------------------------------------------

def _cmd_asprobit(s: "Session", args: str, cmd: str) -> None:
    if _replay(s, cmd, args):
        return
    ranking = cmd == "asroprobit"
    spec = {"case": 4, "alternatives": 3, "casevars": 5, "basealternative": 3, "correlation": 4, "stddev": 3,
            "structural": 6, "noconstant": 4, "intpoints": 4, "intmethod": 4, "intburn": 4, "intseed": 4,
            "altwise": 3, "scalealternative": 5, "reverse": 3, "noorder": 4}
    d = _asc_design(s, args, spec, cmd, ranking=ranking)
    pr, o = d.pr, d.o
    J = len(d.levels)
    if not np.all(d.nalt == J):
        raise StataError(459, "alternatives are not balanced across cases; altwise is not available yet")   # VERIFY
    corr = str(o.get("correlation", "unstructured")).strip().lower()
    sd = str(o.get("stddev", "heteroskedastic")).strip().lower()
    structure = "identity" if (corr.startswith("ind") and sd.startswith("hom")) else "unstructured"
    if structure == "unstructured" and not (corr.startswith("uns") and sd.startswith("het")) \
            and not (corr.startswith("ind") and sd.startswith("hom")):
        raise StataError(198, "only correlation(unstructured) stddev(heteroskedastic) and "
                              "correlation(independent) stddev(homoskedastic) are available")   # VERIFY
    npts = int(o["intpoints"]) if o.get("intpoints") else 200
    U = _draws(npts, J - 1)
    n_theta = _ntheta(J, structure)
    cases = d.cases
    pos = np.array([idx[np.argsort(d.acode[idx])] for idx in cases])        # (cases x J) row index by alternative
    if ranking:
        rk = d.y[pos]
        if o.get("reverse"):
            rk = -rk
        key = np.argsort(rk, axis=1, kind="stable")        # lowest value = best (VERIFY direction)
        kind = "rank"
    else:
        key = np.argmax(d.y[pos], axis=1)
        kind = "choice"
    first = pos[:, 0]
    w = pr.w

    def case_ll(eta_rows, th):
        eta = eta_rows[pos]
        Om = _omega(th, J, d.base, structure)
        return _case_ll(eta, kind, key, Om, U)

    def lik(etas, th):
        eta = np.zeros(d.n)
        for e in etas:
            eta = eta + e
        th = np.asarray(th, dtype=float)
        base_ll = case_ll(eta, th)
        ll = np.zeros(d.n)
        ll[first] = base_ll
        h = 1e-5
        dd = np.zeros(d.n)
        for j in range(J):
            up, dn = eta.copy(), eta.copy()
            up[pos[:, j]] += h
            dn[pos[:, j]] -= h
            dd[pos[:, j]] = (case_ll(up, th) - case_ll(dn, th)) / (2 * h)
        dth = None
        if n_theta:
            dth = np.zeros((d.n, n_theta))
            for a in range(n_theta):
                tu, td = th.copy(), th.copy()
                tu[a] += h
                td[a] -= h
                dth[first, a] = (case_ll(eta, tu) - case_ll(eta, td)) / (2 * h)
        return ll, [dd] * len(etas), dth

    eng = Engine(lik, d.blocks, w, n_theta)
    b0 = np.zeros(eng.K)
    if n_theta:
        b0[eng.offs[-1]:] = 0.0
    res = pr.maximize(eng.fun, b0)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    names = _theta_names(J) if n_theta else []

    def aux_rows(level):
        if not n_theta:
            return []
        from .results import CoefRow
        rows = [CoefRow("sep")]
        for k, nm in enumerate(names):
            rows.append(aux_row("/" + nm, float(res.b[eng.offs[-1] + k]),
                                float(np.sqrt(max(V[eng.offs[-1] + k, eng.offs[-1] + k], 0))), level))
        return rows

    est = _asc_finish(s, args, cmd,
                      "Alternative-specific rank-ordered probit" if ranking else "Alternative-specific multinomial probit",
                      d, eng, res, V, n_clust, SYS, aux_names=[("", nm) for nm in names],
                      aux_rows_fn=aux_rows, extra_scalars=[("ll_c", res.ll)],
                      extra_macros=[("covariance", "unstructured" if structure != "identity" else "identity"),
                                    ("intmethod", "hammersley")])
    est.coef_title, est.eform_on = "Coef.", False
    est.ll_label = "Log simulated-likelihood"
    est.int_lines = ("Hammersley", str(npts))
    if names:
        lb = d.lev_names
        est.footer = lambda s_, e_: [s_.output.write(f"(alt={lb[k]} is the alternative normalizing {w})\n", "text")
                                     for k, w in ((d.base, "location"), (min(1, J - 1) if d.base == 0 else 0, "scale"))]
    display_eq(s, est)


def cmd_asmprobit(s: "Session", args: str) -> None:
    _cmd_asprobit(s, args, "asmprobit")


def cmd_asroprobit(s: "Session", args: str) -> None:
    _cmd_asprobit(s, args, "asroprobit")


# ---------------------------------------------------------------------------
# mprobit
# ---------------------------------------------------------------------------

def cmd_mprobit(s: "Session", args: str) -> None:
    if _replay(s, "mprobit", args):
        return
    spec = {"baseoutcome": 3, "noconstant": 4, "intpoints": 4, "intmethod": 4, "intburn": 4, "intseed": 4,
            "probitparam": 5, "constraints": 3, "noskip": 4, "or": 2}
    pr = _Prep(s, args, spec, extra_vars="")
    pr.write_notes(s.output)
    y = pr.y
    levels = np.unique(y)
    J = len(levels)
    if J < 2:
        raise StataError(148, "outcome does not vary")     # VERIFY
    if J > 8:
        raise StataError(149, "too many outcomes")         # VERIFY
    code = np.searchsorted(levels, y)
    o = pr.o
    if o.get("baseoutcome"):
        bv = float(o["baseoutcome"])
        if bv not in levels:
            raise StataError(459, f"baseoutcome({o['baseoutcome']}) is not an outcome")  # VERIFY
        base = int(np.flatnonzero(levels == bv)[0])
    else:
        base = int(np.argmax(np.bincount(code)))
    from .mleq import block_from_pr
    n = pr.smp.n
    blocks = []
    lev_names = [g(float(v), "%9.0g") for v in levels]
    for j in range(J):
        if j == base:
            continue
        blk = block_from_pr(pr, lev_names[j])
        blocks.append(blk)
    U = _draws(int(o["intpoints"]) if o.get("intpoints") else 200, J - 1)
    Om = _omega(np.zeros(0), J, base, "identity")
    non_base = [j for j in range(J) if j != base]
    w = pr.w
    h = 1e-5

    def eta_matrix(etas):
        E = np.zeros((n, J))
        for e, j in zip(etas, non_base):
            E[:, j] = e
        return E

    def lik(etas, th):
        E = eta_matrix(etas)
        ll = _case_ll(E, "choice", code, Om, U)
        ds = []
        for j in non_base:
            up, dn = E.copy(), E.copy()
            up[:, j] += h
            dn[:, j] -= h
            ds.append((_case_ll(up, "choice", code, Om, U) - _case_ll(dn, "choice", code, Om, U)) / (2 * h))
        return ll, ds, None

    eng = Engine(lik, blocks, w, 0)
    res = pr.maximize(eng.fun, np.zeros(eng.K))
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    idx = []
    for j, blk in enumerate(blocks):
        idx += list(range(eng.offs[j], eng.offs[j] + blk.nslopes))

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [est.title, "", f"Log likelihood = {g(x['ll'], '%10.0g').strip():>9}"]
        right = [("Number of obs", comma(x["N"])), (f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left, right

    est = finish_model(s, pr, "mprobit", args, title="Multinomial probit regression", blocks=blocks, eng=eng,
                       res=res, V=V, n_clust=n_clust, ll0=SYS, aux_names=[], aux_rows_fn=lambda lev: [],
                       hdr_fn=hdr, wald=True, test_idx=idx, df_m=len(idx), k_eq=len(blocks), heads=True,
                       noeform_after_first=False,
                       extra_scalars=[("k_out", float(J)), ("i_base", float(base + 1)), ("N_cd", float(n))],
                       extra_macros=[("baseoutcome", lev_names[base])])
    nxt = next((lev_names[j] for j in range(base + 1, J)), None)
    add_base_row(est, lev_names[base], "(base outcome)", nxt)
    est.coef_title, est.eform_on = "Coef.", False
    display_eq(s, est)
