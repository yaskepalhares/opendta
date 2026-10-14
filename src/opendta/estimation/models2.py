"""Modelos com mais de uma equação ou parâmetros auxiliares:
ologit, oprobit ([R] ologit, [R] oprobit), mlogit ([R] mlogit),
nbreg ([R] nbreg) e tobit ([R] tobit)."""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np
from scipy import special, stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from . import fvars
from .ml import maximize, print_log
from .regress import Sample, cmdline, comma, wexp_macros
from .results import CoefRow, Estimates, coef_table, crit, post, rows_from_columns

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def g(x: float, fmt: str) -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def num_hess(grad: Callable[[np.ndarray], np.ndarray], b: np.ndarray) -> np.ndarray:
    """Hessiana por diferenças centrais do gradiente analítico, com extrapolação
    de Richardson (erro de truncamento O(h⁴))."""
    k = b.size
    H = np.zeros((k, k))
    for i in range(k):
        h = 1e-4 * (abs(b[i]) + 1)

        def D(hh, i=i):
            a, c = b.copy(), b.copy()
            a[i] += hh
            c[i] -= hh
            return (grad(a) - grad(c)) / (2 * hh)
        H[:, i] = (4 * D(h / 2) - D(h)) / 3
    return (H + H.T) / 2


class _Prep:
    """Amostra, pesos e colunas das variáveis independentes (sem constante)."""

    def __init__(self, s: "Session", args: str, spec: dict[str, int], weights=("fweight", "iweight",
                 "pweight"), constant: bool | None = None):
        full = dict(spec)
        full.update({"nolog": 4, "log": 3, "iterate": 4, "tolerance": 3, "ltolerance": 4,
                     "nrtolerance": 4, "difficult": 4, "technique": 4, "noheader": 6})
        smp = Sample(s, args, opts_spec=full, weights=weights)
        self.smp = smp
        self.o = smp.o
        w = smp.w[smp.mask].astype(np.float64) if smp.w is not None else np.ones(smp.n)
        self.w = w
        self.N = float(w.sum()) if smp.wtype in ("fweight", "iweight") else float(smp.n)
        self.y = smp.y.astype(np.float64)
        self.constant = smp.constant if constant is None else constant
        self.notes = fvars.drop_collinear(smp.cols, self.constant, w)
        self.est_cols = [c for c in smp.cols if c.values is not None]
        self.X = np.column_stack([c.values for c in self.est_cols]) if self.est_cols \
            else np.zeros((smp.n, 0))
        self.robust = smp.vce in ("robust", "cluster")

    def write_notes(self, out) -> None:
        for nm in self.notes:
            out.write(f"note: {nm} omitted because of collinearity\n", "text")

    def maximize(self, fun, b0):
        o = self.o
        return maximize(fun, b0, maxiter=int(o["iterate"]) if o.get("iterate") else 300,
                        tol=float(o["tolerance"]) if o.get("tolerance") else 1e-6,
                        ltol=float(o["ltolerance"]) if o.get("ltolerance") else 1e-7,
                        nrtol=float(o["nrtolerance"]) if o.get("nrtolerance") else 1e-5)

    def vce(self, H: np.ndarray, scores: np.ndarray | None) -> tuple[np.ndarray, int | None]:
        """V do modelo (OIM) ou sanduíche com os escores por observação (n x k)."""
        try:
            V = np.linalg.inv(-H)
        except np.linalg.LinAlgError:
            V = np.linalg.pinv(-H)
        V = (V + V.T) / 2
        if not self.robust or scores is None:
            return V, None
        smp = self.smp
        sc = scores * self.w[:, None] if smp.wtype != "fweight" else scores
        n_clust = None
        if smp.vce == "cluster":
            uniq = np.unique(smp.cluster)
            n_clust = len(uniq)
            meat = np.zeros_like(V)
            for u in uniq:
                sg = sc[smp.cluster == u].sum(axis=0)
                meat += np.outer(sg, sg)
            q = n_clust / (n_clust - 1)
        else:
            if smp.wtype == "fweight":
                meat = (scores * self.w[:, None]).T @ scores
            else:
                meat = sc.T @ sc
            q = self.N / (self.N - 1)
        V = q * V @ meat @ V
        return (V + V.T) / 2, n_clust


def _header(s: "Session", left: list[str], right: list[tuple[str, str]]) -> None:
    out = s.output
    out.write("\n", "text")
    n = max(len(left), len(right))
    left = left + [""] * (n - len(left))
    right = right + [("", "")] * (n - len(right))
    for lft, (rl, rv) in zip(left, right):
        if rl:
            out.write(f"{lft:<48}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        else:
            out.write(lft + "\n", "text")
    out.write("\n", "text")


def _std_right(N: float, chi2type: str, df: int, chi2: float, p: float, r2p: float | None):
    right = [("Number of obs", comma(N)), (f"{chi2type} chi2({df})", g(chi2, "%10.2f")),
             ("Prob > chi2", g(p, "%10.4f"))]
    if r2p is not None:
        right.append(("Pseudo R2", g(r2p, "%10.4f")))
    return right


def _model_test(pr: _Prep, b: np.ndarray, V: np.ndarray, idx: list[int], ll: float, ll0: float):
    if pr.robust:
        try:
            chi2 = float(b[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], b[idx])) if idx else 0.0
        except np.linalg.LinAlgError:
            chi2 = SYS
        kind = "Wald"
    else:
        chi2 = 2 * (ll - ll0) if ll0 < SYS else SYS
        kind = "LR"
    p = float(st.chi2.sf(chi2, len(idx))) if chi2 < SYS and idx else SYS
    return kind, chi2, p


def _slope_rows(pr: _Prep, b: np.ndarray, V: np.ndarray, offset: int = 0):
    """Linhas e nomes das inclinações (sem constante), a partir de b/V do bloco."""
    smp = pr.smp
    bmap, semap = {}, {}
    names = []
    j = 0
    for c in smp.cols:
        if c.values is not None:
            bmap[c.name] = b[offset + j]
            v = V[offset + j, offset + j]
            semap[c.name] = np.sqrt(v) if v > 0 else 0.0
            j += 1
        else:
            bmap[c.name], semap[c.name] = 0.0, 0.0
        names.append(c.name if not c.omitted else fvars.omitted_name(c.name))
    rows = rows_from_columns(smp.cols, bmap, semap, s_data(pr), constant=False)
    return names, rows


def s_data(pr: _Prep):
    return pr.smp.s.data


def _expand(pr: _Prep, b_est: np.ndarray, V_est: np.ndarray, layout: list[tuple[str, bool]]):
    """Coloca b/V estimados (só colunas livres) nas posições de e(b), que inclui
    bases e omitidas com zero. `layout` = [(nome, livre?)]."""
    kk = len(layout)
    pos = [i for i, (_, free) in enumerate(layout) if free]
    b = np.zeros(kk)
    V = np.zeros((kk, kk))
    b[pos] = b_est
    V[np.ix_(pos, pos)] = V_est
    return b, V


# ordem de e() observada no Stata 14 (compat 0508, 0509); para os comandos sem
# log de referência, a ordem segue o padrão dos parecidos (VERIFICAR)
_E_ORDER = {
    "ologit": (["rank", "N", "ic", "k", "k_eq", "k_dv", "converged", "rc", "ll", "N_clust", "k_eq_model",
                "ll_0", "df_m", "chi2", "p", "N_cd", "k_cat", "k_aux", "r2_p"],
               ["cmdline", "cmd", "predict", "marginsdefault", "title", "chi2type", "opt", "vce", "vcetype",
                "clustvar", "user", "ml_method", "technique", "which", "wtype", "wexp", "depvar", "properties"],
               ["cat", "ilog", "gradient"]),
    "mlogit": (["rank", "N", "ic", "k", "k_dv", "converged", "rc", "ll", "N_clust", "ll_0", "df_m", "chi2", "p",
                "N_cd", "ibaseout", "baseout", "k_out", "r2_p", "k_eq", "k_eq_model", "k_eq_base"],
               ["cmdline", "cmd", "marginsdefault", "marginsnotok", "predict", "title", "eqnames", "baselab",
                "chi2type", "opt", "vce", "vcetype", "clustvar", "user", "ml_method", "technique", "which",
                "wtype", "wexp", "depvar", "properties"],
               ["out", "ilog", "gradient"]),
    "nbreg": (["rank", "N", "ic", "k", "k_eq", "k_dv", "converged", "rc", "ll", "N_clust", "k_eq_model",
               "ll_0", "df_m", "chi2", "p", "k_aux", "alpha", "ll_c", "chi2_c", "p_c", "r2_p"],
              ["cmdline", "cmd", "predict", "dispers", "diparm1", "title", "chi2type", "chi2_ct", "opt", "vce",
               "vcetype", "clustvar", "user", "ml_method", "technique", "which", "wtype", "wexp", "depvar",
               "properties"],
              ["ilog", "gradient"]),
    "tobit": (["rank", "N", "ic", "k", "k_eq", "k_dv", "converged", "rc", "ll", "N_clust", "k_eq_model",
               "ll_0", "df_m", "chi2", "p", "df_r", "N_unc", "N_lc", "N_rc", "llopt", "ulopt", "k_aux", "r2_p"],
              ["cmdline", "cmd", "predict", "title", "chi2type", "opt", "vce", "vcetype", "clustvar", "user",
               "ml_method", "technique", "which", "wtype", "wexp", "depvar", "properties"],
              ["ilog", "gradient"]),
}
_E_ORDER["oprobit"] = _E_ORDER["ologit"]


def _finish(s: "Session", pr: _Prep, cmd: str, args: str, *, title: str, names: list[str], eqnames,
            b: np.ndarray, V: np.ndarray, rows: list[CoefRow], res, ll0: float, df_m: int, chi2kind,
            chi2, p, r2p, display, n_clust, extra_scalars=(), extra_macros=(), k_eq=1, k_aux=0,
            depvar_name: str | None = None, constant: bool = False, mats=(), k_free=None,
            order: tuple | None = None) -> Estimates:
    from ..commands.matrix import Matrix
    from .ml import ilog
    smp = pr.smp
    est = Estimates(cmd, depvar_name or smp.depname, names, b, V, smp.n, "z", None)
    est.rows = rows
    est.eqnames = list(eqnames)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if pr.robust else ""
    est.clustvar = smp.clustvar
    est.n_clust = n_clust
    est.title = title
    est.robust = pr.robust
    est.extra = {"N": pr.N, "ll": res.ll, "ll_0": ll0, "chi2": chi2, "p": p, "df_m": df_m,
                 "r2_p": r2p, "chi2type": chi2kind}
    est.display = display
    est.ml_cmd = cmd
    sv = {"rank": float(np.linalg.matrix_rank(V)) if V.size else 0.0, "N": pr.N, "ic": float(res.iterations),
          "k": float(k_free if k_free is not None else len(names)), "k_eq": float(k_eq), "k_dv": 1.0,
          "converged": 1.0 if res.converged else 0.0, "rc": 0.0, "ll": res.ll, "k_eq_model": 1.0,
          "ll_0": ll0, "df_m": float(df_m), "chi2": chi2, "p": p}
    if n_clust is not None:
        sv["N_clust"] = float(n_clust)
    if r2p is not None:
        sv["r2_p"] = r2p
    if k_aux:
        sv["k_aux"] = float(k_aux)
    sv.update(dict(extra_scalars))
    mv = {"cmdline": cmdline(cmd, args), "cmd": cmd, "title": title, "chi2type": chi2kind,
          "vce": smp.vce if pr.robust else "oim", "depvar": smp.depname, "opt": "moptimize", "which": "max",
          "ml_method": "d2", "technique": "nr", "properties": "b V", "predict": f"{cmd}_p",
          "user": f"mopt__{cmd}_d2()"}
    if pr.robust:
        mv["vcetype"] = "Robust"
    if smp.clustvar:
        mv["clustvar"] = smp.clustvar
    mv.update(dict(wexp_macros(smp)))
    mv.update(dict(extra_macros))
    md = {"ilog": Matrix(ilog(res)), "gradient": Matrix(res.g.reshape(1, -1))}
    md.update(dict(mats))
    so, mo, xo = order or _E_ORDER.get(cmd, (list(sv), list(mv), list(md)))
    scalars = [(k_, sv[k_]) for k_ in so if k_ in sv]
    macros = [(k_, mv[k_]) for k_ in mo if k_ in mv]
    if order is None:
        # sem ordem explícita, o que não está na lista vai para o fim
        scalars += [(k_, v) for k_, v in sv.items() if k_ not in so]
        macros += [(k_, v) for k_, v in mv.items() if k_ not in mo]
    mlist = [(k_, md[k_]) for k_ in xo if k_ in md]
    post(s, est, scalars, macros, smp.mask, extra_mats=mlist)
    return est


def _log(s: "Session", pr: _Prep, res, label: str | None = None) -> None:
    if pr.o.get("nolog"):
        return
    out = s.output
    out.write("\n", "text")
    print_log(out, res, label or ("log pseudolikelihood" if pr.robust else "log likelihood"))


def _cluster_note(s, est) -> None:
    if est.vcetype and est.clustvar:
        s.output.write(f"{'(Std. Err. adjusted for ' + format(int(est.n_clust), ',') + ' clusters in ' + est.clustvar + ')':>78}\n",
                       "text")


# ---------------------------------------------------------------------------
# ologit / oprobit
# ---------------------------------------------------------------------------

def _ordered(s: "Session", cmd: str, args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return
    pr = _Prep(s, args, {"offset": 3}, constant=False)
    out = s.output
    pr.write_notes(out)
    y = pr.y
    levels = sorted(set(y.tolist()))
    if len(levels) < 2:
        raise StataError(2000, "outcome does not vary")   # VERIFICAR
    if len(levels) > 50:
        raise StataError(149, f"{smp_name(pr)}:  too many outcomes")   # VERIFICAR
    J = len(levels)
    yi = np.searchsorted(levels, y)
    X, w = pr.X, pr.w
    k = X.shape[1]
    if cmd == "ologit":
        F = special.expit
        f = lambda z: special.expit(z) * (1 - special.expit(z))   # noqa: E731
        Finv = special.logit
    else:
        F = st.norm.cdf
        f = st.norm.pdf
        Finv = st.norm.ppf

    def probs(b):
        xb = X @ b[:k]
        cuts = b[k:]
        up = np.where(yi < J - 1, cuts[np.minimum(yi, J - 2)] - xb, np.inf)
        lo = np.where(yi > 0, cuts[np.maximum(yi - 1, 0)] - xb, -np.inf)
        return up, lo

    def ll_obs(b):
        up, lo = probs(b)
        p = F(up) - F(lo)
        return np.log(np.maximum(p, 1e-300))

    def scores(b):
        up, lo = probs(b)
        p = np.maximum(F(up) - F(lo), 1e-300)
        fu = np.where(np.isfinite(up), f(np.where(np.isfinite(up), up, 0)), 0.0)
        fl = np.where(np.isfinite(lo), f(np.where(np.isfinite(lo), lo, 0)), 0.0)
        S = np.zeros((len(y), k + J - 1))
        S[:, :k] = -((fu - fl) / p)[:, None] * X
        for j in range(J - 1):
            S[:, k + j] = np.where(yi == j, fu / p, 0.0) - np.where(yi == j + 1, fl / p, 0.0)
        return S

    def grad(b):
        return scores(b).T @ w

    def fun(b, todo):
        cuts = b[k:]
        if np.any(np.diff(cuts) <= 0):
            return -np.inf, None, None
        ll = float(w @ ll_obs(b))
        if todo == 0:
            return ll, None, None
        return ll, grad(b), num_hess(grad, b)

    cum = np.array([np.sum(w[yi <= j]) for j in range(J - 1)]) / w.sum()
    b0 = np.concatenate([np.zeros(k), Finv(cum)])
    ll0 = float(fun(b0, 0)[0])
    res = pr.maximize(fun, b0)
    _log(s, pr, res)
    V, n_clust = pr.vce(res.H, scores(res.b) if pr.robust else None)
    idx = list(range(k))
    kind, chi2, pv = _model_test(pr, res.b, V, idx, res.ll, ll0)
    r2p = 1 - res.ll / ll0 if ll0 else SYS
    names, rows = _slope_rows(pr, res.b, V)
    layout = [(nm, c.values is not None) for nm, c in zip(names, pr.smp.cols)]
    cutnames = [f"cut{j + 1}" for j in range(J - 1)]
    layout += [("_cons", True)] * (J - 1)
    b_full, V_full = _expand(pr, res.b, V, layout)
    eqs = [pr.smp.depname] * len(names) + cutnames
    rows = list(rows)
    rows.append(CoefRow("sep"))
    for j in range(J - 1):
        se = np.sqrt(V[k + j, k + j])
        rows.append(CoefRow("aux", f"/{cutnames[j]}", float(res.b[k + j]), float(se), name=f"cut{j + 1}"))
    title = "Ordered logistic regression" if cmd == "ologit" else "Ordered probit regression"
    est = _finish(s, pr, cmd, args, title=title, names=names + ["_cons"] * (J - 1), eqnames=eqs,
                  b=b_full, V=V_full, rows=rows, res=res, ll0=ll0, df_m=k, chi2kind=kind, chi2=chi2, p=pv,
                  r2p=r2p, display=_display_std, n_clust=n_clust, k_eq=J, k_aux=J - 1,
                  extra_scalars=[("k_cat", float(J)), ("N_cd", 0.0)],
                  extra_macros=[("marginsdefault", " ".join(f"predict(pr outcome({_lv(v)}))" for v in levels))],
                  mats=[("cat", _mat_row(levels))])
    est.cat = levels
    est.predict = _predict_ordered
    est.kslopes = k
    est.Fcdf = F
    _display_std(s, est)


def _lv(v: float) -> str:
    return f"{int(v)}" if v == int(v) else f"{v:g}"


def _mat_row(vals):
    from ..commands.matrix import Matrix
    return Matrix(np.array(vals, dtype=float).reshape(1, -1))


def smp_name(pr):
    return pr.smp.depname


def _display_std(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                 level: float | None = None, **_kw) -> None:
    x = est.extra
    if header:
        lab = "Log pseudolikelihood" if est.robust else "Log likelihood"
        left = [est.title, "", "", f"{lab} = {g(x['ll'], '%10.0g'):>10}"]
        if getattr(est, "dispersion", ""):
            left[2] = f"Dispersion     = {est.dispersion}"
        r2p = x.get("r2_p")
        _header(s, left, _std_right(x["N"], x["chi2type"], int(x["df_m"]), x["chi2"], x["p"],
                                    r2p if r2p is not None else None))
    if table:
        _cluster_note(s, est)
        coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level,
                   vcetype=est.vcetype, coef_title=getattr(est, "coef_title", "Coef."),
                   eform=getattr(est, "eform_on", False))
        foot = getattr(est, "footer", None)
        if foot:
            foot(s, est)


def _predict_ordered(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat in ("", "pr", "p"):
        # sem outcome(): a primeira categoria (VERIFICAR: o Stata pede uma
        # variável por categoria)
        k = est.kslopes
        cut1 = est.b[est.names.index("_cons", k)] if "_cons" in est.names else 0.0
        p = est.Fcdf(cut1 - xb)
        return np.where(ok, p, SYS), f"(option pr assumed; predicted probability)"   # VERIFICAR
    raise StataError(198, f"option {stat} not allowed")


def cmd_ologit(s, args):
    _ordered(s, "ologit", args)


def cmd_oprobit(s, args):
    _ordered(s, "oprobit", args)


# ---------------------------------------------------------------------------
# mlogit
# ---------------------------------------------------------------------------

def cmd_mlogit(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "mlogit", args)
        return
    pr = _Prep(s, args, {"baseoutcome": 2, "rrr": 3})
    out = s.output
    pr.write_notes(out)
    y = pr.y
    levels = sorted(set(y.tolist()))
    J = len(levels)
    if J < 2:
        raise StataError(2000, "outcome does not vary")   # VERIFICAR
    w = pr.w
    yi = np.searchsorted(levels, y)
    freq = np.array([w[yi == j].sum() for j in range(J)])
    if pr.o.get("baseoutcome"):
        bv = float(str(pr.o["baseoutcome"]).strip().lstrip("#"))
        if str(pr.o["baseoutcome"]).strip().startswith("#"):
            base = int(bv) - 1
        else:
            if bv not in levels:
                raise StataError(198, "outcome not found")   # VERIFICAR
            base = levels.index(bv)
    else:
        base = int(np.argmax(freq))      # a categoria mais frequente (a primeira em empate)
    X = np.column_stack([pr.X, np.ones(len(y))]) if pr.constant else pr.X
    kx = X.shape[1]
    others = [j for j in range(J) if j != base]
    Y = np.zeros((len(y), J))
    Y[np.arange(len(y)), yi] = 1

    def P(b):
        B = b.reshape(len(others), kx)
        eta = np.zeros((len(y), J))
        for m, j in enumerate(others):
            eta[:, j] = X @ B[m]
        eta -= eta.max(axis=1, keepdims=True)
        e = np.exp(eta)
        return e / e.sum(axis=1, keepdims=True)

    def fun(b, todo):
        Pm = P(b)
        ll = float(w @ np.log(np.maximum(Pm[np.arange(len(y)), yi], 1e-300)))
        if todo == 0:
            return ll, None, None
        gr = np.concatenate([X.T @ (w * (Y[:, j] - Pm[:, j])) for j in others])
        H = np.zeros((b.size, b.size))
        for a, ja in enumerate(others):
            for c, jc in enumerate(others):
                wt = w * Pm[:, ja] * ((1.0 if ja == jc else 0.0) - Pm[:, jc])
                H[a * kx:(a + 1) * kx, c * kx:(c + 1) * kx] = -(X * wt[:, None]).T @ X
        return ll, gr, H

    def scores(b):
        Pm = P(b)
        return np.column_stack([X * (Y[:, j] - Pm[:, j])[:, None] for j in others])

    b0 = np.zeros(len(others) * kx)
    if pr.constant:
        for m, j in enumerate(others):
            b0[m * kx + kx - 1] = np.log(freq[j] / freq[base])
    ll0 = float(fun(b0, 0)[0]) if pr.constant else SYS
    res = pr.maximize(fun, b0)
    _log(s, pr, res)
    V, n_clust = pr.vce(res.H, scores(res.b) if pr.robust else None)
    slope_idx = [m * kx + i for m in range(len(others)) for i in range(kx - (1 if pr.constant else 0))]
    kind, chi2, pv = _model_test(pr, res.b, V, slope_idx, res.ll, ll0)
    r2p = 1 - res.ll / ll0 if ll0 < SYS and ll0 else None
    ds = s.data
    depv = ds.get(pr.smp.depcomp.var)
    vlab = ds.value_labels.get(depv.value_label, {}) if depv.value_label else {}

    def eqname(j):
        # rótulos com espaços viram nomes de equação com _ (compat 0508)
        lv = levels[j]
        return vlab.get(int(lv), f"{int(lv)}" if lv == int(lv) else f"{lv:g}").replace(" ", "_")
    # e(b) com todas as equações, inclusive a base (zeros); linhas da tabela
    names_all, eqs_all, free = [], [], []
    rows: list[CoefRow] = []
    bvals, Vidx = [], []
    for j in range(J):
        if rows:
            rows.append(CoefRow("sep"))
        if j == base:
            rows.append(CoefRow("eq", eqname(j), eq="(base outcome)"))
        else:
            rows.append(CoefRow("eq", eqname(j)))
        m = others.index(j) if j != base else None
        slot = 0
        bmap, semap = {}, {}
        for c in pr.smp.cols:
            nm = c.name if not c.omitted else fvars.omitted_name(c.name)
            if m is None and not c.omitted and not c.base:
                nm = fvars.omitted_name(c.name)       # equação de base: o.x1
            names_all.append(nm)
            eqs_all.append(eqname(j))
            if c.values is not None and m is not None:
                Vidx.append(m * kx + slot)
                free.append(True)
                bmap[c.name] = res.b[m * kx + slot]
                semap[c.name] = np.sqrt(max(V[m * kx + slot, m * kx + slot], 0))
                slot += 1
            else:
                free.append(False)
                bmap[c.name], semap[c.name] = 0.0, 0.0
        if pr.constant:
            names_all.append("_cons" if m is not None else "o._cons")
            eqs_all.append(eqname(j))
            if m is not None:
                Vidx.append(m * kx + kx - 1)
                free.append(True)
                bmap["_cons"] = res.b[m * kx + kx - 1]
                semap["_cons"] = np.sqrt(max(V[m * kx + kx - 1, m * kx + kx - 1], 0))
            else:
                free.append(False)
        if m is not None:
            sub = rows_from_columns(pr.smp.cols, bmap, semap, ds, constant=pr.constant)
            for r in sub:
                r.eq = ""
            rows.extend(sub)
    kk = len(names_all)
    b_full = np.zeros(kk)
    V_full = np.zeros((kk, kk))
    pos = [i for i, f_ in enumerate(free) if f_]
    order = Vidx
    b_full[pos] = res.b[order]
    V_full[np.ix_(pos, pos)] = V[np.ix_(order, order)]
    est = _finish(s, pr, "mlogit", args, title="Multinomial logistic regression", names=names_all,
                  eqnames=eqs_all, b=b_full, V=V_full, rows=rows, res=res, ll0=ll0,
                  df_m=len(slope_idx), chi2kind=kind, chi2=chi2, p=pv, r2p=r2p, display=_display_std,
                  n_clust=n_clust, k_eq=J, k_free=len(res.b),
                  extra_scalars=[("k_out", float(J)), ("ibaseout", float(base + 1)),
                                 ("baseout", float(levels[base])), ("N_cd", 0.0), ("k_eq_model", float(J)),
                                 ("k_eq_base", float(base + 1))],
                  extra_macros=[("marginsdefault", " ".join(f"predict(pr outcome({_lv(v)}))" for v in levels)),
                                ("marginsnotok", "stdp stddp SCores"),
                                ("eqnames", " ".join(eqname(j) for j in range(J))),
                                ("baselab", eqname(base))],
                  mats=[("out", _mat_row(levels))])
    est.rows = rows
    if pr.o.get("rrr"):
        est.coef_title, est.eform_on = "RRR", True
    _display_std(s, est)
    est.coef_title, est.eform_on = "Coef.", False


# ---------------------------------------------------------------------------
# nbreg (dispersão média, NB2)
# ---------------------------------------------------------------------------

def cmd_nbreg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "nbreg", args)
        return
    pr = _Prep(s, args, {"dispersion": 4, "irr": 3, "exposure": 3, "offset": 3})
    out = s.output
    pr.write_notes(out)
    y = pr.y
    if np.any(y < 0) or np.any(y != np.trunc(y)):
        raise StataError(459, f"{pr.smp.depname} must be nonnegative integers")   # VERIFICAR
    ds = s.data
    off = np.zeros(len(y))
    if pr.o.get("exposure"):
        off = np.log(ds.get(str(pr.o["exposure"]).strip()).data[pr.smp.mask].astype(float))
    elif pr.o.get("offset"):
        off = ds.get(str(pr.o["offset"]).strip()).data[pr.smp.mask].astype(float)
    X = np.column_stack([pr.X, np.ones(len(y))]) if pr.constant else pr.X
    kx = X.shape[1]
    w = pr.w
    lgy = special.gammaln(y + 1)

    def make(Xm):
        k = Xm.shape[1]

        def ll_obs(b):
            mu = np.exp(np.clip(Xm @ b[:k] + off, -700, 700))
            m = np.exp(-b[k])
            return (special.gammaln(y + m) - special.gammaln(m) - lgy + m * np.log(m / (m + mu))
                    + y * np.log(mu / (m + mu)))

        def scores(b):
            mu = np.exp(np.clip(Xm @ b[:k] + off, -700, 700))
            a = np.exp(b[k])
            m = 1 / a
            sx = (y - mu) / (1 + a * mu)
            dm = special.digamma(y + m) - special.digamma(m) + np.log(m / (m + mu)) + 1 - (m + y) / (m + mu)
            return np.column_stack([Xm * sx[:, None], -m * dm])

        def grad(b):
            return scores(b).T @ w

        def fun(b, todo):
            ll = float(w @ ll_obs(b))
            if todo == 0:
                return ll, None, None
            return ll, grad(b), num_hess(grad, b)
        return fun, scores

    # 1) Poisson, 2) só a constante, 3) modelo completo (como o Stata)
    from .models import _poisson
    pfun, _ = _poisson(y, X, w, off)
    ybar = float(np.sum(w * y) / np.sum(w))
    b0p = np.zeros(kx)
    if pr.constant:
        b0p[-1] = np.log(np.sum(w * y) / np.sum(w * np.exp(off)))
    resp = pr.maximize(pfun, b0p)
    if not pr.o.get("nolog"):
        out.write("\nFitting Poisson model:\n", "text")
    _log(s, pr, resp)
    ll_pois = resp.ll
    fun0, _ = make(np.ones((len(y), 1)))
    b00 = np.array([np.log(np.sum(w * y) / np.sum(w * np.exp(off))), 0.0])
    res0 = pr.maximize(fun0, b00)
    if not pr.o.get("nolog"):
        out.write("\nFitting constant-only model:\n", "text")
    _log(s, pr, res0)
    fun, scores = make(X)
    b0 = np.concatenate([resp.b, [res0.b[1]]])
    res = pr.maximize(fun, b0)
    if not pr.o.get("nolog"):
        out.write("\nFitting full model:\n", "text")
    _log(s, pr, res)
    V, n_clust = pr.vce(res.H, scores(res.b) if pr.robust else None)
    slope = list(range(kx - (1 if pr.constant else 0)))
    kind, chi2, pv = _model_test(pr, res.b, V, slope, res.ll, res0.ll)
    r2p = 1 - res.ll / res0.ll
    names, rows = _slope_rows(pr, res.b, V)
    ds = s.data
    bmap = {}
    rows = rows_from_columns(pr.smp.cols, _bmap(pr, res.b, kx), _semap(pr, V, kx), ds,
                             constant=pr.constant)
    rows.append(CoefRow("sep"))
    lna, se_lna = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    c = crit("z", pr.smp.level, None)
    rows.append(CoefRow("aux", "/lnalpha", lna, se_lna, name="_cons"))
    rows.append(CoefRow("sep"))
    alpha = np.exp(lna)
    rows.append(CoefRow("aux", "alpha", alpha, alpha * se_lna,
                        eq=f"{float(np.exp(lna - c * se_lna))!r} {float(np.exp(lna + c * se_lna))!r}"))
    layout = [(nm, cc.values is not None) for nm, cc in zip(names, pr.smp.cols)]
    if pr.constant:
        layout.append(("_cons", True))
        names = names + ["_cons"]
    layout.append(("_cons", True))
    b_full, V_full = _expand(pr, res.b, V, layout)
    eqs = [pr.smp.depname] * len(names) + ["lnalpha"]
    chibar = max(2 * (res.ll - ll_pois), 0.0)
    p_chibar = 0.5 * float(st.chi2.sf(chibar, 1)) if chibar > 0 else 0.5
    est = _finish(s, pr, "nbreg", args, title="Negative binomial regression", names=names + ["_cons"],
                  eqnames=eqs, b=b_full, V=V_full, rows=rows, res=res, ll0=res0.ll, df_m=len(slope),
                  chi2kind=kind, chi2=chi2, p=pv, r2p=r2p, display=_display_std, n_clust=n_clust, k_eq=2,
                  k_aux=1, extra_scalars=[("alpha", alpha), ("chi2_c", chibar), ("p_c", p_chibar),
                                          ("ll_c", ll_pois)],
                  extra_macros=[("dispers", "mean"), ("chi2_ct", "LR"),
                                ("diparm1", 'lnalpha, exp label("alpha")')], constant=pr.constant)
    est.dispersion = "mean"
    est.chibar = (chibar, p_chibar)
    est.footer = _nbreg_footer
    if pr.o.get("irr"):
        est.coef_title, est.eform_on = "IRR", True
    _display_std(s, est)
    est.coef_title, est.eform_on = "Coef.", False


def _bmap(pr, b, kx):
    bmap = {}
    j = 0
    for c in pr.smp.cols:
        if c.values is not None:
            bmap[c.name] = float(b[j])
            j += 1
        else:
            bmap[c.name] = 0.0
    if pr.constant:
        bmap["_cons"] = float(b[kx - 1])
    return bmap


def _semap(pr, V, kx):
    semap = {}
    j = 0
    for c in pr.smp.cols:
        if c.values is not None:
            semap[c.name] = float(np.sqrt(max(V[j, j], 0)))
            j += 1
        else:
            semap[c.name] = 0.0
    if pr.constant:
        semap["_cons"] = float(np.sqrt(max(V[kx - 1, kx - 1], 0)))
    return semap


def _nbreg_footer(s: "Session", est: Estimates) -> None:
    chibar, p = est.chibar
    s.output.write(f"Likelihood-ratio test of alpha=0:  chibar2(01) = {chibar:7.2f} Prob>=chibar2 = {p:5.3f}\n",
                   "text")


# ---------------------------------------------------------------------------
# tobit
# ---------------------------------------------------------------------------

def cmd_tobit(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "tobit", args)
        return
    pr = _Prep(s, args, {"ll": 2, "ul": 2}, weights=("aweight", "fweight", "iweight", "pweight"))
    out = s.output
    pr.write_notes(out)
    y = pr.y
    o = pr.o
    lo = hi = None
    if o.get("ll") is not None and o.get("ll") is not False:
        lo = float(y.min()) if o["ll"] is True or str(o["ll"]).strip() == "" else float(s.eval(str(o["ll"])))
    if o.get("ul") is not None and o.get("ul") is not False:
        hi = float(y.max()) if o["ul"] is True or str(o["ul"]).strip() == "" else float(s.eval(str(o["ul"])))
    if lo is None and hi is None:
        raise StataError(198, "must specify censoring point")   # VERIFICAR
    left = (y <= lo) if lo is not None else np.zeros(len(y), bool)
    right = (y >= hi) if hi is not None else np.zeros(len(y), bool)
    unc = ~(left | right)
    X = np.column_stack([pr.X, np.ones(len(y))]) if pr.constant else pr.X
    kx = X.shape[1]
    w = pr.w

    def make(Xm):
        k = Xm.shape[1]

        def ll_obs(b):
            xb = Xm @ b[:k]
            sg = np.exp(b[k])
            out_ = np.zeros(len(y))
            out_[unc] = st.norm.logpdf((y[unc] - xb[unc]) / sg) - np.log(sg)
            if lo is not None:
                out_[left] = st.norm.logcdf((lo - xb[left]) / sg)
            if hi is not None:
                out_[right] = st.norm.logsf((hi - xb[right]) / sg)
            return out_

        def scores(b):
            xb = Xm @ b[:k]
            sg = np.exp(b[k])
            sx = np.zeros(len(y))
            ss = np.zeros(len(y))
            z = (y - xb) / sg
            sx[unc] = z[unc] / sg
            ss[unc] = z[unc] ** 2 - 1
            if lo is not None:
                zl = (lo - xb[left]) / sg
                lam = np.exp(st.norm.logpdf(zl) - st.norm.logcdf(zl))
                sx[left] = -lam / sg
                ss[left] = -lam * zl
            if hi is not None:
                zh = (hi - xb[right]) / sg
                lam = np.exp(st.norm.logpdf(zh) - st.norm.logsf(zh))
                sx[right] = lam / sg
                ss[right] = lam * zh
            return np.column_stack([Xm * sx[:, None], ss])

        def grad(b):
            return scores(b).T @ w

        def fun(b, todo):
            ll = float(w @ ll_obs(b))
            if todo == 0:
                return ll, None, None
            return ll, grad(b), num_hess(grad, b)
        return fun, scores

    # início: MQO (VERIFICAR: valores iniciais do Stata)
    beta, *_ = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
    resid = y - X @ beta
    sig0 = np.sqrt(np.sum(w * resid ** 2) / w.sum())
    fun0, _ = make(np.ones((len(y), 1)))
    res0 = pr.maximize(fun0, np.array([np.average(y, weights=w), np.log(np.std(y) or 1.0)]))
    fun, scores = make(X)
    res = pr.maximize(fun, np.concatenate([beta, [np.log(sig0)]]))
    # o tobit do Stata 14 não mostra o log de iterações (compat 0510)
    Vl, n_clust = pr.vce(res.H, scores(res.b) if pr.robust else None)
    # e(b) do Stata 14 guarda sigma (não ln sigma): delta-método
    sg = float(np.exp(res.b[-1]))
    Jm = np.eye(kx + 1)
    Jm[-1, -1] = sg
    V = Jm @ Vl @ Jm.T
    b_rep = np.concatenate([res.b[:kx], [sg]])
    slope = list(range(kx - (1 if pr.constant else 0)))
    kind, chi2, pv = _model_test(pr, b_rep, V, slope, res.ll, res0.ll)
    r2p = 1 - res.ll / res0.ll
    names, _ = _slope_rows(pr, b_rep, V)
    ds = s.data
    rows = rows_from_columns(pr.smp.cols, _bmap(pr, b_rep, kx), _semap(pr, V, kx), ds,
                             constant=pr.constant)
    rows.append(CoefRow("sep"))
    c = crit("z", pr.smp.level, None)
    se_s = float(np.sqrt(V[-1, -1]))
    rows.append(CoefRow("aux", "/sigma", sg, se_s, name="_cons"))
    layout = [(nm, cc.values is not None) for nm, cc in zip(names, pr.smp.cols)]
    if pr.constant:
        layout.append(("_cons", True))
        names = names + ["_cons"]
    layout.append(("_cons", True))
    b_full, V_full = _expand(pr, b_rep, V, layout)
    eqs = ["model"] * len(names) + ["sigma"]
    est = _finish(s, pr, "tobit", args, title="Tobit regression", names=names + ["_cons"], eqnames=eqs,
                  b=b_full, V=V_full, rows=rows, res=res, ll0=res0.ll, df_m=len(slope), chi2kind=kind,
                  chi2=chi2, p=pv, r2p=r2p, display=_display_std, n_clust=n_clust, k_eq=2, k_aux=1,
                  extra_scalars=[("N_unc", float(w[unc].sum())), ("N_lc", float(w[left].sum())),
                                 ("N_rc", float(w[right].sum())), ("sigma", sg)]
                  + ([("llopt", lo)] if lo is not None else []) + ([("ulopt", hi)] if hi is not None else []),
                  constant=pr.constant)
    est.cens = (float(w[left].sum()), float(w[unc].sum()), float(w[right].sum()), lo, hi)
    # o tobit do Stata 14 usa t com N - k graus de liberdade (VERIFICAR)
    # graus de liberdade N - df_m (compat 0510: 60 obs, 2 inclinações -> 58)
    est.stat, est.df_r = "t", pr.N - len(slope)
    s.e["df_r"] = float(pr.N - len(slope))
    est.footer = _tobit_footer
    _display_std(s, est)


def _tobit_footer(s: "Session", est: Estimates) -> None:
    nl, nu, nr, lo, hi = est.cens
    out = s.output
    dep = est.depvar
    out.write(f"{int(nl):>14}  left-censored observations"
              + (f" at {dep} <= {g(lo, '%9.0g')}" if lo is not None and nl else "") + "\n", "text")
    out.write(f"{int(nu):>14}     uncensored observations\n", "text")
    out.write(f"{int(nr):>14} right-censored observations"
              + (f" at {dep} >= {g(hi, '%9.0g')}" if hi is not None and nr else "") + "\n", "text")
