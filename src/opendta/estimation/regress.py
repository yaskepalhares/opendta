"""regress ([R] regress): mínimos quadrados com pesos, vce(robust) e
vce(cluster), tabela ANOVA e tabela de coeficientes do Stata 14."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from . import fvars
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def sig5(x: float) -> str:
    """Root MSE: cinco algarismos significativos, sem o zero à esquerda."""
    if not np.isfinite(x) or abs(x) >= SYS:
        return "."
    t = f"{x:.5g}"
    if "e" in t:
        return format_value(x, "%9.0g", pad=False).strip()
    if t.startswith("0."):
        t = t[1:]
    elif t.startswith("-0."):
        t = "-" + t[2:]
    return t


def g(x: float, fmt: str) -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def comma(n: float) -> str:
    return f"{int(round(n)):,}"


# ---------------------------------------------------------------------------
# preparação comum aos comandos de estimação
# ---------------------------------------------------------------------------

class Sample:
    """Dependente, colunas e amostra de um comando de estimação."""

    def __init__(self, s: "Session", args: str, *, opts_spec: dict[str, int], weights=("aweight",
                 "fweight", "iweight", "pweight"), need_dep: bool = True, constant_opt: bool = True,
                 extra_vars: str = ""):
        with fvars.caching():
            self._init(s, args, opts_spec=opts_spec, weights=weights, need_dep=need_dep,
                       extra_vars=extra_vars)

    def _init(self, s: "Session", args: str, *, opts_spec: dict[str, int], weights,
              need_dep: bool, extra_vars: str = "") -> None:
        from ..commands._util import touse
        from ..commands.summarize import weights as get_weights
        self.s = s
        p = parse_standard(args)
        self.p = p
        spec = dict(opts_spec)
        spec.update({"noconstant": 4, "level": 1, "vce": 3, "robust": 1, "cluster": 2})
        self.o = match_options(p.options, spec) if p.options.strip() else {}
        ds = s.data
        toks = fvars._split_tokens(p.varlist)
        if need_dep and not toks:
            raise StataError(100, "varlist required")
        self.depname = ""
        dep_terms = []
        if need_dep:
            dep_terms = fvars.expand_fv(ds, toks[0])
            if len(dep_terms) != 1 or dep_terms[0].is_interaction or dep_terms[0].comps[0].factor:
                raise StataError(198, "depvar may not be a factor variable")   # VERIFICAR
            self.depcomp = dep_terms[0].comps[0]
            self.depname = self.depcomp.tsname
        rest = " ".join(toks[1:] if need_dep else toks)
        self.terms = fvars.expand_fv(ds, rest) if rest.strip() else []
        mask = touse(s, p)
        extra = fvars.expand_fv(ds, extra_vars) if extra_vars.strip() else []
        bad = fvars.term_missing(ds, dep_terms + self.terms + extra)
        mask &= ~bad
        self.w, self.wtype, mask = get_weights(s, p, mask, weights)
        self.mask = mask
        self.constant = not self.o.get("noconstant")
        self.level = float(self.o["level"]) if self.o.get("level") else float(s.settings.get("level", "95"))
        if self.level < 10 or self.level > 99.99:
            raise StataError(198, "level() must be between 10 and 99.99 inclusive")
        # vce
        vce = str(self.o.get("vce", "") or "").strip()
        self.vce = "ols"
        self.clustvar = ""
        if self.o.get("robust"):
            self.vce = "robust"
        if self.o.get("cluster"):
            self.vce, self.clustvar = "cluster", str(self.o["cluster"]).strip()
        if vce:
            w = vce.split()
            kind = w[0].lower()
            if "robust".startswith(kind) and len(kind) >= 1:
                self.vce = "robust"
            elif "cluster".startswith(kind) and len(kind) >= 2:
                if len(w) < 2:
                    raise StataError(198, "option vce() misspecified")   # VERIFICAR
                self.vce, self.clustvar = "cluster", w[1]
            elif kind in ("ols", "oim", "opg"):
                self.vce = kind
            elif kind in ("hc2", "hc3"):
                self.vce = kind
            else:
                raise StataError(198, f"vcetype '{w[0]}' not allowed")   # VERIFICAR
        if self.wtype == "pweight" and self.vce in ("ols", "oim"):
            self.vce = "robust"
        if self.clustvar:
            from ..core.varlist import resolve_name
            self.clustvar = resolve_name(ds, self.clustvar)
            cv = ds.get(self.clustvar).data
            if ds.get(self.clustvar).is_string:
                codes = {v: i for i, v in enumerate(sorted(set(cv[self.mask])))}
                self.cluster = np.array([codes[v] for v in cv[self.mask]], dtype=float)
            else:
                cm = cv >= SYS
                if np.any(cm & self.mask):
                    self.mask &= ~cm
                    # VERIFICAR: o Stata exclui observações com cluster missing
                self.cluster = cv[self.mask]
        mask = self.mask
        n = int(mask.sum())
        if n == 0:
            raise StataError(2000, "no observations")
        self.y = fvars.comp_values(ds, self.depcomp)[mask] if need_dep else None
        self.cols = fvars.build_design(ds, self.terms, mask, constant=self.constant)
        self.n = int(mask.sum())

    def weights_for_fit(self) -> tuple[np.ndarray, float]:
        """(pesos usados nas contas, N reportado)."""
        n = self.n
        if self.w is None:
            return np.ones(n), float(n)
        w = self.w[self.mask].astype(np.float64)
        if self.wtype == "fweight":
            return w, float(w.sum())
        if self.wtype in ("aweight", "pweight"):
            return w * n / w.sum(), float(n)
        return w, float(w.sum())    # iweight: N = soma dos pesos (compat 0505)


def cmdline(name: str, args: str) -> str:
    return f"{name} {args.strip()}".strip()


def wexp_macros(smp: Sample) -> list[tuple[str, str]]:
    if smp.p.weight is None:
        return []
    return [("wtype", smp.wtype), ("wexp", f"= {smp.p.weight[1]}")]


# ---------------------------------------------------------------------------
# regress
# ---------------------------------------------------------------------------

def cmd_regress(s: "Session", args: str) -> None:
    if not args.strip() or args.strip().startswith(","):
        from .postest import replay
        replay(s, "regress", args)
        return
    smp = Sample(s, args, opts_spec={"beta": 4, "noheader": 6, "notable": 5, "hascons": 7,
                                     "tsscons": 7, "depname": 4, "mse1": 4, "plus": 4})
    o = smp.o
    ds = s.data
    w, N = smp.weights_for_fit()
    out = s.output
    if smp.wtype in ("aweight", "pweight"):
        tot = float(smp.w[smp.mask].sum())
        out.write(f"(sum of wgt is {format_value(tot, '%12.4e')})\n", "text")
    notes = fvars.drop_collinear(smp.cols, smp.constant, w)
    for nm in notes:
        out.write(f"note: {nm} omitted because of collinearity\n", "text")
    est_cols = [c for c in smp.cols if c.values is not None]
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(smp.n)] if smp.constant else [])) \
        if (est_cols or smp.constant) else np.zeros((smp.n, 0))
    y = smp.y
    k = X.shape[1]
    beta, XtWXi = _ols(X, y, w, smp.constant)
    yhat = X @ beta if k else np.zeros(smp.n)
    resid = y - yhat
    rss = float(np.sum(w * resid ** 2))
    if smp.constant:
        ybar = float(np.sum(w * y) / np.sum(w))
        tss = float(np.sum(w * (y - ybar) ** 2))
    else:
        tss = float(np.sum(w * y ** 2))
    mss = tss - rss
    rank = k
    df_m = rank - (1 if smp.constant else 0)
    df_r = N - rank
    if df_r <= 0:
        raise StataError(2001, "insufficient observations")
    s2 = rss / df_r
    r2 = mss / tss if tss > 0 else SYS
    r2_a = 1 - (1 - r2) * ((N - 1) if smp.constant else N) / df_r if tss > 0 else SYS
    rmse = np.sqrt(s2)
    Vb = s2 * XtWXi
    vcetype = ""
    n_clust = None
    if smp.vce in ("robust", "hc2", "hc3", "cluster"):
        wr = w if smp.wtype != "fweight" else w
        if smp.vce == "cluster":
            ids = smp.cluster
            uniq = np.unique(ids)
            n_clust = len(uniq)
            meat = np.zeros((k, k))
            sc = X * (wr * resid)[:, None]
            for u in uniq:
                sg = sc[ids == u].sum(axis=0)
                meat += np.outer(sg, sg)
            q = (n_clust / (n_clust - 1)) * ((N - 1) / (N - k))
            Vb = q * XtWXi @ meat @ XtWXi
            df_r = n_clust - 1
        else:
            if smp.vce == "robust":
                sc = X * (wr * resid)[:, None]
                if smp.wtype == "fweight":
                    meat = (X * (w * resid ** 2)[:, None]).T @ X
                else:
                    meat = sc.T @ sc
                Vb = (N / (N - k)) * XtWXi @ meat @ XtWXi
            else:
                h = np.einsum("ij,jk,ik->i", X, XtWXi, X) * w
                e2 = resid ** 2 / ((1 - h) if smp.vce == "hc2" else (1 - h) ** 2)
                meat = (X * (w * w * e2)[:, None]).T @ X
                Vb = XtWXi @ meat @ XtWXi
        vcetype = {"robust": "Robust", "cluster": "Robust", "hc2": "Robust HC2",
                   "hc3": "Robust HC3"}[smp.vce]
    Vb = (Vb + Vb.T) / 2
    # F
    slopes = list(range(len(est_cols)))
    F = SYS
    df_m_robust = len(slopes)
    if smp.vce == "ols":
        F = (mss / df_m) / s2 if df_m > 0 else 0.0
    elif smp.vce == "cluster" and np.linalg.matrix_rank(Vb, tol=1e-12 * max(np.abs(Vb).max(), 1e-300)) < k:
        # poucos clusters: V tem posto menor que k; o Stata mostra F(posto-1, G-1)
        # com valor missing (compat 0501). VERIFICAR o caso geral.
        rk = np.linalg.matrix_rank(Vb, tol=1e-12 * max(np.abs(Vb).max(), 1e-300))
        df_m_robust = rk - (1 if smp.constant else 0)
        F = SYS
    elif slopes:
        bb = beta[slopes]
        VV = Vb[np.ix_(slopes, slopes)]
        try:
            F = float(bb @ np.linalg.solve(VV, bb)) / len(slopes)
        except np.linalg.LinAlgError:
            F = SYS
    ll = -0.5 * N * (np.log(2 * np.pi) + np.log(rss / N) + 1) if rss > 0 else SYS
    ll_0 = -0.5 * N * (np.log(2 * np.pi) + np.log(tss / N) + 1) if tss > 0 and smp.constant else SYS
    # resultados por nome (inclui bases e omitidas com zero)
    names = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols]
    if smp.constant:
        names.append("_cons")
    kk = len(names)
    b_full = np.zeros(kk)
    V_full = np.zeros((kk, kk))
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([kk - 1] if smp.constant else [])
    b_full[pos] = beta
    V_full[np.ix_(pos, pos)] = Vb
    bmap = {}
    semap = {}
    for i, c in enumerate(smp.cols):
        bmap[c.name] = b_full[i]
        semap[c.name] = np.sqrt(V_full[i, i]) if V_full[i, i] > 0 else 0.0
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b_full[-1], np.sqrt(max(V_full[-1, -1], 0))
    df_m_report = df_m if smp.vce == "ols" else df_m_robust
    est = Estimates("regress", smp.depname, names, b_full, V_full, smp.n, "t", df_r)
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    est.cols = smp.cols
    est.constant = smp.constant
    est.terms = smp.terms
    est.depcomp = smp.depcomp
    est.level = smp.level
    est.vcetype = vcetype
    est.clustvar = smp.clustvar
    est.n_clust = n_clust
    est.extra = {"N": N, "df_m": df_m_report, "df_r": df_r, "F": F, "r2": r2, "rmse": rmse,
                 "mss": mss, "rss": rss, "r2_a": r2_a, "ll": ll, "ll_0": ll_0, "rank": rank,
                 "tss": tss, "anova_df_m": df_m}
    est.wtype = smp.wtype
    est.resid_df = df_r
    est.s2 = s2
    est.XtWXi_full = None
    scalars = [("N", N), ("df_m", df_m_report), ("df_r", df_r), ("F", F), ("r2", r2), ("rmse", rmse),
               ("mss", mss), ("rss", rss), ("r2_a", r2_a), ("ll", ll), ("ll_0", ll_0), ("rank", rank)]
    if n_clust is not None:
        scalars.insert(1, ("N_clust", n_clust))   # VERIFICAR posição
    macros = [("cmdline", cmdline("regress", args)), ("title", "Linear regression"),
              ("marginsok", "XB default"), ("vce", smp.vce if smp.vce != "hc2" else "hc2")]
    if vcetype:
        macros.append(("vcetype", vcetype))
    if smp.clustvar:
        macros.append(("clustvar", smp.clustvar))
    macros += wexp_macros(smp)
    macros += [("depvar", smp.depname), ("cmd", "regress"), ("properties", "b V"),
               ("predict", "regres_p"), ("model", "ols"), ("estat_cmd", "regress_estat")]
    post(s, est, scalars, macros, smp.mask)
    if o.get("beta"):
        est.beta = _std_coefs(smp, w, bmap)
    display_regress(s, est, header=not o.get("noheader"), table=not o.get("notable"),
                    beta=bool(o.get("beta")))


def _ols(X: np.ndarray, y: np.ndarray, w: np.ndarray, constant: bool):
    """Coeficientes e (X'WX)⁻¹ (ordem: inclinações, depois a constante).

    set numerics stata: como o regress do Stata, produtos cruzados dos
    desvios em relação às médias (a constante sai de ȳ - x̄'b), invertidos
    sem truncar autovalores pequenos.
    set numerics precise: decomposição QR de W^(1/2)X (erro proporcional a
    κ(X), não a κ(X)²)."""
    from .numerics import precise
    k = X.shape[1]
    if k == 0:
        return np.zeros(0), np.zeros((0, 0))
    if precise():
        sw = np.sqrt(w)
        Q, R = np.linalg.qr(X * sw[:, None])
        beta = np.linalg.solve(R, Q.T @ (sw * y))
        Rinv = np.linalg.solve(R, np.eye(k))
        return beta, Rinv @ Rinv.T
    W = float(w.sum())
    if constant:
        Xs = X[:, :-1]
        m = (w @ Xs) / W
        ybar = float(w @ y) / W
        Xc = Xs - m
        yc = y - ybar
        if Xs.shape[1]:
            A = np.linalg.inv((Xc * w[:, None]).T @ Xc)
            b = A @ (Xc.T @ (w * yc))
        else:
            A = np.zeros((0, 0))
            b = np.zeros(0)
        cons = ybar - float(m @ b)
        full = np.zeros((k, k))
        full[:-1, :-1] = A
        full[:-1, -1] = full[-1, :-1] = -(A @ m)
        full[-1, -1] = 1.0 / W + float(m @ A @ m)
        return np.append(b, cons), full
    A = np.linalg.inv((X * w[:, None]).T @ X)
    return A @ (X.T @ (w * y)), A


def _std_coefs(smp: Sample, w: np.ndarray, bmap: dict) -> dict:
    def sd(x):
        m = np.sum(w * x) / np.sum(w)
        return np.sqrt(np.sum(w * (x - m) ** 2) / (np.sum(w) - 1))
    sy = sd(smp.y)
    return {c.name: bmap[c.name] * sd(c.values) / sy for c in smp.cols if c.values is not None}


def display_regress(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                    level: float | None = None, beta: bool = False) -> None:
    out = s.output
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        out.write("\n", "text")
        if not est.vcetype:
            _anova(s, est)
        else:
            _robust_header(s, est)
        out.write("\n", "text")
    if table:
        if est.vcetype == "Robust" and est.clustvar:
            # VERIFICAR texto e alinhamento
            out.write(f"{'(Std. Err. adjusted for ' + format(int(est.n_clust), ',') + ' clusters in ' + est.clustvar + ')':>78}\n",
                      "text")
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=lev, vcetype=est.vcetype)


def _anova(s: "Session", est: Estimates) -> None:
    out = s.output
    x = est.extra
    df_m, df_r = x["anova_df_m"], x["df_r"]
    N = x["N"]
    msm = x["mss"] / df_m if df_m > 0 else SYS
    msr = x["rss"] / df_r if df_r > 0 else SYS
    mst = x["tss"] / (N - (1 if est.constant else 0))
    right = [
        ("Number of obs", comma(N)),
        (f"F({int(df_m)}, {int(df_r)})", g(x["F"], "%9.2f")),
        ("Prob > F", g(st.f.sf(x["F"], df_m, df_r) if x["F"] < SYS else SYS, "%9.4f")),
        ("R-squared", g(x["r2"], "%9.4f")),
        ("Adj R-squared", g(x["r2_a"], "%9.4f")),
        ("Root MSE", sig5(x["rmse"])),
    ]
    left = [
        f"{'Source':>12} |       SS           df       MS   ",
        "-" * 13 + "+" + "-" * 34,
        _ss_row("Model", x["mss"], df_m, msm),
        _ss_row("Residual", x["rss"], df_r, msr),
        "-" * 13 + "+" + "-" * 34,
        _ss_row("Total", x["tss"], N - (1 if est.constant else 0), mst),
    ]
    for lft, (lab, val) in zip(left, right):
        out.write(f"{lft:<48}   ", "text")
        out.write(f"{lab:<16}= ", "text")
        out.write(f"{val:>9}\n", "result")


def _ss_row(name: str, ss: float, df: float, ms: float) -> str:
    return f"{name:>12} | {g(ss, '%11.0g'):>11}{int(df):>10}{g(ms, '%11.0g'):>12}"


def _robust_header(s: "Session", est: Estimates) -> None:
    out = s.output
    x = est.extra
    df_m, df_r = x["df_m"], x["df_r"]
    right = [
        ("Number of obs", comma(x["N"])),
        (f"F({int(df_m)}, {int(df_r)})", g(x["F"], "%10.2f")),
        ("Prob > F", g(st.f.sf(x["F"], df_m, df_r) if x["F"] < SYS else SYS, "%10.4f")),
        ("R-squared", g(x["r2"], "%10.4f")),
        ("Root MSE", sig5(x["rmse"])),
    ]
    for i, (lab, val) in enumerate(right):
        left = "Linear regression" if i == 0 else ""
        out.write(f"{left:<48}{lab:<18}= ", "text")
        out.write(f"{val:>10}\n", "result")
