"""Output and postestimation after sem: the header and table (also standardized), estat gof and
predict."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from .results import Estimates, crit, g9
from .sem import Model, _standardized

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")
optimize = lazy("scipy.optimize")

SYS = M.SYSMISS


def equality_constraints(model: Model) -> list[str]:
    """Constraint lines of the paths that share a label: [y1]x1 - [y2]x2 = 0 (first path minus
    each of the others; VERIFY)."""
    groups: dict[str, list[tuple[str, str]]] = {}
    for (kind, a, b), lab in model.ties.items():
        groups.setdefault(lab, []).append((a, b))
    out = []
    for lab, items in groups.items():
        a0, b0 = items[0]
        for a1, b1 in items[1:]:
            out.append(f"[{b0}]{a0} - [{b1}]{a1} = 0")
    return out


def display_sem(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    model: Model = est.model
    spec = model.spec
    lv = level or est.level
    if header:
        out.write("\n" + f"{'Structural equation model':<48}{'Number of obs':<18}= ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        out.write(f"Estimation method  = {x.get('method', 'ml')}\n", "text")
        if x.get("method") != "adf":
            out.write("Log likelihood     = ", "text")
            out.write(f"{g9(x['ll'], '%10.0g'):>10}\n", "result")
        else:
            out.write("Discrepancy        = ", "text")
            out.write(f"{g9(x['chi2_ms'] / x['N'], '%10.0g'):>10}\n", "result")
        # endogenous latent variables first (compat 1302: [y1]G before [x1]F)
        anchors = sorted(model.anchor.items(), key=lambda la: la[0] not in model.endog)
        fixed_paths = [(a, b, v) for (a, b), v in spec.fixed.items()]
        k = 0
        eqs_ = equality_constraints(model)
        if anchors or fixed_paths or eqs_:
            out.write("\n", "text")
        for item in eqs_:
            k += 1
            out.write(f" ({k:2d})  {item}\n", "text")
        for L, a in anchors:
            k += 1
            out.write(f" ({k:2d})  [{a}]{L} = 1\n", "text")
        for a, b, v in fixed_paths:
            k += 1
            out.write(f" ({k:2d})  [{b}]{a} = {v:g}\n", "text")
        if not k:
            out.write("\n", "text")
    if not table:
        return
    th = x["th"]
    V = x["V"]
    se = np.sqrt(np.maximum(np.diag(V), 0))
    std = x.get("standardized")
    fixed_std: dict[tuple[str, str], tuple[float, float]] = {}
    if std:
        fixed_keys = [(a, b) for a, b in spec.paths if ("path", a, b) not in
                      {(k_, a_, b_) for k_, a_, b_ in model.free}]

        def f(t):
            extra = [_standardized_fixed(model, t, a, b, spec.fixed.get((a, b), 1.0)) for a, b in fixed_keys]
            return np.concatenate([_standardized(model, t), extra])
        n0 = len(th)
        full = f(th)
        th_s = full[:n0]
        J = np.zeros((len(full), n0))
        for i in range(n0):
            e = np.zeros(n0)
            h = 1e-6 * (abs(th[i]) + 1)
            e[i] = h
            J[:, i] = (f(th + e) - f(th - e)) / (2 * h)
        se_all = np.sqrt(np.maximum(np.diag(J @ V @ J.T), 0))
        se = se_all[:n0]
        for j, key in enumerate(fixed_keys):
            fixed_std[key] = (full[n0 + j], se_all[n0 + j])
        th = th_s
    c = crit("z", lv, None)
    pos = {(kind, a, b): i for i, (kind, a, b) in enumerate(model.free)}

    def elab(n):
        return f"e.{n}" if n in model.endog else n
    # variances: errors of the observed endogenous, then of the latent endogenous, then the
    # exogenous variables (compat 1302)
    var_order = ([n for n in model.names if n in model.endog and n in spec.observed]
                 + [n for n in model.names if n in model.endog and n not in spec.observed]
                 + [n for n in model.names if n not in model.endog and n not in spec.observed]
                 + [n for n in model.names if n not in model.endog and n in spec.observed])
    var_order = [n for n in var_order if ("var", n, n) not in model.hidden]
    covs = [(kind, a, b) for kind, a, b in model.free if kind == "cov" and (kind, a, b) not in model.hidden]
    labels = [f"var({elab(n)})" for n in var_order] + [f"cov({elab(a)},{elab(b)})" for _, a, b in covs]
    # variance rows end at the bar; a long label widens the whole table (compat 1302)
    W = max([12] + [len(lb) - 1 for lb in labels])
    dash = "-" * (W + 1) + "+" + "-" * 64 + "\n"
    out.write("-" * (W + 66) + "\n", "text")
    if x.get("method") != "adf":
        out.write(f"{'':>{W}} |                 OIM\n", "text")
    out.write(f"{'Standardized' if std else '':>{W}} |      Coef.   Std. Err.      z    P>|z|     "
              f"[{lv:g}% Conf. Interval]\n", "text")
    out.write(dash, "text")

    def row(label, b, s_, kind="coef"):
        if kind in ("var", "cov"):
            out.write(f"{label:>{W + 1}}|", "text")
        else:
            out.write(f"{label:>{W}} |", "text")
        if kind == "constrained":
            out.write(f"  {g9(b):>9}  (constrained)\n", "result")
            return
        if s_ > 0 and np.isfinite(s_):
            if kind == "var":
                # variances: interval on the log scale
                lo, hi = b * math.exp(-c * s_ / b), b * math.exp(c * s_ / b)
                out.write(f"  {g9(b):>9}  {g9(s_):>9}" + " " * 17 + f"    {g9(lo):>9}   {g9(hi):>9}\n",
                          "result")
                return
            z = b / s_
            pv = 2 * float(st.norm.sf(abs(z)))
            out.write(f"  {g9(b):>9}  {g9(s_):>9}  {z:>7.2f}  {pv:>6.3f}    {g9(b - c * s_):>9}"
                      f"   {g9(b + c * s_):>9}\n", "result")
        elif kind == "var":
            # a variance fixed by the standardization (var(F) = 1)
            out.write(f"  {g9(b):>9}  {'.':>9}" + " " * 17 + f"    {'.':>9}   {'.':>9}\n", "result")
        else:
            out.write(f"  {g9(b):>9}  {'.':>9}\n", "result")

    eqs = [n for n in model.names if n in model.endog]
    meas = [n for n in eqs if n in spec.observed and all(a in spec.latent for a, b in spec.paths if b == n)]
    struct = [n for n in eqs if n not in meas]
    for section, members in (("Structural", struct), ("Measurement", meas)):
        if not members:
            continue
        out.write(f"{section:<{W}} |\n", "text")
        for j, y in enumerate(members):
            if j:
                out.write("  " + "-" * (W - 1) + "+" + "-" * 64 + "\n", "text")
            out.write(f"  {y + ' <-':<{W - 2}} |\n", "text")
            for a, b in spec.paths:
                if b != y:
                    continue
                key = ("path", a, b)
                if key in pos:
                    row(a, th[pos[key]], se[pos[key]])
                elif std:
                    row(a, *fixed_std[(a, b)])
                else:
                    row(a, spec.fixed.get((a, b), 1.0), 0, "constrained")
            if y in spec.observed:
                k = pos[("mean", y, y)]
                row("_cons", th[k], se[k])
        out.write(dash, "text")
    means = [n for n in model.exog if n in spec.observed and ("mean", n, n) not in model.hidden]
    if means:
        for n in means:
            k = pos[("mean", n, n)]
            row(f"mean({n})", th[k], se[k])
        out.write(dash, "text")
    for n in var_order:
        key = ("var", n, n)
        lab = f"var({elab(n)})"
        if key in pos:
            row(lab, th[pos[key]], se[pos[key]], "var")
        elif n in model.var_fixed:
            row(lab, model.var_fixed[n], 0, "constrained")
    if covs:
        out.write(dash, "text")
        for key in covs:
            _, a, b = key
            row(f"cov({elab(a)},{elab(b)})", th[pos[key]], se[pos[key]], "cov")
    out.write("-" * (W + 66) + "\n", "text")
    df = x["df"]
    pv = "     ." if df == 0 else f"{x['p_ms']:>6.4f}"
    lab = f"chi2({df})"
    if x.get("method") == "adf":
        out.write(f"Discr. test of model vs. saturated: {lab:<10}= {x['chi2_ms']:>8.2f}, Prob > chi2 = {pv}\n",
                  "text")
    else:
        out.write(f"LR test of model vs. saturated: {lab:<10}= {x['chi2_ms']:>9.2f}, Prob > chi2 = {pv}\n",
                  "text")


def _standardized_fixed(model, th, a, b, v):
    *_, (A, S, m, B, Sv, mv) = model.implied(th)
    ia, ib = model.idx[a], model.idx[b]
    return v * math.sqrt(Sv[ia, ia]) / math.sqrt(Sv[ib, ib])


def estat_gof(s: "Session", est: Estimates, args: str) -> None:
    """estat gof [, stats(...)]: likelihood-ratio tests and, with stats(all), RMSEA, AIC, BIC,
    CFI, TLI, SRMR and CD ([SEM] estat gof; formulas VERIFY)."""
    p = parse_standard(args)
    o = match_options(p.options, {"stats": 4})
    want = set(str(o.get("stats") or "chi2").lower().split())
    if "all" in want:
        want = {"chi2", "rmsea", "ic", "indices", "residuals"}
    x = est.extra
    model: Model = est.model
    N = x["N"]
    prob = x.get("prob")
    G = len(prob.data) if prob is not None else 1
    adf = x.get("method") == "adf"
    if adf:
        want.discard("ic")
    out = s.output
    out.write("\n" + "-" * 76 + "\n", "text")
    out.write(f"{'Fit statistic':<21}|      Value   Description\n", "text")
    sep = "-" * 21 + "+" + "-" * 54 + "\n"
    out.write(sep, "text")
    r = {}

    def line(label, value, desc=""):
        out.write(f"{label:>20} |", "text")
        out.write(f"{value:>11}", "result")
        out.write(("   " + desc if desc else "") + "\n", "text")

    blocks = []
    if "chi2" in want:
        blocks.append("chi2")
    if "rmsea" in want:
        blocks.append("rmsea")
    if "ic" in want:
        blocks.append("ic")
    if "indices" in want:
        blocks.append("indices")
    if "residuals" in want:
        blocks.append("residuals")
    df, chi, dfb, chib = x["df"], x["chi2_ms"], x["df_bs"], x["chi2_bs"]
    for i, b in enumerate(blocks):
        if i:
            out.write(sep, "text")
        if b == "chi2":
            out.write(f"{'Likelihood ratio':<21}|\n", "text")
            line(f"chi2_ms({df})", f"{chi:.3f}", "model vs. saturated")
            line("p > chi2", f"{x['p_ms']:.3f}" if df > 0 else ".")
            line(f"chi2_bs({dfb})", f"{chib:.3f}", "baseline vs. saturated")
            line("p > chi2", f"{x['p_bs']:.3f}" if dfb > 0 else ".")
            r.update({"chi2_ms": chi, "df_ms": float(df), "p_ms": x["p_ms"], "chi2_bs": chib,
                      "df_bs": float(dfb), "p_bs": x["p_bs"]})
        elif b == "rmsea":
            # with several groups the formula includes G; the Browne-Cudeck p-value is not computed
            rm = math.sqrt(G * max(chi - df, 0) / (df * N)) if df > 0 else 0.0
            lo, hi = _rmsea_ci(chi, df, N / G)
            pclose = float(st.ncx2.sf(chi, df, 0.05 ** 2 * (N / G) * df)) if df > 0 else 1.0
            out.write(f"{'Population error':<21}|\n", "text")
            line("RMSEA", f"{rm:.3f}", "Root mean squared error of approximation")
            line("90% CI, lower bound", f"{lo:.3f}")
            line("upper bound", f"{hi:.3f}")
            line("pclose", f"{pclose:.3f}" if pclose < SYS else ".", "Probability RMSEA <= 0.05")
            r.update({"rmsea": rm, "lb90_rmsea": lo, "ub90_rmsea": hi, "pclose": pclose})
        elif b == "ic":
            k = prob.nglob if prob is not None else model.nfree
            aic = -2 * x["ll"] + 2 * k
            bic = -2 * x["ll"] + k * math.log(N)
            out.write(f"{'Information criteria':<21}|\n", "text")
            line("AIC", f"{aic:.3f}", "Akaike's information criterion")
            line("BIC", f"{bic:.3f}", "Bayesian information criterion")
            r.update({"aic": aic, "bic": bic})
        elif b == "indices":
            num = max(chi - df, 0)
            den = max(chib - dfb, chi - df, 0)
            cfi = 1 - num / den if den > 0 else 1.0
            tli = ((chib / dfb - chi / df) / (chib / dfb - 1)) if df > 0 and dfb > 0 and chib / dfb != 1 else 1.0
            out.write(f"{'Baseline comparison':<21}|\n", "text")
            line("CFI", f"{cfi:.3f}", "Comparative fit index")
            line("TLI", f"{tli:.3f}", "Tucker-Lewis index")
            r.update({"cfi": cfi, "tli": tli})
        else:
            srmr, cd = _srmr_cd(est)
            if G > 1:
                from .semestat import _equation_fit, _srmr, _units
                us = _units(est)
                srmr = math.sqrt(sum(_srmr(u["model"], u["th"], u["data"]) ** 2 for u in us) / G)
                cd = float(np.mean([_equation_fit(u["model"], u["th"])[1] for u in us]))   # VERIFY
            out.write(f"{'Size of residuals':<21}|\n", "text")
            line("SRMR", f"{srmr:.3f}", "Standardized root mean squared residual")
            line("CD", f"{cd:.3f}", "Coefficient of determination")
            r.update({"srmr": srmr, "cd": cd})
    out.write("-" * 76 + "\n", "text")
    s.r = r


def _rmsea_ci(chi, df, N):
    if df <= 0:
        return 0.0, 0.0

    def nc_for(p):
        f = lambda lam: float(st.ncx2.cdf(chi, df, lam)) - p   # noqa: E731
        if f(0) < 0:
            return 0.0
        hi = max(chi, 1.0)
        while f(hi) > 0:
            hi *= 2
        return float(optimize.brentq(f, 0, hi))
    lo_nc, hi_nc = nc_for(0.95), nc_for(0.05)
    return math.sqrt(lo_nc / (df * N)), math.sqrt(hi_nc / (df * N))


def _srmr_cd(est):
    model: Model = est.model
    x = est.extra
    mu, Sig, (A, S, m, B, Sv, mv) = model.implied(est.extra["th"])
    cov = x["cov"]
    p = model.p
    d = np.sqrt(np.diag(cov))
    ds = np.sqrt(np.diag(Sig))
    # sample minus implied correlations, and the standardized mean residuals, over the
    # p(p+3)/2 moments (compat 1302)
    R = cov / np.outer(d, d) - Sig / np.outer(ds, ds)
    tri = R[np.tril_indices(p)]
    dm = (x["mean"] - mu) / d
    srmr = float(np.sqrt(((tri ** 2).sum() + (dm ** 2).sum()) / (p * (p + 3) / 2)))
    endog = [model.idx[n] for n in model.endog]
    if endog:
        Psi = S[np.ix_(endog, endog)]
        Sy = Sv[np.ix_(endog, endog)]
        cd = 1 - float(np.linalg.det(Psi) / np.linalg.det(Sy))
    else:
        cd = SYS
    return srmr, cd


def _predict(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    """predict newvar, xb(depvar) and predict newvar, latent(name): linear prediction of an
    observed endogenous variable, and factor scores by the regression method."""
    model: Model = est.model
    spec = model.spec
    ds = s.data
    target = (getattr(est, "_pred_opts", None) or {}).get(stat)
    X = np.column_stack([ds.get(n).data.astype(float) for n in spec.observed])
    ok = np.all(X < SYS, axis=1)
    mu, Sig, (A, S, m, B, Sv, mv) = model.implied(est.extra["th"])
    if stat == "latent":
        names = [target] if isinstance(target, str) and target.strip() else spec.latent[:1]
        L = names[0].strip()
        if L not in spec.latent:
            raise StataError(198, f"{L} is not a latent variable")   # VERIFY
        il = model.idx[L]
        p = model.p
        cov_lx = Sv[il, :p]
        w = np.linalg.solve(Sig, cov_lx)
        vals = np.where(ok, (np.where(ok[:, None], X, 0.0) - mu) @ w + mv[il], SYS)
        return vals, ""
    if stat in ("xb", ""):
        y = target.strip() if isinstance(target, str) and target.strip() else (
            [n for n in model.endog if n in spec.observed] or [""])[0]
        if y not in spec.observed or y not in model.endog:
            raise StataError(198, f"{y} is not an observed endogenous variable")   # VERIFY
        iy = model.idx[y]
        # E(y | its predictors): direct paths times the predictors (latent ones at their
        # predicted values)
        val = np.full(ds.nobs, m[iy])
        for a, b in spec.paths:
            if b != y:
                continue
            ia = model.idx[a]
            coef = A[iy, ia]
            if a in spec.observed:
                val = val + coef * np.where(ok, X[:, spec.observed.index(a)], 0.0)
            else:
                lv, _ = _predict_latent(model, est, X, ok, a)
                val = val + coef * lv
        return np.where(ok, val, SYS), ""
    raise StataError(198, f"option {stat} not allowed")


def _predict_latent(model, est, X, ok, L):
    mu, Sig, (A, S, m, B, Sv, mv) = model.implied(est.extra["th"])
    il = model.idx[L]
    w = np.linalg.solve(Sig, Sv[il, :model.p])
    return (np.where(ok[:, None], X, 0.0) - mu) @ w + mv[il], None
