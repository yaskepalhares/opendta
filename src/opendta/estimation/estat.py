"""estat: ic, vce, summarize, hettest, ovtest, classification and gof
([R] estat, [R] regress postestimation, [R] logistic postestimation).

Formulas from the manuals; layouts checked against the Stata 14 logs
(compat/do/0802_anova_estat.do); what could not be checked has the VERIFY mark."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core.stats import chi2_tail as _chi2_tail

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import find_top, match_options
from .results import current

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or x >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def cmd_estat(s: "Session", args: str) -> None:
    words = args.strip().split(None, 1)
    if not words:
        raise StataError(321, "estat requires a subcommand")   # VERIFY
    sub = words[0].rstrip(",")
    rest = words[1] if len(words) > 1 else ""
    if "," in words[0]:
        rest = words[0][words[0].index(","):] + " " + rest
    est = current(s)
    for name, (minimum, fn) in (getattr(est, "estat", None) or {}).items():
        if len(sub) >= minimum and name.startswith(sub):
            fn(s, est, rest)
            return
    table = [("ic", 2), ("vce", 3), ("summarize", 3), ("hettest", 3), ("ovtest", 3),
             ("classification", 5), ("gof", 3), ("esample", 7), ("bootstrap", 4), ("effects", 3)]
    full = next((n for n, m in table if len(sub) >= m and n.startswith(sub)), None)
    if full is None:
        raise StataError(321, f"invalid subcommand {sub}")   # VERIFY
    fn = {"ic": _ic, "vce": _vce, "summarize": _summarize, "hettest": _hettest, "ovtest": _ovtest,
          "classification": _classification, "gof": _gof, "bootstrap": _bootstrap, "effects": _effects}.get(full)
    if fn is None:
        raise StataError(321, f"invalid subcommand {sub}")
    fn(s, est, rest)


def _opts(rest: str, spec: dict) -> tuple[str, dict]:
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    return body.strip(), (match_options(opt, spec) if opt.strip() else {})


def _e(s: "Session", k: str) -> float:
    v = s.e.get(k)
    return float(v) if isinstance(v, (int, float)) else SYS


# ---------------------------------------------------------------------------
# ic
# ---------------------------------------------------------------------------

def _ic(s: "Session", est, rest: str) -> None:
    _, o = _opts(rest, {"n": 1})
    ll = _e(s, "ll")
    if ll >= SYS:
        raise StataError(321, "estat ic not valid")   # VERIFY
    N = float(o["n"]) if isinstance(o.get("n"), str) else _e(s, "N")
    ll0 = _e(s, "ll_0")
    df = _e(s, "rank")
    if df >= SYS:
        df = float(np.count_nonzero(np.diag(est.V)))
    aic = -2 * ll + 2 * df
    bic = -2 * ll + math.log(N) * df
    out = s.output
    out.write("\nAkaike's information criterion and Bayesian information criterion\n\n", "text")
    out.write("-" * 77 + "\n", "text")
    out.write(f"       Model |{'Obs':>11}{'ll(null)':>10}{'ll(model)':>11}{'df':>8}{'AIC':>12}{'BIC':>11}\n", "text")
    out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
    out.write(f"{'.':>12} |", "text")
    out.write(f"{N:>11,.0f}{_g(ll0):>10}{_g(ll):>11}{df:>8.0f}{_g(aic):>12}{_g(bic):>11}\n", "result")
    out.write("-" * 77 + "\n", "text")
    out.write(" " * 15 + "Note: N=Obs used in calculating BIC; see [R] BIC note.\n", "text")
    from ..commands.matrix import Matrix
    S = Matrix(np.array([[N, ll0, ll, df, aic, bic]]), ["."], ["N", "ll0", "ll", "df", "AIC", "BIC"])
    s.r = {"S": S}


# ---------------------------------------------------------------------------
# vce
# ---------------------------------------------------------------------------

def _vce(s: "Session", est, rest: str) -> None:
    _, o = _opts(rest, {"covariance": 3, "correlation": 3, "format": 3, "block": 3, "diag": 3})
    V = np.asarray(est.V, dtype=float)
    names = [n for n in est.names]
    corr = bool(o.get("correlation"))
    M_ = V.copy()
    if corr:
        d = np.sqrt(np.diag(V))
        with np.errstate(invalid="ignore", divide="ignore"):
            M_ = V / np.outer(d, d)
    fmt = o["format"].strip() if isinstance(o.get("format"), str) else ("%8.4f" if corr else "%10.0g")
    title = "Correlation matrix" if corr else "Covariance matrix"
    out = s.output
    out.write(f"\n{title} of coefficients of {est.cmd} model\n", "text")
    k = len(names)
    from ..lang.functions import _abbrev
    eqs = est.eqnames or [""] * k
    # column width = format width + 1 (compat 0802: %10.0g -> 11, %8.4f -> 9), plus one space
    # between columns
    import re as _re
    mw = _re.match(r"%-?(\d+)", fmt)
    cw = (int(mw.group(1)) if mw else 10) + 1
    per = max(1, 66 // (cw + 1))
    for start in range(0, k, per):
        cols = list(range(start, min(k, start + per)))
        out.write("\n", "text")
        head = "e(V)"
        out.write(f"{head:>12} |" + "".join(f"{_abbrev(_short(names[j]), cw - 1):>{cw}} " for j in cols) + "\n",
                  "text")
        out.write("-" * 13 + "+" + "-" * ((cw + 1) * len(cols)) + "\n", "text")
        last_eq = None
        for i in range(start, k):
            if eqs[i] and eqs[i] != "_" and eqs[i] != last_eq and len(set(eqs)) > 1:
                out.write(f"{_abbrev(eqs[i], 12):<12} |\n", "result")
                last_eq = eqs[i]
            out.write(f"{_abbrev(_short(names[i]), 12):>12} |", "text")
            cells = []
            for j in cols:
                if j > i:
                    cells.append(" " * (cw + 1))
                else:
                    v = M_[i, j]
                    t = format_value(float(v), fmt, pad=False).strip() if np.isfinite(v) else "."
                    cells.append(f"{t:>{cw}} ")
            out.write("".join(cells).rstrip() + "\n", "result")
    from ..commands.matrix import Matrix
    s.r = {"V": Matrix(M_.copy(), list(names), list(names))}


def _short(n: str) -> str:
    return n


# ---------------------------------------------------------------------------
# summarize
# ---------------------------------------------------------------------------

def _summarize(s: "Session", est, rest: str) -> None:
    from .postest import esample
    body, _ = _opts(rest, {"labels": 3, "noheader": 5, "noweights": 5, "equation": 2})
    mask = esample(s)
    names = [est.depvar] if est.depvar else []
    for n in est.names:
        base = n.split(":")[-1]
        if base == "_cons" or "." in base or "#" in base:
            continue
        if base not in names and s.data.has(base):
            names.append(base)
    if body:
        names = s.expand_varlist(body)
    out = s.output
    N = int(mask.sum())
    # Stata 14 layout (compat 0802): values in %9.0g, 67 dashes
    out.write(f"\n  Estimation sample {est.cmd}".ljust(44) + f"Number of obs = {N:>10,}\n", "text")
    out.write("\n  " + "-" * 67 + "\n", "text")
    out.write(f"  {'Variable':>12} |{'Mean':>13}{'Std. Dev.':>15}{'Min':>12}{'Max':>13}\n", "text")
    out.write("  " + "-" * 13 + "+" + "-" * 53 + "\n", "text")
    from ..lang.functions import _abbrev
    stats = []
    for n in names:
        x = s.data.get(n).data[mask]
        x = x[x < SYS]
        mean, sd = float(x.mean()), float(x.std(ddof=1)) if len(x) > 1 else SYS
        mn, mx = float(x.min()), float(x.max())
        stats.append([mean, sd, mn, mx])
        out.write(f"  {_abbrev(n, 12):>12} |", "text")
        out.write(f"{_g(mean):>13}{_g(sd):>14}{_g(mn):>13}{_g(mx):>13}\n", "result")
    out.write("  " + "-" * 67 + "\n", "text")
    from ..commands.matrix import Matrix
    s.r = {"stats": Matrix(np.array(stats), names, ["mean", "sd", "min", "max"])}


# ---------------------------------------------------------------------------
# regress: hettest e ovtest
# ---------------------------------------------------------------------------

def _regress_data(s: "Session", est):
    from .postest import design_all, esample
    if est.cmd != "regress":
        raise StataError(321, f"estat {''} not valid")   # VERIFY
    mask = esample(s)
    X, _ = design_all(s, est)
    keep = np.abs(np.diag(est.V)) > 0
    X = X[mask][:, keep]
    y = s.data.get(est.depvar).data[mask].astype(float)
    b = est.b[keep]
    return X, y, X @ b, mask


def _ols(X: np.ndarray, y: np.ndarray):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ beta
    return beta, e


def _hettest(s: "Session", est, rest: str) -> None:
    body, o = _opts(rest, {"rhs": 3, "iid": 3, "fstat": 3, "normal": 3, "mtest": 2})
    if est.cmd != "regress":
        raise StataError(321, "estat hettest not valid")
    X, y, yhat, mask = _regress_data(s, est)
    e = y - yhat
    N = len(y)
    if o.get("rhs"):
        Z = X[:, :-1] if est.names and est.names[-1] == "_cons" else X
        Z = Z[:, [j for j in range(Z.shape[1]) if np.ptp(Z[:, j]) > 0]]
        desc = " ".join(n for n in est.names if n != "_cons" and "o." not in n and "b." not in n)
    elif body:
        names = s.expand_varlist(body)
        Z = np.column_stack([s.data.get(n).data[mask].astype(float) for n in names])
        desc = " ".join(names)
    else:
        Z = yhat[:, None]
        desc = f"fitted values of {est.depvar}"
    Zc = np.column_stack([Z, np.ones(N)])
    df = Z.shape[1]
    if o.get("iid") or o.get("fstat"):
        u = e * e
        _, r = _ols(Zc, u)
        r2 = 1 - (r @ r) / ((u - u.mean()) @ (u - u.mean()))
        stat = N * r2
        kind = "iid"
    else:
        sig2 = (e @ e) / N
        u = e * e / sig2
        # ESS with centered regressors, without subtracting large sums (compat 0802: chi2 close to
        # zero)
        uc = u - u.mean()
        Zm = Z - Z.mean(axis=0)
        beta = np.linalg.lstsq(Zm, uc, rcond=None)[0]
        stat = float(uc @ (Zm @ beta)) / 2
        kind = "normal"
    p = float(_chi2_tail(stat, df))
    out = s.output
    out.write("\nBreusch-Pagan / Cook-Weisberg test for heteroskedasticity \n", "text")
    out.write("         Ho: Constant variance\n", "text")
    out.write(f"         Variables: {desc}\n", "text")
    out.write("\n", "text")
    out.write(f"         chi2({df})      =  {stat:>7.2f}\n", "text")
    out.write(f"         Prob > chi2  =   {p:.4f}\n", "text")
    s.r = dict(reversed([("p", p), ("df", float(df)), ("chi2", stat)]))


def _ovtest(s: "Session", est, rest: str) -> None:
    _, o = _opts(rest, {"rhs": 3})
    if est.cmd != "regress":
        raise StataError(321, "estat ovtest not valid")
    X, y, yhat, _ = _regress_data(s, est)
    N, k = X.shape
    if o.get("rhs"):
        base = [j for j in range(k) if np.ptp(X[:, j]) > 0]
        extra = np.column_stack([X[:, j] ** p for j in base for p in (2, 3, 4)])
        desc = "powers of the independent variables"   # VERIFY
    else:
        # [0, 1] scale and powers stored as float, as Stata does (compat 0802: F =
        # .2525909782551101)
        z = (yhat - yhat.min()) / (yhat.max() - yhat.min())
        f32 = lambda a: a.astype(np.float32).astype(float)   # noqa: E731
        extra = np.column_stack([f32(z ** 2), f32(z ** 3), f32(z ** 4)])
        desc = "powers of the fitted values of " + est.depvar
    _, e0 = _ols(X, y)
    Xa = np.column_stack([X, extra])
    _, e1 = _ols(Xa, y)
    q = extra.shape[1]
    rank = np.linalg.matrix_rank(Xa)
    df_r = N - rank
    F = ((e0 @ e0 - e1 @ e1) / q) / ((e1 @ e1) / df_r)
    p = float(st.f.sf(F, q, df_r))
    out = s.output
    out.write(f"\nRamsey RESET test using {desc}\n", "text")
    out.write("       Ho:  model has no omitted variables\n", "text")
    out.write(f"{'F(' + str(q) + ', ' + str(df_r) + ')':>26} = {F:>9.2f}\n", "text")
    out.write(f"{'Prob > F':>26} = {p:>11.4f}\n", "text")
    s.r = dict(reversed([("p", p), ("F", F), ("df", float(q)), ("df_r", float(df_r))]))


# ---------------------------------------------------------------------------
# logit/probit: classification e gof
# ---------------------------------------------------------------------------

_BINARY = ("logit", "logistic", "probit")


def _binary(s: "Session", est):
    from .postest import esample, xb_all
    if est.cmd not in _BINARY:
        raise StataError(321, f"estat not valid after {est.cmd}")   # VERIFY
    mask = esample(s)
    xb, _ = xb_all(s, est)
    xb = xb[mask]
    p = st.norm.cdf(xb) if est.cmd == "probit" else 1 / (1 + np.exp(-xb))
    y = s.data.get(est.depvar).data[mask] != 0
    return y, p, mask


def _model_title(est) -> str:
    return ("Probit" if est.cmd == "probit" else "Logistic") + f" model for {est.depvar}"


def _classification(s: "Session", est, rest: str) -> None:
    _, o = _opts(rest, {"cutoff": 3, "all": 3})
    cut = float(o["cutoff"]) if isinstance(o.get("cutoff"), str) else 0.5
    y, p, _ = _binary(s, est)
    pos = p >= cut
    a = int((pos & y).sum())
    b = int((pos & ~y).sum())
    c = int((~pos & y).sum())
    d = int((~pos & ~y).sum())
    n = a + b + c + d
    out = s.output
    out.write(f"\n{_model_title(est)}\n\n", "text")
    out.write("              -------- True --------\n", "text")
    out.write("Classified |         D            ~D  |      Total\n", "text")
    out.write("-----------+--------------------------+-----------\n", "text")
    for lab, u, v in (("+", a, b), ("-", c, d)):
        out.write(f"     {lab}     |", "text")
        out.write(f"{u:>10}{v:>14}  |{u + v:>11}\n", "result")
    out.write("-----------+--------------------------+-----------\n", "text")
    out.write("   Total   |", "text")
    out.write(f"{a + c:>10}{b + d:>14}  |{n:>11}\n", "result")
    cut_s = format_value(cut, "%9.0g", pad=False).strip()
    out.write(f"\nClassified + if predicted Pr(D) >= {cut_s}\n", "text")
    out.write(f"True D defined as {est.depvar} != 0\n", "text")
    out.write("-" * 50 + "\n", "text")

    def pct(num, den):
        return 100 * num / den if den else SYS

    rows1 = [("Sensitivity", "Pr( +| D)", pct(a, a + c)), ("Specificity", "Pr( -|~D)", pct(d, b + d)),
             ("Positive predictive value", "Pr( D| +)", pct(a, a + b)),
             ("Negative predictive value", "Pr(~D| -)", pct(d, c + d))]
    rows2 = [("False + rate for true ~D", "Pr( +|~D)", pct(b, b + d)),
             ("False - rate for true D", "Pr( -| D)", pct(c, a + c)),
             ("False + rate for classified +", "Pr(~D| +)", pct(b, a + b)),
             ("False - rate for classified -", "Pr( D| -)", pct(c, c + d))]
    for rows in (rows1, rows2):
        for lab, pr, v in rows:
            val = ".%" if v >= SYS else f"{v:.2f}%"
            out.write(f"{lab:<32}{pr}", "text")
            out.write(f"{val:>9}\n", "result")
        out.write("-" * 50 + "\n", "text")
    cc = pct(a + d, n)
    out.write(f"{'Correctly classified':<41}", "text")
    out.write(f"{cc:>8.2f}%\n", "result")
    out.write("-" * 50 + "\n", "text")
    # r(P_xy): Pr(x|y) with p/n = classified +/-, 1/0 = true D/~D (VERIFY the return list order)
    s.r = dict(reversed([("P_1n", pct(c, c + d)), ("P_0p", pct(b, a + b)), ("P_0n", pct(d, c + d)),
                         ("P_1p", pct(a, a + b)), ("P_n1", pct(c, a + c)), ("P_p0", pct(b, b + d)),
                         ("P_n0", pct(d, b + d)), ("P_p1", pct(a, a + c)), ("P_corr", cc)]))


def _gof(s: "Session", est, rest: str) -> None:
    _, o = _opts(rest, {"group": 1, "all": 3, "table": 1, "outsample": 3})
    y, p, mask = _binary(s, est)
    out = s.output
    N = len(y)
    out.write(f"\n{_model_title(est)}, goodness-of-fit test\n\n", "text")
    if isinstance(o.get("group"), str):
        G = int(o["group"])
        # groups by the quantiles of the predicted probabilities (like xtile)
        order = np.argsort(p, kind="stable")
        ps = p[order]
        cuts = [ps[min(N - 1, int(math.ceil(N * q / G)) - 1)] for q in range(1, G)]
        grp = np.searchsorted(np.array(cuts), p, side="left")
        chi = 0.0
        ng = 0
        for g in range(G):
            sel = grp == g
            if not sel.any():
                continue
            ng += 1
            obs1, exp1 = float(y[sel].sum()), float(p[sel].sum())
            n_g = float(sel.sum())
            if exp1 > 0 and n_g - exp1 > 0:
                chi += (obs1 - exp1) ** 2 / (exp1 * (1 - exp1 / n_g))
        df = ng - 2
        pv = float(_chi2_tail(chi, df))
        out.write("  (Table collapsed on quantiles of estimated probabilities)\n\n", "text")
        out.write(f"{'number of observations':>29} = {N:>9}\n", "text")
        out.write(f"{'number of groups':>29} = {ng:>9}\n", "text")
        out.write(f"{'Hosmer-Lemeshow chi2(' + str(df) + ')':>29} = {chi:>12.2f}\n", "text")
        out.write(f"{'Prob > chi2':>29} = {pv:>14.4f}\n", "text")
        s.r = dict(reversed([("df", float(df)), ("chi2", chi), ("m", float(ng)), ("N", float(N))]))   # VERIFY
        return
    # covariate patterns: identical rows of X
    from .postest import design_all
    X, _ = design_all(s, est)
    X = X[mask]
    _, inv = np.unique(np.round(X, 12), axis=0, return_inverse=True)
    inv = inv.ravel()
    J = int(inv.max()) + 1
    chi = 0.0
    for j in range(J):
        sel = inv == j
        m = float(sel.sum())
        pj = float(p[sel][0])
        yj = float(y[sel].sum())
        chi += (yj - m * pj) ** 2 / (m * pj * (1 - pj))
    k = int(np.count_nonzero(np.abs(np.diag(est.V)) > 0))
    df = J - k
    pv = float(_chi2_tail(chi, df))
    out.write(f"{'number of observations':>29} = {N:>9}\n", "text")
    out.write(f"{'number of covariate patterns':>29} = {J:>9}\n", "text")
    out.write(f"{'Pearson chi2(' + str(df) + ')':>29} = {chi:>12.2f}\n", "text")
    out.write(f"{'Prob > chi2':>29} = {pv:>14.4f}\n", "text")
    s.r = dict(reversed([("df", float(df)), ("chi2", chi), ("m", float(J)), ("N", float(N))]))


def _bootstrap(s: "Session", est, rest: str) -> None:
    from .bootstrap import estat_bootstrap
    estat_bootstrap(s, est, rest)


def _effects(s: "Session", est, rest: str) -> None:
    from .svy import estat_effects
    estat_effects(s, est, rest)
