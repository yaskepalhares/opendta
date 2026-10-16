"""xtgee and the population-averaged estimator of xtreg, xtlogit and xtpoisson (option pa)
([XT] xtgee).

Generalized estimating equations of Liang and Zeger (1986): the mean follows a GLM (family
and link), observations of the same panel are correlated through a working correlation matrix
(independent, exchangeable, ar #, stationary #, nonstationary #, unstructured), and the
coefficients solve sum_i D_i' V_i^-1 (y_i - mu_i) = 0 by Fisher scoring. Model-based and robust
(sandwich, clustered on the panels) variances.

Undocumented details are marked VERIFY: the divisors of the scale and correlation estimators
(N by default, N-p with nmp), the convergence measure printed in the iteration log and the
robust variance without a small-sample factor."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from . import fvars
from .glm import _FAMILY_ALIASES, _FAMILY_DEFAULT_LINK, _link
from .models2 import _Prep, _bmap, _expand, _semap
from .results import Estimates, coef_table, post, rows_from_columns

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_CORR_NAMES = {"independent": ("independent", 3), "exchangeable": ("exchangeable", 3),
               "ar": ("ar", 2), "stationary": ("stationary", 3), "nonstationary": ("nonstationary", 3),
               "unstructured": ("unstructured", 3)}
_FAMILY_TEXT = {"gaussian": "Gaussian", "binomial": "binomial", "poisson": "Poisson", "gamma": "gamma",
                "igaussian": "inverse Gaussian", "nbinomial": "negative binomial"}
_LINK_TEXT = {"Identity": "identity", "Log": "log", "Logit": "logit", "Probit": "probit",
              "Complementary log-log": "cloglog"}


def _g(x, fmt="%9.0g"):
    if x is None or x >= SYS or not np.isfinite(x):
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _parse_corr(text: str) -> tuple[str, int]:
    words = text.replace("(", " ").replace(")", " ").split()
    if not words:
        return "exchangeable", 0
    head = words[0].lower()
    m = re.fullmatch(r"([a-z]+?)(\d*)", head)
    name, num = (m.group(1), m.group(2)) if m else (head, "")
    for full, (label, mn) in _CORR_NAMES.items():
        if len(name) >= mn and full.startswith(name) or (full == "ar" and name == "ar"):
            k = int(num or (words[1] if len(words) > 1 else 1)) if full in (
                "ar", "stationary", "nonstationary") else 0
            return label, k
    raise StataError(198, f"invalid correlation structure {text}")   # VERIFY


def has_option(args: str, name: str) -> bool:
    from ..lang.syntax import find_top, parse_options
    c = find_top(args, ",")
    if c == -1:
        return False
    return any(n.lower() == name for n, _ in parse_options(args[c + 1:]))


def run_pa(s: "Session", cmd: str, args: str) -> bool:
    """xtreg, pa / xtlogit, pa / xtpoisson, pa: GEE with an exchangeable correlation by
    default. Returns False when the option is absent."""
    if not has_option(args, "pa"):
        return False
    family, link = {"xtreg": ("gaussian", "identity"), "xtlogit": ("binomial", "logit"),
                    "xtpoisson": ("poisson", "log"), "xtprobit": ("binomial", "probit")}[cmd]
    cmd_xtgee(s, pa_args(args, family, link), cmd2=cmd)
    return True


def pa_args(args: str, family: str, link: str) -> str:
    """Options of xtreg/xtlogit/xtpoisson, pa as xtgee options: drops pa (and re/fe) and adds
    family() and link()."""
    from ..lang.syntax import find_top
    c = find_top(args, ",")
    if c == -1:
        return args
    head, opts = args[:c], args[c + 1:]
    from ..lang.syntax import parse_options
    keep = []
    for name, val in parse_options(opts):
        low = name.lower()
        if low == "pa":
            continue
        keep.append(name if val is None else f"{name}({val})")
    keep += [f"family({family})", f"link({link})"]
    return head + ", " + " ".join(keep)


def cmd_xtgee(s: "Session", args: str, *, cmd2: str = "") -> None:
    from .postest import replay
    from .tsops import panel_info, ts_info
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtgee", args)
        return
    pr = _Prep(s, args, {"family": 1, "link": 1, "corr": 1, "i": 1, "t": 1, "eform": 2, "nmp": 3,
                         "scale": 3, "exposure": 3, "offset": 3, "force": 5, "or": 2, "irr": 3,
                         "rgf": 3, "nodisplay": 5},
               weights=("fweight", "iweight", "pweight"))
    o, smp = pr.o, pr.smp
    out = s.output
    ds = s.data
    pr.write_notes(out)
    pvar = str(o["i"]).strip() if isinstance(o.get("i"), str) else panel_info(ds)
    if not pvar:
        raise StataError(459, "panel variable not set; use xtset varname ...")
    pvar = s.expand_varlist(pvar)[0]
    info = ts_info(ds)
    tvar = str(o["t"]).strip() if isinstance(o.get("t"), str) else (info.tvar if info else "")
    if tvar:
        tvar = s.expand_varlist(tvar)[0]
    fam_words = str(o.get("family") or "gaussian").split()
    fam = _FAMILY_ALIASES.get(fam_words[0].lower())
    if fam is None:
        raise StataError(198, f"unknown family {fam_words[0]}")   # VERIFY
    lk = str(o.get("link") or "").split()
    lname, lpow = (lk[0], float(lk[1]) if len(lk) > 1 else None) if lk else _FAMILY_DEFAULT_LINK[fam]
    gfun, ginv0, dmu0, _ltext, link_name = _link(lname, lpow)
    corr, corr_k = _parse_corr(str(o.get("corr") or "exchangeable"))
    mask = smp.mask
    y = pr.y
    n = len(y)
    m_trials = np.ones(n)
    if fam == "binomial" and len(fam_words) > 1:
        tok = fam_words[1]
        try:
            m_trials = np.full(n, float(tok))
        except ValueError:
            m_trials = ds.get(tok).data[mask].astype(float)
    nb_k = float(fam_words[1]) if fam == "nbinomial" and len(fam_words) > 1 else 1.0

    # binomial with trials: the mean is m * p
    def ginv(e):
        return m_trials * ginv0(e) if fam == "binomial" else ginv0(e)

    def dmu(e):
        return m_trials * dmu0(e) if fam == "binomial" else dmu0(e)
    if fam == "binomial" and (np.any(y < 0) or np.any(y > m_trials)):
        raise StataError(2000, f"{smp.depname} has values outside [0, n]")   # VERIFY
    if fam == "poisson" and np.any(y < 0):
        raise StataError(459, f"{smp.depname} must be greater than or equal to zero")   # VERIFY
    off = np.zeros(n)
    if isinstance(o.get("exposure"), str):
        off = np.log(ds.get(s.expand_varlist(o["exposure"])[0]).data[mask].astype(float))
    elif isinstance(o.get("offset"), str):
        off = ds.get(s.expand_varlist(o["offset"])[0]).data[mask].astype(float)
    ids = ds.get(pvar).data[mask]
    if np.any(ids >= SYS):
        raise StataError(459, f"{pvar} missing in sample")   # VERIFY
    if corr in ("ar", "stationary", "nonstationary", "unstructured"):
        if not tvar:
            raise StataError(459, "time variable not set; use xtset idvar timevar or the t() "
                                  "option")   # VERIFY
    tt = ds.get(tvar).data[mask].astype(float) if tvar else np.arange(n, dtype=float)
    order = np.lexsort((tt, ids))
    uniq, inv = np.unique(ids, return_inverse=True)
    G = len(uniq)
    groups = [order[inv[order] == g] for g in range(G)]
    # time positions: 0, 1, ... from the first period of the sample (gaps keep their place)
    delta = info.delta if (info and tvar == info.tvar) else 1.0
    pos = np.zeros(n, dtype=np.int64)
    if tvar:
        t0 = tt.min()
        pos = np.round((tt - t0) / delta).astype(np.int64)
    else:
        for idx in groups:
            pos[idx] = np.arange(len(idx))
    if corr in ("ar", "stationary", "nonstationary", "unstructured"):
        for idx in groups:
            if len(set(pos[idx].tolist())) != len(idx):
                raise StataError(451, "repeated time values within panel")   # VERIFY
    w_obs = pr.w if smp.wtype else np.ones(n)
    wg = np.ones(G)
    if smp.wtype:
        for g, idx in enumerate(groups):
            if np.ptp(w_obs[idx]) > 0:
                raise StataError(199, f"weight must be constant within {pvar}")   # VERIFY
            wg[g] = w_obs[idx[0]]
    X = np.column_stack([pr.X, np.ones(n)]) if pr.constant else pr.X
    p = X.shape[1]
    nmp = bool(o.get("nmp"))
    scale_fixed = fam in ("binomial", "poisson", "nbinomial")
    sc_opt = str(o.get("scale") or "").strip().lower()

    def var(mu):
        if fam == "gaussian":
            return np.ones_like(mu)
        if fam == "binomial":
            return mu * (1 - mu / m_trials)
        if fam == "poisson":
            return mu
        if fam == "gamma":
            return mu ** 2
        if fam == "igaussian":
            return mu ** 3
        return mu + nb_k * mu ** 2

    T = int(pos.max()) + 1 if n else 1

    def working(r, phi):
        """Working correlation for each panel (callable on the positions) and its parameter
        matrix R (T x T)."""
        if corr == "independent":
            return np.eye(T)
        if corr == "exchangeable":
            num = den = 0.0
            for g, idx in enumerate(groups):
                ri = r[idx]
                sr = ri.sum()
                num += wg[g] * (sr * sr - (ri * ri).sum())
                den += wg[g] * len(idx) * (len(idx) - 1)
            a = num / (phi * (den - (p if nmp else 0))) if den > 0 else 0.0
            R = np.full((T, T), a)
            np.fill_diagonal(R, 1.0)
            return R
        # sums of r_s r_t over pairs of positions
        S = np.zeros((T, T))
        C = np.zeros((T, T))
        for g, idx in enumerate(groups):
            pi, ri = pos[idx], r[idx]
            S[np.ix_(pi, pi)] += wg[g] * np.outer(ri, ri)
            C[np.ix_(pi, pi)] += wg[g]
        if corr in ("unstructured", "nonstationary"):
            with np.errstate(invalid="ignore", divide="ignore"):
                R = np.where(C > 0, S / (phi * np.maximum(C - (p if nmp else 0), 1)), 0.0)
            np.fill_diagonal(R, 1.0)
            if corr == "nonstationary":
                lag = np.abs(np.subtract.outer(np.arange(T), np.arange(T)))
                R[lag > corr_k] = 0.0
            return (R + R.T) / 2
        lags = np.abs(np.subtract.outer(np.arange(T), np.arange(T)))
        rho = np.zeros(T)
        rho[0] = 1.0
        for L in range(1, min(T, max(corr_k, 1) + 1)):
            sel = lags == L
            cnt = C[sel].sum() / 2
            rho[L] = (S[sel].sum() / 2) / (phi * cnt) if cnt > 0 else 0.0
        if corr == "stationary":
            rho[corr_k + 1:] = 0.0
            return rho[lags]
        # ar(k): Yule-Walker on the first k autocorrelations, then the implied correlations
        k = corr_k
        Rk = np.array([[rho[abs(i - j)] for j in range(k)] for i in range(k)])
        try:
            phis = np.linalg.solve(Rk, rho[1:k + 1])
        except np.linalg.LinAlgError:
            phis = np.zeros(k)
        full = np.zeros(T)
        full[:k + 1] = rho[:k + 1]
        for L in range(k + 1, T):
            full[L] = float(phis @ full[L - k:L][::-1])
        return full[lags]

    def fit_step(b, R):
        eta = X @ b + off
        mu = ginv(eta)
        d = dmu(eta)
        v = var(mu)
        H = np.zeros((p, p))
        U = np.zeros(p)
        parts = []
        for g, idx in enumerate(groups):
            pi = pos[idx]
            A = np.sqrt(v[idx])
            Vi = A[:, None] * R[np.ix_(pi, pi)] * A[None, :]
            Di = X[idx] * d[idx][:, None]
            Wi = np.linalg.pinv(Vi)
            DW = Di.T @ Wi
            ei = y[idx] - mu[idx]
            H += wg[g] * DW @ Di
            U += wg[g] * DW @ ei
            parts.append(DW @ ei)
        return H, U, mu, v, parts

    # start: independence (GLM by Fisher scoring) from the mean of the outcome
    b = np.zeros(p)
    if pr.constant:
        ybar = float(np.mean(y / m_trials)) if fam == "binomial" else float(np.mean(y))
        if fam == "binomial":
            ybar = min(max(ybar, 1e-4), 1 - 1e-4)
        with np.errstate(all="ignore"):
            start = float(gfun(np.array([max(ybar, 1e-8) if fam != "gaussian" else ybar]))[0])
        b[-1] = start if np.isfinite(start) else 0.0
    R = np.eye(T)
    for _ in range(50):
        H, U, mu, v, _parts = fit_step(b, R)
        step = np.linalg.lstsq(H, U, rcond=None)[0]
        b = b + step
        if np.max(np.abs(step)) < 1e-10:
            break
    maxit = int(float(o["iterate"])) if isinstance(o.get("iterate"), str) else 100
    tol = float(o["tolerance"]) if isinstance(o.get("tolerance"), str) else 1e-6
    log = []
    converged = False
    phi = 1.0
    for it in range(1, maxit + 1):
        eta = X @ b + off
        mu = ginv(eta)
        r = (y - mu) / np.sqrt(var(mu))
        wobs = wg[inv]
        phi = float((wobs * r * r).sum()) / ((wobs.sum()) - (p if nmp else 0))
        R = working(r, phi)
        H, U, mu, v, parts = fit_step(b, R)
        step = np.linalg.lstsq(H, U, rcond=None)[0]
        bnew = b + step
        dif = float(np.max(np.abs(step) / (np.abs(b) + 1)))   # VERIFY: Stata's tolerance measure
        log.append(dif)
        b = bnew
        if dif < tol:
            converged = True
            break
    if not converged:
        out.write("convergence not achieved\n", "error")   # VERIFY
    eta = X @ b + off
    mu = ginv(eta)
    r = (y - mu) / np.sqrt(var(mu))
    wobs = wg[inv]
    pearson = float((wobs * r * r).sum())
    phi_hat = pearson / (wobs.sum() - (p if nmp else 0))
    H, U, mu, v, parts = fit_step(b, R)
    Hinv = np.linalg.pinv(H)
    if sc_opt in ("x2", "chi2"):
        phi_use = phi_hat
    elif sc_opt in ("dev",):
        phi_use = phi_hat
    elif sc_opt and sc_opt not in ("phi",):
        try:
            phi_use = float(sc_opt)
        except ValueError:
            raise StataError(198, "option scale() invalid")   # VERIFY
    else:
        phi_use = 1.0 if scale_fixed else phi_hat
    robust = smp.vce in ("robust", "cluster")
    if robust:
        B = np.zeros((p, p))
        for g, part in enumerate(parts):
            B += wg[g] * np.outer(part, part)
        V = Hinv @ B @ Hinv
    else:
        V = phi_use * Hinv
    V = (V + V.T) / 2

    if not o.get("nolog") and log:
        out.write("\n", "text")
        for k, d_ in enumerate(log, start=1):
            out.write(f"Iteration {k}: tolerance = {_g(d_)}\n", "text")
    k_slopes = pr.X.shape[1]
    bs, Vs = b[:k_slopes], V[:k_slopes, :k_slopes]
    try:
        chi2 = float(bs @ np.linalg.solve(Vs, bs)) if k_slopes else SYS
    except np.linalg.LinAlgError:
        chi2 = SYS
    pchi = chi2_tail(chi2, k_slopes) if chi2 < SYS else SYS
    kx = p
    rows = rows_from_columns(smp.cols, _bmap(pr, b, kx), _semap(pr, V, kx), ds, constant=pr.constant)
    layout = [(c.name if not c.omitted else fvars.omitted_name(c.name), c.values is not None)
              for c in smp.cols]
    if pr.constant:
        layout.append(("_cons", True))
    b_full, V_full = _expand(pr, b, V, layout)
    est = Estimates("xtgee", smp.depname, [nm for nm, _ in layout], b_full, V_full, n, "z", None)
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, pr.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    eform = ""
    if o.get("eform") or o.get("or") or o.get("irr"):
        eform = ("Odds Ratio" if link_name == "Logit" else "IRR" if link_name == "Log" else "exp(b)")
        if o.get("or"):
            eform = "Odds Ratio"
        if o.get("irr"):
            eform = "IRR"
    est.coef_title = eform or "Coef."
    est.eform_on = bool(eform)
    cnt = np.array([len(ix) for ix in groups])
    x = {"N": n, "G": G, "pvar": pvar, "tmin": int(cnt.min()), "tavg": float(cnt.mean()),
         "tmax": int(cnt.max()), "chi2": chi2, "p": pchi, "df": k_slopes, "phi": phi_use,
         "family": _FAMILY_TEXT[fam], "link": _LINK_TEXT.get(link_name, link_name.lower()),
         "corr": corr if corr in ("independent", "exchangeable", "unstructured")
         else f"{corr}({corr_k})" if corr != "ar" else f"AR({corr_k})", "robust": robust}
    est.extra = x
    est.glm_ginv = ginv
    est.predict = _predict
    est.display = display_xtgee

    title = "GEE population-averaged model"
    scalars = [("N", float(n)), ("N_g", float(G)), ("df_m", float(k_slopes)), ("chi2", chi2),
               ("p", pchi), ("df_pear", float(n - p)), ("deviance", SYS), ("dispers", SYS),
               ("phi", phi_use), ("tol", tol), ("dif", log[-1] if log else 0.0),
               ("g_min", float(cnt.min())), ("g_avg", float(cnt.mean())), ("g_max", float(cnt.max())),
               ("rank", float(np.linalg.matrix_rank(V))), ("rc", 0.0 if converged else 430.0)]
    macros = [("cmdline", f"{cmd2 or 'xtgee'} {args.strip()}"), ("cmd", "xtgee"), ("cmd2", cmd2),
              ("depvar", smp.depname), ("family", _FAMILY_TEXT[fam]),
              ("link", x["link"]), ("corr", x["corr"]), ("scale", sc_opt or ("1" if scale_fixed else "phi")),
              ("ivar", pvar), ("tvar", tvar), ("vce", "robust" if robust else "conventional"),
              ("vcetype", "Robust" if robust else ""), ("chi2type", "Wald"), ("model", "pa"),
              ("title", title), ("properties", "b V"), ("predict", "xtgee_p")]
    from ..commands.matrix import Matrix
    Rshow = R if corr != "exchangeable" else np.array([[1.0, R[0, 1] if T > 1 else 0.0],
                                                       [R[0, 1] if T > 1 else 0.0, 1.0]])
    Rmat = Matrix(Rshow, [f"r{i + 1}" for i in range(Rshow.shape[0])],
                  [f"c{i + 1}" for i in range(Rshow.shape[1])])
    post(s, est, list(reversed(scalars)), list(reversed(macros)), mask, extra_mats=[("R", Rmat)])
    if not o.get("nodisplay"):
        display_xtgee(s, est)


def display_xtgee(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                  level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        left = [("GEE population-averaged model", ""), ("Group variable:", x["pvar"]),
                ("Link:", x["link"]), ("Family:", x["family"]), ("Correlation:", x["corr"]),
                ("", ""), ("", ""), ("Scale parameter:", _g(x["phi"]))]
        right = [("Number of obs", f"{x['N']:,}"), ("Number of groups", f"{x['G']:,}"),
                 ("Obs per group:", None), ("min", f"{x['tmin']:,}"), ("avg", f"{x['tavg']:.1f}"),
                 ("max", f"{x['tmax']:,}"), (f"Wald chi2({x['df']})", _g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", _g(x["p"], "%10.4f"))]
        for (ll, lv), (rl, rv) in zip(left, right):
            if lv:
                lt = f"{ll:<16}{lv:>27}"
            else:
                lt = ll
            out.write(f"{lt:<48}", "text")
            if rv is None:
                out.write(f"{rl}\n", "text")
            elif rl in ("min", "avg", "max"):
                out.write(f"{'':<14}{rl} = ", "text")
                out.write(f"{rv:>10}\n", "result")
            else:
                out.write(f"{rl:<18}= ", "text")
                out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if table:
        if x.get("robust"):
            out.write(f"{'(Std. Err. adjusted for clustering on ' + x['pvar'] + ')':>78}\n", "text")
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level,
                   vcetype=est.vcetype, coef_title=getattr(est, "coef_title", "Coef."),
                   eform=getattr(est, "eform_on", False))


def _predict(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import design_all, xb_all
    xb, ok = xb_all(s, est)
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat == "stdp":
        X, _ = design_all(s, est)
        sp = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", X, est.V, X), 0))
        return np.where(ok, sp, SYS), ""
    if stat in ("", "mu", "pr", "n"):
        note = "(option mu assumed; predicted mean " + est.depvar + ")" if not stat else ""
        return np.where(ok, est.glm_ginv(xb), SYS), note
    raise StataError(198, f"option {stat} not allowed")

