"""ivregress ([R] ivregress): two-stage least squares, LIML and GMM with the weighting matrices
and variance estimators of the manual."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail as _chi2_tail
from ..lazy import lazy
from . import fvars
from .panel import _predict_xt, g
from .regress import Sample, cmdline, comma, sig5
from .results import Estimates, coef_table, post, rows_from_columns

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

SYS = M.SYSMISS


def cmd_ivregress(s: "Session", args: str) -> None:
    from .postest import replay
    t = args.strip()
    if not t or t.startswith(","):
        replay(s, "ivregress", args)
        return
    m = re.match(r"^(\w+)\s+(.*)$", t, re.S)
    if not m or m.group(1) not in ("2sls", "liml", "gmm"):
        raise StataError(198, "estimator must be 2sls, liml, or gmm")   # VERIFY
    estimator = m.group(1)
    rest = m.group(2)
    mp = re.search(r"\(([^=()]+)=([^()]*)\)", rest)
    if not mp:
        raise StataError(198, "equation not specified")   # VERIFY
    endog_txt, inst_txt = mp.group(1), mp.group(2)
    exog_txt = rest[:mp.start()] + " " + rest[mp.end():]
    # sample: dependent + exogenous + endogenous + instruments without missing
    smp = Sample(s, exog_txt.strip(), opts_spec={"small": 2, "first": 3, "wmatrix": 2, "center": 6, "igmm": 4},
                 weights=("aweight", "fweight", "pweight"), extra_vars=endog_txt + " " + inst_txt)
    if estimator != "gmm" and smp.o.get("wmatrix"):
        raise StataError(198, "option wmatrix() is allowed only with gmm")   # VERIFY
    for opt_ in ("center", "igmm"):
        if estimator != "gmm" and smp.o.get(opt_):
            raise StataError(198, f"option {opt_} is allowed only with gmm")   # VERIFY
    ds = s.data
    endog_terms = fvars.expand_fv(ds, endog_txt)
    inst_terms = fvars.expand_fv(ds, inst_txt)
    mask = smp.mask
    y = smp.y.astype(float)
    N = len(y)
    endog_cols = fvars.build_design(ds, endog_terms, mask)
    inst_cols = fvars.build_design(ds, inst_terms, mask)
    ex_cols = [c for c in smp.cols if c.values is not None and not c.base]
    Xn = [c for c in endog_cols if c.values is not None and not c.base]
    Zc = [c for c in inst_cols if c.values is not None and not c.base]
    X = np.column_stack([c.values for c in Xn] + [c.values for c in ex_cols] + [np.ones(N)])
    Z = np.column_stack([c.values for c in Zc] + [c.values for c in ex_cols] + [np.ones(N)])
    if Z.shape[1] < X.shape[1]:
        raise StataError(481, "equation not identified; must have at least as many instruments not in\n"
                         "the regression as there are instrumented variables")
    k = X.shape[1]
    wfit = None
    if smp.wtype:
        # weighted estimation: every row scaled by the square root of its (normalized) weight,
        # which gives weighted 2SLS, LIML and GMM (case 1909: svy: ivregress)
        wfit, N = smp.weights_for_fit()
        sq = np.sqrt(wfit)
        y, X, Z = y * sq, X * sq[:, None], Z * sq[:, None]
    PzX = Z @ np.linalg.lstsq(Z, X, rcond=None)[0]
    kappa = 1.0
    QX = PzX
    if estimator == "liml":
        kappa = _liml_kappa(y, Xn, ex_cols, X, Z)
        QX = X - kappa * (X - PzX)         # (I - kappa M_Z) X
    A = np.linalg.inv(QX.T @ X)
    b = A @ (QX.T @ y)
    e = y - X @ b
    rss = float(e @ e)
    small = bool(smp.o.get("small"))
    s2 = rss / (N - k) if small else rss / N
    V = s2 * A
    vcetype = ""
    n_clust = None
    wm_label = ""
    if estimator == "gmm":
        b, e, V, wm_label, n_clust_w, ginfo = _gmm_two_step(smp, y, X, Z, A @ (QX.T @ y), small)
        if smp.o.get("igmm") and not smp.o.get("nolog"):
            for i_, (db, dw) in enumerate(ginfo["trace"], 1):
                s.output.write(f"Iteration {i_}: change in beta = {db:10.3e}   change in W = {dw:10.3e}\n", "text")
        rss = float(e @ e)
        s2 = rss / (N - k) if small else rss / N
        if smp.vce in ("robust", "cluster") or wm_label != "Unadjusted":
            vcetype = "HAC" if wm_label.startswith("HAC") else "Robust"   # VERIFY the HAC label
            n_clust = n_clust_w if (smp.vce == "cluster" or wm_label.startswith("Cluster")) else None
    elif smp.vce in ("robust", "cluster") and getattr(s, "_svy_ctx", None) is not None:
        # svy: the linearized variance from the scores (QX e, already weighted) carried to the
        # design observations, as the other estimators under svy_generic (case 1909)
        ctx = s._svy_ctx
        dd = ctx["dd"]
        Zl = (QX * e[:, None]) @ A.T
        pos = ctx["pos"][np.flatnonzero(mask)]
        Zf = np.zeros((dd.n, Zl.shape[1]))
        Zf[pos] = Zl
        ctx["Z"] = Zf
        V = dd.variance(Zf)
        vcetype = "Robust"
    elif smp.vce in ("robust", "cluster"):
        sc = QX * e[:, None]
        if smp.vce == "cluster":
            uniq = np.unique(smp.cluster)
            n_clust = len(uniq)
            meat = sum(np.outer(sc[smp.cluster == u].sum(0), sc[smp.cluster == u].sum(0)) for u in uniq)
            q = n_clust / (n_clust - 1) * ((N - 1) / (N - k) if small else 1)
        else:
            meat = sc.T @ sc
            q = N / (N - k) if small else 1.0
        V = q * A @ meat @ A
        vcetype = "Robust"
    V = (V + V.T) / 2
    slopes = list(range(k - 1))
    W = float(b[slopes] @ np.linalg.solve(V[np.ix_(slopes, slopes)], b[slopes]))
    if wfit is None:
        ybar = y.mean()
        tss = float(np.sum((y - ybar) ** 2))
    else:
        yo = y / np.sqrt(wfit)
        ybar = float(np.sum(wfit * yo) / np.sum(wfit))
        tss = float(np.sum(wfit * (yo - ybar) ** 2))
    r2 = 1 - rss / tss
    rmse = np.sqrt(rss / (N - k) if small else rss / N)
    names = [c.name for c in Xn] + [c.name for c in ex_cols] + ["_cons"]
    bmap = dict(zip(names, b))
    semap = {nm: float(np.sqrt(max(V[i, i], 0))) for i, nm in enumerate(names)}
    rows = rows_from_columns(Xn + ex_cols, bmap, semap, ds, constant=True)
    stat = "t" if small else "z"
    df_r = N - k if small else None
    est = Estimates("ivregress", smp.depname, names, b, V, N, stat, df_r)
    est.rows = rows
    est.cols = Xn + ex_cols
    est.terms = endog_terms + smp.terms
    _retarget(est, Xn, ex_cols, endog_terms, smp.terms)
    est.constant, est.depcomp = True, smp.depcomp
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust = vcetype, smp.clustvar, n_clust
    est.extra = dict(N=N, W=W, df_m=len(slopes), r2=r2, rmse=rmse, small=small, estimator=estimator,
                     wm_label=wm_label,
                     F=W / len(slopes), instd=" ".join(c.name for c in Xn),
                     insts=" ".join([c.name for c in ex_cols] + [c.name for c in Zc]))
    est.display = display_iv
    est.predict = _predict_xt
    est.predict_opts = {"u": 1, "e": 1, "ue": 2, "xbu": 3}
    r2_a = 1 - (1 - r2) * (N - 1) / (N - k)
    if tss - rss < 0:
        # a negative model sum of squares: Stata reports R-squared as missing (compat 1920)
        r2 = r2_a = SYS
        est.extra["r2"] = SYS
    est.extra.update(mss=tss - rss, rss=rss, tss=tss, r2_a=r2_a, k=k)
    # e(): order observed in Stata 14 with small (compat 0510); without small, VERIFY
    gmm = estimator == "gmm"
    # gmm: J in place of kappa, the iterations, the HAC lags and S and W (compat 1918; the other
    # weight matrices follow the HAC case, VERIFY)
    scalars = [("J", ginfo["J"]) if gmm else ("kappa", kappa), ("rss", rss), ("N", float(N)),
               ("df_m", float(len(slopes)))]
    scalars += [("df_r", float(N - k))] if small else []
    scalars += [("rmse", rmse), ("mss", tss - rss), ("r2", r2), ("r2_a", r2_a)]
    scalars += [("F", W / len(slopes))] if small else [("chi2", W)]
    if n_clust is not None:
        scalars.append(("N_clust", float(n_clust)))
    scalars += [("iterations", float(ginfo["iterations"]) if gmm else 0.0)]
    if gmm and ginfo.get("opt") is not None:
        scalars += [("wlagopt", float(ginfo["opt"][0])), ("vcelagopt", float(ginfo["opt"][1]))]
    if gmm and ginfo["hac"] is not None:
        scalars.append(("hac_lag", float(ginfo["opt"][1] if ginfo.get("opt") else ginfo["hac"][1])))
    scalars.append(("rank", float(k)))
    exogr = " ".join(c.name for c in ex_cols)
    macros = [("cmdline", cmdline("ivregress", args)), ("cmd", "ivregress"), ("estat_cmd", "ivregress_estat")]
    if gmm:
        wtext = " ".join(str(smp.o.get("wmatrix") or "robust").split())
        if ginfo["hac"] is not None:
            macros.append(("hac_kernel", {"quadraticspectral": "Quadratic spectral"}.get(ginfo["hac"][0], ginfo["hac"][0].capitalize())))
            if ginfo.get("opt") is None:
                wtext = f"hac {ginfo['hac'][0]} {ginfo['hac'][1]}"
            # with opt Stata keeps the option as typed: "hac bartlett opt" (compat 1920)
        macros += [("wmatrix", wtext), ("vce", wtext if ginfo["hac"] is not None else (smp.vce if vcetype else "unadjusted"))]
    else:
        macros.append(("vce", smp.vce if vcetype else "unadjusted"))
    if vcetype:
        macros.append(("vcetype", vcetype))
    if smp.clustvar:
        macros.append(("clustvar", smp.clustvar))
    macros += [("estimator", estimator), ("footnote", "ivreg_footnote"), ("marginsnotok", "Residuals SCores"),
               ("marginsok", "XB default"), ("predict", "ivregress_p"), ("depvar", smp.depname)]
    if exogr:
        macros.append(("exogr", exogr))
    macros += [("insts", est.extra["insts"]), ("instd", est.extra["instd"])]
    if small:
        macros.append(("small", "small"))
    title = f"Instrumental variables ({estimator.upper()}) regression"
    est.extra["title"] = title
    macros += [("title", title), ("properties", "b V")]
    mats = []
    if gmm:
        from ..commands.matrix import Matrix
        zn = [c.name for c in Zc] + [c.name for c in ex_cols] + ["_cons"]
        mats = [("S", Matrix(ginfo["S"], zn, zn)), ("W", Matrix(ginfo["W"], zn, zn))]
        est.extra["hac"] = ginfo["hac"]
        est.extra["opt"] = ginfo.get("opt")
    post(s, est, scalars, macros, mask, extra_mats=mats)
    display_iv(s, est)


def _retarget(est, Xn, ex_cols, endog_terms, ex_terms):
    """The endogenous and exogenous columns use term indices of the joined list."""
    off = len(endog_terms)
    for c in ex_cols:
        c.term = c.term + off if not getattr(c, "_retargeted", False) else c.term
        c._retargeted = True


def display_iv(s, est, *, header=True, table=True, level=None, **_kw):
    out = s.output
    x = est.extra
    if header and x["small"] and not est.vcetype and not x["wm_label"]:
        # with small, ANOVA table in the old regress layout (compat 0510)
        out.write(f"\n\n{x['title']}\n\n", "text")
        df_m, df_r, N = int(x["df_m"]), int(est.df_r), x["N"]
        right = [("Number of obs", comma(N)), (f"F({df_m:>3},{df_r:>6})", g(x["F"], "%10.2f")),
                 ("Prob > F", g(st.f.sf(x["F"], df_m, df_r), "%10.4f")), ("R-squared", g(x["r2"], "%10.4f")),
                 ("Adj R-squared", g(x["r2_a"], "%10.4f")), ("Root MSE", g(x["rmse"], "%10.4f"))]

        def row(nm, ss, df, ms):
            return f"{nm:>12} | {g(ss, '%11.0g'):>11}{int(df):>6}{g(ms, '%11.0g'):>12}"
        left = [f"{'Source':>12} |       SS       df       MS", "-" * 13 + "+" + "-" * 30,
                row("Model", x["mss"], df_m, x["mss"] / df_m), row("Residual", x["rss"], df_r, x["rss"] / df_r),
                "-" * 13 + "+" + "-" * 30, row("Total", x["tss"], N - 1, x["tss"] / (N - 1))]
        for lft, (rl, rv) in zip(left, right):
            out.write(f"{lft:<50}{rl:<16}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    elif header:
        out.write("\n", "text")
        if x["small"]:
            right = [("Number of obs", comma(x["N"])), (f"F({int(x['df_m'])}, {int(est.df_r)})", g(x["F"], "%10.2f")),
                     ("Prob > F", g(st.f.sf(x["F"], x["df_m"], est.df_r), "%10.4f"))]
        else:
            right = [("Number of obs", comma(x["N"])), (f"Wald chi2({int(x['df_m'])})", g(x["W"], "%10.2f")),
                     ("Prob > chi2", g(_chi2_tail(x["W"], x["df_m"]), "%10.4f"))]
        right += [("R-squared", g(x["r2"], "%10.4f")), ("Root MSE", sig5(x["rmse"]))]
        for i, (rl, rv) in enumerate(right):
            lft = x["title"] if i == 0 else ""
            if x["wm_label"] and i == len(right) - 1:
                lft = f"GMM weight matrix: {x['wm_label']}"
            out.write(f"{lft:<50}{rl:<16}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level,
               vcetype=est.vcetype)
    out.write(f"Instrumented:  {x['instd']}\n", "text")
    out.write(f"Instruments:   {x['insts']}\n", "text")
    if x.get("hac") is not None:
        kern, lags = x["hac"][0], x["hac"][1]
        kname = {"quadraticspectral": "Quadratic spectral"}.get(kern, kern.capitalize())
        opt = x.get("opt")
        if opt is not None:
            lags = opt[1]
        lags = int(lags) if float(lags) == int(lags) else g(float(lags), "%9.0g").strip()   # VERIFY a fractional lag
        out.write(f"HAC VCE:       {kname} kernel with {lags} lags\n", "text")
        if opt is not None:
            out.write("HAC Lags:      Chosen by Newey-West method\n", "text")
            out.write("WMat Lags:     Chosen by Newey-West method\n", "text")


def _liml_kappa(y, Xn, ex_cols, X, Z) -> float:
    """Smallest eigenvalue of (Y' M_Z Y)^-1 (Y' M_X1 Y), Y = [y, endogenous regressors], X1 =
    included exogenous variables with the constant (kappa >= 1)."""
    from scipy.linalg import eigh
    g = len(Xn)
    Y = np.column_stack([y, X[:, :g]])
    X1 = X[:, g:]

    def resid(M, B):
        return B - M @ np.linalg.lstsq(M, B, rcond=None)[0]
    A_ = resid(X1, Y)
    B_ = resid(Z, Y)
    return float(eigh(A_.T @ A_, B_.T @ B_, eigvals_only=True)[0])


def _hac_weights(kernel: str, lags: int, N: int) -> list[float]:
    """Weights of the autocovariances 1..N-1 for wmatrix(hac kernel #) ([R] ivregress)."""
    out = []
    for j in range(1, N):
        if kernel == "quadraticspectral":
            x = j / (lags + 1)                                # bandwidth # + 1, as the other kernels (compat 1918)
            if not np.isfinite(x):
                out.append(0.0)
                continue
            u = 6 * np.pi * x / 5
            out.append(float(25 / (12 * np.pi ** 2 * x ** 2) * (np.sin(u) / u - np.cos(u))))
            continue
        x = j / (lags + 1)
        if x > 1:
            out.append(0.0)
        elif kernel == "bartlett":
            out.append(1 - x)
        else:                                             # parzen
            out.append(1 - 6 * x * x + 6 * x ** 3 if x <= .5 else 2 * (1 - x) ** 3)
    return out


def _nw94_lags(G: np.ndarray, kernel: str, tune: float = 20.0) -> float:
    """Newey and West's (1994) lag selection as [R] ivregress (Methods and formulas) gives it:
    f_i = (u_i z_i) h with h = 1 except a 0 for the constant (the last column of G), sigma_j the
    autocovariances of f up to m* = int(tune (T/100)^a), s(0) and s(q), gamma = c (s(q)/s(0))^(2/(2q+1))
    and m = gamma T^(1/(2q+1)). The lag is min(int(m), m*) for Bartlett and Parzen and min(m, m*)
    for the quadratic spectral (compat 1920: 27, 24 and 22 with T = 400, all at m*)."""
    n = G.shape[0]
    f = G[:, :-1].sum(axis=1) if G.shape[1] > 1 else G[:, 0]
    q, c, a = {"bartlett": (1, 1.1447, 2 / 9), "parzen": (2, 2.6614, 4 / 25),
               "quadraticspectral": (2, 1.3221, 2 / 25)}[kernel]
    mstar = int(tune * (n / 100) ** a)
    mstar = max(min(mstar, n - 1), 1)
    sig = np.array([f[j:] @ f[:n - j] / n for j in range(mstar + 1)])
    jj = np.arange(1, mstar + 1)
    s0 = sig[0] + 2 * np.sum(sig[1:])
    sq = 2 * np.sum(jj ** q * sig[1:])
    if s0 == 0:
        return 0
    gamma = c * ((sq / s0) ** 2) ** (1 / (2 * q + 1))
    m = gamma * n ** (1 / (2 * q + 1))
    if kernel == "quadraticspectral":
        return float(min(m, mstar))
    return int(min(int(m), mstar))


def _gmm_two_step(smp, y, X, Z, b1, small):
    """Two-step (or, with igmm, iterated) efficient GMM: the weight matrix comes from the 2SLS
    residuals (robust, cluster, HAC or unadjusted; center centers the moments); returns b,
    residuals, V (unadjusted: the efficient variance; otherwise the sandwich with the same kind of
    meat), the weight matrix label and the number of clusters."""
    N = len(y)
    o = smp.o
    wm = str(o.get("wmatrix") or ("cluster " + smp.clustvar if smp.vce == "cluster" else "robust")).split()
    kind = wm[0].lower()
    center = bool(o.get("center"))
    e1 = y - X @ b1
    n_cl = None
    hac = None
    if kind == "hac":
        names = {"bartlett": "bartlett", "nwest": "bartlett", "parzen": "parzen", "gallant": "parzen",
                 "quadraticspectral": "quadraticspectral", "andrews": "quadraticspectral"}
        if len(wm) < 2 or wm[1].lower() not in names:
            raise StataError(198, "wmatrix(hac) requires a kernel")   # VERIFY
        kern = names[wm[1].lower()]
        opt = len(wm) > 2 and wm[2].lower() == "opt"
        if opt:
            tune = float(wm[3]) if len(wm) > 3 else 20.0
            Gc = Z * e1[:, None]
            if center:
                Gc = Gc - Gc.mean(axis=0)
            lags = _nw94_lags(Gc, kern, tune)
        else:
            lags = int(float(wm[2])) if len(wm) > 2 else N - 2          # VERIFY the default
        hac = (kern, lags, _hac_weights(kern, lags, N))

    def meat(e, kind_, weights=None):
        nonlocal n_cl
        G = Z * e[:, None]
        if center and not kind_.startswith("una"):
            G = G - G.mean(axis=0)
        if kind_ == "hac":
            S = G.T @ G
            for j, wj in enumerate(hac[2] if weights is None else weights, start=1):
                if wj != 0:
                    Gj = G[j:].T @ G[:-j]
                    S += wj * (Gj + Gj.T)
            return S
        if kind_.startswith("rob"):
            return G.T @ G
        if kind_.startswith("una"):
            return float(e @ e / N) * (Z.T @ Z)
        if kind_.startswith("clu"):
            from ..core.varlist import resolve_name
            cv = smp.cluster if smp.vce == "cluster" and len(wm) < 2 else \
                smp.s.data.get(resolve_name(smp.s.data, wm[1])).data[smp.mask]
            u = np.unique(cv)
            n_cl = len(u)
            T_ = np.array([G[cv == v].sum(axis=0) for v in u])
            if center:
                T_ = T_ - T_.mean(axis=0)
            return T_.T @ T_
        raise StataError(198, f"wmatrix({kind_}) is not available yet in OpenDTA")
    ZX = Z.T @ X
    W = np.linalg.inv(meat(e1, kind))
    b = None
    trace = []

    def mreldif(a_, b_):
        return float(np.max(np.abs(a_ - b_) / (np.abs(b_) + 1)))
    for _it in range(16000 if o.get("igmm") else 1):          # VERIFY igmm's limits
        Ainv = np.linalg.inv(ZX.T @ W @ ZX)
        b_new = Ainv @ (ZX.T @ W @ (Z.T @ y))
        if b is not None:
            W_new = W
            trace.append((mreldif(b_new, b), mreldif(W_new * N, W_old * N)))
        b_prev, b = b, b_new
        if not o.get("igmm"):
            break
        if trace and trace[-1][0] < 1e-6 and trace[-1][1] < 1e-6:      # VERIFY Stata's rule
            break
        W_old = W
        W = np.linalg.inv(meat(y - X @ b, kind))
    e = y - X @ b
    V = Ainv
    if smp.vce in ("robust", "cluster") or kind[:3] in ("rob", "clu", "hac"):
        # a robust or cluster weight matrix gives the sandwich variance (compat 1907)
        k2 = "hac" if kind == "hac" else ("cluster" if (smp.vce == "cluster" or kind.startswith("clu")) else "robust")
        vce_w = None
        if hac is not None and opt:
            # the VCE lag is chosen again with the GMM residuals ([R] ivregress: S uses b_GMM)
            Gv = Z * e[:, None]
            if center:
                Gv = Gv - Gv.mean(axis=0)
            vcelag = _nw94_lags(Gv, hac[0], tune)
            vce_w = _hac_weights(hac[0], vcelag, N)
        S2 = meat(e, k2, vce_w)
        V = Ainv @ (ZX.T @ W @ S2 @ W @ ZX) @ Ainv
    lagtxt = (lambda v: f"{v:.3f}" if isinstance(v, float) else str(v))   # quadratic spectral opt: 22.000 (compat 1920)
    if hac is not None:
        label = f"HAC {hac[0].capitalize() if hac[0] != 'quadraticspectral' else 'quadratic spectral'} {lagtxt(hac[1])}"   # compat 1918
    else:
        label = {"rob": "Robust", "una": "Unadjusted", "clu": f"Cluster ({wm[1] if len(wm) > 1 else smp.clustvar})"}[kind[:3]]
    Ze = Z.T @ e
    J = float(Ze @ W @ Ze)              # with the weight matrix of the estimates (compat 1918)
    info = {"J": J, "S": meat(e, kind) / N, "W": W * N, "trace": trace, "iterations": max(len(trace), 1),
            "hac": hac, "opt": None}
    if hac is not None and opt:
        info["opt"] = (hac[1], vcelag if vce_w is not None else hac[1])
    return b, e, (V + V.T) / 2, label, n_cl, info
