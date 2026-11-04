"""mi estimate ([MI] mi estimate): fits the model on each imputation and combines the results by
Rubin's rules."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import match_options
from ..lazy import lazy
from .mi import _with_data, dataset_m, nimp, require_mi
from .results import g9

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

SYS = M.SYSMISS


def mi_estimate(s: "Session", opts: str, cmd: str) -> None:
    from .results import Estimates, post
    require_mi(s)
    ds = s.data
    o = match_options(opts, {"level": 1, "dots": 3, "noisily": 1, "nosmall": 5, "vartable": 4,
                             "imputations": 3, "cmdok": 5, "post": 4, "mcerror": 3, "dftable": 3}) \
        if opts.strip() else {}
    M_ = nimp(ds)
    if M_ < 2:
        raise StataError(2000, "insufficient imputations")   # VERIFY
    ms = list(range(1, M_ + 1))
    if isinstance(o.get("imputations"), str):
        from ..lang.words import parse_numlist
        ms = [int(v) for v in parse_numlist(o["imputations"])]
    out = s.output
    bs, Vs, Ns, dfr = [], [], [], []
    names = eqs = None
    title, cmdname, vce = "", "", ""
    for m in ms:
        sub = dataset_m(ds, m)
        with _with_data(s, sub):
            if not o.get("noisily"):
                out.quiet_depth += 1
            try:
                rc = s.run(cmd, echo=False)
            finally:
                if not o.get("noisily"):
                    out.quiet_depth -= 1
            if rc:
                raise StataError(rc, "")
            e = s.e
            if "b" not in e:
                raise StataError(301, "estimation command required")   # VERIFY
            b = e["b"].data.ravel().copy()
            V = e["V"].data.copy()
            if names is None:
                names = list(e["b"].colnames)
                eqs = list(getattr(e["b"], "coleq", []) or ["_"] * len(names))
                title = str(e.get("title", "") or "")
                cmdname = str(e.get("cmd", "") or "")
                vce = str(e.get("vce", "") or "")
            bs.append(b)
            Vs.append(V)
            Ns.append(float(e.get("N", SYS)))
            dfr.append(float(e["df_r"]) if isinstance(e.get("df_r"), float) else None)
    Mn = len(ms)
    Bm = np.array(bs)
    Q = Bm.mean(axis=0)
    W = np.mean(Vs, axis=0)
    dev = Bm - Q
    B = dev.T @ dev / (Mn - 1)
    T = W + (1 + 1 / Mn) * B
    k = Q.size
    wdiag = np.diag(W)
    bdiag = np.diag(B)
    tdiag = np.diag(T)
    keep = wdiag > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(keep, (1 + 1 / Mn) * bdiag / wdiag, 0.0)
        gam = np.where(tdiag > 0, (1 + 1 / Mn) * bdiag / tdiag, 0.0)
        df_old = np.where(gam > 0, (Mn - 1) / gam ** 2, np.inf)
    small = dfr[0] is not None and not o.get("nosmall")
    if small:
        nu_c = dfr[0]
        df_obs = (nu_c + 1) / (nu_c + 3) * nu_c * (1 - gam)
        df = 1 / (1 / df_old + 1 / df_obs)
    else:
        df = df_old
    fmi = np.where(keep, (r + 2 / (df + 3)) / (r + 1), 0.0)
    # model test of the slopes (equal FMI)
    sl = [i for i in range(k) if names[i] != "_cons" and keep[i] and not names[i].startswith(("o.", "b."))
          and "b." not in names[i] and eqs[i] in (eqs[0], "_")]
    if cmdname in ("mlogit",):
        sl = [i for i in range(k) if names[i] != "_cons" and keep[i] and "b." not in names[i]]
    kk = len(sl)
    F = df1 = df2 = pF = SYS
    if kk:
        Ws, Bs = W[np.ix_(sl, sl)], B[np.ix_(sl, sl)]
        rbar = (1 + 1 / Mn) * np.trace(Bs @ np.linalg.inv(Ws)) / kk
        Tt = (1 + rbar) * Ws
        F = float(Q[sl] @ np.linalg.solve(Tt, Q[sl])) / kk
        v1 = kk * (Mn - 1)
        if v1 > 4:
            df2 = 4 + (v1 - 4) * (1 + (1 - 2 / v1) / rbar) ** 2
        else:
            df2 = v1 * (1 + 1 / kk) * (1 + 1 / rbar) ** 2 / 2
        if small:
            # VERIFY: Stata's small-sample df for the model test (Reiter 2007)
            df2 = min(df2, float(dfr[0]))
        df1 = kk
        pF = float(st.f.sf(F, df1, df2))
    N = float(np.mean(Ns))
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    est = Estimates(cmdname, names and "", names, Q, T, int(N), "t", None)
    est.eqnames = eqs
    est.level = level
    est.vcetype = ""
    est.mi = {"df": df, "r": r, "fmi": fmi, "M": Mn, "small": small, "nu_c": dfr[0], "F": F, "df1": df1,
              "df2": df2, "pF": pF, "N": N, "title": title, "vce": vce, "W": W, "B": B}
    est.display = display_mi
    # the depvar comes from the command line (first word after the command name)
    est.depvar = cmd.split()[1] if len(cmd.split()) > 1 else ""
    from ..commands.matrix import Matrix
    scalars = [("N", N), ("M_mi", float(Mn)), ("df_m_mi", float(df1) if df1 < SYS else SYS),
               ("df_r_mi", float(df2) if df2 < SYS else SYS), ("F_mi", F), ("p_mi", pF),
               ("rvi_avg_mi", float(np.mean(r[keep])) if keep.any() else SYS),
               ("fmi_max_mi", float(np.max(fmi[keep])) if keep.any() else SYS)]
    if small:
        scalars.append(("df_c_mi", float(dfr[0])))
    macros = [("mi", "mi"), ("cmd_mi", cmdname), ("cmd", "mi estimate"), ("title_mi", "Multiple-imputation estimates"),
              ("title", title), ("dfadjust_mi", "Small sample" if small else "Large sample"),
              ("modeltest_mi", "Equal FMI"), ("depvar", est.depvar), ("properties", "b V")]
    mats = [("W_mi", Matrix(W, list(names), list(names))), ("B_mi", Matrix(B, list(names), list(names))),
            ("df_mi", Matrix(df.reshape(1, -1), ["df"], list(names))),
            ("fmi_mi", Matrix(fmi.reshape(1, -1), ["fmi"], list(names)))]
    post(s, est, scalars, macros, np.zeros(ds.nobs, dtype=bool), extra_mats=mats)
    display_mi(s, est)


def display_mi(s: "Session", est, *, header: bool = True, table: bool = True, level=None, **_kw) -> None:
    out = s.output
    X = est.mi
    lev = level or est.level
    if header:
        def right(left, lab, val):
            out.write(f"{left:<48}{lab:<18}= ", "text")
            out.write(f"{val:>10}\n", "result")
        keep = np.diag(X["W"]) > 0
        out.write("\n", "text")
        right("Multiple-imputation estimates", "Imputations", f"{X['M']:,}")
        right(X["title"], "Number of obs", f"{int(round(X['N'])):,}")
        right("", "Average RVI", f"{np.mean(X['r'][keep]):.4f}")
        right("", "Largest FMI", f"{np.max(X['fmi'][keep]):.4f}")
        if X["small"]:
            right("", "Complete DF", f"{int(X['nu_c']):,}")
        dfs = X["df"][keep]

        def dfv(v):
            return f"{v:.2f}" if np.isfinite(v) else "."
        left = "DF adjustment:" + f"{'Small sample' if X['small'] else 'Large sample':>15}"
        right(left, "DF:     min", dfv(np.min(dfs)))
        right("", "        avg", dfv(np.mean(dfs)))
        right("", "        max", dfv(np.max(dfs)))
        vce = {"ols": "OLS", "oim": "OIM", "robust": "Robust", "cluster": "Robust", "opg": "OPG"}.get(
            X["vce"].lower(), X["vce"].upper() if X["vce"] else "OLS")
        if X["F"] < SYS:
            right("Model F test:" + f"{'Equal FMI':>16}", f"F({int(X['df1']):>4},{X['df2']:>7.1f})",
                  f"{X['F']:.2f}")
            right("Within VCE type:" + f"{vce:>13}", "Prob > F", f"{X['pF']:.4f}")
        else:
            out.write("Within VCE type:" + f"{vce:>13}\n", "text")
        out.write("\n", "text")
    if not table:
        return
    se = np.sqrt(np.maximum(np.diag(est.V), 0))
    ci = f"[{lev:g}% Conf. Interval]"
    out.write("-" * 78 + "\n", "text")
    out.write(f"{est.depvar[:12]:>12} |      Coef.   Std. Err.      t    P>|t|     {ci:>20}\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    a = 1 - lev / 100
    eq_prev = None
    multi = len(set(est.eqnames)) > 1
    for i, nm in enumerate(est.names):
        if multi and est.eqnames[i] != eq_prev:
            if eq_prev is not None:
                out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
            out.write(f"{est.eqnames[i][:12]:<12} |\n", "result")
            eq_prev = est.eqnames[i]
        lab = nm
        out.write(f"{lab[:12]:>12} |", "text")
        b, sd, d = est.b[i], se[i], X["df"][i]
        if sd <= 0 or not np.isfinite(sd):
            out.write("          0  (omitted)\n" if "o." in nm or "b." in nm else
                      f"  {g9(b):>9}          .        .       .            .           .\n", "result")
            continue
        tv = b / sd
        dd = d if np.isfinite(d) else 1e12
        p = float(2 * st.t.sf(abs(tv), dd))
        c = float(st.t.ppf(1 - a / 2, dd))
        out.write(f"  {g9(b):>9}  {g9(sd):>9}  {tv:>7.2f}  {p:>6.3f}    {g9(b - c * sd):>9}   {g9(b + c * sd):>9}\n",
                  "result")
    out.write("-" * 78 + "\n", "text")
