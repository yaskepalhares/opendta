"""svar: short-run structural vector autoregression, A e_t = B u_t ([TS] svar). A and B are given
as matrices with numbers for fixed entries and missing values for the parameters to estimate."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import match_options, parse_standard
from .results import CoefRow, Estimates, coef_table, post

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def svar_loglik(Sigma: np.ndarray, T: int, A: np.ndarray, B: np.ndarray) -> float:
    """Log likelihood of the A-B model for the residual covariance Sigma:
    -(T/2) [K ln(2 pi) - ln|A|^2 + ln|B|^2 + tr(A'(B^-1)'(B^-1) A Sigma)]."""
    K = Sigma.shape[0]
    sa, la = np.linalg.slogdet(A)
    sb, lb = np.linalg.slogdet(B)
    if sa == 0 or sb == 0:
        return -np.inf
    P = np.linalg.solve(B, A)
    return -0.5 * T * (K * math.log(2 * math.pi) - 2 * la + 2 * lb + float(np.trace(P.T @ P @ Sigma)))


def unpack(phi: np.ndarray, Afix: np.ndarray, Bfix: np.ndarray):
    """A and B from the vector of free parameters (A first, then B, column-major as in Stata)."""
    A, B = Afix.copy(), Bfix.copy()
    fa = np.argwhere(np.isnan(Afix.T))      # row-major over the transpose = column-major
    fb = np.argwhere(np.isnan(Bfix.T))
    q = 0
    for c, r in fa:
        A[r, c] = phi[q]
        q += 1
    for c, r in fb:
        B[r, c] = phi[q]
        q += 1
    return A, B, [(r, c) for c, r in fa], [(r, c) for c, r in fb]


def start_values(Sigma: np.ndarray, Afix: np.ndarray, Bfix: np.ndarray) -> np.ndarray:
    K = Sigma.shape[0]
    A0, B0 = Afix.copy(), Bfix.copy()
    for i in range(K):
        for j in range(K):
            if np.isnan(A0[i, j]):
                A0[i, j] = 1.0 if i == j else 0.0
            if np.isnan(B0[i, j]):
                B0[i, j] = math.sqrt(Sigma[i, i]) if i == j else 0.0
    out = []
    for c in range(K):
        for r in range(K):
            if np.isnan(Afix[r, c]):
                out.append(A0[r, c])
    for c in range(K):
        for r in range(K):
            if np.isnan(Bfix[r, c]):
                out.append(B0[r, c])
    return np.array(out)


def fit_svar(Sigma: np.ndarray, T: int, Afix: np.ndarray, Bfix: np.ndarray, *, log=None, maxiter: int = 200,
             cns=None):
    """ML estimates of the free entries of A and B. Returns (phi, ll, V, converged). cns: (Tm, a) of
    linear constraints on the free entries, phi = Tm phic + a (core/constraints.py)."""
    from .ml import maximize
    from .mixed import _numderiv

    def f(phi):
        A, B, *_ = unpack(phi, Afix, Bfix)
        v = svar_loglik(Sigma, T, A, B)
        return v if np.isfinite(v) else -1e300

    def fun(phi, todo):
        if todo == 0:
            return f(phi), None, None
        ll, g, H = _numderiv(f, phi)
        return ll, g, H

    phi0 = start_values(Sigma, Afix, Bfix)
    if cns is not None:
        Tm, a = cns
        phi0 = a + Tm @ (Tm.T @ (phi0 - a))
        res = maximize(lambda pc, todo: _cfun(fun, Tm, a, pc, todo), Tm.T @ (phi0 - a), maxiter=maxiter)
        Hc = _cfun(fun, Tm, a, res.b, 2)[2]
        res.b = Tm @ res.b + a
        try:
            V = Tm @ np.linalg.inv(-Hc) @ Tm.T
        except np.linalg.LinAlgError:
            V = np.full((len(a), len(a)), np.nan)
        return res, V
    res = maximize(fun, phi0, maxiter=maxiter)
    H = fun(res.b, 2)[2]
    try:
        V = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        V = np.full_like(H, np.nan)
    return res, V


def _cfun(fun, Tm, a, pc, todo):
    ll, g, H = fun(Tm @ pc + a, todo)
    if g is None:
        return ll, g, H
    return ll, Tm.T @ g, Tm.T @ H @ Tm


def _matrix(s: "Session", name: str, K: int, what: str) -> np.ndarray:
    from ..commands.matrix import get_matrix
    m = get_matrix(s, name).data.copy()
    if m.shape != (K, K):
        raise StataError(503, f"{what} matrix must be {K} x {K}")    # VERIFY
    m[m >= SYS] = np.nan
    return m


def cmd_svar(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from .ml import print_log
    from .regress import cmdline
    from .tsmodels import _check_gaps, _lags, _sample_text, _var_fit, _var_stats
    from .tsops import require_ts
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "svar", args)
        return
    p = parse_standard(args)
    o = match_options(p.options, {"lags": 3, "exog": 3, "noconstant": 4, "dfk": 3, "small": 2,
                                  "level": 1, "aeq": 2, "beq": 2, "aconstraints": 2,
                                  "bconstraints": 2, "lrconstraints": 2, "cholesky": 4,
                                  "nolog": 3, "iterate": 4}) if p.options.strip() else {}
    if o.get("lrconstraints"):
        raise StataError(198, "option lrconstraints() not implemented in OpenDTA")   # OpenDTA
    deps = s.expand_varlist(p.varlist)
    K = len(deps)
    lags = _lags(o.get("lags")) if o.get("lags") else [1, 2]
    exog = s.expand_varlist(o["exog"]) if isinstance(o.get("exog"), str) else []
    constant = not o.get("noconstant")
    acn, bcn = o.get("aconstraints"), o.get("bconstraints")
    if not (o.get("aeq") or o.get("beq") or o.get("cholesky") or acn or bcn):
        raise StataError(198, "option aeq() or beq() required")   # VERIFY
    # a matrix given only through constraints has every entry free (VERIFY)
    Afix = _matrix(s, str(o["aeq"]), K, "A") if isinstance(o.get("aeq"), str) else \
        (np.full((K, K), np.nan) if acn else np.where(np.eye(K) == 1, 1.0, 0.0))
    Bfix = _matrix(s, str(o["beq"]), K, "B") if isinstance(o.get("beq"), str) else \
        (np.full((K, K), np.nan) if bcn else np.where(np.eye(K) == 1, 1.0, 0.0))
    if o.get("cholesky"):
        Afix = np.eye(K)
        Bfix = np.where(np.tri(K) == 1, np.nan, 0.0)
    nfree = int(np.isnan(Afix).sum() + np.isnan(Bfix).sum())
    cns, Cmat = None, None
    if acn or bcn:
        from types import SimpleNamespace
        from ..core.constraints import constraint_transform
        from .postest import _constraint
        from .sysmodels import _constraints, _numlist
        _, _, fa0, fb0 = unpack(np.zeros(nfree), Afix, Bfix)
        pnames = [f"a_{i + 1}_{j + 1}" for i, j in fa0] + [f"b_{i + 1}_{j + 1}" for i, j in fb0]
        fake = SimpleNamespace(names=["_cons"] * nfree, eqnames=pnames)
        defined = _constraints(s)
        rows_c = []
        for txt in [t for t in (acn, bcn) if t]:
            for cn_ in _numlist(str(txt)):
                if cn_ not in defined:
                    raise StataError(111, f"constraint {cn_} not found")
                coef, rhs = _constraint(s, fake, defined[cn_])
                row = np.zeros(nfree + 1)
                for j, v in coef.items():
                    row[j] += v
                row[-1] = rhs
                rows_c.append(row)
        Cmat = np.array(rows_c)
        Tm, a_ = constraint_transform(Cmat)
        cns = (Tm, a_.ravel())
    n_eff = nfree - (0 if Cmat is None else int(np.linalg.matrix_rank(Cmat[:, :-1])))
    if n_eff > K * (K + 1) // 2:
        raise StataError(498, "model not identified: too many free parameters")   # VERIFY
    Y, X, ok, xnames = _var_fit(s, deps, lags, exog, constant, touse(s, p))
    t0, t1 = _check_gaps(s, ok)
    Yo, Xo = Y[ok], X[ok]
    T = Yo.shape[0]
    mpar = Xo.shape[1]
    B_, U, Sml, ll_var, aic, hqic, sbic, fpe, det = _var_stats(Yo, Xo)
    Sig = U.T @ U / (T - mpar) if o.get("dfk") else Sml
    out = s.output
    if not o.get("nolog"):
        out.write("\nEstimating short-run parameters\n", "text")      # VERIFY
    res, Vphi = fit_svar(Sig, T, Afix, Bfix, cns=cns)
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, res)
    phi = res.b
    A, B, fa, fb = unpack(phi, Afix, Bfix)
    ll = res.ll
    ll_u = -0.5 * T * (K * (math.log(2 * math.pi) + 1) + np.linalg.slogdet(Sig)[1])
    df_oid = K * (K + 1) // 2 - n_eff
    chi_oid = 2 * (ll_u - ll)
    XtXi = np.linalg.inv(Xo.T @ Xo)
    info = require_ts(s.data)
    names = [f"a_{i + 1}_{j + 1}" for i, j in fa] + [f"b_{i + 1}_{j + 1}" for i, j in fb]
    est = Estimates("svar", " ".join(deps), names, phi, Vphi, T, "z", None)
    est.level = float(o["level"]) if o.get("level") and o["level"] is not True else float(s.settings.get("level", "95"))
    est.vcetype = ""
    rows = []
    se_all = np.sqrt(np.maximum(np.diag(Vphi), 0)) if nfree else np.zeros(0)
    q = 0
    for tag, M_, fix, free in (("a", A, Afix, fa), ("b", B, Bfix, fb)):
        if tag == "b":
            rows.append(CoefRow("sep"))
        for c in range(K):
            for r in range(K):
                nm = f"/{tag}_{r + 1}_{c + 1}"
                if np.isnan(fix[r, c]):
                    # an entry fixed by aconstraints()/bconstraints() is shown as constrained
                    kind = "cons" if Cmat is not None and se_all[q] < 1e-12 else "coef"
                    rows.append(CoefRow(kind, nm, float(M_[r, c]), float(se_all[q])))
                    q += 1
                else:
                    rows.append(CoefRow("cons", nm, float(M_[r, c])))
    est.rows = rows
    est.extra = {"N": T, "ll": ll, "sample": _sample_text(s, t0, t1), "df_oid": df_oid, "chi_oid": chi_oid,
                 "A": A, "B": B, "Afix": Afix, "Bfix": Bfix, "cns_text": []}
    if Cmat is not None:
        from .postest import _ctext
        est.extra["cns_text"] = [_ctext(fake, {j: float(v) for j, v in enumerate(row[:-1]) if v}, float(row[-1]))
                                 for row in Cmat]
    est.display = display_svar
    est.var = {"deps": deps, "lags": lags, "exog": exog, "xnames": xnames, "B": B_, "Sig": Sig,
               "XtXi": XtXi, "T": T, "U": U, "X": Xo, "Y": Yo, "Sml": Sml, "t0": t0, "t1": t1,
               "constant": constant, "small": bool(o.get("small")), "dfk": bool(o.get("dfk"))}
    est.svar = {"A": A, "B": B, "phi": phi, "V": Vphi, "Afix": Afix, "Bfix": Bfix}
    from ..commands.matrix import Matrix
    # e(): every entry of A and B in e(b) (the fixed ones with zero variance), e(Cns) with the
    # constraints in the order shown, and the results of the underlying var with the suffix _var,
    # in the order of Stata's ereturn list (compat 1915)
    allnames, b_all, free_pos = [], [], []
    q = 0
    for tag, M_, fix in (("a", A, Afix), ("b", B, Bfix)):
        for c in range(K):
            for r in range(K):
                allnames.append(f"{tag}_{r + 1}_{c + 1}")
                b_all.append(float(M_[r, c]))
                if np.isnan(fix[r, c]):
                    free_pos.append((len(allnames) - 1, q))
                    q += 1
    nall = len(allnames)
    V_all = np.zeros((nall, nall))
    for i, qi in free_pos:
        for j, qj in free_pos:
            V_all[i, j] = Vphi[qi, qj]
    est.names, est.b, est.V = ["_cons"] * nall, np.array(b_all), V_all
    est.eqnames = allnames
    texts = _cns_texts(est.extra)
    crow = []
    for t_ in texts:
        m_ = re.match(r"^\[(\w+)\]_cons = (.+)$", t_)
        row = np.zeros(nall + 1)
        if m_ and m_.group(1) in allnames:
            row[allnames.index(m_.group(1))] = 1.0
            row[-1] = float(m_.group(2))
        crow.append(row)
    Cm = Matrix(np.array(crow).reshape(-1, nall + 1), [f"r{i}" for i in range(1, len(crow) + 1)],
                ["_cons"] * nall + ["_r"])
    Cm.coleq = allnames + ["_Cns"]
    T_var = Yo.shape[0]
    stats_eq = []
    for j in range(K):
        e_ = U[:, j]
        ssr = float(e_ @ e_)
        yc = Yo[:, j] - Yo[:, j].mean()
        r2 = 1 - ssr / float(yc @ yc) if constant else 1 - ssr / float(Yo[:, j] @ Yo[:, j])
        rmse = math.sqrt(ssr / (T_var - mpar))
        sl = [i for i, (vn, _l) in enumerate(xnames) if vn != "_cons"]
        bj = B_[sl, j]
        Vj = Sig[j, j] * XtXi[np.ix_(sl, sl)]
        chi2_j = float(bj @ np.linalg.solve(Vj, bj)) if sl else SYS
        ll_j = -0.5 * T_var * (math.log(2 * math.pi) + 1 + math.log(ssr / T_var))
        stats_eq.append((chi2_j, float(len(sl)), ll_j, r2, rmse, float(mpar), float(T_var)))
    scalars = [("mlag_var", float(max(lags))), ("fpe_var", fpe), ("aic_var", aic), ("hqic_var", hqic),
               ("sbic_var", sbic), ("tparms_var", float(K * mpar)), ("k_var", float(K * mpar)),
               ("df_eq_var", float(mpar)), ("k_aux", float(nall)), ("k_eq", float(nall)), ("k_eq_var", float(K)),
               ("k_dv_var", float(K)), ("k_dv", float(K)), ("N_cns", float(len(crow))),
               ("ic_ml", float(res.iterations)), ("rc_ml", 0.0 if res.converged else 430.0)]
    if df_oid > 0:
        scalars.append(("chi2_oid", chi_oid))
    scalars += [("oid_df", float(df_oid)), ("N", float(T)), ("rank", float(n_eff))]
    for j in range(K, 0, -1):
        c2, dfm, llj, r2, rmse, kj, nj = stats_eq[j - 1]
        scalars += [(f"chi2_{j}_var", c2), (f"df_m{j}_var", dfm), (f"ll_{j}_var", llj), (f"r2_{j}_var", r2),
                    (f"rmse_{j}_var", rmse), (f"k_{j}_var", kj), (f"obs_{j}_var", nj)]
    scalars += [("ll", ll), ("ll_var", ll_var), ("detsig_ml_var", det), ("detsig_var", float(np.linalg.det(Sig))),
                ("N_gaps_var", 0.0), ("tmin", t0), ("tmax", t1)]
    tfmt = s.data.get(info.tvar).fmt
    macros = [("cmd", "svar"), ("cmdline", cmdline("svar", args)), ("predict", "svar_p"),
              ("lags_var", " ".join(map(str, lags))), ("depvar_var", " ".join(deps)), ("eqnames_var", " ".join(deps)),
              ("endog_var", " ".join(deps))]
    if exog:
        macros.append(("exog_var", " ".join(exog)))
    macros += [("timevar", info.tvar), ("tsfmt", tfmt), ("title", "Structural vector autoregression"),
               ("cns_b", ":".join(t_ for t_ in texts if t_.startswith("[b_"))),
               ("cns_a", ":".join(t_ for t_ in texts if t_.startswith("[a_"))), ("properties", "b V")]
    bf = B_.T.ravel()
    xn = [("_cons" if vn == "_cons" else f"L{l}.{vn}" if l > 0 else vn) for vn, l in xnames]
    vcols = [f"{dv}:{c}" for dv in deps for c in xn]
    mats = [("Cns", Cm), ("Sigma", Matrix(Sig, deps, deps)),
            ("G_var", Matrix(Xo.T @ Xo / T_var, xn, xn)),
            ("bf_var", Matrix(bf.reshape(1, -1), ["y1"], vcols)), ("B", Matrix(B, deps, deps)),
            ("A", Matrix(A, deps, deps)), ("V_var", Matrix(np.kron(Sig, XtXi), vcols, vcols)),
            ("b_var", Matrix(bf.reshape(1, -1), ["y1"], vcols))]
    for nm_, m_ in mats:
        if nm_ != "Cns":
            m_.keep_eq = True
    post(s, est, scalars, macros, ok, extra_mats=mats)
    display_svar(s, est)


def _cns_texts(x: dict) -> list[str]:
    """The constraints in the order Stata lists them: those of aconstraints() and bconstraints(),
    then the fixed entries of aeq() and beq() (compat 1915)."""
    from .tsmodels import _g
    out = list(x.get("cns_text", []))
    K = x["A"].shape[0]
    for tag, fix in (("a", x["Afix"]), ("b", x["Bfix"])):
        for r in range(K):
            for c in range(K):
                if not np.isnan(fix[r, c]):
                    out.append(f"[{tag}_{r + 1}_{c + 1}]_cons = {_g(fix[r, c])}")
    return out


def display_svar(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                 level: float | None = None, **_kw) -> None:
    from .tsmodels import _g
    out = s.output
    x = est.extra
    if header:
        out.write("\nStructural vector autoregression\n\n", "text")
        # the constraints, the fixed entries of A and B among them
        for n, t_ in enumerate(_cns_texts(x), 1):
            out.write(f" ({n:>2})  {t_}\n", "text")
        out.write("\n", "text")
        kind = "Overidentified model" if x["df_oid"] > 0 else "Exactly identified model"
        out.write(f"{x['sample']:<48}{'Number of obs':<18}= {int(x['N']):>10}\n", "text")
        out.write(f"{kind:<48}{'Log likelihood':<18}= {_g(x['ll'], '%9.0g'):>10}\n\n", "text")
    if table:
        coef_table(s, "", est.rows, stat="z", df=None, level=level or est.level)
        if x["df_oid"] > 0:
            out.write(f"LR test of identifying restrictions: chi2({x['df_oid']:>2})= {x['chi_oid']:>8.3f}"
                      f"  Prob > chi2 = {float(chi2_tail(x['chi_oid'], x['df_oid'])):>6.3f}\n", "text")
