"""estat subcommands after sem ([SEM] sem postestimation): gof, eqgof, ggof, eqtest, framework,
residuals, stable, summarize, stdize, teffects, mindices, scoretests and ginvariant.

Every subcommand works on the "units" of the fit: the single group, or each group of a group()
model, with that group's model, parameter vector and covariance matrix. The score (Lagrange
multiplier) tests are computed in an extended parameterization: the constraint that is tested is
released, and the statistic is s' I^-1 s with s the gradient and I minus the Hessian of the
objective at the constrained estimates ([SEM] methods and formulas for sem). The layout of every
table is VERIFY: the manual pages give no output, and compat case 1909 settles them."""

from __future__ import annotations

import copy
import math
import re

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lazy import lazy
from .results import crit
from .sem import Model, _dmoments, _g, _standardized, estat_gof

st = lazy("scipy.stats")
SYS = M.SYSMISS
MINCHI2 = 3.8414588


def _opts(args: str, spec: dict):
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    return body.strip(), (match_options(opt, spec) if opt.strip() else {})


def _units(est) -> list[dict]:
    x = est.extra
    prob = x["prob"]
    if not x.get("groups"):
        return [dict(label=None, model=est.model, th=x["th"], V=x["V"], data=prob.data[0], g=0)]
    out = []
    for g, m in enumerate(prob.models):
        ids = prob.gmap[g]
        out.append(dict(label=x["levels"][g], model=m, th=x["th"][ids],
                        V=x["V"][np.ix_(ids, ids)], data=prob.data[g], g=g))
    return out


def _gtitle(s, est, u) -> None:
    if u["label"] is not None:
        s.output.write(f"\nGroup: {est.extra['gvar']} = {u['label']:g}\n", "text")   # VERIFY


def _fmt(opt: dict, default: str = "%9.0g") -> str:
    f = opt.get("format")
    return f.strip() if isinstance(f, str) and f.strip() else default


def _nobs(s, est):
    from ..commands.matrix import Matrix
    n = [float(d.N) for d in est.extra["prob"].data]
    return Matrix(np.array(n).reshape(-1, 1), [f"{i}" for i in range(1, len(n) + 1)], ["N"])


def _r_name(base: str, u: dict, est) -> str:
    return base if u["label"] is None else f"{base}{u['g'] + 1}"


def _nameof(model: Model, n: str) -> str:
    return n


# ---------------------------------------------------------------------------
# eqgof
# ---------------------------------------------------------------------------

def _equation_fit(model: Model, th: np.ndarray):
    mu, Sig, (A, S, m, B, Sv, mv) = model.implied(th)
    rows = []
    BS = B @ S
    for n in model.endog:
        i = model.idx[n]
        fitted, resid = Sv[i, i], S[i, i]
        cov_ye = BS[i, i]
        var_hat = fitted - 2 * cov_ye + resid
        pred = fitted - resid
        r2 = pred / fitted if fitted else SYS
        cyh = fitted - cov_ye
        mc = cyh / math.sqrt(fitted * var_hat) if fitted > 0 and var_hat > 0 else SYS
        rows.append((n, fitted, pred, resid, r2, mc, mc * mc if mc < SYS else SYS))
    endog = [model.idx[n] for n in model.endog]
    cd = SYS
    if endog:
        cd = 1 - float(np.linalg.det(S[np.ix_(endog, endog)]) / np.linalg.det(Sv[np.ix_(endog, endog)]))
    return rows, cd


def estat_eqgof(s, est, args):
    from ..commands.matrix import Matrix
    _, o = _opts(args, {"format": 3})
    fmt = _fmt(o)
    from ..core.formats import format_value
    out = s.output
    r = {"N_groups": float(len(est.extra["prob"].data))}
    mats = {}
    for u in _units(est):
        _gtitle(s, est, u)
        rows, cd = _equation_fit(u["model"], u["th"])
        spec = u["model"].spec
        out.write("\nEquation-level goodness of fit\n\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(f"{'':>15} |{'Variance':^36}|\n", "text")
        out.write(f"{'depvars':>15} |{'fitted':>12}{'predicted':>12}{'residual':>12}|"
                  f"{'R-squared':>11}{'mc':>9}{'mc2':>9}\n", "text")
        out.write("-" * 16 + "+" + "-" * 61 + "\n", "text")
        for sec, members in (("observed", [x for x in rows if x[0] in spec.observed]),
                             ("latent", [x for x in rows if x[0] in spec.latent])):
            if not members:
                continue
            out.write(f"{sec:<15} |\n", "text")
            for n, *vals in members:
                out.write(f"{n:>15} |", "text")
                out.write("".join(f"{_fv(v, fmt):>12}" for v in vals[:3]) + "|", "result")
                out.write("".join(f"{_fv(v, fmt):>10}" for v in vals[3:]) + "\n", "result")
        out.write("-" * 16 + "+" + "-" * 61 + "\n", "text")
        out.write(f"{'overall':<15} |{'':36}|{_fv(cd, fmt):>11}\n", "result")
        out.write("-" * 78 + "\n", "text")
        out.write("mc  = correlation between depvar and its prediction\n"
                  "mc2 = mc^2 is the Bentler-Raykov squared multiple correlation coefficient\n", "text")
        r[_r_name("CD", u, est)] = cd
        names = [x[0] for x in rows] + ["overall"]
        data = [list(x[1:]) for x in rows] + [[SYS] * 3 + [cd, SYS, SYS]]
        mats[_r_name("eqfit", u, est)] = Matrix(np.array(data, dtype=float), names,
                                                ["fitted", "predicted", "residual", "R2", "mc", "mc2"])
    r["nobs"] = _nobs(s, est)
    r.update(mats)
    s.r = r


def _fv(v, fmt):
    from ..core.formats import format_value
    return "." if v is None or v >= SYS else format_value(float(v), fmt, pad=False).strip()


# ---------------------------------------------------------------------------
# ggof
# ---------------------------------------------------------------------------

def estat_ggof(s, est, args):
    from ..commands.matrix import Matrix
    if not est.extra.get("groups"):
        raise StataError(322, "estat ggof is only available after sem with group()")   # VERIFY
    _, o = _opts(args, {"format": 3})
    fmt = _fmt(o)
    out = s.output
    prob = est.extra["prob"]
    out.write("\nGroup-level goodness of fit\n\n", "text")
    out.write("-" * 60 + "\n", "text")
    out.write(f"{'Group':>12} |{'N':>10}{'SRMR':>12}{'CD':>12}\n", "text")
    out.write("-" * 13 + "+" + "-" * 46 + "\n", "text")
    rows = []
    for u in _units(est):
        srmr = _srmr(u["model"], u["th"], u["data"])
        _, cd = _equation_fit(u["model"], u["th"])
        rows.append([float(u["data"].N), srmr, cd])
        out.write(f"{u['label']:>12g} |{u['data'].N:>10,}", "text")
        out.write(f"{_fv(srmr, fmt):>12}{_fv(cd, fmt):>12}\n", "result")
    out.write("-" * 60 + "\n", "text")
    s.r = {"N_groups": float(len(prob.data)), "nobs": _nobs(s, est),
           "gfit": Matrix(np.array(rows), [f"{u['label']:g}" for u in _units(est)], ["N", "SRMR", "CD"])}


def _srmr(model: Model, th: np.ndarray, data) -> float:
    mu, Sig, _ = model.implied(th)
    d = np.sqrt(np.diag(data.cov))
    ds = np.sqrt(np.diag(Sig))
    p = model.p
    R = data.cov / np.outer(d, d) - Sig / np.outer(ds, ds)
    tri = R[np.tril_indices(p)]
    dm = (data.mean - mu) / d
    return float(np.sqrt(((tri ** 2).sum() + (dm ** 2).sum()) / (p * (p + 3) / 2)))


# ---------------------------------------------------------------------------
# eqtest
# ---------------------------------------------------------------------------

def estat_eqtest(s, est, args):
    from ..commands.matrix import Matrix
    _, o = _opts(args, {"total": 3, "nosvyadjust": 5})
    out = s.output
    r = {"N_groups": float(len(est.extra["prob"].data))}
    mats = {}
    for u in _units(est):
        _gtitle(s, est, u)
        model = u["model"]
        spec = model.spec
        out.write("\nWald tests for equations\n\n", "text")
        out.write("-" * 66 + "\n", "text")
        out.write(f"{'':>16} |{'chi2':>12}{'df':>8}{'p':>12}\n", "text")
        out.write("-" * 17 + "+" + "-" * 48 + "\n", "text")
        data, names = [], []
        for sec, pool in (("observed", spec.observed), ("latent", spec.latent)):
            eqs = [y for y in model.endog if y in pool]
            body = []
            for y in eqs:
                idx = [i for i, it in enumerate(model.free) if it[0] == "path" and it[2] == y]
                if not idx:
                    continue
                b = u["th"][idx]
                Vb = u["V"][np.ix_(idx, idx)]
                w = float(b @ np.linalg.pinv(Vb) @ b)
                df = int(np.linalg.matrix_rank(Vb))
                p = float(st.chi2.sf(w, df))
                body.append((y, w, df, p))
            if body:
                out.write(f"{sec:<16} |\n", "text")
            for y, w, df, p in body:
                out.write(f"{y:>16} |", "text")
                out.write(f"{w:>12.2f}{df:>8}{p:>12.4f}\n", "result")
                data.append([w, df, p])
                names.append(y)
        out.write("-" * 66 + "\n", "text")
        out.write("H0: all coefficients of the equation, other than the intercept, are 0\n", "text")
        if data:
            mats[_r_name("test", u, est)] = Matrix(np.array(data), names, ["chi2", "df", "p"])
    r["nobs"] = _nobs(s, est)
    r.update(mats)
    s.r = r


# ---------------------------------------------------------------------------
# framework
# ---------------------------------------------------------------------------

def estat_framework(s, est, args):
    from ..commands.matrix import Matrix, list_matrix
    _, o = _opts(args, {"standardized": 5, "compact": 4, "fitted": 3, "format": 3})
    if o.get("standardized") or o.get("compact"):
        raise StataError(198, "option standardized and compact are not available yet in OpenDTA")
    fmt = _fmt(o)
    out = s.output
    r = {"N_groups": float(len(est.extra["prob"].data)), "standardized": 0.0}
    mats = {}
    for u in _units(est):
        _gtitle(s, est, u)
        model = u["model"]
        mu, Sig, (A, S, m, B, Sv, mv) = model.implied(u["th"])
        Y = [n for n in model.names if n in model.endog]
        X = [n for n in model.names if n in model.exog]
        iy = [model.idx[n] for n in Y]
        ix = [model.idx[n] for n in X]
        ey = [f"e.{n}" for n in Y]
        blocks = []
        if Y:
            blocks.append(("Beta", "Coefficients of endogenous variables on endogenous variables",
                           Matrix(A[np.ix_(iy, iy)], Y, Y)))
        if Y and X:
            blocks.append(("Gamma", "Coefficients of endogenous variables on exogenous variables",
                           Matrix(A[np.ix_(iy, ix)], Y, X)))
        if Y:
            blocks.append(("Psi", "Covariances of error variables", Matrix(S[np.ix_(iy, iy)], ey, ey)))
            blocks.append(("alpha", "Intercepts of endogenous variables",
                           Matrix(m[iy].reshape(-1, 1), Y, ["_cons"])))
        if X:
            blocks.append(("Phi", "Covariances of exogenous variables", Matrix(S[np.ix_(ix, ix)], X, X)))
            blocks.append(("kappa", "Means of exogenous variables", Matrix(m[ix].reshape(-1, 1), X, ["_cons"])))
        if o.get("fitted"):
            obs = model.spec.observed
            blocks.append(("Sigma", "Fitted covariance matrix of observed variables", Matrix(Sig, obs, obs)))
            blocks.append(("mu", "Fitted means of observed variables", Matrix(mu.reshape(-1, 1), obs, ["_cons"])))
        for name, title, mat in blocks:
            out.write(f"\n{title}\n", "text")
            list_matrix(s, name, mat, fmt, header=False)
            mats[_r_name(name, u, est)] = mat
    r["nobs"] = _nobs(s, est)
    r.update(mats)
    s.r = r


# ---------------------------------------------------------------------------
# residuals
# ---------------------------------------------------------------------------

def estat_residuals(s, est, args):
    from ..commands.matrix import Matrix, list_matrix
    _, o = _opts(args, {"normalized": 4, "standardized": 5, "sample": 3, "zerotolerance": 4, "format": 3})
    fmt = _fmt(o, "%9.3f")
    out = s.output
    kind = "standardized" if o.get("standardized") else "normalized" if o.get("normalized") else "raw"
    zt = float(o["zerotolerance"]) if o.get("zerotolerance") else 0.0
    r = {"N_groups": float(len(est.extra["prob"].data))}
    mats = {}
    for u in _units(est):
        _gtitle(s, est, u)
        model, th, data = u["model"], u["th"], u["data"]
        obs = model.spec.observed
        p = model.p
        mu, Sig, dm = _dmoments(model, th)
        S_ = data.cov * (data.N / (data.N - 1) if o.get("sample") else 1.0)
        rc = S_ - Sig
        rm = data.mean - mu
        if kind != "raw":
            N = data.N
            Vth = u["V"]
            J_mu = np.zeros((p, len(th)))
            J_sig = np.zeros((p, p, len(th)))
            for k, (dmv, dSv) in enumerate(dm):
                if dmv is not None:
                    J_mu[:, k] = dmv
                if dSv is not None:
                    J_sig[:, :, k] = dSv
            zc = np.zeros((p, p))
            zm = np.zeros(p)
            for i in range(p):
                zm[i] = rm[i] / math.sqrt(max(Sig[i, i] / N - (J_mu[i] @ Vth @ J_mu[i] if kind == "standardized" else 0), 1e-300))
                for j in range(i + 1):
                    v = (Sig[i, i] * Sig[j, j] + Sig[i, j] ** 2) / N
                    if kind == "standardized":
                        v = v - float(J_sig[i, j] @ Vth @ J_sig[i, j])
                    zc[i, j] = zc[j, i] = rc[i, j] / math.sqrt(v) if v > 0 else 0.0
            rc, rm = zc, zm
        rc = np.where(np.abs(rc) < zt, 0.0, rc)
        word = {"raw": "Raw", "normalized": "Normalized", "standardized": "Standardized"}[kind]
        out.write(f"\n{word} residuals\n", "text")
        mm = Matrix(rm.reshape(1, -1), ["r1"], obs)
        out.write("\nMean residuals\n", "text")
        list_matrix(s, "", mm, fmt, header=False)
        out.write("\nCovariance residuals\n", "text")
        list_matrix(s, "", Matrix(rc, obs, obs), fmt, header=False)
        mats[_r_name("res_mean", u, est)] = mm
        mats[_r_name("res_cov", u, est)] = Matrix(rc, obs, obs)
    r["nobs"] = _nobs(s, est)
    r.update(mats)
    s.r = r


# ---------------------------------------------------------------------------
# stable
# ---------------------------------------------------------------------------

def estat_stable(s, est, args):
    from ..commands.matrix import Matrix
    _, o = _opts(args, {"detail": 3, "format": 3})
    out = s.output
    r = {"N_groups": float(len(est.extra["prob"].data))}
    for u in _units(est):
        _gtitle(s, est, u)
        model = u["model"]
        Y = [n for n in model.names if n in model.endog]
        iy = [model.idx[n] for n in Y]
        A = model.implied(u["th"])[2][0][np.ix_(iy, iy)]
        ev = np.linalg.eigvals(A)
        mod = np.abs(ev)
        order = np.argsort(-mod, kind="stable")
        ev, mod = ev[order], mod[order]
        out.write("\nStability analysis of simultaneous equation systems\n\n", "text")
        out.write("Eigenvalue stability condition\n", "text")
        out.write("-" * 41 + "\n", "text")
        out.write(f"{'Eigenvalue':>24} |{'Modulus':>14}\n", "text")
        out.write("-" * 25 + "+" + "-" * 15 + "\n", "text")
        for e, md in zip(ev, mod):
            sign = "+" if e.imag >= 0 else "-"
            out.write(f"{e.real:>12.4f} {sign} {abs(e.imag):.4f}i |", "text")
            out.write(f"{md:>14.4f}\n", "result")
        out.write("-" * 41 + "\n", "text")
        idx = float(mod.max()) if len(mod) else 0.0
        out.write(f"\nStability index = {idx:.4f}\n", "text")
        out.write("The stability index is the largest modulus; the system is stable if it is "
                  "less than 1.\n", "text")
        r["stindex" if u["label"] is None else f"stindex{u['g'] + 1}"] = idx
        r[_r_name("Re", u, est)] = Matrix(ev.real.reshape(-1, 1), [f"r{i}" for i in range(1, len(ev) + 1)], ["Re"])
        r[_r_name("Im", u, est)] = Matrix(ev.imag.reshape(-1, 1), [f"r{i}" for i in range(1, len(ev) + 1)], ["Im"])
        r[_r_name("Modulus", u, est)] = Matrix(mod.reshape(-1, 1), [f"r{i}" for i in range(1, len(mod) + 1)], ["Modulus"])
    s.r = r


# ---------------------------------------------------------------------------
# summarize
# ---------------------------------------------------------------------------

def estat_summarize(s, est, args):
    from .estat import _summarize
    body, o = _opts(args, {"group": 3, "labels": 3, "noheader": 5, "noweights": 5})
    names = s.expand_varlist(body) if body else list(est.model.spec.observed)
    _summarize(s, est, " ".join(names))


# ---------------------------------------------------------------------------
# stdize
# ---------------------------------------------------------------------------

def estat_stdize(s, est, args):
    """estat stdize: command -- runs the test, lincom or similar on the standardized parameters
    (b and V are replaced by the standardized estimates and their delta-method covariance for
    the duration of the command)."""
    cmd = args.strip()
    if cmd.startswith(":"):
        cmd = cmd[1:].strip()
    if not cmd:
        raise StataError(198, "command required")   # VERIFY
    if est.extra.get("groups"):
        raise StataError(198, "estat stdize after sem with group() is not available yet in OpenDTA")
    model = est.model
    th, V = est.extra["th"], est.extra["V"]
    bstd = _standardized(model, th)
    J = np.zeros((len(th), len(th)))
    for i in range(len(th)):
        e = np.zeros(len(th))
        h = 1e-6 * (abs(th[i]) + 1)
        e[i] = h
        J[:, i] = (_standardized(model, th + e) - _standardized(model, th - e)) / (2 * h)
    Vs = J @ V @ J.T
    keep = [i for i, it in enumerate(model.free) if it not in model.hidden]
    saved = (est.b, est.V, s.e["b"].data.copy(), s.e["V"].data.copy())
    try:
        est.b = bstd[keep]
        est.V = Vs[np.ix_(keep, keep)]
        s.e["b"].data = est.b.reshape(1, -1).copy()
        s.e["V"].data = est.V.copy()
        s.interp.execute(cmd)
    finally:
        est.b, est.V = saved[0], saved[1]
        s.e["b"].data, s.e["V"].data = saved[2], saved[3]


def _std_names(est):
    return est.names


# ---------------------------------------------------------------------------
# teffects
# ---------------------------------------------------------------------------

def estat_teffects(s, est, args):
    _, o = _opts(args, {"compact": 4, "standardized": 5, "nodirect": 5, "noindirect": 5, "nototal": 5})
    if o.get("standardized"):
        raise StataError(198, "option standardized is not available yet in OpenDTA")
    out = s.output
    lv = est.level
    c = crit("z", lv, None)
    for u in _units(est):
        _gtitle(s, est, u)
        model, th, V = u["model"], u["th"], u["V"]
        spec = model.spec

        def effects(t):
            A = model.implied(t)[2][0]
            B = model.implied(t)[2][3]
            tot = B - np.eye(model.k)
            return A, tot, tot - A
        A0, T0, I0 = effects(th)
        k = len(th)
        Jd = {}
        for key, which in (("d", 0), ("t", 1), ("i", 2)):
            pass
        pairs = [(a, b) for b in model.names for a in model.names if a != b]

        def vec(t):
            A, T, I = effects(t)
            return A, T, I

        base = vec(th)
        J = {0: np.zeros((model.k, model.k, k)), 1: np.zeros((model.k, model.k, k)),
             2: np.zeros((model.k, model.k, k))}
        for i in range(k):
            e = np.zeros(k)
            h = 1e-6 * (abs(th[i]) + 1)
            e[i] = h
            up, dn = vec(th + e), vec(th - e)
            for w in range(3):
                J[w][:, :, i] = (up[w] - dn[w]) / (2 * h)
        tables = []
        if not o.get("nodirect"):
            tables.append(("Direct effects", 0))
        if not o.get("noindirect"):
            tables.append(("Indirect effects", 2))
        if not o.get("nototal"):
            tables.append(("Total effects", 1))
        for title, w in tables:
            M_ = base[w]
            out.write(f"\n{title}\n\n", "text")
            W = 13
            out.write("-" * (W + 66) + "\n", "text")
            out.write(f"{'':>{W}} |                 OIM\n", "text")
            out.write(f"{'':>{W}} |      Coef.   Std. Err.      z    P>|z|     [{lv:g}% Conf. Interval]\n", "text")
            out.write("-" * (W + 1) + "+" + "-" * 64 + "\n", "text")
            eq_names = [n for n in model.names if n in model.endog]
            meas = [n for n in eq_names if n in spec.observed
                    and all(a in spec.latent for a, b in spec.paths if b == n)]
            for sec, members in (("Structural", [n for n in eq_names if n not in meas]),
                                 ("Measurement", meas)):
                rows = []
                for y in members:
                    iy = model.idx[y]
                    srcs = [a for a in model.names if a != y and abs(M_[iy, model.idx[a]]) > 1e-12]
                    if srcs:
                        rows.append((y, srcs))
                if not rows:
                    continue
                out.write(f"{sec:<{W}} |\n", "text")
                for y, srcs in rows:
                    out.write(f"  {y + ' <-':<{W - 2}} |\n", "text")
                    for a in srcs:
                        iy, ia = model.idx[y], model.idx[a]
                        b = M_[iy, ia]
                        se = math.sqrt(max(float(J[w][iy, ia] @ V @ J[w][iy, ia]), 0))
                        if se > 0:
                            z = b / se
                            pv = 2 * float(st.norm.sf(abs(z)))
                            out.write(f"{a:>{W}} |  {_g(b):>9}  {_g(se):>9}  {z:>7.2f}  {pv:>6.3f}    "
                                      f"{_g(b - c * se):>9}   {_g(b + c * se):>9}\n", "result")
                        else:
                            out.write(f"{a:>{W}} |  {_g(b):>9}  {'.':>9}\n", "result")
            out.write("-" * (W + 66) + "\n", "text")


# ---------------------------------------------------------------------------
# score tests
# ---------------------------------------------------------------------------

def _extended(prob, split=frozenset(), split_all=frozenset()):
    from .semfit import Problem
    return Problem(prob.models, prob.data, prob.method, prob.inv, split=set(split), split_all=set(split_all))


def _theta_in(prob, th_glob, ext):
    t = np.zeros(ext.nglob)
    for g in range(len(prob.models)):
        t[ext.gmap[g]] = th_glob[prob.gmap[g]]
    return t


def _lm(prob, th_glob, ext) -> tuple[float, int]:
    t = _theta_in(prob, th_glob, ext)
    sc = ext.grad(t)
    info = -ext.hess(t)
    stat = float(sc @ np.linalg.pinv(info) @ sc)
    return max(stat, 0.0), ext.nglob - prob.nglob


def _th_glob(est):
    x = est.extra
    return x.get("th_glob", x["th"])


def _label(model: Model, it) -> str:
    kind, a, b = it
    if kind == "path":
        return f"{b} <- {a}"
    if kind == "mean":
        return f"mean({a})" if a in model.exog else f"{a} <- _cons"
    if kind == "var":
        return f"var(e.{a})" if a in model.endog else f"var({a})"
    ea = f"e.{a}" if a in model.endog else a
    eb = f"e.{b}" if b in model.endog else b
    return f"cov({ea},{eb})"


def estat_scoretests(s, est, args):
    from ..commands.matrix import Matrix
    _, o = _opts(args, {"minchi2": 3})
    minchi2 = float(o["minchi2"]) if o.get("minchi2") else MINCHI2
    prob = est.extra["prob"]
    th = _th_glob(est)
    out = s.output
    cons = []      # (text, ext problem)
    for g, m in enumerate(prob.models):
        groups: dict[str, list] = {}
        for it, lab in m.ties.items():
            groups.setdefault(lab, []).append(it)
        for lab, items in groups.items():
            for it in items[1:]:
                text = f"[{items[0][2]}]{items[0][1]} - [{it[2]}]{it[1]} = 0"
                cons.append((text, _extended(prob, split={(g, it)})))
    out.write("\nScore tests for linear constraints\n\n", "text")
    rows = []
    for k, (text, ext) in enumerate(cons, 1):
        out.write(f" ({k:2d}) {text}\n", "text")
        stat, df = _lm(prob, th, ext)
        rows.append((k, stat, max(df, 1)))
    shown = [r for r in rows if r[1] >= minchi2]
    if not shown:
        out.write(f"\n(no score tests to report; all chi2 values less than {minchi2})\n", "text")
    else:
        out.write(f"\n{'':>8}{'chi2':>10}{'df':>6}{'P>chi2':>9}\n", "text")
        for k, stat, df in shown:
            out.write(f" ({k:2d}) ", "text")
            out.write(f"{stat:>9.3f}{df:>6}{float(st.chi2.sf(stat, df)):>9.2f}\n", "result")
    s.r = {"N_groups": float(len(prob.data)), "nobs": _nobs(s, est)}
    if rows:
        s.r["Cns_sctest"] = Matrix(np.array([[r[1], r[2], float(st.chi2.sf(r[1], r[2]))] for r in rows]),
                                   [f"c{r[0]}" for r in rows], ["chi2", "df", "p"])


# ---------------------------------------------------------------------------
# modification indices
# ---------------------------------------------------------------------------

def _augment(model: Model, item):
    kind, a, b = item
    spec = copy.deepcopy(model.spec)
    covs = list(model.covs_in)
    if kind == "path":
        spec.paths.append((a, b))
    else:
        covs.append((a, b))
    return Model(spec, covs, dict(model.var_fixed))


def estat_mindices(s, est, args):
    from ..commands.matrix import Matrix
    from .semfit import Problem, item_class, parse_ginvariant
    _, o = _opts(args, {"showpclass": 5, "minchi2": 3})
    minchi2 = float(o["minchi2"]) if o.get("minchi2") else MINCHI2
    want = parse_ginvariant("all" if str(o.get("showpclass") or "all").strip() == "all"
                            else str(o["showpclass"]))
    if est.extra.get("groups"):
        raise StataError(198, "estat mindices after sem with group() is not available yet in OpenDTA")
    prob = est.extra["prob"]
    model = est.model
    spec = model.spec
    th = est.extra["th"]
    data = prob.data[0]
    mu, Sig, (A, S, m, B, Sv, mv) = model.implied(th)
    sd = np.sqrt(np.maximum(np.diag(Sv), 1e-300))
    cands = []
    exo_obs = {n for n in model.exog if n in spec.observed}
    for to in model.names:
        if to not in model.endog:
            continue
        for fr in model.names:
            if fr == to or (fr, to) in spec.paths:
                continue
            cands.append(("path", fr, to))
    err = [n for n in model.names if n in model.endog]
    have = {tuple(sorted(c, key=lambda n: model.idx[n])) for c in model.covs_in}
    for i, a in enumerate(err):
        for b in err[i + 1:]:
            key = tuple(sorted((a, b), key=lambda n: model.idx[n]))
            if key in have:
                continue
            cands.append(("cov", key[0], key[1]))
    rows = []
    for it in cands:
        try:
            aug = _augment(model, it)
        except StataError:
            continue
        cls = item_class(aug, it)
        if cls not in want:
            continue
        th_aug = np.array([th[model.free.index(x)] if x in model.free else 0.0 for x in aug.free])
        sub = Problem([aug], [data], prob.method, set())
        # the observed exogenous variables stay unconstrained as in the fitted model
        gr = sub.grad(th_aug)
        info = -sub.hess(th_aug)
        p = aug.free.index(it)
        c = [i for i in range(len(th_aug)) if i != p]
        ieff = info[p, p] - info[p, c] @ np.linalg.pinv(info[np.ix_(c, c)]) @ info[c, p]
        if ieff <= 1e-12:
            continue
        mi = gr[p] ** 2 / ieff
        epc = gr[p] / ieff
        if it[0] == "path":
            sepc = epc * sd[model.idx[it[1]]] / sd[model.idx[it[2]]]
        else:
            va, vb = S[model.idx[it[1]], model.idx[it[1]]], S[model.idx[it[2]], model.idx[it[2]]]
            sepc = epc / math.sqrt(va * vb) if va > 0 and vb > 0 else SYS
        rows.append((it, cls, mi, epc, sepc))
    shown = sorted([r for r in rows if r[2] >= minchi2], key=lambda r: (r[0][0] != "path", -r[2]))
    out = s.output
    out.write("\nModification indices\n\n", "text")
    if not shown:
        out.write(f"(no modification indices to report; all chi2 values less than {minchi2})\n", "text")
    else:
        out.write("-" * 78 + "\n", "text")
        out.write(f"{'':>26} |{'':>10}{'':>6}{'':>9}{'':>12}{'Standard':>12}\n", "text")
        out.write(f"{'':>26} |{'MI':>10}{'df':>6}{'P>MI':>9}{'EPC':>12}{'EPC':>12}\n", "text")
        out.write("-" * 27 + "+" + "-" * 50 + "\n", "text")
        for it, cls, mi, epc, sepc in shown:
            lab = _label(model, it)
            out.write(f"{lab:<26} |", "text")
            out.write(f"{mi:>10.3f}{1:>6}{float(st.chi2.sf(mi, 1)):>9.2f}{epc:>12.4f}{_fv(sepc, '%9.0g'):>12}\n",
                      "result")
        out.write("-" * 78 + "\n", "text")
        out.write("EPC = expected parameter change\n", "text")
    classes = ["scoef", "scons", "mcoef", "mcons", "serrvar", "merrvar", "smerrcov", "meanex", "covex"]
    r = {"N_groups": 1.0, "nobs": _nobs(s, est)}
    if shown:
        r["mindices_pclass"] = Matrix(np.array([[classes.index(x[1]) + 1] for x in shown], dtype=float),
                                      [_label(model, x[0]) for x in shown], ["pclass"])
        r["mindices"] = Matrix(np.array([[x[2], 1, float(st.chi2.sf(x[2], 1)), x[3], x[4]] for x in shown]),
                               [_label(model, x[0]) for x in shown], ["MI", "df", "P>MI", "EPC", "stdEPC"])
    s.r = r


# ---------------------------------------------------------------------------
# ginvariant
# ---------------------------------------------------------------------------

def estat_ginvariant(s, est, args):
    from ..commands.matrix import Matrix
    from .semfit import PCLASSES, parse_ginvariant
    if not est.extra.get("groups"):
        raise StataError(322, "estat ginvariant is only available after sem with group()")   # VERIFY
    _, o = _opts(args, {"showpclass": 5, "class": 3, "legend": 3})
    show = parse_ginvariant("all" if str(o.get("showpclass") or "all").strip() == "all"
                            else str(o["showpclass"]))
    prob = est.extra["prob"]
    th = _th_glob(est)
    Vg = est.extra["V"]
    G = len(prob.models)
    m0 = prob.models[0]
    out = s.output
    # items present in every group, by class
    items = [it for it in m0.free if all(it in m.free for m in prob.models)
             and not (it in est.extra["hidden"])]
    res = []
    for it in items:
        i0 = m0.free.index(it)
        cls = prob.classes[0][i0]
        if cls not in show:
            continue
        ids = [prob.gmap[g][prob.models[g].free.index(it)] for g in range(G)]
        wald = (np.nan, 0, np.nan)
        score = (np.nan, 0, np.nan)
        if len(set(ids)) == 1:
            ext = _extended(prob, split_all={it})
            stat, df = _lm(prob, th, ext)
            score = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        else:
            R = np.zeros((G - 1, prob.nglob))
            for g in range(1, G):
                R[g - 1, ids[g]] += 1
                R[g - 1, ids[0]] -= 1
            d = R @ th
            Vd = R @ Vg @ R.T
            stat = float(d @ np.linalg.pinv(Vd) @ d)
            df = int(np.linalg.matrix_rank(Vd))
            wald = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        res.append((it, cls, wald, score))
    out.write("\nTests for group invariance of parameters\n\n", "text")
    if o.get("class"):
        _class_tests(s, est, prob, th, Vg, show)
    else:
        W = 28
        out.write("-" * (W + 52) + "\n", "text")
        out.write(f"{'':>{W}} |{'Wald Test':^24}|{'Score Test':^24}\n", "text")
        out.write(f"{'':>{W}} |{'chi2':>8}{'df':>6}{'p':>9} |{'chi2':>8}{'df':>6}{'p':>9}\n", "text")
        out.write("-" * (W + 1) + "+" + "-" * 51 + "\n", "text")
        for it, cls, wald, score in res:
            lab = _label(m0, it)
            out.write(f"{lab:<{W}} |", "text")
            out.write(_tcell(wald) + " |" + _tcell(score) + "\n", "result")
        out.write("-" * (W + 52) + "\n", "text")
        out.write("Wald tests are for parameters that vary across groups; score tests for those "
                  "constrained equal.\n", "text")
    r = {"N_groups": float(G), "nobs": _nobs(s, est)}
    if res:
        r["test"] = Matrix(np.array([[x[2][0], x[2][1], x[2][2], x[3][0], x[3][1], x[3][2]] for x in res],
                                    dtype=float), [_label(m0, x[0]) for x in res],
                           ["Wald_chi2", "Wald_df", "Wald_p", "score_chi2", "score_df", "score_p"])
        r["test_pclass"] = Matrix(np.array([[PCLASSES.index(x[1]) + 1] for x in res], dtype=float),
                                  [_label(m0, x[0]) for x in res], ["pclass"])
    s.r = r


def _tcell(t) -> str:
    stat, df, p = t
    if df == 0 or not np.isfinite(stat):
        return f"{'.':>8}{'.':>6}{'.':>9}"
    return f"{stat:>8.3f}{df:>6}{p:>9.3f}"


def _class_tests(s, est, prob, th, Vg, show):
    from .semfit import PCLASSES
    out = s.output
    G = len(prob.models)
    m0 = prob.models[0]
    out.write("-" * 76 + "\n", "text")
    out.write(f"{'Class':>14} |{'Wald Test':^30}|{'Score Test':^30}\n", "text")
    out.write(f"{'':>14} |{'chi2':>10}{'df':>6}{'p':>9}   |{'chi2':>10}{'df':>6}{'p':>9}\n", "text")
    out.write("-" * 15 + "+" + "-" * 60 + "\n", "text")
    for cls in PCLASSES:
        if cls not in show:
            continue
        its = [it for i, it in enumerate(m0.free) if prob.classes[0][i] == cls
               and all(it in m.free for m in prob.models) and it not in est.extra["hidden"]]
        if not its:
            continue
        constrained = cls in prob.inv
        wald = score = (np.nan, 0, np.nan)
        if constrained:
            ext = _extended(prob, split_all=set(its))
            stat, df = _lm(prob, th, ext)
            score = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        else:
            rows = []
            for it in its:
                ids = [prob.gmap[g][prob.models[g].free.index(it)] for g in range(G)]
                for g in range(1, G):
                    r = np.zeros(prob.nglob)
                    r[ids[g]] += 1
                    r[ids[0]] -= 1
                    rows.append(r)
            R = np.array(rows)
            d = R @ th
            Vd = R @ Vg @ R.T
            stat = float(d @ np.linalg.pinv(Vd) @ d)
            df = int(np.linalg.matrix_rank(Vd))
            wald = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        out.write(f"{cls:>14} |", "text")
        out.write(_tcell(wald).replace("", "") + "   |" + _tcell(score) + "\n", "result")
    out.write("-" * 76 + "\n", "text")


ESTAT = {
    "gof": (3, estat_gof),
    "eqgof": (4, estat_eqgof),
    "ggof": (2, estat_ggof),
    "eqtest": (4, estat_eqtest),
    "framework": (3, estat_framework),
    "residuals": (3, estat_residuals),
    "stable": (3, estat_stable),
    "summarize": (3, estat_summarize),
    "stdize": (4, estat_stdize),
    "teffects": (2, estat_teffects),
    "mindices": (2, estat_mindices),
    "scoretests": (2, estat_scoretests),
    "ginvariant": (2, estat_ginvariant),
}
