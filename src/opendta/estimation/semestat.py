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
        out.write(f"{'':>12} |{'':13}Variance{'':12}|\n", "text")
        out.write(f"{'depvars':>12} |{'fitted':>10}{'predicted':>11}{'residual':>11} |"
                  f"{'R-squared':>10}{'mc':>10}{'mc2':>9}\n", "text")
        out.write("-" * 13 + "+" + "-" * 33 + "+" + "-" * 30 + "\n", "text")
        for sec, members in (("observed", [x for x in rows if x[0] in spec.observed]),
                             ("latent", [x for x in rows if x[0] in spec.latent])):
            if not members:
                continue
            out.write(f"{sec:<12} |{'':33}|\n", "text")
            for n, *vals in members:
                out.write(f"{n[-12:]:>12} |", "text")
                out.write(f"{_fv(vals[0], fmt):>10}" + "".join(f"{_fv(v, fmt):>11}" for v in vals[1:3]) + " |", "result")
                out.write("".join(f"{_fv(v, fmt):>10}" for v in vals[3:]) + "\n", "result")
        out.write("-" * 13 + "+" + "-" * 33 + "+" + "-" * 30 + "\n", "text")
        out.write(f"{'overall':>12} |{'':33}|{_fv(cd, fmt):>10}\n", "result")
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
    out.write("\nGroup-level fit statistics\n", "text")
    out.write("-" * 44 + "\n", "text")
    out.write(f"{'':12} |{'N':>8}{'SRMR':>11}{'CD':>11}\n", "text")
    out.write("-" * 13 + "+" + "-" * 30 + "\n", "text")
    out.write(f"{est.extra['gvar'][-12:]:>12} |\n", "text")
    rows = []
    for u in _units(est):
        srmr = _srmr(u["model"], u["th"], u["data"])
        _, cd = _equation_fit(u["model"], u["th"])
        rows.append([float(u["data"].N), srmr, cd])
        out.write(f"{u['label']:>11g}  |", "text")
        out.write(f"{u['data'].N:>8}{_fv(srmr, '%9.3f'):>11}{_fv(cd, '%9.3f'):>11}\n", "result")
    out.write("-" * 44 + "\n", "text")
    if prob.inv:
        # VERIFY: the chi-squared columns of a model without constraints between groups
        out.write("Note: Group-level chi-squared statistics\n      are not reported because of\n"
                  "      constraints between groups.\n", "text")
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
        out.write("-" * 38 + "\n", "text")
        out.write(f"{'':>12} |{'chi2':>9}{'df':>5}{'p':>10}\n", "text")
        out.write("-" * 13 + "+" + "-" * 24 + "\n", "text")
        data, names = [], []
        for sec, pool in (("observed", spec.observed), ("latent", spec.latent)):
            eqs = [y for y in model.endog if y in pool]
            body = []
            for y in eqs:
                idx = [i for i, it in enumerate(model.free) if it[0] == "path" and it[2] == y]
                if idx:
                    b = u["th"][idx]
                    Vb = u["V"][np.ix_(idx, idx)]
                    w = float(b @ np.linalg.pinv(Vb) @ b)
                    df = int(np.linalg.matrix_rank(Vb))
                    p = float(st.chi2.sf(w, df))
                else:
                    # an equation whose coefficients are all constrained: chi2 0 on 0 df
                    w, df, p = 0.0, 0, SYS
                body.append((y, w, df, p))
            if body:
                out.write(f"{sec:<12} |\n", "text")
            for y, w, df, p in body:
                out.write(f"{y[-12:]:>12} |", "text")
                ps = "." if p >= SYS else f"{p:.4f}"
                out.write(f"{w:>9.2f}{df:>5}{ps:>10}\n", "result")
                data.append([w, df, p])
                names.append(y)
        out.write("-" * 38 + "\n", "text")
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
        obs = set(model.spec.observed)
        grp = lambda n: "observed" if n.split(".")[-1] in obs else "latent"     # noqa: E731
        blocks = []
        if Y:
            blocks.append(("Beta", "Endogenous variables on endogenous variables",
                           Matrix(A[np.ix_(iy, iy)], Y, Y), Y, Y))
        if Y and X:
            blocks.append(("Gamma", "Exogenous variables on endogenous variables",
                           Matrix(A[np.ix_(iy, ix)], Y, X), Y, X))
        if Y:
            blocks.append(("Psi", "Covariances of error variables", Matrix(S[np.ix_(iy, iy)], ey, ey), ey, ey))
            blocks.append(("alpha", "Intercepts of endogenous variables",
                           Matrix(m[iy].reshape(1, -1), ["_cons"], Y), None, Y))
        if X:
            blocks.append(("Phi", "Covariances of exogenous variables", Matrix(S[np.ix_(ix, ix)], X, X), X, X))
            blocks.append(("kappa", "Means of exogenous variables", Matrix(m[ix].reshape(1, -1), ["mean"], X),
                           None, X))
        if o.get("fitted"):
            ob = model.spec.observed
            blocks.append(("Sigma", "Fitted covariances of observed variables", Matrix(Sig, ob, ob), ob, ob))
            blocks.append(("mu", "Fitted means of observed variables", Matrix(mu.reshape(1, -1), ["mu"], ob),
                           None, ob))
        for name, title, mat, rows, cols in blocks:
            out.write(f"\n{title}\n\n", "text")
            _fw_table(out, name, mat, rows, cols, grp, fmt)
            mats[_r_name(name, u, est)] = mat
    r["nobs"] = _nobs(s, est)
    r.update(mats)
    s.r = r


def _fw_table(out, name, mat, rows, cols, grp, fmt, lower=None):
    """One matrix of estat framework: column and row groups (observed, latent), lower triangle
    when the matrix is lower triangular (VERIFY: a nonrecursive Beta)."""
    D = np.asarray(mat.data, float)
    k = D.shape[1]
    width = 10 + 11 * (k - 1) + 1
    if lower is None:
        lower = D.shape[0] == k and rows is not None and (name in ("Psi", "Phi", "Sigma") or np.allclose(np.triu(D, 1), 0))
    if grp is not None:
        cg = [grp(c) for c in cols]
        head = ""
        for j, c in enumerate(cg):
            if j == 0 or c != cg[j - 1]:
                pos = 0 if j == 0 else 10 + 11 * (j - 1)
                head = head.ljust(pos) + " " + c
        out.write(f"{'':16} |{head:<{width}}\n", "text")
    out.write(f"{name:>16} |" + f"{cols[0][-10:]:>10}" + "".join(f"{c[-10:]:>11}" for c in cols[1:]) + " \n", "text")
    out.write("    " + "-" * 13 + "+" + "-" * width + "\n", "text")
    labels = rows if rows is not None else mat.rownames
    last = None
    for i, lab in enumerate(labels):
        if rows is not None and grp is not None and grp(lab) != last:
            last = grp(lab)
            out.write(f"    {last:<12} |{'':{width}}\n", "text")
        vals = [_fv(D[i, j], fmt) for j in range(k)]
        if lower:
            vals = [v if j <= i else "" for j, v in enumerate(vals)]
        line = f"{vals[0]:>10}" + "".join(f"{v:>11}" for v in vals[1:])
        out.write(f"{lab[-12:]:>16} |", "text")
        out.write(f"{line:<{width}}\n", "result")
    out.write("    " + "-" * (14 + width) + "\n", "text")


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
        raw_c, raw_m = rc.copy(), rm.copy()
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
                vm = Sig[i, i] / N - (J_mu[i] @ Vth @ J_mu[i] if kind == "standardized" else 0)
                zm[i] = rm[i] / math.sqrt(vm) if vm > 1e-12 * Sig[i, i] / N else SYS
                for j in range(i + 1):
                    v = (Sig[i, i] * Sig[j, j] + Sig[i, j] ** 2) / N
                    if kind == "standardized":
                        v = v - float(J_sig[i, j] @ Vth @ J_sig[i, j])
                    zc[i, j] = zc[j, i] = rc[i, j] / math.sqrt(v) if v > 0 else SYS
            rc, rm = zc, zm
        rc = np.where(np.abs(rc) < zt, 0.0, rc)
        word = {"raw": "Raw", "normalized": "Normalized", "standardized": "Standardized"}[kind]
        # VERIFY layout of normalized; the raw tables are always shown (case 1909)
        out.write("\nResiduals of observed variables\n\n", "text")
        mrows = [raw_m] + ([rm] if kind != "raw" else [])
        mm = Matrix(np.vstack(mrows), ["raw"] + ([kind] if kind != "raw" else []), obs)
        out.write("  Mean residuals\n\n", "text")
        _fw_table(out, "", mm, None, obs, None, fmt, lower=False)
        out.write("\n  Covariance residuals\n\n", "text")
        _fw_table(out, "", Matrix(raw_c, obs, obs), obs, obs, None, fmt, lower=True)
        if kind != "raw":
            out.write(f"\n  {word} covariance residuals\n\n", "text")
            _fw_table(out, "", Matrix(rc, obs, obs), obs, obs, None, fmt, lower=True)
        mats[_r_name("res_mean", u, est)] = Matrix(raw_m.reshape(1, -1), ["r1"], obs)
        mats[_r_name("res_cov", u, est)] = Matrix(raw_c, obs, obs)
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
        # VERIFY: ties in modulus (+-r) go negative first (compat 1909); the rounding keeps the
        # order from depending on the platform's last bits
        order = np.lexsort((ev.imag, ev.real, -np.round(mod, 10)))
        ev, mod = ev[order], mod[order]
        from .varpost import _cg, _cg7
        out.write("\nStability analysis of simultaneous equation systems\n\n", "text")
        out.write("   Eigenvalue stability condition\n", "text")
        out.write("  +" + "-" * 40 + "+\n", "text")
        out.write("  |        Eigenvalue        |   Modulus   |\n", "text")
        out.write("  |" + "-" * 26 + "+" + "-" * 13 + "|\n", "text")
        for z, md in zip(ev, mod):
            if abs(z.imag) < 1e-12:
                lab = f" {_cg7(z.real):>10}" + " " * 15
            else:
                lab = f" {_cg7(z.real):>10} {'+' if z.imag > 0 else '-'} {_cg7(abs(z.imag)):>9}i  "
            out.write(f"  |{lab}| {_cg(round(md, 6), 9):>9}   |\n", "result")
        out.write("  +" + "-" * 40 + "+\n", "text")
        idx = float(mod.max()) if len(mod) else 0.0
        out.write(f"   stability index = {_cg(round(idx, 6), 9):>9}\n", "text")
        if idx < 1:
            out.write("   All the eigenvalues lie inside the unit circle.\n   SEM satisfies stability condition.\n", "text")
        else:
            # VERIFY
            out.write("   At least one eigenvalue is outside the unit circle.\n"
                      "   SEM does not satisfy stability condition.\n", "text")
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
    # the moments of the observed exogenous variables are conditioned on (fixed) under xconditional
    hid = [i for i, it in enumerate(model.free) if it in model.hidden]
    V = V.copy()
    V[hid, :] = 0.0
    V[:, hid] = 0.0
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
        T_ = base[1]
        for title, w in tables:
            M_ = base[w]
            out.write(f"\n\n{title}\n", "text")
            W = 12
            out.write("-" * 78 + "\n", "text")
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
                    # every pair with an effect is listed in all three tables (case 1909)
                    srcs = [a for a in model.names if a != y and abs(T_[iy, model.idx[a]]) > 1e-12]
                    if srcs:
                        rows.append((y, srcs))
                if not rows:
                    continue
                out.write(f"{sec:<{W}} |\n", "text")
                for k_, (y, srcs) in enumerate(rows):
                    if k_:
                        out.write("  " + "-" * (W - 1) + "+" + "-" * 64 + "\n", "text")
                    out.write(f"  {y + ' <-':<{W - 2}} |\n", "text")
                    for a in srcs:
                        iy, ia = model.idx[y], model.idx[a]
                        b = M_[iy, ia]
                        se = math.sqrt(max(float(J[w][iy, ia] @ V @ J[w][iy, ia]), 0))
                        if abs(b) <= 1e-12 and se == 0:
                            out.write(f"{a[-W:]:>{W}} |", "text")
                            out.write(f"{'0':>11}  (no path)\n", "result")
                        elif se > 0:
                            z = b / se
                            pv = 2 * float(st.norm.sf(abs(z)))
                            out.write(f"{a[-W:]:>{W}} |", "text")
                            out.write(f"  {_g(b):>9}  {_g(se):>9}  {z:>7.2f}  {pv:>6.3f}    "
                                      f"{_g(b - c * se):>9}   {_g(b + c * se):>9}\n", "result")
                        else:
                            out.write(f"{a[-W:]:>{W}} |", "text")
                            out.write(f"  {_g(b):>9}  {'.':>9}\n", "result")
            out.write("-" * 78 + "\n", "text")


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


def _info(ext, t):
    """Information of the extended problem: expected information for ml (as Stata's score tests
    and modification indices, case 1909), minus the Hessian otherwise."""
    if ext.method != "ml":
        return -ext.hess(t)
    I = np.zeros((ext.nglob, ext.nglob))
    for g, (m, d) in enumerate(zip(ext.models, ext.data)):
        ids = ext.gmap[g]
        Ig = _eim(m, t[ids], d.N)
        np.add.at(I, (ids[:, None], ids[None, :]), Ig)
    return I


def _lm(prob, th_glob, ext, t=None) -> tuple[float, int]:
    t = _theta_in(prob, th_glob, ext) if t is None else t
    sc = ext.grad(t)
    info = _info(ext, t)
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
    rows = []
    for k, (text, ext) in enumerate(cons, 1):
        stat, df = _lm(prob, th, ext)
        rows.append((k, stat, max(df, 1)))
    shown = [r for r in rows if r[1] >= minchi2]
    if not shown:
        from ..core.formats import format_value
        # Stata's own value of invchi2tail(1, .05) (case 1909)
        lim = "3.841458820694123" if not o.get("minchi2") else format_value(minchi2, "%18.0g").strip()
        out.write(f"(no score tests to report; all chi2 values less than {lim})\n", "text")
    else:
        # VERIFY: the layout when some test is shown
        out.write("\nScore tests for linear constraints\n\n", "text")
        for k, (text, _ext) in enumerate(cons, 1):
            out.write(f" ({k:2d}) {text}\n", "text")
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


def _eim(model: Model, th: np.ndarray, N: int) -> np.ndarray:
    """Expected information of the normal likelihood: N (dmu' Si dmu + tr(Si dS Si dS) / 2)."""
    mu, Sig, dm = _dmoments(model, th)
    Si = np.linalg.inv(Sig)
    k = len(dm)
    I = np.zeros((k, k))
    for a in range(k):
        ma, Sa = dm[a]
        for b in range(a + 1):
            mb, Sb = dm[b]
            v = 0.0
            if ma is not None and mb is not None:
                v += float(ma @ Si @ mb)
            if Sa is not None and Sb is not None:
                v += 0.5 * float(np.sum((Si @ Sa) * (Si @ Sb).T))
            I[a, b] = I[b, a] = N * v
    return I


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
            # VERIFY (case 1909): Stata does not consider paths out of observed endogenous variables
            if fr == to or (fr, to) in spec.paths or (fr in spec.observed and fr in model.endog):
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
        info = _eim(aug, th_aug, data.N) if prob.method == "ml" else -sub.hess(th_aug)
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
    order = {n: i for i, n in enumerate(model.names)}
    shown = sorted([r for r in rows if r[2] >= minchi2],
                   key=lambda r: (r[0][0] != "path", order[r[0][2]], order[r[0][1]]))
    out = s.output
    out.write("\nModification indices\n\n", "text")
    if not shown:
        out.write(f"(no modification indices to report; all chi2 values less than {minchi2})\n", "text")
    else:
        line = "-" * 14 + "+" + "-" * 46
        out.write("-" * 61 + "\n", "text")
        out.write(f"{'':14}|{'Standard':>46}\n", "text")
        out.write(f"{'':14}|{'MI':>10}{'df':>7}{'P>MI':>7}{'EPC':>11}{'EPC':>11}\n", "text")
        out.write(line + "\n", "text")

        def vals(mi, epc, sepc):
            return (f"{mi:>10.3f}{1:>7}{float(st.chi2.sf(mi, 1)):>7.2f}{_fv(epc, '%9.0g'):>11}"
                    f"{_fv(sepc, '%9.0g'):>11}\n")
        paths = [x for x in shown if x[0][0] == "path"]
        covs = [x for x in shown if x[0][0] != "path"]
        # VERIFY: the Structural section (paths into latent variables) is not in case 1909
        for sec, members in (("Structural", [x for x in paths if x[0][2] in spec.latent]),
                             ("Measurement", [x for x in paths if x[0][2] in spec.observed])):
            if not members:
                continue
            out.write(f"{sec:<14}|\n", "text")
            cur = None
            for it, cls, mi, epc, sepc in members:
                if it[2] != cur:
                    cur = it[2]
                    out.write(f"  {cur + ' <-':<12}|\n", "text")
                out.write(f"{it[1][-12:]:>13} |", "text")
                out.write(vals(mi, epc, sepc), "result")
            out.write(line + "\n", "text")
        for it, cls, mi, epc, sepc in covs:
            out.write(f"{_label(model, it):<14}|", "text")
            out.write(vals(mi, epc, sepc), "result")
        out.write("-" * 61 + "\n", "text")
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
    from .semfit import PCLASSES, item_class, parse_ginvariant
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
    # loadings fixed for identification (same value in every group) are tested by freeing them
    # in the groups after the first (Stata reports a score test for them, case 1909)
    anchors = [("path", m0.names[j], m0.names[i]) for (i, j), v in m0.fixed_A.items()
               if all(m.fixed_A.get((i, j)) == v for m in prob.models)]
    res = []
    for it in anchors:
        stat, df = _lm_anchor(prob, th, it)
        cls = item_class(m0, it)
        if cls in show:
            res.append((it, cls, (np.nan, 0, np.nan), (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)))
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
        res = _ginv_order(m0, res)
        out.write("-" * 78 + "\n", "text")
        out.write(f"{'':12} |{'':12}Wald Test{'':23}Score Test\n", "text")
        out.write(f"{'':12} |{'chi2':>10}{'df':>11}{'p>chi2':>10}{'chi2':>11}{'df':>12}{'p>chi2':>10}\n", "text")
        out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
        eqs = [x for x in res if x[0][0] in ("path", "mean")]
        rest = [x for x in res if x[0][0] not in ("path", "mean")]
        spec = m0.spec
        for sec, pool in (("Structural", spec.latent), ("Measurement", spec.observed)):
            # VERIFY: the Structural section is not in case 1909
            members = [x for x in eqs if x[0][2 if x[0][0] == "path" else 1] in pool]
            if not members:
                continue
            out.write(f"{sec:<12} |\n", "text")
            cur = None
            for it, cls, wald, score in members:
                y = it[2] if it[0] == "path" else it[1]
                if y != cur:
                    if cur is not None:
                        out.write("  " + "-" * 11 + "+" + "-" * 64 + "\n", "text")
                    cur = y
                    out.write(f"  {y + ' <-':<10} |\n", "text")
                lab = it[1] if it[0] == "path" else "_cons"
                out.write(f"{lab[-12:]:>12} |", "text")
                out.write(_tcell(wald) + _tcell(score, wide=True) + "\n", "result")
        if rest:
            out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
            for it, cls, wald, score in rest:
                out.write(f"{_label(m0, it)[-13:]:>13}|", "text")
                out.write(_tcell(wald) + _tcell(score, wide=True) + "\n", "result")
        out.write("-" * 78 + "\n", "text")
    r = {"N_groups": float(G), "nobs": _nobs(s, est)}
    if res:
        r["test"] = Matrix(np.array([[x[2][0], x[2][1], x[2][2], x[3][0], x[3][1], x[3][2]] for x in res],
                                    dtype=float), [_label(m0, x[0]) for x in res],
                           ["Wald_chi2", "Wald_df", "Wald_p", "score_chi2", "score_df", "score_p"])
        r["test_pclass"] = Matrix(np.array([[PCLASSES.index(x[1]) + 1] for x in res], dtype=float),
                                  [_label(m0, x[0]) for x in res], ["pclass"])
    s.r = r


def _tcell(t, wide=False) -> str:
    """Wald cells are 10, 11 and 10 wide; score cells 11, 12 and 10 (case 1909)."""
    stat, df, p = t
    w1, w2 = (11, 12) if wide else (10, 11)
    if df == 0 or not np.isfinite(stat):
        return f"{'.':>{w1}}{'.':>{w2}}{'.':>10}"
    return f"{stat:>{w1}.3f}{df:>{w2}}{p:>10.4f}"


def _ginv_order(m0, res):
    """Equations in the order of the variables, each with its paths and then its intercept;
    variances and covariances after them."""
    pos = {n: i for i, n in enumerate(m0.names)}

    def key(x):
        kind, a, b = x[0]
        if kind == "path":
            return (0, pos[b], 0, pos[a])
        if kind == "mean":
            return (0, pos[a], 1, 0)
        return (1, 0 if kind == "var" else 1, pos[a], pos[b])
    return sorted(res, key=key)


def _lm_anchor(prob, th_glob, it):
    """Score test of a path fixed in every group: the path becomes free (at its fixed value) in
    the groups after the first, which leaves the model identified."""
    from .semfit import Problem
    kind, a, b = it
    models = []
    for g, m in enumerate(prob.models):
        if g == 0:
            models.append(m)
            continue
        mm = copy.copy(m)
        mm.fixed_A = dict(m.fixed_A)
        key = (m.idx[b], m.idx[a])
        v = mm.fixed_A.pop(key)
        mm.free = list(m.free) + [it]
        mm.nfree = len(mm.free)
        mm._anchor_value = v
        models.append(mm)
    ext = Problem(models, prob.data, prob.method, prob.inv,
                  split={(g, it) for g in range(1, len(models))})
    t = np.zeros(ext.nglob)
    for g, m in enumerate(models):
        local = prob.group_theta(th_glob, g)
        vals = np.concatenate([local, [m._anchor_value]]) if g else local
        t[ext.gmap[g]] = vals
    return _lm(prob, th_glob, ext, t=t)


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
