"""streg with distribution(ggamma) and with frailty()/shared() ([ST] streg, Methods and formulas).

Both models have parameters that the one-ancillary-parameter machinery of streg.py does not cover,
so they are fitted here. The generalized gamma log likelihood is differentiated numerically (as
Stata does for it, VERIFY); the frailty models combine the analytic derivatives of the
distribution (through Dist.loglik with d = 1 and d = 0) with a numerical derivative of the
scalar group likelihood. Layout details carry the VERIFY mark (compat/do/1901_survival.do)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lazy import lazy
from .ml import maximize, print_log
from .ml import options as ml_options
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .stfit import STSample, cluster_note, fmt, place
from .vce import bread

special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# generalized gamma
# ---------------------------------------------------------------------------

def _gg_logs(x, mu, lns, kappa, density: bool):
    """ln f(x) (density) or ln S(x) of the generalized gamma with location mu, ln sigma and kappa,
    for x > 0."""
    sig = np.exp(lns)
    lx = np.log(x)
    if abs(kappa) < 1e-6:
        z = (lx - mu) / sig
        if density:
            return -0.5 * z * z - 0.5 * np.log(2 * np.pi) - lns - lx
        return special.log_ndtr(-z)
    g = kappa ** -2
    z = np.sign(kappa) * (lx - mu) / sig
    u = g * np.exp(abs(kappa) * z)
    if density:
        return g * np.log(g) + z * np.sqrt(g) - u - lns - lx - 0.5 * np.log(g) - special.gammaln(g)
    s = special.gammaincc(g, u) if kappa > 0 else special.gammainc(g, u)
    return np.log(np.maximum(s, 1e-300))


class Ggamma:
    """Generalized gamma, in the accelerated failure-time metric only. th = (ln sigma, kappa)."""

    name, title, cmd, predict_sub, user = "ggamma", "Generalized gamma", "gamma", "gamma_p", "gamma_lf0"
    aux = ""
    metrics = ("time",)
    fits_time_metric = True

    def surv(self, x, xb, th, aft):
        lns, kappa = th
        x = np.asarray(x, dtype=float)
        out = np.ones_like(x)
        pos = x > 0
        out[pos] = np.exp(_gg_logs(x[pos], xb[pos] if np.ndim(xb) else xb, lns, kappa, False))
        return out

    def hazard(self, t, xb, th, aft):
        lns, kappa = th
        t = np.asarray(t, dtype=float)
        return np.exp(_gg_logs(t, xb, lns, kappa, True) - _gg_logs(t, xb, lns, kappa, False))

    def median(self, xb, th, aft):
        raise StataError(198, "option median is not available yet in OpenDTA after distribution(ggamma)")


GGAMMA = Ggamma()


def _gg_loglik(b, X, t, t0, d, lnt_c):
    kx = X.shape[1]
    mu = X @ b[:kx]
    lns, kappa = b[kx], b[kx + 1]
    ll = 0.0
    with np.errstate(all="ignore"):
        fail = d > 0
        ll += float(np.sum(_gg_logs(t[fail], mu[fail], lns, kappa, True)))
        ll += float(np.sum(_gg_logs(t[~fail], mu[~fail], lns, kappa, False)))
        has0 = t0 > 0
        ll -= float(np.sum(_gg_logs(t0[has0], mu[has0], lns, kappa, False)))
    return ll + lnt_c


def _numeric(f, b, todo, h0=2e-3):
    """Value, gradient and Hessian of f by central differences (5-point stencils; the Hessian
    differentiates the gradient)."""
    f0 = f(b)
    if todo == 0 or not np.isfinite(f0):
        return f0, None, None
    k = len(b)
    h = h0 * np.maximum(1.0, np.abs(b))

    def grad(x):
        g = np.zeros(k)
        for i in range(k):
            e = np.zeros(k)
            e[i] = h[i]
            g[i] = (-f(x + 2 * e) + 8 * f(x + e) - 8 * f(x - e) + f(x - 2 * e)) / (12 * h[i])
        return g

    g0 = grad(b)
    H = np.zeros((k, k))
    for j in range(k):
        e = np.zeros(k)
        e[j] = h[j]
        H[:, j] = (-grad(b + 2 * e) + 8 * grad(b + e) - 8 * grad(b - e) + grad(b - 2 * e)) / (12 * h[j])
    return f0, g0, (H + H.T) / 2


# ---------------------------------------------------------------------------
# frailty
# ---------------------------------------------------------------------------

def _group_loglik(kind: str, lam, D, vt):
    """Log likelihood of a group with D failures and cumulative hazard lam, frailty with variance
    theta = exp(vt) (the part that does not depend on the log hazards)."""
    th = np.exp(vt)
    if kind == "gamma":
        # D*vt + gammaln(1/th + D) - gammaln(1/th) = sum_{j<D} log(1 + j*th), stable as th -> 0
        dmax = int(np.max(D)) if np.size(D) else 0
        tail = sum(np.where(D > j, np.log1p(j * th), 0.0) for j in range(1, dmax))
        return -(1 / th + D) * np.log1p(th * lam) + tail
    c = (1 + 2 * th * lam) ** -0.5
    x = 1 / (th * c)
    logk = np.log(special.kve(D - 0.5, x)) - x
    return 1 / th + np.log(2) - 0.5 * np.log(2 * np.pi * th) + (D - 0.5) * np.log(c) + logk


def _d1(f, x, h):
    return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h)


def _d2(f, x, h):
    return (-f(x + 2 * h) + 16 * f(x + h) - 30 * f(x) + 16 * f(x - h) - f(x - 2 * h)) / (12 * h * h)


def _group_partials(kind, lam, D, vt, h=2e-3):
    """G and its partial derivatives in (lam, vt), numerically."""
    hl = h / max(1.0, float(np.exp(vt)))
    G = lambda x, v: _group_loglik(kind, x, D, v)           # noqa: E731
    GL = lambda v: _d1(lambda x: G(x, v), lam, hl)          # noqa: E731
    return (G(lam, vt), GL(vt), _d1(lambda v: G(lam, v), vt, h), _d2(lambda x: G(x, vt), lam, hl),
            _d1(GL, vt, h), _d2(lambda v: G(lam, v), vt, h))


def _frailty_fun(dist, X, t, t0, d, gid, kind, lnt_c, fit_aft, k_aux):
    """maximize() function of (beta, [ancillary], ln theta)."""
    kx = X.shape[1]
    lf0 = dist.loglik(t, t0, np.zeros_like(d), fit_aft)
    lf1 = dist.loglik(t, t0, np.ones_like(d), fit_aft)
    ng = int(gid.max()) + 1
    fail = d > 0

    def fun(b, todo):
        th = b[kx] if k_aux else 0.0
        vt = float(b[-1])
        with np.errstate(all="ignore"):
            l0 = lf0(X @ b[:kx], th)
            l1 = lf1(X @ b[:kx], th)
        lam = np.bincount(gid, weights=-l0[0], minlength=ng)
        D = np.bincount(gid, weights=d, minlength=ng)
        with np.errstate(all="ignore"):
            logh = np.where(fail, l1[0] - l0[0], 0.0)
        with np.errstate(all="ignore"):
            G, GL, Gv, GLL, GLv, Gvv = _group_partials(kind, lam, D, vt)
        ll = float(np.sum(logh)) + float(np.sum(G)) + lnt_c
        if not np.isfinite(ll):
            return -np.inf, None, None
        if todo == 0:
            return ll, None, None
        # per observation: lam_x = d lam / d x, h_x = d log h / d x (zero for censored rows)
        lu, lt_, luu, lut, ltt = -l0[1], -l0[2], -l0[3], -l0[4], -l0[5]
        hu, ht = np.where(fail, l1[1] - l0[1], 0.0), np.where(fail, l1[2] - l0[2], 0.0)
        huu, hut, htt = (np.where(fail, l1[3] - l0[3], 0.0), np.where(fail, l1[4] - l0[4], 0.0),
                         np.where(fail, l1[5] - l0[5], 0.0))
        gl = GL[gid]
        p = kx + k_aux + 1
        g = np.zeros(p)
        H = np.zeros((p, p))
        ru = hu + gl * lu                                    # d ll / d u_j
        g[:kx] = X.T @ ru
        Au = np.zeros((ng, kx))
        np.add.at(Au, gid, lu[:, None] * X)                  # sum_j lam_u x_j by group
        H[:kx, :kx] = (X * (huu + gl * luu)[:, None]).T @ X + (Au * GLL[:, None]).T @ Au
        g[-1] = float(np.sum(Gv))
        H[:kx, -1] = H[-1, :kx] = (Au * GLv[:, None]).sum(axis=0)
        H[-1, -1] = float(np.sum(Gvv))
        if k_aux:
            g[kx] = float(np.sum(ht + gl * lt_))
            At = np.bincount(gid, weights=lt_, minlength=ng)
            H[:kx, kx] = H[kx, :kx] = X.T @ (hut + gl * lut) + (Au * (GLL * At)[:, None]).sum(axis=0)
            H[kx, kx] = float(np.sum(htt + gl * ltt)) + float(np.sum(GLL * At * At))
            H[kx, -1] = H[-1, kx] = float(np.sum(GLv * At))
        return ll, g, H
    return fun


def _fit_name(title: str) -> str:
    return title if title in ("Weibull", "Gompertz") else title.lower()


def _fl_title(kind: str, shared: bool) -> str:
    return ("Gamma" if kind == "gamma" else "Inverse-Gaussian") + (" shared frailty" if shared else " frailty")


# ---------------------------------------------------------------------------
# the command
# ---------------------------------------------------------------------------

def fit_extended(s: "Session", stx: STSample, dist, o, args: str, aft: bool) -> None:
    from .streg import DISTS, _model, _pos_log, _to_time_metric
    ds = s.data
    t, t0, d = stx.t, stx.t0, stx.d
    if stx.wtype:
        raise StataError(198, "weights are not available yet in OpenDTA with ggamma or frailty()")
    if stx.robust:
        raise StataError(198, "vce(robust) and vce(cluster) are not available yet in OpenDTA with "
                              "ggamma or frailty()")   # VERIFY: Stata allows them for unshared frailty
    gg = dist.name == "ggamma"
    kind = str(o["frailty"]).strip().lower() if isinstance(o.get("frailty"), str) else ""
    shared = str(o["shared"]).strip() if isinstance(o.get("shared"), str) else ""
    if gg and (kind or shared):
        raise StataError(198, "frailty() and shared() are not allowed with distribution(ggamma)")   # VERIFY
    if shared and not kind:
        kind = "gamma"
    if kind and not gg:
        for full, k in (("gamma", 1), ("invgaussian", 1)):
            if len(kind) >= k and full.startswith(kind):
                kind = full
                break
        else:
            raise StataError(198, f"frailty({kind}) invalid")   # VERIFY
    if shared and np.any(t0 > 0):
        raise StataError(198, "shared() is not allowed with delayed entry or gaps")   # VERIFY
    X = np.column_stack([stx.X, np.ones(stx.n)])
    kx = X.shape[1]
    lnt_c = float(np.sum(d * _pos_log(t)))
    kw = ml_options(o)
    out = s.output
    lab = "log likelihood"
    nslope = kx - 1
    tail: list[CoefRow] = []
    if gg:
        DISTS["ggamma"] = GGAMMA
        mu0 = float(np.mean(_pos_log(t)))
        sd0 = float(np.std(_pos_log(t))) or 1.0
        x1 = np.ones((stx.n, 1))

        def fit(Xm, start):
            return maximize(lambda b, todo: _numeric(lambda z: _gg_loglik(z, Xm, t, t0, d, lnt_c), b, todo),
                            np.asarray(start, dtype=float), **kw)
        res0 = fit(x1, [mu0, float(np.log(sd0)), 0.1])
        b0 = np.zeros(kx + 2)
        b0[kx - 1:] = res0.b
        res = fit(X, b0)
        k_aux = 2
        if not o.get("nolog"):
            out.write("Fitting constant-only model:\n\n", "text")
            print_log(out, res0, lab)
            out.write("\nFitting full model:\n\n", "text")
            print_log(out, res, lab)
        V = bread(res.H)
        b = res.b
        ll_nf = None
    else:
        k_aux = 1 if dist.aux else 0
        fit_aft = aft and dist.fits_time_metric
        rate = float(np.sum(d) / np.sum(t - t0))
        w = np.ones(stx.n)
        lf = dist.loglik(t, t0, d, fit_aft)
        p0 = np.array(dist.start(rate, fit_aft))
        r_nf0 = maximize(_model(lf, np.ones((stx.n, 1)), w, lnt_c, k_aux), p0, **kw)
        b_nf = np.zeros(kx + k_aux)
        b_nf[kx - 1:] = r_nf0.b
        r_nf = maximize(_model(lf, X, w, lnt_c, k_aux), b_nf, **kw)
        gid = (np.unique(_shared_codes(ds, shared, stx.mask), return_inverse=True)[1]
               if shared else np.arange(stx.n))
        gmin = gavg = gmax = ng = None
        if shared:
            cnt = np.bincount(gid)
            ng, gmin, gavg, gmax = len(cnt), int(cnt.min()), float(cnt.mean()), int(cnt.max())
        f0 = _frailty_fun(dist, np.ones((stx.n, 1)), t, t0, d, gid, kind, lnt_c, fit_aft, k_aux)
        f1 = _frailty_fun(dist, X, t, t0, d, gid, kind, lnt_c, fit_aft, k_aux)
        res0 = maximize(f0, np.concatenate([r_nf0.b, [0.0]]), **kw)
        b0 = np.concatenate([r_nf.b, [0.0]])
        if not o.get("nolog"):
            out.write(f"Fitting {_fit_name(dist.title)} model:\n\nFitting constant-only model:\n\n", "text")
            print_log(out, res0, lab)
            out.write("\nFitting full model:\n\n", "text")
        res = maximize(f1, b0, **kw)
        if not o.get("nolog"):
            print_log(out, res, lab)
        V = bread(res.H)
        b = res.b
        if aft and not fit_aft:
            b, V = _to_time_metric(b, V, kx)
        ll_nf = r_nf.ll
    chi2 = 2 * (res.ll - res0.ll)
    pval = float(chi2_tail(chi2, nslope)) if nslope else SYS
    smp = stx.smp
    names, bf, Vf, full_pos, bmap, semap = place(smp, b, V, constant=True, aux=k_aux + (0 if gg else 1))
    est = Estimates("streg", "_t", names, bf, Vf, stx.n, "z", None)
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=True)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, True, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust, est.robust = "", "", None, False
    z = float(st.norm.ppf(1 - (1 - est.level / 100) / 2))
    se = np.sqrt(np.maximum(np.diag(Vf), 0))
    kk = len(bf)
    if gg:
        est.eqnames = ["_t"] * kx + ["ln_sig", "kappa"]
        tail += [CoefRow("sep"),
                 CoefRow("coef", "/ln_sig", bf[kx], se[kx], name="ln_sig:_cons", noeform=True),
                 CoefRow("coef", "/kappa", bf[kx + 1], se[kx + 1], name="kappa:_cons", noeform=True),
                 CoefRow("sep"),
                 _exp_row("sigma", bf[kx], se[kx], z)]
        est.aux = None
        est.th = (bf[kx], bf[kx + 1])
        est.frailty = ""
    else:
        names_aux = ([dist.aux] if dist.aux else []) + ["ln_the"]
        est.eqnames = ["_t"] * kx + names_aux
        if dist.aux:
            tail += [CoefRow("sep"), CoefRow("coef", f"/{dist.aux}", bf[kx], se[kx], name=f"{dist.aux}:_cons",
                                              noeform=True)]
        tail += [CoefRow("coef", "/ln_the", bf[-1], se[-1], name="ln_the:_cons", noeform=True)] \
            if dist.aux else [CoefRow("sep"), CoefRow("coef", "/ln_the", bf[-1], se[-1], name="ln_the:_cons",
                                                      noeform=True)]
        tail += [CoefRow("sep")]
        if dist.aux:
            tail += dist.aux_rows(bf[kx], se[kx], bf[kx] - z * se[kx], bf[kx] + z * se[kx])
        tail += [_exp_row("theta", bf[-1], se[-1], z)]
        est.aux = (float(bf[kx]), float(se[kx])) if dist.aux else None
        est.th = float(bf[kx]) if dist.aux else 0.0
        est.frailty = kind
        est.shared = shared
        est.theta = float(np.exp(bf[-1]))
        lr = 2 * (res.ll - ll_nf)
        est.lr_theta = (max(lr, 0.0), 0.5 * float(chi2_tail(lr, 1)) if lr > 0 else 1.0)
        est.groups = (ng, gmin, gavg, gmax) if shared else None
        est.k_aux, est.ll_nf = k_aux + 1, ll_nf
    est.tail_rows = tail
    x = stx.scalars()
    x.update({"ll": res.ll, "ll_0": res0.ll, "chi2": chi2, "p": pval, "df_m": nslope, "chi2type": "LR"})
    est.extra = x
    est.dist, est.aft = dist.name, aft
    est.distobj = dist
    est.eform_kind = "tr" if (aft and o.get("tr")) else ("hr" if (not aft and not o.get("nohr")) else "")
    est.display = _display
    est.predict = _predict
    est.predict_opts = {"xb": 2}
    est.cmd = "streg"
    _post(s, est, args, dist, res, res0, kk, stx.mask)
    _display(s, est, header=not o.get("noheader"))


def _shared_codes(ds, name: str, mask: np.ndarray) -> np.ndarray:
    from ..core.varlist import resolve_name
    v = ds.get(resolve_name(ds, name))
    x = v.raw[mask]
    return np.array([str(z) for z in x]) if v.is_string else v.data[mask]


def _exp_row(label: str, b: float, se: float, z: float) -> CoefRow:
    with np.errstate(over="ignore"):
        v = float(np.exp(b))
        return CoefRow("aux", label, v, v * se, eq=f"{float(np.exp(b - z * se))!r} {float(np.exp(b + z * se))!r}")


def _post(s: "Session", est: Estimates, args: str, dist, res, res0, k: int, mask) -> None:
    from .regress import cmdline
    x = est.extra
    fr = getattr(est, "frailty", "")
    scalars = [("rank", float(k)), ("N", x["N"]), ("ic", float(res.iterations)), ("k", float(k)),
               ("k_eq", float(len(set(est.eqnames)))),
               ("k_dv", 1.0), ("converged", 1.0 if res.converged else 0.0), ("rc", 0.0),
               ("ll", res.ll), ("k_eq_model", 1.0), ("ll_0", res0.ll), ("rank0", float(len(res0.b))),
               ("df_m", float(x["df_m"])), ("chi2", x["chi2"]), ("p", x["p"])]
    if fr:
        scalars.append(("theta", est.theta))
        if est.groups is not None:
            ng, gmin, gavg, gmax = est.groups
            scalars += [("N_g", float(ng)), ("g_min", float(gmin)), ("g_max", float(gmax))]
        lr, pc = est.lr_theta
        scalars += [("risk", x["risk"]), ("N_fail", x["N_fail"]), ("N_sub", x["N_sub"]),
                    ("k_aux", float(est.k_aux)), ("p_c", pc), ("chi2_c", lr), ("ll_c", est.ll_nf)]
        if est.groups is not None:
            scalars.append(("g_avg", gavg))
        het = dist.cmd + "het"   # VERIFY: names of the evaluators for the other distributions
        user = f"{het}_{'i' if fr == 'invgaussian' else 'g'}lf" + ("_sh" if est.shared else "")
        macros = [("cmdline", cmdline("streg", args)), ("predict", "streg_p"), ("predict_sub", het + "_ip"),
                  ("cmd2", "streg"), ("title", dist.title + " regression"),
                  ("title2", "accelerated failure-time form" if est.aft else "log relative-hazard form"),
                  ("cmd", het), ("footnote", "streghet_footnote"), ("stcurve", "stcurve"),
                  ("shared", est.shared), ("dead", "_d"), ("t0", "_t0"),
                  ("fr_title", _fl_title(fr, bool(est.shared))),
                  ("frm2", "time" if est.aft else "hazard"), ("chi2type", "LR"), ("opt", "moptimize"),
                  ("vce", "oim"), ("user", user), ("ml_method", "d2"), ("technique", "nr"),
                  ("which", "max"), ("depvar", "_t"), ("properties", "b V")]
    else:
        scalars += [("k_aux", 2.0), ("N_sub", x["N_sub"]), ("sigma", float(np.exp(est.th[0]))),
                    ("kappa", float(est.th[1])), ("risk", x["risk"]), ("N_fail", x["N_fail"])]
        macros = [("cmdline", cmdline("streg", args)), ("stcurve", "stcurve"), ("predict", "streg_p"),
                  ("predict_sub", dist.predict_sub), ("cmd2", "streg"), ("frm2", "time"),
                  ("title", dist.title + " regression"), ("cmd", dist.cmd), ("dead", "_d"), ("t0", "_t0"),
                  ("title2", "accelerated failure-time form"), ("chi2type", "LR"), ("opt", "moptimize"),
                  ("vce", "oim"), ("user", "gamma_d2"), ("ml_method", "e2"), ("technique", "nr"),
                  ("which", "max"), ("depvar", "_t"), ("properties", "b V")]   # VERIFY: the e() lists
    post(s, est, scalars, macros, mask)


def _display(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
             level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    form = "accelerated failure-time form" if est.aft else "log relative-hazard form"
    if header:
        shared = bool(est.frailty and est.shared)
        subj = [f"{'No. of subjects':<15} ={fmt(x['N_sub'], '%12.0gc'):>13}",
                f"{'No. of failures':<15} ={fmt(x['N_fail'], '%12.0gc'):>13}",
                f"{'Time at risk':<15} ={fmt(x['risk'], '%12.0g'):>13}"]
        ll = f"{'Log likelihood':<15} ={fmt(x['ll'], '%10.0g'):>13}"
        lr = (f"LR chi2({int(x['df_m'])})", f"{x['chi2']:.2f}")
        pr = ("Prob > chi2", f"{x['p']:.4f}" if x["p"] < SYS else ".")
        if shared:
            ng, gmin, gavg, gmax = est.groups
            out.write(f"\n{est.distobj.title} regression --\n", "text")
            for lft, lab, v in ((form, "Number of obs", x["N"]), (_fl_title(est.frailty, True), "Number of groups", ng)):
                out.write(f"{'':9}{lft:<39}{lab:<18}=", "text")
                out.write(f"{fmt(v, '%12.0gc'):>11}\n", "result")
            out.write(f"Group variable: {est.shared}\n{'':48}Obs per group:\n", "text")
            for line, lab, v in zip(subj, ("min", "avg", "max"), (gmin, gavg, gmax)):
                out.write(f"{line:<48}", "text")
                out.write(f"{lab:>17} =" + f"{fmt(v, '%12.0g'):>11}\n", "result")
            out.write("\n", "text")
            for line, (lab, v) in ((" ", lr), (ll, pr)):
                out.write(f"{line.rstrip() if line.strip() else '':<48}{lab:<18}=", "text")
                out.write(f"{v:>10}\n", "result")
            out.write("\n", "text")
        else:
            out.write(f"\n{est.distobj.title} regression -- {form}\n", "text")
            if est.frailty:
                out.write(f"{'':22}{_fl_title(est.frailty, False)}\n", "text")
            out.write("\n", "text")
            right = [("Number of obs", fmt(x["N"], "%12.0gc")), None, None, lr, pr]
            for lft, rgt in zip(subj + [None, ll], right):
                line = lft if lft is not None else ""
                if rgt is not None:
                    out.write(f"{line:<48}{rgt[0]:<17}=", "text")
                    out.write(f"{rgt[1]:>12}\n", "result")
                else:
                    out.write(line + "\n", "text")
            out.write("\n", "text")
    if not table:
        return
    rows = list(est.rows) + list(est.tail_rows)
    title = {"hr": "Haz. Ratio", "tr": "Tm. Ratio"}.get(est.eform_kind, "Coef.")
    coef_table(s, "_t", rows, stat="z", df=None, level=level or est.level, vcetype="",
               coef_title=title, eform=bool(est.eform_kind))
    if est.frailty and est.lr_theta is not None:
        lr, p = est.lr_theta
        lrs = f"{lr:>8.2f}" if lr == 0 or lr >= 0.005 else f"{lr:>8.2g}"   # VERIFY: tiny values
        out.write(f"Likelihood-ratio test of theta=0: chibar2(01) = {lrs} Prob>=chibar2 = {p:.3f}\n", "text")


def _predict(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    if stat in ("", "xb"):
        return np.where(ok, xb, SYS), "" if stat else "(option xb assumed; fitted values)"
    raise StataError(198, f"option {stat} is not available yet in OpenDTA after distribution(ggamma) or "
                          "frailty()")
