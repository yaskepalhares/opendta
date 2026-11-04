"""mi impute ([MI] mi impute): univariate methods (regress, pmm, logit, ologit, mlogit, poisson,
nbreg, truncreg, intreg), mvn and chained, with the posterior draws of the parameters."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lazy import lazy
from .mi import _METHODS, _ch, _list, _num, _set_list, _touch, _update_miss, nimp, require_mi

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

SYS = M.SYSMISS


def _rng_normal(n: int) -> np.ndarray:
    from ..core import rng
    return special.ndtri(rng.RNG.uniform(n))


def _rng_uniform(n: int) -> np.ndarray:
    from ..core import rng
    return rng.RNG.uniform(n)


def _rng_chi2(df: float) -> float:
    return float(st.chi2.ppf(_rng_uniform(1)[0], df))


def _design(ds, names: list[str], rows: np.ndarray, constant: bool = True) -> np.ndarray:
    from . import fvars
    terms = fvars.expand_fv(ds, " ".join(names)) if names else []
    mask = np.zeros(ds.nobs, dtype=bool)
    mask[rows] = True
    cols = fvars.build_design(ds, terms, mask, constant=constant) if terms else []
    fvars.drop_collinear(cols, constant)
    X = [c.values for c in cols if c.values is not None]
    if constant:
        X.append(np.ones(len(rows)))
    return np.column_stack(X) if X else np.zeros((len(rows), 0))


def _fit_logit_like(kind: str, y, X):
    """ML estimates and variance for logit, poisson (and the linear index for ologit/mlogit)."""
    from .ml import maximize
    k = X.shape[1]
    if kind == "logit":
        def fun(b, todo):
            xb = X @ b
            ll = float(np.sum(y * xb - np.logaddexp(0, xb)))
            if todo == 0:
                return ll, None, None
            p = special.expit(xb)
            return ll, X.T @ (y - p), -(X * (p * (1 - p))[:, None]).T @ X
    else:
        def fun(b, todo):
            xb = X @ b
            mu = np.exp(np.clip(xb, -700, 700))
            ll = float(np.sum(y * xb - mu))
            if todo == 0:
                return ll, None, None
            return ll, X.T @ (y - mu), -(X * mu[:, None]).T @ X
    r = maximize(fun, np.zeros(k))
    return r.b, np.linalg.inv(-r.H)


def _impute_one(method: str, y_obs, X_obs, X_mis, opts, levels=None):
    """Imputed values for the missing rows (one imputation)."""
    n, k = X_obs.shape
    if method in ("regress", "pmm"):
        XtXi = np.linalg.pinv(X_obs.T @ X_obs)
        b = XtXi @ X_obs.T @ y_obs
        e = y_obs - X_obs @ b
        df = n - k
        s2 = float(e @ e) / df
        sig2 = s2 * df / _rng_chi2(df)
        L = np.linalg.cholesky(XtXi + 1e-14 * np.eye(k))
        bstar = b + math.sqrt(sig2) * (L @ _rng_normal(k))
        if method == "regress":
            return X_mis @ bstar + math.sqrt(sig2) * _rng_normal(X_mis.shape[0])
        knn = int(opts.get("knn", 1) or 1)
        yhat_obs = X_obs @ b
        yhat_mis = X_mis @ bstar
        out = np.empty(X_mis.shape[0])
        u = _rng_uniform(X_mis.shape[0])
        for i, v in enumerate(yhat_mis):
            d = np.abs(yhat_obs - v)
            nn = np.argsort(d, kind="stable")[:knn]
            out[i] = y_obs[nn[min(int(u[i] * len(nn)), len(nn) - 1)]]
        return out
    if method in ("logit", "poisson"):
        lev0 = None
        if method == "logit":
            lev = np.unique(y_obs)
            lev0, lev1 = (lev[0], lev[-1]) if len(lev) == 2 else (0.0, 1.0)
            yb = (y_obs == lev1).astype(float)
        else:
            yb = y_obs
        b, V = _fit_logit_like(method, yb, X_obs)
        L = np.linalg.cholesky((V + V.T) / 2 + 1e-14 * np.eye(k))
        bstar = b + L @ _rng_normal(k)
        eta = X_mis @ bstar
        if method == "logit":
            p = special.expit(eta)
            return np.where(_rng_uniform(len(p)) < p, lev1, lev0)
        mu = np.exp(eta)
        return st.poisson.ppf(_rng_uniform(len(mu)), mu)
    if method in ("ologit", "mlogit"):
        levs = np.unique(y_obs)
        J = len(levs)
        idx = np.searchsorted(levs, y_obs)
        from .ml import maximize
        Xn = X_obs[:, :-1]          # no constant: cut points or one constant per outcome
        if method == "ologit":
            p_ = Xn.shape[1]

            def fun(par, todo):
                b, cut = par[:p_], np.cumsum(np.concatenate([[par[p_]], np.exp(par[p_ + 1:])]))
                xb = Xn @ b
                c = np.concatenate([[-np.inf], cut, [np.inf]])
                pr = special.expit(c[idx + 1] - xb) - special.expit(c[idx] - xb)
                return float(np.sum(np.log(np.maximum(pr, 1e-300)))), None, None
        else:
            p_ = X_obs.shape[1]

            def fun(par, todo):
                B = np.column_stack([np.zeros(p_), par.reshape(J - 1, p_).T])
                eta = X_obs @ B
                ll = float(np.sum(eta[np.arange(len(idx)), idx] - special.logsumexp(eta, axis=1)))
                return ll, None, None
        from .mixed import _numderiv

        def full(par, todo):
            if todo == 0:
                return fun(par, 0)
            return _numderiv(lambda q: fun(q, 0)[0], par)
        npar = (Xn.shape[1] + J - 1) if method == "ologit" else (J - 1) * X_obs.shape[1]
        start = np.zeros(npar)
        if method == "ologit":
            start[Xn.shape[1]] = -1.0
        r = maximize(full, start, maxiter=200)
        V = np.linalg.inv(-r.H)
        L = np.linalg.cholesky((V + V.T) / 2 + 1e-12 * np.eye(npar))
        par = r.b + L @ _rng_normal(npar)
        if method == "ologit":
            b, cut = par[:Xn.shape[1]], np.cumsum(np.concatenate([[par[Xn.shape[1]]], np.exp(par[Xn.shape[1] + 1:])]))
            xb = X_mis[:, :-1] @ b
            cum = special.expit(cut[None, :] - xb[:, None])
            u = _rng_uniform(len(xb))
            k_ = np.sum(u[:, None] > cum, axis=1)
            return levs[k_]
        B = np.column_stack([np.zeros(p_), par.reshape(J - 1, p_).T])
        pr = special.softmax(X_mis @ B, axis=1)
        u = _rng_uniform(len(pr))
        k_ = np.minimum(np.sum(u[:, None] > np.cumsum(pr, axis=1), axis=1), J - 1)
        return levs[k_]
    if method == "nbreg":
        return _impute_nbreg(y_obs, X_obs, X_mis)
    if method in ("truncreg", "intreg"):
        return _impute_censored(method, y_obs, X_obs, X_mis, opts)
    raise StataError(198, f"mi impute {method}: not available in OpenDTA")   # OpenDTA


def _draw_params(ll, start: np.ndarray) -> np.ndarray:
    """ML estimates of ll(par) and a draw from the normal approximation of the posterior."""
    from .ml import maximize
    from .mixed import _numderiv

    def full(par, todo):
        return (ll(par), None, None) if todo == 0 else _numderiv(ll, par)
    r = maximize(full, start, maxiter=200)
    V = np.linalg.pinv(-r.H)    # a flat direction (alpha -> 0 in nbreg) must not stop the draw
    L = np.linalg.cholesky((V + V.T) / 2 + 1e-12 * np.eye(len(start)))
    return r.b + L @ _rng_normal(len(start))


def _impute_nbreg(y, X, Xm):
    """Negative binomial (NB2) regression: parameters drawn from their asymptotic posterior,
    then the counts from the gamma-Poisson mixture (VERIFY: the draws are not Stata's)."""
    k = X.shape[1]
    lg = special.gammaln

    def ll(par):
        mu = np.exp(np.clip(X @ par[:k], -30, 30))
        a = math.exp(min(max(par[k], -20), 20))
        r = 1.0 / a
        return float(np.sum(lg(y + r) - lg(r) - lg(y + 1) + y * np.log(a * mu / (1 + a * mu))
                            - r * np.log1p(a * mu)))
    b0, _ = _fit_logit_like("poisson", y, X)
    par = _draw_params(ll, np.concatenate([b0, [0.0]]))
    mu = np.exp(np.clip(Xm @ par[:k], -30, 30))
    a = math.exp(min(max(par[k], -20), 20))
    r = 1.0 / a
    return st.nbinom.ppf(np.minimum(_rng_uniform(len(mu)), 1 - 1e-12), r, 1.0 / (1.0 + a * mu))


def _impute_censored(method, y, X, Xm, opts):
    """truncreg: the observed values come from a normal truncated to [ll, ul]; intreg: the
    observed values are exact and the missing ones are interval-censored in [ll_i, ul_i]. In
    both the parameters are drawn from the normal approximation of the posterior and the
    missing values from the truncated normal (VERIFY against Stata's draws)."""
    k = X.shape[1]
    ndtr, ndtri = special.ndtr, special.ndtri
    lo_o, hi_o = opts["ll_obs"], opts["ul_obs"]
    lo_m, hi_m = opts["ll_mis"], opts["ul_mis"]
    b_ols = np.linalg.pinv(X.T @ X) @ X.T @ y
    s0 = math.sqrt(max(float(np.mean((y - X @ b_ols) ** 2)), 1e-12))

    def mass(z_lo, z_hi):
        return np.maximum(ndtr(z_hi) - ndtr(z_lo), 1e-300)

    def ll(par):
        s = math.exp(min(max(par[k], -20), 20))
        z = (y - X @ par[:k]) / s
        v = float(np.sum(-0.5 * z * z - 0.5 * math.log(2 * math.pi) - par[k]))
        if method == "truncreg":
            xb = X @ par[:k]
            v -= float(np.sum(np.log(mass((lo_o - xb) / s, (hi_o - xb) / s))))
        else:
            xm = Xm @ par[:k]
            v += float(np.sum(np.log(mass((lo_m - xm) / s, (hi_m - xm) / s))))
        return v
    par = _draw_params(ll, np.concatenate([b_ols, [math.log(s0)]]))
    s = math.exp(min(max(par[k], -20), 20))
    xm = Xm @ par[:k]
    zl, zh = (lo_m - xm) / s, (hi_m - xm) / s
    pl, ph = ndtr(zl), ndtr(zh)
    u = _rng_uniform(len(xm))
    p = np.clip(pl + u * (ph - pl), 1e-15, 1 - 1e-15)
    return xm + s * ndtri(p)


def _parse_impute(text: str):
    c = find_top(text, ",")
    body, opts = (text, "") if c == -1 else (text[:c], text[c + 1:])
    return body.strip(), opts


def mi_impute(s: "Session", args: str) -> None:
    style = require_mi(s)
    ds = s.data
    words = args.split(None, 1)
    if not words:
        raise StataError(198, "method required")   # VERIFY
    meth = next((m for m, k in _METHODS.items() if len(words[0]) >= k and m.startswith(words[0])), None)
    if meth is None:
        raise StataError(198, f"mi impute {words[0]}: invalid method")   # VERIFY
    rest = words[1] if len(words) > 1 else ""
    body, opts = _parse_impute(rest)
    o = match_options(opts, {"add": 3, "replace": 3, "rseed": 2, "noisily": 1, "force": 5, "knn": 3,
                             "burnin": 4, "dots": 3, "nolegend": 5, "augment": 3, "orderasis": 6,
                             "noisure": 5, "ll": 2, "ul": 2, "burnbetween": 9}) if opts.strip() else {}
    if o.get("rseed"):
        from ..core import rng
        rng.reset(int(float(o["rseed"])))
    M_old = nimp(ds)
    add = int(o["add"]) if o.get("add") else 0
    if not add and not o.get("replace"):
        raise StataError(198, "either add() or replace must be specified")
    if o.get("replace") and not add:
        ms = list(range(1, M_old + 1))
        M_new = M_old
    else:
        ms = list(range(M_old + 1, M_old + add + 1))
        M_new = M_old + add
    ivars_all = _list(ds, "_mi_ivars")
    # specification: [(method, target, opts)] and the covariates
    if meth in ("chained", "monotone"):
        specs, covars = _parse_chained(s, body)
        if meth == "monotone":
            _check_monotone(ds, specs)
    elif meth == "mvn":
        lhs, _eq, rhs = body.partition("=")
        tg = s.expand_varlist(lhs.strip())
        if len(tg) < 2:
            raise StataError(102, "too few variables specified")
        specs = [("mvn", v, {}) for v in tg]
        covars = s.expand_varlist(rhs) if rhs.strip() else []
    else:
        if "=" in body:
            lhs, rhs = body.split("=", 1)
        else:
            lhs, rhs = body, ""
        tgt = s.expand_varlist(lhs.strip())
        if len(tgt) != 1:
            raise StataError(103, "too many variables specified")
        specs = [(meth, tgt[0], {"knn": o.get("knn"), "ll": o.get("ll"), "ul": o.get("ul")})]
        if meth in ("truncreg", "intreg") and (o.get("ll") is None and o.get("ul") is None):
            raise StataError(198, "option ll() or ul() required")   # VERIFY
        covars = s.expand_varlist(rhs) if rhs.strip() else []
    for _m, v, _o in specs:
        if v not in ivars_all:
            raise StataError(498, f"variable {v} not registered as imputed")   # VERIFY
    mm = _num(ds, "_mi_m")
    # make room for the new imputations (mlong: incomplete rows; flong: all rows)
    if M_new > M_old:
        _add_rows(ds, style, M_old, M_new)
        mm = _num(ds, "_mi_m")
    ids = _num(ds, "_mi_id")
    m0 = np.nonzero(mm == 0)[0]
    out = s.output
    counts = {}
    mvn_draws = None
    for km, m in enumerate(ms):
        rows_m = _rows_of(ds, style, m)         # m=0 order, with the rows of m where they exist
        work = {v: _num(ds, v)[rows_m] for _m, v, _o in specs}
        missing = {v: _num(ds, v)[m0] >= SYS for _m, v, _o in specs}
        cov = _design_rows(ds, covars, rows_m)
        if meth == "chained":
            burn = int(o["burnin"]) if o.get("burnin") else 10
            _chained(specs, work, missing, cov, burn)
        elif meth == "monotone":
            done: list[str] = []
            for mth, v, so in specs:
                _cond(mth, v, so, work, missing, cov, done, o, ds, rows_m)
                done.append(v)
        elif meth == "mvn":
            if mvn_draws is None:
                mvn_draws = _mvn_draws(ds, specs, covars, rows_m, len(ms), o)
                ll_em, it_em, n_omit = o["_em"]
                out.write("\nPerforming EM optimization:\n", "text")
                if n_omit:
                    out.write(f"note: {n_omit} observations omitted from EM estimation because of all "
                              "imputation variables missing\n", "text")
                out.write(f"  observed log likelihood = {_g10(ll_em)} at iteration {it_em}\n", "text")
                out.write("\nPerforming MCMC data augmentation ... \n", "text")
            for v, arr in zip([sp[1] for sp in specs], mvn_draws[km].T):
                work[v] = np.where(missing[v], np.where(np.isfinite(arr), arr, SYS), work[v])
        else:
            mth, v, so = specs[0]
            mis = missing[v]
            ok = ~mis & np.all(np.isfinite(cov), axis=1)
            if mis.any():
                Xo = np.column_stack([cov[ok], np.ones(int(ok.sum()))])
                Xm = np.column_stack([cov[mis], np.ones(int(mis.sum()))])
                so = _limits(ds, so, rows_m, ok, mis)
                if mth == "intreg":
                    # Stata registers the variables of ll() and ul() as passive (case 1909; VERIFY why)
                    lim = [str(so_[k]).strip() for so_ in (specs[0][2],) for k in ("ll", "ul")
                           if so_.get(k) is not None and ds.has(str(so_[k]).strip())]
                    regd = set(_list(ds, "_mi_ivars") + _list(ds, "_mi_pvars") + _list(ds, "_mi_rvars"))
                    new = [v_ for v_ in lim if v_ not in regd]
                    if new:
                        _set_list(ds, "_mi_pvars", _list(ds, "_mi_pvars") + new)
                vals = _impute_one(mth, work[v][ok], Xo, Xm, so)
                w = work[v].copy()
                w[mis] = vals
                work[v] = w
        for _m2, v, _o in specs:
            arr = ds.get(v).data.astype(np.float64)
            arr[rows_m] = np.where(missing[v], work[v], arr[rows_m])
            ds.get(v).data = arr
            counts.setdefault(v, []).append(int(np.sum(missing[v] & (work[v] < SYS))))
    _ch(ds)["_mi_M"] = str(M_new)
    _update_miss(ds)
    _touch(ds)
    s.notify_state()
    # report
    n0 = len(m0)
    title = "Multivariate imputation" if meth in ("chained", "monotone", "mvn") else "Univariate imputation"
    mname = {"regress": "Linear regression", "pmm": "Predictive mean matching", "logit": "Logistic regression",
             "ologit": "Ordered logistic regression", "mlogit": "Multinomial logistic regression",
             "poisson": "Poisson regression", "chained": "Chained equations",
             "nbreg": "Negative binomial regression", "truncreg": "Truncated regression",
             "intreg": "Interval regression", "monotone": "Monotone method",
             "mvn": "Multivariate normal regression"}[meth]
    if meth == "monotone":
        out.write("\nConditional models:\n", "text")
        prev: list[str] = []
        for mth, v, _o in specs:
            out.write(f"{v:>18}: " + " ".join([mth, v] + prev + covars) + "\n", "text")
            prev.append(v)
        out.write("\n", "text")
    out.write("\n", "text")
    out.write(f"{title:<44}Imputations = ", "text")
    out.write(f"{M_new:>8,}\n", "result")
    out.write(f"{mname:<50}added = ", "text")
    out.write(f"{(M_new - M_old):>8,}\n", "result")
    out.write(f"{'Imputed: m=' + str(ms[0]) + ' through m=' + str(ms[-1]):<48}updated = ", "text")
    out.write(f"{(len(ms) if not add else 0):>8,}\n", "result")
    if meth == "chained":
        out.write(f"\nInitialization: monotone{'':20}Iterations = ", "text")
        out.write(f"{int(o['burnin']) * len(ms) if o.get('burnin') else 10 * len(ms):>8,}\n", "result")
        out.write(f"{'burn-in = ':>62}", "text")
        out.write(f"{int(o['burnin']) if o.get('burnin') else 10:>8}\n", "result")
        out.write("\n", "text")
        for mth, v, _o in specs:
            out.write(f"{v:>17}: {mth} {v} {' '.join([w for _m, w, _o in specs if w != v] + covars)}\n", "text")
    if meth == "mvn":
        burn = int(o["burnin"]) if o.get("burnin") else 100
        btw = int(o["burnbetween"]) if o.get("burnbetween") else 100
        out.write(f"\n{'Prior: uniform':<45}Iterations = {burn + (len(ms) - 1) * btw:>8}\n", "text")
        out.write(f"{'burn-in = ':>58}{burn:>8}\n{'between = ':>58}{btw:>8}\n", "text")
    elif meth == "monotone":
        names = {"regress": "linear regression", "pmm": "predictive mean matching", "logit": "logistic regression",
                 "ologit": "ordered logistic regression", "mlogit": "multinomial logistic regression",
                 "poisson": "Poisson regression", "nbreg": "negative binomial regression",
                 "truncreg": "truncated regression", "intreg": "interval regression"}     # VERIFY all but regress
        out.write("\n", "text")
        for mth, v, _o in specs:
            out.write(f"{v:>18}: {names.get(mth, mth)}\n", "text")
    elif meth == "nbreg":
        out.write("\nDispersion: mean\n", "text")
    elif meth in ("truncreg", "intreg"):
        out.write(_limit_lines(ds, specs[0], m0, o), "text")
    out.write("\n" + "-" * 66 + "\n", "text")
    out.write(f"{'':19}|{'':15}Observations per m{'':13}\n", "text")
    out.write(f"{'':19}|" + "-" * 46 + "\n", "text")
    out.write(f"{'Variable':>18} |   Complete   Incomplete   Imputed |     Total\n", "text")
    out.write("-" * 19 + "+" + "-" * 35 + "+" + "-" * 10 + "\n", "text")
    for _m, v, _o in specs:
        inc = int(np.sum(_num(ds, v)[m0] >= SYS))
        out.write(f"{v[:17]:>18} |", "text")
        out.write(f" {n0 - inc:>10,} {inc:>12,} {min(counts[v]):>9,} |{n0:>10,}\n", "result")
    out.write("-" * 66 + "\n", "text")
    out.write("(complete + incomplete = total; imputed is the minimum across m\n"
              " of the number of filled-in observations.)\n", "text")


def _g10(x: float) -> str:
    from ..core.formats import format_value
    return format_value(float(x), "%10.0g", pad=False).strip()


def _limit_lines(ds, spec, m0, o) -> str:
    """The Limit: block of mi impute truncreg and intreg (case 1909)."""
    meth, v, so = spec
    lo, hi = so.get("ll"), so.get("ul")

    def show(x, inf):
        if x is None:
            return inf
        try:
            return _g10(float(x))
        except ValueError:
            return str(x)

    def vals(x):
        if x is None:
            return None
        try:
            return np.full(len(m0), float(x))
        except ValueError:
            return _num(ds, str(x))[m0]
    lv, hv = vals(lo), vals(hi)
    y = _num(ds, v)[m0]
    if meth == "truncreg":
        left = int(np.sum((y < SYS) & (lv is not None) & (y <= (lv if lv is not None else -np.inf)))) if lv is not None else 0
        right = int(np.sum((y < SYS) & (y >= hv))) if hv is not None else 0
        return (f"\nLimit: lower = {show(lo, '-inf'):>12}            Number truncated = {left + right:>8}\n"
                f"       upper = {show(hi, '+inf'):>12}{'left = ':>31}{left:>8}\n"
                f"{'right = ':>58}{right:>8}\n")
    mis = y >= SYS
    lmis = lv >= SYS if lv is not None else np.ones(len(y), bool)
    hmis = hv >= SYS if hv is not None else np.ones(len(y), bool)
    both = mis & lmis & hmis
    interval = int(np.sum(mis & ~lmis & ~hmis))
    left = int(np.sum(mis & lmis & ~hmis))
    right = int(np.sum(mis & ~lmis & hmis))
    cens = interval + left + right
    return (f"\nLimit: lower = {show(lo, '-inf'):>12}              Number missing = {int(both.sum()):>8}\n"
            f"       upper = {show(hi, '+inf'):>12}{'Number censored = ':>31}{cens:>8}\n"
            f"{'interval = ':>58}{interval:>8}\n{'left = ':>58}{left:>8}\n{'right = ':>58}{right:>8}\n")


def _limits(ds, so: dict, rows_m, ok, mis) -> dict:
    """truncreg and intreg: the limits ll() and ul() (a number or a variable) of the observed
    and of the missing rows."""
    if so.get("ll") is None and so.get("ul") is None:
        return so

    def get(spec, sign):
        if spec is None:
            return np.full(len(rows_m), sign * np.inf)
        try:
            return np.full(len(rows_m), float(spec))
        except ValueError:
            v = _num(ds, str(spec).strip())[rows_m]
            return np.where(v >= SYS, sign * np.inf, v)
    lo, hi = get(so.get("ll"), -1), get(so.get("ul"), 1)
    out = dict(so)
    out.update(ll_obs=lo[ok], ul_obs=hi[ok], ll_mis=lo[mis], ul_mis=hi[mis])
    return out


def _check_monotone(ds, specs) -> None:
    """The missing values must be nested in the order given: v1 observed whenever v2 is."""
    prev = None
    for _m, v, _o in specs:
        cur = _num(ds, v)[_num(ds, "_mi_m") == 0] >= SYS
        if prev is not None and np.any(prev & ~cur):
            raise StataError(459, "nonmonotone missing-value pattern")   # VERIFY
        prev = cur


def _mvn_draws(ds, specs, covars, rows_m, nimp, o):
    """Completed data of mi impute mvn: data augmentation started from the observed data, with
    burnin() iterations (default 100) before the first imputation and burnbetween() (default
    100) between the others (VERIFY: the sampler is not Stata's)."""
    from ..core import mvnmiss
    Y = np.column_stack([_num(ds, v)[rows_m] for _m, v, _o in specs])
    Y = np.where(Y >= SYS, np.nan, Y)
    cov = _design_rows(ds, covars, rows_m)
    if cov.shape[1] == 0 or not np.any(np.all(cov == 1.0, axis=0)):
        cov = np.column_stack([cov, np.ones(len(rows_m))])
    good = np.all(np.isfinite(cov), axis=1)
    # the EM estimates (the starting values of Stata's data augmentation) and their report
    usable = good & ~np.all(np.isnan(Y), axis=1)
    _B, _S, ll_em, it_em = mvnmiss.em_mvreg(Y[usable], cov[usable])
    o["_em"] = (ll_em, it_em, int(np.sum(good & ~usable)))
    burn = int(o["burnin"]) if o.get("burnin") else 100
    btw = int(o["burnbetween"]) if o.get("burnbetween") else 100
    comp = mvnmiss.da_mvn(Y[good], cov[good], _rng_normal, _rng_chi2, nimp, burn, btw)
    out = []
    for Yc in comp:
        full = Y.copy()
        full[good] = Yc
        out.append(full)
    return out


def _parse_chained(s, body: str):
    """(method [, opts]) varlist (method) varlist ... [= covariates]"""
    covars = []
    if "=" in body:
        body, rhs = body.rsplit("=", 1)
        covars = s.expand_varlist(rhs) if rhs.strip() else []
    specs = []
    i = 0
    text = body.strip()
    while i < len(text):
        if text[i].isspace():
            i += 1
            continue
        if text[i] != "(":
            raise StataError(198, "mi impute chained: (method) expected")   # VERIFY
        j = text.index(")", i)
        inner = text[i + 1:j].strip()
        mopt = {}
        if "," in inner:
            inner, op = inner.split(",", 1)
            mopt = match_options(op, {"knn": 3, "augment": 3, "noconstant": 4, "omit": 4, "include": 4})
        mth = next((m for m, k in _METHODS.items() if len(inner.strip()) >= k and m.startswith(inner.strip())), None)
        if mth is None or mth == "chained":
            raise StataError(198, f"{inner}: invalid method")   # VERIFY
        k = text.find("(", j + 1)
        names = s.expand_varlist(text[j + 1:k if k != -1 else len(text)])
        for v in names:
            specs.append((mth, v, mopt))
        i = k if k != -1 else len(text)
    return specs, covars


def _design_rows(ds, covars, rows):
    if not covars:
        return np.zeros((len(rows), 0))
    from . import fvars
    terms = fvars.expand_fv(ds, " ".join(covars))
    mask = np.zeros(ds.nobs, dtype=bool)
    mask[rows] = True
    # build_design works on the mask in row order: map back to the order of rows
    cols = fvars.build_design(ds, terms, mask, constant=True)
    X = np.column_stack([c.values for c in cols if c.values is not None]) if cols else np.zeros((len(rows), 0))
    X = np.where(X >= SYS, np.nan, X)
    order = np.argsort(rows, kind="stable")
    inv = np.empty_like(order)
    inv[order] = np.arange(len(rows))
    return X[inv]


def _chained(specs, work, missing, cov, burn):
    names = [v for _m, v, _o in specs]
    # initialization: monotone order (fewest missing first), each with the variables before it
    order = sorted(range(len(specs)), key=lambda i: int(missing[specs[i][1]].sum()))
    done: list[str] = []
    for i in order:
        mth, v, so = specs[i]
        _cond(mth, v, so, work, missing, cov, done)
        done.append(v)
    for _ in range(burn):
        for i in order:
            mth, v, so = specs[i]
            _cond(mth, v, so, work, missing, cov, [w for w in names if w != v])


def _cond(mth, v, so, work, missing, cov, preds, o=None, ds=None, rows_m=None):
    mis = missing[v]
    if not mis.any():
        return
    cols = [cov] + [work[w][:, None] for w in preds]
    X = np.column_stack(cols + [np.ones((len(mis), 1))])
    ok = ~mis & np.all(np.isfinite(X), axis=1) & np.all(X < SYS, axis=1)
    fill = mis & np.all(np.isfinite(X), axis=1) & np.all(X < SYS, axis=1)
    if not fill.any():
        return
    if ds is not None:
        so = _limits(ds, so, rows_m, ok, fill)
    vals = _impute_one(mth, work[v][ok], X[ok], X[fill], so)
    w = work[v].copy()
    w[fill] = vals
    work[v] = w


def _rows_of(ds, style, m) -> np.ndarray:
    mm = _num(ds, "_mi_m")
    if style == "flong" or m == 0:
        return np.nonzero(mm == m)[0]
    ids = _num(ds, "_mi_id")
    base = np.nonzero(mm == 0)[0]
    imp = np.nonzero(mm == m)[0]
    repl = dict(zip(ids[imp].tolist(), imp.tolist()))
    return np.array([repl.get(i, r) for i, r in zip(ids[base].tolist(), base.tolist())], dtype=int)


def _add_rows(ds, style, M_old, M_new) -> None:
    mm = _num(ds, "_mi_m")
    base = np.nonzero(mm == 0)[0]
    if style == "flong":
        src = base
    else:
        # incomplete observations in any registered imputed variable
        miss = np.zeros(len(base), dtype=bool)
        for v in _list(ds, "_mi_ivars"):
            x = ds.get(v)
            miss |= (x.data[base] == "") if x.is_string else (x.data[base] >= SYS)
        src = base[miss]
    reps = [src] * (M_new - M_old)
    if not reps or not len(src):
        return
    take = np.concatenate([np.arange(ds.nobs)] + reps)
    n_old = ds.nobs
    for v in ds.vars:
        v.raw = v.raw[take]
    ds.nobs = len(take)
    mvals = _num(ds, "_mi_m")
    k = n_old
    for m in range(M_old + 1, M_new + 1):
        mvals[k:k + len(src)] = m
        k += len(src)
    ds.get("_mi_m").data = mvals
    if ds.esample is not None:
        ds.esample = None
    ds.sortlist = []
