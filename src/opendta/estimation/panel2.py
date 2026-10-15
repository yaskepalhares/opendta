"""xtpoisson (fe, re) and xtlogit (fe, re) ([XT] xtpoisson, [XT] xtlogit).

- xtpoisson, fe: conditional likelihood (multinomial within the group);
- xtpoisson, re: gamma effect, closed-form likelihood;
- xtlogit, fe: conditional logit (sum over the subsets with the same
  total, by dynamic programming);
- xtlogit, re: normal effect integrated by adaptive Gauss-Hermite with 12
  points (Stata's default is mvaghermite with 12 points; details of the
  adaptation VERIFY).

Layouts checked against case 0804."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lazy import lazy
from . import fvars
from .models2 import _Prep, num_hess
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x, fmt="%10.0g"):
    if x is None or x >= SYS or not np.isfinite(x):
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _num_grad(f, b, h=1e-5):
    g = np.zeros_like(b)
    for i in range(len(b)):
        e = np.zeros_like(b)
        hh = h * (abs(b[i]) + 1)
        e[i] = hh
        g[i] = (f(b + e) - f(b - e)) / (2 * hh)
    return g


class _Panel:
    def __init__(self, s: "Session", cmd: str, args: str):
        from .tsops import panel_info
        self.s = s
        self.pr = _Prep(s, args, {"fe": 2, "re": 2, "pa": 2, "i": 1, "or": 2, "irr": 3, "exposure": 3,
                                  "offset": 3, "normal": 4, "intpoints": 4, "intmethod": 4})
        o = self.pr.o
        if o.get("pa"):
            raise StataError(198, f"{cmd}, pa not implemented")   # VERIFY
        ds = s.data
        pvar = str(o["i"]).strip() if isinstance(o.get("i"), str) else panel_info(ds)
        if not pvar:
            raise StataError(459, "panel variable not set; use xtset varname ...")
        self.pvar = pvar
        ids = ds.get(pvar).data[self.pr.smp.mask]
        self.uniq, self.inv = np.unique(ids, return_inverse=True)
        self.fe = bool(o.get("fe"))
        self.y = self.pr.y.astype(float)


def _obs_per_group(inv: np.ndarray):
    cnt = np.bincount(inv)
    return int(cnt.min()), float(cnt.mean()), int(cnt.max())


def _header(s, title: str, sub: str, x: dict) -> None:
    out = s.output
    out.write("\n", "text")
    rows = [(title, ("Number of obs", f"{x['N']:,}")), (f"Group variable: {x['pvar']}", ("Number of groups", f"{x['G']:,}"))]
    for lft, (rl, rv) in rows:
        out.write(f"{lft:<48}{rl:<18}= ", "text")
        out.write(f"{rv:>10}\n", "result")
    out.write("\n", "text")
    out.write(f"{sub:<48}Obs per group:\n", "text")
    for lab, val in (("min", f"{x['tmin']:,}"), ("avg", f"{x['tavg']:.1f}"), ("max", f"{x['tmax']:,}")):
        out.write(f"{'':<62}{lab} = ", "text")
        out.write(f"{val:>10}\n", "result")
    out.write("\n", "text")
    if x.get("intpoints"):
        # xtlogit, re: integration method and points (compat 0804)
        out.write(f"{'Integration method: mvaghermite':<48}{'Integration pts.':<18}= ", "text")
        out.write(f"{x['intpoints']:>10}\n", "result")
        out.write("\n", "text")
    chi = f"{x['chi2type']} chi2({x['df_m']})"
    out.write(f"{'':<48}{chi:<18}= ", "text")
    out.write(f"{_g(x['chi2'], '%10.2f'):>10}\n", "result")
    lab = "Log likelihood  = "
    ll_txt = format(_g(x["ll"]).strip(), ">10")     # "Log likelihood  =  -394.4428" (compat 0804)
    out.write(f"{lab + ll_txt:<48}{'Prob > chi2':<18}= ", "text")
    out.write(f"{_g(x['p'], '%10.4f'):>10}\n", "result")
    out.write("\n", "text")


def _finish(s, P: _Panel, cmd: str, args: str, b, V, res, names_extra, title, sub, keep_groups, ll0=None,
            chibar=None, aux_rows=(), eform="", lr0=None, intpoints=None):
    pr = P.pr
    smp = pr.smp
    k = pr.X.shape[1] + (0 if P.fe else (1 if pr.constant else 0))
    slopes = list(range(pr.X.shape[1]))
    bs = b[slopes]
    Vs = V[np.ix_(slopes, slopes)]
    try:
        chi2 = float(bs @ np.linalg.solve(Vs, bs)) if slopes else SYS
    except np.linalg.LinAlgError:
        chi2 = SYS
    chi2type = "Wald"
    if lr0 is not None:
        # xtlogit, fe: LR test against the model without covariates (compat 0804)
        chi2, chi2type = max(2 * (res.ll - lr0), 0.0), "LR"
    p = chi2_tail(chi2, len(slopes)) if chi2 < SYS else SYS
    names = [c.name for c in smp.cols] + ([] if P.fe or not pr.constant else ["_cons"])
    from .models2 import _bmap, _expand, _semap
    kx = pr.X.shape[1] + (0 if P.fe or not pr.constant else 1)
    if P.fe:
        bmap = {}
        j = 0
        for c in smp.cols:
            bmap[c.name] = float(b[j]) if c.values is not None else 0.0
            j += 1 if c.values is not None else 0
        semap = {}
        j = 0
        for c in smp.cols:
            semap[c.name] = float(np.sqrt(max(V[j, j], 0))) if c.values is not None else 0.0
            j += 1 if c.values is not None else 0
        rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=False)
    else:
        rows = rows_from_columns(smp.cols, _bmap(pr, b, kx), _semap(pr, V, kx), s.data, constant=pr.constant)
    rows += list(aux_rows)
    layout = [(c.name if not c.omitted else fvars.omitted_name(c.name), c.values is not None) for c in smp.cols]
    if not P.fe and pr.constant:
        layout.append(("_cons", True))
    for nm in names_extra:
        layout.append((nm, True))
    b_full, V_full = _expand(pr, b, V, layout)
    all_names = [nm for nm, _ in layout]
    eqs = [smp.depname] * (len(all_names) - len(names_extra)) + [nm.strip("/") for nm in names_extra]
    est = Estimates(cmd, smp.depname, [nm if not nm.startswith("/") else "_cons" for nm in all_names],
                    b_full, V_full, int(keep_groups.sum()), "z", None)
    est.rows = rows
    est.eqnames = eqs
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, (not P.fe and pr.constant), smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = ""
    inv = P.inv[keep_groups]
    tmin, tavg, tmax = _obs_per_group(np.unique(inv, return_inverse=True)[1])
    x = {"N": int(keep_groups.sum()), "G": len(np.unique(inv)), "pvar": P.pvar, "tmin": tmin, "tavg": tavg,
         "tmax": tmax, "chi2": chi2, "p": p, "df_m": len(slopes), "ll": res.ll, "chi2type": chi2type,
         "intpoints": intpoints}
    est.extra = x
    est.coef_title = eform or "Coef."
    est.eform_on = bool(eform)

    def disp(ss, ee, *, header=True, table=True, level=None, **_kw):
        if header:
            _header(ss, title, sub, ee.extra)
        if table:
            coef_table(ss, ee.depvar, ee.rows, stat="z", df=None, level=level or ee.level,
                       coef_title=ee.coef_title, eform=ee.eform_on)
            if chibar is not None:
                lab, val, pv = chibar
                # "LR test of rho=0: chibar2(01) = 8.21" with Prob at column 56 (compat 0804)
                ss.output.write(f"{lab}: chibar2(01) = {val:.2f}".ljust(55)
                                + f"Prob >= chibar2 = {pv:5.3f}\n", "text")
    est.display = disp
    scalars = [("N", float(x["N"])), ("N_g", float(x["G"])), ("g_min", float(tmin)), ("g_avg", tavg),
               ("g_max", float(tmax)), ("df_m", float(len(slopes))), ("ll", res.ll), ("chi2", chi2), ("p", p),
               ("rank", float(np.linalg.matrix_rank(V))), ("converged", 1.0 if res.converged else 0.0),
               ("ic", float(res.iterations))]
    if ll0 is not None:
        scalars.append(("ll_0", ll0))
    elif lr0 is not None:
        scalars.append(("ll_0", lr0))
    if chibar is not None:
        scalars += [("chi2_c", chibar[1]), ("p_c", chibar[2])]
    macros = [("cmdline", f"{cmd} {args.strip()}"), ("cmd", cmd), ("ivar", P.pvar), ("model", "fe" if P.fe else "re"),
              ("depvar", smp.depname), ("chi2type", chi2type), ("vce", "oim"), ("title", title),
              ("properties", "b V"), ("predict", f"{cmd}_p")]
    sample = np.zeros(s.data.nobs, dtype=bool)
    sample[np.flatnonzero(smp.mask)[keep_groups]] = True
    post(s, est, list(reversed(scalars)), list(reversed(macros)), sample)
    disp(s, est)


# ---------------------------------------------------------------------------
# xtpoisson
# ---------------------------------------------------------------------------

def cmd_xtpoisson(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtpoisson", args)
        return
    P = _Panel(s, "xtpoisson", args)
    pr = P.pr
    y, X, inv = P.y, pr.X, P.inv
    if np.any(y < 0) or np.any(y != np.trunc(y)):
        raise StataError(459, f"{pr.smp.depname} must be nonnegative integers")   # VERIFY
    out = s.output
    eform = "IRR" if pr.o.get("irr") else ""
    if P.fe:
        tot = np.bincount(inv, weights=y)
        keep = tot[inv] > 0
        drop_g = int((tot == 0).sum())
        if drop_g:
            out.write(f"note: {drop_g} group{'s' if drop_g > 1 else ''} ({int((~keep).sum())} obs) dropped because "
                      "of all zero outcomes\n", "text")
        yk, Xk = y[keep], X[keep]
        _, invk = np.unique(inv[keep], return_inverse=True)
        G = invk.max() + 1
        nk = np.bincount(invk, weights=yk)
        const = float(special.gammaln(nk + 1).sum() - special.gammaln(yk + 1).sum())

        def fun(b, todo):
            eta = Xk @ b
            mx = np.full(G, -np.inf)
            np.maximum.at(mx, invk, eta)
            e = np.exp(eta - mx[invk])
            den = np.bincount(invk, weights=e, minlength=G)
            lse = np.log(den) + mx
            ll = float(yk @ eta - nk @ lse) + const
            if todo == 0:
                return ll, None, None
            pi = e / den[invk]
            xbar = np.zeros((G, Xk.shape[1]))
            np.add.at(xbar, invk, Xk * pi[:, None])
            g = Xk.T @ yk - nk @ xbar
            Xc = Xk - xbar[invk]
            H = -(Xc * (pi * nk[invk])[:, None]).T @ Xc
            return ll, g, H
        res = pr.maximize(fun, np.zeros(X.shape[1]))
        if not pr.o.get("nolog"):
            from .ml import print_log
            out.write("\n", "text")
            print_log(out, res, "log likelihood")
        V = np.linalg.inv(-res.H)
        _finish(s, P, "xtpoisson", args, res.b, V, res, [], "Conditional fixed-effects Poisson regression", "",
                keep, eform=eform)
        return
    # gamma random effects
    Xc = np.column_stack([X, np.ones(len(y))]) if pr.constant else X
    k = Xc.shape[1]
    lgy = special.gammaln(y + 1)
    G = inv.max() + 1
    n_i = np.bincount(inv, weights=y, minlength=G)

    def ll_re(par):
        b, lna = par[:k], par[k]
        eta = Xc @ b
        u = np.exp(np.clip(eta, -700, 700))
        th = math.exp(-lna)
        U = np.bincount(inv, weights=u, minlength=G)
        return float((y * eta - lgy).sum() + (special.gammaln(th + n_i) - special.gammaln(th)
                                               + th * math.log(th) - (th + n_i) * np.log(th + U)).sum())

    def fun(par, todo):
        ll = ll_re(par)
        if todo == 0:
            return ll, None, None
        g = _num_grad(ll_re, par)
        H = num_hess(lambda q: _num_grad(ll_re, q), par)
        return ll, g, H
    # start: pooled Poisson and ln(alpha) = 0
    from .models import _poisson
    pf, _ = _poisson(y, Xc, np.ones(len(y)), np.zeros(len(y)))
    from .ml import maximize, print_log
    b0 = np.zeros(k)
    if pr.constant:
        b0[-1] = math.log(max(y.mean(), 1e-8))
    rp = maximize(pf, b0)
    # Stata shows the Poisson starting already at the optimum (iterations 0 and 1 equal; compat
    # 0804): the exact fit stays out of the log
    rp = maximize(pf, rp.b)
    res = pr.maximize(fun, np.concatenate([rp.b, [0.0]]))
    if not pr.o.get("nolog"):
        out.write("\nFitting Poisson model:\n\n", "text")
        print_log(out, rp, "log likelihood")
        out.write("\nFitting full model:\n\n", "text")
        print_log(out, res, "log likelihood")
    V = np.linalg.inv(-res.H)
    lna, se_lna = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    alpha = math.exp(lna)
    from .results import crit
    c = crit("z", pr.smp.level, None)
    aux = [CoefRow("sep"), CoefRow("aux", "/lnalpha", lna, se_lna, name="_cons"), CoefRow("sep"),
           CoefRow("aux", "alpha", alpha, alpha * se_lna,
                   eq=f"{math.exp(lna - c * se_lna)} {math.exp(lna + c * se_lna)}")]
    chib = max(2 * (res.ll - rp.ll), 0.0)
    pchib = 0.5 * chi2_tail(chib, 1) if chib > 0 else 0.5
    _finish(s, P, "xtpoisson", args, res.b, V, res, ["/lnalpha"], "Random-effects Poisson regression",
            "Random effects u_i ~ Gamma", np.ones(len(y), bool), chibar=("LR test of alpha=0", chib, pchib),
            aux_rows=aux, eform=eform)


# ---------------------------------------------------------------------------
# xtlogit
# ---------------------------------------------------------------------------

def _clogit_ll(eta: np.ndarray, y: np.ndarray) -> float:
    """log P(y | Σy) of the conditional logit for one group."""
    k = int(y.sum())
    T = len(eta)
    # f[j] = log Σ_{|S|=j} exp(Σ_S eta), over the observations seen so far
    f = np.full(k + 1, -np.inf)
    f[0] = 0.0
    for t in range(T):
        upto = min(t + 1, k)
        for j in range(upto, 0, -1):
            f[j] = np.logaddexp(f[j], f[j - 1] + eta[t])
    return float(eta @ y - f[k])


_GH_X, _GH_W = np.polynomial.hermite.hermgauss(12)


def cmd_xtlogit(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtlogit", args)
        return
    P = _Panel(s, "xtlogit", args)
    pr = P.pr
    y = (P.y != 0).astype(float)
    X, inv = pr.X, P.inv
    out = s.output
    eform = "OR" if pr.o.get("or") else ""
    from .ml import print_log
    if P.fe:
        tot = np.bincount(inv, weights=y)
        cnt = np.bincount(inv)
        good = (tot > 0) & (tot < cnt)
        keep = good[inv]
        dg = int((~good).sum())
        if np.any(tot[good] > 1):
            # this note comes before the one about dropped groups (compat 0804)
            out.write("note: multiple positive outcomes within groups encountered.\n", "text")
        if dg:
            out.write(f"note: {dg} group{'s' if dg > 1 else ''} ({int((~keep).sum())} obs) dropped because of "
                      "all positive or\n      all negative outcomes.\n", "text")
        groups = [np.flatnonzero((inv == g)) for g in np.flatnonzero(good)]

        def ll_fe(b):
            eta = X @ b
            return sum(_clogit_ll(eta[ix], y[ix]) for ix in groups)

        def fun(b, todo):
            ll = ll_fe(b)
            if todo == 0:
                return ll, None, None
            g = _num_grad(ll_fe, b)
            return ll, g, num_hess(lambda q: _num_grad(ll_fe, q), b)
        # start: pooled logit (with a constant) on the groups kept (compat 0804)
        from .ml import maximize
        from .models import _logit
        Xk = np.column_stack([X[keep], np.ones(int(keep.sum()))])
        lk, _ = _logit(y[keep], Xk, np.ones(int(keep.sum())))
        b_start = maximize(lk, np.zeros(Xk.shape[1])).b[:-1]
        res = pr.maximize(fun, b_start)
        if not pr.o.get("nolog"):
            out.write("\n", "text")
            print_log(out, res, "log likelihood")
        V = np.linalg.inv(-res.H)
        _finish(s, P, "xtlogit", args, res.b, V, res, [], "Conditional fixed-effects logistic regression", "",
                keep, eform="Odds Ratio" if eform else "", lr0=ll_fe(np.zeros(X.shape[1])))
        return
    Xc = np.column_stack([X, np.ones(len(y))]) if pr.constant else X
    k = Xc.shape[1]
    G = inv.max() + 1
    sign = 2 * y - 1
    order = np.argsort(inv, kind="stable")
    bounds = np.searchsorted(inv[order], np.arange(G + 1))

    def ll_re(par):
        b, lns2 = par[:k], par[k]
        sig = math.exp(0.5 * lns2)
        eta = Xc @ b
        tot = 0.0
        for g in range(G):
            ix = order[bounds[g]:bounds[g + 1]]
            e, sg = eta[ix], sign[ix]
            # mode and curvature of the integrand in nu (adaptive)
            mu, tau = 0.0, sig
            for _ in range(8):
                p = special.expit(sg * (e + mu))
                d1 = float((sg * (1 - p)).sum()) - mu / sig ** 2
                d2 = -float((p * (1 - p)).sum()) - 1 / sig ** 2
                mu -= d1 / d2
                tau = 1 / math.sqrt(-d2)
                if abs(d1) < 1e-10:
                    break
            nodes = mu + math.sqrt(2) * tau * _GH_X
            vals = (special.log_expit(sg[None, :] * (e[None, :] + nodes[:, None]))).sum(axis=1) \
                - 0.5 * (nodes / sig) ** 2 - math.log(sig * math.sqrt(2 * math.pi))
            terms = np.log(_GH_W) + _GH_X ** 2 + vals
            tot += math.log(math.sqrt(2) * tau) + float(special.logsumexp(terms))
        return tot

    def fun(par, todo):
        ll = ll_re(par)
        if todo == 0:
            return ll, None, None
        g = _num_grad(ll_re, par)
        return ll, g, num_hess(lambda q: _num_grad(ll_re, q), par)
    from .ml import maximize
    from .models import _logit
    lf, _ = _logit(y, Xc, np.ones(len(y)))
    b0 = np.zeros(k)
    if pr.constant:
        p0 = y.mean()
        b0[-1] = math.log(p0 / (1 - p0))
    rl = maximize(lf, b0)
    # search over tau = sigma2/(1 + sigma2) with the pooled logit coefficients and plain
    # Gauss-Hermite quadrature, until the likelihood falls; the full model starts from the best tau
    # (compat 0804)

    def ll_grid(t):
        if t == 0:
            return rl.ll
        sig = math.sqrt(t / (1 - t))
        eta = Xc @ rl.b
        tot = 0.0
        for g in range(G):
            ix = order[bounds[g]:bounds[g + 1]]
            e, sg = eta[ix], sign[ix]
            nodes = math.sqrt(2) * sig * _GH_X
            vals = special.log_expit(sg[None, :] * (e[None, :] + nodes[:, None])).sum(axis=1)
            tot += float(special.logsumexp(np.log(_GH_W) + vals)) - 0.5 * math.log(math.pi)
        return tot
    grid = []
    t = 0.0
    while t < 0.95:
        grid.append((t, ll_grid(t)))
        if len(grid) > 1 and grid[-1][1] < grid[-2][1]:
            break
        t = round(t + 0.1, 1)
    best_t = max(grid, key=lambda tv: tv[1])[0]
    start_l2 = math.log(best_t / (1 - best_t)) if best_t > 0 else 0.0
    res = pr.maximize(fun, np.concatenate([rl.b, [start_l2]]))
    if not pr.o.get("nolog"):
        out.write("\nFitting comparison model:\n\n", "text")
        print_log(out, rl, "log likelihood")
        out.write("\nFitting full model:\n\n", "text")
        for tv, lv in grid:
            out.write(f"tau = {tv:4.1f}     log likelihood = {_g(lv).strip()}\n", "text")
        out.write("\n", "text")
        print_log(out, res, "log likelihood")
    V = np.linalg.inv(-res.H)
    l2, se2 = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    sig = math.exp(0.5 * l2)
    rho = sig ** 2 / (sig ** 2 + math.pi ** 2 / 3)
    from .results import crit
    c = crit("z", pr.smp.level, None)
    lo2, hi2 = l2 - c * se2, l2 + c * se2

    def rho_of(v):
        s2 = math.exp(v)
        return s2 / (s2 + math.pi ** 2 / 3)
    aux = [CoefRow("sep"), CoefRow("aux", "/lnsig2u", l2, se2, name="_cons"), CoefRow("sep"),
           CoefRow("aux", "sigma_u", sig, 0.5 * sig * se2,
                   eq=f"{math.exp(0.5 * lo2)} {math.exp(0.5 * hi2)}"),
           CoefRow("aux", "rho", rho, rho * (1 - rho) * se2, eq=f"{rho_of(lo2)} {rho_of(hi2)}")]
    chib = max(2 * (res.ll - rl.ll), 0.0)
    pchib = 0.5 * chi2_tail(chib, 1) if chib > 0 else 0.5
    _finish(s, P, "xtlogit", args, res.b, V, res, ["/lnsig2u"], "Random-effects logistic regression",
            "Random effects u_i ~ Gaussian", np.ones(len(y), bool), chibar=("LR test of rho=0", chib, pchib),
            aux_rows=aux, eform="Odds Ratio" if eform else "", intpoints=12)
