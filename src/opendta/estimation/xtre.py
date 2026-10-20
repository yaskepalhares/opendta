"""Panel models with a random effect integrated by quadrature, and xtnbreg: xtprobit, xtcloglog,
xtologit, xtoprobit, xttobit, xtintreg, xtnbreg and quadchk ([XT] xtprobit, [XT] xtcloglog,
[XT] xtologit, [XT] xtoprobit, [XT] xttobit, [XT] xtintreg, [XT] xtnbreg, [XT] quadchk).

The likelihoods are in core/panelquad.py; this module reads the data, maximizes, writes the
output and stores e(). Layouts follow xtlogit (compat 0804); where the manual does not show
the layout the code carries a VERIFY mark and the case is in compat 1905."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import panelquad as pq
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lazy import lazy
from .ml import maximize, print_log
from .models2 import num_hess
from .panel2 import _finish, _num_grad, _Panel
from .results import CoefRow, crit

special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

# the part of rho that is not the random effect: logistic, normal and extreme-value errors
_RHO_C = {"logit": math.pi ** 2 / 3, "probit": 1.0, "cloglog": math.pi ** 2 / 6}


# ---------------------------------------------------------------------------
# options and common pieces
# ---------------------------------------------------------------------------

def _quad_options(P: _Panel) -> tuple[int, bool, str]:
    o = P.pr.o
    nq = 12
    if o.get("intpoints") not in (None, False, True):
        try:
            nq = int(str(o["intpoints"]).strip())
        except ValueError:
            raise StataError(198, "option intpoints() incorrectly specified")
        if nq < 1 or nq > 195:
            raise StataError(198, "intpoints() must be between 1 and 195")      # VERIFY
    meth = str(o["intmethod"]).strip().lower() if o.get("intmethod") not in (None, False, True) else "mvaghermite"
    if meth and "mvaghermite".startswith(meth):
        return nq, True, "mvaghermite"
    if meth and "ghermite".startswith(meth):
        return nq, False, "ghermite"
    raise StataError(198, "option intmethod() incorrectly specified")


def _exp(v: float) -> float:
    return math.exp(min(max(v, -700.0), 700.0))


def _no_fe(P: _Panel, cmd: str) -> None:
    o = P.pr.o
    if o.get("fe"):
        raise StataError(198, f"option fe not allowed")      # VERIFY


def _numeric_fun(ll):
    """fun(par, todo) for maximize(), with numerical derivatives."""
    def fun(par, todo):
        v = ll(par)
        if not np.isfinite(v):
            return -np.inf, None, None
        if todo == 0:
            return v, None, None
        g = _num_grad(ll, par)
        return v, g, num_hess(lambda q: _num_grad(ll, q), par)
    return fun


def _robust_vce(llg, par, H, grp, cluster):
    """Sandwich for the clusters (default: the panels): scores of the group likelihoods by
    central differences, G/(G-1) as the small-sample factor. Returns V and the number of clusters."""
    p = len(par)
    S = np.zeros((grp.G, p))
    for j in range(p):
        h = 1e-5 * (abs(par[j]) + 1)
        a, c = par.copy(), par.copy()
        a[j] += h
        c[j] -= h
        S[:, j] = (llg(a) - llg(c)) / (2 * h)
    if cluster is not None:
        codes = np.unique(cluster, return_inverse=True)[1]
        first = np.zeros(grp.G, dtype=int)
        first[grp.inv[::-1]] = np.arange(len(grp.inv))[::-1]
        gc = codes[first]
        M_ = int(gc.max()) + 1
        S2 = np.zeros((M_, p))
        np.add.at(S2, gc, S)
        S = S2
    else:
        M_ = grp.G
    A = np.linalg.inv(-H)
    meat = S.T @ S * (M_ / (M_ - 1) if M_ > 1 else 1.0)
    return A @ meat @ A, M_


def _vce(P: _Panel, llg, res, grp, plain=None):
    """(V, robust) with robust = (clustvar, N clusters) or None."""
    smp = P.pr.smp
    V = np.linalg.inv(-res.H) if plain is None else plain
    if not P.pr.robust:
        return V, None
    cl = smp.cluster if smp.vce == "cluster" else None
    Vr, M_ = _robust_vce(llg, res.b, res.H, grp, cl)
    return Vr, (smp.clustvar if smp.vce == "cluster" else P.pvar, M_)


def _drop_note(s, P, y_is_const: bool) -> None:
    pass


def _remember(s: "Session", cmd: str, args: str, nq: int, adaptive: bool, ll: float, est_b) -> None:
    s._quadchk = {"cmd": cmd, "args": args, "nq": nq, "adaptive": adaptive, "ll": ll}


def _store(s, **kv):
    for k, v in kv.items():
        s.e[k] = v


# ---------------------------------------------------------------------------
# xtprobit, xtcloglog
# ---------------------------------------------------------------------------

def _pooled_binary(kind: str, y: np.ndarray, Xc: np.ndarray, pr):
    if kind == "probit":
        from .models import _probit
        fn, _ = _probit(y, Xc, np.ones(len(y)))
        return fn
    # cloglog: Pr(y=1) = 1 - exp(-exp(xb))
    def grad(b):
        e = np.exp(np.clip(Xc @ b, -700, 700))
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(y != 0, e * np.exp(-e) / -np.expm1(-e), -e)
        return Xc.T @ np.nan_to_num(r)

    def fn(b, todo):
        v = float(pq.ll_cloglog(Xc @ b, y).sum())
        if todo == 0:
            return v, None, None
        return v, grad(b), num_hess(grad, b)
    return fn


def _binary_re(s: "Session", cmd: str, args: str, kind: str, title: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return
    from .xtgee import run_pa
    if run_pa(s, cmd, args):
        return
    P = _Panel(s, cmd, args)
    _no_fe(P, cmd)
    pr = P.pr
    nq, adaptive, meth = _quad_options(P)
    y = (P.y != 0).astype(float)
    if y.min() == y.max():
        raise StataError(2000, f"{pr.smp.depname} does not vary")         # VERIFY
    Xc = np.column_stack([pr.X, np.ones(len(y))]) if pr.constant else pr.X
    k = Xc.shape[1]
    grp = pq.Groups(P.inv)
    sign = 2 * y - 1
    if kind == "probit":
        def f(u):
            return pq.ll_probit(u, sign if u.ndim == 1 else sign[:, None])
    else:
        def f(u):
            return pq.ll_cloglog(u, y if u.ndim == 1 else y[:, None])

    def llg(par):
        return pq.re_loglik_groups(f, Xc @ par[:k], grp, math.exp(0.5 * par[k]), nq, adaptive)

    def ll(par):
        return float(llg(par).sum())
    fun = _numeric_fun(ll)

    out = s.output
    pooled = _pooled_binary(kind, y, Xc, pr)
    b0 = np.zeros(k)
    if pr.constant:
        p0 = y.mean()
        b0[-1] = st.norm.ppf(p0) if kind == "probit" else math.log(-math.log(1 - p0))
    rl = maximize(pooled, b0)
    grid = pq.tau_grid(f, Xc @ rl.b, grp, rl.ll, nq)
    best_t = max(grid, key=lambda tv: tv[1])[0]
    start_l2 = math.log(best_t / (1 - best_t)) if best_t > 0 else 0.0
    res = pr.maximize(fun, np.concatenate([rl.b, [start_l2]]))
    if not pr.o.get("nolog"):
        out.write("\nFitting comparison model:\n\n", "text")
        print_log(out, rl, "log likelihood")
        out.write("\nFitting full model:\n\n", "text")
        lab = "rho" if kind == "probit" else "tau"           # VERIFY (xtcloglog)
        for tv, lv in grid:
            out.write(f"{lab} = {tv:4.1f}     log likelihood = {_g(lv)}\n", "text")
        out.write("\n", "text")
        print_log(out, res, "log likelihood")
    V, robust = _vce(P, llg, res, grp)
    l2, se2 = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    sig = math.exp(0.5 * l2)
    c = _RHO_C[kind]
    rho = sig ** 2 / (sig ** 2 + c)
    z = crit("z", pr.smp.level, None)
    lo2, hi2 = l2 - z * se2, l2 + z * se2

    def rho_of(v):
        s2 = _exp(v)
        return s2 / (s2 + c)
    aux = [CoefRow("sep"), CoefRow("aux", "/lnsig2u", l2, se2, name="_cons"), CoefRow("sep"),
           CoefRow("aux", "sigma_u", sig, 0.5 * sig * se2, eq=f"{_exp(0.5 * lo2)} {_exp(0.5 * hi2)}"),
           CoefRow("aux", "rho", rho, rho * (1 - rho) * se2, eq=f"{rho_of(lo2)} {rho_of(hi2)}")]
    chib = max(2 * (res.ll - rl.ll), 0.0)
    pchib = 0.5 * chi2_tail(chib, 1) if chib > 0 else 0.5
    _finish(s, P, cmd, args, res.b, V, res, ["/lnsig2u"], title, "Random effects u_i ~ Gaussian",
            np.ones(len(y), bool), chibar=("LR test of rho=0", chib, pchib), aux_rows=aux, intpoints=nq,
            intmethod=meth, robust=robust,
            extra_scalars=[("sigma_u", sig), ("rho", rho), ("ll_c", rl.ll), ("chi2_c", chib), ("p_c", pchib)],
            extra_macros=[("distrib", "Gaussian")])
    _remember(s, cmd, args, nq, adaptive, res.ll, res.b)


def _g(x):
    from ..core.formats import format_value
    return format_value(float(x), "%10.0g", pad=False).strip()


def cmd_xtprobit(s: "Session", args: str) -> None:
    _binary_re(s, "xtprobit", args, "probit", "Random-effects probit regression")


def cmd_xtcloglog(s: "Session", args: str) -> None:
    _binary_re(s, "xtcloglog", args, "cloglog", "Random-effects complementary log-log regression")


# ---------------------------------------------------------------------------
# xtologit, xtoprobit
# ---------------------------------------------------------------------------

def _ordered_re(s: "Session", cmd: str, args: str, logistic: bool) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return
    P = _Panel(s, cmd, args, constant=False)
    _no_fe(P, cmd)
    pr = P.pr
    out = s.output
    pr.write_notes(out)
    nq, adaptive, meth = _quad_options(P)
    y = P.y
    levels = sorted(set(y.tolist()))
    J = len(levels)
    if J < 2:
        raise StataError(2000, "outcome does not vary")                   # VERIFY
    yi = np.searchsorted(levels, y)
    X = pr.X
    k = X.shape[1]
    grp = pq.Groups(P.inv)
    Finv = special.logit if logistic else st.norm.ppf

    def llg(par):
        cuts = par[k:k + J - 1]
        if np.any(np.diff(cuts) <= 0):
            return np.full(grp.G, -np.inf)
        f = lambda u: pq.ll_ordered(u, yi, cuts, logistic)         # noqa: E731
        return pq.re_loglik_groups(f, X @ par[:k], grp, math.exp(0.5 * par[-1]), nq, adaptive)

    def ll(par):
        return float(llg(par).sum())

    def ll_pooled(b):
        cuts = b[k:]
        if np.any(np.diff(cuts) <= 0):
            return -np.inf
        return float(pq.ll_ordered(X @ b[:k], yi, cuts, logistic).sum())
    cum = np.array([np.mean(yi <= j) for j in range(J - 1)])
    rl = maximize(_numeric_fun(ll_pooled), np.concatenate([np.zeros(k), Finv(cum)]))
    # one-dimensional search for the variance with the pooled coefficients
    base = rl.b

    def ll_var(v):
        return ll(np.concatenate([base, v]))
    rv = maximize(_numeric_fun(ll_var), np.array([-1.0]))
    res = pr.maximize(_numeric_fun(ll), np.concatenate([base, rv.b]))
    if not pr.o.get("nolog"):
        out.write("\nFitting comparison model:\n\n", "text")
        print_log(out, rl, "log likelihood")
        out.write("\nRefining starting values:\n\n", "text")                   # VERIFY
        print_log(out, rv, "log likelihood")
        out.write("\nFitting full model:\n\n", "text")
        print_log(out, res, "log likelihood")
    V, robust = _vce(P, llg, res, grp)
    # sigma2_u is reported in levels (Stata estimates it directly)
    l2 = float(res.b[-1])
    s2 = math.exp(l2)
    Jm = np.eye(len(res.b))
    Jm[-1, -1] = s2
    Vl = Jm @ V @ Jm.T
    b_rep = res.b.copy()
    b_rep[-1] = s2
    se_s2 = float(np.sqrt(Vl[-1, -1]))
    z = crit("z", pr.smp.level, None)
    aux = [CoefRow("sep")]
    for j in range(J - 1):
        aux.append(CoefRow("aux", f"/cut{j + 1}", float(res.b[k + j]), float(np.sqrt(V[k + j, k + j])),
                           name=f"cut{j + 1}"))
    aux += [CoefRow("sep"), CoefRow("aux", "/sigma2_u", s2, se_s2, name="sigma2_u")]
    chib = max(2 * (res.ll - rl.ll), 0.0)
    pchib = 0.5 * chi2_tail(chib, 1) if chib > 0 else 0.5
    name = "ologit" if logistic else "oprobit"
    title = "Random-effects ordered logistic regression" if logistic else "Random-effects ordered probit regression"
    _finish(s, P, cmd, args, b_rep, Vl, res, [f"/cut{j + 1}" for j in range(J - 1)] + ["/sigma2_u"], title,
            "Random effects u_i ~ Gaussian", np.ones(len(y), bool),
            chibar=(f"LR test vs. {name} model", chib, pchib), aux_rows=aux, intpoints=nq, intmethod=meth,
            robust=robust, extra_scalars=[("sigma2_u", s2), ("ll_c", rl.ll), ("chi2_c", chib), ("p_c", pchib),
                                          ("k_cat", float(J))],
            extra_macros=[("distrib", "Gaussian")])
    _remember(s, cmd, args, nq, adaptive, res.ll, res.b)


def cmd_xtologit(s: "Session", args: str) -> None:
    _ordered_re(s, "xtologit", args, True)


def cmd_xtoprobit(s: "Session", args: str) -> None:
    _ordered_re(s, "xtoprobit", args, False)


# ---------------------------------------------------------------------------
# xttobit, xtintreg
# ---------------------------------------------------------------------------

def _limit(s, v, default):
    if v is None or v is False:
        return None
    if v is True or str(v).strip() == "":
        return default
    return float(s.eval(str(v)))


def _censored_re(s: "Session", cmd: str, args: str, interval: bool) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return
    ds = s.data
    opts = {"ll": 2, "ul": 2} if not interval else {}
    if interval:
        parts = args.split(None, 2)
        if len(parts) < 2 or parts[1].startswith(","):
            raise StataError(102, "too few variables specified")
        dep1, dep2 = parts[0], parts[1]
        rest = parts[2] if len(parts) > 2 else ""
        P = _Panel(s, cmd, rest if rest.strip() and not rest.lstrip().startswith(",") else " " + rest,
                   extra_opts=opts, need_dep=False)
        from ..core.varlist import resolve_name
        n1, n2 = resolve_name(ds, dep1), resolve_name(ds, dep2)
        mask = P.pr.smp.mask
        y1 = ds.get(n1).data[mask].astype(float)
        y2 = ds.get(n2).data[mask].astype(float)
        y1 = np.where(y1 >= SYS, np.nan, y1)
        y2 = np.where(y2 >= SYS, np.nan, y2)
        P.pr.smp.depname = n1
        P.pr.smp.depcomp = None
    else:
        P = _Panel(s, cmd, args, extra_opts=opts)
        y = P.y
    _no_fe(P, cmd)
    pr = P.pr
    out = s.output
    pr.write_notes(out)
    nq, adaptive, meth = _quad_options(P)
    if interval:
        keep = ~(np.isnan(y1) & np.isnan(y2))
        if np.any(np.isnan(y1) & np.isnan(y2)):
            pass
        lo = hi = None
    else:
        lo = _limit(s, pr.o.get("ll"), float(y.min()))
        hi = _limit(s, pr.o.get("ul"), float(y.max()))
        if lo is None and hi is None and not interval:
            pass
        keep = np.ones(len(y), bool)
    Xc = np.column_stack([pr.X, np.ones(len(keep))]) if pr.constant else pr.X
    Xc = Xc[keep]
    inv_all = P.inv[keep]
    _, inv = np.unique(inv_all, return_inverse=True)
    grp = pq.Groups(inv)
    k = Xc.shape[1]
    if interval:
        a1, a2 = y1[keep], y2[keep]
        left = np.isnan(a1)
        right = np.isnan(a2)
        point = ~left & ~right & (a1 == a2)
        intv = ~left & ~right & (a1 != a2)
        unc = point
        yy = None
    else:
        yy = y[keep]
        left = (yy <= lo) if lo is not None else np.zeros(len(yy), bool)
        right = (yy >= hi) if hi is not None else np.zeros(len(yy), bool)
        unc = ~(left | right)
        intv = np.zeros(len(yy), bool)

    def obs_f(sig_e):
        if interval:
            return lambda u: pq.ll_interval(u, a1, a2, sig_e)
        return lambda u: pq.ll_tobit(u, yy, sig_e, lo, hi)

    def llg(par):
        return pq.re_loglik_groups(obs_f(math.exp(par[k + 1])), Xc @ par[:k], grp, math.exp(par[k]), nq,
                                   adaptive)

    def ll(par):
        return float(llg(par).sum())

    def ll_pooled(b):
        return float(obs_f(math.exp(b[k]))(Xc @ b[:k]).sum())
    # starting values: least squares on the points (or the midpoints), then the pooled tobit
    if interval:
        mid = np.where(left, a2, np.where(right, a1, (a1 + a2) / 2))
    else:
        mid = yy
    beta0, *_ = np.linalg.lstsq(Xc, mid, rcond=None)
    sd0 = float(np.std(mid - Xc @ beta0)) or 1.0
    rl = maximize(_numeric_fun(ll_pooled), np.concatenate([beta0, [math.log(sd0)]]))
    se0 = math.exp(rl.b[k])
    # standard deviation of the effect: the best of a few values with the pooled coefficients
    cand = [se0 * m for m in (0.1, 0.25, 0.5, 1.0, 2.0)]
    vals = [ll(np.concatenate([rl.b[:k], [math.log(c), rl.b[k]]])) for c in cand]
    su0 = cand[int(np.argmax(vals))]
    start = np.concatenate([rl.b[:k], [math.log(su0), rl.b[k]]])
    res = pr.maximize(_numeric_fun(ll), start)
    if not pr.o.get("nolog"):
        out.write("\nFitting full model:\n\n", "text")                          # VERIFY
        print_log(out, res, "log likelihood")
    Vl, robust = _vce(P, llg, res, grp)
    su, se_ = math.exp(res.b[k]), math.exp(res.b[k + 1])
    Jm = np.eye(len(res.b))
    Jm[k, k], Jm[k + 1, k + 1] = su, se_
    V = Jm @ Vl @ Jm.T
    b_rep = res.b.copy()
    b_rep[k], b_rep[k + 1] = su, se_
    rho = su ** 2 / (su ** 2 + se_ ** 2)
    g_ = np.array([2 * rho * (1 - rho), -2 * rho * (1 - rho)])
    se_rho = float(math.sqrt(g_ @ Vl[k:k + 2, k:k + 2] @ g_))
    lr = math.log(rho / (1 - rho))
    se_lr = se_rho / (rho * (1 - rho))
    z = crit("z", pr.smp.level, None)
    aux = [CoefRow("sep"),
           CoefRow("aux", "/sigma_u", su, float(np.sqrt(V[k, k])), name="sigma_u"),
           CoefRow("aux", "/sigma_e", se_, float(np.sqrt(V[k + 1, k + 1])), name="sigma_e"),
           CoefRow("aux", "rho", rho, se_rho,
                   eq=f"{1 / (1 + math.exp(-(lr - z * se_lr)))} {1 / (1 + math.exp(-(lr + z * se_lr)))}")]
    chib = max(2 * (res.ll - rl.ll), 0.0)
    pchib = 0.5 * chi2_tail(chib, 1) if chib > 0 else 0.5
    counts = (int(left.sum()), int(unc.sum()), int(right.sum()), int(intv.sum()))

    def footer(ss, ee):
        nl, nu, nr, ni = counts
        w = ss.output
        w.write(f"Observation summary: {nl:>12,}  left-censored observations\n", "text")
        w.write(f"{nu:>33,}     uncensored observations\n", "text")
        w.write(f"{nr:>33,} right-censored observations\n", "text")
        if interval:
            w.write(f"{ni:>33,}     interval observations\n", "text")
    title = "Random-effects interval-data regression" if interval else "Random-effects tobit regression"
    sc = [("sigma_u", su), ("sigma_e", se_), ("rho", rho), ("ll_c", rl.ll), ("chi2_c", chib), ("p_c", pchib),
          ("N_lc", float(left.sum())), ("N_rc", float(right.sum())), ("N_unc", float(unc.sum()))]
    if interval:
        sc.append(("N_int", float(intv.sum())))
    else:
        if lo is not None:
            sc.append(("llopt", lo))
        if hi is not None:
            sc.append(("ulopt", hi))
    _finish(s, P, cmd, args, b_rep, V, res, ["/sigma_u", "/sigma_e"], title, "Random effects u_i ~ Gaussian",
            keep, chibar=("LR test of sigma_u=0", chib, pchib), aux_rows=aux, intpoints=nq, intmethod=meth,
            robust=robust, extra_scalars=sc, extra_macros=[("distrib", "Gaussian")], footer=footer)
    _remember(s, cmd, args, nq, adaptive, res.ll, res.b)


def cmd_xttobit(s: "Session", args: str) -> None:
    _censored_re(s, "xttobit", args, False)


def cmd_xtintreg(s: "Session", args: str) -> None:
    _censored_re(s, "xtintreg", args, True)


# ---------------------------------------------------------------------------
# xtnbreg
# ---------------------------------------------------------------------------

def cmd_xtnbreg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtnbreg", args)
        return
    from .xtgee import run_pa
    if run_pa(s, "xtnbreg", args):
        return
    P = _Panel(s, "xtnbreg", args)
    pr = P.pr
    y, inv = P.y, P.inv
    if np.any(y < 0) or np.any(y != np.trunc(y)):
        raise StataError(459, f"{pr.smp.depname} must be nonnegative integers")   # VERIFY
    out = s.output
    pr.write_notes(out)
    ds = s.data
    off = np.zeros(len(y))
    if pr.o.get("exposure"):
        off = np.log(ds.get(str(pr.o["exposure"]).strip()).data[pr.smp.mask].astype(float))
    elif pr.o.get("offset"):
        off = ds.get(str(pr.o["offset"]).strip()).data[pr.smp.mask].astype(float)
    eform = "IRR" if pr.o.get("irr") else ""
    X = pr.X
    from .ml import maximize as _max
    from .models import _poisson
    if P.fe:
        tot = np.bincount(inv, weights=y)
        keep = tot[inv] > 0
        drop_g = int((tot == 0).sum())
        if drop_g:
            out.write(f"note: {drop_g} group{'s' if drop_g > 1 else ''} ({int((~keep).sum())} obs) dropped because "
                      "of all zero outcomes\n", "text")
        Xk, yk, ok = X[keep], y[keep], off[keep]
        _, invk = np.unique(inv[keep], return_inverse=True)
        grp = pq.Groups(invk)

        def llg(b):
            return pq.nbreg_fe_group_ll(np.exp(np.clip(Xk @ b + ok, -700, 700)), yk, grp)

        def ll(b):
            return float(llg(b).sum())
        pf, _ = _poisson(yk, np.column_stack([Xk, np.ones(len(yk))]), np.ones(len(yk)), ok)
        rp = _max(pf, np.zeros(Xk.shape[1] + 1))
        res = pr.maximize(_numeric_fun(ll), rp.b[:-1])
        if not pr.o.get("nolog"):
            out.write("\n", "text")
            print_log(out, res, "log likelihood")
        V, robust = _vce(P, llg, res, grp)
        _finish(s, P, "xtnbreg", args, res.b, V, res, [], "Conditional FE negative binomial regression", "", keep,
                eform=eform, robust=robust)
        return
    Xc = np.column_stack([X, np.ones(len(y))]) if pr.constant else X
    k = Xc.shape[1]
    grp = pq.Groups(inv)

    def llg(par):
        mu = np.exp(np.clip(Xc @ par[:k] + off, -700, 700))
        return pq.nbreg_re_group_ll(mu, y, grp, par[k], par[k + 1])

    def ll(par):
        return float(llg(par).sum())
    lg = special.gammaln
    lgy = lg(y + 1)

    def ll_c(b):
        # nbreg, dispersion(constant): variance mu (1 + delta); delta = exp(b[k])
        mu = np.exp(np.clip(Xc @ b[:k] + off, -700, 700))
        d = math.exp(b[k])
        lam = mu / d
        return float((lg(lam + y) - lg(lam) - lgy - lam * math.log1p(d) + y * (b[k] - math.log1p(d))).sum())
    pf, _ = _poisson(y, Xc, np.ones(len(y)), off)
    b0 = np.zeros(k)
    if pr.constant:
        b0[-1] = math.log(max(y.mean(), 1e-8))
    rp = _max(pf, b0)
    rc = _max(_numeric_fun(ll_c), np.concatenate([rp.b, [0.0]]))
    # the constant of the dispersion-constant model is ln(mu/lambda): the effect enters the other way
    bstart = rc.b[:k].copy()
    if pr.constant:
        bstart[-1] -= rc.b[k]
    best = None
    for lr_ in range(-1, 6):
        for ls_ in range(-1, 6):
            v = ll(np.concatenate([bstart, [lr_, ls_]]))
            if np.isfinite(v) and (best is None or v > best[0]):
                best = (v, lr_, ls_)
    res = pr.maximize(_numeric_fun(ll), np.concatenate([bstart, [best[1], best[2]]]))
    if not pr.o.get("nolog"):
        out.write("\nFitting negative binomial (constant dispersion) model:\n\n", "text")
        print_log(out, rc, "log likelihood")
        out.write("\nFitting full model:\n\n", "text")
        print_log(out, res, "log likelihood")
    V, robust = _vce(P, llg, res, grp)
    lr_, ls_ = float(res.b[k]), float(res.b[k + 1])
    se_r, se_s = float(np.sqrt(V[k, k])), float(np.sqrt(V[k + 1, k + 1]))
    z = crit("z", pr.smp.level, None)
    r_, s_ = math.exp(lr_), math.exp(ls_)
    aux = [CoefRow("sep"), CoefRow("aux", "/ln_r", lr_, se_r, name="_cons"),
           CoefRow("aux", "/ln_s", ls_, se_s, name="_cons"), CoefRow("sep"),
           CoefRow("aux", "r", r_, r_ * se_r, eq=f"{_exp(lr_ - z * se_r)} {_exp(lr_ + z * se_r)}"),
           CoefRow("aux", "s", s_, s_ * se_s, eq=f"{_exp(ls_ - z * se_s)} {_exp(ls_ + z * se_s)}")]
    chib = max(2 * (res.ll - rc.ll), 0.0)
    pchib = 0.5 * chi2_tail(chib, 1) if chib > 0 else 0.5
    _finish(s, P, "xtnbreg", args, res.b, V, res, ["/ln_r", "/ln_s"], "Random-effects negative binomial regression",
            "Random effects u_i ~ Beta", np.ones(len(y), bool), chibar=("LR test vs. pooled", chib, pchib),
            aux_rows=aux, eform=eform, robust=robust,
            extra_scalars=[("r", r_), ("s", s_), ("ll_c", rc.ll), ("chi2_c", chib), ("p_c", pchib)])


# ---------------------------------------------------------------------------
# quadchk
# ---------------------------------------------------------------------------

def cmd_quadchk(s: "Session", args: str) -> None:
    from .postest import _restore, _snapshot
    info = getattr(s, "_quadchk", None)
    est = getattr(s, "_estimates", None)
    if info is None or est is None or est.cmd != info["cmd"]:
        raise StataError(301, "last estimates not found")
    head, _, opts = args.partition(",")
    nums = head.split()
    nq = info["nq"]
    if nums:
        try:
            n1, n2 = (int(v) for v in nums[:2]) if len(nums) >= 2 else (int(nums[0]), None)
        except ValueError:
            raise StataError(198, "invalid syntax")
        if n2 is None:
            raise StataError(198, "invalid syntax")
    else:
        n1, n2 = max(int(round(2 * nq / 3)), 1), int(round(4 * nq / 3))
    show = "nooutput" not in opts.replace(" ", "").lower()
    snap = _snapshot(s)
    saved_info = dict(info)
    base_args = _drop_intpoints(info["args"])
    runs = [(nq, info["ll"], est.b.copy(), list(est.names), list(est.eqnames))]
    from ..commands.registry import lookup
    cmd_fn = lookup(info["cmd"]).fn
    out = s.output
    try:
        for n in (n1, n2):
            if show:
                out.write(f"\nRefitting model intpoints() = {n}\n", "text")
            run_args = f"{base_args}, intpoints({n})" if "," not in base_args else f"{base_args} intpoints({n})"
            if show:
                cmd_fn(s, run_args)
            else:
                with out.quietly():
                    cmd_fn(s, run_args)
            e2 = s._estimates
            runs.append((n, float(s.e["ll"]), e2.b.copy(), list(e2.names), list(e2.eqnames)))
    finally:
        _restore(s, snap)
        s._quadchk = saved_info
    _quad_table(s, runs, est)


def _drop_intpoints(args: str) -> str:
    import re
    return re.sub(r"\bintpoints\s*\(\s*\d+\s*\)", "", args).rstrip()


def _quad_table(s, runs, est) -> None:
    from ..core.formats import format_value
    out = s.output

    def f(v):
        return format_value(float(v), "%11.0g", pad=False).strip()

    def e3(v):
        return format_value(float(v), "%10.4e", pad=False).strip()
    (n0, ll0, b0, names, eqs), (n1, ll1, b1, _, _), (n2, ll2, b2, _, _) = runs
    out.write("\n", "text")
    out.write(f"{'':>25}Quadrature check\n", "text")
    out.write(f"{'Fitted':>26}{'Comparison':>16}{'Comparison':>16}\n", "text")
    out.write(f"{'quadrature':>26}{'quadrature':>16}{'quadrature':>16}\n", "text")
    out.write(f"{str(n0) + ' points':>26}{str(n1) + ' points':>16}{str(n2) + ' points':>16}\n", "text")
    out.write("-" * 79 + "\n", "text")
    out.write(f"{'Log':<14}{f(ll0):>12}{f(ll1):>16}{f(ll2):>16}\n", "text")
    out.write(f"{'likelihood':<14}{'':>12}{e3(ll1 - ll0):>16}{e3(ll2 - ll0):>16}  Difference\n", "text")
    out.write(f"{'':<14}{'':>12}{e3((ll1 - ll0) / abs(ll0)):>16}{e3((ll2 - ll0) / abs(ll0)):>16}  Relative difference\n",
              "text")
    last_eq = None
    for i, nm in enumerate(names):
        eq = eqs[i] if i < len(eqs) else ""
        if eq != last_eq:
            out.write("-" * 79 + "\n", "text")
            out.write(f"{eq + ':':<14}\n", "text")
            last_eq = eq
        lab = {"/lnsig2u": "lnsig2u"}.get(nm, nm)
        out.write(f"{lab:<14}{f(b0[i]):>12}{f(b1[i]):>16}{f(b2[i]):>16}\n", "text")
        d1, d2 = b1[i] - b0[i], b2[i] - b0[i]
        out.write(f"{'':<14}{'':>12}{e3(d1):>16}{e3(d2):>16}  Difference\n", "text")
        rel1 = d1 / abs(b0[i]) if b0[i] else float("nan")
        rel2 = d2 / abs(b0[i]) if b0[i] else float("nan")
        out.write(f"{'':<14}{'':>12}{e3(rel1):>16}{e3(rel2):>16}  Relative difference\n", "text")
    out.write("-" * 79 + "\n", "text")
    out.write("\n", "text")
