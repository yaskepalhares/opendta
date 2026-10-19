"""Plain functions for survival-time data manipulation ([ST] stsplit, stjoin, stfill, stgen, stvary,
snapspan): they take arrays and return arrays, without a Session, variables or printing."""

from __future__ import annotations

import numpy as np


def factorize(keys: np.ndarray) -> tuple[np.ndarray, int]:
    """Subject number (0, 1, ...) of each element of `keys` (numbers or strings), by sorted value."""
    _u, inv = np.unique(keys, return_inverse=True)
    return inv.astype(np.int64), len(_u)


def time_order(gid: np.ndarray, t0: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Positions that sort the records by subject, then entry time, then exit time (stable)."""
    return np.lexsort((t, t0, gid))


def group_bounds(gid_sorted: np.ndarray) -> np.ndarray:
    """Start position of each run of equal subject numbers in a sorted vector (plus the end)."""
    if len(gid_sorted) == 0:
        return np.zeros(1, dtype=np.int64)
    cut = np.flatnonzero(np.diff(gid_sorted)) + 1
    return np.concatenate([[0], cut, [len(gid_sorted)]]).astype(np.int64)


# ---------------------------------------------------------------------------
# stsplit
# ---------------------------------------------------------------------------

def split_pieces(t0: np.ndarray, t: np.ndarray, cuts: list[np.ndarray]):
    """Splits each record (t0, t] at the cut points of that record that fall strictly inside it.

    `cuts[i]` is the sorted vector of cut points of record i. Returns (parent, a, b, last): the
    record each piece comes from, its limits (a, b] and whether it is the final piece of the record
    (the one that keeps the failure indicator)."""
    parent: list[int] = []
    lo: list[float] = []
    hi: list[float] = []
    last: list[bool] = []
    for i in range(len(t)):
        c = cuts[i]
        inside = c[(c > t0[i]) & (c < t[i])] if len(c) else c
        edges = np.concatenate([[t0[i]], inside, [t[i]]])
        k = len(edges) - 1
        parent.extend([i] * k)
        lo.extend(edges[:-1].tolist())
        hi.extend(edges[1:].tolist())
        last.extend([False] * (k - 1) + [True])
    return (np.array(parent, dtype=np.int64), np.array(lo), np.array(hi), np.array(last, dtype=bool))


def category_labels(b: np.ndarray, cuts: list[np.ndarray], parent: np.ndarray, at: np.ndarray
                    ) -> np.ndarray:
    """The lower bound of the interval of `at` that contains each piece: pieces end at or before
    the next cut point, so the number of cut points strictly below the end of the piece picks the
    interval. The first interval is labelled 0, or one less than the smallest at() value when some
    value is negative ([ST] stsplit)."""
    first = 0.0 if not len(at) or at.min() >= 0 else float(at.min()) - 1
    labels = np.concatenate([[first], at])
    out = np.empty(len(b))
    for j in range(len(b)):
        c = cuts[parent[j]]
        out[j] = labels[int(np.sum(c < b[j]))] if len(c) else labels[0]
    return out


# ---------------------------------------------------------------------------
# stjoin
# ---------------------------------------------------------------------------

def join_chains(order: np.ndarray, gid: np.ndarray, t0: np.ndarray, t: np.ndarray,
                censored: np.ndarray, same: np.ndarray) -> np.ndarray:
    """For records in `order` (sorted by subject and time), True where the record is joined to the
    one before it: same subject, it starts where the previous one ended, the previous one ended
    without an event and `same` (record i has the same values as the one before it in `order`)."""
    n = len(order)
    join = np.zeros(n, dtype=bool)
    if n < 2:
        return join
    g, a, b, c = gid[order], t0[order], t[order], censored[order]
    join[1:] = (g[1:] == g[:-1]) & (a[1:] == b[:-1]) & c[:-1] & same
    return join


# ---------------------------------------------------------------------------
# stfill, stvary, stgen
# ---------------------------------------------------------------------------

def is_missing_values(x: np.ndarray, is_string: bool, sysmiss: float) -> np.ndarray:
    if is_string:
        return np.array([v == "" for v in x], dtype=bool)
    return x >= sysmiss


def fill_forward(values: np.ndarray, miss: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    """Positions of the source value of each element after carrying the last nonmissing value
    forward within each subject (-1 where nothing has been observed yet)."""
    n = len(values)
    src = np.full(n, -1, dtype=np.int64)
    for a, b in zip(bounds[:-1], bounds[1:]):
        cur = -1
        for k in range(a, b):
            if not miss[k]:
                cur = k
            src[k] = cur
    return src


def subject_summary(values: np.ndarray, miss: np.ndarray, bounds: np.ndarray):
    """Per subject: (constant, varying, never missing, always missing, sometimes missing) as in
    stvary. A subject whose values are all missing is neither constant nor varying."""
    cons = vary = never = always = some = 0
    for a, b in zip(bounds[:-1], bounds[1:]):
        m = miss[a:b]
        ok = values[a:b][~m]
        if m.all():
            always += 1
        elif m.any():
            some += 1
        else:
            never += 1
        if len(ok):
            if len(set(ok.tolist())) == 1:
                cons += 1
            else:
                vary += 1
    return cons, vary, never, always, some


def cumulative_before(flag: np.ndarray, bounds: np.ndarray, include_current: bool) -> np.ndarray:
    """Number of records with flag set up to the previous record (or up to this one) of the same
    subject: count() and count0() of stgen."""
    out = np.zeros(len(flag))
    for a, b in zip(bounds[:-1], bounds[1:]):
        c = np.cumsum(flag[a:b].astype(float))
        out[a:b] = c if include_current else c - flag[a:b]
    return out


def gaps(t0: np.ndarray, t: np.ndarray, bounds: np.ndarray):
    """For records sorted by subject and time: the number of gaps before each record and their
    total length. Delayed entry (the first record) is not a gap."""
    n = len(t)
    ng = np.zeros(n)
    gl = np.zeros(n)
    for a, b in zip(bounds[:-1], bounds[1:]):
        cnt, tot = 0.0, 0.0
        for k in range(a, b):
            if k > a and t0[k] > t[k - 1]:
                cnt += 1
                tot += t0[k] - t[k - 1]
            ng[k], gl[k] = cnt, tot
    return ng, gl


# ---------------------------------------------------------------------------
# rates: strate, stptime
# ---------------------------------------------------------------------------

def _z(level: float) -> float:
    from scipy import stats
    return float(stats.norm.ppf(1 - (1 - level / 100) / 2))


def poisson_rate_ci(d: float, rate: float, level: float) -> tuple[float, float]:
    """Limits from a quadratic approximation to the Poisson log likelihood of ln(rate): the
    standard error of ln(rate) is 1/sqrt(failures) ([ST] strate). Missing without failures."""
    if not d > 0:
        return float("nan"), float("nan")
    f = float(np.exp(_z(level) / np.sqrt(d)))
    return rate / f, rate * f


def jackknife_log_ci(num: np.ndarray, den: np.ndarray, cluster: np.ndarray, level: float
                     ) -> tuple[float, float, float]:
    """Rate = sum(num)/sum(den) with a delete-one-cluster jackknife on ln(rate): returns
    (rate, lower, upper). VERIFY: Stata's jackknife intervals."""
    ids, inv = np.unique(cluster, return_inverse=True)
    g = len(ids)
    N, D = float(num.sum()), float(den.sum())
    rate = N / D
    if g < 2 or not rate > 0:
        return rate, float("nan"), float("nan")
    cn = np.bincount(inv, weights=num, minlength=g)
    cd = np.bincount(inv, weights=den, minlength=g)
    with np.errstate(divide="ignore", invalid="ignore"):
        loo = np.log((N - cn) / (D - cd))
    theta = np.log(rate)
    pseudo = g * theta - (g - 1) * loo
    if not np.all(np.isfinite(pseudo)):
        return rate, float("nan"), float("nan")
    se = float(np.sqrt(np.var(pseudo, ddof=1) / g))
    z = _z(level)
    return rate, float(np.exp(theta - z * se)), float(np.exp(theta + z * se))


# ---------------------------------------------------------------------------
# stci
# ---------------------------------------------------------------------------

def km_percentile_se(times: np.ndarray, S: np.ndarray, gw: np.ndarray, tp: float, p: float,
                     eps: float = 0.05) -> float:
    """Large-sample standard error of the p-th percentile, sqrt(Var S(tp)) / f(tp), with the
    density estimated as in Collett (2003): (S(u) - S(l)) / (l - u) with u the largest time with
    S >= 1 - p + eps and l the smallest time with S <= 1 - p - eps (VERIFY: Stata's epsilon)."""
    q = 1 - p
    up = np.flatnonzero(S >= q + eps)
    lo = np.flatnonzero(S <= q - eps)
    if not len(lo):
        return float("nan")
    u = float(times[up[-1]]) if len(up) else 0.0
    Su = float(S[up[-1]]) if len(up) else 1.0
    l_ = float(times[lo[0]])
    dens = (Su - float(S[lo[0]])) / (l_ - u) if l_ > u else float("nan")
    j = np.searchsorted(times, tp, "right") - 1
    var = (S[j] ** 2) * gw[j] if j >= 0 else 0.0
    return float(np.sqrt(var) / dens) if dens and np.isfinite(dens) else float("nan")


def km_rmean(times: np.ndarray, S: np.ndarray, d: np.ndarray, n: np.ndarray, tmax: float,
             ccorr: bool = False):
    """Area under the Kaplan-Meier curve up to tmax and its standard error ([ST] stci). Returns
    (mean, se). The continuity correction multiplies the variance by m/(m-1) (VERIFY)."""
    keep = times <= tmax
    t, Sx = times[keep], S[keep]
    edges = np.concatenate([[0.0], t, [tmax]]) if (len(t) == 0 or t[-1] < tmax) else np.concatenate([[0.0], t])
    vals = np.concatenate([[1.0], Sx])
    if len(edges) > len(vals):
        vals = np.concatenate([vals, [vals[-1]]])
    seg = np.diff(edges)
    area_after = np.concatenate([np.cumsum((seg * vals[:len(seg)])[::-1])[::-1], [0.0]])
    mean = float(area_after[0])
    # A_i: area from the i-th failure time to tmax
    idx = np.searchsorted(edges, t) if len(t) else np.array([], dtype=int)
    A = area_after[np.minimum(idx, len(area_after) - 1)] if len(t) else np.array([])
    dk, nk = d[keep], n[keep]
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where((dk > 0) & (nk > dk), dk / (nk * (nk - dk)), 0.0)
    var = float(np.sum(A * A * term))
    if ccorr:
        m = float(dk.sum())
        var *= m / (m - 1) if m > 1 else 1.0
    return mean, float(np.sqrt(var))


def km_emean(times: np.ndarray, S: np.ndarray, rmean: float, tmax: float) -> float:
    """Extended mean: the restricted mean plus the area under an exponential tail that starts at
    the last time with the last survivor value (VERIFY: Stata's extension)."""
    if not len(S) or S[-1] <= 0:
        return rmean
    return rmean + float(S[-1]) * tmax / -float(np.log(S[-1]))


# ---------------------------------------------------------------------------
# ltable
# ---------------------------------------------------------------------------


def life_table(t: np.ndarray, dead: np.ndarray, w: np.ndarray, spec, adjust: bool, level: float):
    """The rows of [ST] ltable: for each interval with observations, a dict with lo, hi (None for the
    open last interval), N (beginning total), d, m, n (adjusted at risk), S, se, lb, ub, hazard,
    hazard se and interval. `spec` is a width or the list of cut points."""
    if np.isscalar(spec):
        wd = float(spec)
        lo_all = np.floor(t / wd) * wd
        starts = np.unique(lo_all)
        his = {s_: s_ + wd for s_ in starts}
        key = lo_all
    else:
        cuts = np.asarray(sorted(spec), dtype=float)
        idx = np.searchsorted(cuts, t, "right") - 1
        key = np.where(idx >= 0, cuts[np.maximum(idx, 0)], np.nan)
        starts = np.unique(key[~np.isnan(key)])
        his = {c: (cuts[i + 1] if i + 1 < len(cuts) else None) for i, c in enumerate(cuts)}
    z = _z(level)
    rows = []
    alive = float(w[~np.isnan(key)].sum())
    S, gsum = 1.0, 0.0
    zero = False
    for c in starts:
        sel = key == c
        d = float(w[sel & (dead != 0)].sum())
        m = float(w[sel & (dead == 0)].sum())
        N = alive
        n = N - m / 2 if adjust else N
        row = {"lo": float(c), "hi": his[c], "N": N, "d": d, "m": m, "n": n}
        if not zero:
            q = 1 - d / n if n > 0 else 1.0
            if n > 0 and n - d > 0:
                gsum += d / (n * (n - d))
            S *= q
            if n > 0 and n - d <= 0:
                zero = True
        row["S"] = S
        if zero or not (S > 0):
            row.update(se=float("nan"), lb=float("nan"), ub=float("nan"))
        elif S >= 1:
            row.update(se=0.0, lb=float("nan"), ub=float("nan"))
        else:
            se = S * float(np.sqrt(gsum))
            sig = float(np.sqrt(gsum)) / abs(float(np.log(S)))
            row.update(se=se, lb=S ** float(np.exp(z * sig)), ub=S ** float(np.exp(-z * sig)))
        # hazard
        if row["hi"] is None:
            row.update(h=float("nan"), hse=float("nan"), hlb=float("nan"), hub=float("nan"))
        else:
            width = row["hi"] - row["lo"]
            if adjust:
                f = d / n if n > 0 else float("nan")
                h = f / ((1 - f / 2) * width)
                hse = h * float(np.sqrt(max(1 - (width * h / 2) ** 2, 0)) / np.sqrt(d)) if d > 0 else float("nan")
                hlb, hub = h - z * hse, h + z * hse
            else:
                f = d / N if N > 0 else float("nan")
                h = f / width
                hse = h / float(np.sqrt(d)) if d > 0 else float("nan")
                from scipy import stats
                if d > 0:
                    hlb = h / (2 * d) * float(stats.chi2.ppf((1 - level / 100) / 2, 2 * d))
                    hub = h / (2 * d) * float(stats.chi2.ppf(1 - (1 - level / 100) / 2, 2 * d))
                else:
                    hlb = hub = float("nan")
            row.update(h=h, hse=hse, hlb=hlb, hub=hub)
        rows.append(row)
        alive -= d + m
    return rows


def lr_exponential(d: np.ndarray, T: np.ndarray) -> tuple[float, int]:
    """Likelihood-ratio test of equal exponential hazards: 2[sum d_g ln(d_g/T_g) - D ln(D/T)]."""
    d, T = np.asarray(d, float), np.asarray(T, float)
    ok = d > 0
    D, TT = d.sum(), T.sum()
    stat = 2 * (float(np.sum(d[ok] * np.log(d[ok] / T[ok]))) - D * float(np.log(D / TT)))
    return stat, len(d) - 1


def logrank_chi2(t: np.ndarray, dead: np.ndarray, w: np.ndarray, g: np.ndarray, G: int) -> tuple[float, int]:
    """Log-rank statistic for G groups (g = 0..G-1) with the hypergeometric variance."""
    O = np.zeros(G)
    E = np.zeros(G)
    V = np.zeros((G, G))
    for x in np.unique(t[dead != 0]):
        risk = t >= x
        nt = float(w[risk].sum())
        ev = (t == x) & (dead != 0)
        dt = float(w[ev].sum())
        ng = np.bincount(g[risk], weights=w[risk], minlength=G)
        dg = np.bincount(g[ev], weights=w[ev], minlength=G)
        O += dg
        E += dt * ng / nt
        if nt > 1:
            V += dt * (nt - dt) / (nt - 1) / nt * (np.diag(ng) - np.outer(ng, ng) / nt)
    k = G - 1
    diff = (O - E)[:k]
    try:
        chi2 = float(diff @ np.linalg.solve(V[:k, :k], diff))
    except np.linalg.LinAlgError:
        chi2 = float(diff @ np.linalg.pinv(V[:k, :k]) @ diff)
    return chi2, k


# ---------------------------------------------------------------------------
# exact ties of stcox
# ---------------------------------------------------------------------------

_GL = None


def _gauss_laguerre():
    global _GL
    if _GL is None:
        _GL = np.polynomial.laguerre.laggauss(15)
    return _GL


def exactm_term(xd: np.ndarray, xn: np.ndarray, b: np.ndarray):
    """Log of the exact marginal-likelihood term of one failure time and its gradient and Hessian.
    xd: covariates of the d tied failures; xn: those of the other subjects at risk. The term is
    the integral over (0, inf) of exp(-t) times the product over the failures of
    1 - exp(-t e_l / s), with s the sum of exp(xb) over the others, by 15-point Gauss-Laguerre
    quadrature ([ST] stcox, Methods and formulas)."""
    k = len(b)
    rd, rn = np.exp(xd @ b), np.exp(xn @ b)
    A = float(rn.sum())
    if A <= 0:
        return 0.0, np.zeros(k), np.zeros((k, k))
    z, w = _gauss_laguerre()
    abar = (rn[:, None] * xn).sum(axis=0) / A
    A2 = np.einsum("j,jk,jl->kl", rn, xn, xn) / A
    xc = xd - abar
    y = (z[:, None] / A) * rd[None, :]            # K x d
    E = np.exp(-y)
    F = np.prod(1 - E, axis=1)
    L = float(np.sum(w * F))
    gm = y * E / (1 - E)
    gp = E / (1 - E) - y * E / (1 - E) ** 2
    v = np.einsum("ki,ic->kc", gm, xc)
    pi = w * F / L
    vbar = pi @ v
    quad = np.einsum("k,kc,kl->cl", pi, v, v)
    quad += np.einsum("k,ki,ic,il->cl", pi, gp * y, xc, xc)
    quad -= np.einsum("k,k->", pi, gm.sum(axis=1)) * (A2 - np.outer(abar, abar))
    return float(np.log(L)), vbar, quad - np.outer(vbar, vbar)


def exactp_term(xr: np.ndarray, is_fail: np.ndarray, b: np.ndarray):
    """Log of the exact partial-likelihood (conditional logistic) term of one failure time and its
    gradient and Hessian: sum of xb over the failures minus log f(r, d), with the recursion
    f(r, d) = f(r-1, d) + f(r-1, d-1) exp(xb_r)."""
    k = len(b)
    eta = xr @ b
    c = float(eta.max())
    e = np.exp(eta - c)
    d = int(is_fail.sum())
    f = np.zeros(d + 1)
    f[0] = 1.0
    g = np.zeros((d + 1, k))
    H = np.zeros((d + 1, k, k))
    for i in range(len(e)):
        x, ei = xr[i], e[i]
        xx = np.outer(x, x)
        for m in range(min(d, i + 1), 0, -1):
            H[m] += ei * (xx * f[m - 1] + np.outer(x, g[m - 1]) + np.outer(g[m - 1], x) + H[m - 1])
            g[m] += ei * (x * f[m - 1] + g[m - 1])
            f[m] += ei * f[m - 1]
    ll = float(eta[is_fail].sum()) - d * c - float(np.log(f[d]))
    grad = xr[is_fail].sum(axis=0) - g[d] / f[d]
    hess = -(H[d] / f[d] - np.outer(g[d], g[d]) / f[d] ** 2)
    return ll, grad, hess
