"""Treatment-effects estimators on arrays ([TE] teffects, eteffects): regression adjustment,
inverse-probability weighting, their combinations, the control-function estimator and nearest-
neighbour matching. Binary treatment. Each estimator is a set of stacked moment conditions whose
solution is found in closed form from sub-fits; the variance is the robust M-estimation sandwich
with a numerical Jacobian. No Session and no output."""

from __future__ import annotations

import numpy as np
from scipy import stats


# ---------------------------------------------------------------------------
# generalized linear fits: mean, derivative of the mean and variance functions
# ---------------------------------------------------------------------------

def _mean(kind, eta):
    if kind == "linear":
        return eta, np.ones_like(eta), np.ones_like(eta)
    if kind == "logit":
        m = 1 / (1 + np.exp(-eta))
        return m, m * (1 - m), m * (1 - m)
    if kind == "probit":
        m = np.clip(stats.norm.cdf(eta), 1e-15, 1 - 1e-15)
        return m, stats.norm.pdf(eta), m * (1 - m)
    if kind == "poisson":
        m = np.exp(np.clip(eta, -50, 50))
        return m, m, m
    raise ValueError(kind)


def glm_mean(kind, X, b):
    return _mean(kind, X @ b)[0]


def glm_score(kind, X, y, b):
    """Score contributions (n x k) of a GLM with the canonical-or-probit link."""
    m, dm, v = _mean(kind, X @ b)
    return X * ((y - m) * dm / v)[:, None]


def glm_fit(kind, X, y, w=None, tol=1e-12, maxit=100):
    """Weighted maximum likelihood by Fisher scoring."""
    n, k = X.shape
    w = np.ones(n) if w is None else w
    b = np.zeros(k)
    if kind == "linear":
        sw = np.sqrt(w)
        return np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    if kind == "poisson":
        b[-1] = np.log(max(np.average(y, weights=w), 1e-8)) if np.all(X[:, -1] == 1) else 0.0
    for _ in range(maxit):
        m, dm, v = _mean(kind, X @ b)
        sc = X.T @ (w * (y - m) * dm / v)
        info = (X * (w * dm * dm / v)[:, None]).T @ X
        step = np.linalg.solve(info + 1e-12 * np.eye(k), sc)
        b = b + step
        if np.max(np.abs(step)) < tol:
            break
    return b


class GLMOutcome:
    """Outcome model of a GLM: fit with case weights, mean and score of the linear index."""

    def __init__(self, kind, y):
        self.kind, self.y = kind, y

    def fit(self, X, w):
        m = w > 0
        return glm_fit(self.kind, X[m], self.y[m], w[m])

    def mean(self, X, b):
        return glm_mean(self.kind, X, b)

    def score(self, X, b):
        return glm_score(self.kind, X, self.y, b)


class LikOutcome:
    """Outcome model given by a log likelihood lf(u, th) -> (l, l_u, l_th, ...) per observation
    (a parametric survival model) with an optional ancillary parameter th, and the mean of the
    outcome mean_fn(u, th). The parameters are (beta, th)."""

    def __init__(self, lf, mean_fn, kaux, start):
        self.lf, self.mean_fn, self.kaux, self.start = lf, mean_fn, kaux, start

    def _split(self, X, b):
        k = X.shape[1]
        return X @ b[:k], (b[k] if self.kaux else 0.0)

    def fit(self, X, w):
        from scipy import optimize
        k = X.shape[1]
        sw = w.sum()

        def fun(b):
            u, th = self._split(X, b)
            with np.errstate(all="ignore"):
                l, lu, lt = self.lf(u, th)[:3]
            val = -np.sum(w * l) / sw
            if not np.isfinite(val):
                return 1e300, np.zeros_like(b)
            g = np.concatenate([X.T @ (w * lu), [np.sum(w * lt)] if self.kaux else []])
            return val, -g / sw

        b0 = np.zeros(k + self.kaux)
        b0[k - 1] = self.start[0]
        if self.kaux and len(self.start) > 1:
            b0[k] = self.start[1]
        r = optimize.minimize(fun, b0, jac=True, method="BFGS", options={"gtol": 1e-11, "maxiter": 500})
        return r.x

    def mean(self, X, b):
        u, th = self._split(X, b)
        return self.mean_fn(u, th)

    def score(self, X, b):
        u, th = self._split(X, b)
        _l, lu, lt = self.lf(u, th)[:3]
        return np.column_stack([X * lu[:, None]] + ([lt[:, None]] if self.kaux else []))


# ---------------------------------------------------------------------------
# M-estimation variance
# ---------------------------------------------------------------------------

def sandwich_m(gfun, theta, n, clusters=None):
    """V = A^-1 B A^-T / n with A the numerical Jacobian of the mean moments and B their outer
    product; `gfun(theta)` returns the (n x p) per-observation moments."""
    p = len(theta)
    A = np.zeros((p, p))
    for j in range(p):
        h = 1e-6 * max(1.0, abs(theta[j]))
        up, dn = theta.copy(), theta.copy()
        up[j] += h
        dn[j] -= h
        A[:, j] = (gfun(up).mean(axis=0) - gfun(dn).mean(axis=0)) / (2 * h)
    g = gfun(theta)
    if clusters is None:
        B = g.T @ g / n
    else:
        _, inv = np.unique(clusters, return_inverse=True)
        gc = np.zeros((inv.max() + 1, p))
        np.add.at(gc, inv, g)
        B = gc.T @ gc / n
    Ai = np.linalg.inv(A)
    V = Ai @ B @ Ai.T / n
    return (V + V.T) / 2


# ---------------------------------------------------------------------------
# regression adjustment, IPW, IPWRA, AIPW
# ---------------------------------------------------------------------------

class Fit:
    """Result of one estimator: theta = [mu0, mu1, beta0, beta1, gamma] (blocks that the estimator
    does not use are empty) and its robust variance."""

    def __init__(self, theta, V, mu_names, n):
        self.theta, self.V, self.n = theta, V, n


def teffects(est, y, t, Xo, Z, omodel, tmodel, stat, clusters=None):
    """est: ra, ipw, ipwra, aipw. t is 0/1; Xo (outcome design, with constant), Z (treatment design,
    with constant); stat: ate, atet or pomeans. Returns (theta, V, layout) where layout gives the
    slices of mu (2), beta0, beta1, gamma."""
    n, ko = Xo.shape
    om = GLMOutcome(omodel, y) if isinstance(omodel, str) else omodel
    kz = Z.shape[1] if est != "ra" else 0
    d = np.ones(n) if stat != "atet" else t.copy()
    use_out = est != "ipw"
    use_trt = est != "ra"
    sl_mu = slice(0, 2)
    kb = ko + getattr(om, "kaux", 0)
    sl_b0 = slice(2, 2 + kb) if use_out else slice(2, 2)
    sl_b1 = slice(sl_b0.stop, sl_b0.stop + kb) if use_out else slice(2, 2)
    sl_g = slice(sl_b1.stop, sl_b1.stop + kz)
    t0 = 1 - t

    def prob(g):
        return _mean(tmodel, Z @ g)[0]

    def weights(p):
        if stat == "atet":
            return t / 1.0, t0 * p / (1 - p)
        return t / p, t0 / (1 - p)

    # ---- point estimates
    g = glm_fit(tmodel, Z, t) if use_trt else np.zeros(0)
    p = prob(g) if use_trt else None
    if use_trt:
        w1, w0 = weights(p)
    b0 = b1 = np.zeros(0)
    if use_out:
        if est in ("ra", "aipw"):               # the aipw outcome models are fitted arm by arm, unweighted
            b0, b1 = om.fit(Xo, t0), om.fit(Xo, t)
        else:
            b0, b1 = om.fit(Xo, w0), om.fit(Xo, w1)
    if est == "ipw":
        mu = np.array([np.sum(w0 * y) / np.sum(w0), np.sum(w1 * y) / np.sum(w1)])
    else:
        m0, m1 = om.mean(Xo, b0), om.mean(Xo, b1)
        if est == "aipw":
            a0 = m0 + w0 * (y - m0)
            a1 = m1 + w1 * (y - m1)
            mu = np.array([a0.mean(), a1.mean()])
        else:
            mu = np.array([np.sum(d * m0) / d.sum(), np.sum(d * m1) / d.sum()])
    theta = np.concatenate([mu, b0, b1, g])

    # ---- moments
    def gfun(th):
        mu_, b0_, b1_, g_ = th[sl_mu], th[sl_b0], th[sl_b1], th[sl_g]
        cols = []
        if est == "ra":
            cols.append(t0[:, None] * om.score(Xo, b0_))
            cols.append(t[:, None] * om.score(Xo, b1_))
            m0_, m1_ = om.mean(Xo, b0_), om.mean(Xo, b1_)
            mom = np.column_stack([d * (m0_ - mu_[0]), d * (m1_ - mu_[1])])
        else:
            pp = prob(g_)
            cols.append(glm_score(tmodel, Z, t, g_))
            ww1, ww0 = weights(pp)
            if est == "ipw":
                mom = np.column_stack([ww0 * (y - mu_[0]), ww1 * (y - mu_[1])])
            else:
                a0_, a1_ = (t0, t) if est == "aipw" else (ww0, ww1)
                cols.append(a0_[:, None] * om.score(Xo, b0_))
                cols.append(a1_[:, None] * om.score(Xo, b1_))
                m0_, m1_ = om.mean(Xo, b0_), om.mean(Xo, b1_)
                if est == "ipwra":
                    mom = np.column_stack([d * (m0_ - mu_[0]), d * (m1_ - mu_[1])])
                else:
                    mom = np.column_stack([m0_ + ww0 * (y - m0_) - mu_[0], m1_ + ww1 * (y - m1_) - mu_[1]])
        # order: mu moments first, then the sub-model scores in the order of theta
        if est == "ra":
            return np.column_stack([mom, cols[0], cols[1]])
        if est == "ipw":
            return np.column_stack([mom, cols[0]])
        return np.column_stack([mom, cols[1], cols[2], cols[0]])

    V = sandwich_m(gfun, theta, n, clusters)
    return theta, V, (sl_mu, sl_b0, sl_b1, sl_g)


# ---------------------------------------------------------------------------
# control function (eteffects)
# ---------------------------------------------------------------------------

def eteffects(y, t, Xo, Z, stat, clusters=None, omodel="linear"):
    """Control-function estimator of [TE] eteffects (Methods and formulas): a probit treatment
    model, the residual nu = t - Phi(z'pi), and in each arm an outcome model in x and nu (linear,
    probit, fractional probit or exponential mean), all fitted together with the potential-outcome
    means as an exactly identified GMM problem. POM_j is the mean over the sample (the treated with
    atet) of m(x'b_1j + nu b_2j). theta = [mu0, mu1, beta0, rho0, beta1, rho1, pi]."""
    n, ko = Xo.shape
    d = np.ones(n) if stat != "atet" else t.copy()
    t0 = 1 - t
    kind = {"linear": "linear", "probit": "probit", "fractional": "probit", "exponential": "poisson"}[omodel]

    def resid_t(g):
        return t - stats.norm.cdf(Z @ g)

    def mean(W, a):
        return _mean(kind, W @ a)[0] if kind != "linear" else W @ a

    def score(W, yv, a):
        if kind == "linear":
            return W * (yv - W @ a)[:, None]
        return glm_score(kind, W, yv, a)

    g = glm_fit("probit", Z, t)
    W = np.column_stack([Xo, resid_t(g)])
    arms = []
    for tj in (t0, t):
        m = tj > 0
        if kind == "linear":
            arms.append(np.linalg.lstsq(W[m], y[m], rcond=None)[0])
        else:
            arms.append(glm_fit(kind, W[m], y[m]))
    c0, c1 = arms
    mu = np.array([np.sum(d * mean(W, c0)) / d.sum(), np.sum(d * mean(W, c1)) / d.sum()])
    theta = np.concatenate([mu, c0, c1, g])
    k1 = ko + 1

    # the variance takes the probit coefficients as known: the sandwich of the outcome and POM
    # moments at the estimated residual reproduces Stata's standard errors (compat 1908)
    def gfun(th):
        mu_, a0, a1 = th[:2], th[2:2 + k1], th[2 + k1:2 + 2 * k1]
        s0 = t0[:, None] * score(W, y, a0)
        s1 = t[:, None] * score(W, y, a1)
        mom = np.column_stack([d * (mean(W, a0) - mu_[0]), d * (mean(W, a1) - mu_[1])])
        return np.column_stack([mom, s0, s1])

    V = sandwich_m(gfun, theta[:2 + 2 * k1], n, clusters)
    return theta, V, k1


# ---------------------------------------------------------------------------
# nearest-neighbour matching
# ---------------------------------------------------------------------------

def _dist_matrix(A, B, metric, S=None):
    if metric == "euclidean":
        M = np.eye(A.shape[1])
    elif metric == "ivariance":
        M = np.diag(1 / np.diag(S))
    else:
        M = np.linalg.inv(S)
    D = np.zeros((A.shape[0], B.shape[0]))
    for i in range(A.shape[0]):
        diff = B - A[i]
        D[i] = np.einsum("ij,jk,ik->i", diff, M, diff)
    return D


def _knn(D, k):
    """For each row of D the indices of its k nearest columns, ties included."""
    out = []
    for row in D:
        order = np.argsort(row, kind="stable")
        kth = row[order[min(k, len(row)) - 1]]
        out.append(np.flatnonzero(row <= kth + 1e-10 * max(1.0, abs(kth))))
    return out


def nnmatch(y, t, X, stat="ate", nn=1, metric="maha", vce_nn=2, S=None):
    """Matching estimator of the ATE or ATET with Abadie-Imbens robust standard errors. Returns
    (tau, se, mu, (min matches, max matches)) with mu the potential-outcome means of the target
    group (point estimates only). VERIFY: ties and the variance."""
    n = len(y)
    if S is None:
        S = np.cov(X.T, bias=False).reshape(X.shape[1], X.shape[1])
    one = np.flatnonzero(t > 0)
    zero = np.flatnonzero(t == 0)
    yh = np.zeros((n, 2))                # imputed outcomes (control, treated)
    yh[one, 1] = y[one]
    yh[zero, 0] = y[zero]
    K = np.zeros(n)
    counts = []
    D10 = _dist_matrix(X[one], X[zero], metric, S)    # treated -> controls
    D01 = _dist_matrix(X[zero], X[one], metric, S)    # controls -> treated
    for rows, cols, tgt in ((one, zero, 0), (zero, one, 1)):
        D = D10 if tgt == 0 else D01
        for ii, nb in zip(rows, _knn(D, nn)):
            w = 1.0 / len(nb)
            counts.append(len(nb))
            yh[ii, tgt] = y[cols[nb]].mean()
            K[cols[nb]] += w
    # conditional variances from nearest same-group neighbours
    s2 = np.zeros(n)
    for grp in (one, zero):
        Dg = _dist_matrix(X[grp], X[grp], metric, S)
        np.fill_diagonal(Dg, np.inf)
        J = vce_nn
        for ii, nb in zip(grp, _knn(Dg, J)):
            s2[ii] = J / (J + 1) * (y[ii] - y[grp[nb]].mean()) ** 2
    diff = yh[:, 1] - yh[:, 0]
    M = float(nn)
    if stat == "ate":
        tau = diff.mean()
        var = np.sum((diff - tau) ** 2 + (K ** 2 + ((2 * M - 1) / M) * K) * s2) / n ** 2
        mu = (yh[:, 0].mean(), yh[:, 1].mean())
    else:
        n1 = len(one)
        tau = diff[one].mean()
        var = (np.sum((diff[one] - tau) ** 2) + np.sum(K[zero] * (K[zero] - 1 / M) * s2[zero])) / n1 ** 2
        mu = (yh[one, 0].mean(), y[one].mean())
    return float(tau), float(np.sqrt(max(var, 0.0))), mu, (min(counts), max(counts))


# ---------------------------------------------------------------------------
# multivalued treatments (multinomial logit treatment model)
# ---------------------------------------------------------------------------

def mlogit_prob(Z, g, L):
    """P(t = l), l = 0..L-1, with level 0 as the base of the multinomial logit; g stacks the
    (L - 1) coefficient vectors."""
    k = Z.shape[1]
    E = np.column_stack([np.zeros(len(Z))] + [Z @ g[(l - 1) * k:l * k] for l in range(1, L)])
    E -= E.max(axis=1, keepdims=True)
    P = np.exp(E)
    return P / P.sum(axis=1, keepdims=True)


def mlogit_score(Z, tcode, g, L):
    P = mlogit_prob(Z, g, L)
    return np.column_stack([Z * ((tcode == l) - P[:, l])[:, None] for l in range(1, L)])


def mlogit_fit(Z, tcode, L, tol=1e-12, maxit=100):
    """Multinomial logit by Newton-Raphson."""
    n, k = Z.shape
    g = np.zeros(k * (L - 1))
    for _ in range(maxit):
        P = mlogit_prob(Z, g, L)
        grad = mlogit_score(Z, tcode, g, L).sum(axis=0)
        H = np.zeros((len(g), len(g)))
        for a in range(1, L):
            for b in range(1, L):
                wab = P[:, a] * ((a == b) - P[:, b])
                H[(a - 1) * k:a * k, (b - 1) * k:b * k] = -(Z * wab[:, None]).T @ Z
        step = np.linalg.solve(H, grad)
        g = g - step
        if np.max(np.abs(step)) < tol:
            break
    return g


def teffects_multi(est, y, tcode, L, Xo, Z, omodel, stat, tl, clusters=None):
    """teffects ra, ipw, ipwra and aipw with L >= 2 treatment levels coded 0..L-1 (multinomial
    logit treatment model). stat: ate, atet (treated level tl) or pomeans. theta = [mu_0 ..
    mu_{L-1}, beta_0 .. beta_{L-1}, gamma]; returns (theta, V)."""
    n, ko = Xo.shape
    om = GLMOutcome(omodel, y) if isinstance(omodel, str) else omodel
    use_out, use_trt = est != "ipw", est != "ra"
    kb = (ko + getattr(om, "kaux", 0)) if use_out else 0
    kz = Z.shape[1] * (L - 1) if use_trt else 0
    T_ = np.column_stack([(tcode == l).astype(float) for l in range(L)])
    d = np.ones(n) if stat != "atet" else T_[:, tl]

    def weights(P):
        if stat == "atet":
            return T_ * (P[:, [tl]] / P)
        return T_ / P

    g = mlogit_fit(Z, tcode, L) if use_trt else np.zeros(0)
    W = weights(mlogit_prob(Z, g, L)) if use_trt else None
    bs = []
    if use_out:
        for l in range(L):
            bs.append(om.fit(Xo, T_[:, l] if est in ("ra", "aipw") else W[:, l]))
    if est == "ipw":
        mu = np.array([np.sum(W[:, l] * y) / np.sum(W[:, l]) for l in range(L)])
    else:
        M_ = np.column_stack([om.mean(Xo, b) for b in bs])
        if est == "aipw":
            mu = (M_ + W * (y[:, None] - M_)).mean(axis=0)
        else:
            mu = (d[:, None] * M_).sum(axis=0) / d.sum()
    theta = np.concatenate([mu] + bs + [g])

    def gfun(th):
        mu_ = th[:L]
        bs_ = [th[L + l * kb:L + (l + 1) * kb] for l in range(L)] if use_out else []
        g_ = th[L + L * kb:]
        cols = []
        if use_trt:
            W_ = weights(mlogit_prob(Z, g_, L))
        if est == "ipw":
            mom = W_ * (y[:, None] - mu_[None, :])
        else:
            M2 = np.column_stack([om.mean(Xo, b) for b in bs_])
            a = T_ if est in ("ra", "aipw") else W_
            cols = [a[:, [l]] * om.score(Xo, bs_[l]) for l in range(L)]
            if est == "aipw":
                mom = M2 + W_ * (y[:, None] - M2) - mu_[None, :]
            else:
                mom = d[:, None] * (M2 - mu_[None, :])
        parts = [mom] + cols + ([mlogit_score(Z, tcode, g_, L)] if use_trt else [])
        return np.column_stack(parts)

    V = sandwich_m(gfun, theta, n, clusters)
    return theta, V
