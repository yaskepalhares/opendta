"""Fitting machinery of sem: groups, missing values and the ADF discrepancy.

A `Problem` holds one `Model` per group, the data summary of each group, and the map from each
group's free parameters to the global parameter vector (parameters constrained equal across
groups share one global entry). The objective is the log likelihood summed over groups for
method(ml) and method(mlmv) (the data summary is one pattern for ml and one per missing-value
pattern for mlmv), and minus the sum over groups of N_g F_WLS for method(adf), where
F_WLS = (v - tau)' W^-1 (v - tau), v = [z; vech(S)], tau = [mu; vech(Sigma)] and W is the sample
covariance of the observation-level moment contributions ([SEM] methods and formulas for sem).
The details Stata does not document (starting values, the iteration log) are marked VERIFY."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from ..core import mvnmiss
from ..core.errors import StataError
from .sem import (Model, _dmoments, _g, _grad_blocks, _loglik_blocks, _start, baseline_ex_correlated)

PCLASSES = ("scoef", "scons", "mcoef", "mcons", "serrvar", "merrvar", "smerrcov", "meanex", "covex")


@dataclass
class GroupData:
    label: float
    N: int
    pats: list                       # data summary: one pattern (ml) or one per missing pattern
    mean: np.ndarray                 # saturated mean
    cov: np.ndarray                  # saturated covariance (divided by N)
    ll_s: float                      # saturated log likelihood
    ll_b: float                      # baseline log likelihood
    X: np.ndarray | None = None      # complete raw data (adf)
    adf: dict | None = None          # v, Wi and the baseline discrepancy (adf)


# ---------------------------------------------------------------------------
# data summaries
# ---------------------------------------------------------------------------

def _structure_blocks(model: Model):
    """Column blocks of the baseline model: the observed exogenous variables may be correlated
    with each other when the model has observed endogenous variables."""
    obs = model.spec.observed
    ex = [i for i, n in enumerate(obs) if n in model.exog]
    return [ex] if (ex and baseline_ex_correlated(model)) else []


def build_group_data(model: Model, X: np.ndarray, method: str, label: float = 0.0) -> GroupData:
    """Summarize the rows of X (columns = model.spec.observed; NaN = missing) for the method."""
    N, p = X.shape
    if method == "mlmv":
        pats = mvnmiss.patterns(X)
        N = sum(q.n for q in pats)
        mean, cov, ll_s, _ = mvnmiss.em_mvn(X)
        _, _, ll_b, _ = mvnmiss.em_mvn(X, blocks=_structure_blocks(model))
        return GroupData(label, N, pats, mean, cov, ll_s, ll_b)
    mean = X.mean(axis=0)
    D = X - mean
    cov = D.T @ D / N
    pats = [mvnmiss.Pattern(np.arange(p), np.zeros(0, dtype=int), N, mean, cov)]
    ll_s = float(-0.5 * N * (p * math.log(2 * math.pi) + np.linalg.slogdet(cov)[1] + p))
    Sig = np.diag(np.diag(cov)).copy()
    for blk in _structure_blocks(model):
        Sig[np.ix_(blk, blk)] = cov[np.ix_(blk, blk)]
    ll_b = float(-0.5 * N * (p * math.log(2 * math.pi) + np.linalg.slogdet(Sig)[1]
                             + np.trace(cov @ np.linalg.inv(Sig))))
    gd = GroupData(label, N, pats, mean, cov, ll_s, ll_b)
    if method == "adf":
        gd.X = X
        gd.adf = _adf_summary(model, X, mean, cov)
    return gd


def _vech(M: np.ndarray) -> np.ndarray:
    return M[np.tril_indices(M.shape[0])]


def _adf_summary(model: Model, X: np.ndarray, mean, cov) -> dict:
    N, p = X.shape
    D = X - mean
    il = np.tril_indices(p)
    H = np.column_stack([D, D[:, il[0]] * D[:, il[1]] - _vech(cov)])
    W = H.T @ H / N
    q = W.shape[0]
    if np.linalg.matrix_rank(W) < q:
        raise StataError(498, "the weight matrix of method(adf) is singular; more observations "
                              "than moments are required")   # VERIFY
    Wi = np.linalg.inv(W)
    v = np.concatenate([mean, _vech(cov)])
    # baseline model: tau = A phi is linear in the free means, variances and covariances
    cols = []
    for i in range(p):
        e = np.zeros(q)
        e[i] = 1
        cols.append(e)
    allowed = {(i, i) for i in range(p)}
    for blk in _structure_blocks(model):
        allowed |= {(i, j) for i in blk for j in blk if i > j}
    for k, (i, j) in enumerate(zip(*il)):
        if (int(i), int(j)) in allowed:
            e = np.zeros(q)
            e[p + k] = 1
            cols.append(e)
    A = np.column_stack(cols)
    phi = np.linalg.solve(A.T @ Wi @ A, A.T @ Wi @ v)
    r = v - A @ phi
    return dict(v=v, Wi=Wi, F_b=float(r @ Wi @ r))


# ---------------------------------------------------------------------------
# parameter classes and the map to the global vector
# ---------------------------------------------------------------------------

def item_class(model: Model, item: tuple[str, str, str]) -> str:
    """The parameter class of a free item (for ginvariant())."""
    kind, a, b = item
    spec = model.spec
    endog_obs = [n for n in model.endog if n in spec.observed]
    meas = {n for n in endog_obs if all(x in spec.latent for x, y in spec.paths if y == n)}

    def error_kind(n):
        return "m" if n in meas else "s"
    if kind == "path":
        return "mcoef" if (a in spec.latent and b in meas) else "scoef"
    if kind == "mean":
        if a in model.exog:
            return "meanex"
        return "mcons" if a in meas else "scons"
    if kind == "var":
        if a in model.exog:
            return "covex"
        return "merrvar" if a in meas else "serrvar"
    if a in model.exog and b in model.exog:
        return "covex"
    ka, kb = error_kind(a) if a in model.endog else "x", error_kind(b) if b in model.endog else "x"
    if ka == kb == "m":
        return "merrvar"
    if ka == kb == "s":
        return "serrvar"
    return "smerrcov"


def parse_ginvariant(text: str | None) -> set[str]:
    if text is None or not str(text).strip():
        return {"mcoef", "mcons"}
    out: set[str] = set()
    for w in str(text).split():
        w = w.lower()
        if w == "none":
            continue
        if w == "all":
            out |= set(PCLASSES)
            continue
        hit = [c for c in PCLASSES if c == w]
        if not hit:
            raise StataError(198, f"invalid ginvariant() class {w}")   # VERIFY
        out.add(hit[0])
    return out


@dataclass
class Problem:
    models: list[Model]
    data: list[GroupData]
    method: str
    inv: set[str] = field(default_factory=set)
    gmap: list[np.ndarray] = field(default_factory=list)
    nglob: int = 0
    classes: list[list[str]] = field(default_factory=list)
    ilog: list[float] = field(default_factory=list)
    split: set = field(default_factory=set)        # (group, item) pairs freed from ties and invariance
    split_all: set = field(default_factory=set)    # items freed in every group

    def __post_init__(self):
        G = len(self.models)
        self.classes = [[item_class(m, it) for it in m.free] for m in self.models]
        keys: dict = {}
        self.gmap = []
        for g, m in enumerate(self.models):
            ids = []
            for it, cl in zip(m.free, self.classes[g]):
                if it in m.ties and (g, it) not in self.split and it not in self.split_all:
                    key = ("label", m.ties[it]) if G == 1 else ("label", g, m.ties[it])
                elif it in self.split_all or (g, it) in self.split:
                    key = ("split", g, it)
                else:
                    key = (cl, it) if (G > 1 and cl in self.inv) or G == 1 else (g, it)
                if key not in keys:
                    keys[key] = len(keys)
                ids.append(keys[key])
            self.gmap.append(np.asarray(ids, dtype=int))
        self.nglob = len(keys)

    # global <-> group parameters
    def group_theta(self, th: np.ndarray, g: int) -> np.ndarray:
        return th[self.gmap[g]]

    def expand(self, th: np.ndarray) -> np.ndarray:
        return np.concatenate([self.group_theta(th, g) for g in range(len(self.models))])

    def expand_matrix(self) -> np.ndarray:
        """T such that the expanded vector is T @ (global vector)."""
        rows = sum(len(m) for m in self.gmap)
        T = np.zeros((rows, self.nglob))
        r = 0
        for ids in self.gmap:
            for i in ids:
                T[r, i] = 1.0
                r += 1
        return T

    # objective --------------------------------------------------------------
    def f(self, th: np.ndarray) -> float:
        tot = 0.0
        for g, (m, d) in enumerate(zip(self.models, self.data)):
            t = self.group_theta(th, g)
            if self.method == "adf":
                tot -= d.N * _adf_F(m, t, d)
            else:
                v = _loglik_blocks(m, t, d.pats)
                if not np.isfinite(v):
                    return -np.inf
                tot += v
        return tot

    def grad(self, th: np.ndarray) -> np.ndarray:
        out = np.zeros(self.nglob)
        for g, (m, d) in enumerate(zip(self.models, self.data)):
            t = self.group_theta(th, g)
            if self.method == "adf":
                gg = _adf_grad(m, t, d)
            else:
                gg = _grad_blocks(m, t, d.pats)
            np.add.at(out, self.gmap[g], gg)
        return out

    def hess(self, th: np.ndarray) -> np.ndarray:
        if self.method == "adf":
            H = np.zeros((self.nglob, self.nglob))
            for g, (m, d) in enumerate(zip(self.models, self.data)):
                J = _adf_jacobian(m, self.group_theta(th, g))
                Hg = -2 * d.N * J.T @ d.adf["Wi"] @ J
                ids = self.gmap[g]
                np.add.at(H, (ids[:, None], ids[None, :]), Hg)
            return (H + H.T) / 2
        k = len(th)
        H = np.zeros((k, k))
        for i in range(k):
            e = np.zeros(k)
            hh = 1e-5 * (abs(th[i]) + 1)
            e[i] = hh
            H[:, i] = (self.grad(th + e) - self.grad(th - e)) / (2 * hh)
        return (H + H.T) / 2

    def start(self) -> np.ndarray:
        acc = np.zeros(self.nglob)
        cnt = np.zeros(self.nglob)
        for g, (m, d) in enumerate(zip(self.models, self.data)):
            t0 = _start(m, d.mean, d.cov)
            np.add.at(acc, self.gmap[g], t0)
            np.add.at(cnt, self.gmap[g], 1)
        return acc / np.maximum(cnt, 1)


# ---------------------------------------------------------------------------
# ADF pieces
# ---------------------------------------------------------------------------

def _tau_and_jacobian(model: Model, th: np.ndarray):
    mu, Sig, dm = _dmoments(model, th)
    p = model.p
    il = np.tril_indices(p)
    tau = np.concatenate([mu, _vech(Sig)])
    J = np.zeros((len(tau), len(th)))
    for k, (dmv, dSv) in enumerate(dm):
        if dmv is not None:
            J[:p, k] = dmv
        if dSv is not None:
            J[p:, k] = dSv[il]
    return tau, J


def _adf_jacobian(model: Model, th: np.ndarray) -> np.ndarray:
    return _tau_and_jacobian(model, th)[1]


def _adf_F(model: Model, th: np.ndarray, d: GroupData) -> float:
    try:
        mu, Sig, _ = model.implied(th)
    except np.linalg.LinAlgError:
        return np.inf
    tau = np.concatenate([mu, _vech(Sig)])
    r = d.adf["v"] - tau
    return float(r @ d.adf["Wi"] @ r)


def _adf_grad(model: Model, th: np.ndarray, d: GroupData) -> np.ndarray:
    """Gradient of -N F_WLS: 2 N J' W^-1 (v - tau)."""
    tau, J = _tau_and_jacobian(model, th)
    r = d.adf["v"] - tau
    return 2 * d.N * J.T @ d.adf["Wi"] @ r


# ---------------------------------------------------------------------------
# maximization
# ---------------------------------------------------------------------------

def fit_problem(prob: Problem, s, log: bool, maxiter: int = 300):
    """Quasi-Newton to approach the optimum, then Newton steps (Gauss-Newton for ADF) to finish.
    Returns (theta, objective, V, converged); V is the inverse of minus the Hessian for ml and
    mlmv, and (N J' W^-1 J)^-1 summed over groups for adf."""
    from scipy import optimize
    f, gr = prob.f, prob.grad
    th = prob.start()
    ll = f(th)
    if not np.isfinite(ll):
        raise StataError(430, "initial values not feasible")   # VERIFY
    lines = [ll]
    # the quasi-Newton phase works on the log of the variances, which keeps them positive (a
    # model with simple structure can otherwise run off to a spurious optimum with negative
    # variances); the Newton steps that finish are unconstrained
    vi = np.zeros(prob.nglob, dtype=bool)
    for g, m in enumerate(prob.models):
        for k, it in enumerate(m.free):
            if it[0] == "var":
                vi[prob.gmap[g][k]] = True
    u0 = th.copy()
    u0[vi] = np.log(np.maximum(th[vi], 1e-8))

    def to_theta(u):
        t = u.copy()
        t[vi] = np.exp(np.clip(u[vi], -50, 50))
        return t

    def fu(u):
        v = f(to_theta(u))
        return -v if np.isfinite(v) else 1e100

    def gu(u):
        t = to_theta(u)
        gg = -gr(t)
        gg[vi] = gg[vi] * t[vi]
        return gg
    res = optimize.minimize(fu, u0, method="BFGS", jac=gu, options={"maxiter": maxiter, "gtol": 1e-7})
    th = to_theta(res.x)
    lines.append(f(th))
    converged = False
    for _ in range(50):
        g = gr(th)
        H = prob.hess(th)
        try:
            step = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            break
        lam = 1.0
        base = f(th)
        while lam > 1e-6 and not (f(th + lam * step) >= base - 1e-12):
            lam /= 2
        th = th + lam * step
        lines.append(f(th))
        if abs(lines[-1] - lines[-2]) < 1e-10 and np.max(np.abs(lam * step)) < 1e-8:
            converged = True
            break
    shown = [lines[0]] + [x for i, x in enumerate(lines[1:], 1)
                          if abs(x - lines[i - 1]) > 1e-9 or i == len(lines) - 1]
    if log:
        out = s.output
        out.write("\nFitting target model:\n\n", "text")
        word = "discrepancy" if prob.method == "adf" else "log likelihood"
        for i, x in enumerate(shown):
            val = -x if prob.method == "adf" else x
            out.write(f"{'Iteration ' + str(i) + ':':<14} {word} = {_g(val, '%10.0g'):>10}  \n", "text")
    H = prob.hess(th)
    try:
        V = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        V = np.linalg.pinv(-H)
    if prob.method == "adf":
        V = 2 * V        # H is -2 N J'WiJ; the covariance is (N J'WiJ)^-1
    prob.ilog = shown
    return th, f(th), (V + V.T) / 2, converged
