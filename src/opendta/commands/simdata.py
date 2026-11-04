"""drawnorm and corr2data ([D] drawnorm, corr2data): data drawn from, or built to have, a given
mean vector and correlation or covariance matrix."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Dataset, Variable
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options

from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _moments(s: "Session", o: dict, k: int):
    """means, standard deviations and correlation (or covariance) from the options."""
    def vec(name):
        v = o.get(name)
        if not isinstance(v, str):
            return None
        return _matrix_arg(s, v)
    means = vec("means")
    means = np.zeros(k) if means is None else np.ravel(means)
    if isinstance(o.get("cov"), str):
        cov = _square(s, str(o["cov"]), k, o)
    else:
        sds = vec("sds")
        sds = np.ones(k) if sds is None else np.ravel(sds)
        corr = _square(s, str(o["corr"]), k, o) if isinstance(o.get("corr"), str) else np.eye(k)
        cov = corr * np.outer(sds, sds)
    if len(means) != k or cov.shape != (k, k):
        raise StataError(503, "conformability error")
    return means, cov


def _matrix_arg(s: "Session", text: str) -> np.ndarray:
    t = text.strip()
    mats = getattr(s, "matrices", {})
    if t in mats:
        return np.asarray(mats[t].data, float)
    words = t.replace(",", " ").split()
    try:
        return np.array([[float(w) for w in words]], float)      # a list of numbers
    except ValueError:
        pass
    if t.startswith("("):
        rows = [r.split(",") if "," in r else r.split() for r in t.strip("()").split("\\")]
        return np.array([[float(s.eval(x)) for x in r if x.strip()] for r in rows], float)
    raise StataError(111, f"matrix {t} not found")


def _square(s: "Session", text: str, k: int, o: dict) -> np.ndarray:
    m = _matrix_arg(s, text)
    store = str(o.get("cstorage", "full")) if isinstance(o.get("cstorage"), str) else "full"
    if store != "full" or m.shape != (k, k):
        vals = np.ravel(m)
        full = np.zeros((k, k))
        idx = np.tril_indices(k) if store in ("lower", "full") else np.triu_indices(k)
        if len(vals) != len(idx[0]):
            raise StataError(503, "conformability error")
        full[idx] = vals
        m = np.where(full != 0, full, full.T) if store == "lower" else full + full.T - np.diag(np.diag(full))
        if store != "upper":
            m = full + full.T - np.diag(np.diag(full))
    return m


def _chol(cov: np.ndarray, force: bool) -> np.ndarray:
    try:
        return np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        if not force:
            raise StataError(506, "matrix not positive definite")
        w, V = np.linalg.eigh(cov)
        return V @ np.diag(np.sqrt(np.maximum(w, 0)))


def _prepare_obs(s: "Session", o: dict) -> int:
    ds = s.data
    if isinstance(o.get("n"), str):
        n = int(float(s.eval(str(o["n"]))))
        if ds.nobs and n != ds.nobs:
            if not o.get("clear"):
                raise StataError(4, "no; data in memory would be lost")   # VERIFY
            s.data = Dataset()
            ds = s.data
        if n > ds.nobs:
            ds.set_obs(n)
            s.output.write(f"(obs {n:,})\n", "text")      # compat 1501
        return n
    if o.get("clear"):
        s.data = Dataset()
    if not s.data.nobs:
        raise StataError(2000, "no observations")   # VERIFY: n() required with empty data
    return s.data.nobs


_GEN_OPTS = {"clear": 5, "double": 3, "n": 1, "sds": 3, "corr": 4, "cov": 3, "cstorage": 2,
             "forcepsd": 5, "means": 2, "seed": 4}


@command("drawnorm", "drawnorm")
def cmd_drawnorm(s: "Session", args: str) -> None:
    from ..core import rng
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    names = body.split()
    if not names:
        raise StataError(198, "invalid syntax")
    o = match_options(opt, _GEN_OPTS)
    if isinstance(o.get("seed"), str):
        s.run(f"quietly set seed {o['seed']}")
    k = len(names)
    means, cov = _moments(s, o, k)
    L = _chol(cov, bool(o.get("forcepsd")))
    n = _prepare_obs(s, o)
    # standard normal draws, one variable after the other (VERIFY: Stata's order of draws)
    Z = np.column_stack([rng.draw("rnormal", [0.0, 1.0], n) for _ in range(k)])
    X = Z @ L.T + means
    ds = s.data
    vt = "double" if o.get("double") else s.default_type()
    for j, nm in enumerate(names):
        col = np.full(ds.nobs, SYS)
        col[:n] = X[:, j]
        ds.add(Variable(nm, vt, col))
    s.notify_state()


@command("corr2data", "corr2data")
def cmd_corr2data(s: "Session", args: str) -> None:
    from ..core import rng
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    names = body.split()
    if not names:
        raise StataError(198, "invalid syntax")
    o = match_options(opt, _GEN_OPTS)
    k = len(names)
    means, cov = _moments(s, o, k)
    L = _chol(cov, bool(o.get("forcepsd")))
    n = _prepare_obs(s, o)
    if n <= k:
        raise StataError(2001, "insufficient observations")   # VERIFY
    Z = np.column_stack([rng.draw("rnormal", [0.0, 1.0], n) for _ in range(k)])
    Z = Z - Z.mean(0)
    S = np.cov(Z, rowvar=False).reshape(k, k)
    Z = Z @ np.linalg.inv(np.linalg.cholesky(S)).T      # exactly identity sample covariance
    X = Z @ L.T + means
    ds = s.data
    vt = "double" if o.get("double") else s.default_type()
    for j, nm in enumerate(names):
        col = np.full(ds.nobs, SYS)
        col[:n] = X[:, j]
        ds.add(Variable(nm, vt, col))
    s.notify_state()
