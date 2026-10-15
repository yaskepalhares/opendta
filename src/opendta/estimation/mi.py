"""Multiple imputation: mi set, mi register, mi impute, mi estimate,
mi passive, mi xeq, mi extract, mi describe, mi misstable ([MI] manual).

Storage follows Stata's mlong and flong styles: _mi_m (0 for the original
data), _mi_id (observation identifier) and _mi_miss (1 for observations
incomplete in m=0), with the settings in the _dta characteristics.

mi impute follows [MI] mi impute regress/pmm/logit/ologit/mlogit/poisson,
Methods and formulas: the parameters are drawn from their posterior (scaled
chi-square for the variance, normal for the coefficients) and then the
missing values. The draws come from OpenDTA's random-number generator in
an order Stata does not document, so imputed values differ from Stata's
(like bootstrap; compat/known_differences.txt). mi estimate combines the
estimates with Rubin's rules, with the small-sample degrees of freedom
of Barnard and Rubin (1999) when the command has residual degrees of
freedom, and the model F test of Li, Raghunathan and Rubin (1991).
"""

from __future__ import annotations

import math
import time
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import find_top, match_options
from ..lazy import lazy

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..core.dataset import Dataset
    from ..session import Session

SYS = M.SYSMISS
SYSV = ("_mi_m", "_mi_id", "_mi_miss")


def _g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


# ---------------------------------------------------------------------------
# settings
# ---------------------------------------------------------------------------

def _ch(ds: "Dataset") -> dict:
    return ds.chars.setdefault("_dta", {})


def mi_style(ds: "Dataset") -> str:
    ch = ds.chars.get("_dta", {})
    return ch.get("_mi_style", "") if ch.get("_mi_marker") else ""


def require_mi(s: "Session") -> str:
    st_ = mi_style(s.data)
    if not st_:
        raise StataError(119, "no; data not mi set")   # VERIFY
    return st_


def _list(ds, key: str) -> list[str]:
    return [v for v in _ch(ds).get(key, "").split() if ds.has(v)]


def _set_list(ds, key: str, names: list[str]) -> None:
    ch = _ch(ds)
    if names:
        ch[key] = " ".join(dict.fromkeys(names))
    else:
        ch.pop(key, None)


def nimp(ds) -> int:
    return int(_ch(ds).get("_mi_M", "0") or 0)


def _touch(ds) -> None:
    _ch(ds)["_mi_update"] = str(time.time())


def _num(ds, name) -> np.ndarray:
    return ds.get(name).data.astype(np.float64)


def _update_miss(ds) -> None:
    """_mi_miss: 1 in m=0 for observations with a missing imputed variable."""
    m = _num(ds, "_mi_m")
    ivars = _list(ds, "_mi_ivars")
    miss = np.zeros(ds.nobs, dtype=bool)
    for v in ivars:
        x = ds.get(v)
        miss |= (x.data == "") if x.is_string else (x.data >= SYS)
    vals = np.where(m == 0, miss.astype(float), SYS)
    ds.get("_mi_miss").data = vals


# ---------------------------------------------------------------------------
# the dataset of imputation m
# ---------------------------------------------------------------------------

def dataset_m(ds: "Dataset", m: int) -> "Dataset":
    """A copy of the data of imputation m (without the system variables)."""
    style = mi_style(ds)
    mm = _num(ds, "_mi_m")
    if style == "flong" or m == 0:
        rows = np.nonzero(mm == m)[0]
    else:
        ids = _num(ds, "_mi_id")
        base = np.nonzero(mm == 0)[0]
        imp = np.nonzero(mm == m)[0]
        repl = dict(zip(ids[imp].tolist(), imp.tolist()))
        rows = np.array([repl.get(i, r) for i, r in zip(ids[base].tolist(), base.tolist())], dtype=int)
    new = ds.copy()
    keep = np.zeros(ds.nobs, dtype=bool)
    keep[rows] = True
    # keep the rows in the order of m=0 (rows may come from the end of the data)
    order = rows
    for v in new.vars:
        v.raw = v.raw[order]
    new.nobs = len(order)
    new.drop_vars([v for v in SYSV if new.has(v)])
    new.esample = None
    ch = new.chars.get("_dta", {})
    for k in [k for k in ch if k.startswith("_mi_")]:
        ch.pop(k)
    return new


class _with_data:
    """Runs commands on another dataset and restores the original one."""

    def __init__(self, s: "Session", ds):
        self.s, self.ds = s, ds

    def __enter__(self):
        self.old = self.s.data
        self.s.data = self.ds
        return self.ds

    def __exit__(self, *exc):
        self.s.data = self.old
        return False


# ---------------------------------------------------------------------------
# mi dispatcher
# ---------------------------------------------------------------------------

def cmd_mi(s: "Session", args: str) -> None:
    text = args.strip()
    if not text:
        raise StataError(198, "subcommand required")   # VERIFY
    head = text.split(None, 1)[0].rstrip(":,")
    rest = text[len(text.split(None, 1)[0]):].strip() if " " in text else ""
    if text.split(None, 1)[0].endswith(":"):
        rest = ": " + rest
    table = [("set", 3, mi_set), ("register", 3, mi_register), ("unregister", 5, mi_unregister),
             ("describe", 1, mi_describe), ("impute", 3, mi_impute), ("estimate", 3, mi_estimate),
             ("passive", 3, mi_passive), ("xeq", 3, mi_xeq), ("extract", 3, mi_extract),
             ("misstable", 5, mi_misstable), ("query", 1, mi_query)]
    for name, k, fn in table:
        if len(head) >= k and name.startswith(head):
            if name in ("estimate", "passive", "xeq") and ":" in text:
                i = text.index(":")
                fn(s, text[len(text.split(None, 1)[0]):i].strip().lstrip(",").strip()
                   if len(text.split(None, 1)[0]) <= i else "", text[i + 1:].strip())
            elif name in ("estimate", "passive", "xeq"):
                raise StataError(198, "colon required")   # VERIFY
            else:
                fn(s, rest)
            return
    raise StataError(198, f"mi {head}: invalid subcommand")   # VERIFY


def mi_set(s: "Session", args: str) -> None:
    from ..core.dataset import Variable
    ds = s.data
    style = args.strip().lower()
    if style in ("m", "ml", "mlo", "mlon", "mlong"):
        style = "mlong"
    elif style in ("flong",):
        style = "flong"
    elif style in ("wide", "w", "flongsep"):
        raise StataError(198, f"mi set {style}: OpenDTA supports the mlong and flong styles")   # OpenDTA
    else:
        raise StataError(198, "mi set: style required (mlong or flong)")   # VERIFY
    if mi_style(ds):
        raise StataError(119, "data already mi set")   # VERIFY
    n = ds.nobs
    ds.add(Variable("_mi_m", "int", np.zeros(n)))
    ds.add(Variable("_mi_id", "long" if n > 32000 else "int", np.arange(1, n + 1, dtype=float)))
    ds.add(Variable("_mi_miss", "byte", np.zeros(n)))
    ch = _ch(ds)
    ch.update({"_mi_marker": "_mi_ds_1", "_mi_style": style, "_mi_M": "0"})
    _touch(ds)
    s.notify_state()


def _reg(s, args, key):
    ds = s.data
    require_mi(s)
    words = args.split()
    if not words:
        raise StataError(198, "invalid syntax")
    kind = words[0]
    names = s.expand_varlist(" ".join(words[1:]))
    full = next((k for k in ("imputed", "passive", "regular") if len(kind) >= 3 and k.startswith(kind)), None)
    if full is None:
        raise StataError(198, f"mi register {kind}: invalid")   # VERIFY
    return full, names


def mi_register(s: "Session", args: str) -> None:
    full, names = _reg(s, args, "")
    ds = s.data
    for nm in names:
        if nm in SYSV:
            raise StataError(198, f"{nm} is a system variable")   # VERIFY
    for key in ("_mi_ivars", "_mi_pvars", "_mi_rvars"):
        _set_list(ds, key, [v for v in _list(ds, key) if v not in names])
    key = {"imputed": "_mi_ivars", "passive": "_mi_pvars", "regular": "_mi_rvars"}[full]
    _set_list(ds, key, _list(ds, key) + names)
    if full == "imputed":
        _update_miss(ds)
    _touch(ds)


def mi_unregister(s: "Session", args: str) -> None:
    require_mi(s)
    ds = s.data
    names = s.expand_varlist(args)
    for key in ("_mi_ivars", "_mi_pvars", "_mi_rvars"):
        _set_list(ds, key, [v for v in _list(ds, key) if v not in names])
    _update_miss(ds)


def mi_query(s: "Session", args: str) -> None:
    ds = s.data
    style = mi_style(ds)
    if not style:
        s.output.write("(data not mi set)\n", "text")
        return
    s.output.write(f"data mi set {style}, M = {nimp(ds)}\n", "text")
    s.r = {"M": float(nimp(ds)), "style": style}


def mi_describe(s: "Session", args: str) -> None:
    style = require_mi(s)
    ds = s.data
    out = s.output
    m = _num(ds, "_mi_m")
    miss = _num(ds, "_mi_miss")
    m0 = m == 0
    n_inc = int(np.sum(miss[m0] == 1))
    n0 = int(m0.sum())
    M_ = nimp(ds)
    ivars, pvars, rvars = _list(ds, "_mi_ivars"), _list(ds, "_mi_pvars"), _list(ds, "_mi_rvars")
    out.write(f"\n  Style:  {style}\n", "text")
    ago = time.time() - float(_ch(ds).get("_mi_update", time.time()))
    stamp = time.strftime("%d%b%Y %H:%M:%S", time.localtime(time.time() - ago)).lower()
    out.write(f"          last mi update {stamp}, {int(ago)} seconds ago\n\n", "text")
    out.write(f"  Obs.:   complete    {n0 - n_inc:>9,}\n", "text")
    out.write(f"          incomplete  {n_inc:>9,}  (M = {M_} imputation{'s' if M_ != 1 else ''})\n", "text")
    out.write("          ---------------------\n", "text")
    out.write(f"          total       {n0:>9,}\n\n", "text")

    def lst(title, names, count=False):
        parts = []
        for v in names:
            if count:
                x = ds.get(v)
                k = int(np.sum(((x.data == "") if x.is_string else (x.data >= SYS))[m0]))
                parts.append(f"{v}({k})")
            else:
                parts.append(v)
        out.write(f"  {'Vars.:  ' if title == 'imputed' else '        '}{title + ':':<9}{len(names):>2}"
                  f"{'; ' + ' '.join(parts) if parts else ''}\n\n", "text")
    lst("imputed", ivars, True)
    lst("passive", pvars)
    lst("regular", rvars)
    out.write(f"          system:   3; {' '.join(SYSV)}\n\n", "text")
    unreg = [v for v in ds.names if v not in ivars + pvars + rvars and v not in SYSV]
    if unreg:
        out.write(f"         (there {'are' if len(unreg) != 1 else 'is'} {len(unreg)} unregistered "
                  f"variable{'s' if len(unreg) != 1 else ''})\n", "text")
    s.r = {"M": float(M_), "N": float(n0), "N_incomplete": float(n_inc), "N_complete": float(n0 - n_inc),
           "ivars": " ".join(ivars)}


def mi_extract(s: "Session", args: str) -> None:
    require_mi(s)
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    m = int(float(body.strip()))
    if m < 0 or m > nimp(s.data):
        raise StataError(198, f"m={m} does not exist")   # VERIFY
    new = dataset_m(s.data, m)
    s.data = new
    new.changed = True
    s.notify_state()


def mi_misstable(s: "Session", args: str) -> None:
    require_mi(s)
    with _with_data(s, dataset_m(s.data, 0)):
        s.run("misstable " + args, echo=False)


def mi_xeq(s: "Session", mlist: str, cmd: str) -> None:
    require_mi(s)
    ds = s.data
    from ..lang.words import parse_numlist
    ms = [int(v) for v in parse_numlist(mlist)] if mlist.strip() else list(range(nimp(ds) + 1))
    out = s.output
    for m in ms:
        out.write(f"\nm={m} data:\n", "text")
        for c in cmd.split(";"):
            if not c.strip():
                continue
            out.write(f"-> {c.strip()}\n", "text")
            with _with_data(s, dataset_m(ds, m)):
                rc = s.run(c.strip(), echo=False)
            if rc:
                raise StataError(rc, "")


def mi_passive(s: "Session", opts: str, cmd: str) -> None:
    """mi passive: generate|replace|egen ... : runs the command on every imputation."""
    require_mi(s)
    ds = s.data
    words = cmd.split()
    if not words or words[0] not in ("generate", "gen", "g", "replace", "egen"):
        raise StataError(198, "mi passive: generate, replace or egen required")   # VERIFY
    target = words[1].split("=")[0]
    if target in ("byte", "int", "long", "float", "double") or target.startswith("str"):
        target = words[2].split("=")[0]
    out = s.output
    mm = _num(ds, "_mi_m")
    M_ = nimp(ds)
    style = mi_style(ds)
    rc = 0
    vals = None
    for m in range(M_ + 1):
        out.write(f"m={m}:\n", "text")
        sub = dataset_m(ds, m)
        if words[0] != "replace" and sub.has(target):
            sub.drop_vars([target])
        with _with_data(s, sub):
            rc = s.run(cmd, echo=False)
        if rc:
            raise StataError(rc, "")
        v = sub.get(target)
        if vals is None:
            from ..core.dataset import Variable
            if not ds.has(target):
                ds.add(Variable(target, v.vtype, np.full(ds.nobs, "" if v.is_string else SYS,
                                                         dtype=object if v.is_string else float)))
            vals = ds.get(target)
        # write back the rows of m
        if m == 0 or style == "flong":
            rows = np.nonzero(mm == m)[0]
            vals_arr = np.array(vals.data, copy=True)
            vals_arr[rows] = v.data
            vals.data = vals_arr
        else:
            ids = _num(ds, "_mi_id")
            imp = np.nonzero(mm == m)[0]
            base_ids = _num(sub, "_mi_id") if sub.has("_mi_id") else None
            pos0 = {int(i): k for k, i in enumerate(ids[mm == 0].tolist())}
            arr = np.array(vals.data, copy=True)
            for r in imp:
                arr[r] = v.data[pos0[int(ids[r])]]
            vals.data = arr
            del base_ids
    _set_list(ds, "_mi_pvars", _list(ds, "_mi_pvars") + [target])
    _touch(ds)
    s.notify_state()


# ---------------------------------------------------------------------------
# mi impute
# ---------------------------------------------------------------------------

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
    raise StataError(198, f"mi impute {method}: not available in OpenDTA")   # OpenDTA


_METHODS = {"regress": 3, "pmm": 3, "logit": 4, "ologit": 5, "mlogit": 5, "poisson": 4, "chained": 5,
            "monotone": 4, "mvn": 3, "truncreg": 5, "intreg": 4, "nbreg": 5}


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
    if meth in ("mvn", "monotone", "truncreg", "intreg", "nbreg"):
        raise StataError(198, f"mi impute {meth} is not available in OpenDTA")   # OpenDTA
    rest = words[1] if len(words) > 1 else ""
    body, opts = _parse_impute(rest)
    o = match_options(opts, {"add": 3, "replace": 3, "rseed": 2, "noisily": 1, "force": 5, "knn": 3,
                             "burnin": 4, "dots": 3, "nolegend": 5, "augment": 3, "orderasis": 6,
                             "noisure": 5}) if opts.strip() else {}
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
    if meth == "chained":
        specs, covars = _parse_chained(s, body)
    else:
        if "=" in body:
            lhs, rhs = body.split("=", 1)
        else:
            lhs, rhs = body, ""
        tgt = s.expand_varlist(lhs.strip())
        if len(tgt) != 1:
            raise StataError(103, "too many variables specified")
        specs = [(meth, tgt[0], {"knn": o.get("knn")})]
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
    for m in ms:
        rows_m = _rows_of(ds, style, m)         # m=0 order, with the rows of m where they exist
        work = {v: _num(ds, v)[rows_m] for _m, v, _o in specs}
        missing = {v: _num(ds, v)[m0] >= SYS for _m, v, _o in specs}
        cov = _design_rows(ds, covars, rows_m)
        if meth == "chained":
            burn = int(o["burnin"]) if o.get("burnin") else 10
            _chained(specs, work, missing, cov, burn)
        else:
            mth, v, so = specs[0]
            mis = missing[v]
            ok = ~mis & np.all(np.isfinite(cov), axis=1)
            if mis.any():
                Xo = np.column_stack([cov[ok], np.ones(int(ok.sum()))])
                Xm = np.column_stack([cov[mis], np.ones(int(mis.sum()))])
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
    title = "Multivariate imputation" if meth == "chained" else "Univariate imputation"
    mname = {"regress": "Linear regression", "pmm": "Predictive mean matching", "logit": "Logistic regression",
             "ologit": "Ordered logistic regression", "mlogit": "Multinomial logistic regression",
             "poisson": "Poisson regression", "chained": "Chained equations"}[meth]
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
    out.write("\n" + "-" * 66 + "\n", "text")
    out.write(f"{'':19}|{'Observations per m':^46}\n", "text")
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


def _cond(mth, v, so, work, missing, cov, preds):
    mis = missing[v]
    if not mis.any():
        return
    cols = [cov] + [work[w][:, None] for w in preds]
    X = np.column_stack(cols + [np.ones((len(mis), 1))])
    ok = ~mis & np.all(np.isfinite(X), axis=1) & np.all(X < SYS, axis=1)
    fill = mis & np.all(np.isfinite(X), axis=1) & np.all(X < SYS, axis=1)
    if not fill.any():
        return
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


# ---------------------------------------------------------------------------
# mi estimate
# ---------------------------------------------------------------------------

def mi_estimate(s: "Session", opts: str, cmd: str) -> None:
    from .results import Estimates, post
    require_mi(s)
    ds = s.data
    o = match_options(opts, {"level": 1, "dots": 3, "noisily": 1, "nosmall": 5, "vartable": 4,
                             "imputations": 3, "cmdok": 5, "post": 4, "mcerror": 3, "dftable": 3}) \
        if opts.strip() else {}
    M_ = nimp(ds)
    if M_ < 2:
        raise StataError(2000, "insufficient imputations")   # VERIFY
    ms = list(range(1, M_ + 1))
    if isinstance(o.get("imputations"), str):
        from ..lang.words import parse_numlist
        ms = [int(v) for v in parse_numlist(o["imputations"])]
    out = s.output
    bs, Vs, Ns, dfr = [], [], [], []
    names = eqs = None
    title, cmdname, vce = "", "", ""
    for m in ms:
        sub = dataset_m(ds, m)
        with _with_data(s, sub):
            if not o.get("noisily"):
                out.quiet_depth += 1
            try:
                rc = s.run(cmd, echo=False)
            finally:
                if not o.get("noisily"):
                    out.quiet_depth -= 1
            if rc:
                raise StataError(rc, "")
            e = s.e
            if "b" not in e:
                raise StataError(301, "estimation command required")   # VERIFY
            b = e["b"].data.ravel().copy()
            V = e["V"].data.copy()
            if names is None:
                names = list(e["b"].colnames)
                eqs = list(getattr(e["b"], "coleq", []) or ["_"] * len(names))
                title = str(e.get("title", "") or "")
                cmdname = str(e.get("cmd", "") or "")
                vce = str(e.get("vce", "") or "")
            bs.append(b)
            Vs.append(V)
            Ns.append(float(e.get("N", SYS)))
            dfr.append(float(e["df_r"]) if isinstance(e.get("df_r"), float) else None)
    Mn = len(ms)
    Bm = np.array(bs)
    Q = Bm.mean(axis=0)
    W = np.mean(Vs, axis=0)
    dev = Bm - Q
    B = dev.T @ dev / (Mn - 1)
    T = W + (1 + 1 / Mn) * B
    k = Q.size
    wdiag = np.diag(W)
    bdiag = np.diag(B)
    tdiag = np.diag(T)
    keep = wdiag > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(keep, (1 + 1 / Mn) * bdiag / wdiag, 0.0)
        gam = np.where(tdiag > 0, (1 + 1 / Mn) * bdiag / tdiag, 0.0)
        df_old = np.where(gam > 0, (Mn - 1) / gam ** 2, np.inf)
    small = dfr[0] is not None and not o.get("nosmall")
    if small:
        nu_c = dfr[0]
        df_obs = (nu_c + 1) / (nu_c + 3) * nu_c * (1 - gam)
        df = 1 / (1 / df_old + 1 / df_obs)
    else:
        df = df_old
    fmi = np.where(keep, (r + 2 / (df + 3)) / (r + 1), 0.0)
    # model test of the slopes (equal FMI)
    sl = [i for i in range(k) if names[i] != "_cons" and keep[i] and not names[i].startswith(("o.", "b."))
          and "b." not in names[i] and eqs[i] in (eqs[0], "_")]
    if cmdname in ("mlogit",):
        sl = [i for i in range(k) if names[i] != "_cons" and keep[i] and "b." not in names[i]]
    kk = len(sl)
    F = df1 = df2 = pF = SYS
    if kk:
        Ws, Bs = W[np.ix_(sl, sl)], B[np.ix_(sl, sl)]
        rbar = (1 + 1 / Mn) * np.trace(Bs @ np.linalg.inv(Ws)) / kk
        Tt = (1 + rbar) * Ws
        F = float(Q[sl] @ np.linalg.solve(Tt, Q[sl])) / kk
        v1 = kk * (Mn - 1)
        if v1 > 4:
            df2 = 4 + (v1 - 4) * (1 + (1 - 2 / v1) / rbar) ** 2
        else:
            df2 = v1 * (1 + 1 / kk) * (1 + 1 / rbar) ** 2 / 2
        if small:
            # VERIFY: Stata's small-sample df for the model test (Reiter 2007)
            df2 = min(df2, float(dfr[0]))
        df1 = kk
        pF = float(st.f.sf(F, df1, df2))
    N = float(np.mean(Ns))
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    est = Estimates(cmdname, names and "", names, Q, T, int(N), "t", None)
    est.eqnames = eqs
    est.level = level
    est.vcetype = ""
    est.mi = {"df": df, "r": r, "fmi": fmi, "M": Mn, "small": small, "nu_c": dfr[0], "F": F, "df1": df1,
              "df2": df2, "pF": pF, "N": N, "title": title, "vce": vce, "W": W, "B": B}
    est.display = display_mi
    # the depvar comes from the command line (first word after the command name)
    est.depvar = cmd.split()[1] if len(cmd.split()) > 1 else ""
    from ..commands.matrix import Matrix
    scalars = [("N", N), ("M_mi", float(Mn)), ("df_m_mi", float(df1) if df1 < SYS else SYS),
               ("df_r_mi", float(df2) if df2 < SYS else SYS), ("F_mi", F), ("p_mi", pF),
               ("rvi_avg_mi", float(np.mean(r[keep])) if keep.any() else SYS),
               ("fmi_max_mi", float(np.max(fmi[keep])) if keep.any() else SYS)]
    if small:
        scalars.append(("df_c_mi", float(dfr[0])))
    macros = [("mi", "mi"), ("cmd_mi", cmdname), ("cmd", "mi estimate"), ("title_mi", "Multiple-imputation estimates"),
              ("title", title), ("dfadjust_mi", "Small sample" if small else "Large sample"),
              ("modeltest_mi", "Equal FMI"), ("depvar", est.depvar), ("properties", "b V")]
    mats = [("W_mi", Matrix(W, list(names), list(names))), ("B_mi", Matrix(B, list(names), list(names))),
            ("df_mi", Matrix(df.reshape(1, -1), ["df"], list(names))),
            ("fmi_mi", Matrix(fmi.reshape(1, -1), ["fmi"], list(names)))]
    post(s, est, scalars, macros, np.zeros(ds.nobs, dtype=bool), extra_mats=mats)
    display_mi(s, est)


def display_mi(s: "Session", est, *, header: bool = True, table: bool = True, level=None, **_kw) -> None:
    out = s.output
    X = est.mi
    lev = level or est.level
    if header:
        def right(left, lab, val):
            out.write(f"{left:<48}{lab:<18}= ", "text")
            out.write(f"{val:>10}\n", "result")
        keep = np.diag(X["W"]) > 0
        out.write("\n", "text")
        right("Multiple-imputation estimates", "Imputations", f"{X['M']:,}")
        right(X["title"], "Number of obs", f"{int(round(X['N'])):,}")
        right("", "Average RVI", f"{np.mean(X['r'][keep]):.4f}")
        right("", "Largest FMI", f"{np.max(X['fmi'][keep]):.4f}")
        if X["small"]:
            right("", "Complete DF", f"{int(X['nu_c']):,}")
        dfs = X["df"][keep]

        def dfv(v):
            return f"{v:.2f}" if np.isfinite(v) else "."
        left = "DF adjustment:" + f"{'Small sample' if X['small'] else 'Large sample':>15}"
        right(left, "DF:     min", dfv(np.min(dfs)))
        right("", "        avg", dfv(np.mean(dfs)))
        right("", "        max", dfv(np.max(dfs)))
        vce = {"ols": "OLS", "oim": "OIM", "robust": "Robust", "cluster": "Robust", "opg": "OPG"}.get(
            X["vce"].lower(), X["vce"].upper() if X["vce"] else "OLS")
        if X["F"] < SYS:
            right("Model F test:" + f"{'Equal FMI':>16}", f"F({int(X['df1']):>4},{X['df2']:>7.1f})",
                  f"{X['F']:.2f}")
            right("Within VCE type:" + f"{vce:>13}", "Prob > F", f"{X['pF']:.4f}")
        else:
            out.write("Within VCE type:" + f"{vce:>13}\n", "text")
        out.write("\n", "text")
    if not table:
        return
    se = np.sqrt(np.maximum(np.diag(est.V), 0))
    ci = f"[{lev:g}% Conf. Interval]"
    out.write("-" * 78 + "\n", "text")
    out.write(f"{est.depvar[:12]:>12} |      Coef.   Std. Err.      t    P>|t|     {ci:>20}\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    a = 1 - lev / 100
    eq_prev = None
    multi = len(set(est.eqnames)) > 1
    for i, nm in enumerate(est.names):
        if multi and est.eqnames[i] != eq_prev:
            if eq_prev is not None:
                out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
            out.write(f"{est.eqnames[i][:12]:<12} |\n", "result")
            eq_prev = est.eqnames[i]
        lab = nm
        out.write(f"{lab[:12]:>12} |", "text")
        b, sd, d = est.b[i], se[i], X["df"][i]
        if sd <= 0 or not np.isfinite(sd):
            out.write("          0  (omitted)\n" if "o." in nm or "b." in nm else
                      f"  {_g(b):>9}          .        .       .            .           .\n", "result")
            continue
        tv = b / sd
        dd = d if np.isfinite(d) else 1e12
        p = float(2 * st.t.sf(abs(tv), dd))
        c = float(st.t.ppf(1 - a / 2, dd))
        out.write(f"  {_g(b):>9}  {_g(sd):>9}  {tv:>7.2f}  {p:>6.3f}    {_g(b - c * sd):>9}   {_g(b + c * sd):>9}\n",
                  "result")
    out.write("-" * 78 + "\n", "text")
