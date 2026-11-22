"""bsqreg, iqreg and sqreg ([R] qreg): quantile regression with bootstrap standard errors,
interquantile-range regression and simultaneous quantile regression.

The coefficients are those of qreg (estimation/extra.py, _qfit); the variance is the covariance
of the bootstrap replicates. The draws use OpenDTA's generator, so the standard errors differ
from Stata's for the same seed (the coefficients do not)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.words import parse_numlist
from .extra import _qfit
from .resample import _dot, _dots
from .results import CoefRow, Estimates, coef_table, g9, post

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# computations
# ---------------------------------------------------------------------------

def _rho(u: np.ndarray, q: float) -> np.ndarray:
    return u * (q - (u < 0))


def deviations(X: np.ndarray, y: np.ndarray, q: float, b: np.ndarray) -> tuple[float, float, float]:
    """(raw sum of deviations about the q quantile of y, min sum of deviations, the quantile);
    the quantile is the observation of order ceil(n q), as qreg reports it."""
    n = len(y)
    ys = np.sort(y, kind="stable")
    yq = float(ys[min(max(int(np.ceil(n * q - 1e-9)) - 1, 0), n - 1)])
    raw = float(_rho(y - yq, q).sum())
    mn = float(_rho(y - X @ b, q).sum())
    return raw, mn, yq


def boot_coefs(X: np.ndarray, y: np.ndarray, qs: list[float], reps: int, uniform, on_rep=None) -> np.ndarray:
    """reps x (k * len(qs)) matrix of the coefficients refitted on bootstrap samples (rows of
    NaN when a fit fails)."""
    n, k = X.shape
    out = np.full((reps, k * len(qs)), np.nan)
    for r in range(reps):
        idx = np.minimum((uniform(n) * n).astype(np.int64), n - 1)
        ok = True
        try:
            out[r] = np.concatenate([_qfit(X[idx], y[idx], q) for q in qs])
        except StataError:
            ok = False
        if on_rep is not None:
            on_rep(ok, r)
    return out


def boot_vce(R: np.ndarray) -> np.ndarray:
    R = R[~np.isnan(R).any(axis=1)]
    if len(R) < 2:
        return np.full((R.shape[1], R.shape[1]), np.nan)
    return np.atleast_2d(np.cov(R, rowvar=False, ddof=1))


# ---------------------------------------------------------------------------
# common setup
# ---------------------------------------------------------------------------

def _quantile(text: str) -> float:
    q = float(text)
    if q >= 1:
        q /= 100
    if not 0 < q < 1:
        raise StataError(198, "quantiles must be between 0 and 1")   # VERIFY
    return q


def _setup(s: "Session", cmd: str, args: str, extra: dict):
    from .postest import replay
    from .regress import Sample
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return None
    spec = {"reps": 1, "level": 1, "nodots": 5, "wlsiter": 3, **extra}
    smp = Sample(s, args, opts_spec=spec, weights=())
    from . import fvars
    fvars.drop_collinear(smp.cols, smp.constant)
    cols = [c for c in smp.cols if c.values is not None]
    y = smp.y.astype(float)
    n = len(y)
    X = np.column_stack([c.values for c in cols] + ([np.ones(n)] if smp.constant else []))
    names = [c.name for c in cols] + (["_cons"] if smp.constant else [])
    reps = int(smp.o["reps"]) if smp.o.get("reps") else 20
    if reps < 2:
        raise StataError(198, "reps() must be an integer greater than 1")   # VERIFY
    return smp, X, y, names, reps


def _bootstrap(s: "Session", X, y, qs, reps, o) -> np.ndarray:
    from ..core import rng
    if not o.get("nodots"):
        s.output.write("(fitting base model)\n", "text")   # not with nodots (compat 2008)
    show = _dots(s, "Bootstrap replications", reps, o)
    cb = (lambda ok, r: _dot(s, ok, r, reps)) if show else None
    return boot_coefs(X, y, qs, reps, rng.RNG.uniform, cb)


def _eqnames(qs: list[float]) -> list[str]:
    """q and the quantile with as many decimals as the most precise one needs, at least two
    (compat 2008: .25 .5 .75 give q25 q50 q75; .1 .333 give q100 q333)."""
    d = 2
    while d < 8 and any(abs(q * 10 ** d - round(q * 10 ** d)) > 1e-9 for q in qs):
        d += 1
    return [f"q{round(q * 10 ** d)}" for q in qs]


def _qlabel(q: float) -> str:
    return format_value(q, "%9.0g", pad=False).strip()


def _rows(names: list[str], b: np.ndarray, V: np.ndarray, eq: str = "") -> list[CoefRow]:
    se = np.sqrt(np.maximum(np.diag(V), 0))
    return [CoefRow("coef", n, float(bi), float(si), eq=eq, name=n) for n, bi, si in zip(names, b, se)]


def _macros(cmd: str, args: str, depvar: str, seed: str, extra=()) -> list[tuple[str, str]]:
    """e() macros in the order ereturn list shows them (compat 2008); e(seed) is the state of
    the random-number generator before the replications."""
    return [("cmd", cmd), ("cmdline", f"{cmd} {args.strip()}"),
            ("predict", "sqreg_p" if cmd == "sqreg" else "qreg_p"),
            ("marginsnotok", "stdp stddp Residuals"),
            *([("vcetype", "Bootstrap")] if cmd != "bsqreg" else []), *extra,
            ("depvar", depvar), ("seed", seed), ("properties", "b V")]


def _seed(s: "Session") -> str:
    return str(s.creturn("rngstate"))


# ---------------------------------------------------------------------------
# bsqreg
# ---------------------------------------------------------------------------

def cmd_bsqreg(s: "Session", args: str) -> None:
    """bsqreg depvar indepvars [, quantile(#) reps(#) level(#)]. VERIFY: the header and the
    e() results."""
    setup = _setup(s, "bsqreg", args, {"quantile": 1})
    if setup is None:
        return
    smp, X, y, names, reps = setup
    q = _quantile(smp.o["quantile"]) if smp.o.get("quantile") else 0.5
    b = _qfit(X, y, q)
    raw, mn, yq = deviations(X, y, q, b)
    seed = _seed(s)
    V = boot_vce(_bootstrap(s, X, y, [q], reps, smp.o))
    n, k = X.shape
    est = Estimates("bsqreg", smp.depname, names, b, V, n, "t", float(n - k))
    est.rows = _rows(names, b, V)
    est.level = smp.level
    est.extra = {"q": q, "raw": raw, "min": mn, "yq": yq, "r2_p": 1 - mn / raw if raw > 0 else SYS,
                 "N": n, "reps": reps}
    est.display = _display_bsqreg
    scalars = [("reps", float(reps)), ("N", float(n)), ("df_r", float(n - k)), ("q_v", yq),
               ("sum_rdev", raw), ("q", q), ("sum_adev", mn), ("convcode", 0.0), ("rank", float(k))]
    post(s, est, scalars, _macros("bsqreg", args, smp.depname, seed), smp.mask)
    _display_bsqreg(s, est)


def _display_bsqreg(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        q = x["q"]
        title = "Median regression" if abs(q - 0.5) < 1e-12 else f"{_qlabel(q)} Quantile regression"
        title += f", bootstrap({x['reps']}) SEs"
        out.write(f"\n{title:<52}Number of obs = ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        out.write(f"  Raw sum of deviations {g9(x['raw'], '%9.0g'):>8} (about {g9(x['yq'], '%10.0g').strip()})\n", "text")
        out.write(f"  Min sum of deviations {g9(x['min'], '%9.0g'):>8}{'':<20}Pseudo R2     = ", "text")
        out.write(f"{x['r2_p']:>10.4f}\n\n", "result")
    if table:
        # no Bootstrap heading over the table (compat 2008)
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level)


# ---------------------------------------------------------------------------
# iqreg
# ---------------------------------------------------------------------------

def cmd_iqreg(s: "Session", args: str) -> None:
    """iqreg depvar indepvars [, quantiles(# #) reps(#) level(#) nodots]: the difference of the
    coefficients at two quantiles. VERIFY: the header and the e() results."""
    setup = _setup(s, "iqreg", args, {"quantiles": 1})
    if setup is None:
        return
    smp, X, y, names, reps = setup
    qs = [_quantile(str(v)) for v in parse_numlist(str(smp.o["quantiles"]))] if smp.o.get("quantiles") \
        else [0.25, 0.75]
    if len(qs) != 2:
        raise StataError(198, "quantiles() must contain exactly two values")   # VERIFY
    lo, hi = sorted(qs)
    b0, b1 = _qfit(X, y, lo), _qfit(X, y, hi)
    b = b1 - b0
    k = X.shape[1]
    seed = _seed(s)
    R = _bootstrap(s, X, y, [lo, hi], reps, smp.o)
    V = boot_vce(R[:, k:] - R[:, :k])
    n = len(y)
    r2 = []
    for q, bq in ((lo, b0), (hi, b1)):
        raw, mn, _ = deviations(X, y, q, bq)
        r2.append(1 - mn / raw if raw > 0 else SYS)
    est = Estimates("iqreg", smp.depname, names, b, V, n, "t", float(n - k))
    est.rows = _rows(names, b, V)
    est.level = smp.level
    est.extra = {"q0": lo, "q1": hi, "r2": r2, "N": n, "reps": reps}
    est.display = _display_iqreg
    # compat 2008: the sums with index 0 are those of the upper quantile
    d_hi, d_lo = deviations(X, y, hi, b1), deviations(X, y, lo, b0)
    scalars = [("reps", float(reps)), ("N", float(n)), ("df_r", float(n - k)), ("q0", lo), ("q1", hi),
               ("sumrdev0", d_hi[0]), ("sumrdev1", d_lo[0]), ("sumadev0", d_hi[1]),
               ("sumadev1", d_lo[1]), ("convcode", 0.0), ("rank", float(k))]
    post(s, est, scalars, _macros("iqreg", args, smp.depname, seed), smp.mask)
    _display_iqreg(s, est)


def _display_iqreg(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        lo, hi = _qlabel(x["q0"]), _qlabel(x["q1"])
        out.write(f"\n{hi + '-' + lo + ' Interquantile regression':<52}Number of obs = ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        out.write(f"{'  bootstrap(' + str(x['reps']) + ') SEs':<52}{hi + ' Pseudo R2':>13} = ", "text")
        out.write(f"{x['r2'][1]:>10.4f}\n", "result")
        out.write(f"{'':<52}{lo + ' Pseudo R2':>13} = ", "text")
        out.write(f"{x['r2'][0]:>10.4f}\n\n", "result")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level,
                   vcetype="Bootstrap")


# ---------------------------------------------------------------------------
# sqreg
# ---------------------------------------------------------------------------

def cmd_sqreg(s: "Session", args: str) -> None:
    """sqreg depvar indepvars [, quantiles(# ...) reps(#) level(#) nodots]: one equation per
    quantile (q25, q50...), with the bootstrap covariance across equations. VERIFY: the header,
    the equation names of quantiles that are not whole percents and the e() results."""
    setup = _setup(s, "sqreg", args, {"quantiles": 1})
    if setup is None:
        return
    smp, X, y, names, reps = setup
    qs = sorted({_quantile(str(v)) for v in parse_numlist(str(smp.o["quantiles"]))}) \
        if smp.o.get("quantiles") else [0.5]
    bs = [_qfit(X, y, q) for q in qs]
    b = np.concatenate(bs)
    seed = _seed(s)
    V = boot_vce(_bootstrap(s, X, y, qs, reps, smp.o))
    n, k = X.shape
    eqs = _eqnames(qs)
    rows: list[CoefRow] = []
    r2 = []
    for j, (q, eq, bq) in enumerate(zip(qs, eqs, bs)):
        if j:
            rows.append(CoefRow("sep", ""))
        rows.append(CoefRow("eq", eq))
        sl = slice(j * k, (j + 1) * k)
        rows += _rows(names, bq, V[sl, sl], eq=eq)
        raw, mn, _ = deviations(X, y, q, bq)
        r2.append(1 - mn / raw if raw > 0 else SYS)
    full = [nm for _ in eqs for nm in names]
    est = Estimates("sqreg", smp.depname, full, b, V, n, "t", float(n - k))
    est.eqnames = [eq for eq in eqs for _ in names]
    est.rows = rows
    est.level = smp.level
    est.extra = {"qs": qs, "r2": r2, "N": n, "reps": reps}
    est.display = _display_sqreg
    scalars = [("reps", float(reps)), ("N", float(n)), ("df_r", float(n - k))]
    for j, q in enumerate(qs, 1):
        raw, mn, _ = deviations(X, y, q, bs[j - 1])
        scalars += [(f"q{j}", q), (f"sumrdv{j}", raw), (f"sumadv{j}", mn)]
    scalars += [("n_q", float(len(qs))), ("convcode", 0.0), ("rank", float(k * len(qs)))]
    post(s, est, scalars, _macros("sqreg", args, smp.depname, seed, [("eqnames", " ".join(eqs))]),
         smp.mask)
    _display_sqreg(s, est)


def _display_sqreg(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        out.write(f"\n{'Simultaneous quantile regression':<52}Number of obs = ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        for j, (q, r2) in enumerate(zip(x["qs"], x["r2"])):
            left = f"  bootstrap({x['reps']}) SEs" if j == 0 else ""
            lab = f"{q:.2f}".lstrip("0")   # compat 2008: .10, .33, .50
            out.write(f"{left:<52}{lab + ' Pseudo R2':>13} = ", "text")
            out.write(f"{r2:>10.4f}\n", "result")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level,
                   vcetype="Bootstrap", eq_headers=True)
