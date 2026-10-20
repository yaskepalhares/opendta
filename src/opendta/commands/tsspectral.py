"""Periodogram, spectral and white-noise tools: pergram, cumsp, xcorr, wntestb, wntestq and
psdensity ([TS] pergram, [TS] cumsp, [TS] xcorr, [TS] wntestb, [TS] wntestq, [TS] psdensity).

The spectral formulas are in the manual: omega_k = (k-1)/n, C_k^2 = |sum x(t) e^(2 pi i (t-1)
omega_k)|^2 / n^2 and the sample spectral density f(omega_k) = n C_k^2. The layouts the manual
does not show are marked VERIFY (compat/do/1904_timeseries.do)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..estimation.tsops import require_ts
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from ._util import touse
from .registry import command
from .tshelp import kv, new_variable, panels, to_nan

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _one_series(s: "Session", p, *, dense: bool = True):
    """The single variable of the command over the sample, in time order, with its times."""
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) != 1:
        raise StataError(103 if names else 102,
                         "too many variables specified" if names else "too few variables specified")
    return names[0]


def _values(s: "Session", tok: str) -> np.ndarray:
    from .tsfilters import _series
    return _series(s, tok)


def _sample_series(s: "Session", p, tok: str, *, what: str = "") -> np.ndarray:
    """Values of tok over the sample in time order; the sample must be dense."""
    info = require_ts(s.data)
    x = _values(s, tok)
    mask = touse(s, p)
    groups = panels(s, mask)
    if info.panel and len(groups) > 1:
        raise StataError(498, "panel data not allowed; use a single time series")      # VERIFY
    _pv, rows, t = groups[0]
    if not len(rows):
        raise StataError(2000, "no observations")
    if np.isnan(x[rows]).any() or np.any(np.abs(np.diff(t) - info.delta) > 1e-9):
        raise StataError(498, "sample may not contain gaps" if not np.isnan(x[rows]).any()
                         else "sample may not contain missing values")      # VERIFY
    return x[rows], rows


# ---------------------------------------------------------------------------
# periodogram
# ---------------------------------------------------------------------------

def periodogram(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """(frequencies, n C_k^2, sample variance) at the natural frequencies k = 1..n/2 + 1; the
    series is demeaned."""
    n = len(x)
    xc = x - x.mean()
    f = np.fft.rfft(xc)
    c2 = (np.abs(f) ** 2) / n ** 2
    om = np.arange(len(f)) / n
    return om, n * c2, float((xc @ xc) / n)


def _tsplot_graph(s: "Session", p, graw: dict, layers, xtitle: str, ytitle: str, title: str,
                  note: str | None = None):
    from ..graphics.baseplots import _graph
    return _graph(s, graw, layers, xtitle, ytitle, True, title, note)


@command("pergram", "pergram")
def cmd_pergram(s: "Session", args: str) -> None:
    from ..graphics.baseplots import _own, _split
    from ..graphics.model import Layer
    from .graphcmd import finish
    p, lopts, graw = _split(s, args, {"generate": 1, "nograph": 4})
    tok = _one_series(s, p)
    x, rows = _sample_series(s, p, tok)
    om, pg, var = periodogram(x)
    if isinstance(lopts.get("generate"), str):
        col = np.full(s.data.nobs, np.nan)
        col[:len(pg)] = pg
        new_variable(s, lopts["generate"], col, vtype="double")
    if lopts.get("nograph"):
        return
    with np.errstate(divide="ignore"):
        y = np.log(pg / var)
    y = np.clip(y, -6, 6)
    keep = np.arange(len(om)) >= 1
    lab = s.data.get(tok).label if s.data.has(tok) else ""
    layers = [Layer("connected", om[keep].tolist(), y[keep].tolist(), {}, _own(lopts, {"generate", "nograph"}), "")]
    g = _tsplot_graph(s, p, graw, layers, "Frequency", "Log Periodogram", lab or "Sample spectral density function",
                      "Evaluated at the natural frequencies")
    g.opts["_ylim"] = (-6, 6)
    finish(s, g, graw)


@command("cumsp", "cumsp")
def cmd_cumsp(s: "Session", args: str) -> None:
    from ..graphics.baseplots import _own, _split
    from ..graphics.model import Layer
    from .graphcmd import finish
    p, lopts, graw = _split(s, args, {"generate": 1})
    tok = _one_series(s, p)
    x, rows = _sample_series(s, p, tok)
    om, pg, _v = periodogram(x)
    F = np.cumsum(pg) / pg.sum()
    if isinstance(lopts.get("generate"), str):
        col = np.full(s.data.nobs, np.nan)
        col[:len(F)] = F
        new_variable(s, lopts["generate"], col, vtype="double")
    layers = [Layer("connected", om.tolist(), F.tolist(), {}, _own(lopts, {"generate"}), "")]
    g = _tsplot_graph(s, p, graw, layers, "Frequency", "Cumulative spectral distribution",
                      "Sample spectral distribution function", "Points evaluated at the natural frequencies")
    finish(s, g, graw)


# ---------------------------------------------------------------------------
# xcorr
# ---------------------------------------------------------------------------

def crosscorr(x1: np.ndarray, x2: np.ndarray, Q: int) -> np.ndarray:
    """rho_12(k) = Cov[x1(t), x2(t+k)] / sqrt(R11(0) R22(0)), k = -Q..Q."""
    n = len(x1)
    a, b = x1 - x1.mean(), x2 - x2.mean()
    den = math.sqrt(float(a @ a) * float(b @ b))
    out = []
    for k in range(-Q, Q + 1):
        if k >= 0:
            out.append(float(a[:n - k] @ b[k:]) / den)
        else:
            out.append(float(a[-k:] @ b[:n + k]) / den)
    return np.array(out)


def _bar(v: float) -> str:
    k = min(int(math.floor(abs(v) * 8 + 1e-12)), 8)
    return " " * 8 + "|" + "-" * k if v >= 0 else " " * (8 - k) + "-" * k + "|"


@command("xcorr", "xc")
def cmd_xcorr(s: "Session", args: str) -> None:
    from ..graphics.baseplots import _own, _split
    from ..graphics.model import Layer
    from .graphcmd import finish
    p, lopts, graw = _split(s, args, {"generate": 1, "lags": 1, "table": 1, "noplot": 3, "base": 4})
    names = p.varlist.split()
    if len(names) != 2:
        raise StataError(102 if len(names) < 2 else 103,
                         "too few variables specified" if len(names) < 2 else "too many variables specified")
    x1, rows = _sample_series(s, p, names[0])
    x2, _r = _sample_series(s, p, names[1])
    n = len(x1)
    Q = int(float(s.eval(str(lopts["lags"])))) if isinstance(lopts.get("lags"), str) else min(n // 2 - 2, 20)
    if Q < 1 or Q >= n:
        raise StataError(198, "lags() out of range")      # VERIFY
    r = crosscorr(x1, x2, Q)
    if isinstance(lopts.get("generate"), str):
        col = np.full(s.data.nobs, np.nan)
        col[:len(r)] = r
        new_variable(s, lopts["generate"], col, vtype="double")
    lags = np.arange(-Q, Q + 1)
    if lopts.get("table"):
        out = s.output
        noplot = bool(lopts.get("noplot"))
        out.write("\n", "text")
        if not noplot:
            out.write(f"{'':28}-1       0       1\n", "text")
            out.write(" LAG       CORR     [Cross-correlation]\n", "text")
        else:
            out.write(" LAG      CORR\n", "text")
        out.write("-" * (47 if not noplot else 15) + "\n", "text")
        for k, v in zip(lags, r):
            if noplot:
                out.write(f"{k:<8}{v:>7.4f}   \n", "result")      # VERIFY: Stata 14 layout
                continue
            line = f"{k:<5}{v:>10.4f}"
            line += f"  {_bar(v)}"
            out.write(line.rstrip() + "\n", "result")
        return
    layers = [Layer("spike", lags.tolist(), r.tolist(), {}, _own(lopts, {"generate", "lags", "table", "noplot", "base"}), "")]
    g = _tsplot_graph(s, p, graw, layers, "Lag", "Cross-correlogram", "")
    finish(s, g, graw)


# ---------------------------------------------------------------------------
# wntestb, wntestq
# ---------------------------------------------------------------------------

def kolmogorov_sf(a: float) -> float:
    """1 - G(a), G(a) = sum_{j=-inf}^{inf} (-1)^j exp(-2 a^2 j^2)."""
    if a <= 0:
        return 1.0
    j = np.arange(1, 200)
    v = 2 * float(np.sum((-1.0) ** (j - 1) * np.exp(-2 * a * a * j * j)))
    return min(max(v, 0.0), 1.0)


def bartlett_b(x: np.ndarray):
    """Bartlett's B statistic, its p-value and the cumulative periodogram."""
    n = len(x)
    om, pg, _v = periodogram(x)
    q = n // 2 + 1
    F = np.cumsum(pg[:q]) / pg[:q].sum()
    k = np.arange(1, q + 1)
    dev = np.abs(F - k / q)
    B = math.sqrt(n / 2) * float(dev.max())
    return B, kolmogorov_sf(B), om[:q], F, q


@command("wntestb", "wntestb")
def cmd_wntestb(s: "Session", args: str) -> None:
    from ..graphics.baseplots import _own, _split
    from ..graphics.model import Layer
    from .graphcmd import finish
    p, lopts, graw = _split(s, args, {"table": 1, "level": 1, "noplot": 3})
    tok = _one_series(s, p)
    x, _rows = _sample_series(s, p, tok)
    B, pv, om, F, q = bartlett_b(x)
    kv(s, [("stat", B), ("p", pv)])
    out = s.output
    if lopts.get("table"):
        out.write("\nCumulative periodogram white-noise test\n", "text")
        out.write("-" * 39 + "\n", "text")
        out.write(f" Bartlett's (B) statistic  = {B:>10.4f}\n", "text")
        out.write(f" Prob > B".ljust(27) + f"= {pv:>10.4f}\n", "text")
        return
    if lopts.get("noplot"):
        return
    n = len(x)
    level = float(s.eval(str(lopts["level"]))) if isinstance(lopts.get("level"), str) else float(s.settings.get("level", "95"))
    lo_, hi_ = 0.0, 5.0
    target = 1 - level / 100
    for _ in range(80):
        mid = (lo_ + hi_) / 2
        if kolmogorov_sf(mid) > target:
            lo_ = mid
        else:
            hi_ = mid
    a = (lo_ + hi_) / 2 / math.sqrt(n / 2)
    k = np.arange(1, q + 1) / q
    layers = [Layer("connected", om.tolist(), F.tolist(), {}, _own(lopts, {"table", "level", "noplot"}), ""),
              Layer("line", om.tolist(), k.tolist(), {}, {"lcolor": "gs8"}, ""),
              Layer("line", om.tolist(), (k + a).tolist(), {}, {"lcolor": "gs8", "lpattern": "dash"}, ""),
              Layer("line", om.tolist(), (k - a).tolist(), {}, {"lcolor": "gs8", "lpattern": "dash"}, "")]
    g = _tsplot_graph(s, p, graw, layers, "Frequency", "Cumulative periodogram",
                      f"Cumulative periodogram for {tok}",
                      f"Bartlett's (B) statistic = {B:.2f}  Prob > B = {pv:.4f}")
    finish(s, g, graw)


@command("wntestq", "wntestq")
def cmd_wntestq(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"lags": 1}) if p.options.strip() else {}
    tok = _one_series(s, p)
    x, _rows = _sample_series(s, p, tok)
    n = len(x)
    m = int(float(s.eval(str(o["lags"])))) if isinstance(o.get("lags"), str) else min(n // 2 - 2, 40)
    if m < 1 or m >= n:
        raise StataError(198, "lags() out of range")      # VERIFY
    xc = x - x.mean()
    c0 = float(xc @ xc)
    ac = np.array([float(xc[k:] @ xc[:-k]) / c0 for k in range(1, m + 1)])
    Q = n * (n + 2) * float(np.sum(ac ** 2 / (n - np.arange(1, m + 1))))
    pv = float(chi2_tail(Q, m))
    out = s.output
    out.write("\nPortmanteau test for white noise\n", "text")
    out.write("-" * 39 + "\n", "text")
    out.write(f" Portmanteau (Q) statistic = {Q:>10.4f}\n", "text")
    out.write(f" Prob > chi2({m})".ljust(27) + f"= {pv:>10.4f}\n", "text")
    kv(s, [("stat", Q), ("df", float(m)), ("p", pv)])


# ---------------------------------------------------------------------------
# psdensity
# ---------------------------------------------------------------------------

def arma_density(ar: dict, ma: dict, sigma2: float, omega: np.ndarray, d: float = 0.0):
    """sigma^2 / (2 pi) |theta(e^-iw)|^2 / |phi(e^-iw)|^2 |1 - e^iw|^(-2d); ar and ma map the lag
    to its coefficient."""
    z = np.exp(-1j * omega)
    phi = 1 - sum(c * z ** k for k, c in ar.items()) if ar else np.ones_like(z)
    theta = 1 + sum(c * z ** k for k, c in ma.items()) if ma else np.ones_like(z)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = sigma2 / (2 * math.pi) * np.abs(theta) ** 2 / np.abs(phi) ** 2
        if d:
            f = f * np.abs(1 - np.exp(1j * omega)) ** (-2 * d)
    return f


def arma_variance(ar: dict, ma: dict, sigma2: float, d: float = 0.0) -> float:
    """Variance of the ARMA process (psi weights summed until they vanish)."""
    from scipy import signal
    P = max(ar) if ar else 0
    Q = max(ma) if ma else 0
    a = np.zeros(P + 1)
    a[0] = 1
    for k, c in ar.items():
        a[k] = -c
    b = np.zeros(Q + 1)
    b[0] = 1
    for k, c in ma.items():
        b[k] = c
    imp = np.zeros(20000)
    imp[0] = 1
    psi = signal.lfilter(b, a, imp)
    return float(sigma2 * np.sum(psi ** 2))


@command("psdensity", "psdensity")
def cmd_psdensity(s: "Session", args: str) -> None:
    from ..estimation.results import current
    p = parse_standard(args)
    o = match_options(p.options, {"pspectrum": 2, "range": 3, "cycle": 3, "smemory": 2}) if p.options.strip() else {}
    names = p.varlist.split()
    if len(names) != 2:
        raise StataError(102 if len(names) < 2 else 103,
                         "too few variables specified" if len(names) < 2 else "too many variables specified")
    est = current(s)
    sample = touse(s, p)
    rows = np.flatnonzero(sample)
    N = len(rows)
    if est.cmd == "arima":
        A = est.arima
        b = est.b
        k = A["k"]
        na, nm = len(A["ar"]), len(A["ma"])
        ar = {lag: float(b[k + i]) for i, lag in enumerate(A["ar"])}
        ma = {lag: float(b[k + na + i]) for i, lag in enumerate(A["ma"])}
        sigma2 = float(b[-1]) ** 2
        d = 0.0
    elif est.cmd == "arfima":
        A = est.arfima
        ar, ma, sigma2, d = A["ar"], A["ma"], A["sigma2"], (0.0 if o.get("smemory") else A["d"])
    else:
        raise StataError(321, "psdensity only works after arima, arfima and ucm")      # VERIFY
    lo, hi = 0.0, math.pi
    if isinstance(o.get("range"), str):
        v = [float(w) for w in o["range"].split()]
        if len(v) != 2 or not 0 <= v[0] < v[1] < math.pi:
            raise StataError(198, "range() invalid")
        lo, hi = v
    omega = math.pi * np.arange(N) / (N - 1) if N > 1 else np.array([0.0])
    omega = omega[(omega >= lo) & (omega < hi)] if isinstance(o.get("range"), str) else omega
    f = arma_density(ar, ma, sigma2, omega, d)
    if not o.get("pspectrum"):
        f = f / arma_variance(ar, ma, sigma2, d) if d == 0 else f / A["gamma0"]
    for nm in names:
        if s.data.has(nm):
            raise StataError(110, f"{nm} already defined")
    for nm, vals in ((names[0], f), (names[1], omega)):
        col = np.full(s.data.nobs, np.nan)
        col[rows[:len(vals)]] = vals
        new_variable(s, nm, col)
