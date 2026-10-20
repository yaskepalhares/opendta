"""Graph commands of the base-statistics block: quantile, symplot, serrbar, dotplot, sunflower,
lpoly, cusum and gladder ([R] quantile, [R] diagnostic plots, [G-2] graph twoway, [R] dotplot,
[R] sunflower, [R] lpoly, [R] cusum, [R] ladder).

The drawing is OpenDTA's engine; the default titles and axes follow the manuals and the details
are marked VERIFY (compat/do/1903_base_stats.do)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from ..lang.syntax import parse_standard
from ..lazy import lazy
from . import data as D
from .model import Graph, Layer, Panel
from .opts import graph_options, split_layer_graph
from .twoway import layer_spec

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session


def _split(s: "Session", text: str, extra: dict):
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, {**layer_spec(), **extra})
    return p, lopts, graw


def _own(lopts: dict, extra: dict) -> dict:
    """The layer options that are styles (not this command's own)."""
    return {k: v for k, v in lopts.items() if k not in extra}


def _vars(s: "Session", p, n: int) -> list[str]:
    names = s.expand_varlist(p.varlist)
    if len(names) < n:
        raise StataError(102, "too few variables specified")
    if len(names) > n:
        raise StataError(103, "too many variables specified")
    return names


def _graph(s, graw, layers, xtitle, ytitle, legend_off=True, title=None, note=None):
    opts = graph_options(graw)
    opts["_xtitle"], opts["_ytitle"] = xtitle, ytitle
    if legend_off:
        opts.setdefault("legend", {"off": True})
    if title:
        opts.setdefault("title", {"lines": [title]})
    if note:
        opts.setdefault("note", {"lines": [note]})
    return Graph(family="twoway", opts=opts, panels=[Panel("", layers)])


# ---------------------------------------------------------------------------
# quantile, symplot
# ---------------------------------------------------------------------------

def quantile(s: "Session", text: str):
    """Ordered values against the fraction of the data, (i - .5)/N ([R] diagnostic plots)."""
    p, lopts, graw = _split(s, text, {"grid": 4})
    v = _vars(s, p, 1)[0]
    x = D.numeric(s, v)
    mask = D.base_mask(s, p) & D.valid(x)
    xs = np.sort(x[mask])
    if len(xs) < 1:
        raise StataError(2000, "no observations")
    n = len(xs)
    frac = (np.arange(1, n + 1) - 0.5) / n
    layers = [Layer("scatter", frac.tolist(), xs.tolist(), {}, _own(lopts, {"grid"}), D.var_title(s, v))]
    return _graph(s, graw, layers, "Fraction of the data", D.var_title(s, v)), graw


def symplot(s: "Session", text: str):
    """Distance above the median against distance below, for the pairs of ordered values
    equidistant from the median ([R] diagnostic plots)."""
    p, lopts, graw = _split(s, text, {"grid": 4})
    v = _vars(s, p, 1)[0]
    x = D.numeric(s, v)
    mask = D.base_mask(s, p) & D.valid(x)
    xs = np.sort(x[mask])
    n = len(xs)
    if n < 2:
        raise StataError(2000, "no observations")
    med = float(np.median(xs))
    k = n // 2
    below = med - xs[:k]
    above = xs[::-1][:k] - med
    hi = float(max(below.max(), above.max()))
    layers = [Layer("scatter", below.tolist(), above.tolist(), {}, _own(lopts, {"grid"}), ""),
              Layer("line", [0.0, hi], [0.0, hi], {}, {"lcolor": "gs8"}, "")]
    g = _graph(s, graw, layers, "Distance below median", "Distance above median")
    return g, graw


# ---------------------------------------------------------------------------
# serrbar
# ---------------------------------------------------------------------------

def serrbar(s: "Session", text: str):
    """mvar +- scale*svar against xvar, with caps ([R] serrbar)."""
    p, lopts, graw = _split(s, text, {"scale": 5, "rcapopts": 4})
    mv, sv, xv = _vars(s, p, 3)
    m, sd, x = D.numeric(s, mv), D.numeric(s, sv), D.numeric(s, xv)
    mask = D.base_mask(s, p) & D.valid(m) & D.valid(sd) & D.valid(x)
    scale = float(lopts["scale"]) if isinstance(lopts.get("scale"), str) else 1.0
    o = np.argsort(x[mask], kind="stable")
    xx, mm, ss = x[mask][o], m[mask][o], sd[mask][o]
    style = _own(lopts, {"scale", "rcapopts"})
    layers = [Layer("rcap", xx.tolist(), (mm - scale * ss).tolist(),
                    {"y2": (mm + scale * ss).tolist()}, {}, ""),
              Layer("scatter", xx.tolist(), mm.tolist(), {}, style, D.var_title(s, mv))]
    return _graph(s, graw, layers, D.var_title(s, xv), D.var_title(s, mv)), graw


# ---------------------------------------------------------------------------
# dotplot
# ---------------------------------------------------------------------------

def dotplot(s: "Session", text: str):
    """Dot plot of one or more variables, or of one variable over the groups of over().
    The vertical axis is binned into ny bins and the dots of a bin are stacked horizontally
    (VERIFY the default ny and the spacing against Stata)."""
    p, lopts, graw = _split(s, text, {"over": 2, "by": 2, "nx": 2, "ny": 2, "center": 2,
                                      "nogroup": 3, "nogroupvars": 3, "bar": 3, "mean": 4,
                                      "median": 3, "xlabel": 2})
    names = s.expand_varlist(p.varlist)
    if not names:
        raise StataError(102, "too few variables specified")
    base = D.base_mask(s, p)
    over = lopts.get("over") or lopts.get("by")
    cols: list[tuple[str, np.ndarray]] = []
    if over and isinstance(over, str):
        if len(names) != 1:
            raise StataError(198, "over() requires a single variable")
        y = D.numeric(s, names[0])
        ov = D.numeric(s, s.expand_varlist(over)[0])
        m = base & D.valid(y) & D.valid(ov)
        for val in sorted(set(ov[m].tolist())):
            cols.append((D.value_text(s, over, val), y[m & (ov == val)]))
    else:
        for nm in names:
            y = D.numeric(s, nm)
            cols.append((nm, y[base & D.valid(y)]))
    allv = np.concatenate([c[1] for c in cols])
    if len(allv) == 0:
        raise StataError(2000, "no observations")
    ny = int(lopts["ny"]) if isinstance(lopts.get("ny"), str) else max(5, min(35, int(round(math.sqrt(len(allv)) * 2))))
    lo, hi = float(allv.min()), float(allv.max())
    step = (hi - lo) / ny if hi > lo else 1.0
    xs, ys = [], []
    for j, (_lab, y) in enumerate(cols, start=1):
        b = np.floor((y - lo) / step + 1e-9)
        b = np.minimum(b, ny)             # the maximum goes in the last bin
        yy = lo + b * step
        for val in np.unique(yy):
            cnt = int((yy == val).sum())
            off = (np.arange(cnt) - (cnt - 1) / 2) * (0.1 if len(cols) > 1 else 1.0)
            if len(cols) == 1:
                xs += (np.arange(1, cnt + 1)).tolist()
            else:
                xs += (j + off).tolist()
            ys += [float(val)] * cnt
    layers = [Layer("scatter", xs, ys, {}, _own(lopts, {"over", "by", "nx", "ny", "center", "nogroup",
                                                         "nogroupvars", "bar", "mean", "median", "xlabel"}), "")]
    if len(cols) == 1:
        g = _graph(s, graw, layers, "Frequency", D.var_title(s, names[0]))
    else:
        g = _graph(s, graw, layers, "", "")
        g.opts["xlabel"] = {"values": list(range(1, len(cols) + 1)), "labels": [c[0] for c in cols]}
    return g, graw


# ---------------------------------------------------------------------------
# sunflower
# ---------------------------------------------------------------------------

def sunflower(s: "Session", text: str):
    """Density-distribution sunflower plot: the plane is cut in square bins; bins with fewer
    than `light` observations show the observations, with at least `dark` a dark flower,
    otherwise a light one with one petal per `petalweight` observations (VERIFY the drawing)."""
    p, lopts, graw = _split(s, text, {"binwidth": 4, "light": 3, "dark": 3, "petalweight": 6,
                                      "flowersonly": 4, "nograph": 4, "notable": 4, "marker": 3})
    yv, xv = _vars(s, p, 2)
    y, x = D.numeric(s, yv), D.numeric(s, xv)
    mask = D.base_mask(s, p) & D.valid(x) & D.valid(y)
    xx, yy = x[mask], y[mask]
    n = len(xx)
    if n == 0:
        raise StataError(2000, "no observations")
    bw = float(lopts["binwidth"]) if isinstance(lopts.get("binwidth"), str) else \
        max((xx.max() - xx.min()), (yy.max() - yy.min()), 1e-9) / max(8.0, math.sqrt(n) / 1.0)
    light = int(lopts["light"]) if isinstance(lopts.get("light"), str) else 3
    dark = int(lopts["dark"]) if isinstance(lopts.get("dark"), str) else 13
    pw = float(lopts["petalweight"]) if isinstance(lopts.get("petalweight"), str) else 1.0
    ix = np.floor((xx - xx.min()) / bw).astype(int)
    iy = np.floor((yy - yy.min()) / bw).astype(int)
    key = ix * (iy.max() + 1) + iy
    cells: dict[int, list[int]] = {}
    for i, kk in enumerate(key):
        cells.setdefault(int(kk), []).append(i)
    single_x, single_y, lx, ly, dx, dy = [], [], [], [], [], []
    for idx in cells.values():
        c = len(idx)
        cx, cy = float(xx[idx].mean()), float(yy[idx].mean())
        if c < light:
            if not lopts.get("flowersonly"):
                single_x += xx[idx].tolist()
                single_y += yy[idx].tolist()
        elif c >= dark:
            dx.append(cx)
            dy.append(cy)
        else:
            lx.append(cx)
            ly.append(cy)
    layers = []
    if single_x:
        layers.append(Layer("scatter", single_x, single_y, {}, {"msymbol": "oh"}, ""))
    if lx:
        layers.append(Layer("scatter", lx, ly, {}, {"msymbol": "T", "mcolor": "gs10"}, ""))
    if dx:
        layers.append(Layer("scatter", dx, dy, {}, {"msymbol": "T", "mcolor": "black"}, ""))
    g = _graph(s, graw, layers, D.var_title(s, xv), D.var_title(s, yv),
               note=f"bin width = {bw:g}; bin height = {bw:g}")
    if lopts.get("nograph"):
        return None, graw
    return g, graw


# ---------------------------------------------------------------------------
# lpoly (graph) -- the smoother is lpoly_fit, shared with the command-line output
# ---------------------------------------------------------------------------

_KERNELS = {
    "epanechnikov": lambda u: np.where(np.abs(u) < 1, 0.75 * (1 - u * u), 0.0),
    "epan2": lambda u: np.where(np.abs(u) < 1, 0.75 * (1 - u * u), 0.0),
    "biweight": lambda u: np.where(np.abs(u) < 1, 15 / 16 * (1 - u * u) ** 2, 0.0),
    "cosine": lambda u: np.where(np.abs(u) < 1, 1 + np.cos(2 * math.pi * u), 0.0),
    "gaussian": lambda u: np.exp(-0.5 * u * u) / math.sqrt(2 * math.pi),
    "parzen": lambda u: np.where(np.abs(u) <= 0.5, 4 / 3 - 8 * u * u + 8 * np.abs(u) ** 3,
                                 np.where(np.abs(u) <= 1, 8 / 3 * (1 - np.abs(u)) ** 3, 0.0)),
    "rectangle": lambda u: np.where(np.abs(u) < 1, 0.5, 0.0),
    "triangle": lambda u: np.where(np.abs(u) < 1, 1 - np.abs(u), 0.0),
}


def lpoly_fit(x: np.ndarray, y: np.ndarray, grid: np.ndarray, h: float, degree: int = 0,
              kernel: str = "epanechnikov") -> np.ndarray:
    """Local polynomial smooth of y on x at the points of `grid` ([R] lpoly, Methods and
    formulas): weighted least squares of degree `degree` with kernel weights K((x - g)/h)."""
    K = _KERNELS[kernel]
    out = np.full(len(grid), np.nan)
    for i, g0 in enumerate(grid):
        w = K((x - g0) / h)
        keep = w > 0
        if keep.sum() <= degree:
            continue
        d = x[keep] - g0
        Z = np.vander(d, degree + 1, increasing=True)
        sw = np.sqrt(w[keep])
        try:
            b, *_ = np.linalg.lstsq(Z * sw[:, None], y[keep] * sw, rcond=None)
        except np.linalg.LinAlgError:
            continue
        out[i] = b[0]
    return out


def rot_bandwidth(x: np.ndarray) -> float:
    """Default bandwidth: a reference-rule value, 0.9 min(sd, IQR/1.349) n^(-1/5). VERIFY: Stata
    uses the Fan-Gijbels rule-of-thumb plug-in estimator; the difference is checked in the
    1903 compat case."""
    n = len(x)
    sd = float(x.std(ddof=1))
    iqr = float(np.subtract(*np.percentile(x, [75, 25]))) / 1.349
    return 0.9 * min(sd, iqr if iqr > 0 else sd) * n ** -0.2


def lpoly(s: "Session", text: str):
    p, lopts, graw = _split(s, text, {"kernel": 3, "degree": 3, "bwidth": 2, "n": 1, "at": 2,
                                      "generate": 3, "nograph": 3, "noscatter": 4, "ci": 2,
                                      "level": 3, "msymbol": 3})
    yv, xv = _vars(s, p, 2)
    y, x = D.numeric(s, yv), D.numeric(s, xv)
    mask = D.base_mask(s, p) & D.valid(x) & D.valid(y)
    xx, yy = x[mask], y[mask]
    if len(xx) < 2:
        raise StataError(2001, "insufficient observations")
    kn = str(lopts.get("kernel", "epanechnikov"))
    full = [k for k in _KERNELS if k.startswith(kn)]
    if not full:
        raise StataError(198, f"kernel {kn} not allowed")
    kn = full[0]
    deg = int(lopts["degree"]) if isinstance(lopts.get("degree"), str) else 0
    h = float(lopts["bwidth"]) if isinstance(lopts.get("bwidth"), str) else rot_bandwidth(xx)
    npts = int(lopts["n"]) if isinstance(lopts.get("n"), str) else min(50, len(xx))
    if isinstance(lopts.get("at"), str):
        grid = D.numeric(s, s.expand_varlist(lopts["at"])[0])[mask]
        grid_all = D.numeric(s, s.expand_varlist(lopts["at"])[0])
        gm = D.valid(grid_all) & D.base_mask(s, p)
        grid = grid_all[gm]
    else:
        grid = np.linspace(xx.min(), xx.max(), npts)
    sm = lpoly_fit(xx, yy, grid, h, deg, kn)
    if isinstance(lopts.get("generate"), str):
        from ..core.dataset import Variable
        names = lopts["generate"].split()
        nobs = s.data.nobs
        vals = np.full(nobs, D.SYS)
        if len(names) == 1:
            if isinstance(lopts.get("at"), str):
                gm = D.valid(D.numeric(s, s.expand_varlist(lopts["at"])[0]))
                vals[gm & D.base_mask(s, p)] = np.where(np.isnan(sm), D.SYS, sm)
                s.data.add(Variable(names[0], "float", vals))
        else:
            if len(grid) > nobs:
                raise StataError(198, "number of grid points exceeds the number of observations")
            gx = np.full(nobs, D.SYS)
            gx[:len(grid)] = grid
            vals[:len(grid)] = np.where(np.isnan(sm), D.SYS, sm)
            s.data.add(Variable(names[0], "float", gx))
            s.data.add(Variable(names[1], "float", vals))
    if lopts.get("nograph"):
        return None, graw
    layers = []
    if not lopts.get("noscatter"):
        layers.append(Layer("scatter", xx.tolist(), yy.tolist(), {}, {"msymbol": "oh"}, D.var_title(s, yv)))
    ok = ~np.isnan(sm)
    layers.append(Layer("line", grid[ok].tolist(), sm[ok].tolist(), {}, {}, "lpoly smooth"))
    g = _graph(s, graw, layers, D.var_title(s, xv), D.var_title(s, yv),
               title="Local polynomial smooth",
               note=f"kernel = {kn}, degree = {deg}, bandwidth = {h:.4g}")
    g.opts["_ytitle"] = D.var_title(s, yv)
    return g, graw


# ---------------------------------------------------------------------------
# cusum
# ---------------------------------------------------------------------------

def cusum(s: "Session", text: str):
    """Cumulative sum of the standardized residuals of the regression of yvar on xvar, in the
    order of xvar ([R] cusum; VERIFY the test statistics, which are not displayed yet)."""
    p, lopts, graw = _split(s, text, {"generate": 3, "yfit": 4, "cshighlight": 4})
    yv, xv = _vars(s, p, 2)
    y, x = D.numeric(s, yv), D.numeric(s, xv)
    mask = D.base_mask(s, p) & D.valid(x) & D.valid(y)
    xx, yy = x[mask], y[mask]
    n = len(xx)
    if n < 3:
        raise StataError(2001, "insufficient observations")
    o = np.argsort(xx, kind="stable")
    xs, ys = xx[o], yy[o]
    Z = np.column_stack([np.ones(n), xs])
    b, *_ = np.linalg.lstsq(Z, ys, rcond=None)
    e = ys - Z @ b
    sd = math.sqrt(float(e @ e) / (n - 2))
    cs = np.cumsum(e / sd)
    if isinstance(lopts.get("generate"), str):
        from ..core.dataset import Variable
        full = np.full(s.data.nobs, D.SYS)
        idx = np.flatnonzero(mask)[o]
        full[idx] = cs
        s.data.add(Variable(lopts["generate"], "float", full))
    layers = [Layer("line", list(range(1, n + 1)), cs.tolist(), {}, {}, "")]
    g = _graph(s, graw, layers, "Order of " + D.var_title(s, xv), "Cusum")
    return g, graw


# ---------------------------------------------------------------------------
# gladder
# ---------------------------------------------------------------------------

def gladder(s: "Session", text: str):
    """Histograms of the variable under the ladder-of-powers transformations, each with the
    normal density overlaid, combined in a 3x3 graph ([R] ladder)."""
    from ..commands.basestat2 import _LADDER
    p, lopts, graw = _split(s, text, {"bin": 3, "start": 3, "width": 3, "percent": 3,
                                      "fraction": 3, "frequency": 3, "density": 3,
                                      "normal": 3, "l1title": 4, "ycommon": 4, "xcommon": 4})
    v = _vars(s, p, 1)[0]
    x = D.numeric(s, v)
    mask = D.base_mask(s, p) & D.valid(x)
    xs = x[mask]
    if len(xs) < 2:
        raise StataError(2000, "no observations")
    kids = []
    for lab, form, f in _LADDER:
        with np.errstate(all="ignore"):
            t = f(xs)
        t = t[np.isfinite(t)]
        if len(t) < 2:
            kids.append(Graph(family="twoway", opts={"legend": {"off": True}, "title": {"lines": [lab]}},
                              panels=[Panel("", [])]))
            continue
        h = D.hist_bins(t, None, {})
        mu, sd = float(t.mean()), float(t.std(ddof=1))
        lo, hi = float(t.min()), float(t.max())
        gx = np.linspace(lo, hi, 100)
        dens = np.exp(-0.5 * ((gx - mu) / sd) ** 2) / (sd * math.sqrt(2 * math.pi)) if sd > 0 else gx * 0
        layers = [Layer("histogram", h["centers"], h["heights"], {"width": h["width"]}, {}, ""),
                  Layer("line", gx.tolist(), dens.tolist(), {}, {"lcolor": "maroon"}, "")]
        kid = _graph(s, graw, layers, form if form != "x" else v, "Density")
        kid.opts = {k: val for k, val in kid.opts.items() if k not in ("name", "saving", "nodraw")}
        kid.opts["_xtitle"], kid.opts["_ytitle"] = form if form != "x" else D.var_title(s, v), "Density"
        kid.opts["title"] = {"lines": [lab]}
        kids.append(kid)
    opts = graph_options(graw)
    opts["cols"] = 3
    return Graph(family="combine", opts=opts, children=kids), graw
