"""Dados dos gráficos: amostra (if/in), rótulos, grupos de by() e as contas
de histogram, kdensity, ajustes (lfit/qfit) e caixas (graph box).

Fórmulas tiradas dos manuais [G-2] e [R] (histogram, kdensity, graph box);
pontos não documentados têm a marca VERIFICAR."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Iterator

import numpy as np

from ..core import missing as M
from ..core import stats as S
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import Parsed

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def base_mask(s: "Session", p: Parsed) -> np.ndarray:
    from ..commands._util import touse
    return touse(s, p)


def numeric(s: "Session", name: str) -> np.ndarray:
    v = s.data.get(name)
    if v.is_string:
        raise StataError(109, "type mismatch")
    return np.asarray(v.data, dtype=np.float64)


def valid(x: np.ndarray) -> np.ndarray:
    return x < SYS


def nan(x: np.ndarray) -> np.ndarray:
    y = np.asarray(x, dtype=np.float64).copy()
    y[y >= SYS] = np.nan
    return y


def var_title(s: "Session", name: str) -> str:
    v = s.data.get(name)
    return v.label or v.name


def value_labels(s: "Session", name: str) -> dict:
    v = s.data.get(name)
    if not v.value_label:
        return {}
    return dict(s.data.value_labels.get(v.value_label, {}))


def value_text(s: "Session", name: str, x) -> str:
    v = s.data.get(name)
    if v.is_string:
        return str(x)
    if x >= SYS:
        return M.missing_name(x)
    labs = value_labels(s, name)
    if float(x).is_integer() and int(x) in labs:
        return labs[int(x)]
    return format_value(float(x), v.fmt, pad=False).strip()


def groups(s: "Session", mask: np.ndarray, byvars: list[str], *, missing: bool = False
           ) -> Iterator[tuple[list, np.ndarray]]:
    """Grupos de by()/over(): (valores, máscara), em ordem crescente."""
    if not byvars:
        yield [], mask
        return
    cols = []
    for b in byvars:
        v = s.data.get(b)
        cols.append(np.asarray(v.raw if v.is_string else v.data))
    keep = mask.copy()
    if not missing:
        for b, c in zip(byvars, cols):
            if not s.data.get(b).is_string:
                keep &= c < SYS
    idx = np.flatnonzero(keep)
    keys = sorted({tuple(c[i] for c in cols) for i in idx},
                  key=lambda t: tuple((0, x) if not isinstance(x, str) else (1, x) for x in t))
    for key in keys:
        m = keep.copy()
        for c, k in zip(cols, key):
            m &= (c == k)
        yield list(key), m


def group_title(s: "Session", byvars: list[str], key: list) -> str:
    return ", ".join(value_text(s, b, k) for b, k in zip(byvars, key))


def weights(s: "Session", p: Parsed, mask: np.ndarray, allowed=("fweight", "aweight")):
    if p.weight is None:
        return None, mask
    from ..commands.summarize import weights as _w
    w, _, mask = _w(s, p, mask, allowed)
    return w, mask


# ---------------------------------------------------------------------------
# histogram
# ---------------------------------------------------------------------------

def hist_bins(x: np.ndarray, w: np.ndarray | None, opts: dict) -> dict:
    """Barras do histograma. opts: bin, width, start, discrete, e a escala
    (density, fraction, frequency, percent). Devolve centros, alturas,
    largura, início e número de barras."""
    n = len(x)
    if n == 0:
        raise StataError(2000, "no observations")
    ww = np.ones(n) if w is None else w
    W = float(ww.sum())
    xmin, xmax = float(x.min()), float(x.max())
    if opts.get("discrete"):
        vals = np.unique(x)
        if "width" in opts:
            width = float(opts["width"])
        else:
            # VERIFICAR: largura = menor diferença entre valores distintos
            width = float(np.min(np.diff(vals))) if len(vals) > 1 else 1.0
        start = float(opts.get("start", xmin))
        k = int(math.floor((xmax - start) / width + 1e-9)) + 1
        idx = np.floor((x - start) / width + 0.5).astype(int)
        idx = np.clip(idx, 0, k - 1)
        centers = start + width * np.arange(k)
    else:
        start = float(opts["start"]) if "start" in opts else xmin
        if start > xmin:
            raise StataError(198, "start() must be less than or equal to the minimum")   # VERIFICAR
        if "width" in opts:
            width = float(opts["width"])
            # o máximo entra na última barra (compat 0701: width(2) start(0)
            # com máximo 10 dá 5 barras)
            k = max(1, int(math.ceil((xmax - start) / width - 1e-9)))
        else:
            if "bin" in opts:
                k = int(opts["bin"])
            else:
                # [R] histogram: k = min(sqrt(N), 10 ln(N)/ln(10)), truncado
                k = int(min(math.sqrt(W), 10 * math.log10(W)))
                k = max(k, 1)
            width = (xmax - start) / k if xmax > start else 1.0
        # VERIFICAR: intervalos [a, b), o último fechado em xmax
        idx = np.floor((x - start) / width + 1e-12 * max(1.0, abs(xmax))).astype(int)
        idx = np.clip(idx, 0, k - 1)
        centers = start + width * (np.arange(k) + 0.5)
    counts = np.bincount(idx, weights=ww, minlength=k)[:k]
    scale = opts.get("scale", "density")
    if scale == "frequency":
        h = counts
    elif scale == "fraction":
        h = counts / W
    elif scale == "percent":
        h = 100 * counts / W
    else:
        h = counts / (W * width)
    return {"centers": centers.tolist(), "heights": h.tolist(), "width": width, "start": start,
            "bins": k, "N": W, "min": xmin, "max": xmax}


# ---------------------------------------------------------------------------
# kdensity
# ---------------------------------------------------------------------------

_SQ5 = math.sqrt(5.0)


def _kernel(name: str):
    if name == "epanechnikov":
        return lambda z: np.where(np.abs(z) < _SQ5, 3 / (4 * _SQ5) * (1 - z * z / 5), 0.0)
    if name == "epan2":
        return lambda z: np.where(np.abs(z) < 1, 0.75 * (1 - z * z), 0.0)
    if name == "biweight":
        return lambda z: np.where(np.abs(z) < 1, 15 / 16 * (1 - z * z) ** 2, 0.0)
    if name == "cosine":
        return lambda z: np.where(np.abs(z) < 0.5, 1 + np.cos(2 * np.pi * z), 0.0)
    if name == "gaussian":
        return lambda z: np.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)
    if name == "parzen":
        def k(z):
            a = np.abs(z)
            return np.where(a <= 0.5, 4 / 3 - 8 * z * z + 8 * a ** 3,
                            np.where(a <= 1, 8 * (1 - a) ** 3 / 3, 0.0))
        return k
    if name == "rectangle":
        return lambda z: np.where(np.abs(z) < 1, 0.5, 0.0)
    if name == "triangle":
        return lambda z: np.where(np.abs(z) < 1, 1 - np.abs(z), 0.0)
    raise StataError(198, f"kernel {name} not allowed")   # VERIFICAR mensagem


KERNELS = ("epanechnikov", "epan2", "biweight", "cosine", "gaussian", "parzen", "rectangle", "triangle")


def kernel_name(text: str | None) -> str:
    if not text:
        return "epanechnikov"
    t = text.strip().lower()
    if t == "epan2":
        return "epan2"
    for k in KERNELS:
        if k != "epan2" and len(t) >= 2 and k.startswith(t):
            return k
    raise StataError(198, f"kernel({text}) invalid")   # VERIFICAR


def default_bandwidth(x: np.ndarray, w: np.ndarray | None) -> float:
    """[R] kdensity: h = 0.9 m / n^(1/5), m = min(sd, IQR/1.349)."""
    order = np.argsort(x, kind="stable")
    xs = x[order]
    ws = None if w is None else w[order]
    mom = S.moments(xs, ws, "aweight" if w is not None else "")
    n = float(mom.N if w is None else ws.sum())
    q1 = S.percentile(xs, 25, ws)
    q3 = S.percentile(xs, 75, ws)
    m = min(mom.sd, (q3 - q1) / 1.349)
    if m <= 0:
        m = mom.sd
    return 0.9 * m / n ** 0.2


def kdensity(x: np.ndarray, w: np.ndarray | None, *, kernel: str = "epanechnikov",
             bw: float | None = None, n: int | None = None, at: np.ndarray | None = None) -> dict:
    N = len(x)
    if N == 0:
        raise StataError(2000, "no observations")
    h = bw if bw is not None else default_bandwidth(x, w)
    if at is None:
        npts = n if n is not None else min(N, 50)
        # grade de mínimo - h a máximo + h (compat 0701)
        lo, hi = float(x.min()) - h, float(x.max()) + h
        grid = np.linspace(lo, hi, npts) if npts > 1 else np.array([lo])
    else:
        grid = at
    K = _kernel(kernel)
    ww = np.ones(N) if w is None else w / w.sum() * N
    z = (grid[:, None] - x[None, :]) / h
    dens = (K(z) * ww[None, :]).sum(axis=1) / (N * h)
    return {"x": grid, "d": dens, "bwidth": h, "n": len(grid), "kernel": kernel}


# ---------------------------------------------------------------------------
# ajustes
# ---------------------------------------------------------------------------

def poly_fit(x: np.ndarray, y: np.ndarray, degree: int, grid: np.ndarray, ci: bool = False,
             level: float = 95.0, w: np.ndarray | None = None):
    X = np.vander(x - x.mean(), degree + 1, increasing=True)
    Xg = np.vander(grid - x.mean(), degree + 1, increasing=True)
    ww = np.ones(len(x)) if w is None else w
    XtX = X.T @ (X * ww[:, None])
    beta = np.linalg.solve(XtX, X.T @ (ww * y))
    fit = Xg @ beta
    if not ci:
        return fit, None, None
    resid = y - X @ beta
    df = len(x) - (degree + 1)
    s2 = float((ww * resid * resid).sum() / df)
    V = s2 * np.linalg.inv(XtX)
    se = np.sqrt(np.einsum("ij,jk,ik->i", Xg, V, Xg))
    from scipy import stats as _st
    t = float(_st.t.ppf(0.5 + level / 200, df))
    return fit, fit - t * se, fit + t * se


# ---------------------------------------------------------------------------
# graph box
# ---------------------------------------------------------------------------

def box_stats(x: np.ndarray, w: np.ndarray | None = None) -> dict:
    """[G-2] graph box: quartis pela definição de summarize; bigodes até os
    valores adjacentes (dentro de 1,5 IQR); o resto é ponto externo."""
    order = np.argsort(x, kind="stable")
    xs = x[order]
    ws = None if w is None else w[order]
    q1 = S.percentile(xs, 25, ws)
    med = S.percentile(xs, 50, ws)
    q3 = S.percentile(xs, 75, ws)
    iqr = q3 - q1
    lo_lim, hi_lim = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    inside = xs[(xs >= lo_lim) & (xs <= hi_lim)]
    lo = float(inside.min()) if len(inside) else q1
    hi = float(inside.max()) if len(inside) else q3
    out = xs[(xs < lo_lim) | (xs > hi_lim)]
    return {"q1": q1, "med": med, "q3": q3, "lo": lo, "hi": hi, "out": out.tolist()}


# ---------------------------------------------------------------------------
# graph bar: estatísticas
# ---------------------------------------------------------------------------

STATS = ("mean", "median", "p", "sum", "rawsum", "count", "percent", "max", "min", "sd",
         "first", "last", "firstnm", "lastnm", "iqr", "semean")


def stat_value(stat: str, x: np.ndarray, w: np.ndarray | None) -> float:
    ok = x < SYS
    x = x[ok]
    if w is not None:
        w = w[ok]
    if stat == "count":
        return float(len(x) if w is None else w.sum())
    if len(x) == 0:
        return SYS
    ww = np.ones(len(x)) if w is None else w
    if stat == "mean":
        return float((ww * x).sum() / ww.sum())
    if stat in ("sum", "rawsum"):
        return float((ww * x).sum()) if stat == "sum" else float(x.sum())
    if stat == "max":
        return float(x.max())
    if stat == "min":
        return float(x.min())
    if stat == "sd":
        if len(x) < 2:
            return SYS
        m = (ww * x).sum() / ww.sum()
        return float(math.sqrt((ww * (x - m) ** 2).sum() / (ww.sum() - 1)))
    if stat == "semean":
        sd = stat_value("sd", x, w)
        return sd / math.sqrt(ww.sum()) if sd < SYS else SYS
    if stat in ("first", "firstnm"):
        return float(x[0])
    if stat in ("last", "lastnm"):
        return float(x[-1])
    order = np.argsort(x, kind="stable")
    xs, ws = x[order], (None if w is None else w[order])
    if stat == "median":
        return S.percentile(xs, 50, ws)
    if stat == "iqr":
        return S.percentile(xs, 75, ws) - S.percentile(xs, 25, ws)
    if stat.startswith("p") and stat[1:].isdigit():
        return S.percentile(xs, float(stat[1:]), ws)
    raise StataError(198, f"({stat}) invalid statistic")   # VERIFICAR
