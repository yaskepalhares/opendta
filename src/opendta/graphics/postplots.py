"""Gráficos de diagnóstico e de pós-estimação: marginsplot, lowess, qnorm,
pnorm, rvfplot, rvpplot ([R] marginsplot, [R] lowess, [R] diagnostic
plots, [R] regress postestimation diagnostic plots).

O desenho é o do motor do OpenDTA; os títulos e eixos seguem o padrão
descrito nos manuais (VERIFICAR detalhes de texto)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core.errors import StataError
from ..lang.syntax import parse_standard
from . import data as D
from .model import Graph, Layer, Panel
from .opts import graph_options, split_layer_graph
from .twoway import layer_spec

if TYPE_CHECKING:
    from ..session import Session


# ---------------------------------------------------------------------------
# lowess (Cleveland 1979, sem as iterações robustas, como o lowess do Stata)
# ---------------------------------------------------------------------------

def lowess_fit(x: np.ndarray, y: np.ndarray, bwidth: float = 0.8, mean: bool = False,
               noweight: bool = False) -> np.ndarray:
    """Valor suavizado em cada x ([R] lowess, Methods and formulas): com os
    dados ordenados por x, a janela de i vai de max(1, i-k) a min(i+k, N),
    k = floor((N·bwidth - 0.5)/2), cortada nas pontas (não deslocada); pesos
    tricúbicos com Δ = 1.0001·max(x(i+) - x(i), x(i) - x(i-)); regressão
    local ponderada (ou média ponderada com `mean`). Conferido no compat 0804."""
    n = len(x)
    order = np.argsort(x, kind="stable")
    xs, ys = x[order], y[order]
    half = max(int(math.floor((n * bwidth - 0.5) / 2)), 0)
    out = np.empty(n)
    for i in range(n):
        d = np.abs(xs - xs[i])
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        sel = np.arange(lo, hi)
        delta = 1.0001 * max(xs[hi - 1] - xs[i], xs[i] - xs[lo])
        if noweight or delta <= 0:
            wts = np.ones(len(sel))
        else:
            u = d[sel] / delta
            wts = np.where(u < 1, (1 - u ** 3) ** 3, 0.0)
        if mean:
            out[i] = float(wts @ ys[sel] / wts.sum()) if wts.sum() > 0 else ys[i]
            continue
        W = wts.sum()
        if W <= 0:
            out[i] = ys[i]
            continue
        xm = wts @ xs[sel] / W
        ym = wts @ ys[sel] / W
        sxx = wts @ (xs[sel] - xm) ** 2
        b = (wts @ ((xs[sel] - xm) * (ys[sel] - ym))) / sxx if sxx > 0 else 0.0
        out[i] = ym + b * (xs[i] - xm)
    res = np.empty(n)
    res[order] = out
    return res


LOWESS_OPTS = {"mean": 4, "noweight": 3, "bwidth": 2, "logit": 3, "adjust": 3, "generate": 3, "nograph": 6,
               "lineopts": 5}


def lowess(s: "Session", text: str):
    from ..core.dataset import Variable
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, {**layer_spec(), **LOWESS_OPTS})
    names = s.expand_varlist(p.varlist)
    if len(names) != 2:
        raise StataError(102 if len(names) < 2 else 103,
                         "too few variables specified" if len(names) < 2 else "too many variables specified")
    yv, xv = names
    y, x = D.numeric(s, yv), D.numeric(s, xv)
    mask = D.base_mask(s, p) & D.valid(x) & D.valid(y)
    bw = float(lopts["bwidth"]) if isinstance(lopts.get("bwidth"), str) else 0.8
    fit = lowess_fit(x[mask], y[mask], bw, bool(lopts.get("mean")), bool(lopts.get("noweight")))
    if lopts.get("adjust"):
        fit = fit * (y[mask].mean() / fit.mean())
    if lopts.get("logit"):
        with np.errstate(divide="ignore", invalid="ignore"):
            fit = np.log(fit / (1 - fit))
    if isinstance(lopts.get("generate"), str):
        nm = lopts["generate"].split()[0]
        vals = np.full(s.data.nobs, D.SYS)
        vals[mask] = fit.astype(np.float32).astype(float)     # variável float, como no Stata
        s.data.add(Variable(nm, "float", vals))
    if lopts.get("nograph"):
        return None, graw
    o = np.argsort(x[mask], kind="stable")
    opts = graph_options(graw)
    layers = [Layer("scatter", x[mask].tolist(), y[mask].tolist(), {}, {"msymbol": "oh"}, D.var_title(s, yv)),
              Layer("line", x[mask][o].tolist(), fit[o].tolist(), {}, {}, "lowess")]
    opts.setdefault("title", {"lines": ["Lowess smoother"]})
    opts.setdefault("note", {"lines": [f"bandwidth = {bw:g}"]})
    opts["_xtitle"], opts["_ytitle"] = D.var_title(s, xv), D.var_title(s, yv)
    opts.setdefault("legend", {"off": True})
    return Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), graw


# ---------------------------------------------------------------------------
# qnorm, pnorm
# ---------------------------------------------------------------------------

def _one_var(s: "Session", text: str, extra: dict):
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, {**layer_spec(), **extra})
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(102 if not names else 103,
                         "too few variables specified" if not names else "too many variables specified")
    v = names[0]
    x = D.numeric(s, v)
    mask = D.base_mask(s, p) & D.valid(x)
    xs = np.sort(x[mask])
    if len(xs) < 2:
        raise StataError(2000, "no observations")
    return v, xs, lopts, graw


def qnorm(s: "Session", text: str):
    """Quantis de x contra os da normal com a mesma média e desvio-padrão;
    posições (i - .5)/n (VERIFICAR: o Stata usa i/(n+1)?)."""
    v, xs, lopts, graw = _one_var(s, text, {"grid": 4})
    n = len(xs)
    q = st.norm.ppf((np.arange(1, n + 1) - 0.5) / n) * xs.std(ddof=1) + xs.mean()
    lo, hi = float(min(q.min(), xs.min())), float(max(q.max(), xs.max()))
    opts = graph_options(graw)
    layers = [Layer("scatter", q.tolist(), xs.tolist(), {}, {}, D.var_title(s, v)),
              Layer("line", [lo, hi], [lo, hi], {}, {"lcolor": "gs8"}, "")]
    opts["_xtitle"], opts["_ytitle"] = "Inverse Normal", D.var_title(s, v)
    opts.setdefault("legend", {"off": True})
    return Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), graw


def pnorm(s: "Session", text: str):
    """P_i = Φ((x_i - média)/dp) contra i/(n+1) ([R] diagnostic plots)."""
    v, xs, lopts, graw = _one_var(s, text, {"grid": 4})
    n = len(xs)
    emp = np.arange(1, n + 1) / (n + 1)
    th = st.norm.cdf((xs - xs.mean()) / xs.std(ddof=1))
    opts = graph_options(graw)
    layers = [Layer("scatter", emp.tolist(), th.tolist(), {}, {}, D.var_title(s, v)),
              Layer("line", [0, 1], [0, 1], {}, {"lcolor": "gs8"}, "")]
    opts["_xtitle"], opts["_ytitle"] = "Empirical P[i] = i/(N+1)", f"Normal F[({v}-m)/s]"
    opts.setdefault("legend", {"off": True})
    return Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), graw


# ---------------------------------------------------------------------------
# rvfplot, rvpplot (depois do regress)
# ---------------------------------------------------------------------------

def _resid_fit(s: "Session"):
    from ..estimation.postest import esample, xb_all
    from ..estimation.results import current
    est = current(s)
    if est.cmd not in ("regress", "anova"):
        raise StataError(301, "last estimates not found")   # VERIFICAR
    mask = esample(s)
    xb, ok = xb_all(s, est)
    y = s.data.get(est.depvar).data.astype(float)
    m = mask & ok & (y < D.SYS)
    return est, m, xb, y


def rvfplot(s: "Session", text: str):
    _p = parse_standard(text)
    lopts, graw = split_layer_graph(_p.options, layer_spec())
    est, m, xb, y = _resid_fit(s)
    opts = graph_options(graw)
    layers = [Layer("scatter", xb[m].tolist(), (y - xb)[m].tolist(), {}, {}, "Residuals")]
    opts["_xtitle"], opts["_ytitle"] = "Fitted values", "Residuals"
    opts.setdefault("yline", {"values": [0.0]})
    return Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), graw


def rvpplot(s: "Session", text: str):
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, layer_spec())
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(102, "too few variables specified")
    est, m, xb, y = _resid_fit(s)
    x = D.numeric(s, names[0])
    m = m & D.valid(x)
    opts = graph_options(graw)
    layers = [Layer("scatter", x[m].tolist(), (y - xb)[m].tolist(), {}, {}, "Residuals")]
    opts["_xtitle"], opts["_ytitle"] = D.var_title(s, names[0]), "Residuals"
    opts.setdefault("yline", {"values": [0.0]})
    return Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), graw


# ---------------------------------------------------------------------------
# marginsplot
# ---------------------------------------------------------------------------

def marginsplot(s: "Session", text: str):
    """Margens e intervalos do último margins. Com at() variando uma
    variável, ela vai no eixo x; senão, os níveis dos fatores."""
    info = getattr(s, "_margins_plot", None)
    if not info:
        raise StataError(301, "previous command was not margins")   # VERIFICAR
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, {**layer_spec(), "noci": 4, "recast": 6, "recastci": 7,
                                                "xdimension": 4, "plotdimension": 4, "allxlabels": 4})
    res = info["results"]
    vals, lo, hi = info["vals"], info["lo"], info["hi"]
    ats = info["ats"]
    varying = [k for k in (ats[0] if ats else {}) if len({a.get(k) for a in ats}) > 1] if len(ats) > 1 else []
    series: dict[str, list] = {}
    xlabels: dict = {}
    if varying:
        xv = varying[0]
        for i, (head, lab, at_lab) in enumerate(res):
            key = f"{head}={lab}" if head else (lab if info["dydx"] else "")
            x = ats[int(at_lab) - 1][xv] if at_lab else 0.0
            series.setdefault(key, []).append((x, vals[i], lo[i], hi[i]))
        xtitle = xv
    else:
        heads = [h for h, _l, _a in res]
        two = [h for h in heads if "#" in h]
        positions: dict = {}
        for i, (head, lab, at_lab) in enumerate(res):
            if two and "#" in head:
                a, b_ = lab.split("#", 1)
                key = f"{head.split('#')[1]}={b_}"
                xkey = a
            else:
                key = ""
                xkey = lab if not head else lab
            if xkey not in positions:
                positions[xkey] = len(positions)
                xlabels[positions[xkey]] = xkey
            series.setdefault(key, []).append((float(positions[xkey]), vals[i], lo[i], hi[i]))
        xtitle = (two[0].split("#")[0] if two else (heads[0] if heads and heads[0] else ""))
    layers = []
    pal = 0
    for key, pts in series.items():
        pts.sort()
        xs = [q[0] for q in pts]
        layers.append(Layer("connected", xs, [q[1] for q in pts], {}, {"_pal": pal}, key or info["title"]))
        if not lopts.get("noci"):
            layers.append(Layer("rcap", xs, [q[2] for q in pts], {"y2": [q[3] for q in pts]}, {"_pal": pal}, ""))
        pal += 1
    opts = graph_options(graw)
    lev = info["level"]
    ttl = info["title"].title().replace("Of", "of")
    opts.setdefault("title", {"lines": [f"{ttl} with {lev:g}% CIs"]})
    ylab = info["expr"].split(",")[0]
    if info["dydx"]:
        ylab = f"Effects on {ylab}"
    opts["_ytitle"] = ylab
    opts["_xtitle"] = xtitle
    if xlabels:
        opts.setdefault("xlabel", {"ticks": [[float(k), v] for k, v in xlabels.items()]})
    if len(series) <= 1:
        opts.setdefault("legend", {"off": True})
    return Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), graw
