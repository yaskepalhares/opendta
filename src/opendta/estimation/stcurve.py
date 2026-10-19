"""stcurve: curves of the survivor, hazard, cumulative hazard and cumulative incidence functions
after stcox, streg and stcrreg ([ST] stcurve). Covariates not given in at() are held at their
means in the estimation sample. Layout details carry the VERIFY mark
(compat/do/1901_survival.do)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import parse_standard
from .results import current
from .st import header, require_st

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_AT = {f"at{i}": 3 for i in range(1, 11)}
_NOT_YET = ("outfile", "width", "kernel", "noboundary", "fixedonly", "unconditional",
            "marginal", "addplot")


def _values(est, spec: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for tok in spec.split():
        name, eq, val = tok.partition("=")
        if not eq:
            raise StataError(198, f"at({spec}) invalid")
        try:
            out[name] = float(val)
        except ValueError:
            raise StataError(198, f"at({spec}) invalid") from None
    return out


def _xb(s: "Session", est, given: dict[str, float]) -> float:
    """Linear index at the covariate values: the means of the estimation sample, except for the
    variables in `given`."""
    cols = list(est.cols)
    known = {c.name for c in cols if c.values is not None}
    for k in given:
        if k not in known:
            raise StataError(198, f"{k} is not a covariate of the model")   # VERIFY
    xb = 0.0
    for i, c in enumerate(cols):
        if c.values is None:
            continue
        xb += float(est.b[i]) * (given[c.name] if c.name in given else float(np.mean(c.values)))
    if getattr(est, "constant", False):
        xb += float(est.b[len(cols)])
    return xb


def _step(times: np.ndarray, vals: np.ndarray, x0: float, x1: float, y0: float):
    xs, ys = [x0], [y0]
    cur = y0
    for t, v in zip(times, vals):
        if t > x1:
            break
        if t < x0:
            cur = v
            xs[0], ys[0] = x0, cur
            continue
        xs += [float(t), float(t)]
        ys += [cur, float(v)]
        cur = float(v)
    xs.append(float(x1))
    ys.append(cur)
    return xs, ys


def cmd_stcurve(s: "Session", args: str) -> None:
    from ..commands.graphcmd import finish
    from ..graphics.model import Graph, Layer, Panel
    from ..graphics.opts import axis_label, graph_options, split_layer_graph
    p = parse_standard(args)
    spec = {"survival": 1, "hazard": 1, "cumhaz": 2, "cif": 3, "at": 2, "range": 2, "noshow": 6,
            "outfile": 3, "width": 2, "kernel": 2, "noboundary": 3, "alpha1": 3, "fixedonly": 3,
            "unconditional": 3, "marginal": 3, "addplot": 3, **_AT}
    lopts, graw = split_layer_graph(p.options, spec)
    est = current(s)
    kind = ("cox" if est.cmd == "cox" else "crreg" if est.cmd == "stcrreg"
            else "streg" if hasattr(est, "dist") else "")
    if not kind:
        raise StataError(301, "last estimates not found")   # VERIFY: stcurve after another command
    info = require_st(s)
    for k in _NOT_YET:
        if lopts.get(k):
            raise StataError(198, f"option {k} is not available yet in OpenDTA")
    what = [k for k in ("survival", "hazard", "cumhaz", "cif") if lopts.get(k)]
    if len(what) != 1:
        raise StataError(198, "exactly one of survival, hazard, cumhaz, or cif must be specified")   # VERIFY
    what = what[0]
    if kind == "crreg" and what in ("survival", "hazard"):
        raise StataError(198, f"option {what} not allowed after stcrreg")   # VERIFY
    if kind != "crreg" and what == "cif":
        raise StataError(198, "option cif is allowed only after stcrreg")   # VERIFY
    if kind == "cox" and what == "hazard":
        raise StataError(198, "option hazard is not available yet in OpenDTA after stcox")
    if kind == "cox" and getattr(est, "strata_vars", None):
        raise StataError(198, "stcurve after a stratified stcox is not available yet in OpenDTA")
    if not lopts.get("noshow") and not lopts.get("nodraw"):   # VERIFY: nodraw prints no header
        header(s, info)
    ds = s.data
    tt = ds.get("_t").data
    mask = ds.esample if ds.esample is not None else np.ones(ds.nobs, dtype=bool)
    tt = tt[mask & (tt < SYS)]
    lo, hi = (float(tt.min()), float(tt.max())) if len(tt) else (0.0, 0.0)
    if lopts.get("range"):
        rg = str(lopts["range"]).split()
        if len(rg) != 2:
            raise StataError(198, "range() invalid")
        lo, hi = float(rg[0]), float(rg[1])
    sets = [(lab, _values(est, str(lopts[lab]))) for lab in ["at"] + list(_AT) if lopts.get(lab)]
    curves = sets or [("", {})]
    base = _baseline(s, est, kind)
    layers = []
    for i, (lab, given) in enumerate(curves):
        xb = _xb(s, est, given)
        xs, ys = _curve(est, kind, what, xb, base, lo, hi)
        text = str(lopts[lab]).replace("  ", " ") if lab else ""
        layers.append(Layer("line", xs, ys, {}, {"_pal": i} if len(curves) > 1 else {}, text))
    opts = graph_options(graw)
    sub = {"cox": "Cox proportional hazards regression", "crreg": "Competing-risks regression"}.get(kind)
    if sub is None:
        from .streg import DISTS
        sub = DISTS[est.dist].title + " regression"
    opts.setdefault("subtitle", {"lines": [sub]})
    opts["_xtitle"] = "analysis time"
    opts["_ytitle"] = {"survival": "Survival", "hazard": "Hazard function", "cumhaz": "Cumulative hazard",
                       "cif": "Cumulative Incidence"}[what]   # VERIFY
    if what == "survival" and "ylabel" not in graw:
        try:
            opts["ylabel"] = axis_label("0(.25)1, format(%4.2f)")
        except Exception:   # noqa: BLE001
            pass
    if len(curves) == 1:
        opts.setdefault("legend", {"off": True})
    g = Graph(family="twoway", opts=opts, panels=[Panel("", layers)])
    finish(s, g, graw)


def _baseline(s: "Session", est, kind: str):
    """(times, baseline cumulative hazard, baseline survivor) as step functions."""
    if kind == "crreg":
        return est.fg["ft"], est.fg["H0"], None
    if kind == "cox":
        from .stcox import _baseline as cox_baseline
        H_t, _H_t0, _hc, S_t, _info = cox_baseline(s, est)
        t = s.data.get("_t").data
        ok = (H_t < SYS) & (t < SYS)
        order = np.argsort(t[ok], kind="stable")
        tt, hh, ss = t[ok][order], H_t[ok][order], S_t[ok][order]
        keep = np.r_[tt[1:] != tt[:-1], True]
        return tt[keep], hh[keep], ss[keep]
    return None


def _curve(est, kind: str, what: str, xb: float, base, lo: float, hi: float):
    if kind == "streg":
        from .streg import DISTS
        dist, aft = DISTS[est.dist], est.aft
        th = est.th if dist.name == "ggamma" else (est.aux[0] if est.aux else 0.0)
        # after frailty(): the curves are conditional on a frailty of one (alpha1; VERIFY default)
        t = np.linspace(lo, hi, 200)
        t = np.where(t > 0, t, 1e-8) if what == "hazard" else t
        xbv = np.full_like(t, xb)
        if what == "survival":
            y = dist.surv(t, xbv, th, aft)
        elif what == "hazard":
            y = dist.hazard(t, xbv, th, aft)
        else:
            y = -np.log(dist.surv(t, xbv, th, aft))
        return [float(v) for v in t], [float(v) for v in y]
    times, H0, S0 = base
    r = float(np.exp(xb))
    if what == "survival":
        return _step(times, S0 ** r, lo, hi, 1.0)
    if what == "cif":
        return _step(times, 1 - np.exp(-H0 * r), lo, hi, 0.0)
    return _step(times, H0 * r, lo, hi, 0.0)
