"""histogram, kdensity, graph bar/hbar, graph box/hbox, graph pie,
graph matrix and graph combine ([R] histogram, [R] kdensity, [G-2] graph ...).

Each function builds a Graph and also returns the graph options as text
(name(), saving(), nodraw stay with the calling command)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, parse_standard
from ..lang.words import strip_outer_quotes
from . import data as D
from .model import Graph, Layer, Panel
from .opts import graph_options, split_layer_graph, split_sub
from .twoway import Builder, HIST_LAYER, finish, layer_spec

if TYPE_CHECKING:
    from ..session import Session


def _gfmt(x: float) -> str:
    """Numbers in the histogram message: %9.0g without spaces (VERIFY)."""
    from ..core.formats import format_value
    return format_value(float(x), "%9.0g", pad=False).strip()


# ---------------------------------------------------------------------------
# histogram
# ---------------------------------------------------------------------------

HIST_OPTS = {**HIST_LAYER, "normal": 4, "kdensity": 4, "addlabels": 4, "horizontal": 3,
             "normopts": 8, "kdenopts": 7, "addlabopts": 8}


def histogram(s: "Session", text: str) -> tuple[Graph, dict, str]:
    p = parse_standard(text)
    spec = {**layer_spec(), **HIST_OPTS}
    lopts, graw = split_layer_graph(p.options, spec)
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103 if len(names) > 1 else 102,
                         "too many variables specified" if len(names) > 1 else "too few variables specified")
    v = names[0]
    if s.data.get(v).is_string:
        raise StataError(109, "type mismatch")
    x = D.numeric(s, v)
    mask = D.base_mask(s, p) & D.valid(x)
    w, mask = D.weights(s, p, mask, ("fweight",))
    if not mask.any():
        raise StataError(2000, "no observations")
    from .twoway import hist_opts, hist_ytitle
    ho = hist_opts(lopts)
    total = D.hist_bins(x[mask], None if w is None else w[mask], ho)
    if ho.get("discrete"):
        msg = f"(start={_gfmt(total['start'])}, width={_gfmt(total['width'])})"
    else:
        msg = f"(bin={total['bins']}, start={_gfmt(total['start'])}, width={_gfmt(total['width'])})"
    style = {k: (strip_outer_quotes(val) if isinstance(val, str) else val) for k, val in lopts.items()
             if k not in HIST_OPTS}
    if lopts.get("horizontal"):
        style["horizontal"] = True
    opts = graph_options(graw)
    g = Graph(family="twoway", opts=opts)
    scale = ho.get("scale", "density")

    def panel_layers(m: np.ndarray) -> list[Layer]:
        xx = x[m]
        ww = None if w is None else w[m]
        if len(xx) == 0:
            return []
        ho2 = dict(ho)
        ho2["start"] = total["start"]
        if ho.get("discrete"):
            ho2["width"] = total["width"]
        else:
            ho2["bin"] = total["bins"]
            ho2.pop("width", None)
            ho2["_xmax"] = total["max"]
        h = _bins_fixed(xx, ww, ho2, total)
        out = [Layer("histogram", h["centers"], h["heights"], {"width": h["width"]}, dict(style), "")]
        W = float(len(xx) if ww is None else ww.sum())
        factor = {"density": 1.0, "fraction": h["width"], "frequency": h["width"] * W,
                  "percent": 100 * h["width"]}[scale]
        if ho.get("discrete"):
            lo, hi = total["start"] - total["width"] / 2, total["max"] + total["width"] / 2
        else:
            lo, hi = total["start"], total["start"] + total["width"] * total["bins"]
        if lopts.get("normal"):
            mu = float(np.average(xx, weights=ww))
            sd = float(math.sqrt(np.average((xx - mu) ** 2, weights=ww) * W / (W - 1))) if W > 1 else 0.0
            grid = np.linspace(lo, hi, 100)
            dens = np.exp(-0.5 * ((grid - mu) / sd) ** 2) / (sd * math.sqrt(2 * math.pi)) if sd > 0 else grid * 0
            out.append(Layer("line", grid.tolist(), (dens * factor).tolist(), {}, {"lcolor": "maroon"}, ""))
        if lopts.get("kdensity"):
            kd = D.kdensity(xx, ww)
            out.append(Layer("line", kd["x"].tolist(), (kd["d"] * factor).tolist(), {}, {"lcolor": "navy"}, ""))
        if lopts.get("addlabels"):
            out[0].extra["blabel"] = True
        return out

    byo = opts.get("by")
    if byo:
        byvars = s.expand_varlist(" ".join(byo["vars"]))
        for key, gm in D.groups(s, mask, byvars, missing=bool(byo.get("missing"))):
            lys = panel_layers(gm)
            if lys:
                g.panels.append(Panel(D.group_title(s, byvars, key), lys))
        if byo.get("total"):
            g.panels.append(Panel("Total", panel_layers(mask)))
        opts.setdefault("note", {"lines": ["Graphs by " + ", ".join(D.var_title(s, b) for b in byvars)]})
    else:
        g.panels.append(Panel("", panel_layers(mask)))
    if lopts.get("horizontal"):
        opts["_xtitle"], opts["_ytitle"] = hist_ytitle(lopts), D.var_title(s, v)
    else:
        opts["_xtitle"], opts["_ytitle"] = D.var_title(s, v), hist_ytitle(lopts)
    opts.setdefault("legend", {"off": True})
    return g, graw, msg


def _bins_fixed(x: np.ndarray, w, ho: dict, total: dict) -> dict:
    """Bins with the start, width and number of bins of the whole sample
    (with by(), the panels use the same bins)."""
    k = total["bins"]
    width = total["width"]
    start = total["start"]
    ww = np.ones(len(x)) if w is None else w
    W = float(ww.sum())
    if ho.get("discrete"):
        idx = np.floor((x - start) / width + 0.5).astype(int)
        centers = start + width * np.arange(k)
    else:
        idx = np.floor((x - start) / width + 1e-12 * max(1.0, abs(total["max"]))).astype(int)
        centers = start + width * (np.arange(k) + 0.5)
    idx = np.clip(idx, 0, k - 1)
    counts = np.bincount(idx, weights=ww, minlength=k)[:k]
    scale = ho.get("scale", "density")
    h = {"frequency": counts, "fraction": counts / W, "percent": 100 * counts / W}.get(scale)
    if h is None:
        h = counts / (W * width)
    return {"centers": centers.tolist(), "heights": h.tolist(), "width": width}


# ---------------------------------------------------------------------------
# kdensity
# ---------------------------------------------------------------------------

KDEN_OPTS = {"kernel": 1, "bwidth": 2, "generate": 3, "n": 1, "at": 2, "nograph": 6, "normal": 4,
             "student": 3, "normopts": 8, "stopts": 6, "range": 1}


def kdensity(s: "Session", text: str) -> tuple[Graph | None, dict]:
    p = parse_standard(text)
    spec = {**layer_spec(), **KDEN_OPTS}
    lopts, graw = split_layer_graph(p.options, spec)
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103 if len(names) > 1 else 102,
                         "too many variables specified" if len(names) > 1 else "too few variables specified")
    v = names[0]
    if s.data.get(v).is_string:
        raise StataError(109, "type mismatch")
    x = D.numeric(s, v)
    mask = D.base_mask(s, p) & D.valid(x)
    w, mask = D.weights(s, p, mask, ("fweight", "aweight"))
    if not mask.any():
        raise StataError(2000, "no observations")
    kname = D.kernel_name(lopts.get("kernel") if isinstance(lopts.get("kernel"), str) else None)
    bw = float(lopts["bwidth"]) if isinstance(lopts.get("bwidth"), str) else None
    npts = int(lopts["n"]) if isinstance(lopts.get("n"), str) else None
    at = None
    at_rows = None
    if isinstance(lopts.get("at"), str):
        av = D.numeric(s, s.expand_varlist(lopts["at"])[0])
        at_rows = np.flatnonzero(D.valid(av))
        at = av[at_rows]
    xx = x[mask]
    ww = None if w is None else w[mask]
    kd = D.kdensity(xx, ww, kernel=kname, bw=bw, n=npts, at=at)
    gen = lopts.get("generate")
    if isinstance(gen, str):
        _generate(s, gen, kd, at_rows)
    s.r = {"scale": 1.0, "n": float(kd["n"]), "bwidth": float(kd["bwidth"]), "kernel": kname}   # VERIFY
    if lopts.get("nograph"):
        return None, graw
    style = {k: (strip_outer_quotes(val) if isinstance(val, str) else val) for k, val in lopts.items()
             if k not in KDEN_OPTS}
    opts = graph_options(graw)
    layers = [Layer("kdensity", kd["x"].tolist(), kd["d"].tolist(), {}, style, "Kernel density estimate")]
    if lopts.get("normal"):
        mu = float(np.average(xx, weights=ww))
        n = float(len(xx) if ww is None else ww.sum())
        sd = float(math.sqrt(np.average((xx - mu) ** 2, weights=ww) * n / (n - 1))) if n > 1 else 0.0
        gx = np.asarray(kd["x"])
        nd = np.exp(-0.5 * ((gx - mu) / sd) ** 2) / (sd * math.sqrt(2 * math.pi))
        layers.append(Layer("line", gx.tolist(), nd.tolist(), {}, {}, "Normal density"))
    g = Graph(family="twoway", opts=opts, panels=[Panel("", layers)])
    opts.setdefault("title", {"lines": ["Kernel density estimate"]})
    opts.setdefault("note", {"lines": [f"kernel = {kname}, bandwidth = {_fmt4(kd['bwidth'])}"]})
    opts["_xtitle"] = D.var_title(s, v)
    opts["_ytitle"] = "Density"
    return g, graw


def _fmt4(x: float) -> str:
    # VERIFY: Stata's note uses %9.4f without spaces
    t = f"{x:.4f}"
    return t[1:] if t.startswith("0.") else t


def _generate(s: "Session", spec: str, kd: dict, at_rows) -> None:
    from ..core.dataset import Variable
    names = spec.split()
    ds = s.data
    if at_rows is not None:
        if len(names) != 1:
            raise StataError(198, "only one variable may be specified in generate() with at()")   # VERIFY
        dname = names[0]
        vals = np.full(ds.nobs, M.SYSMISS)
        vals[at_rows] = kd["d"]
        ds.add(Variable(dname, "double", vals, fmt="%10.0g"))
        return
    if len(names) != 2:
        raise StataError(198, "generate() requires two new variable names")   # VERIFY
    for nm in names:
        if ds.has(nm):
            raise StataError(110, f"{nm} already defined")
    n = len(kd["x"])
    if n > ds.nobs:
        ds.set_obs(n)   # VERIFY: does Stata increase _N when n() > _N?
    xs = np.full(ds.nobs, M.SYSMISS)
    ys = np.full(ds.nobs, M.SYSMISS)
    xs[:n] = kd["x"]
    ys[:n] = kd["d"]
    # double with format %10.0g (compat 0701: list shows -1.0146058)
    ds.add(Variable(names[0], "double", xs, fmt="%10.0g"))
    ds.add(Variable(names[1], "double", ys, fmt="%10.0g"))


# ---------------------------------------------------------------------------
# graph bar / hbar
# ---------------------------------------------------------------------------

BAR_OPTS = {"over": 4, "asyvars": 4, "stack": 5, "percentages": 7, "blabel": 3, "bargap": 6,
            "outergap": 6, "nofill": 6, "missing": 4, "cw": 2, "yvaroptions": 6, "showyvars": 5, "bar": 3, "nolabel": 7, "exclude0": 8, "yalternate": 4, "xalternate": 4}


def _yvars(s: "Session", text: str, default: str = "mean") -> list[tuple[str, str, str]]:
    """'(mean) a b (sum) c d=e' -> [(stat, var, label)]."""
    out = []
    stat = default
    t = text.strip()
    i = 0
    while i < len(t):
        if t[i].isspace():
            i += 1
            continue
        if t[i] == "(":
            j = t.index(")", i)
            stat = t[i + 1:j].strip().lower()
            if stat in ("p50",):
                stat = "median"
            if stat not in D.STATS and not (stat.startswith("p") and stat[1:].isdigit()):
                raise StataError(198, f"({stat}) invalid statistic")   # VERIFY
            i = j + 1
            continue
        j = i
        while j < len(t) and not t[j].isspace() and t[j] != "(":
            j += 1
        tok = t[i:j]
        i = j
        if "=" in tok:
            new, old = tok.split("=", 1)
            out.append((stat, s.expand_varlist(old)[0], new))
        else:
            for nm in s.expand_varlist(tok):
                out.append((stat, nm, ""))
    return out


_STAT_TITLE = {"mean": "mean", "median": "median", "sum": "sum", "rawsum": "sum", "count": "count",
               "percent": "percent", "max": "max", "min": "min", "sd": "sd", "first": "first",
               "last": "last", "firstnm": "first", "lastnm": "last", "iqr": "iqr", "semean": "semean"}


def _over_list(graw: dict, key: str = "over") -> list[dict]:
    """over() may be repeated; parse_options keeps only the last one, so
    we reopen the raw text (graw['_raw'])."""
    out = []
    from ..lang.syntax import parse_options
    for n, a in parse_options(graw.get("_raw", "")):
        if len(n) >= 4 and "over".startswith(n):
            main, subs = split_sub(a or "")
            out.append({"var": main.strip(), "subs": subs})
    return out


def bar(s: "Session", text: str, family: str) -> tuple[Graph, dict]:
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, BAR_OPTS)
    graw["_raw"] = p.options
    overs = _over_list(graw)
    if len(overs) > 3:
        raise StataError(198, "too many over() options")   # VERIFY
    ys = _yvars(s, p.varlist) if p.varlist.strip() else []
    mask = D.base_mask(s, p)
    w, mask = D.weights(s, p, mask, ("fweight", "aweight", "pweight", "iweight"))
    if lopts.get("cw"):
        for _, nm, _ in ys:
            mask &= D.valid(D.numeric(s, nm))
    overvars = [s.expand_varlist(o["var"])[0] for o in overs]
    if not ys and not overvars:
        raise StataError(102, "too few variables specified")   # VERIFY
    opts = graph_options({k: v for k, v in graw.items() if k != "_raw"})

    def series_for(m: np.ndarray) -> tuple[list[str], list[dict]]:
        keys = list(D.groups(s, m, overvars, missing=bool(lopts.get("missing")))) if overvars else [([], m)]
        cats = [D.group_title(s, overvars, k) if overvars else "" for k, _ in keys]
        series = []
        if not ys:
            tot = float(m.sum() if w is None else w[m].sum())
            vals = [float(gm.sum() if w is None else w[gm].sum()) / tot * 100 if tot else M.SYSMISS
                    for _, gm in keys]
            series.append({"label": "percent", "values": vals})
        for stat, nm, lab in ys:
            xv = D.numeric(s, nm)
            vals = []
            for _, gm in keys:
                if stat == "percent":
                    tot = float(xv[m & D.valid(xv)].sum())
                    sub = float(xv[gm & D.valid(xv)].sum())
                    vals.append(100 * sub / tot if tot else M.SYSMISS)
                else:
                    vals.append(D.stat_value(stat, xv[gm], None if w is None else w[gm]))
            label = lab or f"{_STAT_TITLE.get(stat, stat)} of {nm}"
            series.append({"label": label, "values": [None if v >= M.SYSMISS else v for v in vals]})
        if lopts.get("asyvars") and overvars and len(series) == 1:
            # the first over() becomes the series (colors and legend)
            vals = series[0]["values"]
            series = [{"label": c, "values": [v]} for c, v in zip(cats, vals)]
            cats = [""]
        if lopts.get("percentages") and len(series) > 1:
            for k in range(len(cats)):
                tot = sum((ser["values"][k] or 0) for ser in series)
                for ser in series:
                    if ser["values"][k] is not None and tot:
                        ser["values"][k] = 100 * ser["values"][k] / tot
        return cats, series

    g = Graph(family=family, opts=opts)
    blabel = isinstance(lopts.get("blabel"), str) and lopts["blabel"].split(",")[0].strip() in ("bar", "total")
    byo = opts.get("by")
    if byo:
        byvars = s.expand_varlist(" ".join(byo["vars"]))
        for key, gm in D.groups(s, mask, byvars):
            cats, series = series_for(gm)
            g.panels.append(Panel(D.group_title(s, byvars, key), [Layer("bars", extra={
                "cats": cats, "series": series, "stack": bool(lopts.get("stack")), "blabel": blabel})]))
        opts.setdefault("note", {"lines": ["Graphs by " + ", ".join(D.var_title(s, b) for b in byvars)]})
    else:
        cats, series = series_for(mask)
        g.panels.append(Panel("", [Layer("bars", extra={"cats": cats, "series": series,
                                                       "stack": bool(lopts.get("stack")), "blabel": blabel})]))
    first = g.panels[0].layers[0].extra["series"]
    opts["_ytitle"] = first[0]["label"] if len(first) == 1 else ""
    opts["_xtitle"] = ""
    return g, graw


# ---------------------------------------------------------------------------
# graph box / hbox
# ---------------------------------------------------------------------------

BOX_OPTS = {"over": 4, "asyvars": 4, "nooutsides": 5, "missing": 4, "cw": 2, "medtype": 4,
            "box": 3, "marker": 3, "yvaroptions": 6, "showyvars": 5, "nolabel": 7}


def box(s: "Session", text: str, family: str) -> tuple[Graph, dict]:
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, BOX_OPTS)
    graw["_raw"] = p.options
    overs = _over_list(graw)
    names = s.expand_varlist(p.varlist)
    if not names:
        raise StataError(102, "too few variables specified")
    mask = D.base_mask(s, p)
    w, mask = D.weights(s, p, mask, ("fweight", "aweight", "pweight"))
    overvars = [s.expand_varlist(o["var"])[0] for o in overs]
    opts = graph_options({k: v for k, v in graw.items() if k != "_raw"})

    def series_for(m: np.ndarray):
        keys = list(D.groups(s, m, overvars, missing=bool(lopts.get("missing")))) if overvars else [([], m)]
        cats = [D.group_title(s, overvars, k) if overvars else "" for k, _ in keys]
        series = []
        for nm in names:
            xv = D.numeric(s, nm)
            st = []
            for _, gm in keys:
                mm = gm & D.valid(xv)
                if not mm.any():
                    st.append(None)
                    continue
                b = D.box_stats(xv[mm], None if w is None else w[mm])
                if lopts.get("nooutsides"):
                    b["out"] = []
                st.append(b)
            series.append({"label": D.var_title(s, nm), "stats": st})
        return cats, series

    g = Graph(family=family, opts=opts)
    byo = opts.get("by")
    if byo:
        byvars = s.expand_varlist(" ".join(byo["vars"]))
        for key, gm in D.groups(s, mask, byvars):
            cats, series = series_for(gm)
            g.panels.append(Panel(D.group_title(s, byvars, key),
                                  [Layer("boxes", extra={"cats": cats, "series": series})]))
        opts.setdefault("note", {"lines": ["Graphs by " + ", ".join(D.var_title(s, b) for b in byvars)]})
    else:
        cats, series = series_for(mask)
        g.panels.append(Panel("", [Layer("boxes", extra={"cats": cats, "series": series})]))
    opts["_ytitle"] = D.var_title(s, names[0]) if len(names) == 1 else ""
    opts["_xtitle"] = ""
    return g, graw


# ---------------------------------------------------------------------------
# graph pie
# ---------------------------------------------------------------------------

PIE_OPTS = {"over": 4, "missing": 4, "allcategories": 3, "cw": 2, "noclockwise": 9, "angle0": 6,
            "plabel": 3, "ptext": 3, "intensity": 5, "line": 4, "pie": 3, "sort": 4, "descending": 4,
            "nolabel": 7}


def pie(s: "Session", text: str) -> tuple[Graph, dict]:
    p = parse_standard(text)
    from ..lang.syntax import parse_options
    lopts, graw = split_layer_graph(p.options, PIE_OPTS)
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    mask = D.base_mask(s, p)
    w, mask = D.weights(s, p, mask, ("fweight", "aweight", "pweight"))
    over = None
    if isinstance(lopts.get("over"), str):
        over = s.expand_varlist(split_sub(lopts["over"])[0])[0]
    if not names and over is None:
        raise StataError(102, "too few variables specified")
    if len(names) > 1 and over is not None:
        raise StataError(198, "only one variable allowed with over()")   # VERIFY
    plabel = None
    for n, a in parse_options(p.options):
        if len(n) >= 3 and "plabel".startswith(n) and a:
            words = a.replace(",", " ").split()
            for k in ("percent", "sum", "name"):
                if any(k.startswith(x) and len(x) >= 3 for x in words[1:2]):
                    plabel = k
    opts = graph_options(graw)

    def slices(m: np.ndarray):
        if over is None:
            vals = [float((D.numeric(s, nm)[m & D.valid(D.numeric(s, nm))] *
                           (1 if w is None else w[m & D.valid(D.numeric(s, nm))])).sum()) for nm in names]
            return [D.var_title(s, nm) for nm in names], vals
        labels, vals = [], []
        for key, gm in D.groups(s, m, [over], missing=bool(lopts.get("missing"))):
            labels.append(D.group_title(s, [over], key))
            if names:
                xv = D.numeric(s, names[0])
                mm = gm & D.valid(xv)
                vals.append(float((xv[mm] * (1 if w is None else w[mm])).sum()))
            else:
                vals.append(float(gm.sum() if w is None else w[gm].sum()))
        return labels, vals

    g = Graph(family="pie", opts=opts)
    byo = opts.get("by")
    if byo:
        byvars = s.expand_varlist(" ".join(byo["vars"]))
        for key, gm in D.groups(s, mask, byvars):
            labels, vals = slices(gm)
            g.panels.append(Panel(D.group_title(s, byvars, key),
                                  [Layer("pie", extra={"labels": labels, "values": vals, "plabel": plabel})]))
    else:
        labels, vals = slices(mask)
        g.panels.append(Panel("", [Layer("pie", extra={"labels": labels, "values": vals, "plabel": plabel})]))
    return g, graw


# ---------------------------------------------------------------------------
# graph matrix
# ---------------------------------------------------------------------------

def matrix(s: "Session", text: str) -> tuple[Graph, dict]:
    p = parse_standard(text)
    lopts, graw = split_layer_graph(p.options, {**layer_spec(), "half": 4, "diagonal": 4, "maxes": 5,
                                                 "iscale": 6, "jitter": 6})
    names = s.expand_varlist(p.varlist)
    if len(names) < 2:
        raise StataError(102, "too few variables specified")
    mask = D.base_mask(s, p)
    cols = []
    for nm in names:
        if s.data.get(nm).is_string:
            raise StataError(109, "type mismatch")
    # VERIFY: Stata uses, for each pair, the obs. without missing in that pair
    for nm in names:
        mask &= D.valid(D.numeric(s, nm))
    for nm in names:
        cols.append(D.numeric(s, nm)[mask].tolist())
    labels = [D.var_title(s, nm) for nm in names]
    g = Graph(family="matrix", opts=graph_options(graw),
              panels=[Panel("", [Layer("matrix", extra={"names": labels, "cols": cols,
                                                         "half": bool(lopts.get("half"))})])])
    return g, graw


# ---------------------------------------------------------------------------
# graph combine
# ---------------------------------------------------------------------------

def combine(s: "Session", text: str) -> tuple[Graph, dict]:
    from .store import gph_path, load_gph, store
    c = find_top(text, ",")
    body, opt = (text, "") if c == -1 else (text[:c], text[c + 1:])
    _, graw = split_layer_graph(opt, {"cols": 3, "rows": 3, "colfirst": 4, "holes": 4, "iscale": 6,
                                      "imargin": 6, "altshrink": 6, "commonscheme": 6, "ycommon": 7,
                                      "xcommon": 7, "nocopies": 8})
    st = store(s)
    kids = []
    from ..lang.words import split_words
    for nm in split_words(body):
        if nm in st.graphs:
            kids.append(Graph.from_json(st.graphs[nm].to_json()))
        elif nm.endswith(".gph") or gph_path(nm).exists():
            kids.append(load_gph(gph_path(nm)))
        else:
            raise StataError(111, f"graph {nm} not found")   # VERIFY
    if not kids:
        raise StataError(198, "nothing to combine")   # VERIFY
    opts = graph_options(graw)
    for k in ("cols", "rows"):
        if isinstance(graw.get(k), str):
            opts[k] = int(graw[k])
    g = Graph(family="combine", opts=opts, children=kids)
    return g, graw


# ---------------------------------------------------------------------------
# twoway with an already built builder (used by tsline etc.)
# ---------------------------------------------------------------------------

def twoway_from_builder(s: "Session", b: Builder):
    return finish(s, b)
