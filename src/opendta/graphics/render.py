"""Drawing graphs with matplotlib.

Takes a Graph (model.py) and returns a matplotlib Figure, used both
in the Graph window of the interface and in `graph export`. The look follows
OpenDTA's own schemes (scheme.py); the general layout (titles at the
top, legend and notes at the bottom, by() panels in a grid) follows
Stata's as described in [G-3], without copying StataCorp graphical elements."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..lazy import lazy
from . import scheme as SC
from .model import Graph, Layer, Panel

if TYPE_CHECKING:
    from matplotlib.figure import Figure

# matplotlib is imported when the first graph is drawn (start-up time; see opendta/lazy.py)
_mfigure = lazy("matplotlib.figure")
_mticker = lazy("matplotlib.ticker")

# Stata's default size: 5.5 × 4 inches ([G-3] region_options)
DEFAULT_SIZE = (5.5, 4.0)

_HOLLOW = {"oh", "dh", "sh", "th", "circle_hollow", "diamond_hollow", "square_hollow", "triangle_hollow", "smcircle_hollow"}
_MSIZE = {"vtiny": 4, "tiny": 6, "vsmall": 9, "small": 14, "medsmall": 20, "medium": 26, "medlarge": 34,
          "large": 44, "vlarge": 60, "huge": 80, "vhuge": 110, "ehuge": 150}
_LWIDTH = {"none": 0, "vvthin": 0.3, "vthin": 0.5, "thin": 0.8, "medthin": 1.1, "medium": 1.4,
           "medthick": 1.8, "thick": 2.3, "vthick": 3.0, "vvthick": 4.0, "vvvthick": 5.0}
_TSIZE = {"vtiny": 5, "tiny": 6, "vsmall": 7, "small": 8, "medsmall": 9, "medium": 10, "medlarge": 11,
          "large": 13, "vlarge": 15, "huge": 18, "vhuge": 22}

FORMATS = {".png": ("png", "PNG"), ".pdf": ("pdf", "PDF"), ".svg": ("svg", "SVG"),
           ".eps": ("eps", "EPS"), ".ps": ("ps", "PS"), ".tif": ("tiff", "TIFF"),
           ".tiff": ("tiff", "TIFF"), ".jpg": ("jpg", "JPEG"), ".jpeg": ("jpg", "JPEG")}


# ---------------------------------------------------------------------------
# estilos
# ---------------------------------------------------------------------------

def _msymbol(style: dict, i: int, sch: dict) -> tuple[str, bool]:
    sym = str(style.get("msymbol", "")).strip().lower()
    if not sym:
        return sch["markers"][i % len(sch["markers"])], False
    return SC.SYMBOLS.get(sym, "o"), sym in _HOLLOW


def _color(style: dict, keys: tuple[str, ...], default: str) -> str:
    for k in keys:
        if style.get(k):
            return SC.color(str(style[k]), default)
    return default


def _lw(style: dict, default: float = 1.4) -> float:
    v = str(style.get("lwidth", "")).strip().lower()
    if not v:
        return default
    if v in _LWIDTH:
        return _LWIDTH[v]
    try:
        return float(v) * 1.4
    except ValueError:
        return default


def _ls(style: dict):
    v = str(style.get("lpattern", "")).strip().lower()
    return SC.PATTERNS.get(v, "-") if v else "-"


def _ms(style: dict, default: float = 22.0) -> float:
    v = str(style.get("msize", "")).strip().lower()
    if not v:
        return default
    if v in _MSIZE:
        return _MSIZE[v]
    try:
        return float(v) * 8
    except ValueError:
        return default


def _tsize(spec: dict, default: float) -> float:
    v = str(spec.get("size", "")).strip().lower() if isinstance(spec, dict) else ""
    return _TSIZE.get(v, default)


def _arr(v) -> np.ndarray:
    return np.asarray([np.nan if x is None else x for x in v], dtype=float)


# ---------------------------------------------------------------------------
# twoway layers
# ---------------------------------------------------------------------------

def draw_layer(ax, ly: Layer, i: int, sch: dict) -> list:
    """Draws a layer; returns the artists that go into the legend."""
    st = ly.style
    pal = sch["colors"]
    base_color = pal[int(st.get("_pal", i)) % len(pal)]
    col = _color(st, ("color", "mcolor", "lcolor"), base_color)
    lcol = _color(st, ("lcolor", "color"), base_color)
    fcol = _color(st, ("fcolor", "bcolor", "color"), base_color)
    x, y = _arr(ly.x), _arr(ly.y)
    horiz = bool(st.get("horizontal"))
    kind = ly.kind
    label = ly.label or None
    if kind in ("scatter", "scatteri"):
        mk, hollow = _msymbol(st, i, sch)
        if mk == "":
            return []
        sizes = ly.extra.get("sizes")
        s = _ms(st) if sizes is None else _arr(sizes)
        h = ax.scatter(x, y, s=s, marker=mk, label=label, linewidths=0.8,
                       facecolors="none" if hollow else col, edgecolors=col, zorder=3)
        labs = ly.extra.get("mlabel")
        if labs:
            for xi, yi, t in zip(x, y, labs):
                if np.isfinite(xi) and np.isfinite(yi):
                    ax.annotate(str(t), (xi, yi), xytext=(3, 3), textcoords="offset points",
                                fontsize=7, color=sch["text"])
        return [h]
    if kind in ("line", "connected", "lfit", "qfit", "function", "kdensity", "tsline", "mspline"):
        mk = ""
        if kind == "connected":
            mk, _ = _msymbol(st, i, sch)
        (h,) = ax.plot(y if horiz else x, x if horiz else y, color=lcol, linestyle=_ls(st), linewidth=_lw(st),
                       marker=mk or None, markersize=math.sqrt(_ms(st)), label=label, zorder=2)
        return [h]
    if kind == "area":
        base = float(st.get("base", 0) or 0)
        h = ax.fill_between(x, base, y, color=fcol, alpha=0.85, linewidth=0, label=label, zorder=1)
        return [h]
    if kind in ("bar", "histogram"):
        w = float(ly.extra.get("width", st.get("barwidth", 1) or 1))
        base = float(st.get("base", 0) or 0)
        edge = sch["background"] if kind == "histogram" else fcol
        if horiz:
            h = ax.barh(x, y - base, height=w, left=base, color=fcol, edgecolor=edge, linewidth=0.6,
                        label=label, zorder=1)
        else:
            h = ax.bar(x, y - base, width=w, bottom=base, color=fcol, edgecolor=edge, linewidth=0.6,
                       label=label, zorder=1)
        return [h]
    if kind in ("spike", "dropline"):
        base = float(st.get("base", 0) or 0)
        h = ax.vlines(x, base, y, colors=lcol, linewidth=_lw(st, 1.0), label=label)
        if kind == "dropline":
            mk, hollow = _msymbol(st, i, sch)
            ax.scatter(x, y, s=_ms(st), marker=mk, facecolors="none" if hollow else col, edgecolors=col)
        return [h]
    y2 = _arr(ly.extra.get("y2", []))
    if kind == "rarea":
        h = ax.fill_between(x, y, y2, color=fcol, alpha=0.35, linewidth=0, label=label, zorder=1)
        return [h]
    if kind == "rbar":
        w = float(st.get("barwidth", 0.6) or 0.6)
        h = ax.bar(x, y2 - y, bottom=y, width=w, color=fcol, label=label, zorder=1)
        return [h]
    if kind in ("rcap", "rspike", "rcapsym"):
        if horiz:
            h = ax.hlines(x, y, y2, colors=lcol, linewidth=_lw(st, 1.0), label=label)
        else:
            h = ax.vlines(x, y, y2, colors=lcol, linewidth=_lw(st, 1.0), label=label)
        if kind == "rcap":
            span = np.nanmax(x) - np.nanmin(x) if len(x) > 1 else 1.0
            cap = 0.012 * (span if span > 0 else 1.0)
            fn = ax.vlines if horiz else ax.hlines
            for yy in (y, y2):
                fn(yy, x - cap, x + cap, colors=lcol, linewidth=_lw(st, 1.0))
        if kind == "rcapsym":
            mk, _ = _msymbol(st, i, sch)
            ax.scatter(np.r_[x, x], np.r_[y, y2], s=_ms(st), marker=mk, color=lcol)
        return [h]
    if kind in ("rline", "rconnected"):
        mk = _msymbol(st, i, sch)[0] if kind == "rconnected" else None
        (h,) = ax.plot(x, y, color=lcol, linestyle=_ls(st), linewidth=_lw(st), marker=mk, label=label)
        ax.plot(x, y2, color=lcol, linestyle=_ls(st), linewidth=_lw(st), marker=mk)
        return [h]
    if kind == "rscatter":
        mk, hollow = _msymbol(st, i, sch)
        h = ax.scatter(np.r_[x, x], np.r_[y, y2], s=_ms(st), marker=mk, label=label,
                       facecolors="none" if hollow else col, edgecolors=col)
        return [h]
    return []


# ---------------------------------------------------------------------------
# eixos
# ---------------------------------------------------------------------------

def _apply_axis(ax, which: str, o: dict, sch: dict, value_labels: dict | None) -> None:
    lab = o.get(which + "label", {})
    sc = o.get(which + "scale", {})
    axis = ax.xaxis if which == "x" else ax.yaxis
    if sc.get("log"):
        (ax.set_xscale if which == "x" else ax.set_yscale)("log")
    if sc.get("range"):
        lo, hi = sc["range"]
        cur = ax.get_xlim() if which == "x" else ax.get_ylim()
        lo, hi = min(lo, cur[0]), max(hi, cur[1])
        (ax.set_xlim if which == "x" else ax.set_ylim)(lo, hi)
    if sc.get("reverse"):
        (ax.invert_xaxis if which == "x" else ax.invert_yaxis)()
    if sc.get("off"):
        axis.set_visible(False)
    if lab.get("none"):
        axis.set_ticks([])
    elif lab.get("ticks"):
        vals = [t[0] for t in lab["ticks"]]
        texts = []
        for v, t in lab["ticks"]:
            if t is not None:
                texts.append(t)
            elif lab.get("valuelabel") and value_labels and float(v).is_integer() and int(v) in value_labels:
                texts.append(value_labels[int(v)])
            else:
                texts.append(_num(v))
        axis.set_ticks(vals)
        axis.set_ticklabels(texts)
    elif lab.get("auto"):
        axis.set_major_locator(_mticker.MaxNLocator(nbins=max(1, lab["auto"] - 1)))
    elif lab.get("valuelabel") and value_labels:
        keys = sorted(value_labels)
        axis.set_ticks(keys)
        axis.set_ticklabels([value_labels[k] for k in keys])
    if "angle" in lab:
        for t in axis.get_ticklabels():
            t.set_rotation(lab["angle"])
    grid = lab.get("grid", which == "y")
    if grid:
        ax.grid(True, axis=which, color=sch["grid"], linewidth=0.7, zorder=0)
    else:
        ax.grid(False, axis=which)
    for v in o.get(which + "line", {}).get("values", []) if which in ("x", "y") else []:
        spec = o.get(which + "line", {})
        c = SC.color(spec.get("lcolor"), sch["axis"])
        (ax.axvline if which == "x" else ax.axhline)(v, color=c, linestyle=SC.PATTERNS.get(spec.get("lpattern", ""), "-"),
                                                     linewidth=0.9, zorder=1)


def _num(v: float) -> str:
    if float(v).is_integer():
        return str(int(v))
    t = f"{v:.6g}"
    return t[1:] if t.startswith("0.") else ("-" + t[2:] if t.startswith("-0.") else t)


def _style_axes(ax, sch: dict) -> None:
    ax.set_facecolor(sch["plotregion"])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(sch["axis"])
    ax.tick_params(colors=sch["axis"], labelcolor=sch["text"], labelsize=8)
    ax.set_axisbelow(True)


def _join(spec) -> str:
    if not spec:
        return ""
    return "\n".join(spec.get("lines", [])) if isinstance(spec, dict) else str(spec)


# ---------------------------------------------------------------------------
# families
# ---------------------------------------------------------------------------

def _draw_twoway_panel(ax, panel: Panel, g: Graph, sch: dict) -> list:
    handles = []
    for i, ly in enumerate(panel.layers):
        hs = draw_layer(ax, ly, i, sch)
        if hs and ly.label:
            handles.append((hs[0], ly.label))
    if g.opts.get("text"):
        for yv, xv, t in g.opts["text"]:
            ax.text(xv, yv, t, fontsize=8, color=sch["text"], ha="center", va="center")
    return handles


def _draw_bars(ax, ly: Layer, sch: dict, horizontal: bool) -> list:
    cats = ly.extra.get("cats", [])
    series = ly.extra.get("series", [])
    stack = ly.extra.get("stack", False)
    nser = max(1, len(series))
    pos = np.arange(len(cats), dtype=float)
    width = 0.8 / (1 if stack else nser)
    bottom = np.zeros(len(cats))
    handles = []
    for j, ser in enumerate(series):
        vals = _arr(ser["values"])
        c = sch["colors"][j % len(sch["colors"])]
        offs = pos if stack else pos - 0.4 + width * (j + 0.5)
        if horizontal:
            h = ax.barh(offs, np.nan_to_num(vals), height=width, left=bottom if stack else 0, color=c, zorder=2)
        else:
            h = ax.bar(offs, np.nan_to_num(vals), width=width, bottom=bottom if stack else 0, color=c, zorder=2)
        if ly.extra.get("blabel"):
            ax.bar_label(h, labels=[_num(round(v, 4)) if np.isfinite(v) else "" for v in vals], fontsize=7,
                         color=sch["text"])
        if stack:
            bottom = bottom + np.nan_to_num(vals)
        handles.append((h, ser.get("label", "")))
    if horizontal:
        ax.set_yticks(pos)
        ax.set_yticklabels(cats)
        ax.invert_yaxis()
        ax.grid(True, axis="x", color=sch["grid"], linewidth=0.7, zorder=0)
    else:
        ax.set_xticks(pos)
        ax.set_xticklabels(cats)
        ax.grid(True, axis="y", color=sch["grid"], linewidth=0.7, zorder=0)
    return handles if len(series) > 1 else []


def _orientation_kw(horizontal: bool) -> dict:
    """matplotlib 3.10 replaced vert= with orientation= in bxp."""
    import matplotlib
    major, minor = (int(p) for p in matplotlib.__version__.split(".")[:2])
    if (major, minor) >= (3, 10):
        return {"orientation": "horizontal" if horizontal else "vertical"}
    return {"vert": not horizontal}


def _draw_dots(ax, ly: Layer, sch: dict) -> list:
    """graph dot: one dotted line per category and one marker per series."""
    cats = ly.extra.get("cats", [])
    series = ly.extra.get("series", [])
    pos = np.arange(len(cats), dtype=float)
    for p_ in pos:
        ax.axhline(p_, color=sch["grid"], linestyle=":", linewidth=0.9, zorder=0)
    handles = []
    for j, ser in enumerate(series):
        vals = _arr(ser["values"])
        c = sch["colors"][j % len(sch["colors"])]
        h = ax.scatter(vals, pos, s=30, color=c, marker=sch["markers"][j % len(sch["markers"])], zorder=3)
        handles.append((h, ser.get("label", "")))
    ax.set_yticks(pos)
    ax.set_yticklabels(cats)
    ax.invert_yaxis()
    return handles if len(series) > 1 else []


def _draw_boxes(ax, ly: Layer, sch: dict, horizontal: bool) -> list:
    cats = ly.extra.get("cats", [])
    series = ly.extra.get("series", [])
    nser = max(1, len(series))
    width = 0.7 / nser
    handles = []
    for j, ser in enumerate(series):
        c = sch["colors"][j % len(sch["colors"])]
        stats = []
        positions = []
        for k, b in enumerate(ser["stats"]):
            if b is None:
                continue
            stats.append({"med": b["med"], "q1": b["q1"], "q3": b["q3"], "whislo": b["lo"], "whishi": b["hi"],
                          "fliers": b["out"]})
            positions.append(k - 0.35 + width * (j + 0.5))
        if not stats:
            continue
        orient = _orientation_kw(horizontal)
        r = ax.bxp(stats, positions=positions, widths=width * 0.85, patch_artist=True, **orient,
                   showfliers=True, manage_ticks=False,
                   boxprops={"facecolor": c, "edgecolor": c},
                   medianprops={"color": sch["background"], "linewidth": 1.5},
                   whiskerprops={"color": c}, capprops={"color": c},
                   flierprops={"marker": "o", "markerfacecolor": "none", "markeredgecolor": c, "markersize": 4})
        handles.append((r["boxes"][0], ser.get("label", "")))
    pos = np.arange(len(cats))
    if horizontal:
        ax.set_yticks(pos)
        ax.set_yticklabels(cats)
        ax.grid(True, axis="x", color=sch["grid"], linewidth=0.7, zorder=0)
    else:
        ax.set_xticks(pos)
        ax.set_xticklabels(cats)
        ax.grid(True, axis="y", color=sch["grid"], linewidth=0.7, zorder=0)
    if len(cats) == 1 and cats[0] == "":
        (ax.set_yticks if horizontal else ax.set_xticks)([])
    return handles if len(series) > 1 else []


def _draw_pie(ax, ly: Layer, sch: dict) -> list:
    vals = _arr(ly.extra.get("values", []))
    labels = ly.extra.get("labels", [])
    ok = np.isfinite(vals) & (vals > 0)
    vals = vals[ok]
    labels = [lab for lab, k in zip(labels, ok) if k]
    pl = ly.extra.get("plabel")
    auto = None
    if pl == "percent":
        auto = "%1.1f%%"
    elif pl == "sum":
        total = vals.sum()
        auto = lambda p: _num(round(p * total / 100, 4))   # noqa: E731
    wedges, *_ = ax.pie(vals, colors=[sch["colors"][k % len(sch["colors"])] for k in range(len(vals))],
                        autopct=auto, startangle=90, counterclock=False,
                        wedgeprops={"edgecolor": sch["background"], "linewidth": 1},
                        textprops={"fontsize": 8, "color": sch["text"]})
    ax.set_aspect("equal")
    return list(zip(wedges, labels))


def _draw_matrix(fig, ly: Layer, sch: dict) -> None:
    names = ly.extra.get("names", [])
    cols = [_arr(c) for c in ly.extra.get("cols", [])]
    half = ly.extra.get("half", False)
    k = len(names)
    gs = fig.add_gridspec(k, k, hspace=0.08, wspace=0.08)
    for r in range(k):
        for c in range(k):
            if half and c > r:
                continue
            ax = fig.add_subplot(gs[r, c])
            _style_axes(ax, sch)
            ax.tick_params(labelsize=6)
            if r == c:
                ax.text(0.5, 0.5, names[r], ha="center", va="center", fontsize=9, color=sch["text"],
                        transform=ax.transAxes)
                ax.set_xticks([])
                ax.set_yticks([])
                for sp in ax.spines.values():
                    sp.set_visible(True)
                    sp.set_color(sch["grid"])
                continue
            ax.scatter(cols[c], cols[r], s=6, color=sch["colors"][0], linewidths=0)
            if r != k - 1:
                ax.set_xticklabels([])
            if c != 0:
                ax.set_yticklabels([])


# ---------------------------------------------------------------------------
# montagem
# ---------------------------------------------------------------------------

def draw(container, g: Graph, sch: dict) -> None:
    """Draws `g` in a Figure or SubFigure."""
    o = g.opts
    container.set_facecolor(sch["background"])
    tcolor = sch["text"]
    title = _join(o.get("title"))
    subtitle = _join(o.get("subtitle"))
    head = title + (("\n" + subtitle) if subtitle and title else subtitle)
    if head:
        container.suptitle(head, fontsize=_tsize(o.get("title", {}), 12), color=tcolor)
    if g.family == "combine":
        kids = g.children
        n = len(kids)
        cols = int(o.get("cols") or 0) or (math.ceil(math.sqrt(n)) if n else 1)
        rows = int(o.get("rows") or 0) or max(1, math.ceil(n / cols))
        if o.get("rows") and not o.get("cols"):
            cols = max(1, math.ceil(n / rows))
        subs = container.subfigures(rows, cols, squeeze=False)
        for k, child in enumerate(kids):
            sub = subs[k // cols][k % cols]
            draw(sub, child, SC.get(child.opts.get("scheme")) if child.opts.get("scheme") else sch)
        for k in range(n, rows * cols):
            subs[k // cols][k % cols].set_facecolor(sch["background"])
        _notes(container, o, sch)
        return
    panels = g.panels or [Panel()]
    if g.family == "matrix":
        _draw_matrix(container, panels[0].layers[0], sch)
        _notes(container, o, sch)
        return
    byo = o.get("by")
    n = len(panels)
    if byo or n > 1:
        cols = int(byo.get("cols") or 0) if byo else 0
        rows = int(byo.get("rows") or 0) if byo else 0
        if not cols:
            cols = math.ceil(n / rows) if rows else (math.ceil(math.sqrt(n)) if n > 1 else 1)
        rows = max(1, math.ceil(n / cols))
        share = g.family not in ("pie",)
        axes = container.subplots(rows, cols, squeeze=False, sharex=share, sharey=share)
    else:
        rows = cols = 1
        axes = [[container.add_subplot(1, 1, 1)]]
    handles: list = []
    horizontal = g.family in ("hbar", "hbox", "dot")
    for k in range(rows * cols):
        ax = axes[k // cols][k % cols]
        if k >= n:
            ax.set_visible(False)
            continue
        p = panels[k]
        _style_axes(ax, sch)
        if g.family == "pie":
            ax.axis("off")
        if p.title and (byo or n > 1):
            ax.set_title(p.title, fontsize=9, color=tcolor, backgroundcolor=sch["grid"])
        elif subtitle and not title and n == 1:
            pass
        hs: list = []
        for ly in p.layers:
            if ly.kind == "bars" and g.family == "dot":
                hs += _draw_dots(ax, ly, sch)
            elif ly.kind == "bars":
                hs += _draw_bars(ax, ly, sch, horizontal)
            elif ly.kind == "boxes":
                hs += _draw_boxes(ax, ly, sch, horizontal)
            elif ly.kind == "pie":
                hs += _draw_pie(ax, ly, sch)
        if g.family == "twoway":
            hs = _draw_twoway_panel(ax, p, g, sch)
            vl = o.get("_xvaluelabels")
            _apply_axis(ax, "x", o, sch, {int(a): b for a, b in vl.items()} if vl else None)
            vl = o.get("_yvaluelabels")
            _apply_axis(ax, "y", o, sch, {int(a): b for a, b in vl.items()} if vl else None)
        elif g.family in ("bar", "hbar", "box", "hbox", "dot"):
            vax = "x" if horizontal else "y"
            _apply_axis(ax, vax, {**o, vax + "label": {**o.get("ylabel", {}), "grid": True}}, sch, None)
        if not handles:
            handles = hs
        if g.family != "pie":
            last_row = (k // cols == rows - 1) or (k + cols >= n)
            if last_row or n == 1:
                xt = _join(o.get("xtitle")) if "xtitle" in o else o.get("_xtitle", "")
                if horizontal:
                    xt = _join(o.get("ytitle")) if "ytitle" in o else o.get("_ytitle", "")
                ax.set_xlabel(xt, fontsize=9, color=tcolor)
            if k % cols == 0:
                yt = _join(o.get("ytitle")) if "ytitle" in o else o.get("_ytitle", "")
                if horizontal:
                    yt = ""
                ax.set_ylabel(yt, fontsize=9, color=tcolor)
    leg = o.get("legend", {})
    show = (len(handles) > 1 or (g.family == "pie" and handles)) and not leg.get("off")
    if leg.get("on") and handles:
        show = True
    if show:
        order = leg.get("order")
        labels_over = leg.get("labels", {})
        items = []
        if order:
            for idx, txt in order:
                if 1 <= idx <= len(handles):
                    h, lab = handles[idx - 1]
                    items.append((h, txt if txt is not None else labels_over.get(str(idx), lab)))
        else:
            for idx, (h, lab) in enumerate(handles, start=1):
                items.append((h, labels_over.get(str(idx), lab)))
        ncol = leg.get("cols") or (math.ceil(len(items) / leg["rows"]) if leg.get("rows") else min(2, len(items)))
        L = container.legend([h for h, _ in items], [t for _, t in items], loc="outside lower center",
                             ncol=ncol, fontsize=8, frameon=True, title=leg.get("title"))
        L.get_frame().set_edgecolor(sch["grid"])
        L.get_frame().set_facecolor(sch["background"])
        for t in L.get_texts():
            t.set_color(tcolor)
    _notes(container, o, sch)


def _notes(container, o: dict, sch: dict) -> None:
    parts = [_join(o.get("note")), _join(o.get("caption"))]
    text = "\n".join(p for p in parts if p)
    if text:
        container.supxlabel(text, x=0.01, ha="left", fontsize=7, color=sch["text"])


def figure(g: Graph, default_scheme: str = "opendta") -> Figure:
    sch = SC.get(g.opts.get("scheme") or default_scheme)
    w = float(g.opts.get("xsize", DEFAULT_SIZE[0]))
    h = float(g.opts.get("ysize", DEFAULT_SIZE[1]))
    fig = _mfigure.Figure(figsize=(w, h), layout="constrained")
    fig.patch.set_facecolor(sch["background"])
    draw(fig, g, sch)
    return fig


def export(g: Graph, path, fmt: str, *, width: int | None = None, height: int | None = None,
           default_scheme: str = "opendta") -> None:
    fig = figure(g, default_scheme)
    dpi = 150
    w, h = fig.get_size_inches()
    if width:
        dpi = width / w
    elif height:
        dpi = height / h
    fig.savefig(path, format=fmt, dpi=dpi, facecolor=fig.get_facecolor())
