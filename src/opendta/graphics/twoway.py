"""graph twoway: syntax parsing and assembly of the layers ([G-2] graph twoway).

Accepted forms:
    twoway (scatter y x) (lfit y x), title("...")
    twoway scatter y x || lfit y x ||, title("...")
    scatter y x, by(foreign)
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from ..lang.syntax import parse_standard
from ..lang.words import parse_numlist, split_words, strip_outer_quotes
from . import data as D
from .model import Graph, Layer, Panel
from .opts import graph_options, split_layer_graph

if TYPE_CHECKING:
    from ..session import Session

PLOTTYPES = {
    "scatter": 2, "line": 4, "connected": 4, "lfit": 4, "qfit": 4, "lfitci": 6, "qfitci": 6,
    "area": 4, "bar": 3, "spike": 5, "dropline": 8, "rcap": 4, "rcapsym": 7, "rarea": 5, "rbar": 4,
    "rspike": 6, "rline": 5, "rconnected": 5, "rscatter": 4, "function": 4, "histogram": 4,
    "kdensity": 4, "scatteri": 8, "tsline": 6, "lowess": 3,
}
RANGE_PLOTS = ("rcap", "rcapsym", "rarea", "rbar", "rspike", "rline", "rconnected", "rscatter")


def plottype(word: str) -> str | None:
    w = word.lower()
    if w in PLOTTYPES:
        return w
    for full, m in PLOTTYPES.items():
        if len(w) >= m and full.startswith(w):
            return full
    return None


# ---------------------------------------------------------------------------
# splitting the command into plots
# ---------------------------------------------------------------------------

def split_plots(text: str) -> tuple[list[str], str]:
    """Splits into plots (by `||` or parentheses) and also returns the final
    global options (after `)` or an empty `||`)."""
    pieces: list[str] = []
    glob = ""
    buf: list[str] = []
    i, n = 0, len(text)
    depth = 0

    def push() -> None:
        t = "".join(buf).strip()
        buf.clear()
        if t:
            pieces.append(t)

    while i < n:
        c = text[i]
        if c == '"' or text.startswith('`"', i):
            if c == '"':
                j = text.find('"', i + 1)
                j = n if j == -1 else j + 1
            else:
                d, j = 1, i + 2
                while j < n and d:
                    if text.startswith('`"', j):
                        d += 1
                        j += 2
                    elif text.startswith("\"'", j):
                        d -= 1
                        j += 2
                    else:
                        j += 1
            buf.append(text[i:j])
            i = j
            continue
        if depth == 0 and text.startswith("||", i):
            push()
            i += 2
            continue
        if depth == 0 and c == "(" and not "".join(buf).strip():
            depth, j = 1, i + 1
            while j < n and depth:
                if text[j] == '"':
                    k = text.find('"', j + 1)
                    j = n if k == -1 else k + 1
                    continue
                if text[j] == "(":
                    depth += 1
                elif text[j] == ")":
                    depth -= 1
                j += 1
            if depth:
                raise StataError(198, "unmatched parenthesis")   # VERIFY
            pieces.append(text[i + 1:j - 1].strip())
            i = j
            continue
        if depth == 0 and c == "," and not "".join(buf).strip():
            glob = text[i + 1:]
            break
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        buf.append(c)
        i += 1
    push()
    # a piece that starts with a comma is a global option (scatter y x ||, title())
    out = []
    for p in pieces:
        if p.startswith(","):
            glob = (glob + " " + p[1:]).strip()
        else:
            out.append(p)
    return out, glob.strip()


# ---------------------------------------------------------------------------
# montagem
# ---------------------------------------------------------------------------

class Builder:
    """Joins the layers of all plots; with by(), builds one panel per group."""

    def __init__(self, s: "Session"):
        self.s = s
        self.specs: list[tuple[str, object, dict, np.ndarray]] = []   # (type, parsed, layer opts, mask)
        self.gopts: dict = {}
        self.ytitles: list[str] = []
        self.xtitles: list[str] = []
        self.xvl: dict = {}
        self.yvl: dict = {}

    def add(self, spec: str) -> None:
        words = spec.split(None, 1)
        if not words:
            return
        pt = plottype(words[0])
        if pt is None:
            raise StataError(198, f"{words[0]} is not a twoway plot type")   # VERIFY
        rest = words[1] if len(words) > 1 else ""
        self.add_typed(pt, rest)

    def add_typed(self, pt: str, rest: str) -> None:
        s = self.s
        if pt == "scatteri":
            from ..lang.syntax import find_top
            c = find_top(rest, ",")
            body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
            lopts, g = split_layer_graph(opt, layer_spec())
            self.gopts.update(g)
            nums = [w for w in split_words(body) if not (w.startswith("(") and w.endswith(")"))]
            vals = parse_numlist(" ".join(nums))
            if len(vals) % 2:
                raise StataError(198, "invalid syntax")
            self.specs.append(("scatteri", (vals[0::2], vals[1::2]), lopts, None))
            return
        if pt == "function":
            self._function(rest)
            return
        p = parse_standard(rest)
        lopts, g = split_layer_graph(p.options, layer_spec())
        self.gopts.update(g)
        mask = D.base_mask(s, p)
        names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
        for nm in names:
            if s.data.get(nm).is_string:
                raise StataError(109, "type mismatch")   # VERIFY
        if pt == "tsline":
            from ..estimation.tsops import ts_info
            ti = ts_info(s.data)
            if ti is None:
                raise StataError(111, "time variable not set, use tsset varname ...")
            names = names + [ti.tvar]
        need = {"histogram": 1, "kdensity": 1}.get(pt, 3 if pt in RANGE_PLOTS else 2)
        if len(names) < need:
            raise StataError(102, "too few variables specified")
        if pt in ("histogram", "kdensity", "lfit", "qfit", "lfitci", "qfitci", "lowess") or pt in RANGE_PLOTS:
            maxv = {"histogram": 1, "kdensity": 1}.get(pt, 3 if pt in RANGE_PLOTS else 2)
            if len(names) > maxv:
                raise StataError(103, "too many variables specified")
        w, mask = D.weights(s, p, mask, ("fweight", "aweight", "pweight", "iweight"))
        self.specs.append((pt, (names, w, p), lopts, mask))

    def _function(self, rest: str) -> None:
        from ..lang.syntax import find_top
        c = find_top(rest, ",")
        body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
        lopts, g = split_layer_graph(opt, layer_spec())
        self.gopts.update(g)
        body = body.strip()
        m = re.match(r"^([A-Za-z_]\w*)\s*=(?!=)\s*(.*)$", body)
        yname, expr = (m.group(1), m.group(2)) if m else ("y", body)
        self.specs.append(("function", (yname, expr), lopts, None))

    # -- computing the layers for a group mask -------------------------
    def layers(self, group: np.ndarray | None) -> list[Layer]:
        s = self.s
        out: list[Layer] = []
        for pt, payload, lopts, mask in self.specs:
            style = {k: (strip_outer_quotes(v) if isinstance(v, str) else v) for k, v in lopts.items()}
            if pt == "scatteri":
                ys, xs = payload
                out.append(Layer("scatteri", list(xs), list(ys), {}, style, ""))
                continue
            if pt == "function":
                out.append(self._function_layer(payload, style))
                continue
            names, w, p = payload
            m = mask if group is None else (mask & group)
            out.extend(self._layers_for(pt, names, w, m, style))
        return out

    def _function_layer(self, payload, style) -> Layer:
        yname, expr = payload
        rng = style.get("range")
        lo, hi = 0.0, 1.0
        if isinstance(rng, str) and rng.strip():
            vals = rng.split()
            if len(vals) == 2 and all(D_isnum(v) for v in vals):
                lo, hi = float(vals[0]), float(vals[1])
            elif len(vals) == 1 and self.s.data.has(vals[0]):
                x = D.numeric(self.s, vals[0])
                x = x[D.valid(x)]
                lo, hi = float(x.min()), float(x.max())
        n = int(style.get("n", 300) or 300)
        grid = np.linspace(lo, hi, n)
        y = eval_on_grid(self.s, expr, grid)
        self.xtitles.append("x")
        self.ytitles.append(yname)
        return Layer("function", grid.tolist(), y.tolist(), {}, style, f"{yname} = {expr}")

    def _layers_for(self, pt: str, names: list[str], w, m: np.ndarray, style: dict) -> list[Layer]:
        s = self.s
        out: list[Layer] = []
        if pt in ("histogram", "kdensity"):
            v = names[0]
            x = D.numeric(s, v)
            mm = m & D.valid(x)
            ww = None if w is None else w[mm]
            self.xtitles.append(D.var_title(s, v))
            if pt == "histogram":
                h = D.hist_bins(x[mm], ww, hist_opts(style))
                self.ytitles.append(hist_ytitle(style))
                for k in ("bin", "width", "start", "discrete", "density", "fraction", "frequency", "percent"):
                    style.pop(k, None)
                out.append(Layer("histogram", h["centers"], h["heights"], {"width": h["width"]}, style,
                                 D.var_title(s, v)))
            else:
                kd = D.kdensity(x[mm], ww, kernel=D.kernel_name(style.get("kernel")),
                                bw=float(style["bwidth"]) if style.get("bwidth") else None,
                                n=int(style["n"]) if style.get("n") else None)
                self.ytitles.append("Density")   # VERIFY
                out.append(Layer("kdensity", kd["x"].tolist(), kd["d"].tolist(), {}, style, "kdensity " + v))
            return out
        xname = names[-1] if pt not in RANGE_PLOTS else names[2]
        x = D.numeric(s, xname)
        self.xtitles.append(D.var_title(s, xname))
        vl = D.value_labels(s, xname)
        if vl:
            self.xvl.update(vl)
        if pt in RANGE_PLOTS:
            y1, y2 = D.numeric(s, names[0]), D.numeric(s, names[1])
            mm = m & D.valid(x) & D.valid(y1) & D.valid(y2)
            order = np.argsort(x[mm], kind="stable") if pt in ("rarea", "rline", "rconnected") else slice(None)
            out.append(Layer(pt, x[mm][order].tolist(), y1[mm][order].tolist(),
                             {"y2": y2[mm][order].tolist()}, style, f"{names[0]}/{names[1]}"))
            self.ytitles.append("")
            return out
        if pt == "lowess":
            from .postplots import lowess_fit
            y = D.numeric(s, names[0])
            mm = m & D.valid(x) & D.valid(y)
            bw = float(style.pop("bwidth")) if isinstance(style.get("bwidth"), str) else 0.8
            fit = lowess_fit(x[mm], y[mm], bw, bool(style.pop("mean", False)))
            o = np.argsort(x[mm], kind="stable")
            out.append(Layer("line", x[mm][o].tolist(), fit[o].tolist(), {}, style,
                             f"lowess {names[0]} {names[1]}"))
            self.ytitles.append(D.var_title(s, names[0]))
            return out
        if pt in ("lfit", "qfit", "lfitci", "qfitci"):
            y = D.numeric(s, names[0])
            mm = m & D.valid(x) & D.valid(y)
            xv, yv = x[mm], y[mm]
            if len(xv) < 2:
                raise StataError(2000, "no observations")
            deg = 1 if pt.startswith("l") else 2
            npts = int(style.get("n", 0) or 0) or (2 if deg == 1 and not pt.endswith("ci") else 100)
            grid = np.linspace(float(xv.min()), float(xv.max()), max(npts, 2))
            ci = pt.endswith("ci")
            lev = float(style.get("level", self.s.settings.get("level", 95)))
            fit, lo, hi = D.poly_fit(xv, yv, deg, grid, ci, lev, None if w is None else w[mm])
            if ci:
                band = dict(style)
                band.setdefault("color", "gs12")
                out.append(Layer("rarea", grid.tolist(), lo.tolist(), {"y2": hi.tolist()}, band,
                                 f"{lev:g}% CI"))
            out.append(Layer("lfit" if deg == 1 else "qfit", grid.tolist(), fit.tolist(), {}, style,
                             "Fitted values"))
            self.ytitles.append("Fitted values" if not ci else "")
            return out
        ys = names[:-1]
        for yn in ys:
            y = D.numeric(s, yn)
            mm = m & D.valid(x) & D.valid(y)
            xv, yv = x[mm], y[mm]
            if style.get("sort") or pt == "tsline":
                o = np.argsort(xv, kind="stable")
                xv, yv = xv[o], yv[o]
            extra: dict = {}
            if style.get("mlabel") and isinstance(style["mlabel"], str):
                lv = style["mlabel"].strip()
                var = s.data.get(s.expand_varlist(lv)[0])
                extra["mlabel"] = [D.value_text(s, var.name, var.raw[i] if var.is_string else var.data[i])
                                   for i in np.flatnonzero(mm)]
            if w is not None and pt == "scatter":
                ww = w[mm]
                extra["sizes"] = (22 * ww / ww.mean()).tolist() if ww.mean() > 0 else None
            kind = "line" if pt == "tsline" else pt
            out.append(Layer(kind, xv.tolist(), yv.tolist(), extra, style, D.var_title(s, yn)))
            self.ytitles.append(D.var_title(s, yn))
            yl = D.value_labels(s, yn)
            if yl:
                self.yvl.update(yl)
        return out


def D_isnum(t: str) -> bool:
    try:
        float(t)
        return True
    except ValueError:
        return False


def eval_on_grid(s: "Session", expr: str, grid: np.ndarray, var: str = "x") -> np.ndarray:
    """Evaluates `expr` with `x` taking each point of the grid (twoway function)."""
    from ..commands._util import eval_vector
    from ..core.dataset import Dataset, Variable
    saved, saved_by = s.data, s.by_groups
    tmp = Dataset()
    tmp.set_obs(len(grid))
    tmp.add(Variable(var, "double", grid.astype(np.float64)))
    s.data, s.by_groups = tmp, None
    try:
        v = eval_vector(s, expr)
    finally:
        s.data, s.by_groups = saved, saved_by
    y = np.broadcast_to(np.asarray(v, dtype=np.float64), grid.shape).copy()
    y[y >= D.SYS] = np.nan
    return y


def hist_opts(style: dict) -> dict:
    o: dict = {}
    for k in ("bin", "width", "start"):
        if k in style and style[k] is not True:
            o[k] = float(style[k])
    if style.get("discrete"):
        o["discrete"] = True
    for k in ("fraction", "frequency", "percent", "density"):
        if style.get(k):
            o["scale"] = k
    return o


def hist_ytitle(style: dict) -> str:
    for k, t in (("fraction", "Fraction"), ("frequency", "Frequency"), ("percent", "Percent")):
        if style.get(k):
            return t
    return "Density"


HIST_LAYER = {"mean": 4, "bin": 3, "width": 1, "start": 5, "discrete": 4, "density": 4, "fraction": 4,
              "frequency": 4, "percent": 4, "kernel": 1, "bwidth": 2, "gap": 3}


def layer_spec() -> dict:
    from .opts import LAYER_OPTS
    return {**LAYER_OPTS, **HIST_LAYER}


def build(s: "Session", text: str, *, default_type: str | None = None) -> tuple[Graph, dict]:
    """Builds the Graph; also returns the graph options as text (name,
    saving, nodraw stay with the calling command)."""
    if default_type:
        text = default_type + " " + text
    pieces, glob = split_plots(text)
    if not pieces:
        raise StataError(198, "invalid syntax")
    b = Builder(s)
    for piece in pieces:
        b.add(piece)
    if glob:
        lopts, g = split_layer_graph(glob, layer_spec())
        if lopts:
            # layer options after the global `,`: apply to all plots
            for k, spec in enumerate(b.specs):
                b.specs[k] = (spec[0], spec[1], {**lopts, **spec[2]}, spec[3])
        b.gopts.update(g)
    return finish(s, b)


def finish(s: "Session", b: Builder) -> tuple[Graph, dict]:
    graph_raw = b.gopts
    opts = graph_options(graph_raw)
    g = Graph(family="twoway", opts=opts)
    byo = opts.get("by")
    if byo:
        byvars = s.expand_varlist(" ".join(byo["vars"]))
        byo["vars"] = byvars
        # groups with no observations in the sample of any plot do not become a panel
        allmask = np.ones(s.data.nobs, dtype=bool)
        for key, gm in D.groups(s, allmask, byvars, missing=bool(byo.get("missing"))):
            layers = b.layers(gm)
            if not any(len(ly.x) for ly in layers):
                continue
            g.panels.append(Panel(D.group_title(s, byvars, key), layers))
        if byo.get("total"):
            g.panels.append(Panel("Total", b.layers(None)))
        if not g.panels:
            raise StataError(2000, "no observations")
        labs = ", ".join(D.var_title(s, v) for v in byvars)
        note = opts.get("note", {}).get("lines") if "note" in opts else None
        if note is None and "note" not in byo:
            opts["note"] = {"lines": [f"Graphs by {labs}"]}
        elif byo.get("note"):
            opts["note"] = {"lines": byo["note"]}
        if byo.get("title"):
            opts["title"] = {"lines": byo["title"]}
    else:
        layers = b.layers(None)
        if not any(len(ly.x) for ly in layers):
            raise StataError(2000, "no observations")
        g.panels.append(Panel("", layers))
    # automatic axis titles: a single y → its label; several → blank
    yt = [t for t in dict.fromkeys(b.ytitles) if t is not None]
    xt = list(dict.fromkeys(b.xtitles))
    opts["_ytitle"] = yt[0] if len(yt) == 1 else ""
    opts["_xtitle"] = xt[0] if xt else ""
    if b.xvl:
        opts["_xvaluelabels"] = {str(k): v for k, v in b.xvl.items()}
    if b.yvl:
        opts["_yvaluelabels"] = {str(k): v for k, v in b.yvl.items()}
    # with a single layer, no legend (Stata only shows a legend with 2 or more)
    return g, graph_raw
