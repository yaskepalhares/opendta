"""Leitura das opções de gráfico ([G-3] options).

As opções são separadas em dois grupos:
  - opções da camada (cores, símbolos, padrões de linha...), que valem só
    para o plot em que aparecem;
  - opções do gráfico (títulos, eixos, legenda, by(), name()...), que valem
    para o gráfico inteiro, em qualquer plot onde estejam escritas (como no
    Stata, onde `twoway (scatter y x, title("a"))` dá título ao gráfico).

O resultado é sempre um dicionário simples (JSON), para caber no .gph.
"""

from __future__ import annotations

import re

from ..core.errors import StataError
from ..lang.syntax import find_top, parse_options
from ..lang.words import parse_numlist, split_words, strip_outer_quotes

# nome completo -> abreviação mínima
LAYER_OPTS = {
    "mcolor": 3, "msymbol": 2, "msize": 3, "mfcolor": 3, "mlcolor": 3, "mlabel": 4, "mlabposition": 6,
    "mlabcolor": 6, "mlabsize": 6, "lcolor": 3, "lpattern": 3, "lwidth": 3, "lstyle": 3,
    "color": 5, "fcolor": 3, "fintensity": 5, "bcolor": 3, "barwidth": 4, "base": 4,
    "sort": 4, "connect": 1, "cmissing": 5, "jitter": 6, "horizontal": 5, "vertical": 4,
    "range": 1, "n": 1, "droplines": 5, "recast": 6, "yaxis": 5, "xaxis": 5,
    "estopts": 7, "predopts": 8, "level": 3, "noweights": 8,
}

GRAPH_OPTS = {
    "title": 2, "subtitle": 3, "note": 4, "caption": 3,
    "xtitle": 6, "ytitle": 6, "xlabel": 4, "ylabel": 4, "xtick": 5, "ytick": 5,
    "xmlabel": 5, "ymlabel": 5, "xmtick": 5, "ymtick": 5,
    "xscale": 3, "yscale": 3, "legend": 3, "name": 4, "saving": 3, "scheme": 3,
    "by": 2, "yline": 5, "xline": 5, "aspectratio": 6, "xsize": 5, "ysize": 5,
    "nodraw": 6, "graphregion": 6, "plotregion": 5, "text": 4, "xcommon": 7, "ycommon": 7,
}

_PLAIN = re.compile(r"^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$")


def _match(name: str, spec: dict[str, int]) -> str | None:
    if name in spec:
        return name
    for full, minimum in spec.items():
        if len(name) >= minimum and full.startswith(name):
            return full
    return None


def split_layer_graph(text: str, layer_spec: dict[str, int] | None = None,
                      extra_graph: dict[str, int] | None = None) -> tuple[dict, dict]:
    """Separa as opções em (camada, gráfico). Devolve dicionários
    {nome completo: argumento (str) ou True}. `no...` é aceito nas opções do
    gráfico como valor False (por exemplo nodraw → nodraw=True)."""
    layer_spec = layer_spec if layer_spec is not None else LAYER_OPTS
    gspec = dict(GRAPH_OPTS)
    if extra_graph:
        gspec.update(extra_graph)
    layer: dict = {}
    graph: dict = {}
    for name, arg in parse_options(text):
        full = _match(name, layer_spec)
        if full is not None:
            layer[full] = arg if arg is not None else True
            continue
        full = _match(name, gspec)
        if full is not None:
            graph[full] = arg if arg is not None else True
            continue
        raise StataError(198, f"option {name} not allowed")
    return layer, graph


# ---------------------------------------------------------------------------
# argumentos
# ---------------------------------------------------------------------------

def split_sub(arg: str) -> tuple[str, dict]:
    """'texto, sub(1) outra' -> ('texto', {'sub': '1', 'outra': True})."""
    if arg is None or arg is True:
        return "", {}
    c = find_top(arg, ",")
    if c == -1:
        return arg.strip(), {}
    subs: dict = {}
    for n, a in parse_options(arg[c + 1:]):
        subs[n] = a if a is not None else True
    return arg[:c].strip(), subs


def text_lines(arg) -> list[str]:
    """Argumento de title(), note()...: várias strings entre aspas são
    linhas; texto sem aspas é uma linha só."""
    if arg is None or arg is True:
        return []
    main, _ = split_sub(arg)
    if not main:
        return []
    if main.startswith('"') or main.startswith('`"'):
        return [w for w in split_words(main)]
    return [main]


def text_opt(arg) -> dict:
    main, subs = split_sub(arg if isinstance(arg, str) else "")
    lines = text_lines(main) if main else []
    out: dict = {"lines": lines}
    for k in ("size", "color", "position", "ring", "angle"):
        for n, v in subs.items():
            if k.startswith(n) and len(n) >= 3:
                out[k] = strip_outer_quotes(v) if isinstance(v, str) else v
    return out


def axis_label(arg) -> dict:
    """xlabel()/ylabel(): valores (numlist), pares valor "rótulo", #k,
    none, minmax e subopções angle(), grid/nogrid, valuelabel, format()."""
    out: dict = {}
    if arg is None or arg is True:
        return out
    main, subs = split_sub(arg)
    for n, v in subs.items():
        nl = n.lower()
        if nl.startswith("ang") or nl == "a":
            try:
                out["angle"] = float(strip_outer_quotes(v)) if v not in ("vertical", "v") else 90.0
            except (TypeError, ValueError):
                out["angle"] = 90.0 if str(v).startswith("v") else 0.0
        elif nl == "grid":
            out["grid"] = True
        elif nl == "nogrid":
            out["grid"] = False
        elif nl.startswith("val"):
            out["valuelabel"] = True
        elif nl.startswith("form"):
            out["format"] = strip_outer_quotes(v or "")
        elif nl.startswith("labs"):
            out["size"] = v
        elif nl == "nolabels":
            out["none"] = True
    words = split_words(main, keep_quotes=True)
    ticks: list = []
    nums: list[str] = []

    def flush() -> None:
        if nums:
            for x in parse_numlist(" ".join(nums)):
                ticks.append([x, None])
            nums.clear()

    for w in words:
        if w.startswith('"') or w.startswith('`"'):
            flush()
            if not ticks:
                raise StataError(198, "invalid label specifier")
            ticks[-1][1] = strip_outer_quotes(w)
        elif w.startswith("##"):
            flush()
            out["minor"] = int(w[2:])
        elif w.startswith("#"):
            flush()
            try:
                out["auto"] = int(w[1:])
            except ValueError:
                raise StataError(198, "invalid label specifier")
        elif w in ("none", "minmax"):
            flush()
            out[w] = True
        else:
            nums.append(w)
    flush()
    if ticks:
        out["ticks"] = ticks
    return out


def axis_scale(arg) -> dict:
    out: dict = {}
    if arg is None or arg is True:
        return out
    for n, a in parse_options(arg):
        if n == "log":
            out["log"] = True
        elif n.startswith("rev"):
            out["reverse"] = True
        elif n.startswith("r") and "range".startswith(n):
            vals = parse_numlist(a or "")
            out["range"] = [min(vals), max(vals)] if vals else None
        elif n in ("off", "noline"):
            out[n] = True
        elif n.startswith("alt"):
            out["alt"] = True
        else:
            raise StataError(198, f"option {n} not allowed")
    return out


def legend_opt(arg) -> dict:
    out: dict = {}
    if arg is None or arg is True:
        return out
    for n, a in parse_options(arg):
        if n == "off":
            out["off"] = True
        elif n == "on":
            out["on"] = True
        elif n.startswith("ord"):
            order: list = []
            for w in split_words(a or "", keep_quotes=True):
                if w.startswith('"') or w.startswith('`"'):
                    if order:
                        order[-1][1] = strip_outer_quotes(w)
                elif w == "-":
                    continue
                else:
                    order.append([int(w), None])
            out["order"] = order
        elif n.startswith("lab"):
            ws = split_words(a or "")
            if ws:
                out.setdefault("labels", {})[str(int(ws[0]))] = " ".join(ws[1:])
        elif n in ("rows", "row", "r"):
            out["rows"] = int(a)
        elif n in ("cols", "col", "c"):
            out["cols"] = int(a)
        elif n.startswith("pos") or n == "position":
            out["position"] = int(a)
        elif n == "ring":
            out["ring"] = int(a)
        elif n.startswith("ti"):
            out["title"] = " ".join(text_lines(a))
        else:
            pass     # VERIFICAR: subopções de aparência ignoradas (size, region...)
    return out


def name_opt(arg) -> tuple[str, bool]:
    main, subs = split_sub(arg or "")
    name = main.strip()
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,31}", name):
        raise StataError(198, f"{name} invalid name")
    return name, any("replace".startswith(k) and len(k) >= 3 for k in subs)


def saving_opt(arg) -> tuple[str, bool]:
    main, subs = split_sub(arg or "")
    return strip_outer_quotes(main), any("replace".startswith(k) and len(k) >= 3 for k in subs)


def lines_opt(arg) -> dict:
    main, subs = split_sub(arg or "")
    out: dict = {"values": parse_numlist(main) if main else []}
    for n, v in subs.items():
        if n.startswith("lc"):
            out["lcolor"] = strip_outer_quotes(v)
        elif n.startswith("lp"):
            out["lpattern"] = strip_outer_quotes(v)
        elif n.startswith("lw"):
            out["lwidth"] = strip_outer_quotes(v)
    return out


def by_opt(arg) -> dict:
    main, subs = split_sub(arg or "")
    out: dict = {"vars": main.split()}
    for n, v in subs.items():
        if n.startswith("tot"):
            out["total"] = True
        elif n in ("cols", "col", "c"):
            out["cols"] = int(v)
        elif n in ("rows", "row", "r"):
            out["rows"] = int(v)
        elif n.startswith("ti"):
            out["title"] = text_lines(v)
        elif n.startswith("note"):
            out["note"] = text_lines(v)
        elif n.startswith("leg"):
            out["legend"] = legend_opt(v)
        elif n in ("compact", "missing", "ixaxes", "iyaxes", "iscale", "noiyaxes", "noixaxes"):
            out[n] = True
    return out


def graph_options(graph: dict) -> dict:
    """Converte as opções do gráfico (texto) para a forma guardada no Graph."""
    o: dict = {}
    for key in ("title", "subtitle", "note", "caption", "xtitle", "ytitle"):
        if key in graph:
            o[key] = text_opt(graph[key])
    for key in ("xlabel", "ylabel"):
        if key in graph:
            o[key] = axis_label(graph[key])
    for key in ("xscale", "yscale"):
        if key in graph:
            o[key] = axis_scale(graph[key])
    if "legend" in graph:
        o["legend"] = legend_opt(graph["legend"])
    for key in ("yline", "xline"):
        if key in graph:
            o[key] = lines_opt(graph[key])
    if "by" in graph:
        o["by"] = by_opt(graph["by"])
    if "scheme" in graph:
        o["scheme"] = strip_outer_quotes(str(graph["scheme"]))
    if "aspectratio" in graph:
        o["aspectratio"] = float(split_sub(graph["aspectratio"])[0])
    for key in ("xsize", "ysize"):
        if key in graph:
            o[key] = float(graph[key])
    if "text" in graph:
        ws = split_words(graph["text"])
        o["text"] = [[float(ws[i]), float(ws[i + 1]), ws[i + 2]] for i in range(0, len(ws) - 2, 3)]
    return o


def is_number(text: str) -> bool:
    return bool(_PLAIN.match(text.strip()))
