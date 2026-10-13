"""Comandos de gráfico: graph (e subcomandos), twoway, scatter, line,
connected, tsline, histogram e kdensity ([G-2], [R] histogram, [R] kdensity).

Os gráficos ficam guardados na sessão (graph dir/drop/rename) e aparecem
na janela Graph da interface; sem interface (modo batch), são só
guardados e podem ser gravados com graph export/graph save. Mensagens
seguem as do Stata 14; o que não é documentado tem a marca VERIFICAR e
casos em compat/do/07xx."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..graphics import families as FAM
from ..graphics import render as R
from ..graphics import twoway as TW
from ..graphics.model import Graph
from ..graphics.opts import name_opt, saving_opt
from ..graphics.store import gph_path, load_gph, save_gph, store
from ..lang.syntax import find_top, match_options
from ..lang.words import split_words, strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


# ---------------------------------------------------------------------------
# parte comum: name(), saving(), nodraw, esquema e exibição
# ---------------------------------------------------------------------------

def finish(s: "Session", g: Graph, graw: dict) -> None:
    name, replace = ("Graph", True)
    if "name" in graw and isinstance(graw["name"], str):
        name, replace = name_opt(graw["name"])
    g.name = name
    if "scheme" not in g.opts:
        g.opts["scheme"] = s.settings.get("scheme", "s2color")
    store(s).put(g, replace=replace)
    if "saving" in graw and isinstance(graw["saving"], str):
        fname, rep = saving_opt(graw["saving"])
        path = gph_path(fname)
        save_gph(g, path, replace=rep, out=s.output)
        s.output.write(f"(file {_shown(fname, path)} saved)\n", "text")
    if not graw.get("nodraw"):
        show(s, g)


def show(s: "Session", g: Graph) -> None:
    if s.settings.get("graphics", "on") != "on":
        return
    hook = s.ui_hooks.get("graph")
    if hook is not None:
        hook(g)


def _shown(text: str, path: Path) -> str:
    """Nome do arquivo como o usuário escreveu, com a extensão acrescentada."""
    t = strip_outer_quotes(text.strip())
    return t if Path(t).suffix else t + path.suffix


# ---------------------------------------------------------------------------
# twoway e atalhos
# ---------------------------------------------------------------------------

def _twoway(s: "Session", text: str, default_type: str | None = None) -> None:
    g, graw = TW.build(s, text, default_type=default_type)
    finish(s, g, graw)


@command("twoway", "tw")
def cmd_twoway(s: "Session", args: str) -> None:
    _twoway(s, args)


@command("scatter", "sc")
def cmd_scatter(s: "Session", args: str) -> None:
    _twoway(s, args, "scatter")


@command("line")
def cmd_line(s: "Session", args: str) -> None:
    _twoway(s, args, "line")


@command("connected")
def cmd_connected(s: "Session", args: str) -> None:
    _twoway(s, args, "connected")


@command("tsline")
def cmd_tsline(s: "Session", args: str) -> None:
    _twoway(s, args, "tsline")


@command("histogram", "hist")
def cmd_histogram(s: "Session", args: str) -> None:
    g, graw, msg = FAM.histogram(s, args)
    # com by() o Stata não mostra a mensagem das barras (compat 0701)
    if "by" not in graw:
        s.output.write(msg + "\n", "text")
    finish(s, g, graw)


@command("kdensity", "kdens")
def cmd_kdensity(s: "Session", args: str) -> None:
    g, graw = FAM.kdensity(s, args)
    if g is not None:
        finish(s, g, graw)


# ---------------------------------------------------------------------------
# diagnósticos e pós-estimação
# ---------------------------------------------------------------------------

def _post(fn, s: "Session", args: str) -> None:
    g, graw = fn(s, args)
    if g is not None:
        finish(s, g, graw)


@command("lowess", "lowess")
def cmd_lowess(s: "Session", args: str) -> None:
    from ..graphics.postplots import lowess
    _post(lowess, s, args)


@command("qnorm", "qnorm")
def cmd_qnorm(s: "Session", args: str) -> None:
    from ..graphics.postplots import qnorm
    _post(qnorm, s, args)


@command("pnorm", "pnorm")
def cmd_pnorm(s: "Session", args: str) -> None:
    from ..graphics.postplots import pnorm
    _post(pnorm, s, args)


@command("rvfplot", "rvfplot")
def cmd_rvfplot(s: "Session", args: str) -> None:
    from ..graphics.postplots import rvfplot
    _post(rvfplot, s, args)


@command("rvpplot", "rvpplot")
def cmd_rvpplot(s: "Session", args: str) -> None:
    from ..graphics.postplots import rvpplot
    _post(rvpplot, s, args)


@command("marginsplot", "marginsplot")
def cmd_marginsplot(s: "Session", args: str) -> None:
    from ..graphics.postplots import marginsplot
    _post(marginsplot, s, args)


# ---------------------------------------------------------------------------
# graph
# ---------------------------------------------------------------------------

_SUB = [("twoway", 2), ("bar", 3), ("hbar", 4), ("box", 3), ("hbox", 4), ("pie", 3), ("matrix", 3),
        ("combine", 4), ("export", 3), ("save", 4), ("use", 3), ("dir", 3), ("drop", 4),
        ("display", 2), ("rename", 3), ("copy", 4), ("describe", 1), ("query", 1), ("close", 5),
        ("dot", 3)]


@command("graph", "gr")
def cmd_graph(s: "Session", args: str) -> None:
    import re
    m = re.match(r"\s*([A-Za-z_]\w*)(.*)$", args, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    head, rest = m.group(1), m.group(2).strip()
    sub = None
    for name, m in _SUB:
        if len(head) >= m and name.startswith(head):
            sub = name
            break
    if sub is None:
        # graph scatter y x == twoway scatter y x
        if TW.plottype(head) is not None:
            _twoway(s, args)
            return
        raise StataError(198, f"{head} is not a graph command")   # VERIFICAR
    if sub == "twoway":
        _twoway(s, rest)
    elif sub in ("bar", "hbar"):
        g, graw = FAM.bar(s, rest, sub)
        finish(s, g, graw)
    elif sub == "dot":
        g, graw = FAM.bar(s, rest, "dot")
        finish(s, g, graw)
    elif sub in ("box", "hbox"):
        g, graw = FAM.box(s, rest, sub)
        finish(s, g, graw)
    elif sub == "pie":
        g, graw = FAM.pie(s, rest)
        finish(s, g, graw)
    elif sub == "matrix":
        g, graw = FAM.matrix(s, rest)
        finish(s, g, graw)
    elif sub == "combine":
        g, graw = FAM.combine(s, rest)
        finish(s, g, graw)
    elif sub == "export":
        _export(s, rest)
    elif sub == "save":
        _save(s, rest)
    elif sub == "use":
        _use(s, rest)
    elif sub == "dir":
        _dir(s, rest)
    elif sub == "drop":
        _drop(s, rest)
    elif sub == "display":
        _display(s, rest)
    elif sub == "rename":
        _rename(s, rest)
    elif sub == "copy":
        _copy(s, rest)
    elif sub == "describe":
        _describe(s, rest)
    elif sub == "query":
        _query(s, rest)
    elif sub == "close":
        hook = s.ui_hooks.get("graph_close")
        if hook is not None:
            hook(rest.strip())


def _split(text: str) -> tuple[str, str]:
    c = find_top(text, ",")
    return (text, "") if c == -1 else (text[:c], text[c + 1:])


def _export(s: "Session", text: str) -> None:
    body, opt = _split(text)
    o = match_options(opt, {"as": 2, "replace": 3, "name": 4, "width": 5, "height": 6})
    fname = strip_outer_quotes(body.strip())
    if not fname:
        raise StataError(198, "invalid syntax")
    path = Path(fname).expanduser()
    if isinstance(o.get("as"), str):
        ext = "." + o["as"].strip().lower()
        if ext not in R.FORMATS:
            raise StataError(198, f"as({o['as']}) invalid")   # VERIFICAR
        if path.suffix == "":
            path = path.with_suffix(ext)
            fname += ext
    else:
        ext = path.suffix.lower()
        if ext == "":
            raise StataError(198, "you must specify a file extension or the as() option")   # VERIFICAR
    if ext not in R.FORMATS:
        raise StataError(198, f"{ext[1:]} is not a supported graph export format")   # VERIFICAR
    g = store(s).get(o["name"].strip() if isinstance(o.get("name"), str) else None)
    if path.exists() and not o.get("replace"):
        raise StataError(602, f"file {fname} already exists")
    fmt, label = R.FORMATS[ext]
    width = int(o["width"]) if isinstance(o.get("width"), str) else None
    height = int(o["height"]) if isinstance(o.get("height"), str) else None
    try:
        R.export(g, path, fmt, width=width, height=height)
    except OSError:
        raise StataError(603, f"file {fname} could not be opened")
    # no Mac, o Stata 14 mostra o caminho completo só no PDF (compat 0701)
    shown = str(path.resolve()) if label == "PDF" else fname
    s.output.write(f"(file {shown} written in {label} format)\n", "text")


def _save(s: "Session", text: str) -> None:
    body, opt = _split(text)
    o = match_options(opt, {"replace": 3, "asis": 4})
    words = split_words(body)
    st = store(s)
    if len(words) == 1:
        g = st.get(None)
        fname = words[0]
    elif len(words) == 2:
        g = st.get(words[0])
        fname = words[1]
    else:
        raise StataError(198, "invalid syntax")
    path = gph_path(fname)
    save_gph(g, path, replace=bool(o.get("replace")), out=s.output)
    s.output.write(f"(file {_shown(fname, path)} saved)\n", "text")


def _use(s: "Session", text: str) -> None:
    body, opt = _split(text)
    o = match_options(opt, {"name": 4, "nodraw": 6, "scheme": 3})
    words = split_words(body)
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    path = gph_path(words[0])
    g = load_gph(path)
    # VERIFICAR: o nome do gráfico carregado é o nome do arquivo sem extensão
    name, replace = (path.stem, True)
    if isinstance(o.get("name"), str):
        name, replace = name_opt(o["name"])
    if isinstance(o.get("scheme"), str):
        g.opts["scheme"] = strip_outer_quotes(o["scheme"])
    g.name = name
    store(s).put(g, replace=replace)
    if not o.get("nodraw"):
        show(s, g)


def _dir(s: "Session", text: str) -> None:
    body, opt = _split(text)
    o = match_options(opt, {"memory": 3, "gph": 3, "detail": 3})
    names = list(store(s).graphs)
    pattern = body.strip()
    if pattern:
        import fnmatch
        names = [n for n in names if fnmatch.fnmatchcase(n, pattern)]
    files: list[str] = []
    if not o.get("memory"):
        import glob
        files = sorted(Path(f).name for f in glob.glob("*.gph"))
    items = names + files
    if not items:
        return
    # em ordem alfabética, colunas da largura do maior nome + 2, recuo de 4
    # (compat 0701: "    d       e       f       odta_a  p")
    items = sorted(items)
    L = max(len(x) for x in items) + 2
    per = max(1, (int(s.settings.get("linesize", 80)) - 4) // L)
    for k in range(0, len(items), per):
        line = "".join(f"{x:<{L}}" for x in items[k:k + per]).rstrip()
        s.output.write("    " + line + "\n", "result")
    s.r = {"list": " ".join(names)}


def _drop(s: "Session", text: str) -> None:
    st = store(s)
    words = split_words(text)
    if not words:
        raise StataError(198, "invalid syntax")
    if words == ["_all"]:
        st.clear()
        hook = s.ui_hooks.get("graph_close")
        if hook is not None:
            hook("_all")
        return
    import fnmatch
    for w in words:
        hits = [n for n in st.graphs if fnmatch.fnmatchcase(n, w)] if any(c in w for c in "*?") else [w]
        for n in hits:
            if n not in st.graphs:
                raise StataError(111, f"graph {n} not found\nno graphs dropped")
            st.drop(n)
            hook = s.ui_hooks.get("graph_close")
            if hook is not None:
                hook(n)


def _display(s: "Session", text: str) -> None:
    body, opt = _split(text)
    o = match_options(opt, {"scheme": 3, "xsize": 5, "ysize": 5, "margins": 3})
    st = store(s)
    g = st.get(body.strip() or None)
    if isinstance(o.get("scheme"), str):
        g.opts["scheme"] = strip_outer_quotes(o["scheme"])
    for k in ("xsize", "ysize"):
        if isinstance(o.get(k), str):
            g.opts[k] = float(o[k])
    st.current = g.name
    show(s, g)


def _rename(s: "Session", text: str) -> None:
    body, opt = _split(text)
    o = match_options(opt, {"replace": 3})
    words = split_words(body)
    st = store(s)
    if len(words) == 1:
        old, new = st.current, words[0]
    elif len(words) == 2:
        old, new = words
    else:
        raise StataError(198, "invalid syntax")
    g = st.get(old)
    name_opt(new)
    if new in st.graphs and not o.get("replace"):
        raise StataError(110, f"graph {new} already exists")   # VERIFICAR
    st.drop(old)
    g.name = new
    st.put(g, replace=True)


def _copy(s: "Session", text: str) -> None:
    body, opt = _split(text)
    o = match_options(opt, {"replace": 3})
    words = split_words(body)
    st = store(s)
    if len(words) == 1:
        src, new = st.current, words[0]
    elif len(words) == 2:
        src, new = words
    else:
        raise StataError(198, "invalid syntax")
    g = Graph.from_json(st.get(src).to_json())
    g.name = new
    st.put(g, replace=bool(o.get("replace")))


def _describe(s: "Session", text: str) -> None:
    st = store(s)
    t = text.strip()
    if t.endswith(".gph") or (t and t not in st.graphs and gph_path(t).exists()):
        g = load_gph(gph_path(t))
        where = str(gph_path(t))
    else:
        g = st.get(t or None)
        where = "memory"
    out = s.output
    # VERIFICAR: layout de graph describe
    out.write(f"\n{g.name} stored in {where}\n\n", "text")
    out.write(f"{'family:':>12}  {g.family}\n", "text")
    out.write(f"{'scheme:':>12}  {g.opts.get('scheme', '')}\n", "text")
    out.write(f"{'panels:':>12}  {len(g.panels)}\n", "text")


def _query(s: "Session", text: str) -> None:
    t = text.strip()
    if t.startswith("sch"):
        from ..graphics.scheme import ALIASES, SCHEMES
        out = s.output
        out.write("Available schemes are\n\n", "text")
        for name in list(SCHEMES) + sorted(ALIASES):
            note = "" if name in SCHEMES else f"   (= {ALIASES[name]})"
            out.write(f"    {name}{note}\n", "text")
        return
    s.output.write("Graph families: twoway bar hbar box hbox pie matrix combine\n", "text")
