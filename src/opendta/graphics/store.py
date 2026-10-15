"""Graphs in memory (graph dir/drop/rename/display) and .gph files.

OpenDTA's .gph is its own JSON (not StataCorp's format): it keeps the
description of the graph (layers and options), not the image. `graph use` of a
.gph saved by Stata is not supported."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from .model import Graph

if TYPE_CHECKING:
    from ..session import Session

MAGIC = "OpenDTA graph"
VERSION = 1


class GraphStore:
    def __init__(self) -> None:
        self.graphs: dict[str, Graph] = {}
        self.current: str = ""

    def put(self, g: Graph, *, replace: bool) -> None:
        if g.name in self.graphs and not replace and g.name != "Graph":
            # VERIFY: message and code when name() already exists without replace
            raise StataError(110, f"graph {g.name} already exists")
        # replacing keeps the position in the list (VERIFY the order of graph dir)
        self.graphs[g.name] = g
        self.current = g.name

    def get(self, name: str | None = None) -> Graph:
        n = name or self.current
        if not n or n not in self.graphs:
            if not self.graphs:
                raise StataError(111, "no graphs in memory")   # VERIFY
            raise StataError(111, f"graph {n} not found")      # VERIFY
        return self.graphs[n]

    def drop(self, name: str) -> None:
        if name not in self.graphs:
            raise StataError(111, f"graph {name} not found")   # VERIFY
        del self.graphs[name]
        if self.current == name:
            self.current = next(reversed(self.graphs), "") if self.graphs else ""

    def clear(self) -> None:
        self.graphs.clear()
        self.current = ""


def store(s: "Session") -> GraphStore:
    st = getattr(s, "_graphs", None)
    if st is None:
        st = GraphStore()
        s._graphs = st
    return st


def gph_path(text: str) -> Path:
    from ..lang.words import strip_outer_quotes
    t = strip_outer_quotes(text.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    p = Path(t).expanduser()
    if p.suffix == "":
        p = p.with_suffix(".gph")
    return p


def save_gph(g: Graph, path: Path, *, replace: bool, out=None) -> None:
    if path.exists() and not replace:
        raise StataError(602, f"file {path} already exists")
    if replace and not path.exists() and out is not None:
        # like save, replace: a note when the file does not exist yet (compat 0701)
        out.write(f"(note: file {path} not found)\n", "text")
    payload = {"format": MAGIC, "version": VERSION, "graph": g.to_json()}
    try:
        path.write_text(json.dumps(payload), encoding="utf-8")
    except OSError:
        raise StataError(603, f"file {path} could not be opened")


def load_gph(path: Path) -> Graph:
    if not path.exists():
        raise StataError(601, f"file {path} not found")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, UnicodeDecodeError):
        payload = None
    if not isinstance(payload, dict) or payload.get("format") != MAGIC:
        # Stata .gph (its own binary format) or a corrupt file
        raise StataError(610, f"file {path} not an OpenDTA graph file")
    return Graph.from_json(payload["graph"])
