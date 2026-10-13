"""Gráficos em memória (graph dir/drop/rename/display) e arquivos .gph.

O .gph do OpenDTA é JSON próprio (não é o formato da StataCorp): guarda a
descrição do gráfico (camadas e opções), não a imagem. `graph use` de um
.gph gravado pelo Stata não é suportado."""

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
            # VERIFICAR: mensagem e código quando name() já existe sem replace
            raise StataError(110, f"graph {g.name} already exists")
        # substituir mantém a posição na lista (VERIFICAR a ordem de graph dir)
        self.graphs[g.name] = g
        self.current = g.name

    def get(self, name: str | None = None) -> Graph:
        n = name or self.current
        if not n or n not in self.graphs:
            if not self.graphs:
                raise StataError(111, "no graphs in memory")   # VERIFICAR
            raise StataError(111, f"graph {n} not found")      # VERIFICAR
        return self.graphs[n]

    def drop(self, name: str) -> None:
        if name not in self.graphs:
            raise StataError(111, f"graph {name} not found")   # VERIFICAR
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


def save_gph(g: Graph, path: Path, *, replace: bool) -> None:
    if path.exists() and not replace:
        raise StataError(602, f"file {path} already exists")
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
        # .gph do Stata (binário próprio) ou arquivo corrompido
        raise StataError(610, f"file {path} not an OpenDTA graph file")
    return Graph.from_json(payload["graph"])
