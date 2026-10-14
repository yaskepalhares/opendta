"""Descrição de um gráfico, independente de como é desenhado.

Um Graph guarda as camadas (scatter, line, barras...), as opções gerais
(títulos, eixos, legenda) e, com by(), os painéis. O desenho fica em
render.py; a gravação em .gph (formato próprio do OpenDTA, JSON) em store.py.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class Layer:
    kind: str                                   # scatter, line, connected, lfit, qfit, area, bar, ...
    x: list = field(default_factory=list)
    y: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)   # y2 (rcap), pesos, rótulos...
    style: dict = field(default_factory=dict)   # mcolor, msymbol, lcolor, lpattern, lwidth...
    label: str = ""


@dataclass
class Panel:
    title: str = ""
    layers: list[Layer] = field(default_factory=list)


@dataclass
class Graph:
    family: str = "twoway"                      # twoway, bar, hbar, box, hbox, pie, matrix, combine
    panels: list[Panel] = field(default_factory=list)
    opts: dict = field(default_factory=dict)    # title, subtitle, note, caption, xtitle, ytitle, ...
    name: str = "Graph"
    data: dict = field(default_factory=dict)    # para bar/box/pie/matrix: estrutura própria
    children: list = field(default_factory=list)   # graph combine

    def to_json(self) -> dict[str, Any]:
        d = asdict(self)
        d["children"] = [c.to_json() if isinstance(c, Graph) else c for c in self.children]
        return d

    @staticmethod
    def from_json(d: dict) -> "Graph":
        g = Graph(family=d.get("family", "twoway"), opts=d.get("opts", {}), name=d.get("name", "Graph"),
                  data=d.get("data", {}))
        for p in d.get("panels", []):
            g.panels.append(Panel(p.get("title", ""), [Layer(**ly) for ly in p.get("layers", [])]))
        g.children = [Graph.from_json(c) for c in d.get("children", [])]
        return g
