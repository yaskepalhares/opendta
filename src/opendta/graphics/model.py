"""Description of a graph, independent of how it is drawn.

A Graph keeps the layers (scatter, line, bars...), the general options
(titles, axes, legend) and, with by(), the panels. Drawing lives in
render.py; saving to .gph (OpenDTA's own JSON format) in store.py.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass
class Layer:
    kind: str                                   # scatter, line, connected, lfit, qfit, area, bar, ...
    x: list = field(default_factory=list)
    y: list = field(default_factory=list)
    extra: dict = field(default_factory=dict)   # y2 (rcap), weights, labels...
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
    data: dict = field(default_factory=dict)    # for bar/box/pie/matrix: their own structure
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
