"""Simple Graph Editor: side panel of the Graph window.

Edits titles, legend, scheme, size, grid and log scales of the
graph in memory and redraws it at once. As in Stata, the edits change the
graph (and whatever is saved with graph save/export afterward), not the command
that created it. "Revert" goes back to the graph as it was when the
editor was opened.

The editing itself (apply_edits) does not depend on Qt, so that it can be tested.
"""

from __future__ import annotations

import copy

from ..graphics import scheme as SC
from ..graphics.model import Graph
from ..graphics.render import DEFAULT_SIZE
from ..i18n import tr

TEXT_KEYS = ("title", "subtitle", "note", "caption", "xtitle", "ytitle")


def _join(spec) -> str:
    if not spec:
        return ""
    return "\n".join(spec.get("lines", [])) if isinstance(spec, dict) else str(spec)


def current_values(g: Graph) -> dict:
    """Values the editor shows for the graph."""
    o = g.opts
    vals: dict = {}
    for k in TEXT_KEYS:
        if k in o:
            vals[k] = _join(o[k])
        else:
            vals[k] = str(o.get("_" + k, "")) if k in ("xtitle", "ytitle") else ""
    vals["legend"] = not o.get("legend", {}).get("off", False)
    vals["scheme"] = o.get("scheme") or ""
    vals["xsize"] = float(o.get("xsize", DEFAULT_SIZE[0]))
    vals["ysize"] = float(o.get("ysize", DEFAULT_SIZE[1]))
    vals["ygrid"] = o.get("ylabel", {}).get("grid", True)
    vals["xlog"] = bool(o.get("xscale", {}).get("log"))
    vals["ylog"] = bool(o.get("yscale", {}).get("log"))
    return vals


def apply_edits(g: Graph, new: dict) -> None:
    """Writes into the graph only what changed relative to current_values(g)."""
    old = current_values(g)
    o = g.opts
    for k in TEXT_KEYS:
        if k in new and new[k] != old[k]:
            lines = [ln for ln in str(new[k]).split("\n")] if new[k] else [""]
            spec = dict(o.get(k, {})) if isinstance(o.get(k), dict) else {}
            spec["lines"] = lines if new[k] else []
            if not new[k] and k not in ("xtitle", "ytitle"):
                o.pop(k, None)
            else:
                o[k] = spec
    if "legend" in new and new["legend"] != old["legend"]:
        leg = dict(o.get("legend", {}))
        if new["legend"]:
            leg.pop("off", None)
            leg["on"] = True
        else:
            leg["off"] = True
        o["legend"] = leg
    if "scheme" in new and new["scheme"] != old["scheme"]:
        if new["scheme"]:
            o["scheme"] = new["scheme"]
        else:
            o.pop("scheme", None)
    for k in ("xsize", "ysize"):
        if k in new and abs(float(new[k]) - old[k]) > 1e-9:
            o[k] = float(new[k])
    if "ygrid" in new and new["ygrid"] != old["ygrid"]:
        o["ylabel"] = {**o.get("ylabel", {}), "grid": bool(new["ygrid"])}
    for which in ("x", "y"):
        key = which + "log"
        if key in new and new[key] != old[key]:
            sc = dict(o.get(which + "scale", {}))
            if new[key]:
                sc["log"] = True
            else:
                sc.pop("log", None)
            o[which + "scale"] = sc


def scheme_names() -> list[str]:
    return list(SC.SCHEMES)


def make_editor(parent, on_apply):
    """Editor panel (QWidget). on_apply(values) receives the values."""
    from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout,
                                   QHBoxLayout, QLineEdit, QPushButton, QVBoxLayout, QWidget)

    class Editor(QWidget):
        def __init__(self):
            super().__init__(parent)
            lay = QVBoxLayout(self)
            form = QFormLayout()
            self.fields: dict[str, object] = {}
            for k, label in (("title", "Title"), ("subtitle", "Subtitle"), ("note", "Note"),
                             ("caption", "Caption"), ("xtitle", "X axis title"),
                             ("ytitle", "Y axis title")):
                e = QLineEdit()
                self.fields[k] = e
                form.addRow(tr(label), e)
            self.legend = QCheckBox(tr("Show legend"))
            self.ygrid = QCheckBox(tr("Grid lines on the y axis"))
            self.xlog = QCheckBox(tr("Log scale (x)"))
            self.ylog = QCheckBox(tr("Log scale (y)"))
            for w in (self.legend, self.ygrid, self.xlog, self.ylog):
                form.addRow("", w)
            self.scheme = QComboBox()
            self.scheme.addItem("(default)", "")
            for n in scheme_names():
                self.scheme.addItem(n, n)
            form.addRow(tr("Scheme"), self.scheme)
            self.xsize = QDoubleSpinBox()
            self.ysize = QDoubleSpinBox()
            for sp in (self.xsize, self.ysize):
                sp.setRange(1.0, 20.0)
                sp.setSingleStep(0.5)
                sp.setDecimals(1)
            form.addRow(tr("Width (in)"), self.xsize)
            form.addRow(tr("Height (in)"), self.ysize)
            lay.addLayout(form)
            row = QHBoxLayout()
            self.btn_revert = QPushButton(tr("Revert"))
            self.btn_apply = QPushButton(tr("Apply"))
            self.btn_apply.setDefault(True)
            row.addWidget(self.btn_revert)
            row.addStretch(1)
            row.addWidget(self.btn_apply)
            lay.addLayout(row)
            lay.addStretch(1)
            self.btn_apply.clicked.connect(lambda: on_apply(self.values()))
            for e in self.fields.values():
                e.returnPressed.connect(lambda: on_apply(self.values()))

        def load(self, vals: dict) -> None:
            for k, e in self.fields.items():
                e.setText(str(vals.get(k, "")).replace("\n", " / "))
            self.legend.setChecked(bool(vals.get("legend", True)))
            self.ygrid.setChecked(bool(vals.get("ygrid", True)))
            self.xlog.setChecked(bool(vals.get("xlog")))
            self.ylog.setChecked(bool(vals.get("ylog")))
            sch = vals.get("scheme", "")
            if sch and self.scheme.findData(sch) < 0:
                self.scheme.addItem(sch, sch)         # Stata names (s2color...) are aliases
            self.scheme.setCurrentIndex(max(0, self.scheme.findData(sch)))
            self.xsize.setValue(float(vals.get("xsize", DEFAULT_SIZE[0])))
            self.ysize.setValue(float(vals.get("ysize", DEFAULT_SIZE[1])))

        def values(self) -> dict:
            out = {k: e.text().replace(" / ", "\n") for k, e in self.fields.items()}
            out.update(legend=self.legend.isChecked(), ygrid=self.ygrid.isChecked(),
                       xlog=self.xlog.isChecked(), ylog=self.ylog.isChecked(),
                       scheme=self.scheme.currentData() or "",
                       xsize=self.xsize.value(), ysize=self.ysize.value())
            return out

    return Editor()


def snapshot(g: Graph) -> dict:
    return copy.deepcopy(g.opts)
