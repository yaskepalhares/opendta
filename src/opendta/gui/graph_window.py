"""Graph window: one tab per named graph (name()), as in Stata 14.

The drawing comes from graphics/render.py; the window only hosts the matplotlib
Figure and offers saving (.gph), exporting (PNG, PDF, SVG...), copying and
a simple Graph Editor (gui/graph_editor.py) in a side panel."""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QFileDialog,
    QMainWindow,
    QTabWidget,
    QToolBar,
)

from ..graphics import render as R
from ..graphics.model import Graph
from ..i18n import tr
from ..lazy import lazy
from . import graph_editor as GE
from .style import apply_rounded

# imported with the first graph (start-up time; see opendta/lazy.py)
_qtagg = lazy("matplotlib.backends.backend_qtagg")

if TYPE_CHECKING:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

    from ..session import Session


class GraphWindow(QMainWindow):
    def __init__(self, session: "Session", parent=None, *, run: Callable[[str], int] | None = None):
        super().__init__(parent)
        self.session = session
        self.run = run
        self.setWindowTitle(tr("Graph"))
        self.resize(720, 560)
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self._close_tab)
        self.tabs.currentChanged.connect(self._tab_changed)
        self.setCentralWidget(self.tabs)
        apply_rounded(self)
        self.canvases: dict[str, FigureCanvasQTAgg] = {}

        tb = QToolBar("Graph")
        tb.setMovable(False)
        self._tb_actions: list[tuple[QAction, str]] = []
        for text, slot, key in (("Save...", self.save_graph, "Ctrl+S"), ("Export...", self.export_graph, None),
                                ("Copy", self.copy_graph, "Ctrl+C")):
            act = QAction(tr(text), self)
            act.triggered.connect(slot)
            if key:
                act.setShortcut(QKeySequence(key))
            tb.addAction(act)
            self._tb_actions.append((act, text))
        tb.addSeparator()
        self.act_edit = QAction(tr("Edit graph"), self, checkable=True)
        self.act_edit.setShortcut(QKeySequence("Ctrl+E"))
        self.act_edit.toggled.connect(self._toggle_editor)
        tb.addAction(self.act_edit)
        self._tb_actions.append((self.act_edit, "Edit graph"))
        self.addToolBar(tb)

        # Graph Editor (side panel, hidden until requested)
        self.editor = GE.make_editor(self, self._apply_edits)
        self.editor.btn_revert.clicked.connect(self._revert_edits)
        self.editor_dock = QDockWidget(tr("Graph Editor"), self)
        self.editor_dock.setObjectName("GraphEditor")
        self.editor_dock.setWidget(self.editor)
        self.editor_dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable
                                     | QDockWidget.DockWidgetFeature.DockWidgetMovable)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.editor_dock)
        self.editor_dock.hide()
        self.editor_dock.visibilityChanged.connect(
            lambda vis: self.act_edit.setChecked(vis) if not vis and self.act_edit.isChecked() else None)
        self._originals: dict[str, dict] = {}
        QShortcut(QKeySequence.StandardKey.Close, self, self.close)

    # -- display ------------------------------------------------------------------
    def show_graph(self, g: Graph, *, edited: bool = False) -> None:
        if not edited:
            self._originals.pop(g.name, None)     # new graph: earlier edits no longer apply
        fig = R.figure(g)
        canvas = _qtagg.FigureCanvasQTAgg(fig)
        old = self.canvases.get(g.name)
        if old is not None:
            idx = self.tabs.indexOf(old)
            self.tabs.removeTab(idx)
            old.deleteLater()
            self.tabs.insertTab(idx, canvas, g.name)
        else:
            self.tabs.addTab(canvas, g.name)
        self.canvases[g.name] = canvas
        self.tabs.setCurrentWidget(canvas)
        self.tabs.tabBar().setVisible(self.tabs.count() > 1)
        self.setWindowTitle(f"{tr('Graph')} - {g.name}")
        canvas.draw_idle()
        if not edited and self.editor_dock.isVisible():
            self._load_editor()
        self.show()
        self.raise_()

    def close_graph(self, name: str) -> None:
        names = list(self.canvases) if name in ("", "_all") else [name]
        for n in names:
            c = self.canvases.pop(n, None)
            if c is not None:
                self.tabs.removeTab(self.tabs.indexOf(c))
                c.deleteLater()
        if not self.canvases:
            self.hide()

    def _close_tab(self, idx: int) -> None:
        w = self.tabs.widget(idx)
        for n, c in list(self.canvases.items()):
            if c is w:
                self.close_graph(n)

    def _tab_changed(self, idx: int) -> None:
        name = self.current_name()
        if name:
            self.setWindowTitle(f"{tr('Graph')} - {name}")
            if self.editor_dock.isVisible():
                self._load_editor()
            from ..graphics.store import store
            st = store(self.session)
            if name in st.graphs:
                st.current = name

    def current_name(self) -> str:
        w = self.tabs.currentWidget()
        for n, c in self.canvases.items():
            if c is w:
                return n
        return ""

    # -- actions ------------------------------------------------------------------------
    def _cmd(self, text: str) -> None:
        if self.run is not None:
            self.run(text)
        else:
            self.session.run(text)

    def save_graph(self) -> None:
        name = self.current_name()
        if not name:
            return
        path, _ = QFileDialog.getSaveFileName(self, tr("Save graph"), f"{name}.gph",
                                              tr("OpenDTA graph (*.gph)"))
        if path:
            self._cmd(f'graph save {name} "{path}", replace')

    def export_graph(self) -> None:
        name = self.current_name()
        if not name:
            return
        filters = "PNG (*.png);;PDF (*.pdf);;SVG (*.svg);;EPS (*.eps);;TIFF (*.tif)"
        path, _ = QFileDialog.getSaveFileName(self, tr("Export graph"), f"{name}.png", filters)
        if path:
            self._cmd(f'graph export "{path}", name({name}) replace')

    def copy_graph(self) -> None:
        w = self.tabs.currentWidget()
        if w is None:
            return
        QApplication.clipboard().setPixmap(w.grab())

    # -- Graph Editor ------------------------------------------------------------------
    def _graph(self) -> Graph | None:
        from ..graphics.store import store
        return store(self.session).graphs.get(self.current_name())

    def _toggle_editor(self, on: bool) -> None:
        if on:
            self._load_editor()
            self.editor_dock.show()
        else:
            self.editor_dock.hide()

    def _load_editor(self) -> None:
        g = self._graph()
        if g is None:
            return
        self._originals.setdefault(g.name, GE.snapshot(g))
        self.editor.load(GE.current_values(g))

    def _apply_edits(self, values: dict) -> None:
        g = self._graph()
        if g is None:
            return
        self._originals.setdefault(g.name, GE.snapshot(g))
        GE.apply_edits(g, values)
        self.show_graph(g, edited=True)

    def _revert_edits(self) -> None:
        g = self._graph()
        if g is None or g.name not in self._originals:
            return
        g.opts = self._originals.pop(g.name)
        self.show_graph(g, edited=True)
        self._load_editor()

    def retranslate(self) -> None:
        for act, key in self._tb_actions:
            act.setText(tr(key))
        self.editor_dock.setWindowTitle(tr("Graph Editor"))
        name = self.current_name()
        self.setWindowTitle(f"{tr('Graph')} - {name}" if name else tr("Graph"))
