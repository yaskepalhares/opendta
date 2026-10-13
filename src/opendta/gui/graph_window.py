"""Janela Graph: uma aba por gráfico com nome (name()), como no Stata 14.

O desenho vem de graphics/render.py; a janela só hospeda a Figure do
matplotlib e oferece gravar (.gph), exportar (PNG, PDF, SVG...) e copiar."""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QKeySequence, QShortcut
from PySide6.QtWidgets import QApplication, QFileDialog, QMainWindow, QTabWidget, QToolBar

from ..graphics import render as R
from ..graphics.model import Graph

if TYPE_CHECKING:
    from ..session import Session


class GraphWindow(QMainWindow):
    def __init__(self, session: "Session", parent=None, *, run: Callable[[str], int] | None = None):
        super().__init__(parent)
        self.session = session
        self.run = run
        self.setWindowTitle("Graph")
        self.resize(720, 560)
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self._close_tab)
        self.tabs.currentChanged.connect(self._tab_changed)
        self.setCentralWidget(self.tabs)
        self.canvases: dict[str, FigureCanvasQTAgg] = {}

        tb = QToolBar("Graph")
        tb.setMovable(False)
        for text, slot, key in (("Save...", self.save_graph, "Ctrl+S"), ("Export...", self.export_graph, None),
                                ("Copy", self.copy_graph, "Ctrl+C")):
            act = QAction(text, self)
            act.triggered.connect(slot)
            if key:
                act.setShortcut(QKeySequence(key))
            tb.addAction(act)
        self.addToolBar(tb)
        QShortcut(QKeySequence.StandardKey.Close, self, self.close)

    # -- exibição ------------------------------------------------------------------
    def show_graph(self, g: Graph) -> None:
        fig = R.figure(g)
        canvas = FigureCanvasQTAgg(fig)
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
        self.setWindowTitle(f"Graph - {g.name}")
        canvas.draw_idle()
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
            self.setWindowTitle(f"Graph - {name}")
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

    # -- ações ------------------------------------------------------------------------
    def _cmd(self, text: str) -> None:
        if self.run is not None:
            self.run(text)
        else:
            self.session.run(text)

    def save_graph(self) -> None:
        name = self.current_name()
        if not name:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save graph", f"{name}.gph", "OpenDTA graph (*.gph)")
        if path:
            self._cmd(f'graph save {name} "{path}", replace')

    def export_graph(self) -> None:
        name = self.current_name()
        if not name:
            return
        filters = "PNG (*.png);;PDF (*.pdf);;SVG (*.svg);;EPS (*.eps);;TIFF (*.tif)"
        path, _ = QFileDialog.getSaveFileName(self, "Export graph", f"{name}.png", filters)
        if path:
            self._cmd(f'graph export "{path}", name({name}) replace')

    def copy_graph(self) -> None:
        w = self.tabs.currentWidget()
        if w is None:
            return
        QApplication.clipboard().setPixmap(w.grab())
