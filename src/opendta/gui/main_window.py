"""Janela principal no layout do Stata 14.

    ┌──────────────────────────────────────────────────────────┐
    │ Menus: File Edit Data Graphics Statistics User Window Help│
    │ Barra de ferramentas                                      │
    ├──────────┬───────────────────────────────┬───────────────┤
    │ Review   │ Results                       │ Variables     │
    │          │                               ├───────────────┤
    │          ├───────────────────────────────┤ Properties    │
    │          │ Command                       │               │
    ├──────────┴───────────────────────────────┴───────────────┤
    │ diretório de trabalho                       CAP  NUM  OVR │
    └──────────────────────────────────────────────────────────┘

Ícones e identidade visual são próprios do OpenDTA. Itens de menu ainda não
implementados aparecem desabilitados com a fase prevista na dica.
"""

from __future__ import annotations

import os
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import (QAction, QColor, QFont, QFontDatabase, QKeySequence,
                           QTextCharFormat, QTextCursor)
from PySide6.QtWidgets import (QDockWidget, QFileDialog, QHeaderView, QLabel,
                               QLineEdit, QMainWindow, QMenu, QMessageBox,
                               QPlainTextEdit, QSplitter, QStyle, QToolBar,
                               QTreeWidget, QTreeWidgetItem, QVBoxLayout,
                               QWidget)

from .. import __version__
from ..core.errors import ExitRequest
from ..session import Session
from .theme import MONOSPACE_FAMILIES, MONOSPACE_SIZE, STANDARD, ResultsScheme


def monospace_font() -> QFont:
    available = set(QFontDatabase.families())
    for fam in MONOSPACE_FAMILIES:
        if fam in available:
            f = QFont(fam, MONOSPACE_SIZE)
            break
    else:
        f = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        f.setPointSize(MONOSPACE_SIZE)
    f.setStyleHint(QFont.StyleHint.Monospace)
    return f


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

class ResultsView(QPlainTextEdit):
    def __init__(self, scheme: ResultsScheme = STANDARD):
        super().__init__()
        self.setReadOnly(True)
        self.setUndoRedoEnabled(False)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setFont(monospace_font())
        self.setObjectName("Results")
        self.setStyleSheet(f"QPlainTextEdit#Results {{ background: {scheme.background}; }}")
        self._formats: dict[str, QTextCharFormat] = {}
        for style, color in scheme.colors.items():
            fmt = QTextCharFormat()
            fmt.setForeground(QColor(color))
            if style in scheme.bold:
                fmt.setFontWeight(QFont.Weight.Bold)
            self._formats[style] = fmt

    def append_styled(self, text: str, style: str) -> None:
        cursor = self.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        cursor.insertText(text, self._formats.get(style, self._formats["text"]))
        self.setTextCursor(cursor)
        self.ensureCursorVisible()


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------

class CommandLine(QLineEdit):
    submitted = Signal(str)

    def __init__(self):
        super().__init__()
        self.setFont(monospace_font())
        self.setObjectName("Command")
        self.history: list[str] = []
        self._pos = 0
        self.returnPressed.connect(self._submit)

    def _submit(self) -> None:
        text = self.text()
        if not text.strip():
            return
        self.history.append(text)
        self._pos = len(self.history)
        self.clear()
        self.submitted.emit(text)

    def keyPressEvent(self, event) -> None:  # noqa: N802
        key = event.key()
        if key in (Qt.Key.Key_PageUp, Qt.Key.Key_Up) and self.history:
            self._pos = max(0, self._pos - 1)
            self.setText(self.history[self._pos])
            return
        if key in (Qt.Key.Key_PageDown, Qt.Key.Key_Down) and self.history:
            self._pos = min(len(self.history), self._pos + 1)
            self.setText(self.history[self._pos] if self._pos < len(self.history) else "")
            return
        if key == Qt.Key.Key_Escape:
            self.clear()
            return
        super().keyPressEvent(event)


def _titled(title: str, widget: QWidget) -> QWidget:
    """Painel com barra de título fina, como os painéis encaixados do Stata."""
    box = QWidget()
    lay = QVBoxLayout(box)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(0)
    label = QLabel(f"  {title}")
    label.setObjectName("PaneTitle")
    lay.addWidget(label)
    lay.addWidget(widget)
    return box


# ---------------------------------------------------------------------------
# Janela principal
# ---------------------------------------------------------------------------

_PENDING = {
    1: "Disponível na fase 1 (dados em memória)",
    2: "Disponível na fase 2 (programação)",
    3: "Disponível na fase 3 (descritivas e manipulação)",
    5: "Disponível na fase 5 (estimação)",
    6: "Disponível na fase 6 (amostras complexas)",
    7: "Disponível na fase 7 (gráficos)",
    8: "Disponível na fase 8 (ferramentas da interface)",
}


class MainWindow(QMainWindow):
    def __init__(self, session: Session | None = None):
        super().__init__()
        self.session = session or Session()
        self.resize(1280, 800)
        self._build_central()
        self._build_docks()
        self._build_menus()
        self._build_toolbar()
        self._build_statusbar()
        self._apply_styles()

        self.session.output.add_listener(self.results.append_styled)
        self.session.add_state_listener(self.refresh_state)
        self.refresh_state()
        self.command.setFocus()
        self.setAcceptDrops(True)

    # -- montagem ------------------------------------------------------------
    def _build_central(self) -> None:
        self.results = ResultsView()
        self.command = CommandLine()
        self.command.submitted.connect(self.run_command)

        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(_titled("Results", self.results))
        split.addWidget(_titled("Command", self.command))
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 0)
        split.setCollapsible(0, False)
        split.setCollapsible(1, False)
        split.setSizes([700, 60])
        self.setCentralWidget(split)

    def _dock(self, title: str, widget: QWidget, area: Qt.DockWidgetArea) -> QDockWidget:
        d = QDockWidget(title, self)
        d.setObjectName(title)
        d.setWidget(widget)
        d.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable
                      | QDockWidget.DockWidgetFeature.DockWidgetMovable
                      | QDockWidget.DockWidgetFeature.DockWidgetFloatable)
        self.addDockWidget(area, d)
        return d

    def _build_docks(self) -> None:
        # Review
        self.review = QTreeWidget()
        self.review.setHeaderLabels(["Command", "_rc"])
        self.review.setRootIsDecorated(False)
        self.review.setFont(monospace_font())
        self.review.header().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.review.header().setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        self.review.header().setStretchLastSection(False)
        self.review.itemClicked.connect(lambda it, _c: self.command.setText(it.text(0)))
        self.review.itemDoubleClicked.connect(lambda it, _c: self.run_command(it.text(0)))
        review_box = QWidget()
        lay = QVBoxLayout(review_box)
        lay.setContentsMargins(0, 0, 0, 0)
        self.review_filter = QLineEdit()
        self.review_filter.setPlaceholderText("Filter commands here")
        self.review_filter.textChanged.connect(self._filter_review)
        lay.addWidget(self.review_filter)
        lay.addWidget(self.review)
        self.dock_review = self._dock("Review", review_box, Qt.DockWidgetArea.LeftDockWidgetArea)

        # Variables
        self.variables = QTreeWidget()
        self.variables.setHeaderLabels(["Name", "Label"])
        self.variables.setRootIsDecorated(False)
        var_box = QWidget()
        lay = QVBoxLayout(var_box)
        lay.setContentsMargins(0, 0, 0, 0)
        self.var_filter = QLineEdit()
        self.var_filter.setPlaceholderText("Filter variables here")
        lay.addWidget(self.var_filter)
        lay.addWidget(self.variables)
        self.dock_variables = self._dock("Variables", var_box, Qt.DockWidgetArea.RightDockWidgetArea)

        # Properties
        self.properties = QTreeWidget()
        self.properties.setHeaderLabels(["Property", "Value"])
        self.properties.header().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self._prop_vars = QTreeWidgetItem(["Variables"])
        self._prop_data = QTreeWidgetItem(["Data"])
        for parent, keys in ((self._prop_vars, ["Name", "Label", "Type", "Format", "Value label", "Notes"]),
                             (self._prop_data, ["Filename", "Label", "Notes", "Variables",
                                                "Observations", "Size", "Memory", "Sorted by"])):
            for k in keys:
                parent.addChild(QTreeWidgetItem([k, ""]))
            self.properties.addTopLevelItem(parent)
            parent.setExpanded(True)
        self.dock_properties = self._dock("Properties", self.properties, Qt.DockWidgetArea.RightDockWidgetArea)

        self.splitDockWidget(self.dock_variables, self.dock_properties, Qt.Orientation.Vertical)
        self.resizeDocks([self.dock_review], [260], Qt.Orientation.Horizontal)
        self.resizeDocks([self.dock_variables], [300], Qt.Orientation.Horizontal)
        self.resizeDocks([self.dock_variables, self.dock_properties], [420, 330], Qt.Orientation.Vertical)

    def _action(self, menu, text: str, slot=None, *, shortcut: str | None = None,
                phase: int | None = None) -> QAction:
        act = QAction(text, self)
        if shortcut:
            act.setShortcut(QKeySequence(shortcut))
        if slot is not None:
            act.triggered.connect(slot)
        if phase is not None:
            act.setEnabled(False)
            act.setToolTip(_PENDING.get(phase, ""))
            act.setStatusTip(_PENDING.get(phase, ""))
        menu.addAction(act)
        return act

    def _build_menus(self) -> None:
        mb = self.menuBar()

        m = mb.addMenu("&File")
        self._action(m, "Open...", shortcut="Ctrl+O", phase=1)
        self._action(m, "Save", shortcut="Ctrl+S", phase=1)
        self._action(m, "Save as...", phase=1)
        m.addSeparator()
        self._action(m, "Do...", self.choose_do_file)
        self._action(m, "Change working directory...", self.choose_directory)
        m.addSeparator()
        log = m.addMenu("Log")
        self._action(log, "Begin...", phase=2)
        self._action(log, "Close", phase=2)
        self._action(m, "Print...", phase=8)
        m.addSeparator()
        self._action(m, "Exit", self.close)

        m = mb.addMenu("&Edit")
        self._action(m, "Copy", self.results.copy, shortcut="Ctrl+C")
        self._action(m, "Copy table", phase=3)
        self._action(m, "Find...", phase=8)
        m.addSeparator()
        self._action(m, "Clear Results", self.results.clear)
        self._action(m, "Preferences", phase=8)

        m = mb.addMenu("&Data")
        self._action(m, "Describe data", phase=1)
        self._action(m, "Data Editor", phase=1)
        self._action(m, "Create or change data", phase=1)
        self._action(m, "Variables Manager", phase=1)
        self._action(m, "Data utilities", phase=3)
        self._action(m, "Sort", phase=1)
        self._action(m, "Combine datasets", phase=3)

        m = mb.addMenu("&Graphics")
        for item in ("Twoway graph (scatter, line, etc.)", "Bar chart", "Dot chart", "Pie chart",
                     "Histogram", "Box plot", "Scatterplot matrix", "Graph Editor"):
            self._action(m, item, phase=7)

        m = mb.addMenu("&Statistics")
        for item, phase in (("Summaries, tables, and tests", 3), ("Linear models and related", 5),
                            ("Binary outcomes", 5), ("Ordinal outcomes", 5),
                            ("Count outcomes", 5), ("Generalized linear models", 5),
                            ("Longitudinal/panel data", 5), ("Survival analysis", 8),
                            ("Survey data analysis", 6), ("Postestimation", 5)):
            self._action(m, item, phase=phase)

        mb.addMenu("&User")

        m = mb.addMenu("&Window")
        for dock in (self.dock_review, self.dock_variables, self.dock_properties):
            m.addAction(dock.toggleViewAction())
        m.addSeparator()
        self._action(m, "Command", lambda: self.command.setFocus())
        self._action(m, "Results", lambda: self.results.setFocus())
        self._action(m, "Graph", phase=7)
        self._action(m, "Viewer", phase=2)
        self._action(m, "Data Editor", phase=1)
        self._action(m, "Do-file Editor", phase=8)

        m = mb.addMenu("&Help")
        self._action(m, "Search...", phase=2)
        self._action(m, "About OpenDTA", self.show_about)

        for menu in mb.findChildren(QMenu):
            menu.setToolTipsVisible(True)

    def _build_toolbar(self) -> None:
        tb = QToolBar("Toolbar")
        tb.setObjectName("Toolbar")
        tb.setMovable(False)
        st = self.style()
        P = QStyle.StandardPixmap
        items = [
            ("Open", P.SP_DialogOpenButton, None, 1),
            ("Save", P.SP_DialogSaveButton, None, 1),
            ("Print", P.SP_FileIcon, None, 8),
            ("Log", P.SP_FileDialogDetailedView, None, 2),
            ("Viewer", P.SP_FileDialogInfoView, None, 2),
            ("Graph", P.SP_DesktopIcon, None, 7),
            ("Do-file Editor", P.SP_FileDialogContentsView, None, 8),
            ("Data Editor (Edit)", P.SP_FileDialogListView, None, 1),
            ("Data Browser (Browse)", P.SP_FileDialogStart, None, 1),
            ("Variables Manager", P.SP_DirIcon, None, 1),
            ("Clear --more-- condition", P.SP_ArrowDown, None, 8),
            ("Break", P.SP_BrowserStop, None, 8),
        ]
        for text, pix, slot, phase in items:
            act = QAction(st.standardIcon(pix), text, self)
            act.setToolTip(text + (f" — {_PENDING[phase]}" if phase else ""))
            act.setEnabled(slot is not None)
            if slot is not None:
                act.triggered.connect(slot)
            tb.addAction(act)
        self.addToolBar(tb)

    def _build_statusbar(self) -> None:
        sb = self.statusBar()
        self.cwd_label = QLabel()
        sb.addWidget(self.cwd_label, 1)
        for flag in ("CAP", "NUM", "OVR"):
            lab = QLabel(flag)
            lab.setEnabled(False)
            sb.addPermanentWidget(lab)

    def _apply_styles(self) -> None:
        self.setStyleSheet("""
            QLabel#PaneTitle { background: #e8e8e8; border-bottom: 1px solid #c8c8c8;
                               padding: 3px 0; font-weight: bold; }
            QLineEdit#Command { border: none; padding: 4px; background: #ffffff; }
        """)

    # -- ações ---------------------------------------------------------------
    def run_command(self, line: str) -> None:
        try:
            rc = self.session.run_command(line)
        except ExitRequest:
            self.close()
            return
        item = QTreeWidgetItem([line, str(rc) if rc else ""])
        if rc:
            for col in (0, 1):
                item.setForeground(col, QColor("#cc0000"))
        self.review.addTopLevelItem(item)
        self.review.scrollToItem(item)
        self.command.setFocus()

    def _filter_review(self, text: str) -> None:
        for k in range(self.review.topLevelItemCount()):
            it = self.review.topLevelItem(k)
            it.setHidden(bool(text) and text.lower() not in it.text(0).lower())

    def choose_do_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Do", os.getcwd(), "Do-files (*.do *.ado);;All files (*)")
        if path:
            self.run_command(f'do "{path}"')

    def choose_directory(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Change working directory", os.getcwd())
        if path:
            self.run_command(f'cd "{path}"')

    def show_about(self) -> None:
        QMessageBox.about(self, "About OpenDTA",
                          f"<b>OpenDTA {__version__}</b><br>"
                          "Interpretador de do-files com sintaxe compatível com o Stata 14.<br>"
                          "Projeto pessoal de estudo. Stata é marca registrada da StataCorp LLC;<br>"
                          "o OpenDTA não é afiliado à StataCorp.")

    def refresh_state(self) -> None:
        cwd = os.getcwd()
        self.cwd_label.setText(cwd)
        self.setWindowTitle(f"OpenDTA {__version__} — {Path(cwd).name or cwd}")
        data_props = {
            "Filename": "", "Variables": str(0), "Observations": str(self.session.nobs),
        }
        for k in range(self._prop_data.childCount()):
            child = self._prop_data.child(k)
            if child.text(0) in data_props:
                child.setText(1, data_props[child.text(0)])

    # -- arrastar e soltar do-files -------------------------------------------
    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:  # noqa: N802
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.lower().endswith((".do", ".ado")):
                self.run_command(f'do "{path}"')
