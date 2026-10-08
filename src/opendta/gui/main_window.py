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

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import (QAction, QColor, QFont, QFontDatabase, QKeySequence,
                           QTextCharFormat, QTextCursor)
from PySide6.QtWidgets import (QApplication, QDockWidget, QFileDialog, QHeaderView, QLabel,
                               QLineEdit, QMainWindow, QMenu, QMessageBox,
                               QPlainTextEdit, QSplitter, QToolBar,
                               QTreeWidget, QTreeWidgetItem, QVBoxLayout,
                               QWidget)

from .. import __version__
from ..core.errors import ExitRequest
from ..session import Session
from .data_browser import DataBrowser
from .variables_manager import VariablesManager
from .viewer import Viewer
from .icons import ACCENT_RED, TOOLBAR_SIZE, app_icon, icon
from .icons import set_theme as set_icon_theme
from .preferences import PreferencesDialog
from .settings import Preferences
from .theme import MONOSPACE_FALLBACKS, STANDARD, ResultsScheme, apply_theme, results_scheme


def _dta_notes(ds) -> list[str]:
    try:
        n = int(ds.chars.get("_dta", {}).get("note0", "0"))
    except ValueError:
        return []
    return [ds.chars["_dta"].get(f"note{k}", "") for k in range(1, n + 1)]


def _human_size(nbytes: int) -> str:
    for unit in ("bytes", "K", "M", "G"):
        if nbytes < 1024 or unit == "G":
            return f"{nbytes:,} {unit}" if unit == "bytes" else f"{nbytes:,.2f}{unit}"
        nbytes /= 1024
    return str(nbytes)


def monospace_font(prefs: Preferences | None = None) -> QFont:
    """Fonte das janelas de texto: a escolhida nas preferências, ou o padrão do
    sistema (Menlo no macOS, Courier New no Windows, DejaVu Sans Mono no Linux)."""
    prefs = prefs or Preferences()
    available = set(QFontDatabase.families())
    candidates = [prefs.font_family] + MONOSPACE_FALLBACKS
    for fam in candidates:
        if fam in available:
            f = QFont(fam, prefs.font_size)
            break
    else:
        f = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        f.setPointSize(prefs.font_size)
    f.setStyleHint(QFont.StyleHint.Monospace)
    return f


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

class ResultsView(QPlainTextEdit):
    MAX_SEGMENTS = 200_000

    def __init__(self, scheme: ResultsScheme = STANDARD):
        super().__init__()
        self.setReadOnly(True)
        self.setUndoRedoEnabled(False)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setFont(monospace_font())
        self.setObjectName("Results")
        self._segments: list[tuple[str, str]] = []   # para redesenhar ao trocar o tema
        self.set_scheme(scheme)

    def set_scheme(self, scheme: ResultsScheme) -> None:
        self._scheme = scheme
        self.setStyleSheet(f"QPlainTextEdit#Results {{ background: {scheme.background}; "
                           f"color: {scheme.colors['text']}; }}")
        self._formats: dict[str, QTextCharFormat] = {}
        for style, color in scheme.colors.items():
            fmt = QTextCharFormat()
            fmt.setForeground(QColor(color))
            if style in scheme.bold:
                fmt.setFontWeight(QFont.Weight.Bold)
            self._formats[style] = fmt
        if self._segments:
            segments, self._segments = self._segments, []
            super().clear()
            for text, style in segments:
                self.append_styled(text, style)

    def clear(self) -> None:  # noqa: D401
        self._segments = []
        super().clear()

    def append_styled(self, text: str, style: str) -> None:
        self._segments.append((text, style))
        if len(self._segments) > self.MAX_SEGMENTS:
            del self._segments[: len(self._segments) // 2]
        cursor = self.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        cursor.insertText(text, self._formats.get(style, self._formats["text"]))
        self.setTextCursor(cursor)
        self.ensureCursorVisible()


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------

class CommandLine(QPlainTextEdit):
    """Janela Command: redimensionável e com várias linhas.

    Enter executa; Shift+Enter quebra a linha. Texto colado com várias
    linhas é executado linha a linha. PgUp/PgDn percorrem o histórico; as
    setas ↑/↓ também, enquanto o comando tiver uma linha só."""

    submitted = Signal(str)

    def __init__(self):
        super().__init__()
        self.setFont(monospace_font())
        self.setObjectName("Command")
        self.setTabChangesFocus(True)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        self.history: list[str] = []
        self._pos = 0
        self._update_min_height()

    def setFont(self, font) -> None:  # noqa: N802
        super().setFont(font)
        self._update_min_height()

    def _update_min_height(self) -> None:
        # pelo menos uma linha visível, qualquer que seja a fonte
        self.setMinimumHeight(self.fontMetrics().lineSpacing() + 14)

    # compatibilidade com a API de QLineEdit usada pela janela
    def text(self) -> str:
        return self.toPlainText()

    def setText(self, text: str) -> None:  # noqa: N802
        self.setPlainText(text)
        self.moveCursor(QTextCursor.MoveOperation.End)

    def insert(self, text: str) -> None:
        self.insertPlainText(text)

    def _submit(self) -> None:
        text = self.toPlainText()
        if not text.strip():
            return
        self.history.append(text)
        self._pos = len(self.history)
        self.clear()
        self.submitted.emit(text)

    def _recall(self, step: int) -> None:
        if not self.history:
            return
        self._pos = min(max(self._pos + step, 0), len(self.history))
        self.setText(self.history[self._pos] if self._pos < len(self.history) else "")

    def keyPressEvent(self, event) -> None:  # noqa: N802
        key = event.key()
        mods = event.modifiers()
        if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
            if mods & Qt.KeyboardModifier.ShiftModifier:
                self.insertPlainText("\n")
            else:
                self._submit()
            return
        single = "\n" not in self.toPlainText()
        if key == Qt.Key.Key_PageUp or (key == Qt.Key.Key_Up and single):
            self._recall(-1)
            return
        if key == Qt.Key.Key_PageDown or (key == Qt.Key.Key_Down and single):
            self._recall(+1)
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
        self.prefs = Preferences()
        self.resize(1280, 800)
        self._build_central()
        self._build_docks()
        self._build_menus()
        self._build_toolbar()
        self._build_statusbar()
        self._apply_styles()

        self.session.settings["hints"] = "on" if self.prefs.hints else "off"
        self.session.ui_hooks["set_permanently"] = self._set_permanently
        self.session.output.add_listener(self.results.append_styled)
        self.session.add_state_listener(self.refresh_state)
        self.refresh_state()
        self.apply_preferences()
        self._restore_layout()
        self.command.setFocus()
        self.setAcceptDrops(True)
        self.browser: DataBrowser | None = None
        self.varmanager: VariablesManager | None = None
        self.session.ui_hooks["browse"] = self.show_browser
        self.session.ui_hooks["varmanage"] = self.show_variables_manager
        self.viewer: Viewer | None = None
        self.session.ui_hooks["help"] = self.show_viewer
        self.session.ui_hooks["view"] = self.show_viewer

    # -- montagem ------------------------------------------------------------
    def _build_central(self) -> None:
        self.results = ResultsView()
        self.command = CommandLine()
        self.command.submitted.connect(self.run_commands)

        # Results e Command separados por uma divisória arrastável nos dois sentidos
        split = QSplitter(Qt.Orientation.Vertical)
        split.setObjectName("CentralSplitter")
        split.setHandleWidth(6)
        split.addWidget(_titled("Results", self.results))
        split.addWidget(_titled("Command", self.command))
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 0)
        split.setCollapsible(0, False)
        split.setCollapsible(1, False)
        split.setSizes([700, 70])
        self.splitter = split
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
        self.review.itemClicked.connect(lambda it, _c: self.command.setText(self._review_text(it)))
        self.review.itemDoubleClicked.connect(lambda it, _c: self.run_commands(self._review_text(it)))
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
        self.var_filter.textChanged.connect(self._filter_variables)
        self.variables.itemSelectionChanged.connect(self._show_variable_properties)
        # duplo clique envia o nome para a janela Command, como no Stata
        self.variables.itemDoubleClicked.connect(
            lambda it, _c: (self.command.insert(it.text(0) + " "), self.command.setFocus()))
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
        self._action(m, "Open...", self.open_dataset, shortcut="Ctrl+O")
        self._action(m, "Save", self.save_dataset, shortcut="Ctrl+S")
        self._action(m, "Save as...", self.save_dataset_as, shortcut="Ctrl+Shift+S")
        m.addSeparator()
        self._action(m, "Do...", self.choose_do_file)
        self._action(m, "Change working directory...", self.choose_directory)
        m.addSeparator()
        log = m.addMenu("Log")
        self._action(log, "Begin...", self.begin_log)
        self._action(log, "Close", lambda: self.run_command("log close"))
        self._action(log, "Suspend", lambda: self.run_command("log off"))
        self._action(log, "Resume", lambda: self.run_command("log on"))
        self._action(log, "View...", self.view_file)
        self._action(m, "Print...", phase=8)
        m.addSeparator()
        self._action(m, "Exit", self.close)

        m = mb.addMenu("&Edit")
        self._action(m, "Copy", self.results.copy, shortcut="Ctrl+C")
        self._action(m, "Copy table", phase=3)
        self._action(m, "Find...", phase=8)
        m.addSeparator()
        self._action(m, "Clear Results", self.results.clear)
        self._action(m, "Preferences...", self.show_preferences, shortcut="Ctrl+,")

        m = mb.addMenu("&Data")
        self._action(m, "Describe data", lambda: self.run_command("describe"))
        self._action(m, "Data Editor", lambda: self.run_command("browse"))
        self._action(m, "Create or change data", phase=1)
        self._action(m, "Variables Manager", lambda: self.run_command("varmanage"))
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
        self._action(m, "Viewer", lambda: self.run_command("help"))
        self._action(m, "Data Editor", lambda: self.run_command("browse"))
        self._action(m, "Do-file Editor", phase=8)

        m = mb.addMenu("&Help")
        self._action(m, "Search...", self.search_help)
        self._action(m, "Contents", lambda: self.run_command("help"))
        m.addSeparator()
        self.act_hints = QAction("Explain Errors (set hints)", self, checkable=True)
        self.act_hints.setToolTip("After an error message, explain what went wrong and how to fix "
                                  "it. Same as typing \"set hints on\" or \"set hints off\".")
        self.act_hints.setStatusTip(self.act_hints.toolTip())
        self.act_hints.triggered.connect(self._toggle_hints)
        m.addAction(self.act_hints)
        self._action(m, "About Error Explanations", self.show_hints_help)
        m.addSeparator()
        self._action(m, "About OpenDTA", self.show_about)

        for menu in mb.findChildren(QMenu):
            menu.setToolTipsVisible(True)

    def _build_toolbar(self) -> None:
        tb = QToolBar("Toolbar")
        tb.setObjectName("Toolbar")
        tb.setMovable(False)
        tb.setIconSize(TOOLBAR_SIZE)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
        # mesma disposição da barra do Stata 14; desenho próprio (gui/icons/)
        groups = [
            [("open", "Open", self.open_dataset, 1), ("save", "Save", self.save_dataset, 1),
             ("print", "Print", None, 8)],
            [("log", "Log (begin or close)", self.toggle_log, 2),
             ("viewer", "Viewer", lambda: self.run_command("help"), 2), ("graph", "Graph", None, 7),
             ("dofile", "Do-file Editor", None, 8)],
            [("dataeditor", "Data Editor (Edit)", lambda: self.run_command("edit"), 1),
             ("databrowser", "Data Browser (Browse)", lambda: self.run_command("browse"), 1),
             ("variables", "Variables Manager", lambda: self.run_command("varmanage"), 1)],
            [("more", "Clear --more-- condition", None, 8), ("break", "Break", None, 8)],
        ]
        self.toolbar_actions: dict[str, QAction] = {}
        for gi, group in enumerate(groups):
            if gi:
                tb.addSeparator()
            for name, text, slot, phase in group:
                accent = ACCENT_RED if name == "break" else None
                act = QAction(icon(name, accent), text, self)
                act.setData((name, accent))
                act.setToolTip(text + (f" — {_PENDING[phase]}" if phase and slot is None else ""))
                act.setEnabled(slot is not None)
                if slot is not None:
                    act.triggered.connect(slot)
                tb.addAction(act)
                self.toolbar_actions[name] = act
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
            QLabel#PaneTitle { background: palette(window); border-bottom: 1px solid palette(mid);
                               padding: 4px 0; font-weight: 600; }
            QPlainTextEdit#Command { border: none; padding: 3px; background: palette(base); }
            QSplitter#CentralSplitter::handle { background: palette(window); }
            QSplitter#CentralSplitter::handle:hover { background: palette(mid); }
            QToolBar#Toolbar { border: none; spacing: 2px; padding: 3px 6px; }
            QToolBar#Toolbar QToolButton { border: none; border-radius: 6px; padding: 4px; }
            QToolBar#Toolbar QToolButton:hover { background: rgba(127, 127, 127, 0.16); }
            QToolBar#Toolbar QToolButton:pressed { background: rgba(127, 127, 127, 0.28); }
            QToolBar#Toolbar::separator { width: 1px; margin: 5px 6px; background: rgba(127, 127, 127, 0.3); }
        """)

    # -- ações ---------------------------------------------------------------
    def run_commands(self, text: str) -> None:
        """Texto da janela Command. Uma linha: comando comum. Várias linhas:
        rodam juntas como um trecho de do-file (blocos { } inteiros)."""
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if len(lines) <= 1:
            self.run_command(lines[0] if lines else text)
            return
        try:
            rc = self.session.run_text("\n".join(lines))
        except ExitRequest:
            self.close()
            return
        self._add_review(text, rc)
        self.command.setFocus()

    def _add_review(self, text: str, rc: int) -> None:
        first = text.strip().splitlines()[0]
        shown = first + (" …" if "\n" in text.strip() else "")
        item = QTreeWidgetItem([shown, str(rc) if rc else ""])
        item.setData(0, Qt.ItemDataRole.UserRole, text)
        item.setToolTip(0, text)
        if rc:
            err = QColor(results_scheme(self.prefs.theme).colors["error"])
            for col in (0, 1):
                item.setForeground(col, err)
        self.review.addTopLevelItem(item)
        self.review.scrollToItem(item)

    def run_command(self, line: str) -> None:
        try:
            rc = self.session.run_command(line)
        except ExitRequest:
            self.close()
            return
        self._add_review(line, rc)
        self.command.setFocus()

    @staticmethod
    def _review_text(item) -> str:
        full = item.data(0, Qt.ItemDataRole.UserRole)
        return full if full else item.text(0)

    def _filter_review(self, text: str) -> None:
        for k in range(self.review.topLevelItemCount()):
            it = self.review.topLevelItem(k)
            it.setHidden(bool(text) and text.lower() not in it.text(0).lower())

    # -- arquivos de dados ---------------------------------------------------
    def _confirm_discard(self) -> bool:
        """Dados alterados: oferece salvar antes de substituir. False = cancelar."""
        ds = self.session.data
        if not (ds.changed and ds.nvars):
            return True
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Icon.Warning)
        box.setText("Data in memory have changed.")
        box.setInformativeText("Do you want to save the changes before opening another dataset?")
        box.setStandardButtons(QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard
                               | QMessageBox.StandardButton.Cancel)
        box.setDefaultButton(QMessageBox.StandardButton.Save)
        answer = box.exec()
        if answer == QMessageBox.StandardButton.Cancel:
            return False
        if answer == QMessageBox.StandardButton.Save:
            return self.save_dataset()
        return True

    def open_dataset(self, path: str | None = None) -> None:
        if not path:
            path, _ = QFileDialog.getOpenFileName(self, "Open", os.getcwd(),
                                                  "Stata data (*.dta);;All files (*)")
        if not path or not self._confirm_discard():
            return
        self.run_command(f'use "{path}", clear')

    def save_dataset(self) -> bool:
        ds = self.session.data
        if not ds.fullpath:
            return self.save_dataset_as()
        self.run_command(f'save "{ds.fullpath}", replace')
        return self.session.rc == 0

    def save_dataset_as(self) -> bool:
        ds = self.session.data
        start = ds.filename or os.path.join(os.getcwd(), "untitled.dta")
        path, _ = QFileDialog.getSaveFileName(self, "Save as", start, "Stata data (*.dta)")
        if not path:
            return False
        if not path.lower().endswith(".dta"):
            path += ".dta"
        # o diálogo já confirmou a substituição, por isso replace
        self.run_command(f'save "{path}", replace')
        return self.session.rc == 0

    def show_browser(self, columns=None, rows=None, *, edit: bool = False,
                     nolabel: bool = False) -> None:
        if self.browser is None:
            self.browser = DataBrowser(self.session, monospace_font(self.prefs), self,
                                       run=self._run_from_editor)
            self.browser.setWindowFlag(Qt.WindowType.Window, True)
            self.browser.set_dark(self.prefs.theme == "dark")
        self.browser.model.nolabel = nolabel
        self.browser.model.set_subset(columns, rows)
        self.browser.set_mode(edit)
        self.browser.show()
        self.browser.raise_()
        self.browser.activateWindow()

    # -- Viewer e log ---------------------------------------------------------------
    def show_viewer(self, path, title: str) -> None:
        if self.viewer is None:
            font = monospace_font(self.prefs)
            self.viewer = Viewer(self.session, self, font_family=font.family(),
                                 font_size=max(9, font.pointSize()), dark=self.prefs.theme == "dark",
                                 run=self._run_from_editor)
            self.viewer.setWindowFlag(Qt.WindowType.Window, True)
        self.viewer.dark = self.prefs.theme == "dark"
        self.viewer.show_page(path, title)
        self.viewer.show()
        self.viewer.raise_()
        self.viewer.activateWindow()

    def search_help(self) -> None:
        self.run_command("help")
        if self.viewer is not None:
            self.viewer.field.setFocus()
            self.viewer.field.selectAll()

    def begin_log(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Begin log", os.path.join(os.getcwd(), "log.smcl"),
                                              "SMCL log (*.smcl);;Text log (*.log)")
        if path:
            # o diálogo já confirmou a substituição
            opts = "replace text" if path.lower().endswith(".log") else "replace"
            self.run_command(f'log using "{path}", {opts}')

    def toggle_log(self) -> None:
        if getattr(self.session, "logs", None):
            self.run_command("log close")
        else:
            self.begin_log()

    def view_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "View", os.getcwd(),
                                              "Logs and help (*.smcl *.log *.sthlp *.txt);;All files (*)")
        if path:
            self.run_command(f'view "{path}"')

    def show_variables_manager(self) -> None:
        if self.varmanager is None:
            self.varmanager = VariablesManager(self.session, monospace_font(self.prefs), self,
                                               run=self._run_from_editor)
            self.varmanager.setWindowFlag(Qt.WindowType.Window, True)
        self.varmanager.refresh()
        self.varmanager.show()
        self.varmanager.raise_()
        self.varmanager.activateWindow()

    def _run_from_editor(self, line: str) -> int:
        """Comando gerado pelo Data Editor ou pelo Variables Manager: ecoa em
        Results e entra em Review, como no Stata."""
        self.run_command(line)
        return self.session.rc

    def choose_do_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Do", os.getcwd(), "Do-files (*.do *.ado);;All files (*)")
        if path:
            self.run_command(f'do "{path}"')

    def choose_directory(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Change working directory", os.getcwd())
        if path:
            self.run_command(f'cd "{path}"')

    def show_preferences(self) -> None:
        dlg = PreferencesDialog(self.prefs, self)
        dlg.applied.connect(self.apply_preferences)
        dlg.exec()

    def apply_preferences(self) -> None:
        """Aplica tema, fonte e ícone escolhidos (também chamado na abertura)."""
        self.apply_theme(self.prefs.theme)
        font = monospace_font(self.prefs)
        for w in (self.results, self.command, self.review):
            w.setFont(font)
        if getattr(self, "browser", None) is not None:
            self.browser.table.setFont(font)
        ico = app_icon(self.prefs.app_icon)
        self.setWindowIcon(ico)
        app = QApplication.instance()
        if app is not None:
            app.setWindowIcon(ico)

    # -- layout salvo entre sessões -------------------------------------------
    def _restore_layout(self) -> None:
        geo, state, split = self.prefs.layout()
        if geo:
            self.restoreGeometry(geo)
        if state:
            self.restoreState(state)
        if split:
            self.splitter.restoreState(split)

    def closeEvent(self, event) -> None:  # noqa: N802
        self.prefs.save_layout(self.saveGeometry(), self.saveState(), self.splitter.saveState())
        super().closeEvent(event)

    def apply_theme(self, theme: str) -> None:
        """Tema claro ou escuro da interface inteira, independente do sistema."""
        app = QApplication.instance()
        if app is not None:
            apply_theme(app, theme)
        set_icon_theme(theme)
        self.results.set_scheme(results_scheme(theme))
        err = QColor(results_scheme(theme).colors["error"])
        for k in range(self.review.topLevelItemCount()):
            it = self.review.topLevelItem(k)
            if it.text(1):
                for col in (0, 1):
                    it.setForeground(col, err)
        for act in getattr(self, "toolbar_actions", {}).values():
            name, accent = act.data()
            act.setIcon(icon(name, accent))
        self._apply_styles()
        if getattr(self, "browser", None) is not None:
            self.browser.set_dark(theme == "dark")

    def _set_permanently(self, name: str, value: str) -> None:
        """set ..., permanently: lembrado na próxima abertura (só hints, por ora)."""
        if name == "hints":
            self.prefs.hints = value == "on"

    def _toggle_hints(self, on: bool) -> None:
        # vira comando, para ficar em Results e Review como qualquer set
        self.run_command(f"set hints {'on' if on else 'off'}, permanently")

    def show_hints_help(self) -> None:
        QMessageBox.information(
            self, "Error explanations",
            "<p>When a command fails, OpenDTA shows the usual error message and return "
            "code <tt>r(#)</tt>, and between them a short explanation of what went wrong "
            "in that command and a hint on how to fix it:</p>"
            "<pre>. replace nota = 9 in 2\nObs. nos. out of range\n"
            "  \u2192 \"in 2\" asks for observations that do not exist; ...\n"
            "    Hint: create observations with \"set obs #\" ...\nr(198);</pre>"
            "<p>The return codes do not change, so <tt>capture</tt> and <tt>_rc</tt> work as "
            "before. Errors raised on purpose with <tt>error #</tt> are not explained.</p>"
            "<p><b>set hints on</b> &nbsp;turns the explanations on (default)<br>"
            "<b>set hints off</b> &nbsp;shows only the error message and return code</p>"
            "<p>Add <b>, permanently</b> to remember the choice the next time OpenDTA "
            "opens. The menu item <i>Help \u2192 Explain Errors</i> does that.</p>")

    def show_about(self) -> None:
        QMessageBox.about(self, "About OpenDTA",
                          f"<b>OpenDTA {__version__}</b><br>"
                          "Interpretador de do-files com sintaxe compatível com o Stata 14.<br>"
                          "Projeto pessoal de estudo. Stata é marca registrada da StataCorp LLC;<br>"
                          "o OpenDTA não é afiliado à StataCorp.")

    def refresh_state(self) -> None:
        if hasattr(self, "act_hints"):
            self.act_hints.setChecked(self.session.settings.get("hints", "on") == "on")
        cwd = os.getcwd()
        self.cwd_label.setText(cwd)
        ds = self.session.data
        title = Path(ds.filename).name if ds.filename else (Path(cwd).name or cwd)
        self.setWindowTitle(f"OpenDTA {__version__} — {title}")

        # Variables
        selected = self._selected_variable()
        self.variables.clear()
        for v in ds.vars:
            item = QTreeWidgetItem([v.name, v.label])
            self.variables.addTopLevelItem(item)
            if v.name == selected:
                item.setSelected(True)
        self._filter_variables(self.var_filter.text())

        # Properties > Data
        size = ds.width() * ds.nobs
        data_props = {
            "Filename": Path(ds.filename).name if ds.filename else "",
            "Label": ds.label,
            "Notes": str(len(_dta_notes(ds))) if _dta_notes(ds) else "",
            "Variables": f"{ds.nvars:,}",
            "Observations": f"{ds.nobs:,}",
            "Size": _human_size(size),
            "Memory": _human_size(sum(v.nbytes for v in ds.vars)),
            "Sorted by": " ".join(ds.sortlist),
        }
        for k in range(self._prop_data.childCount()):
            child = self._prop_data.child(k)
            child.setText(1, data_props.get(child.text(0), ""))
        self._show_variable_properties()
        # as janelas auxiliares se atualizam depois do comando terminar
        # (a edição de uma célula dispara comandos de dentro do próprio modelo)
        QTimer.singleShot(0, self._refresh_aux_windows)

    def _refresh_aux_windows(self) -> None:
        if getattr(self, "browser", None) is not None and self.browser.isVisible():
            self.browser.refresh()
        if getattr(self, "varmanager", None) is not None and self.varmanager.isVisible():
            self.varmanager.refresh()

    def _selected_variable(self) -> str:
        items = self.variables.selectedItems()
        return items[0].text(0) if items else ""

    def _filter_variables(self, text: str) -> None:
        for k in range(self.variables.topLevelItemCount()):
            it = self.variables.topLevelItem(k)
            it.setHidden(bool(text) and text.lower() not in (it.text(0) + " " + it.text(1)).lower())

    def _show_variable_properties(self) -> None:
        name = self._selected_variable()
        ds = self.session.data
        props = {}
        if name and ds.has(name):
            v = ds.get(name)
            props = {"Name": v.name, "Label": v.label, "Type": v.vtype, "Format": v.fmt,
                     "Value label": v.value_label, "Notes": ""}
        for k in range(self._prop_vars.childCount()):
            child = self._prop_vars.child(k)
            child.setText(1, props.get(child.text(0), ""))

    # -- arrastar e soltar do-files -------------------------------------------
    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:  # noqa: N802
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.lower().endswith((".do", ".ado")):
                self.run_command(f'do "{path}"')
            elif path.lower().endswith(".dta"):
                self.open_dataset(path)
                break
