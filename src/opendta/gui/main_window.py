"""Main window in the Stata 14 layout.

    ┌──────────────────────────────────────────────────────────┐
    │ Menus: File Edit Data Graphics Statistics User Window Help│
    │ Toolbar                                                   │
    ├──────────┬───────────────────────────────┬───────────────┤
    │ Review   │ Results                       │ Variables     │
    │          │                               ├───────────────┤
    │          ├───────────────────────────────┤ Properties    │
    │          │ Command                       │               │
    ├──────────┴───────────────────────────────┴───────────────┤
    │ working directory                           CAP  NUM  OVR │
    └──────────────────────────────────────────────────────────┘

Icons and visual identity are OpenDTA's own. Menu items not yet
implemented appear disabled, with a tooltip explaining why.
"""

from __future__ import annotations

import os
import re
import time
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QColor, QFont, QKeySequence
from PySide6.QtWidgets import (QApplication, QDockWidget, QFileDialog, QHeaderView, QLabel,
                               QLineEdit, QMainWindow, QMenu, QMessageBox, QSizePolicy, QSplitter, QToolBar, QToolButton,
                               QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget)
from shiboken6 import isValid

from .. import __version__
from .. import i18n
from . import qt_i18n
from ..i18n import tr
from ..core.errors import ExitRequest
from ..session import Session
from .data_browser import DataBrowser
from .variables_manager import VariablesManager
from .viewer import Viewer
from .graph_window import GraphWindow
from .dofile_editor import DoFileEditor
from .icons import ACCENT_RED, TOOLBAR_SIZE, app_icon, icon
from .icons import set_theme as set_icon_theme
from .menu_templates import DATA_TEMPLATES, GRAPHICS_TEMPLATES, STATISTICS_TEMPLATES
from .panes import CommandLine, ResultsView, monospace_font, titled
from .preferences import PreferencesDialog
from .settings import Preferences
from .theme import apply_theme, results_scheme
from .window_cmds import build_user_menu
from .window_cmds import install as install_window_hooks


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


# short name shown under each toolbar icon, as in Stata (the tooltip keeps the long text)
_TOOLBAR_LABELS = {
    "open": "Open", "save": "Save", "print": "Print", "log": "Log", "viewer": "Viewer", "graph": "Graph",
    "dofile": "Do-file Editor", "dataeditor": "Data Editor", "databrowser": "Data Browser",
    "varmanager": "Manager", "review": "Review", "variables": "Variables", "properties": "Properties",
    "more": "More", "break": "Break",
}

# toolbar buttons that exist in Stata but do nothing here, with the reason in their tooltip
_PENDING = {
    "more": "No effect: OpenDTA always runs with set more off",
    "break": "Not available yet: commands run in the interface thread",
}


def _present(window: QWidget) -> None:
    window.show()
    window.raise_()
    window.activateWindow()


def _paint(item: QTreeWidgetItem, color: QColor) -> None:
    for col in (0, 1):
        item.setForeground(col, color)


class MainWindow(QMainWindow):
    def __init__(self, session: Session | None = None):
        super().__init__()
        self.session = session or Session()
        self.session.input_provider = self._request_input
        self.prefs = Preferences()
        from .update_dialog import UpdateController
        self.updates = UpdateController(self)
        self._pending_update = None
        # windows created on first use
        self.browser: DataBrowser | None = None
        self.varmanager: VariablesManager | None = None
        self.viewer: Viewer | None = None
        self.graph_window: GraphWindow | None = None
        self.dofile_editor: DoFileEditor | None = None
        self._settings: PreferencesDialog | None = None
        self._last_find = ""
        i18n.set_language(self.prefs.language)
        i18n.set_messages(self.prefs.messages)
        qt_i18n.apply(i18n.language())
        self.resize(1280, 800)
        self._build_central()
        self._build_docks()
        self._build_menus()
        self._build_toolbar()
        self._build_statusbar()
        self._apply_styles()

        self.session.settings["hints"] = "on" if self.prefs.hints else "off"
        if not self.prefs.superscript:
            self.session.run("quietly set superscript off")
        if self.prefs.numerics == "precise":
            self.session.run("quietly set numerics precise")
        self.session.ui_hooks["set_permanently"] = self._set_permanently
        self.session.ui_hooks["language"] = lambda: QTimer.singleShot(0, self.retranslate)
        self.session.output.add_listener(self.results.append_styled)
        self.session.add_state_listener(self.refresh_state)
        self._startup_session()
        self.refresh_state()
        self.apply_preferences()
        self._restore_layout()
        self.command.setFocus()
        self.setAcceptDrops(True)
        self.session.ui_hooks.update(browse=self.show_browser, varmanage=self.show_variables_manager,
                                     help=self.show_viewer, view=self.show_viewer, graph=self.show_graph,
                                     graph_close=self.close_graph, doedit=self.show_dofile_editor)
        install_window_hooks(self)

    # -- montagem ------------------------------------------------------------
    def _build_central(self) -> None:
        self.results = ResultsView()
        self.command = CommandLine()
        self.command.submitted.connect(self.run_commands)
        self.results.linkActivated.connect(self.run_command)    # "help topic" of an explanation

        # Results and Command separated by a splitter that can be dragged both ways
        split = QSplitter(Qt.Orientation.Vertical)
        split.setObjectName("CentralSplitter")
        split.setHandleWidth(6)
        split.addWidget(titled("Results", self.results))
        split.addWidget(titled("Command", self.command))
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 0)
        split.setCollapsible(0, False)
        split.setCollapsible(1, False)
        split.setSizes([700, 70])
        self.splitter = split
        self.setCentralWidget(split)

    def _dock(self, title: str, widget: QWidget, area: Qt.DockWidgetArea) -> QDockWidget:
        d = QDockWidget(tr(title), self)
        d.setObjectName(title)          # the object name stays in English (saved layout)
        d.setWidget(widget)
        d.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable
                      | QDockWidget.DockWidgetFeature.DockWidgetMovable
                      | QDockWidget.DockWidgetFeature.DockWidgetFloatable)
        self.addDockWidget(area, d)
        return d

    def _build_docks(self) -> None:
        # Review
        self.review = QTreeWidget()
        self.review.setHeaderLabels([tr("Command"), "_rc"])
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
        self.review_filter.setPlaceholderText(tr("Filter commands here"))
        self.review_filter.textChanged.connect(self._filter_review)
        lay.addWidget(self.review_filter)
        lay.addWidget(self.review)
        self.dock_review = self._dock("Review", review_box, Qt.DockWidgetArea.LeftDockWidgetArea)

        # Variables
        self.variables = QTreeWidget()
        self.variables.setHeaderLabels([tr("Name"), tr("Label")])
        self.variables.setRootIsDecorated(False)
        var_box = QWidget()
        lay = QVBoxLayout(var_box)
        lay.setContentsMargins(0, 0, 0, 0)
        self.var_filter = QLineEdit()
        self.var_filter.setPlaceholderText(tr("Filter variables here"))
        self.var_filter.textChanged.connect(self._filter_variables)
        self.variables.itemSelectionChanged.connect(self._show_variable_properties)
        # double click sends the name to the Command window, as in Stata
        self.variables.itemDoubleClicked.connect(
            lambda it, _c: (self.command.insert(it.text(0) + " "), self.command.setFocus()))
        lay.addWidget(self.var_filter)
        lay.addWidget(self.variables)
        self.dock_variables = self._dock("Variables", var_box, Qt.DockWidgetArea.RightDockWidgetArea)

        # Properties
        self.properties = QTreeWidget()
        self.properties.setHeaderLabels([tr("Property"), tr("Value")])
        self.properties.header().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self._prop_vars = QTreeWidgetItem([tr("Variables")])
        self._prop_data = QTreeWidgetItem([tr("Data")])
        self._prop_vars.setData(0, Qt.ItemDataRole.UserRole, "Variables")
        self._prop_data.setData(0, Qt.ItemDataRole.UserRole, "Data")
        for parent, keys in ((self._prop_vars, ["Name", "Label", "Type", "Format", "Value label", "Notes"]),
                             (self._prop_data, ["Filename", "Label", "Notes", "Variables",
                                                "Observations", "Size", "Memory", "Sorted by"])):
            for k in keys:
                item = QTreeWidgetItem([tr(k), ""])
                item.setData(0, Qt.ItemDataRole.UserRole, k)     # key in English
                parent.addChild(item)
            self.properties.addTopLevelItem(parent)
            parent.setExpanded(True)
        self.dock_properties = self._dock("Properties", self.properties, Qt.DockWidgetArea.RightDockWidgetArea)

        self.splitDockWidget(self.dock_variables, self.dock_properties, Qt.Orientation.Vertical)
        self.resizeDocks([self.dock_review], [260], Qt.Orientation.Horizontal)
        self.resizeDocks([self.dock_variables], [300], Qt.Orientation.Horizontal)
        self.resizeDocks([self.dock_variables, self.dock_properties], [420, 330], Qt.Orientation.Vertical)
        # chosen again in the Window menu, a closed window comes back to its fixed place
        for dock in (self.dock_review, self.dock_variables, self.dock_properties):
            dock.toggleViewAction().triggered.connect(lambda on, d=dock: on and self._place_dock(d))

    def _action(self, menu, text: str, slot, *, shortcut: str | None = None) -> QAction:
        act = QAction(tr(text), self)
        if shortcut:
            act.setShortcut(QKeySequence(shortcut))
        act.triggered.connect(slot)
        menu.addAction(act)
        return act

    def _build_menus(self) -> None:
        mb = self.menuBar()

        m = mb.addMenu(tr("&File"))
        self._action(m, "Open...", self.open_dataset, shortcut="Ctrl+O")
        self._action(m, "Save", self.save_dataset, shortcut="Ctrl+S")
        self._action(m, "Save as...", self.save_dataset_as, shortcut="Ctrl+Shift+S")
        m.addSeparator()
        self._action(m, "Do...", self.choose_do_file)
        self._action(m, "Change working directory...", self.choose_directory)
        m.addSeparator()
        log = m.addMenu(tr("Log"))
        self._action(log, "Begin...", self.begin_log)
        self._action(log, "Close", lambda: self.run_command("log close"))
        self._action(log, "Suspend", lambda: self.run_command("log off"))
        self._action(log, "Resume", lambda: self.run_command("log on"))
        self._action(log, "View...", self.view_file)
        self._action(m, "Print...", self.print_results, shortcut="Ctrl+P")
        m.addSeparator()
        self._action(m, "Exit", self.close)

        m = mb.addMenu(tr("&Edit"))
        self._action(m, "Copy", self.results.copy, shortcut="Ctrl+C")
        self._action(m, "Copy table", self.copy_table, shortcut="Ctrl+Shift+C")
        self._action(m, "Find...", self.find_results, shortcut="Ctrl+F")
        m.addSeparator()
        self._action(m, "Clear Results", self.results.clear)
        self._action(m, "Settings...", lambda: self.show_preferences(), shortcut="Ctrl+,")

        m = mb.addMenu(tr("&Data"))
        self._action(m, "Describe data", lambda: self.run_command("describe"))
        self._action(m, "Data Editor", lambda: self.run_command("browse"))
        for group, items in DATA_TEMPLATES:
            sub = m.addMenu(tr(group))
            sub.setToolTipsVisible(True)
            for item, cmd in items:
                self._template_action(sub, item, cmd)

        # items without dialogs write the command to the Command window, ready to be completed with
        # the variables
        m = mb.addMenu(tr("&Graphics"))
        for item, cmd in GRAPHICS_TEMPLATES:
            if item is None:
                m.addSeparator()
                continue
            self._template_action(m, item, cmd)

        m = mb.addMenu(tr("&Statistics"))
        for group, items in STATISTICS_TEMPLATES:
            sub = m.addMenu(tr(group))
            for item, cmd in items:
                self._template_action(sub, item, cmd)
            sub.setToolTipsVisible(True)

        self.user_menu = mb.addMenu(tr("&User"))
        if getattr(self, "_user_items", None):
            build_user_menu(self.user_menu, self._user_items, self.run_command)

        m = mb.addMenu(tr("&Window"))
        for dock in (self.dock_review, self.dock_variables, self.dock_properties):
            m.addAction(dock.toggleViewAction())
        m.addSeparator()
        self._action(m, "Command", lambda: self.command.setFocus())
        self._action(m, "Results", lambda: self.results.setFocus())
        self._action(m, "Graph", self.show_graph_window)
        self._action(m, "Viewer", lambda: self.run_command("help"))
        self._action(m, "Data Editor", lambda: self.run_command("browse"))
        self._action(m, "Do-file Editor", lambda: self.show_dofile_editor(""), shortcut="Ctrl+9")

        m = mb.addMenu(tr("&Help"))
        self._action(m, "Search...", self.search_help)
        self._action(m, "Contents", lambda: self.run_command("help"))
        m.addSeparator()
        self.act_hints = QAction(tr("Explain Errors (set hints)"), self, checkable=True)
        self.act_hints.setToolTip(tr("After an error message, explain what went wrong and how to fix "
                                     "it. Same as typing \"set hints on\" or \"set hints off\"."))
        self.act_hints.setStatusTip(self.act_hints.toolTip())
        self.act_hints.triggered.connect(self._toggle_hints)
        m.addAction(self.act_hints)
        self._action(m, "About Error Explanations", self.show_hints_help)
        self.act_superscript = QAction(tr("Superscript Exponents (set superscript)"), self, checkable=True)
        self.act_superscript.setToolTip(tr("Show scientific notation as 1.0000\u00d710\u00b9\u2070 "
                                           "instead of 1.00000e+10. Same as typing \"set superscript "
                                           "on\" or \"set superscript off\"."))
        self.act_superscript.setStatusTip(self.act_superscript.toolTip())
        self.act_superscript.triggered.connect(self._toggle_superscript)
        m.addAction(self.act_superscript)
        m.addSeparator()
        self._action(m, "Check for Updates...", lambda: self.updates.check_now())
        self._action(m, "About OpenDTA", self.show_about)

        for menu in mb.findChildren(QMenu):
            menu.setToolTipsVisible(True)

    def _template_action(self, menu, item: str, cmd: str) -> QAction:
        """Menu item for a command: opens the dialog box, when there is one,
        or writes the command to the Command window, ready to be completed."""
        from .dialogs import spec_for
        spec = spec_for(cmd)
        if spec is not None:
            act = self._action(menu, item, lambda _=False, sp=spec, c=cmd: self.open_dialog(sp, c))
            act.setText(act.text() + "...")
            act.setToolTip(tr("Open a dialog to build the command"))
        else:
            act = self._action(menu, item, lambda _=False, c=cmd: self.insert_command(c))
            act.setToolTip(tr('Write "{cmd}" in the Command window').format(cmd=cmd.strip()))
        return act

    def open_dialog(self, spec, template: str = "") -> None:
        from .dialogs import CommandDialog
        dlg = CommandDialog(spec, self.session, self, run=self.run_command,
                            to_command=self.insert_command, template=template)
        dlg.show()
        dlg.raise_()

    def retranslate(self) -> None:
        """Switches the interface language without restarting: rebuilds menus and texts."""
        qt_i18n.apply(i18n.language())
        self.menuBar().clear()
        self._build_menus()
        self.refresh_state()
        for act in self.toolbar_actions.values():
            key = act.property("key")
            act.setText(tr(act.property("label")))
            act.setToolTip(tr(key) + (act.property("suffix") or ""))
        self._fit_toolbar_buttons()
        for d in (self.dock_review, self.dock_variables, self.dock_properties):
            d.setWindowTitle(tr(d.objectName()))
        for lab in self.findChildren(QLabel, "PaneTitle"):
            lab.setText(f"  {tr(lab.property('key'))}")
        self.review.setHeaderLabels([tr("Command"), "_rc"])
        self.variables.setHeaderLabels([tr("Name"), tr("Label")])
        self.properties.setHeaderLabels([tr("Property"), tr("Value")])
        self.review_filter.setPlaceholderText(tr("Filter commands here"))
        self.var_filter.setPlaceholderText(tr("Filter variables here"))
        for parent in (self._prop_vars, self._prop_data):
            parent.setText(0, tr(parent.data(0, Qt.ItemDataRole.UserRole)))
            for k in range(parent.childCount()):
                ch = parent.child(k)
                ch.setText(0, tr(ch.data(0, Qt.ItemDataRole.UserRole)))
        for name in ("dofile_editor", "graph_window", "browser", "varmanager", "viewer"):
            w = getattr(self, name, None)
            if w is not None and hasattr(w, "retranslate"):
                w.retranslate()

    def _build_toolbar(self) -> None:
        tb = QToolBar("Toolbar")
        tb.setObjectName("Toolbar")
        tb.setMovable(False)
        tb.setIconSize(TOOLBAR_SIZE)
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        # same arrangement as Stata 14's toolbar; OpenDTA's own drawing (gui/icons/)
        # buttons without a slot are the _PENDING ones: shown disabled, with the reason
        groups = [
            [("open", "Open", self.open_dataset), ("save", "Save", self.save_dataset),
             ("print", "Print", self.print_results)],
            [("log", "Log (begin or close)", self.toggle_log),
             ("viewer", "Viewer", lambda: self.run_command("help")), ("graph", "Graph", self.show_graph_window),
             ("dofile", "Do-file Editor", lambda: self.show_dofile_editor(""))],
            [("dataeditor", "Data Editor (Edit)", lambda: self.run_command("edit")),
             ("databrowser", "Data Browser (Browse)", lambda: self.run_command("browse")),
             ("varmanager", "Variables Manager", lambda: self.run_command("varmanage"))],
            # the side windows: each button appears only while its window is closed and brings it back
            [(name, f"Show the {title} window", lambda _=False, d=dock: d.toggleViewAction().trigger())
             for name, title, dock in (("review", "Review", self.dock_review),
                                       ("variables", "Variables", self.dock_variables),
                                       ("properties", "Properties", self.dock_properties))],
            [("more", "Clear --more-- condition", None), ("break", "Break", None)],
        ]
        self.toolbar_actions: dict[str, QAction] = {}
        for gi, group in enumerate(groups):
            if gi:
                sep = tb.addSeparator()
                if gi == 3:
                    panel_separator = sep
            for name, text, slot in group:
                accent = ACCENT_RED if name == "break" else None
                label = _TOOLBAR_LABELS[name]
                act = QAction(icon(name, accent), tr(label), self)
                act.setData((name, accent))
                act.setProperty("label", label)
                suffix = f" — {tr(_PENDING[name])}" if slot is None else ""
                act.setProperty("key", text)
                act.setProperty("suffix", suffix)
                act.setToolTip(tr(text) + suffix)
                act.setEnabled(slot is not None)
                if slot is not None:
                    act.triggered.connect(slot)
                tb.addAction(act)
                self.toolbar_actions[name] = act
        # Settings on the far right, away from the Stata-style groups
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        tb.addWidget(spacer)
        act = QAction(icon("settings"), tr("Settings"), self)
        act.setData(("settings", None))
        act.setProperty("label", "Settings")
        act.setProperty("key", "Settings")
        suffix = f"  [{QKeySequence('Ctrl+,').toString(QKeySequence.SequenceFormat.NativeText)}]"
        act.setProperty("suffix", suffix)
        act.setToolTip(tr("Settings") + suffix)
        act.triggered.connect(lambda: self.show_preferences())
        tb.addAction(act)
        self.toolbar_actions["settings"] = act
        self.addToolBar(tb)
        self._fit_toolbar_buttons()
        self._link_dock_buttons(panel_separator)

    def _link_dock_buttons(self, separator: QAction) -> None:
        """The Review, Variables and Properties buttons show only while their window is closed (by its
        X or by the Window menu); open, the window is there and the button leaves the toolbar."""
        buttons = [(self.toolbar_actions[name], dock) for name, dock in
                   (("review", self.dock_review), ("variables", self.dock_variables),
                    ("properties", self.dock_properties))]

        def sync(act: QAction, closed: bool) -> None:
            if not all(isValid(a) for a in (separator, *(b for b, _ in buttons))):
                return                      # the window is being destroyed and takes its toolbar along
            act.setVisible(closed)
            separator.setVisible(any(a.isVisible() for a, _ in buttons))

        for act, dock in buttons:
            dock.toggleViewAction().toggled.connect(lambda on, a=act: sync(a, not on))
            sync(act, dock.isHidden())

    def _fit_toolbar_buttons(self) -> None:
        """Small font and a width that follows the label, the same on every platform style."""
        tb = self.findChild(QToolBar, "Toolbar")
        font = QFont(self.font())
        font.setPixelSize(11)
        for act in tb.actions():
            btn = tb.widgetForAction(act)
            if isinstance(btn, QToolButton):
                btn.setFont(font)
                btn.setFixedWidth(max(TOOLBAR_SIZE.width() + 16, btn.fontMetrics().horizontalAdvance(act.text()) + 16))

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
            QToolBar#Toolbar QToolButton { border: none; border-radius: 6px; padding: 3px 4px; }
            QToolBar#Toolbar QToolButton:hover { background: rgba(127, 127, 127, 0.16); }
            QToolBar#Toolbar QToolButton:pressed { background: rgba(127, 127, 127, 0.28); }
            QToolBar#Toolbar::separator { width: 1px; margin: 5px 6px; background: rgba(127, 127, 127, 0.3); }
        """)

    def _request_input(self) -> str | None:
        """display ..., _request(): the typed line, or None if cancelled."""
        from PySide6.QtWidgets import QInputDialog
        text, ok = QInputDialog.getText(self, "OpenDTA", tr("Input requested by display:"))
        return text if ok else None

    # -- actions ---------------------------------------------------------------
    def run_commands(self, text: str) -> None:
        """Text of the Command window. One line: an ordinary command. Several lines:
        they run together as a piece of a do-file (whole { } blocks)."""
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

    def _add_review(self, text: str, rc: int, children: list | None = None) -> None:
        item = self._review_item(text, rc, children or [])
        if children:
            # Settings > General > Session: the commands of the do-files, under the do line
            label = self._do_label(text)
            if label:
                item.setText(0, label)
            self.review.setRootIsDecorated(True)
        self.review.addTopLevelItem(item)
        self.review.scrollToItem(item)

    def _review_item(self, text: str, rc: int, children: list) -> QTreeWidgetItem:
        first = text.strip().splitlines()[0]
        shown = first + (" …" if "\n" in text.strip() else "")
        item = QTreeWidgetItem([shown, str(rc) if rc else ""])
        item.setData(0, Qt.ItemDataRole.UserRole, text)
        item.setToolTip(0, text)
        if rc:
            _paint(item, QColor(self._scheme().colors["error"]))
        for ch in children:
            item.addChild(self._review_item(ch["text"], ch["rc"], ch["children"]))
        return item

    def _do_label(self, line: str) -> str:
        """`do "<temporary file>"` of the Do-file Editor shown as `do (selection: name.do)`."""
        m = re.match(r'^\s*(\w+)\s+"([^"]+)"\s*$', line)
        if not m or self.dofile_editor is None:
            return ""
        src = self.dofile_editor.tmp_sources.get(m.group(2))
        return f"{m.group(1)} ({src})" if src else ""

    def _startup_session(self) -> None:
        """Initial folder and automatic log (Settings > General > Session)."""
        folder = self.prefs.start_dir
        if folder and os.path.isdir(folder):
            self.session.run(f'quietly cd "{folder}"')
        if self.prefs.auto_log:
            base = Path(folder) if folder and os.path.isdir(folder) else Path.home()
            try:
                (base / "logs").mkdir(exist_ok=True)
            except OSError:
                return
            path = base / "logs" / time.strftime("opendta-%Y-%m-%d-%H%M%S.log")
            self.session.run(f'log using "{path}", text')

    def _confirm(self, line: str) -> bool:
        """Settings > General > Safety: asks before a command that discards much of the data."""
        from .. import guard
        skip = frozenset(c for c in guard.COMMANDS if not self.prefs.confirm_enabled(c))
        if len(skip) == len(guard.COMMANDS):
            return True
        risk = guard.check(self.session, line, skip)
        if risk is None or self._ask_risk(risk, line):
            return True
        self.session.output.write(i18n.L("(not run: cancelled at the confirmation)\n",
                                    "(não executado: cancelado na confirmação)\n"), "text")
        return False

    def _ask_risk(self, risk, line: str) -> bool:
        """The confirmation window. The first time, it also says where to turn it off."""
        first = not self.prefs.confirm_notice_shown
        self.prefs.confirm_notice_shown = True
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Icon.Warning)
        box.setWindowTitle(tr("Confirm command"))
        box.setText(risk.message)
        info = f". {line.strip()}"
        if first:
            info += "\n\n" + tr("You can turn these confirmations off in Settings > General > Safety.")
        box.setInformativeText(info)
        run = box.addButton(tr("Run anyway"), QMessageBox.ButtonRole.AcceptRole)
        cancel = box.addButton(tr("Cancel"), QMessageBox.ButtonRole.RejectRole)
        settings = box.addButton(tr("Cancel and open Settings"), QMessageBox.ButtonRole.ActionRole) if first else None
        box.setDefaultButton(cancel)
        box.setEscapeButton(cancel)
        box.exec()
        if settings is not None and box.clickedButton() is settings:
            self.show_preferences("general", reveal=f"confirm_{risk.command}")
        return box.clickedButton() is run

    def run_command(self, line: str) -> None:
        if not self._confirm(line):
            return
        collect = self.prefs.review_do_commands
        if collect:
            self.session.review_root = {"children": []}
        try:
            rc = self.session.run_command(line)
        except ExitRequest:
            self.close()
            return
        finally:
            root = self.session.review_root
            self.session.review_root = None
            self.session._review_stack.clear()
        self._add_review(line, rc, root["children"] if collect and root else None)
        self.command.setFocus()

    @staticmethod
    def _review_text(item) -> str:
        full = item.data(0, Qt.ItemDataRole.UserRole)
        return full if full else item.text(0)

    def _filter_review(self, text: str) -> None:
        """An item stays visible when it or one of the commands under it matches."""
        t = text.lower()

        def show(it) -> bool:
            kids = [show(it.child(j)) for j in range(it.childCount())]
            ok = not t or t in it.text(0).lower() or any(kids)
            it.setHidden(not ok)
            return ok
        for k in range(self.review.topLevelItemCount()):
            show(self.review.topLevelItem(k))

    # -- data files ---------------------------------------------------
    def _confirm_discard(self) -> bool:
        """Changed data: offers to save before replacing. False = cancel."""
        ds = self.session.data
        if not (ds.changed and ds.nvars):
            return True
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Icon.Warning)
        box.setText(tr("Data in memory have changed."))
        box.setInformativeText(tr("Do you want to save the changes before opening another dataset?"))
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
            path, _ = QFileDialog.getOpenFileName(self, tr("Open"), os.getcwd(),
                                                  tr("Stata data (*.dta);;All files (*)"))
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
        path, _ = QFileDialog.getSaveFileName(self, tr("Save as"), start, tr("Stata data (*.dta)"))
        if not path:
            return False
        if not path.lower().endswith(".dta"):
            path += ".dta"
        # the dialog has already confirmed the replacement, hence replace
        self.run_command(f'save "{path}", replace')
        return self.session.rc == 0

    def show_browser(self, columns=None, rows=None, *, edit: bool = False,
                     nolabel: bool = False) -> None:
        if self.browser is None:
            self.browser = DataBrowser(self.session, monospace_font(self.prefs), self,
                                       run=self._run_from_editor)
            self.browser.setWindowFlag(Qt.WindowType.Window, True)
            self.browser.set_dark(self.prefs.effective_theme == "dark")
        self.browser.model.nolabel = nolabel
        self.browser.model.set_subset(columns, rows)
        self.browser.set_mode(edit)
        _present(self.browser)

    # -- Viewer e log ---------------------------------------------------------------
    def show_viewer(self, path, title: str) -> None:
        if self.viewer is None:
            font = monospace_font(self.prefs)
            self.viewer = Viewer(self.session, self, font_family=font.family(),
                                 font_size=max(9, font.pointSize()), dark=self.prefs.effective_theme == "dark",
                                 run=self._run_from_editor)
            self.viewer.setWindowFlag(Qt.WindowType.Window, True)
        self.viewer.dark = self.prefs.effective_theme == "dark"
        self.viewer.show_page(path, title)
        _present(self.viewer)

    # -- graphs and Do-file Editor ------------------------------------------------
    def show_graph(self, g) -> None:
        if self.graph_window is None:
            self.graph_window = GraphWindow(self.session, self, run=self._run_from_editor)
            self.graph_window.setWindowFlag(Qt.WindowType.Window, True)
        self.graph_window.show_graph(g)

    def close_graph(self, name: str) -> None:
        if self.graph_window is not None:
            self.graph_window.close_graph(name)

    def show_graph_window(self) -> None:
        from ..graphics.store import store
        st = store(self.session)
        if self.graph_window is None or not self.graph_window.canvases:
            if st.current:
                self.show_graph(st.get())
            return
        _present(self.graph_window)

    def show_dofile_editor(self, path: str = "") -> None:
        if self.dofile_editor is None:
            self.dofile_editor = DoFileEditor(self.session, self, run=self.run_command,
                                              font=monospace_font(self.prefs), theme=self.prefs.effective_theme,
                                              palette=self.prefs.editor_palette)
            self.dofile_editor.setWindowFlag(Qt.WindowType.Window, True)
            if self.prefs.reopen_files:
                for f in self.prefs.editor_files:
                    if os.path.isfile(f):
                        self.dofile_editor.open_path(f)
            # connected after the files are restored, so the stored list is not rewritten halfway
            self.dofile_editor.files_changed.connect(
                lambda: setattr(self.prefs, "editor_files", self.dofile_editor.open_paths()))
        if path:
            self.dofile_editor.open_path(path)
        _present(self.dofile_editor)

    def insert_command(self, text: str) -> None:
        self.command.setPlainText(text)
        cur = self.command.textCursor()
        cur.movePosition(cur.MoveOperation.End)
        self.command.setTextCursor(cur)
        self.command.setFocus()

    def copy_table(self) -> None:
        """Copies the Results selection with columns separated by tabs
        (to paste into a spreadsheet), like Stata's Copy table."""
        text = self.results.textCursor().selectedText().replace("\u2029", "\n")
        if not text:
            return
        rows = []
        for line in text.splitlines():
            line = line.replace("|", " ").strip()
            if not line or set(line) <= set("-+= "):
                continue
            rows.append("\t".join(re.split(r"\s{2,}", line)))
        QApplication.clipboard().setText("\n".join(rows))

    def print_results(self) -> None:
        from PySide6.QtPrintSupport import QPrintDialog, QPrinter
        printer = QPrinter()
        dlg = QPrintDialog(printer, self)
        if dlg.exec():
            self.results.print_(printer)

    def find_results(self) -> None:
        from PySide6.QtWidgets import QInputDialog
        text, ok = QInputDialog.getText(self, tr("Find"), tr("Find in Results:"), text=self._last_find)
        if not ok or not text:
            return
        self._last_find = text
        if not self.results.find(text):
            cur = self.results.textCursor()
            cur.movePosition(cur.MoveOperation.Start)
            self.results.setTextCursor(cur)
            if not self.results.find(text):
                self.statusBar().showMessage(f'"{text}" not found', 3000)

    def search_help(self) -> None:
        self.run_command("help")
        if self.viewer is not None:
            self.viewer.field.setFocus()
            self.viewer.field.selectAll()

    def begin_log(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, tr("Begin log"), os.path.join(os.getcwd(), "log.smcl"),
                                              tr("SMCL log (*.smcl);;Text log (*.log)"))
        if path:
            # the dialog has already confirmed the replacement
            opts = "replace text" if path.lower().endswith(".log") else "replace"
            self.run_command(f'log using "{path}", {opts}')

    def toggle_log(self) -> None:
        if getattr(self.session, "logs", None):
            self.run_command("log close")
        else:
            self.begin_log()

    def view_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, tr("View"), os.getcwd(),
                                              tr("Logs and help (*.smcl *.log *.sthlp *.txt);;All files (*)"))
        if path:
            self.run_command(f'view "{path}"')

    def show_variables_manager(self) -> None:
        if self.varmanager is None:
            self.varmanager = VariablesManager(self.session, monospace_font(self.prefs), self,
                                               run=self._run_from_editor)
            self.varmanager.setWindowFlag(Qt.WindowType.Window, True)
        self.varmanager.refresh()
        _present(self.varmanager)

    def _run_from_editor(self, line: str) -> int:
        """Command generated by the Data Editor or the Variables Manager: echoed in
        Results and added to Review, as in Stata."""
        self.run_command(line)
        return self.session.rc

    def choose_do_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, tr("Do"), os.getcwd(),
                                              tr("Do-files (*.do *.ado);;All files (*)"))
        if path:
            self.run_command(f'do "{path}"')

    def choose_directory(self) -> None:
        path = QFileDialog.getExistingDirectory(self, tr("Change working directory"), os.getcwd())
        if path:
            self.run_command(f'cd "{path}"')

    def show_preferences(self, page: str | None = None, reveal: str | None = None) -> None:
        """The Settings window (not modal: the changes show in the windows behind it). `reveal` scrolls
        to a row of the page and flashes it."""
        dlg = self._settings
        if dlg is None or not dlg.isVisible():
            dlg = PreferencesDialog(self.prefs, self)
            dlg.applied.connect(self.apply_preferences)
            dlg.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
            dlg.destroyed.connect(lambda *_: setattr(self, "_settings", None))
            self._settings = dlg
        if page:
            dlg.select(page)
        _present(dlg)
        if reveal:
            QTimer.singleShot(150, lambda: dlg.reveal(reveal))

    def apply_preferences(self) -> None:
        """Applies the chosen theme, font, icon and language (also called at startup)."""
        old = (i18n.language(), i18n.messages())
        i18n.set_language(self.prefs.language)
        i18n.set_messages(self.prefs.messages)
        if (i18n.language(), i18n.messages()) != old:
            self.retranslate()
        self.apply_theme(self.prefs.effective_theme)
        font = monospace_font(self.prefs)
        for w in (self.results, self.command, self.review):
            w.setFont(font)
        self.results.set_solid_rules(self.prefs.solid_table_lines)
        if self.browser is not None:
            self.browser.table.setFont(font)
        if self.dofile_editor is not None:
            self.dofile_editor.set_theme(self.prefs.effective_theme, self.prefs.editor_palette)
        ico = app_icon(self.prefs.app_icon)
        self.setWindowIcon(ico)
        app = QApplication.instance()
        if app is not None:
            app.setWindowIcon(ico)

    # -- layout saved between sessions -------------------------------------------
    def _place_dock(self, dock: QDockWidget) -> None:
        """Review, Variables and Properties always go back to where the application puts them: Review on
        the left, Variables above Properties on the right, docked and at their usual size."""
        dock.setFloating(False)
        if dock is self.dock_review:
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, dock)
            self.resizeDocks([dock], [260], Qt.Orientation.Horizontal)
            return
        right = [d for d in (self.dock_variables, self.dock_properties) if d is dock or d.isVisible()]
        for d in right:
            d.setFloating(False)
            self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, d)
        if len(right) == 2:
            self.splitDockWidget(self.dock_variables, self.dock_properties, Qt.Orientation.Vertical)
            self.resizeDocks(right, [420, 330], Qt.Orientation.Vertical)
        self.resizeDocks([dock], [300], Qt.Orientation.Horizontal)

    def _restore_layout(self) -> None:
        """The size of the window and of the Results/Command split come back; the side windows (Review,
        Variables, Properties) do not: at every start they open at their default place, even if they
        were closed or moved last time."""
        geo, _state, split = self.prefs.layout()
        if geo:
            self.restoreGeometry(geo)
        if split:
            self.splitter.restoreState(split)

    def closeEvent(self, event) -> None:  # noqa: N802
        if self.dofile_editor is not None:
            if not self.dofile_editor.close():
                event.ignore()
                return
            self.dofile_editor.cleanup()
        if self.prefs.auto_log:
            self.session.run("capture quietly log close _all")
        self.prefs.save_layout(self.saveGeometry(), self.saveState(), self.splitter.saveState())
        super().closeEvent(event)
        start = self._pending_update
        if start is not None and event.isAccepted():
            # the installer (or the replacement script) waits for this process to end
            self._pending_update = None
            start()
            QApplication.quit()

    def _scheme(self, theme: str | None = None):
        """Colors of the Results window: the chosen color scheme in the current theme."""
        return results_scheme(theme or self.prefs.effective_theme, self.prefs.results_palette)

    def apply_theme(self, theme: str) -> None:
        """Light or dark theme for the whole interface, independent of the system."""
        app = QApplication.instance()
        if app is not None:
            apply_theme(app, theme)
        set_icon_theme(theme)
        self.results.set_scheme(self._scheme(theme))
        err = QColor(self._scheme(theme).colors["error"])

        def repaint(it) -> None:
            if it.text(1):
                _paint(it, err)
            for j in range(it.childCount()):
                repaint(it.child(j))
        for k in range(self.review.topLevelItemCount()):
            repaint(self.review.topLevelItem(k))
        for act in self.toolbar_actions.values():
            name, accent = act.data()
            act.setIcon(icon(name, accent))
        self._apply_styles()
        if self.browser is not None:
            self.browser.set_dark(theme == "dark")

    def _set_permanently(self, name: str, value: str) -> None:
        """set ..., permanently: remembered at the next startup (hints, superscript, numerics)."""
        if name in ("hints", "superscript"):
            setattr(self.prefs, name, value == "on")
        elif name in ("numerics", "language", "messages"):
            setattr(self.prefs, name, value)

    def _toggle_superscript(self, on: bool) -> None:
        self.run_command(f"set superscript {'on' if on else 'off'}, permanently")

    def _toggle_hints(self, on: bool) -> None:
        # becomes a command, so that it shows in Results and Review like any set
        self.run_command(f"set hints {'on' if on else 'off'}, permanently")

    def show_hints_help(self) -> None:
        QMessageBox.information(self, tr("Error explanations"), i18n.L(
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
            "opens. The menu item <i>Help \u2192 Explain Errors</i> does that.</p>",
            "<p>Quando um comando falha, o OpenDTA mostra a mensagem de erro e o código "
            "<tt>r(#)</tt> de sempre e, entre eles, uma explicação curta do que deu errado "
            "naquele comando e uma dica de como corrigir:</p>"
            "<pre>. replace nota = 9 in 2\nObs. nos. out of range\n"
            "  \u2192 \"in 2\" pede observações que não existem; ...\n"
            "    Dica: crie observações com \"set obs #\" ...\nr(198);</pre>"
            "<p>Os códigos de retorno não mudam: <tt>capture</tt> e <tt>_rc</tt> funcionam "
            "como antes. Erros provocados de propósito com <tt>error #</tt> não são explicados.</p>"
            "<p><b>set hints on</b> &nbsp;liga as explicações (padrão)<br>"
            "<b>set hints off</b> &nbsp;mostra só a mensagem de erro e o código</p>"
            "<p>Acrescente <b>, permanently</b> para guardar a escolha na próxima abertura. "
            "O item <i>Ajuda \u2192 Explicar erros</i> faz isso.</p>"))

    def show_about(self) -> None:
        QMessageBox.about(self, tr("About OpenDTA"),
                          f"<b>OpenDTA {__version__}</b><br>" + i18n.L(
                              "Do-file interpreter with syntax compatible with Stata 14.<br>"
                              "© 2026 OpenDTA contributors — MIT License.<br>"
                              "Stata is a registered trademark of StataCorp LLC;<br>"
                              "OpenDTA is not affiliated with StataCorp.<br>"
                              '<a href="https://github.com/yaskepalhares/opendta">'
                              "github.com/yaskepalhares/opendta</a>",
                              "Interpretador de do-files com sintaxe compatível com o Stata 14.<br>"
                              "© 2026 colaboradores do OpenDTA — licença MIT.<br>"
                              "Stata é marca registrada da StataCorp LLC;<br>"
                              "o OpenDTA não é afiliado à StataCorp.<br>"
                              '<a href="https://github.com/yaskepalhares/opendta">'
                              "github.com/yaskepalhares/opendta</a>"))

    def refresh_state(self) -> None:
        if hasattr(self, "act_hints"):
            self.act_hints.setChecked(self.session.settings.get("hints", "on") == "on")
        if hasattr(self, "act_superscript"):
            self.act_superscript.setChecked(self.session.settings.get("superscript", "on") == "on")
        cwd = os.getcwd()
        self.cwd_label.setText(cwd)
        ds = self.session.data
        title = Path(ds.filename).name if ds.filename else (Path(cwd).name or cwd)
        # window manage maintitle replaces the title until "maintitle reset"
        self.setWindowTitle(getattr(self, "main_title", "") or f"OpenDTA {__version__} — {title}")

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
            child.setText(1, data_props.get(child.data(0, Qt.ItemDataRole.UserRole), ""))
        self._show_variable_properties()
        # the auxiliary windows refresh after the command ends (editing a cell runs commands from
        # inside the model itself)
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
            child.setText(1, props.get(child.data(0, Qt.ItemDataRole.UserRole), ""))

    # -- drag and drop of do-files -------------------------------------------
    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event) -> None:  # noqa: N802
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.lower().endswith((".do", ".ado")):
                if self.prefs.run_dropped_dofiles:
                    self.run_command(f'do "{path}"')
                else:
                    self.show_dofile_editor(path)
            elif path.lower().endswith(".dta"):
                self.open_dataset(path)
                break
