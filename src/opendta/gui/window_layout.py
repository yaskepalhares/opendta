"""How the main window is built: the central splitter with the Results and Command panes, the
docks (Review, Variables, Properties), the menus, the toolbar and the status bar, and the
retranslation of their texts when the interface language changes. MainWindow (main_window.py)
inherits these methods and keeps the behavior: running commands, opening windows and files."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QFont, QKeySequence
from PySide6.QtWidgets import (QDockWidget, QHeaderView, QLabel, QLineEdit, QMenu, QSizePolicy,
                               QSplitter, QToolBar, QToolButton, QTreeWidget, QTreeWidgetItem,
                               QVBoxLayout, QWidget)
from shiboken6 import isValid

from .. import i18n
from . import qt_i18n
from ..i18n import tr
from .icons import ACCENT_RED, TOOLBAR_SIZE, icon
from .menu_templates import DATA_TEMPLATES, GRAPHICS_TEMPLATES, STATISTICS_TEMPLATES
from .panes import CommandLine, ResultsView, monospace_font, titled
from .window_cmds import build_user_menu

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


class WindowLayout:
    """The parts of MainWindow that build its widgets (a mixin of MainWindow)."""

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
        build_user_menu(self.user_menu, getattr(self, "_user_items", None) or [], self.run_command)

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
        self._action(m, "Contents", self.help_contents)
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
                btn.setFixedWidth(max(TOOLBAR_SIZE.width() + 16, btn.fontMetrics().horizontalAdvance(act.text()) + 14))

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
