"""Data Editor: grid of the data in memory, in Browse and Edit modes.

As in Stata 14's Data Editor: one column per variable, rows numbered
by observation, values shown in the variable's format. Strings in
red, labeled values in blue, numbers in the text color.

In Edit mode, each changed cell becomes a command (replace ... in #,
set obs, generate) run by the session, echoed in Results and kept in
Review, so that every change is recorded and reproducible. The empty row
and column after the data create new observations and variables.
"""

from __future__ import annotations

from ..i18n import L, tr
from .icons import icon
from .style import apply_rounded

from typing import TYPE_CHECKING, Callable

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSize, Qt, Signal
from PySide6.QtGui import QAction, QActionGroup, QColor, QFont, QKeySequence, QShortcut
from PySide6.QtWidgets import (QDialog, QHBoxLayout, QHeaderView, QInputDialog, QLabel, QLineEdit, QMainWindow,
                               QMessageBox, QPushButton, QSizePolicy, QSplitter, QTableView, QToolBar,
                               QToolButton, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget)

from ..commands.inspect import cell_text
from ..core.formats import parse_format
from ..io.fixed import to_number

if TYPE_CHECKING:
    from ..session import Session

STRING_COLOR = QColor("#C0392B")
LABEL_COLOR = QColor("#1F5FBF")
STRING_COLOR_DARK = QColor("#FF7B72")
LABEL_COLOR_DARK = QColor("#79B8FF")


def _quote(text: str) -> str:
    """Quoted text for a command; compound quotes if there are quotes."""
    if '"' in text:
        return '`"' + text + '"\''
    return '"' + text + '"'


def _new_var_name(ds) -> str:
    k = 1
    while ds.has(f"var{k}"):
        k += 1
    return f"var{k}"


def edit_commands(ds, var, obs: int, text: str) -> list[str]:
    """Commands that write `text` into observation obs (0-based) of variable var.

    var None = new column; obs == ds.nobs = new observation. Raises
    ValueError with the message for the user if the text does not fit."""
    cmds: list[str] = []
    if obs >= ds.nobs:
        cmds.append(f"set obs {obs + 1}")
    n = obs + 1
    text = text.strip() if var is None or not var.is_string else text
    if var is None:
        name = _new_var_name(ds)
        x = to_number(text)
        if x is not None:
            cmds.append(f"generate {name} = {text or '.'} in {n}")
        else:
            cmds.append(f"generate {name} = {_quote(text)} in {n}")
        return cmds
    if var.is_string:
        cmds.append(f"replace {var.name} = {_quote(text)} in {n}")
        return cmds
    value = text.strip()
    if to_number(value) is not None:
        cmds.append(f"replace {var.name} = {value or '.'} in {n}")
        return cmds
    labels = ds.value_labels.get(var.value_label, {}) if var.value_label else {}
    for code, lab in labels.items():
        if lab == value:
            cmds.append(f"replace {var.name} = {code} in {n}")
            return cmds
    # VERIFY: Stata's message when text is typed into a numeric variable
    raise ValueError(f"'{value}' is not a number and is not a label of {var.name}")


class DatasetModel(QAbstractTableModel):
    def __init__(self, session: "Session", parent=None):
        super().__init__(parent)
        self.session = session
        self.dark = False
        self.columns: list[str] | None = None   # browse varlist
        self.rows = None                         # browse if/in (indices)
        self._cols: list = []
        self._seen = None                        # displayed dataset
        self._shape = (-1, -1)
        self.editable = False
        self.nolabel = False
        self.run: Callable[[str], int] | None = None   # runs a command (returns rc)
        self.on_error: Callable[[str], None] | None = None
        self.fill = (0, 0)                       # rows and columns the view can show
        self._pad = (0, 0)                       # empty rows/columns added to fill the view

    @property
    def ds(self):
        return self.session.data

    def set_subset(self, columns, rows) -> None:
        self._seen = self.ds
        self.columns, self.rows = columns, rows
        self._shape = (-1, -1)
        self.refresh()

    def set_editable(self, on: bool) -> None:
        self.editable = on
        self._shape = (-1, -1)
        self.refresh()

    def _extra(self) -> tuple[int, int]:
        """Empty row and column to create observations and variables (Edit mode)."""
        if not self.editable:
            return 0, 0
        return (1 if self.rows is None else 0, 1 if self.columns is None else 0)

    def _padding(self) -> tuple[int, int]:
        """Empty rows and columns that fill the view while there is no data (like an empty sheet)."""
        ds = self.ds
        if ds.nobs or self.rows is not None or self.columns is not None:
            return 0, 0
        er, ec = self._extra()
        rows = max(0, self.fill[0] - er)
        cols = max(0, self.fill[1] - len(self._cols) - ec) if not ds.nvars else 0
        return rows, cols

    def set_fill(self, rows: int, cols: int) -> None:
        """Size of the view in cells; empty cells are added up to it when the data is empty."""
        self.fill = (rows, cols)
        new = self._padding()
        old = self._pad
        if new == old:
            return
        base_r = self.rowCount() - old[0]
        base_c = self.columnCount() - old[1]
        if new[0] > old[0]:
            self.beginInsertRows(QModelIndex(), base_r + old[0], base_r + new[0] - 1)
            self._pad = (new[0], self._pad[1])
            self.endInsertRows()
        elif new[0] < old[0]:
            self.beginRemoveRows(QModelIndex(), base_r + new[0], base_r + old[0] - 1)
            self._pad = (new[0], self._pad[1])
            self.endRemoveRows()
        if new[1] > old[1]:
            self.beginInsertColumns(QModelIndex(), base_c + old[1], base_c + new[1] - 1)
            self._pad = (self._pad[0], new[1])
            self.endInsertColumns()
        elif new[1] < old[1]:
            self.beginRemoveColumns(QModelIndex(), base_c + new[1], base_c + old[1] - 1)
            self._pad = (self._pad[0], new[1])
            self.endRemoveColumns()

    def refresh(self) -> None:
        ds = self.ds
        same_vars = (ds is self._seen and self.columns is None
                     and [v.name for v in self._cols] == ds.names
                     and all(a is b for a, b in zip(self._cols, ds.vars)))
        shape = (ds.nobs, ds.nvars)
        if same_vars and shape == self._shape and self.rows is None:
            # only values changed: keeps the scroll position and the selection
            if self.rowCount() and self.columnCount():
                self.dataChanged.emit(self.index(0, 0),
                                      self.index(self.rowCount() - 1, self.columnCount() - 1))
            return
        self.beginResetModel()
        if ds is not self._seen:                 # another file: leaves the filter
            self.columns = self.rows = None
            self._seen = ds
        if self.columns is not None and all(ds.has(n) for n in self.columns):
            self._cols = [ds.get(n) for n in self.columns]
        else:
            self.columns = None
            self._cols = list(ds.vars)
        if self.rows is not None and (len(self.rows) and self.rows.max() >= ds.nobs):
            self.rows = None
        self._shape = (ds.nobs, ds.nvars)
        self._pad = self._padding()
        self.endResetModel()

    def obs(self, row: int) -> int:
        return int(self.rows[row]) if self.rows is not None else row

    def var(self, col: int):
        return self._cols[col]

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        if parent.isValid():
            return 0
        n = len(self.rows) if self.rows is not None else self.ds.nobs
        return n + self._extra()[0] + self._pad[0]

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._cols) + self._extra()[1] + self._pad[1]

    def _is_extra(self, index: QModelIndex) -> bool:
        n = len(self.rows) if self.rows is not None else self.ds.nobs
        return index.row() >= n or index.column() >= len(self._cols)

    def flags(self, index: QModelIndex):
        base = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        # only the first empty column creates a variable; empty rows below the data create observations
        if self.editable and index.isValid() and index.column() <= len(self._cols):
            return base | Qt.ItemFlag.ItemIsEditable
        return base

    def setData(self, index: QModelIndex, value, role: int = Qt.ItemDataRole.EditRole) -> bool:  # noqa: N802
        if role != Qt.ItemDataRole.EditRole or not self.editable or self.run is None:
            return False
        var = self._cols[index.column()] if index.column() < len(self._cols) else None
        n = len(self.rows) if self.rows is not None else self.ds.nobs
        obs = self.obs(index.row()) if index.row() < n else (index.row() if self.rows is None else self.ds.nobs)
        text = str(value)
        if var is not None and obs < self.ds.nobs and text == self.data(index, Qt.ItemDataRole.EditRole):
            return False                          # nada mudou
        if var is None and (not text.strip() or index.column() > len(self._cols)):
            return False
        try:
            cmds = edit_commands(self.ds, var, obs, text)
        except ValueError as e:
            if self.on_error:
                self.on_error(str(e))
            return False
        for cmd in cmds:
            if self.run(cmd):
                return False
        return True

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        if self._is_extra(index):
            if role == Qt.ItemDataRole.BackgroundRole and self._extra() != (0, 0) and not self._pad[0] \
                    and not self._pad[1]:
                return QColor(127, 127, 127, 28)      # row/column for new data
            return "" if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole) else None
        var = self._cols[index.column()]
        i = self.obs(index.row())
        if role == Qt.ItemDataRole.DisplayRole:
            return cell_text(self.ds, var, i, use_labels=not self.nolabel)
        if role == Qt.ItemDataRole.EditRole:
            # raw value for editing (no label; strings as they are)
            return var.value(i) if var.is_string else cell_text(self.ds, var, i, use_labels=False)
        if role == Qt.ItemDataRole.TextAlignmentRole:
            left = var.is_string and parse_format(var.fmt).left
            h = Qt.AlignmentFlag.AlignLeft if left else Qt.AlignmentFlag.AlignRight
            return int(h | Qt.AlignmentFlag.AlignVCenter)
        if role == Qt.ItemDataRole.ForegroundRole:
            if self.nolabel and not var.is_string:
                return None
            if var.is_string:
                return STRING_COLOR_DARK if self.dark else STRING_COLOR
            if var.value_label and var.value_label in self.ds.value_labels:
                text = cell_text(self.ds, var, i)
                raw = cell_text(self.ds, var, i, use_labels=False)
                if text != raw:
                    return LABEL_COLOR_DARK if self.dark else LABEL_COLOR
        if role == Qt.ItemDataRole.ToolTipRole and var.value_label:
            return cell_text(self.ds, var, i, use_labels=False)
        return None

    def headerData(self, section: int, orientation, role: int = Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if role != Qt.ItemDataRole.DisplayRole:
            if (role == Qt.ItemDataRole.ToolTipRole and orientation == Qt.Orientation.Horizontal
                    and section < len(self._cols)):
                return self._cols[section].label or None
            return None
        if orientation == Qt.Orientation.Horizontal:
            return self._cols[section].name if section < len(self._cols) else ""
        n = len(self.rows) if self.rows is not None else self.ds.nobs
        return str(self.obs(section) + 1) if section < n else ""




ROW_HEIGHT = 20
COL_WIDTH = 96

_EDITOR_QSS = """
QToolBar#EditorToolbar { border: none; spacing: 2px; padding: 4px 10px; }
QToolBar#EditorToolbar QToolButton { border: none; border-radius: 8px; padding: 4px 4px; }
QToolBar#EditorToolbar QToolButton:hover { background: rgba(127, 127, 127, 0.16); }
QToolBar#EditorToolbar QToolButton:pressed { background: rgba(127, 127, 127, 0.28); }
QToolBar#EditorToolbar QToolButton:checked { background: rgba(10, 132, 255, 0.16); color: #0a84ff; }
QWidget#CellBar QLineEdit { border: 1px solid rgba(127, 127, 127, 0.4); border-radius: 8px; padding: 4px 8px; }
QWidget#CellBar QLineEdit:focus { border: 1px solid #0a84ff; }
QTableView#Grid { border: none; gridline-color: rgba(127, 127, 127, 0.28); }
QTableView#Grid { selection-background-color: rgba(10, 132, 255, 0.22); selection-color: palette(text); }
QTableView#Grid QHeaderView::section { background: rgba(127, 127, 127, 0.10); border: none;
    border-right: 1px solid rgba(127, 127, 127, 0.25); border-bottom: 1px solid rgba(127, 127, 127, 0.25);
    padding: 0 6px; }
QTableView#Grid QHeaderView::section:checked { background: rgba(10, 132, 255, 0.22); }
QLabel#PaneTitle { font-weight: 600; padding: 6px 10px; }
QTreeWidget#SidePane { border: none; background: transparent; outline: none; }
QTreeWidget#SidePane::item { padding: 2px 0; border-radius: 6px; }
QTreeWidget#SidePane::item:selected { background: rgba(10, 132, 255, 0.22); color: palette(text); }
QSplitter::handle { background: rgba(127, 127, 127, 0.18); }
QStatusBar { border-top: 1px solid rgba(127, 127, 127, 0.25); }
QStatusBar QLabel { padding: 0 8px; }
"""


class _Grid(QTableView):
    resized = Signal()

    def resizeEvent(self, e) -> None:  # noqa: N802
        super().resizeEvent(e)
        self.resized.emit()


class SnapshotsDialog(QDialog):
    """List, take, restore and erase snapshots (the snapshot command), as a small modeless window."""

    def __init__(self, session: "Session", run: Callable[[str], int], parent=None):
        super().__init__(parent)
        self.session, self.run = session, run
        self.setWindowTitle(L("Snapshots", "Instantâneos"))
        self.resize(460, 340)
        apply_rounded(self)
        lay = QVBoxLayout(self)
        self.tree = QTreeWidget()
        self.tree.setRootIsDecorated(False)
        self.tree.setHeaderLabels(["#", L("Label", "Rótulo"), L("Created", "Criado em")])
        self.tree.header().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.tree.itemSelectionChanged.connect(self._sync)
        lay.addWidget(self.tree, 1)
        row = QHBoxLayout()
        self.label = QLineEdit()
        self.label.setPlaceholderText(L("Label of the new snapshot", "Rótulo do novo instantâneo"))
        self.take = QPushButton(L("Take snapshot", "Criar instantâneo"))
        self.take.clicked.connect(self._take)
        self.label.returnPressed.connect(self._take)
        row.addWidget(self.label, 1)
        row.addWidget(self.take)
        lay.addLayout(row)
        row = QHBoxLayout()
        self.restore = QPushButton(L("Restore", "Restaurar"))
        self.erase = QPushButton(L("Erase", "Excluir"))
        self.restore.clicked.connect(lambda: self._on_selected("restore"))
        self.erase.clicked.connect(lambda: self._on_selected("erase"))
        row.addWidget(self.restore)
        row.addWidget(self.erase)
        row.addStretch(1)
        close = QPushButton(L("Close", "Fechar"))
        close.clicked.connect(self.close)
        row.addWidget(close)
        lay.addLayout(row)
        self.reload()

    def reload(self) -> None:
        self.tree.clear()
        for k, snap in sorted(getattr(self.session, "_snapshots", {}).items()):
            it = QTreeWidgetItem([str(k), snap["label"], snap["when"]])
            it.setData(0, Qt.ItemDataRole.UserRole, k)
            self.tree.addTopLevelItem(it)
        self._sync()

    def _sync(self) -> None:
        on = bool(self.tree.selectedItems())
        self.restore.setEnabled(on)
        self.erase.setEnabled(on)

    def _take(self) -> None:
        text = self.label.text().replace('"', "'").strip()
        self.run(f'snapshot save, label("{text}")' if text else "snapshot save")
        self.label.clear()
        self.reload()

    def _on_selected(self, verb: str) -> None:
        items = self.tree.selectedItems()
        if items:
            self.run(f"snapshot {verb} {items[0].data(0, Qt.ItemDataRole.UserRole)}")
            self.reload()


class DataBrowser(QMainWindow):
    """Separate window, as in Stata: toolbar, cell bar, grid, Variables and Properties panes and a
    status bar. Stays open and follows the data."""

    def __init__(self, session: "Session", font=None, parent=None,
                 run: Callable[[str], int] | None = None):
        super().__init__(parent)
        self.resize(1100, 640)
        self.session = session
        self._run = run
        self._busy = False
        self._snapshots: SnapshotsDialog | None = None
        self.model = DatasetModel(session, self)
        self.model.run = run
        self.model.on_error = self._error
        self.table = _Grid()
        self.table.setObjectName("Grid")
        self.table.setModel(self.model)
        self.table.setAlternatingRowColors(False)
        self.table.verticalHeader().setDefaultSectionSize(ROW_HEIGHT)
        self.table.horizontalHeader().setDefaultSectionSize(COL_WIDTH)
        self.table.horizontalHeader().setHighlightSections(True)
        self.table.verticalHeader().setHighlightSections(True)
        self.table.resized.connect(self._fit_fill)
        if font is not None:
            self.table.setFont(font)

        # -- toolbar: Edit and Browse on the left; Filter, Variables, Properties and Snapshots on the right
        tb = QToolBar("Data Editor")
        tb.setObjectName("EditorToolbar")
        tb.setMovable(False)
        tb.setIconSize(QSize(22, 22))
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        group = QActionGroup(self)
        self.act_edit = QAction(icon("dataeditor"), tr("Edit"), self, checkable=True)
        self.act_browse = QAction(icon("databrowser"), tr("Browse"), self, checkable=True)
        for act in (self.act_edit, self.act_browse):
            group.addAction(act)
            tb.addAction(act)
        self.act_edit.triggered.connect(lambda: self._switch(True))
        self.act_browse.triggered.connect(lambda: self._switch(False))
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        tb.addWidget(spacer)
        self.act_filter = QAction(icon("filter"), L("Filter", "Filtro"), self, checkable=True)
        self.act_vars = QAction(icon("variables"), tr("Variables"), self, checkable=True)
        self.act_props = QAction(icon("properties"), L("Properties", "Propriedades"), self, checkable=True)
        self.act_snaps = QAction(icon("snapshots"), L("Snapshots", "Instantâneos"), self)
        for act in (self.act_filter, self.act_vars, self.act_props, self.act_snaps):
            tb.addAction(act)
        self.act_vars.setChecked(True)
        self.act_props.setChecked(True)
        self.act_filter.triggered.connect(self._toggle_filter)
        self.act_vars.toggled.connect(lambda on: self.vars_pane.setVisible(on))
        self.act_props.toggled.connect(lambda on: self.props_pane.setVisible(on))
        self.act_snaps.triggered.connect(self._show_snapshots)
        self.addToolBar(tb)
        self.toolbar = tb
        self._fit_toolbar()

        # -- cell bar: name[obs] and the value of the current cell
        bar = QWidget()
        bar.setObjectName("CellBar")
        blay = QHBoxLayout(bar)
        blay.setContentsMargins(10, 8, 10, 8)
        blay.setSpacing(8)
        self.cell_name = QLineEdit()
        self.cell_name.setReadOnly(True)
        self.cell_name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.cell_name.setFixedWidth(150)
        self.cell_name.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.cell_value = QLineEdit()
        self.cell_value.returnPressed.connect(self._commit_cell)
        blay.addWidget(self.cell_name)
        blay.addWidget(self.cell_value, 1)

        # -- side panes
        self.vars_list = QTreeWidget()
        self.vars_list.setObjectName("SidePane")
        self.vars_list.setRootIsDecorated(False)
        self.vars_list.setHeaderLabels([tr("Name"), tr("Label")])
        self.vars_list.itemClicked.connect(self._go_to_variable)
        self.vars_search = QLineEdit()
        self.vars_search.setPlaceholderText(L("Filter variables", "Filtrar variáveis"))
        self.vars_search.textChanged.connect(self._filter_vars)
        self.vars_pane = self._pane(tr("Variables"), self.vars_list, self.vars_search)
        self.props = QTreeWidget()
        self.props.setObjectName("SidePane")
        self.props.setHeaderHidden(True)
        self.props.setColumnCount(2)
        self.props.header().setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        self.props.itemChanged.connect(self._property_edited)
        self.props_pane = self._pane(L("Properties", "Propriedades"), self.props)
        side = QSplitter(Qt.Orientation.Vertical)
        side.addWidget(self.vars_pane)
        side.addWidget(self.props_pane)
        side.setSizes([320, 280])
        left = QWidget()
        llay = QVBoxLayout(left)
        llay.setContentsMargins(0, 0, 0, 0)
        llay.setSpacing(0)
        llay.addWidget(bar)
        llay.addWidget(self.table, 1)
        main = QSplitter(Qt.Orientation.Horizontal)
        main.addWidget(left)
        main.addWidget(side)
        main.setStretchFactor(0, 1)
        main.setSizes([820, 280])
        self.setCentralWidget(main)

        # -- status bar: Vars / Order / Obs / Length / Filter
        self.st = {k: QLabel() for k in ("vars", "obs", "length", "filter")}
        sb = self.statusBar()
        sb.setSizeGripEnabled(False)
        sb.addWidget(self.st["vars"])
        sb.addWidget(self.st["obs"])
        sb.addPermanentWidget(self.st["length"])
        sb.addPermanentWidget(self.st["filter"])

        self.setStyleSheet(_EDITOR_QSS)
        apply_rounded(self)
        QShortcut(QKeySequence.StandardKey.Close, self, self.close)
        self.table.selectionModel().currentChanged.connect(self._show_cell)
        self.set_mode(False)

    # -- construction helpers
    def _pane(self, title: str, tree: QTreeWidget, search: QLineEdit | None = None) -> QWidget:
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        head = QLabel(title)
        head.setObjectName("PaneTitle")
        lay.addWidget(head)
        lay.addWidget(tree, 1)
        if search is not None:
            box = QWidget()
            blay = QHBoxLayout(box)
            blay.setContentsMargins(8, 6, 8, 8)
            blay.addWidget(search)
            lay.addWidget(box)
        return w

    def _fit_toolbar(self) -> None:
        font = QFont(self.font())
        font.setPixelSize(11)
        for act in self.toolbar.actions():
            btn = self.toolbar.widgetForAction(act)
            if isinstance(btn, QToolButton):
                btn.setFont(font)
                btn.setFixedWidth(max(48, btn.fontMetrics().horizontalAdvance(act.text()) + 20))

    def _fit_fill(self) -> None:
        vp = self.table.viewport()
        self.model.set_fill(vp.height() // ROW_HEIGHT + 1, vp.width() // COL_WIDTH + 1)
        # the empty grid fills the view exactly: no scroll bars until there is data
        policy = Qt.ScrollBarPolicy.ScrollBarAlwaysOff if self.model._pad != (0, 0) \
            else Qt.ScrollBarPolicy.ScrollBarAsNeeded
        self.table.setHorizontalScrollBarPolicy(policy)
        self.table.setVerticalScrollBarPolicy(policy)

    # -- modes
    def _switch(self, edit: bool) -> None:
        self.set_mode(edit)

    def set_mode(self, edit: bool) -> None:
        self.model.set_editable(edit and self.model.run is not None)
        triggers = (QTableView.EditTrigger.DoubleClicked | QTableView.EditTrigger.EditKeyPressed
                    | QTableView.EditTrigger.AnyKeyPressed) if self.model.editable \
            else QTableView.EditTrigger.NoEditTriggers
        self.table.setEditTriggers(triggers)
        self.act_edit.setChecked(self.model.editable)
        self.act_browse.setChecked(not self.model.editable)
        self.cell_value.setReadOnly(not self.model.editable)
        self.setWindowTitle(tr("Data Editor (Edit)") if self.model.editable else tr("Data Editor (Browse)"))
        self._fit_fill()
        self.refresh()

    def _error(self, message: str) -> None:
        QMessageBox.warning(self, tr("Data Editor"), message)

    def set_dark(self, dark: bool) -> None:
        self.model.dark = dark
        self.model._shape = (-1, -1)
        self.model.refresh()
        for act, name in ((self.act_edit, "dataeditor"), (self.act_browse, "databrowser"),
                          (self.act_filter, "filter"), (self.act_vars, "variables"),
                          (self.act_props, "properties"), (self.act_snaps, "snapshots")):
            act.setIcon(icon(name))

    # -- filter and snapshots
    def _toggle_filter(self, on: bool) -> None:
        cmd = "edit" if self.model.editable else "browse"
        if not on:
            self._run_cmd(cmd)
            return
        expr, ok = QInputDialog.getText(self, L("Filter", "Filtro"),
                                        L("Show only the observations where (an expression):",
                                          "Mostrar só as observações em que (uma expressão):"))
        if not ok or not expr.strip():
            self.act_filter.setChecked(self.model.rows is not None)
            return
        self._run_cmd(f"{cmd} if {expr.strip()}")
        self.act_filter.setChecked(self.model.rows is not None)

    def _run_cmd(self, cmd: str) -> int:
        return self._run(cmd) if self._run is not None else 1

    def _show_snapshots(self) -> None:
        if self._run is None:
            return
        if self._snapshots is None:
            self._snapshots = SnapshotsDialog(self.session, self._run, self)
            self._snapshots.setWindowFlag(Qt.WindowType.Window, True)
        self._snapshots.reload()
        self._snapshots.show()
        self._snapshots.raise_()

    # -- cell bar
    def _cell_label(self, current: QModelIndex) -> tuple[str, str]:
        """(name[obs], text) of the cell; empty cells name the variable they would create."""
        if not current.isValid():
            return "", ""
        if self.model._is_extra(current):
            return f"var{current.column() + 1}[{current.row() + 1}]", ""
        var = self.model.var(current.column())
        i = self.model.obs(current.row())
        return f"{var.name}[{i + 1}]", cell_text(self.model.ds, var, i, use_labels=False)

    def _show_cell(self, current: QModelIndex, _prev=None) -> None:
        name, text = self._cell_label(current)
        self.cell_name.setText(name)
        self.cell_value.setText(text)
        self.st["length"].setText(L(f"Length: {len(text)}", f"Tamanho: {len(text)}"))
        self._show_properties()

    def _commit_cell(self) -> None:
        cur = self.table.currentIndex()
        if cur.isValid() and self.model.editable:
            self.model.setData(cur, self.cell_value.text())
            self.table.setFocus()

    # -- Variables pane
    def _fill_variables(self) -> None:
        cols = self.model._cols
        sig = [(v.name, v.label) for v in cols]
        if sig == getattr(self, "_vars_sig", None):
            return
        self._vars_sig = sig
        self.vars_list.clear()
        for k, (name, label) in enumerate(sig):
            it = QTreeWidgetItem([name, label])
            it.setData(0, Qt.ItemDataRole.UserRole, k)
            self.vars_list.addTopLevelItem(it)
        self._filter_vars(self.vars_search.text())

    def _filter_vars(self, text: str) -> None:
        text = text.strip().lower()
        for k in range(self.vars_list.topLevelItemCount()):
            it = self.vars_list.topLevelItem(k)
            it.setHidden(bool(text) and text not in it.text(0).lower() and text not in it.text(1).lower())

    def _go_to_variable(self, item: QTreeWidgetItem) -> None:
        col = item.data(0, Qt.ItemDataRole.UserRole)
        row = max(0, self.table.currentIndex().row())
        idx = self.model.index(min(row, max(0, self.model.rowCount() - 1)), col)
        self.table.setCurrentIndex(idx)
        self.table.scrollTo(idx)

    # -- Properties pane
    def _show_properties(self) -> None:
        self._busy = True
        try:
            self.props.clear()
            ds = self.model.ds
            cur = self.table.currentIndex()
            var = self.model.var(cur.column()) if cur.isValid() and cur.column() < len(self.model._cols) else None
            editable = self.model.editable

            def group(title: str, rows: list[tuple[str, str, bool]], key: str) -> None:
                top = QTreeWidgetItem([title, ""])
                top.setFlags(Qt.ItemFlag.ItemIsEnabled)
                self.props.addTopLevelItem(top)
                for label, value, can_edit in rows:
                    child = QTreeWidgetItem([label, value])
                    flags = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
                    if can_edit and editable:
                        flags |= Qt.ItemFlag.ItemIsEditable
                    child.setFlags(flags)
                    child.setData(0, Qt.ItemDataRole.UserRole, (key, label))
                    top.addChild(child)
                top.setExpanded(True)

            if var is not None:
                group(tr("Variables"), [(tr("Name"), var.name, True), (tr("Label"), var.label, True),
                                        (tr("Type"), var.vtype, False), (tr("Format"), var.fmt, True),
                                        (tr("Value label"), var.value_label or "", False)], "var")
            else:
                group(tr("Variables"), [], "var")
            order = ", ".join(ds.sortlist) if ds.sortlist else L("(unsorted)", "(sem ordem)")
            group(tr("Data"), [(tr("Label"), ds.label, True), (tr("Filename"), ds.filename, False),
                               (L("Variables", "Variáveis"), f"{ds.nvars:,}", False),
                               (tr("Observations"), f"{ds.nobs:,}", False), (tr("Sorted by"), order, False)], "data")
        finally:
            self._busy = False

    def _property_edited(self, item: QTreeWidgetItem, column: int) -> None:
        key = item.data(0, Qt.ItemDataRole.UserRole)
        if self._busy or column != 1 or not key or self._run is None:
            return
        kind, label = key
        text = item.text(1).strip().replace('"', "'")
        cur = self.table.currentIndex()
        var = self.model.var(cur.column()) if cur.isValid() and cur.column() < len(self.model._cols) else None
        if kind == "data":
            self._run_cmd(f'label data "{text}"')
        elif var is not None and label == tr("Name") and text:
            self._run_cmd(f"rename {var.name} {text}")
        elif var is not None and label == tr("Label"):
            self._run_cmd(f'label variable {var.name} "{text}"')
        elif var is not None and label == tr("Format") and text:
            self._run_cmd(f"format {var.name} {text}")
        self.refresh()

    # -- refresh
    def refresh(self) -> None:
        self.model.refresh()
        self._fit_fill()
        ds = self.model.ds
        nvars = len(self.model._cols)
        filtered = self.model.rows is not None
        self.act_filter.setEnabled(ds.nvars > 0)
        self.act_filter.setChecked(filtered)
        self.st["vars"].setText(L(f"Vars: {nvars:,}    Order: Dataset", f"Vars: {nvars:,}    Ordem: Dados"))
        self.st["obs"].setText(L(f"Obs: {ds.nobs:,}", f"Obs: {ds.nobs:,}"))
        self.st["filter"].setText(L("Filter: ", "Filtro: ") + (L("On", "Ligado") if filtered else L("Off", "Desligado")))
        if not self.table.currentIndex().isValid() and self.model.rowCount() and self.model.columnCount():
            self.table.setCurrentIndex(self.model.index(0, 0))
        self._fill_variables()
        self._show_cell(self.table.currentIndex())
