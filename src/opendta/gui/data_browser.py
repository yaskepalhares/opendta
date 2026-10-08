"""Data Editor: grade dos dados em memória, nos modos Browse e Edit.

Como no Data Editor do Stata 14: uma coluna por variável, linhas numeradas
pela observação, valores exibidos no formato da variável. Strings em
vermelho, valores com rótulo em azul, números na cor do texto.

No modo Edit, cada célula alterada vira um comando (replace ... in #,
set obs, generate) executado pela sessão, ecoado em Results e guardado em
Review, de modo que toda alteração fica registrada e reproduzível. A linha
e a coluna vazias depois dos dados criam observações e variáveis novas.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtGui import QAction, QActionGroup, QColor, QKeySequence, QShortcut
from PySide6.QtWidgets import (QLabel, QMainWindow, QMessageBox, QTableView, QToolBar,
                               QVBoxLayout, QWidget)

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
    """Texto entre aspas para um comando; aspas compostas se houver aspas."""
    if '"' in text:
        return '`"' + text + '"\''
    return '"' + text + '"'


def _new_var_name(ds) -> str:
    k = 1
    while ds.has(f"var{k}"):
        k += 1
    return f"var{k}"


def edit_commands(ds, var, obs: int, text: str) -> list[str]:
    """Comandos que gravam `text` na observação obs (0-based) da variável var.

    var None = coluna nova; obs == ds.nobs = observação nova. Levanta
    ValueError com a mensagem para o usuário se o texto não servir."""
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
    # VERIFICAR: mensagem do Stata ao digitar texto numa variável numérica
    raise ValueError(f"'{value}' is not a number and is not a label of {var.name}")


class DatasetModel(QAbstractTableModel):
    def __init__(self, session: "Session", parent=None):
        super().__init__(parent)
        self.session = session
        self.dark = False
        self.columns: list[str] | None = None   # browse varlist
        self.rows = None                         # browse if/in (índices)
        self._cols: list = []
        self._seen = None                        # conjunto de dados exibido
        self._shape = (-1, -1)
        self.editable = False
        self.nolabel = False
        self.run: Callable[[str], int] | None = None   # executa um comando (devolve rc)
        self.on_error: Callable[[str], None] | None = None

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
        """Linha e coluna vazias para criar observações e variáveis (modo Edit)."""
        if not self.editable:
            return 0, 0
        return (1 if self.rows is None else 0, 1 if self.columns is None else 0)

    def refresh(self) -> None:
        ds = self.ds
        same_vars = (ds is self._seen and self.columns is None
                     and [v.name for v in self._cols] == ds.names
                     and all(a is b for a, b in zip(self._cols, ds.vars)))
        shape = (ds.nobs, ds.nvars)
        if same_vars and shape == self._shape and self.rows is None:
            # só valores mudaram: mantém posição de rolagem e seleção
            if self.rowCount() and self.columnCount():
                self.dataChanged.emit(self.index(0, 0),
                                      self.index(self.rowCount() - 1, self.columnCount() - 1))
            return
        self.beginResetModel()
        if ds is not self._seen:                 # outro arquivo: sai do filtro
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
        self.endResetModel()

    def obs(self, row: int) -> int:
        return int(self.rows[row]) if self.rows is not None else row

    def var(self, col: int):
        return self._cols[col]

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        if parent.isValid():
            return 0
        n = len(self.rows) if self.rows is not None else self.ds.nobs
        return n + self._extra()[0]

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._cols) + self._extra()[1]

    def _is_extra(self, index: QModelIndex) -> bool:
        n = len(self.rows) if self.rows is not None else self.ds.nobs
        return index.row() >= n or index.column() >= len(self._cols)

    def flags(self, index: QModelIndex):
        base = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        if self.editable and index.isValid():
            return base | Qt.ItemFlag.ItemIsEditable
        return base

    def setData(self, index: QModelIndex, value, role: int = Qt.ItemDataRole.EditRole) -> bool:  # noqa: N802
        if role != Qt.ItemDataRole.EditRole or not self.editable or self.run is None:
            return False
        var = self._cols[index.column()] if index.column() < len(self._cols) else None
        n = len(self.rows) if self.rows is not None else self.ds.nobs
        obs = self.obs(index.row()) if index.row() < n else self.ds.nobs
        text = str(value)
        if var is not None and obs < self.ds.nobs and text == self.data(index, Qt.ItemDataRole.EditRole):
            return False                          # nada mudou
        if var is None and not text.strip():
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
            if role == Qt.ItemDataRole.BackgroundRole:
                return QColor(127, 127, 127, 28)      # linha/coluna para dados novos
            return "" if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole) else None
        var = self._cols[index.column()]
        i = self.obs(index.row())
        if role == Qt.ItemDataRole.DisplayRole:
            return cell_text(self.ds, var, i, use_labels=not self.nolabel)
        if role == Qt.ItemDataRole.EditRole:
            # valor cru para edição (sem rótulo; strings como estão)
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


class DataBrowser(QMainWindow):
    """Janela separada, como no Stata. Fica aberta e acompanha os dados."""

    def __init__(self, session: "Session", font=None, parent=None,
                 run: Callable[[str], int] | None = None):
        super().__init__(parent)
        self.resize(900, 560)
        self.model = DatasetModel(session, self)
        self.model.run = run
        self.model.on_error = self._error
        self.table = QTableView()
        self.table.setModel(self.model)
        self.table.setAlternatingRowColors(False)
        self.table.verticalHeader().setDefaultSectionSize(20)
        self.table.horizontalHeader().setDefaultSectionSize(96)
        if font is not None:
            self.table.setFont(font)
        self.info = QLabel()
        self.info.setContentsMargins(6, 2, 6, 2)

        tb = QToolBar("Mode")
        tb.setMovable(False)
        group = QActionGroup(self)
        self.act_edit = QAction("Edit", self, checkable=True)
        self.act_browse = QAction("Browse", self, checkable=True)
        for act in (self.act_edit, self.act_browse):
            group.addAction(act)
            tb.addAction(act)
        self.act_edit.triggered.connect(lambda: self.set_mode(True))
        self.act_browse.triggered.connect(lambda: self.set_mode(False))
        self.addToolBar(tb)

        box = QWidget()
        lay = QVBoxLayout(box)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.table, 1)
        lay.addWidget(self.info)
        self.setCentralWidget(box)
        QShortcut(QKeySequence.StandardKey.Close, self, self.close)
        self.table.selectionModel().currentChanged.connect(self._show_cell)
        self.set_mode(False)

    def set_mode(self, edit: bool) -> None:
        self.model.set_editable(edit and self.model.run is not None)
        triggers = (QTableView.EditTrigger.DoubleClicked | QTableView.EditTrigger.EditKeyPressed
                    | QTableView.EditTrigger.AnyKeyPressed) if self.model.editable \
            else QTableView.EditTrigger.NoEditTriggers
        self.table.setEditTriggers(triggers)
        self.act_edit.setChecked(self.model.editable)
        self.act_browse.setChecked(not self.model.editable)
        self.setWindowTitle("Data Editor (Edit)" if self.model.editable else "Data Editor (Browse)")
        self.refresh()

    def _error(self, message: str) -> None:
        QMessageBox.warning(self, "Data Editor", message)

    def set_dark(self, dark: bool) -> None:
        self.model.dark = dark
        self.model._shape = (-1, -1)
        self.model.refresh()

    def refresh(self) -> None:
        self.model.refresh()
        ds = self.model.ds
        filt = "On" if self.model.rows is not None else "Off"
        mode = "Edit" if self.model.editable else "Browse"
        nvars = len(self.model._cols)
        self.info.setText(f"Vars: {nvars:,}    Obs: {ds.nobs:,}    Filter: {filt}    Mode: {mode}")

    def _show_cell(self, current: QModelIndex, _prev=None) -> None:
        if not current.isValid() or self.model._is_extra(current):
            self.statusBar().clearMessage()
            return
        ds = self.model.ds
        var = self.model.var(current.column())
        i = self.model.obs(current.row())
        raw = cell_text(ds, var, i, use_labels=False)
        self.statusBar().showMessage(f"[{i + 1},{ds.index(var.name) + 1}]  {var.name} = {raw}")
