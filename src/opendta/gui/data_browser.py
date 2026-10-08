"""Data Browser: grade somente leitura dos dados em memória.

Como no Data Editor (Browse) do Stata 14: uma coluna por variável, linhas
numeradas pela observação, valores exibidos no formato da variável. Strings
em vermelho, valores com rótulo em azul, números na cor do texto. A edição
de células fica para a fase de ferramentas da interface.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import QAbstractTableModel, QModelIndex, Qt
from PySide6.QtGui import QColor, QKeySequence, QShortcut
from PySide6.QtWidgets import QLabel, QMainWindow, QTableView, QVBoxLayout, QWidget

from ..commands.inspect import cell_text
from ..core.formats import parse_format

if TYPE_CHECKING:
    from ..session import Session

STRING_COLOR = QColor("#C0392B")
LABEL_COLOR = QColor("#1F5FBF")
STRING_COLOR_DARK = QColor("#FF7B72")
LABEL_COLOR_DARK = QColor("#79B8FF")


class DatasetModel(QAbstractTableModel):
    def __init__(self, session: "Session", parent=None):
        super().__init__(parent)
        self.session = session
        self.dark = False
        self.columns: list[str] | None = None   # browse varlist
        self.rows = None                         # browse if/in (índices)
        self._cols: list = []
        self._seen = None                        # conjunto de dados exibido

    @property
    def ds(self):
        return self.session.data

    def set_subset(self, columns, rows) -> None:
        self._seen = self.ds
        self.columns, self.rows = columns, rows
        self.refresh()

    def refresh(self) -> None:
        self.beginResetModel()
        ds = self.ds
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
        self.endResetModel()

    def obs(self, row: int) -> int:
        return int(self.rows[row]) if self.rows is not None else row

    def var(self, col: int):
        return self._cols[col]

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        if parent.isValid():
            return 0
        return len(self.rows) if self.rows is not None else self.ds.nobs

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._cols)

    def data(self, index: QModelIndex, role: int = Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        var = self._cols[index.column()]
        i = self.obs(index.row())
        if role == Qt.ItemDataRole.DisplayRole:
            return cell_text(self.ds, var, i)
        if role == Qt.ItemDataRole.TextAlignmentRole:
            left = var.is_string and parse_format(var.fmt).left
            h = Qt.AlignmentFlag.AlignLeft if left else Qt.AlignmentFlag.AlignRight
            return int(h | Qt.AlignmentFlag.AlignVCenter)
        if role == Qt.ItemDataRole.ForegroundRole:
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
            if role == Qt.ItemDataRole.ToolTipRole and orientation == Qt.Orientation.Horizontal:
                return self._cols[section].label or None
            return None
        if orientation == Qt.Orientation.Horizontal:
            return self._cols[section].name
        return str(self.obs(section) + 1)


class DataBrowser(QMainWindow):
    """Janela separada, como no Stata. Fica aberta e acompanha os dados."""

    def __init__(self, session: "Session", font=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Data Editor (Browse)")
        self.resize(900, 560)
        self.model = DatasetModel(session, self)
        self.table = QTableView()
        self.table.setModel(self.model)
        self.table.setEditTriggers(QTableView.EditTrigger.NoEditTriggers)
        self.table.setAlternatingRowColors(False)
        self.table.verticalHeader().setDefaultSectionSize(20)
        self.table.horizontalHeader().setDefaultSectionSize(96)
        if font is not None:
            self.table.setFont(font)
        self.info = QLabel()
        self.info.setContentsMargins(6, 2, 6, 2)
        box = QWidget()
        lay = QVBoxLayout(box)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.table, 1)
        lay.addWidget(self.info)
        self.setCentralWidget(box)
        QShortcut(QKeySequence.StandardKey.Close, self, self.close)
        self.table.selectionModel().currentChanged.connect(self._show_cell)
        self.refresh()

    def set_dark(self, dark: bool) -> None:
        self.model.dark = dark
        self.model.refresh()

    def refresh(self) -> None:
        self.model.refresh()
        ds = self.model.ds
        filt = "On" if self.model.rows is not None else "Off"
        self.info.setText(f"Vars: {self.model.columnCount():,}    Obs: {ds.nobs:,}    "
                          f"Filter: {filt}    Mode: Browse")

    def _show_cell(self, current: QModelIndex, _prev=None) -> None:
        if not current.isValid():
            return
        ds = self.model.ds
        var = self.model.var(current.column())
        i = self.model.obs(current.row())
        raw = cell_text(ds, var, i, use_labels=False)
        self.statusBar().showMessage(f"[{i + 1},{ds.index(var.name) + 1}]  {var.name} = {raw}")
