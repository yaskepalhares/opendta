"""Variables Manager: propriedades de todas as variáveis numa tabela.

Cada alteração vira um comando (rename, label variable, recast, format,
label values) executado pela sessão, ecoado em Results e guardado em
Review, como no Stata.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSortFilterProxyModel, Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import QLineEdit, QMainWindow, QTableView, QVBoxLayout, QWidget

if TYPE_CHECKING:
    from ..session import Session

COLUMNS = ("Name", "Label", "Type", "Format", "Value label")


def _quote(text: str) -> str:
    if '"' in text:
        return '`"' + text + '"\''
    return '"' + text + '"'


def property_command(name: str, column: int, value: str) -> str | None:
    """Comando que muda a propriedade `column` (índice em COLUMNS) da variável."""
    v = value.strip()
    if column == 0:
        return f"rename {name} {v}" if v and v != name else None
    if column == 1:
        return f"label variable {name} {_quote(value)}" if value else f"label variable {name}"
    if column == 2:
        return f"recast {v} {name}" if v else None
    if column == 3:
        return f"format {name} {v}" if v else None
    if column == 4:
        return f"label values {name} {v}".rstrip()
    return None


class VariablesModel(QAbstractTableModel):
    def __init__(self, session: "Session", run: Callable[[str], int] | None, parent=None):
        super().__init__(parent)
        self.session = session
        self.run = run

    def refresh(self) -> None:
        self.beginResetModel()
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else self.session.data.nvars

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if role == Qt.ItemDataRole.DisplayRole and orientation == Qt.Orientation.Horizontal:
            return COLUMNS[section]
        return None

    def _value(self, row: int, col: int) -> str:
        v = self.session.data.vars[row]
        return (v.name, v.label, v.vtype, v.fmt, v.value_label)[col]

    def data(self, index: QModelIndex, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        if role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.EditRole):
            return self._value(index.row(), index.column())
        return None

    def flags(self, index: QModelIndex):
        base = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        return base | Qt.ItemFlag.ItemIsEditable if self.run is not None else base

    def setData(self, index: QModelIndex, value, role=Qt.ItemDataRole.EditRole) -> bool:  # noqa: N802
        if role != Qt.ItemDataRole.EditRole or self.run is None:
            return False
        text = str(value)
        if text == self._value(index.row(), index.column()):
            return False
        name = self.session.data.vars[index.row()].name
        cmd = property_command(name, index.column(), text)
        if cmd is None:
            return False
        return self.run(cmd) == 0


class VariablesManager(QMainWindow):
    def __init__(self, session: "Session", font=None, parent=None,
                 run: Callable[[str], int] | None = None):
        super().__init__(parent)
        self.setWindowTitle("Variables Manager")
        self.resize(760, 460)
        self.model = VariablesModel(session, run, self)
        self.proxy = QSortFilterProxyModel(self)
        self.proxy.setSourceModel(self.model)
        self.proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self.proxy.setFilterKeyColumn(-1)
        self.filter = QLineEdit()
        self.filter.setPlaceholderText("Filter variables here")
        self.filter.textChanged.connect(self.proxy.setFilterFixedString)
        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.verticalHeader().setVisible(False)
        self.table.verticalHeader().setDefaultSectionSize(20)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setColumnWidth(0, 140)
        self.table.setColumnWidth(1, 260)
        if font is not None:
            self.table.setFont(font)
        box = QWidget()
        lay = QVBoxLayout(box)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self.filter)
        lay.addWidget(self.table, 1)
        self.setCentralWidget(box)
        QShortcut(QKeySequence.StandardKey.Close, self, self.close)

    def refresh(self) -> None:
        self.model.refresh()
