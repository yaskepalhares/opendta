"""Diálogo Edit > Preferences (primeira versão: aparência e fonte)."""

from __future__ import annotations

from PySide6.QtCore import QSize, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (QButtonGroup, QDialog, QDialogButtonBox,
                               QFontComboBox, QFormLayout, QGroupBox,
                               QHBoxLayout, QLabel, QRadioButton, QSpinBox,
                               QVBoxLayout)

from .icons import app_icon
from .settings import APP_ICON_VARIANTS, Preferences
from .theme import default_monospace


class PreferencesDialog(QDialog):
    applied = Signal()

    def __init__(self, prefs: Preferences, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.prefs = prefs
        lay = QVBoxLayout(self)

        # ícone do aplicativo
        box = QGroupBox("Ícone do aplicativo")
        row = QHBoxLayout(box)
        self.icon_group = QButtonGroup(self)
        for key, label in APP_ICON_VARIANTS.items():
            rb = QRadioButton(label)
            rb.setIcon(app_icon(key))
            rb.setIconSize(QSize(48, 48))
            rb.setProperty("variant", key)
            rb.setChecked(prefs.app_icon == key)
            self.icon_group.addButton(rb)
            row.addWidget(rb)
        lay.addWidget(box)

        # fonte
        box = QGroupBox("Fonte das janelas Results, Command e Review")
        form = QFormLayout(box)
        self.family = QFontComboBox()
        self.family.setFontFilters(QFontComboBox.FontFilter.MonospacedFonts)
        self.family.setCurrentFont(QFont(prefs.font_family))
        self.size_box = QSpinBox()
        self.size_box.setRange(6, 36)
        self.size_box.setValue(prefs.font_size)
        fam, size = default_monospace()
        form.addRow("Família", self.family)
        form.addRow("Tamanho", self.size_box)
        form.addRow("", QLabel(f"Padrão neste sistema: {fam} {size}"))
        lay.addWidget(box)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok
                                   | QDialogButtonBox.StandardButton.Cancel
                                   | QDialogButtonBox.StandardButton.RestoreDefaults)
        buttons.accepted.connect(self._ok)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(self._defaults)
        lay.addWidget(buttons)

    def _defaults(self) -> None:
        fam, size = default_monospace()
        self.family.setCurrentFont(QFont(fam))
        self.size_box.setValue(size)
        for b in self.icon_group.buttons():
            b.setChecked(b.property("variant") == "light")

    def _ok(self) -> None:
        checked = self.icon_group.checkedButton()
        self.prefs.app_icon = checked.property("variant") if checked else "light"
        self.prefs.font_family = self.family.currentFont().family()
        self.prefs.font_size = self.size_box.value()
        self.prefs.sync()
        self.applied.emit()
        self.accept()
