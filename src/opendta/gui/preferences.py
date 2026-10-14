"""Diálogo Edit > Preferences: aparência, ícone, fonte e idioma."""

from __future__ import annotations

from PySide6.QtCore import QSize, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDialog,
                               QDialogButtonBox, QFontComboBox, QFormLayout,
                               QGroupBox, QHBoxLayout, QLabel, QRadioButton,
                               QSpinBox, QVBoxLayout)

from ..i18n import LANGUAGES, tr
from .icons import app_icon
from .settings import APP_ICON_VARIANTS, Preferences
from .theme import default_monospace


class PreferencesDialog(QDialog):
    applied = Signal()

    def __init__(self, prefs: Preferences, parent=None):
        super().__init__(parent)
        self.setWindowTitle(tr("Preferences"))
        self.prefs = prefs
        lay = QVBoxLayout(self)

        # tema da interface
        box = QGroupBox(tr("Interface appearance"))
        row = QHBoxLayout(box)
        self.theme_group = QButtonGroup(self)
        for key, label in (("light", "Light (default)"), ("dark", "Dark")):
            rb = QRadioButton(tr(label))
            rb.setProperty("theme", key)
            rb.setChecked(prefs.theme == key)
            self.theme_group.addButton(rb)
            row.addWidget(rb)
        row.addStretch(1)
        lay.addWidget(box)

        # ícone do aplicativo
        box = QGroupBox(tr("Application icon"))
        row = QHBoxLayout(box)
        self.icon_group = QButtonGroup(self)
        for key, label in APP_ICON_VARIANTS.items():
            rb = QRadioButton(tr(label))
            rb.setIcon(app_icon(key))
            rb.setIconSize(QSize(48, 48))
            rb.setProperty("variant", key)
            rb.setChecked(prefs.app_icon == key)
            self.icon_group.addButton(rb)
            row.addWidget(rb)
        lay.addWidget(box)

        # fonte
        box = QGroupBox(tr("Font of the Results, Command and Review windows"))
        form = QFormLayout(box)
        self.family = QFontComboBox()
        self.family.setFontFilters(QFontComboBox.FontFilter.MonospacedFonts)
        self.family.setCurrentFont(QFont(prefs.font_family))
        self.size_box = QSpinBox()
        self.size_box.setRange(6, 36)
        self.size_box.setValue(prefs.font_size)
        fam, size = default_monospace()
        form.addRow(tr("Family"), self.family)
        form.addRow(tr("Size"), self.size_box)
        form.addRow("", QLabel(f"{tr('Default on this system')}: {fam} {size}"))
        lay.addWidget(box)

        # idioma
        box = QGroupBox(tr("Language"))
        form = QFormLayout(box)
        self.lang_box = QComboBox()
        for key, label in LANGUAGES.items():
            self.lang_box.addItem(label, key)
        self.lang_box.setCurrentIndex(max(0, self.lang_box.findData(prefs.language)))
        form.addRow(tr("Interface language"), self.lang_box)
        self.translate_msgs = QCheckBox(tr("Also translate Stata-style messages (errors)"))
        self.translate_msgs.setChecked(prefs.messages == "translated")
        self.translate_msgs.setToolTip(tr("Keep them in English to compare with Stata logs and manuals."))
        form.addRow("", self.translate_msgs)
        note = QLabel(tr("Keep them in English to compare with Stata logs and manuals."))
        note.setWordWrap(True)
        note.setEnabled(False)
        form.addRow("", note)
        self.lang_box.currentIndexChanged.connect(self._lang_changed)
        self._lang_changed()
        lay.addWidget(box)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok
                                   | QDialogButtonBox.StandardButton.Cancel
                                   | QDialogButtonBox.StandardButton.RestoreDefaults)
        buttons.accepted.connect(self._ok)
        buttons.rejected.connect(self.reject)
        buttons.button(QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(self._defaults)
        lay.addWidget(buttons)

    def _lang_changed(self) -> None:
        # a tradução das mensagens do Stata só existe em português
        self.translate_msgs.setEnabled(self.lang_box.currentData() == "pt")

    def _defaults(self) -> None:
        fam, size = default_monospace()
        self.family.setCurrentFont(QFont(fam))
        self.size_box.setValue(size)
        for b in self.icon_group.buttons():
            b.setChecked(b.property("variant") == "light")
        for b in self.theme_group.buttons():
            b.setChecked(b.property("theme") == "light")
        self.lang_box.setCurrentIndex(max(0, self.lang_box.findData("en")))
        self.translate_msgs.setChecked(False)

    def _ok(self) -> None:
        theme = self.theme_group.checkedButton()
        self.prefs.theme = theme.property("theme") if theme else "light"
        checked = self.icon_group.checkedButton()
        self.prefs.app_icon = checked.property("variant") if checked else "light"
        self.prefs.font_family = self.family.currentFont().family()
        self.prefs.font_size = self.size_box.value()
        self.prefs.language = self.lang_box.currentData() or "en"
        self.prefs.messages = ("translated" if self.translate_msgs.isChecked()
                               and self.prefs.language == "pt" else "stata")
        self.prefs.sync()
        self.applied.emit()
        self.accept()
