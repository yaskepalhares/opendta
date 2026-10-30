"""Rounded look shared by the secondary windows (tabs, buttons, inputs), matching the main window.

Neutral translucent greys, so it works in the light and dark themes without reading the palette.
"""

from __future__ import annotations

from PySide6.QtWidgets import QWidget

ROUNDED_QSS = """
QTabWidget::pane { border: none; }
QTabBar { background: transparent; }
QTabBar::tab { background: transparent; border: none; border-radius: 8px; padding: 5px 12px; margin: 4px 2px; }
QTabBar::tab:selected { background: rgba(127, 127, 127, 0.22); font-weight: 600; }
QTabBar::tab:hover:!selected { background: rgba(127, 127, 127, 0.12); }
QPushButton { border: 1px solid rgba(127, 127, 127, 0.4); border-radius: 8px; padding: 5px 14px; }
QPushButton:hover { background: rgba(127, 127, 127, 0.16); }
QPushButton:pressed { background: rgba(127, 127, 127, 0.28); }
QPushButton:default { background: #0a84ff; border: 1px solid #0a84ff; color: #ffffff; }
QPushButton:default:hover { background: #3b9bff; }
QPushButton:disabled { color: rgba(127, 127, 127, 0.7); }
QLineEdit, QSpinBox, QDoubleSpinBox { border: 1px solid rgba(127, 127, 127, 0.4); border-radius: 8px;
                                      padding: 4px 8px; }
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus { border: 1px solid #0a84ff; }
QComboBox { border: 1px solid rgba(127, 127, 127, 0.4); border-radius: 8px; padding: 4px 8px; }
QComboBox:hover { background: rgba(127, 127, 127, 0.12); }
QToolBar { border: none; spacing: 2px; padding: 3px 8px; }
QToolBar QToolButton { border: none; border-radius: 8px; padding: 4px 8px; }
QToolBar QToolButton:hover { background: rgba(127, 127, 127, 0.16); }
QToolBar QToolButton:pressed { background: rgba(127, 127, 127, 0.28); }
"""


def apply_rounded(widget: QWidget) -> None:
    """Add the rounded look to the widget's own stylesheet."""
    widget.setStyleSheet(widget.styleSheet() + ROUNDED_QSS)
