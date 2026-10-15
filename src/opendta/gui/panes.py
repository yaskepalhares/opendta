"""The panes of the main window: Results and Command, and the monospace font of the text
windows."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QTextCharFormat, QTextCursor
from PySide6.QtWidgets import QLabel, QPlainTextEdit, QVBoxLayout, QWidget

from ..i18n import tr
from .settings import Preferences
from .theme import MONOSPACE_FALLBACKS, STANDARD, ResultsScheme


def monospace_font(prefs: Preferences | None = None) -> QFont:
    """Font of the text windows: the one chosen in the preferences, or the system
    default (Menlo on macOS, Courier New on Windows, DejaVu Sans Mono on Linux)."""
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
        self._segments: list[tuple[str, str]] = []   # to redraw when the theme changes
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
    """Command window: resizable and multi-line.

    Enter runs; Shift+Enter breaks the line. Pasted text with several
    lines is run line by line. PgUp/PgDn walk through the history; so do the
    ↑/↓ arrows, while the command has a single line."""

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
        # at least one visible line, whatever the font
        self.setMinimumHeight(self.fontMetrics().lineSpacing() + 14)

    # compatibility with the QLineEdit API used by the window
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


def titled(title: str, widget: QWidget) -> QWidget:
    """Panel with a thin title bar, like Stata's docked panels."""
    box = QWidget()
    lay = QVBoxLayout(box)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(0)
    label = QLabel(f"  {tr(title)}")
    label.setObjectName("PaneTitle")
    label.setProperty("key", title)
    lay.addWidget(label)
    lay.addWidget(widget)
    return box
