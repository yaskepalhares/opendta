"""The panes of the main window: Results and Command, and the monospace font of the text
windows."""

from __future__ import annotations

import re

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QPainter, QPen, QTextCharFormat, QTextCursor
from PySide6.QtWidgets import QLabel, QPlainTextEdit, QVBoxLayout, QWidget

from ..i18n import tr
from .rules import DOWN, LEFT, RIGHT, UP, table_cells
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

# "(519,227 observations deleted)": shown in the "alert" style, unless nothing was deleted
_DELETED = re.compile(r"\(([\d,]+) observations? deleted\)")


class ResultsView(QPlainTextEdit):
    """Results window. Text in the "hint_link" style (the "help topic" of an error explanation)
    is a link: clicking it emits linkActivated with that command."""

    MAX_SEGMENTS = 200_000
    linkActivated = Signal(str)

    def __init__(self, scheme: ResultsScheme = STANDARD):
        super().__init__()
        self.viewport().setMouseTracking(True)
        self.setReadOnly(True)
        self.setUndoRedoEnabled(False)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setFont(monospace_font())
        self.setObjectName("Results")
        self._segments: list[tuple[str, str]] = []   # to redraw when the theme changes
        self._solid_rules = True
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
            if style == "hint_link":
                fmt.setFontUnderline(True)
                fmt.setAnchor(True)
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
        fmt = self._formats.get(style, self._formats["text"])
        if style == "hint_link":
            fmt = QTextCharFormat(fmt)
            fmt.setAnchorHref(text)
        pos = 0
        if style != "hint_link" and "deleted)" in text:
            for m in _DELETED.finditer(text):
                if int(m.group(1).replace(",", "")) == 0:
                    continue
                cursor.insertText(text[pos:m.start()], fmt)
                cursor.insertText(m.group(), self._formats["alert"])
                pos = m.end()
        cursor.insertText(text[pos:], fmt)
        self.setTextCursor(cursor)
        self.ensureCursorVisible()

    def set_solid_rules(self, on: bool) -> None:
        """Draw the lines of the tables as continuous lines (off: the text's dashes, as in Stata)."""
        self._solid_rules = on
        self.viewport().update()

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        if self._solid_rules:
            self._paint_rules(event.rect().bottom())

    def _paint_rules(self, last_y: int) -> None:
        """Covers the dashes, bars and crossings of the visible lines with continuous lines. The text
        is not changed: only what is drawn."""
        painter = QPainter(self.viewport())
        painter.setPen(QPen(QColor(self._scheme.colors["text"]), 1))
        block = self.firstVisibleBlock()
        top = self.blockBoundingGeometry(block).translated(self.contentOffset()).top()
        while block.isValid() and top <= last_y:
            height = self.blockBoundingRect(block).height()
            text = block.text()
            if "-" in text or "|" in text:
                before, after = block.previous(), block.next()
                cells = table_cells(before.text() if before.isValid() else None, text,
                                    after.text() if after.isValid() else None)
                if cells:
                    self._paint_cells(painter, block, round(top), round(height), cells)
            top += height
            block = block.next()
        painter.end()

    def _paint_cells(self, painter: QPainter, block, top: int, height: int, cells) -> None:
        cursor = QTextCursor(block)
        edges: dict[int, int] = {}

        def edge(col: int) -> int:
            if col not in edges:
                cursor.setPosition(block.position() + col)
                edges[col] = self.cursorRect(cursor).left()
            return edges[col]

        middle = top + height // 2
        bottom = top + height
        background = QColor(self._scheme.background)
        for col, sides in cells:
            left, right = edge(col), edge(col + 1)
            centre = (left + right) // 2
            painter.fillRect(left, top, right - left, height, background)
            if sides & LEFT:
                painter.drawLine(left, middle, centre, middle)
            if sides & RIGHT:
                painter.drawLine(centre, middle, right, middle)
            if sides & UP:
                painter.drawLine(centre, top, centre, middle)
            if sides & DOWN:
                painter.drawLine(centre, middle, centre, bottom)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        over = bool(self.anchorAt(event.position().toPoint()))
        self.viewport().setCursor(Qt.CursorShape.PointingHandCursor if over else Qt.CursorShape.IBeamCursor)
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        super().mouseReleaseEvent(event)
        href = self.anchorAt(event.position().toPoint())
        if href and event.button() == Qt.MouseButton.LeftButton and not self.textCursor().hasSelection():
            self.linkActivated.emit(href)


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
