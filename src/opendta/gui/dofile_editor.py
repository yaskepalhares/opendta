"""Do-file Editor: tabs, line numbers, syntax highlighting, open and
save, and running the file or the selection (Do = do; Run = run, without echo).

As in Stata, running writes the text to a temporary file and calls
`do` (or `run`) on it; a selection is extended to whole lines."""

from __future__ import annotations

from ..i18n import tr

import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Callable

from PySide6.QtCore import QRect, QSize, Qt, Signal
from PySide6.QtGui import (QAction, QColor, QFont, QKeySequence, QPainter, QShortcut, QSyntaxHighlighter,
                           QTextCharFormat, QTextCursor, QTextDocument)
from PySide6.QtWidgets import (QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox,
                               QPlainTextEdit, QPushButton, QTabWidget, QTextEdit, QToolBar, QWidget)

if TYPE_CHECKING:
    from ..session import Session

# OpenDTA's own colors (light, dark)
PALETTE = {
    "light": {"command": "#1f4f99", "string": "#a0522d", "comment": "#3a7d44", "macro": "#8a3fa0",
              "number": "#b5651d", "gutter": "#f2f3f5", "gutter_text": "#8a9099", "current": "#f7f9fc"},
    "dark": {"command": "#7fb2f0", "string": "#e6a374", "comment": "#7fbf86", "macro": "#d29be0",
             "number": "#e8c07a", "gutter": "#262a30", "gutter_text": "#7a828c", "current": "#2a2f36"},
}

_EXTRA_KEYWORDS = ("if", "in", "using", "else", "foreach", "forvalues", "while", "of", "local", "global",
                   "program", "end", "define", "capture", "quietly", "noisily", "by", "bysort", "version",
                   "mata", "tempvar", "tempname", "tempfile", "args", "syntax", "return", "ereturn",
                   "preserve", "restore", "exit", "continue", "break")


def _keywords() -> list[str]:
    from ..commands.registry import REGISTRY
    return sorted(set(REGISTRY) | set(_EXTRA_KEYWORDS))


class StataHighlighter(QSyntaxHighlighter):
    """Highlighting: commands, strings, comments (*, //, ///, /* */), macros and numbers."""

    def __init__(self, doc: QTextDocument, theme: str = "light"):
        super().__init__(doc)
        words = "|".join(re.escape(w) for w in sorted(_keywords(), key=len, reverse=True))
        self.re_cmd = re.compile(rf"\b(?:{words})\b")
        self.re_str = re.compile(r'`"(?:[^"]|"(?!\'))*"\'|"[^"\n]*"')
        self.re_macro = re.compile(r"`[A-Za-z_0-9:()\s.]*'|\$\{?[A-Za-z_][A-Za-z_0-9]*\}?")
        self.re_num = re.compile(r"(?<![\w.])(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?\b")
        self.set_theme(theme)

    def set_theme(self, theme: str) -> None:
        pal = PALETTE["dark" if theme == "dark" else "light"]
        self.fmt = {}
        for key in ("command", "string", "comment", "macro", "number"):
            f = QTextCharFormat()
            f.setForeground(QColor(pal[key]))
            if key == "command":
                f.setFontWeight(QFont.Weight.DemiBold)
            if key == "comment":
                f.setFontItalic(True)
            self.fmt[key] = f
        self.rehighlight()

    def highlightBlock(self, text: str) -> None:
        self.setCurrentBlockState(0)
        start = 0
        if self.previousBlockState() == 1:
            end = text.find("*/")
            if end == -1:
                self.setFormat(0, len(text), self.fmt["comment"])
                self.setCurrentBlockState(1)
                return
            self.setFormat(0, end + 2, self.fmt["comment"])
            start = end + 2
        stripped = text.lstrip()
        if start == 0 and stripped.startswith("*"):
            self.setFormat(0, len(text), self.fmt["comment"])
            return
        for m in self.re_num.finditer(text, start):
            self.setFormat(m.start(), m.end() - m.start(), self.fmt["number"])
        for m in self.re_cmd.finditer(text, start):
            self.setFormat(m.start(), m.end() - m.start(), self.fmt["command"])
        first = re.match(r"\s*([A-Za-z_]\w*)", text[start:])
        if first:
            from ..commands.registry import lookup
            if lookup(first.group(1)) is not None:
                self.setFormat(start + first.start(1), len(first.group(1)), self.fmt["command"])
        for m in self.re_macro.finditer(text, start):
            self.setFormat(m.start(), m.end() - m.start(), self.fmt["macro"])
        strings = [(m.start(), m.end()) for m in self.re_str.finditer(text, start)]
        for a, b in strings:
            self.setFormat(a, b - a, self.fmt["string"])

        def in_string(pos: int) -> bool:
            return any(a <= pos < b for a, b in strings)

        i = start
        while i < len(text):
            if in_string(i):
                i += 1
                continue
            if text.startswith("/*", i):
                end = text.find("*/", i + 2)
                if end == -1:
                    self.setFormat(i, len(text) - i, self.fmt["comment"])
                    self.setCurrentBlockState(1)
                    return
                self.setFormat(i, end + 2 - i, self.fmt["comment"])
                i = end + 2
                continue
            if text.startswith("//", i) and (i == 0 or text[i - 1].isspace()):
                self.setFormat(i, len(text) - i, self.fmt["comment"])
                return
            i += 1


class _Gutter(QWidget):
    def __init__(self, editor: "CodeEditor"):
        super().__init__(editor)
        self.editor = editor

    def sizeHint(self) -> QSize:
        return QSize(self.editor.gutter_width(), 0)

    def paintEvent(self, event) -> None:
        self.editor.paint_gutter(event)


class CodeEditor(QPlainTextEdit):
    def __init__(self, font: QFont, theme: str = "light"):
        super().__init__()
        self.theme = theme
        self.path: Path | None = None
        self.setFont(font)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setTabStopDistance(self.fontMetrics().horizontalAdvance(" ") * 4)
        self.gutter = _Gutter(self)
        self.highlighter = StataHighlighter(self.document(), theme)
        self.blockCountChanged.connect(lambda _: self._update_margin())
        self.updateRequest.connect(self._update_gutter)
        self.cursorPositionChanged.connect(self._highlight_line)
        self._update_margin()
        self._highlight_line()

    def gutter_width(self) -> int:
        digits = max(3, len(str(max(1, self.blockCount()))))
        return 12 + self.fontMetrics().horizontalAdvance("9") * digits

    def _update_margin(self) -> None:
        self.setViewportMargins(self.gutter_width(), 0, 0, 0)

    def _update_gutter(self, rect, dy) -> None:
        if dy:
            self.gutter.scroll(0, dy)
        else:
            self.gutter.update(0, rect.y(), self.gutter.width(), rect.height())

    def resizeEvent(self, e) -> None:
        super().resizeEvent(e)
        cr = self.contentsRect()
        self.gutter.setGeometry(QRect(cr.left(), cr.top(), self.gutter_width(), cr.height()))

    def paint_gutter(self, event) -> None:
        pal = PALETTE["dark" if self.theme == "dark" else "light"]
        p = QPainter(self.gutter)
        p.fillRect(event.rect(), QColor(pal["gutter"]))
        p.setPen(QColor(pal["gutter_text"]))
        p.setFont(self.font())
        block = self.firstVisibleBlock()
        n = block.blockNumber()
        top = round(self.blockBoundingGeometry(block).translated(self.contentOffset()).top())
        bottom = top + round(self.blockBoundingRect(block).height())
        h = self.fontMetrics().height()
        while block.isValid() and top <= event.rect().bottom():
            if block.isVisible() and bottom >= event.rect().top():
                p.drawText(0, top, self.gutter.width() - 6, h, Qt.AlignmentFlag.AlignRight, str(n + 1))
            block = block.next()
            top = bottom
            bottom = top + round(self.blockBoundingRect(block).height())
            n += 1

    def _highlight_line(self) -> None:
        pal = PALETTE["dark" if self.theme == "dark" else "light"]
        sel = QTextEdit.ExtraSelection()
        sel.format.setBackground(QColor(pal["current"]))
        sel.format.setProperty(QTextCharFormat.Property.FullWidthSelection, True)
        sel.cursor = self.textCursor()
        sel.cursor.clearSelection()
        self.setExtraSelections([sel])

    def set_theme(self, theme: str) -> None:
        self.theme = theme
        self.highlighter.set_theme(theme)
        self._highlight_line()
        self.gutter.update()

    def code_to_run(self) -> str:
        """Selection extended to whole lines; without a selection, the whole file."""
        cur = self.textCursor()
        if not cur.hasSelection():
            return self.toPlainText()
        doc = self.document()
        a = doc.findBlock(cur.selectionStart())
        b = doc.findBlock(cur.selectionEnd())
        if cur.selectionEnd() == b.position() and b.blockNumber() > a.blockNumber():
            b = b.previous()
        lines = []
        blk = a
        while blk.isValid():
            lines.append(blk.text())
            if blk == b:
                break
            blk = blk.next()
        return "\n".join(lines)


class _FindBar(QWidget):
    def __init__(self, owner: "DoFileEditor"):
        super().__init__()
        self.owner = owner
        lay = QHBoxLayout(self)
        lay.setContentsMargins(6, 2, 6, 2)
        lay.addWidget(QLabel(tr("Find:")))
        self.find = QLineEdit()
        self.find.returnPressed.connect(lambda: owner.find_next(self.find.text()))
        lay.addWidget(self.find)
        lay.addWidget(QLabel(tr("Replace:")))
        self.repl = QLineEdit()
        lay.addWidget(self.repl)
        for text, slot in (("Next", lambda: owner.find_next(self.find.text())),
                           ("Replace", lambda: owner.replace_one(self.find.text(), self.repl.text())),
                           ("Replace all", lambda: owner.replace_all(self.find.text(), self.repl.text())),
                           ("✕", self.hide)):
            b = QPushButton(tr(text))
            b.clicked.connect(slot)
            lay.addWidget(b)


class DoFileEditor(QMainWindow):
    """Do-file Editor window."""

    ran = Signal()

    def __init__(self, session: "Session", parent=None, *, run: Callable[[str], None],
                 font: QFont | None = None, theme: str = "light"):
        super().__init__(parent)
        self.session = session
        self.run = run
        self.theme = theme
        self.editor_font = font or QFont("Menlo", 12)
        self.resize(860, 640)
        self.untitled = 0
        self.tmpfiles: list[str] = []

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        self.tabs.currentChanged.connect(lambda _: self._update_title())
        self.findbar = _FindBar(self)
        self.findbar.hide()
        central = QWidget()
        from PySide6.QtWidgets import QVBoxLayout
        lay = QVBoxLayout(central)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.tabs)
        lay.addWidget(self.findbar)
        self.setCentralWidget(central)

        tb = QToolBar("Do-file Editor")
        tb.setMovable(False)
        acts = [("New", self.new_file, QKeySequence.StandardKey.New),
                ("Open...", self.open_dialog, QKeySequence.StandardKey.Open),
                ("Save", self.save, QKeySequence.StandardKey.Save),
                ("Save as...", self.save_as, QKeySequence.StandardKey.SaveAs),
                (None, None, None),
                ("Find", self.show_find, QKeySequence.StandardKey.Find),
                (None, None, None),
                ("Do", self.do_selection, "Ctrl+D"),
                ("Run", self.run_selection, "Ctrl+R")]
        for text, slot, key in acts:
            if text is None:
                tb.addSeparator()
                continue
            act = QAction(tr(text), self)
            act.setProperty("key", text)
            act.triggered.connect(slot)
            if key is not None:
                act.setShortcut(QKeySequence(key))
            tip = {"Do": "Execute the selection, or the whole file (do)",
                   "Run": "Execute the selection, or the whole file, without showing commands (run)"}.get(text)
            if tip:
                act.setToolTip(tr(tip) + f"  [{QKeySequence(key).toString(QKeySequence.SequenceFormat.NativeText)}]")
            tb.addAction(act)
        self.addToolBar(tb)
        self.toolbar = tb
        QShortcut(QKeySequence("Ctrl+Shift+D"), self, self.do_selection)
        QShortcut(QKeySequence("Escape"), self, self.findbar.hide)
        self.new_file()

    # -- abas -----------------------------------------------------------------------------
    def current(self) -> CodeEditor | None:
        w = self.tabs.currentWidget()
        return w if isinstance(w, CodeEditor) else None

    def _add(self, ed: CodeEditor, title: str) -> None:
        idx = self.tabs.addTab(ed, title)
        self.tabs.setCurrentIndex(idx)
        ed.document().modificationChanged.connect(lambda _m, e=ed: self._refresh_tab(e))
        ed.setFocus()
        self._update_title()

    def _tab_title(self, ed: CodeEditor) -> str:
        name = ed.path.name if ed.path else ed.property("untitled")
        return name + ("*" if ed.document().isModified() else "")

    def _refresh_tab(self, ed: CodeEditor) -> None:
        idx = self.tabs.indexOf(ed)
        if idx >= 0:
            self.tabs.setTabText(idx, self._tab_title(ed))
        self._update_title()

    def _update_title(self) -> None:
        ed = self.current()
        if ed is None:
            self.setWindowTitle(tr("Do-file Editor"))
            return
        self.setWindowTitle(f"{tr('Do-file Editor')} - {self._tab_title(ed)}")

    def new_file(self) -> None:
        self.untitled += 1
        ed = CodeEditor(self.editor_font, self.theme)
        ed.setProperty("untitled", f"Untitled{self.untitled}.do")
        self._add(ed, ed.property("untitled"))

    def open_path(self, path: str | Path) -> None:
        p = Path(path).expanduser()
        if p.suffix == "":
            p = p.with_suffix(".do")
        for i in range(self.tabs.count()):
            ed = self.tabs.widget(i)
            if isinstance(ed, CodeEditor) and ed.path and ed.path.resolve() == p.resolve():
                self.tabs.setCurrentIndex(i)
                return
        ed = CodeEditor(self.editor_font, self.theme)
        ed.path = p
        if p.exists():
            ed.setPlainText(p.read_text(encoding="utf-8", errors="replace"))
        ed.document().setModified(False)
        cur = self.current()
        # the initial empty tab is replaced by the opened file
        if cur is not None and cur.path is None and not cur.toPlainText() and not cur.document().isModified():
            self.tabs.removeTab(self.tabs.indexOf(cur))
        self._add(ed, self._tab_title(ed))

    def open_dialog(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(self, tr("Open"), os.getcwd(),
                                                tr("Do-files (*.do *.ado);;All files (*)"))
        for p in paths:
            self.open_path(p)

    def save(self) -> bool:
        ed = self.current()
        if ed is None:
            return False
        if ed.path is None:
            return self.save_as()
        ed.path.write_text(ed.toPlainText(), encoding="utf-8")
        ed.document().setModified(False)
        return True

    def save_as(self) -> bool:
        ed = self.current()
        if ed is None:
            return False
        start = str(ed.path) if ed.path else os.path.join(os.getcwd(), ed.property("untitled"))
        path, _ = QFileDialog.getSaveFileName(self, tr("Save as"), start, "Do-files (*.do);;Ado-files (*.ado)")
        if not path:
            return False
        ed.path = Path(path)
        return self.save()

    def close_tab(self, idx: int) -> bool:
        ed = self.tabs.widget(idx)
        if isinstance(ed, CodeEditor) and ed.document().isModified():
            self.tabs.setCurrentIndex(idx)
            ans = QMessageBox.question(self, tr("Do-file Editor"),
                                       tr("Save changes to {name}?").format(name=self._tab_title(ed).rstrip('*')),
                                       QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard
                                       | QMessageBox.StandardButton.Cancel)
            if ans == QMessageBox.StandardButton.Cancel:
                return False
            if ans == QMessageBox.StandardButton.Save and not self.save():
                return False
        self.tabs.removeTab(idx)
        if self.tabs.count() == 0:
            self.hide()
            self.new_file()
        return True

    def closeEvent(self, e) -> None:
        for i in reversed(range(self.tabs.count())):
            ed = self.tabs.widget(i)
            if isinstance(ed, CodeEditor) and ed.document().isModified():
                if not self.close_tab(i):
                    e.ignore()
                    return
        e.accept()

    # -- busca ------------------------------------------------------------------------------
    def show_find(self) -> None:
        self.findbar.show()
        ed = self.current()
        if ed is not None and ed.textCursor().hasSelection():
            self.findbar.find.setText(ed.textCursor().selectedText())
        self.findbar.find.setFocus()
        self.findbar.find.selectAll()

    def find_next(self, text: str) -> bool:
        ed = self.current()
        if ed is None or not text:
            return False
        if ed.find(text):
            return True
        ed.moveCursor(QTextCursor.MoveOperation.Start)
        return ed.find(text)

    def replace_one(self, text: str, repl: str) -> None:
        ed = self.current()
        if ed is None or not text:
            return
        cur = ed.textCursor()
        if cur.hasSelection() and cur.selectedText() == text:
            cur.insertText(repl)
        self.find_next(text)

    def replace_all(self, text: str, repl: str) -> int:
        ed = self.current()
        if ed is None or not text:
            return 0
        body = ed.toPlainText()
        n = body.count(text)
        if n:
            cur = ed.textCursor()
            cur.beginEditBlock()
            cur.select(QTextCursor.SelectionType.Document)
            cur.insertText(body.replace(text, repl))
            cur.endEditBlock()
        return n

    # -- execution ------------------------------------------------------------------------
    def _execute(self, verb: str) -> None:
        ed = self.current()
        if ed is None:
            return
        code = ed.code_to_run()
        if not code.strip():
            return
        fd, path = tempfile.mkstemp(prefix="SD", suffix=".tmp")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(code if code.endswith("\n") else code + "\n")
        self.tmpfiles.append(path)
        self.run(f'{verb} "{path}"')
        self.ran.emit()

    def do_selection(self) -> None:
        self._execute("do")

    def run_selection(self) -> None:
        self._execute("run")

    def set_theme(self, theme: str) -> None:
        self.theme = theme
        for i in range(self.tabs.count()):
            ed = self.tabs.widget(i)
            if isinstance(ed, CodeEditor):
                ed.set_theme(theme)

    def cleanup(self) -> None:
        for p in self.tmpfiles:
            try:
                os.unlink(p)
            except OSError:
                pass
        self.tmpfiles.clear()

    def retranslate(self) -> None:
        for act in self.toolbar.actions():
            key = act.property("key")
            if key:
                act.setText(tr(key))
        self._update_title()
