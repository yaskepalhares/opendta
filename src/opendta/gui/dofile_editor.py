"""Do-file Editor: tabs, line numbers, syntax highlighting, open and
save, and running the file or the selection (Do = do; Run = run, without echo).

As in Stata, running writes the text to a temporary file and calls
`do` (or `run`) on it; a selection is extended to whole lines."""

from __future__ import annotations

from ..i18n import tr
from .icons import icon

import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Callable

from PySide6.QtCore import QRect, QSize, Qt, Signal
from PySide6.QtGui import (QAction, QColor, QFont, QKeySequence, QPainter, QShortcut, QSyntaxHighlighter,
                           QTextCharFormat, QTextCursor, QTextDocument)
from PySide6.QtWidgets import (QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox,
                               QPlainTextEdit, QPushButton, QSizePolicy, QTabWidget, QTextEdit, QToolBar,
                               QToolButton, QWidget)

if TYPE_CHECKING:
    from ..session import Session

# OpenDTA's own colors (light, dark)
PALETTE = {
    "light": {"command": "#1f4f99", "string": "#a0522d", "comment": "#3a7d44", "macro": "#8a3fa0",
              "number": "#b5651d", "gutter": "#f2f3f5", "gutter_text": "#8a9099", "current": "#f7f9fc"},
    "dark": {"command": "#7fb2f0", "string": "#e6a374", "comment": "#7fbf86", "macro": "#d29be0",
             "number": "#e8c07a", "gutter": "#262a30", "gutter_text": "#7a828c", "current": "#2a2f36"},
}

RUN_GREEN = "#2fa84f"      # Do / Run icons

# short name shown under each toolbar icon, as in the main window (the tooltip keeps the long text)
_SHORT_LABELS = {"Open...": "Open", "Save as...": "Save as", "Comment/uncomment": "Comment"}


def editor_colors(palette: str, theme: str) -> dict:
    """Colors of the editor for a color scheme (palettes.py). "classic" is the original look of the
    theme and has no background of its own (the application's is used)."""
    from .palettes import resolve
    pal = resolve(palette, theme)
    if pal is None:
        return dict(PALETTE["dark" if theme == "dark" else "light"])
    return {"command": pal.command, "string": pal.string, "comment": pal.comment, "macro": pal.macro,
            "number": pal.number, "gutter": pal.gutter, "gutter_text": pal.gutter_text,
            "current": pal.current, "bg": pal.bg, "fg": pal.fg}

_EXTRA_KEYWORDS = ("if", "in", "using", "else", "foreach", "forvalues", "while", "of", "local", "global",
                   "program", "end", "define", "capture", "quietly", "noisily", "by", "bysort", "version",
                   "mata", "tempvar", "tempname", "tempfile", "args", "syntax", "return", "ereturn",
                   "preserve", "restore", "exit", "continue", "break")


def _keywords() -> list[str]:
    from ..commands.registry import REGISTRY
    return sorted(set(REGISTRY) | set(_EXTRA_KEYWORDS))


class StataHighlighter(QSyntaxHighlighter):
    """Highlighting: commands, strings, comments (*, //, ///, /* */), macros and numbers."""

    def __init__(self, doc: QTextDocument, theme: str = "light", palette: str = "classic"):
        super().__init__(doc)
        words = "|".join(re.escape(w) for w in sorted(_keywords(), key=len, reverse=True))
        self.re_cmd = re.compile(rf"\b(?:{words})\b")
        self.re_str = re.compile(r'`"(?:[^"]|"(?!\'))*"\'|"[^"\n]*"')
        self.re_macro = re.compile(r"`[A-Za-z_0-9:()\s.]*'|\$\{?[A-Za-z_][A-Za-z_0-9]*\}?")
        self.re_num = re.compile(r"(?<![\w.])(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?\b")
        self.set_theme(theme, palette)

    def set_theme(self, theme: str, palette: str = "classic") -> None:
        pal = editor_colors(palette, theme)
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


def toggle_comment(lines: list[str]) -> list[str]:
    """Comments the lines out with `// `; if all the non-blank ones already are, uncomments them."""
    body = [ln for ln in lines if ln.strip()]
    if body and all(ln.lstrip().startswith("//") for ln in body):
        out = []
        for ln in lines:
            s = ln.lstrip()
            if s.startswith("//"):
                ln = ln[: len(ln) - len(s)] + s[3 if s.startswith("// ") else 2:]
            out.append(ln)
        return out
    return [ln[: len(ln) - len(ln.lstrip())] + "// " + ln.lstrip() if ln.strip() else ln for ln in lines]


def indent_lines(lines: list[str]) -> list[str]:
    return ["\t" + ln if ln.strip() else ln for ln in lines]


def outdent_lines(lines: list[str]) -> list[str]:
    """Removes one leading tab, or up to four leading spaces."""
    out = []
    for ln in lines:
        if ln.startswith("\t"):
            ln = ln[1:]
        else:
            ln = ln[min(4, len(ln) - len(ln.lstrip(" "))):]
        out.append(ln)
    return out


class _Gutter(QWidget):
    def __init__(self, editor: "CodeEditor"):
        super().__init__(editor)
        self.editor = editor

    def sizeHint(self) -> QSize:
        return QSize(self.editor.gutter_width(), 0)

    def paintEvent(self, event) -> None:
        self.editor.paint_gutter(event)


class CodeEditor(QPlainTextEdit):
    def __init__(self, font: QFont, theme: str = "light", palette: str = "classic"):
        super().__init__()
        self.theme = theme
        self.palette_key = palette
        self.colors = editor_colors(palette, theme)
        self._apply_colors()
        self.path: Path | None = None
        self.setFont(font)
        self.document().setDocumentMargin(8)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.setTabStopDistance(self.fontMetrics().horizontalAdvance(" ") * 4)
        self.gutter = _Gutter(self)
        self.highlighter = StataHighlighter(self.document(), theme, palette)
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
        pal = self.colors
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
        pal = self.colors
        sel = QTextEdit.ExtraSelection()
        sel.format.setBackground(QColor(pal["current"]))
        sel.format.setProperty(QTextCharFormat.Property.FullWidthSelection, True)
        sel.cursor = self.textCursor()
        sel.cursor.clearSelection()
        self.setExtraSelections([sel])

    def _apply_colors(self) -> None:
        bg = self.colors.get("bg")
        self.setStyleSheet(f"QPlainTextEdit {{ background: {bg}; color: {self.colors['fg']}; }}" if bg else "")

    def set_theme(self, theme: str, palette: str | None = None) -> None:
        self.theme = theme
        if palette is not None:
            self.palette_key = palette
        self.colors = editor_colors(self.palette_key, theme)
        self._apply_colors()
        self.highlighter.set_theme(theme, self.palette_key)
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

    def current_line_code(self) -> str:
        """The line of the cursor, with the lines joined to it by `///`."""
        blk = self.textCursor().block()
        while blk.previous().isValid() and blk.previous().text().rstrip().endswith("///"):
            blk = blk.previous()
        lines = [blk.text()]
        while blk.text().rstrip().endswith("///") and blk.next().isValid():
            blk = blk.next()
            lines.append(blk.text())
        return "\n".join(lines)

    def transform_lines(self, fn: Callable[[list[str]], list[str]]) -> None:
        """Applies fn to the lines of the selection (or the cursor line), as a single undo step."""
        cur = self.textCursor()
        had_selection = cur.hasSelection()
        doc = self.document()
        a = doc.findBlock(cur.selectionStart())
        b = doc.findBlock(cur.selectionEnd())
        if had_selection and cur.selectionEnd() == b.position() and b.blockNumber() > a.blockNumber():
            b = b.previous()
        first, start = a.blockNumber(), a.position()
        cur.setPosition(start)
        cur.setPosition(b.position() + b.length() - 1, QTextCursor.MoveMode.KeepAnchor)
        new = fn(cur.selectedText().split("\u2029"))
        cur.beginEditBlock()
        cur.insertText("\n".join(new))
        cur.endEditBlock()
        if had_selection:           # keeps the lines selected, so the button can be pressed again
            last = doc.findBlockByNumber(first + len(new) - 1)
            cur.setPosition(start)
            cur.setPosition(last.position() + last.length() - 1, QTextCursor.MoveMode.KeepAnchor)
            self.setTextCursor(cur)


class _FindBar(QWidget):
    def __init__(self, owner: "DoFileEditor"):
        super().__init__()
        self.owner = owner
        self.setObjectName("FindBar")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(10, 6, 10, 6)
        lay.setSpacing(6)
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
    files_changed = Signal()        # files opened, closed or saved under another name

    def __init__(self, session: "Session", parent=None, *, run: Callable[[str], None],
                 font: QFont | None = None, theme: str = "light", palette: str = "classic"):
        super().__init__(parent)
        self.session = session
        self.run = run
        self.theme = theme
        self.palette = palette
        self._closing = False
        self.editor_font = font or QFont("Menlo", 12)
        self.resize(860, 640)
        self.untitled = 0
        self.tmpfiles: list[str] = []
        self.tmp_sources: dict[str, str] = {}     # temporary file -> what it holds, for Review

        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        self.tabs.currentChanged.connect(lambda _: (self._update_title(), self._update_undo()))
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
        tb.setObjectName("DoToolbar")
        tb.setMovable(False)
        tb.setIconSize(QSize(22, 22))
        tb.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
        key = QKeySequence.StandardKey
        edit = self._edit
        # (text, icon, slot, shortcut, bound?) - undo and redo keep the editor's own shortcuts
        groups = [
            [("New", "new", self.new_file, key.New, True), ("Open...", "open", self.open_dialog, key.Open, True),
             ("Save", "floppy", self.save, key.Save, True),
             ("Save as...", "floppy-saveas", self.save_as, key.SaveAs, True)],
            [("Undo", "undo", lambda: edit(lambda e: e.undo()), key.Undo, False),
             ("Redo", "redo", lambda: edit(lambda e: e.redo()), key.Redo, False)],
            [("Find", "search", self.show_find, key.Find, True)],
            [("Comment/uncomment", "comment", lambda: edit(lambda e: e.transform_lines(toggle_comment)),
              "Ctrl+/", True),
             ("Indent", "indent", lambda: edit(lambda e: e.transform_lines(indent_lines)), "Ctrl+]", True),
             ("Outdent", "outdent", lambda: edit(lambda e: e.transform_lines(outdent_lines)), "Ctrl+[", True)],
            # on the right, as in Stata 14, with Do at the far end
            [("Do line", "play-line", self.do_line, "Ctrl+Return", True),
             ("Run", "play-quiet", self.run_selection, "Ctrl+R", True),
             ("Do", "play", self.do_selection, "Ctrl+D", True)],
        ]
        tips = {"Do": "Execute the selection, or the whole file (do)",
                "Run": "Execute the selection, or the whole file, without showing commands (run)",
                "Do line": "Execute the current line (do)"}
        self.tool_actions: dict[str, QAction] = {}
        for gi, group in enumerate(groups):
            if gi == len(groups) - 1:
                spacer = QWidget()
                spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
                tb.addWidget(spacer)
            elif gi:
                tb.addSeparator()
            for text, name, slot, shortcut, bound in group:
                accent = RUN_GREEN if name.startswith("play") else None
                act = QAction(icon(name, accent), tr(text), self)
                act.setData((name, accent))
                act.setProperty("key", text)
                act.setProperty("label", _SHORT_LABELS.get(text, text))
                act.setProperty("tip", tips.get(text, text))
                act.setProperty("hint", QKeySequence(shortcut).toString(QKeySequence.SequenceFormat.NativeText))
                act.triggered.connect(slot)
                if bound:
                    act.setShortcut(QKeySequence(shortcut))
                self._retitle(act)
                tb.addAction(act)
                self.tool_actions[name] = act
        self.addToolBar(tb)
        self.toolbar = tb
        self._fit_toolbar_buttons()
        self.new_tab_button = QToolButton()          # the + at the end of the tab bar
        self.new_tab_button.setAutoRaise(True)
        self.new_tab_button.setIcon(icon("plus"))
        self.new_tab_button.setIconSize(QSize(16, 16))
        self.new_tab_button.clicked.connect(self.new_file)
        self.tabs.setCornerWidget(self.new_tab_button, Qt.Corner.TopRightCorner)
        self._retitle_tab_button()
        self.setStyleSheet("""
            QToolBar#DoToolbar { border: none; spacing: 2px; padding: 4px 10px; }
            QToolBar#DoToolbar QToolButton { border: none; border-radius: 8px; padding: 4px 4px; }
            QToolBar#DoToolbar QToolButton:hover { background: rgba(127, 127, 127, 0.16); }
            QToolBar#DoToolbar QToolButton:pressed { background: rgba(127, 127, 127, 0.28); }
            QToolBar#DoToolbar::separator { width: 1px; margin: 6px 8px; background: rgba(127, 127, 127, 0.3); }
            QTabWidget::pane { border: none; }
            QTabBar { background: transparent; }
            QTabBar::tab { background: transparent; border: none; border-radius: 8px; padding: 5px 12px;
                           margin: 4px 2px 4px 2px; }
            QTabBar::tab:selected { background: rgba(127, 127, 127, 0.22); font-weight: 600; }
            QTabBar::tab:hover:!selected { background: rgba(127, 127, 127, 0.12); }
            QTabBar::close-button { border-radius: 6px; subcontrol-position: right; }
            QTabBar::close-button:hover { background: rgba(127, 127, 127, 0.3); }
            QPlainTextEdit { border: none; }
            QWidget#FindBar { border-top: 1px solid rgba(127, 127, 127, 0.3); }
            QWidget#FindBar QLineEdit { border: 1px solid rgba(127, 127, 127, 0.4); border-radius: 8px;
                                        padding: 5px 10px; }
            QWidget#FindBar QLineEdit:focus { border: 1px solid #0a84ff; }
            QWidget#FindBar QPushButton { border: 1px solid rgba(127, 127, 127, 0.4); border-radius: 8px;
                                          padding: 5px 14px; }
            QWidget#FindBar QPushButton:hover { background: rgba(127, 127, 127, 0.16); }
        """)
        QShortcut(QKeySequence("Ctrl+Shift+D"), self, self.do_selection)
        QShortcut(QKeySequence("Escape"), self, self.findbar.hide)
        self.new_file()

    def _retitle(self, act: QAction) -> None:
        """Name and tooltip (with the shortcut) in the current language."""
        act.setText(tr(act.property("label") or act.property("key")))
        hint = act.property("hint")
        act.setToolTip(tr(act.property("tip")) + (f"  [{hint}]" if hint else ""))

    def _fit_toolbar_buttons(self) -> None:
        """Small font and a width that follows the label, as in the main window's toolbar."""
        font = QFont(self.font())
        font.setPixelSize(11)
        for act in self.tool_actions.values():
            btn = self.toolbar.widgetForAction(act)
            if isinstance(btn, QToolButton):
                btn.setFont(font)
                btn.setFixedWidth(max(40, btn.fontMetrics().horizontalAdvance(act.text()) + 16))

    def _retitle_tab_button(self) -> None:
        shortcut = QKeySequence(QKeySequence.StandardKey.New).toString(QKeySequence.SequenceFormat.NativeText)
        self.new_tab_button.setToolTip(f"{tr('New tab')}  [{shortcut}]")

    def _edit(self, fn: Callable[[CodeEditor], None]) -> None:
        ed = self.current()
        if ed is not None:
            fn(ed)

    def _update_undo(self) -> None:
        ed = self.current()
        doc = ed.document() if ed is not None else None
        self.tool_actions["undo"].setEnabled(doc is not None and doc.isUndoAvailable())
        self.tool_actions["redo"].setEnabled(doc is not None and doc.isRedoAvailable())

    # -- abas -----------------------------------------------------------------------------
    def current(self) -> CodeEditor | None:
        w = self.tabs.currentWidget()
        return w if isinstance(w, CodeEditor) else None

    def _add(self, ed: CodeEditor, title: str) -> None:
        idx = self.tabs.addTab(ed, title)
        self.tabs.setCurrentIndex(idx)
        ed.document().modificationChanged.connect(lambda _m, e=ed: self._refresh_tab(e))
        ed.undoAvailable.connect(lambda _a: self._update_undo())
        ed.redoAvailable.connect(lambda _a: self._update_undo())
        ed.setFocus()
        self._update_title()
        self._files_changed()

    def _files_changed(self) -> None:
        if not self._closing:
            self.files_changed.emit()

    def open_paths(self) -> list[str]:
        """Saved files of the open tabs, in order."""
        out = []
        for i in range(self.tabs.count()):
            ed = self.tabs.widget(i)
            if isinstance(ed, CodeEditor) and ed.path is not None:
                out.append(str(ed.path))
        return out

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
        ed = CodeEditor(self.editor_font, self.theme, self.palette)
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
        ed = CodeEditor(self.editor_font, self.theme, self.palette)
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
        ok = self.save()
        self._files_changed()
        return ok

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
        self._files_changed()
        if self.tabs.count() == 0:
            self.hide()
            self.new_file()
        return True

    def closeEvent(self, e) -> None:
        self._closing = True        # closing the window keeps the list of files for the next session
        try:
            for i in reversed(range(self.tabs.count())):
                ed = self.tabs.widget(i)
                if isinstance(ed, CodeEditor) and ed.document().isModified():
                    if not self.close_tab(i):
                        e.ignore()
                        return
            e.accept()
        finally:
            self._closing = False

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
    def _execute(self, verb: str, line: bool = False) -> None:
        ed = self.current()
        if ed is None:
            return
        code = ed.current_line_code() if line else ed.code_to_run()
        if not code.strip():
            return
        fd, path = tempfile.mkstemp(prefix="SD", suffix=".tmp")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(code if code.endswith("\n") else code + "\n")
        self.tmpfiles.append(path)
        name = ed.path.name if ed.path else ed.property("untitled")
        if line:
            self.tmp_sources[path] = f"{tr('line')}: {name}"
        elif ed.textCursor().hasSelection():
            self.tmp_sources[path] = f"{tr('selection')}: {name}"
        else:
            self.tmp_sources[path] = name
        self.run(f'{verb} "{path}"')
        self.ran.emit()

    def do_selection(self) -> None:
        self._execute("do")

    def run_selection(self) -> None:
        self._execute("run")

    def do_line(self) -> None:
        self._execute("do", line=True)

    def set_theme(self, theme: str, palette: str | None = None) -> None:
        self.theme = theme
        if palette is not None:
            self.palette = palette
        for i in range(self.tabs.count()):
            ed = self.tabs.widget(i)
            if isinstance(ed, CodeEditor):
                ed.set_theme(theme, self.palette)
        for act in self.tool_actions.values():
            name, accent = act.data()
            act.setIcon(icon(name, accent))
        self.new_tab_button.setIcon(icon("plus"))

    def cleanup(self) -> None:
        for p in self.tmpfiles:
            try:
                os.unlink(p)
            except OSError:
                pass
        self.tmpfiles.clear()

    def retranslate(self) -> None:
        for act in self.tool_actions.values():
            self._retitle(act)
        self._retitle_tab_button()
        self._fit_toolbar_buttons()
        self._update_title()
