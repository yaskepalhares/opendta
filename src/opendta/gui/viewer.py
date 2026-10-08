"""Viewer: mostra páginas de help e arquivos SMCL ou texto.

Os links {help tópico} abrem outra página no mesmo Viewer; links da web
abrem no navegador; Voltar e Avançar percorrem o histórico. A caixa no
topo aceita um tópico ("describe"), "help tópico" ou o caminho de um
arquivo, como a caixa de comando do Viewer do Stata.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Callable

from PySide6.QtCore import QUrl
from PySide6.QtGui import QAction, QDesktopServices, QKeySequence, QShortcut
from PySide6.QtWidgets import QLineEdit, QMainWindow, QTextBrowser, QToolBar

from ..core.smcl import render_html

if TYPE_CHECKING:
    from ..session import Session


class Viewer(QMainWindow):
    def __init__(self, session: "Session", parent=None, *, font_family: str = "Menlo",
                 font_size: int = 12, dark: bool = False,
                 run: Callable[[str], int] | None = None):
        super().__init__(parent)
        self.session = session
        self.run = run
        self.font_family, self.font_size, self.dark = font_family, font_size, dark
        self.resize(820, 640)
        self.history: list[tuple[Path, str]] = []
        self.pos = -1

        tb = QToolBar("Viewer")
        tb.setMovable(False)
        self.act_back = QAction("◀", self)
        self.act_back.setToolTip("Back")
        self.act_back.triggered.connect(lambda: self.go(-1))
        self.act_fwd = QAction("▶", self)
        self.act_fwd.setToolTip("Forward")
        self.act_fwd.triggered.connect(lambda: self.go(1))
        tb.addAction(self.act_back)
        tb.addAction(self.act_fwd)
        self.field = QLineEdit()
        self.field.setPlaceholderText("help topic, or a file to view")
        self.field.returnPressed.connect(self._from_field)
        tb.addWidget(self.field)
        self.addToolBar(tb)

        self.browser = QTextBrowser()
        self.browser.setOpenLinks(False)
        self.browser.anchorClicked.connect(self._link)
        self.setCentralWidget(self.browser)
        QShortcut(QKeySequence.StandardKey.Close, self, self.close)
        QShortcut(QKeySequence.StandardKey.Back, self, lambda: self.go(-1))
        QShortcut(QKeySequence.StandardKey.Forward, self, lambda: self.go(1))

    # -- navegação -----------------------------------------------------------------
    def show_page(self, path: Path, title: str, *, record: bool = True) -> None:
        raw = path.read_bytes()
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            text = raw.decode("latin-1")
        if path.suffix.lower() in (".smcl", ".sthlp", ".hlp") or text.lstrip().startswith("{smcl}"):
            html = render_html(text, dark=self.dark, font=self.font_family, size=self.font_size)
            self.browser.setHtml(html)
        else:
            self.browser.setPlainText(text)
        self.setWindowTitle(f"Viewer — {title}")
        self.field.setText(title if path.suffix == ".sthlp" else str(path))
        if record:
            self.history = self.history[: self.pos + 1] + [(path, title)]
            self.pos = len(self.history) - 1
        self._update_buttons()

    def go(self, step: int) -> None:
        k = self.pos + step
        if 0 <= k < len(self.history):
            self.pos = k
            path, title = self.history[k]
            self.show_page(path, title, record=False)

    def _update_buttons(self) -> None:
        self.act_back.setEnabled(self.pos > 0)
        self.act_fwd.setEnabled(self.pos < len(self.history) - 1)

    def open_topic(self, topic: str) -> bool:
        from ..commands.helpcmd import find_help
        path = find_help(self.session, topic)
        if path is None:
            self.browser.setHtml(f"<p>help for <b>{topic}</b> not found</p>"
                                 '<p>Type a topic above, or see <a href="help:opendta">the index</a>.</p>')
            return False
        self.show_page(path, topic or "opendta")
        return True

    def _from_field(self) -> None:
        text = self.field.text().strip()
        if text.startswith("help "):
            text = text[5:].strip()
        p = Path(text).expanduser()
        if text and p.exists():
            self.show_page(p, p.name)
        else:
            self.open_topic(text)

    def _link(self, url: QUrl) -> None:
        href = url.toString()
        if href.startswith("help:"):
            self.open_topic(href[5:])
        elif href.startswith("view:"):
            p = Path(href[5:]).expanduser()
            if p.exists():
                self.show_page(p, p.name)
        elif href.startswith("stata:") and self.run is not None:
            self.run(href[6:])
        else:
            QDesktopServices.openUrl(url)
