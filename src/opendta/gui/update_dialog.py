"""Update notice, download and installation (Help > Check for Updates and
the daily check at startup). The network work runs in a QThread; the
logic is in opendta.updater."""

from __future__ import annotations

import tempfile
import time
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal, Slot
from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import (QDialog, QFrame, QHBoxLayout, QLabel, QMessageBox, QProgressBar,
                               QPushButton, QTextBrowser, QVBoxLayout)

from .. import __version__, i18n, updater
from ..i18n import L, tr
from .icons import app_icon, is_dark
from .preferences import COLORS, _button, stylesheet

DAY = 24 * 3600


class _Worker(QObject):
    """Runs in its own QThread. Its signals must go to slots of objects that live in the GUI
    thread (never to lambdas, which would run here)."""
    done = Signal(object)
    failed = Signal(str)
    progress = Signal(int, int)

    def __init__(self, job, *args):
        super().__init__()
        self.job, self.args = job, args
        self.cancel = False

    def run(self) -> None:
        try:
            if self.job == "check":
                self.done.emit(updater.fetch_latest())
            else:
                asset, dest = self.args
                path = updater.download(asset, dest, progress=lambda a, b: self.progress.emit(a, b),
                                        cancelled=lambda: self.cancel)
                self.done.emit(path)
        except updater.UpdateError as e:
            self.failed.emit(str(e))
        except Exception as e:      # noqa: BLE001
            self.failed.emit(str(e))


def _start(owner, worker: _Worker) -> QThread:
    th = QThread(owner)
    worker.moveToThread(th)
    th.started.connect(worker.run)
    worker.done.connect(th.quit)
    worker.failed.connect(th.quit)
    th.finished.connect(worker.deleteLater)
    th.start()
    return th


class UpdateController(QObject):
    """Owned by the main window: checks, shows the dialog and keeps the install step."""
    checked = Signal()          # after each check, successful or not (Settings shows the date)

    def __init__(self, window):
        super().__init__(window)
        self.win = window
        self.prefs = window.prefs
        self._thread = None
        self._worker = None
        self._manual = False
        self._orphans: list[tuple[QThread, _Worker]] = []   # download threads still finishing

    # -- threads that outlive the dialog -------------------------------------------
    def adopt(self, thread: QThread, worker: _Worker) -> None:
        """The dialog can close (Cancel, Esc, success) while its QThread is still running. A
        QThread destroyed while running aborts the whole application, so the controller, which
        lives as long as the main window, takes the thread over until it finishes."""
        thread.setParent(self)
        self._orphans.append((thread, worker))
        thread.finished.connect(self._reap)
        if thread.isFinished():
            self._reap()

    @Slot()
    def _reap(self) -> None:
        alive = []
        for th, w in self._orphans:
            if th.isRunning():
                alive.append((th, w))
            else:
                th.wait()
                th.deleteLater()
        self._orphans = alive

    def shutdown(self) -> None:
        """Called when the main window closes: cancels downloads and waits for the threads."""
        for th, w in self._orphans + ([(self._thread, self._worker)] if self._thread is not None else []):
            try:
                w.cancel = True
                th.quit()
                th.wait(5000)
            except RuntimeError:        # already deleted
                pass
        self._orphans = []

    # -- checking --------------------------------------------------------------
    def check_at_startup(self) -> None:
        if not self.prefs.check_updates or updater.install_kind() == "source":
            return
        if time.time() - self.prefs.last_update_check < DAY:
            return
        self._check(manual=False)

    def check_now(self) -> None:
        self._check(manual=True)

    def _check(self, manual: bool) -> None:
        if self._thread is not None and self._thread.isRunning():
            return
        self._manual = manual
        self._worker = _Worker("check")
        # bound slots of this object (which lives in the GUI thread): Qt queues the call to the
        # GUI thread. A lambda would run in the worker thread, and a dialog opened there aborts
        # the application (macOS: NSWindow outside the main thread)
        self._worker.done.connect(self._on_checked)
        self._worker.failed.connect(self._on_check_failed)
        self._thread = _start(self, self._worker)

    @Slot(object)
    def _on_checked(self, rel) -> None:
        self._checked(rel, self._manual)

    @Slot(str)
    def _on_check_failed(self, msg: str) -> None:
        self._check_failed(msg, self._manual)

    def _check_failed(self, msg: str, manual: bool) -> None:
        self.checked.emit()
        if manual:
            QMessageBox.warning(self.win, tr("Check for Updates"),
                                L("Could not check for updates:", "Não foi possível procurar atualizações:")
                                + f"\n{msg}")

    def _checked(self, rel: updater.Release, manual: bool) -> None:
        self.prefs.last_update_check = time.time()
        self.prefs.sync()
        self.checked.emit()
        if not updater.is_newer(rel.version):
            if manual:
                QMessageBox.information(self.win, tr("Check for Updates"),
                                        L(f"OpenDTA {__version__} is the latest version.",
                                          f"O OpenDTA {__version__} é a versão mais recente."))
            return
        if not manual and self.prefs.skipped_version == rel.version:
            return
        UpdateDialog(self, rel, self.win).exec()

    # -- installing (called by the dialog) ----------------------------------------
    def install_after_close(self, start) -> None:
        """Closes OpenDTA (with the usual questions about unsaved work) and starts the installation."""
        self.win._pending_update = start
        if not self.win.close():
            self.win._pending_update = None
            QMessageBox.information(self.win, tr("Update"),
                                    L("The update was not installed because OpenDTA was not closed. "
                                      "It will be offered again next time.",
                                      "A atualização não foi instalada porque o OpenDTA não foi fechado. "
                                      "Ela será oferecida de novo da próxima vez."))


def _style(c: dict) -> str:
    """The Settings window's stylesheet plus what only this dialog has: the notes box, the
    filled main button and the progress bar."""
    return stylesheet(c) + f"""
QTextBrowser#Notes {{ background: transparent; border: none; color: {c['text']}; font-size: 13px;
                      selection-background-color: {c['accent']}; selection-color: #ffffff; }}
QTextBrowser#Notes QScrollBar:vertical {{ background: transparent; width: 10px; margin: 2px; }}
QTextBrowser#Notes QScrollBar::handle:vertical {{ background: {c['border']}; border-radius: 3px;
                                                 min-height: 28px; }}
QTextBrowser#Notes QScrollBar::add-line:vertical, QTextBrowser#Notes QScrollBar::sub-line:vertical {{ height: 0; }}
QTextBrowser#Notes QScrollBar::add-page:vertical, QTextBrowser#Notes QScrollBar::sub-page:vertical {{ background: transparent; }}
QPushButton#Primary {{ background: {c['accent']}; border: 1px solid {c['accent']}; border-radius: 8px;
                       padding: 6px 16px; color: #ffffff; font-size: 13px; font-weight: 600; }}
QPushButton#Primary:hover {{ background: #3b9bff; border-color: #3b9bff; }}
QPushButton#Primary:disabled {{ background: {c['off']}; border-color: {c['off']}; }}
QPushButton#Secondary:disabled, QPushButton#Link:disabled {{ color: {c['off']}; }}
QProgressBar {{ background: {c['control']}; border: none; border-radius: 3px; min-height: 6px;
                max-height: 6px; text-align: center; color: transparent; }}
QProgressBar::chunk {{ background: {c['accent']}; border-radius: 3px; }}
"""


class UpdateDialog(QDialog):
    """Same look as the Settings window. The release notes show only the part in the interface
    language."""

    def __init__(self, ctl: UpdateController, rel: updater.Release, parent=None):
        super().__init__(parent)
        self.ctl, self.rel = ctl, rel
        self.kind = updater.install_kind()
        self.asset = updater.asset_for(rel, self.kind)
        self.setWindowTitle(tr("Update Available"))
        self.setObjectName("SettingsRoot")
        self.resize(600, 560)
        self.setMinimumSize(480, 420)
        c = COLORS["dark" if is_dark() else "light"]
        self.setStyleSheet(_style(c))
        lay = QVBoxLayout(self)
        lay.setContentsMargins(28, 24, 28, 20)
        lay.setSpacing(0)
        head = QHBoxLayout()
        head.setSpacing(16)
        tv = QVBoxLayout()
        tv.setSpacing(2)
        title = QLabel(L(f"OpenDTA {rel.version} is available", f"O OpenDTA {rel.version} está disponível"))
        title.setObjectName("BigTitle")
        have = QLabel(L(f"You have version {__version__}.", f"Você tem a versão {__version__}."))
        have.setObjectName("Muted")
        tv.addWidget(title)
        tv.addWidget(have)
        logo = QLabel()
        logo.setPixmap(app_icon(ctl.prefs.app_icon).pixmap(56, 56))
        head.addWidget(logo)
        head.addLayout(tv, 1)
        lay.addLayout(head)
        lay.addSpacing(18)
        box = QFrame()
        box.setObjectName("Preview")
        bl = QVBoxLayout(box)
        bl.setContentsMargins(6, 6, 6, 6)
        notes = QTextBrowser()
        notes.setObjectName("Notes")
        notes.setOpenExternalLinks(True)
        notes.document().setDocumentMargin(10)
        pal = notes.palette()
        pal.setColor(QPalette.ColorRole.Link, QColor(c["accent"]))
        notes.setPalette(pal)
        text = updater.notes_for(rel.notes, i18n.language())
        notes.setMarkdown(text or L("(no release notes)", "(sem notas da versão)"))
        bl.addWidget(notes)
        lay.addWidget(box, 1)
        lay.addSpacing(14)
        self.status = QLabel("")
        self.status.setObjectName("Muted")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)
        self.bar = QProgressBar()
        self.bar.setTextVisible(False)
        self.bar.hide()
        lay.addWidget(self.bar)
        lay.addSpacing(14)
        row = QHBoxLayout()
        row.setSpacing(10)
        self.skip = _button(tr("Skip This Version"), "Link")
        self.later = _button(tr("Remind Me Later"))
        if self.kind == "source" or self.asset is None:
            self.install = _button(tr("Open Download Page"), "Primary")
            if self.kind == "source":
                self.status.setText(L("OpenDTA is running from source: update it with git pull and "
                                      "pip install -e .",
                                      "O OpenDTA está rodando a partir do código-fonte: atualize com git pull "
                                      "e pip install -e ."))
        else:
            self.install = _button(tr("Install and Restart"), "Primary")
            self.status.setText(L("OpenDTA will download the update, close and open again.",
                                  "O OpenDTA vai baixar a atualização, fechar e abrir de novo."))
        self.install.setDefault(True)
        row.addWidget(self.skip)
        row.addStretch(1)
        row.addWidget(self.later)
        row.addWidget(self.install)
        lay.addLayout(row)
        self.skip.clicked.connect(self._skip)
        self.later.clicked.connect(self.reject)
        self.install.clicked.connect(self._install)
        self._worker = None
        self._thread = None

    def done(self, r: int) -> None:
        """Every way of closing the dialog (buttons, Esc, window close) ends here."""
        self._release_download()
        super().done(r)

    def _release_download(self) -> None:
        th, w = self._thread, self._worker
        self._thread = None
        if th is None or w is None:
            return
        w.cancel = True
        for sig, slot in ((w.progress, self._progress), (w.done, self._downloaded),
                          (w.failed, self._failed)):
            try:
                sig.disconnect(slot)
            except (RuntimeError, TypeError):
                pass
        self.ctl.adopt(th, w)

    def _skip(self) -> None:
        self.ctl.prefs.skipped_version = self.rel.version
        self.ctl.prefs.sync()
        self.reject()

    def _install(self) -> None:
        if self.kind == "source" or self.asset is None:
            updater.open_externally(self.rel.page)
            self.accept()
            return
        self.install.setEnabled(False)
        self.skip.setEnabled(False)
        self.later.setText(tr("Cancel"))
        self.later.clicked.disconnect()
        self.later.clicked.connect(self._cancel)
        self.bar.show()
        self.bar.setRange(0, 0)
        self.status.setText(L("Downloading...", "Baixando..."))
        dest = Path(tempfile.gettempdir()) / "opendta-update"
        self._worker = _Worker("download", self.asset, dest)
        self._worker.progress.connect(self._progress)
        self._worker.done.connect(self._downloaded)
        self._worker.failed.connect(self._failed)
        self._thread = _start(self, self._worker)

    @Slot(int, int)
    def _progress(self, done: int, total: int) -> None:
        if total:
            self.bar.setRange(0, 1000)
            self.bar.setValue(int(1000 * done / total))
            self.status.setText(L(f"Downloading... {done / 1e6:.1f} of {total / 1e6:.1f} MB",
                                  f"Baixando... {done / 1e6:.1f} de {total / 1e6:.1f} MB"))

    def _cancel(self) -> None:
        self.reject()       # done() cancels the download and hands the thread to the controller

    @Slot(str)
    def _failed(self, msg: str) -> None:
        if msg == "cancelled":
            return
        self.bar.hide()
        self.status.setText(L("The download failed: ", "O download falhou: ") + msg)
        self.install.setEnabled(True)
        self.install.setText(tr("Open Download Page"))
        self.install.clicked.disconnect()
        self.install.clicked.connect(lambda: (updater.open_externally(self.rel.page), self.accept()))

    @Slot(object)
    def _downloaded(self, path) -> None:
        self.bar.setRange(0, 1)
        self.bar.setValue(1)
        try:
            start = updater.plan_install(self.kind, Path(path))
        except updater.UpdateError as e:
            # cannot replace in place (e.g. folder not writable): open the downloaded file
            self.status.setText(L(f"OpenDTA cannot update itself here ({e}). The downloaded file will "
                                  "open so that you can install it.",
                                  f"O OpenDTA não consegue se atualizar aqui ({e}). O arquivo baixado vai "
                                  "abrir para você instalar."))
            updater.open_externally(path)
            self.accept()
            return
        self.accept()
        self.ctl.install_after_close(start)
