"""Update notice, download and installation (Help > Check for Updates and
the daily check at startup). The network work runs in a QThread; the
logic is in opendta.updater."""

from __future__ import annotations

import tempfile
import time
from pathlib import Path

from PySide6.QtCore import QObject, QThread, Signal
from PySide6.QtWidgets import (QDialog, QHBoxLayout, QLabel, QMessageBox, QProgressBar,
                               QPushButton, QTextBrowser, QVBoxLayout)

from .. import __version__, updater
from ..i18n import L, tr

DAY = 24 * 3600


class _Worker(QObject):
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

    def __init__(self, window):
        super().__init__(window)
        self.win = window
        self.prefs = window.prefs
        self._thread = None
        self._worker = None

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
        self._worker = _Worker("check")
        self._worker.done.connect(lambda rel: self._checked(rel, manual))
        self._worker.failed.connect(lambda msg: self._check_failed(msg, manual))
        self._thread = _start(self, self._worker)

    def _check_failed(self, msg: str, manual: bool) -> None:
        if manual:
            QMessageBox.warning(self.win, tr("Check for Updates"),
                                L("Could not check for updates:", "Não foi possível procurar atualizações:")
                                + f"\n{msg}")

    def _checked(self, rel: updater.Release, manual: bool) -> None:
        self.prefs.last_update_check = time.time()
        self.prefs.sync()
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


class UpdateDialog(QDialog):
    def __init__(self, ctl: UpdateController, rel: updater.Release, parent=None):
        super().__init__(parent)
        self.ctl, self.rel = ctl, rel
        self.kind = updater.install_kind()
        self.asset = updater.asset_for(rel, self.kind)
        self.setWindowTitle(tr("Update Available"))
        self.resize(560, 440)
        lay = QVBoxLayout(self)
        head = QLabel(L(f"<b>OpenDTA {rel.version}</b> is available. You have version {__version__}.",
                        f"O <b>OpenDTA {rel.version}</b> está disponível. Você tem a versão {__version__}."))
        lay.addWidget(head)
        notes = QTextBrowser()
        notes.setOpenExternalLinks(True)
        notes.setMarkdown(rel.notes or L("(no release notes)", "(sem notas da versão)"))
        lay.addWidget(notes, 1)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        lay.addWidget(self.status)
        self.bar = QProgressBar()
        self.bar.hide()
        lay.addWidget(self.bar)
        row = QHBoxLayout()
        self.skip = QPushButton(tr("Skip This Version"))
        self.later = QPushButton(tr("Remind Me Later"))
        if self.kind == "source" or self.asset is None:
            self.install = QPushButton(tr("Open Download Page"))
            if self.kind == "source":
                self.status.setText(L("OpenDTA is running from source: update it with git pull and "
                                      "pip install -e .",
                                      "O OpenDTA está rodando a partir do código-fonte: atualize com git pull "
                                      "e pip install -e ."))
        else:
            self.install = QPushButton(tr("Install and Restart"))
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

    def _progress(self, done: int, total: int) -> None:
        if total:
            self.bar.setRange(0, 1000)
            self.bar.setValue(int(1000 * done / total))
            self.status.setText(L(f"Downloading... {done / 1e6:.1f} of {total / 1e6:.1f} MB",
                                  f"Baixando... {done / 1e6:.1f} de {total / 1e6:.1f} MB"))

    def _cancel(self) -> None:
        if self._worker is not None:
            self._worker.cancel = True
        self.reject()

    def _failed(self, msg: str) -> None:
        if msg == "cancelled":
            return
        self.bar.hide()
        self.status.setText(L("The download failed: ", "O download falhou: ") + msg)
        self.install.setEnabled(True)
        self.install.setText(tr("Open Download Page"))
        self.install.clicked.disconnect()
        self.install.clicked.connect(lambda: (updater.open_externally(self.rel.page), self.accept()))

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
