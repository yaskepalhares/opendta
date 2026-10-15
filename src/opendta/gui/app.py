"""Entry point of the graphical interface."""

from __future__ import annotations

import sys


def open_paths(win, paths) -> None:
    """Do-files are run (as Stata does with `stata file.do`);
    .dta datasets are loaded; other text files open in the Do-file Editor."""
    for a in paths:
        low = str(a).lower()
        if low.endswith(".do"):
            win.run_command(f'do "{a}"')
        elif low.endswith(".dta"):
            win.open_dataset(str(a))
        elif low.endswith(".ado"):
            win.show_dofile_editor(str(a))


def run_gui(args: list[str] | None = None) -> int:
    from PySide6.QtCore import QEvent, QObject, QTimer
    from PySide6.QtWidgets import QApplication

    from .. import i18n
    from .main_window import MainWindow
    from .settings import Preferences
    from .theme import apply_theme

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("OpenDTA")
    prefs = Preferences()
    i18n.set_language(prefs.language)
    apply_theme(app, prefs.theme)   # before creating the window: no flicker ("system" is resolved inside)
    win = MainWindow()

    # macOS: files opened from the Finder (double click, "Open With", dragging onto the Dock icon)
    # arrive as QFileOpenEvent, not as arguments
    class _FileOpen(QObject):
        def eventFilter(self, obj, event):  # noqa: N802
            if event.type() == QEvent.Type.FileOpen:
                open_paths(win, [event.file()])
                return True
            return False

    opener = _FileOpen(app)
    app.installEventFilter(opener)
    win.show()
    open_paths(win, args or [])
    QTimer.singleShot(4000, win.updates.check_at_startup)   # daily check for a new release
    return app.exec()
