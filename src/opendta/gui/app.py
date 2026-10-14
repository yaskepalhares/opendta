"""Ponto de entrada da interface gráfica."""

from __future__ import annotations

import sys


def open_paths(win, paths) -> None:
    """Do-files são executados (como o Stata faz com `stata arquivo.do`);
    bases .dta são carregadas; outros arquivos de texto abrem no Do-file Editor."""
    for a in paths:
        low = str(a).lower()
        if low.endswith(".do"):
            win.run_command(f'do "{a}"')
        elif low.endswith(".dta"):
            win.open_dataset(str(a))
        elif low.endswith(".ado"):
            win.show_dofile_editor(str(a))


def run_gui(args: list[str] | None = None) -> int:
    from PySide6.QtCore import QEvent, QObject
    from PySide6.QtWidgets import QApplication

    from .. import i18n
    from .main_window import MainWindow
    from .settings import Preferences
    from .theme import apply_theme

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("OpenDTA")
    prefs = Preferences()
    i18n.set_language(prefs.language)
    apply_theme(app, prefs.theme)   # antes de criar a janela: sem piscar
    win = MainWindow()

    # macOS: arquivos abertos pelo Finder (duplo clique, "Abrir com", arrastar
    # para o ícone do Dock) chegam como QFileOpenEvent, não como argumentos
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
    return app.exec()
