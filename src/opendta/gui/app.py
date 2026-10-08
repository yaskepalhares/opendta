"""Ponto de entrada da interface gráfica."""

from __future__ import annotations

import sys


def run_gui(args: list[str] | None = None) -> int:
    from PySide6.QtWidgets import QApplication

    from .main_window import MainWindow

    from .settings import Preferences
    from .theme import apply_theme

    app = QApplication.instance() or QApplication(sys.argv[:1])
    app.setApplicationName("OpenDTA")
    apply_theme(app, Preferences().theme)   # antes de criar a janela: sem piscar
    win = MainWindow()
    win.show()
    # `opendta arquivo.do` abre a janela e executa o do-file, como o Stata faz
    for a in args or []:
        if a.lower().endswith(".do"):
            win.run_command(f'do "{a}"')
    return app.exec()
