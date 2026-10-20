"""Qt's own texts in the language of the interface.

OpenDTA translates its menus and windows (i18n.py), but the texts Qt itself provides (the standard
buttons Save / Don't Save / Cancel / Yes / No / OK, the right-click menu of text fields: Undo, Cut,
Copy, Paste, Select All, and the print and color dialogs) come from Qt's translation files. The
PySide6 package ships them; packaging/opendta.spec copies the Portuguese ones into the app.
"""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import QCoreApplication, QLibraryInfo, QTranslator

_installed: list[QTranslator] = []


def _folders():
    yield Path(QLibraryInfo.path(QLibraryInfo.LibraryPath.TranslationsPath))
    import PySide6
    yield Path(PySide6.__file__).parent / "Qt" / "translations"
    bundle = getattr(sys, "_MEIPASS", None)         # the packaged application
    if bundle:
        yield Path(bundle) / "PySide6" / "Qt" / "translations"


def apply(language: str) -> None:
    """Loads Qt's translations for the interface language ("pt" = Brazilian Portuguese; English is
    Qt's own, so nothing is loaded)."""
    app = QCoreApplication.instance()
    if app is None:
        return
    for translator in _installed:
        app.removeTranslator(translator)
    _installed.clear()
    if language != "pt":
        return
    for base in ("qtbase", "qt"):
        for folder in _folders():
            translator = QTranslator()
            if translator.load(f"{base}_pt_BR", str(folder)):
                app.installTranslator(translator)
                _installed.append(translator)      # kept alive here
                break
