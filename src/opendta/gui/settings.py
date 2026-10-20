"""Persistent interface preferences (QSettings).

They live in each system's standard place: the registry on Windows, a .plist on macOS
and ~/.config on Linux. With the OPENDTA_SETTINGS environment variable pointing
to an .ini file, they live in that file (the tests use this so they neither
read nor change the preferences of whoever runs pytest).
"""

from __future__ import annotations

import json
import os

from PySide6.QtCore import QSettings

from .palettes import KEYS as PALETTE_KEYS
from .theme import default_monospace

ORG, APP = "OpenDTA", "OpenDTA"

APP_ICON_VARIANTS = {
    "light": "Table — light (default)",
    "dark": "Table — dark",
}


def default_store() -> QSettings:
    path = os.environ.get("OPENDTA_SETTINGS")
    if path:
        return QSettings(path, QSettings.Format.IniFormat)
    return QSettings(ORG, APP)


class Preferences:
    def __init__(self, settings: QSettings | None = None):
        self._s = settings or default_store()

    # application icon
    @property
    def app_icon(self) -> str:
        v = str(self._s.value("appearance/app_icon", "light"))
        return v if v in APP_ICON_VARIANTS else "light"

    @app_icon.setter
    def app_icon(self, value: str) -> None:
        self._s.setValue("appearance/app_icon", value if value in APP_ICON_VARIANTS else "light")

    # interface theme: light (default), dark, or "system" (follows the operating system)
    @property
    def theme(self) -> str:
        v = str(self._s.value("appearance/theme", "light"))
        return v if v in ("light", "dark", "system") else "light"

    @theme.setter
    def theme(self, value: str) -> None:
        self._s.setValue("appearance/theme", value if value in ("light", "dark", "system") else "light")

    @property
    def effective_theme(self) -> str:
        """light or dark: the chosen theme, with "system" resolved from the operating system."""
        from .theme import system_theme
        return system_theme() if self.theme == "system" else self.theme

    # font of the Results, Command and Review windows
    @property
    def hints(self) -> bool:
        """Error explanations (set hints) turned on when the program opens."""
        return str(self._s.value("behavior/hints", "on")) != "off"

    @hints.setter
    def hints(self, on: bool) -> None:
        self._s.setValue("behavior/hints", "on" if on else "off")

    @property
    def check_updates(self) -> bool:
        """Look for a new release on GitHub when OpenDTA opens (at most once a day)."""
        return str(self._s.value("updates/check", "on")) != "off"

    @check_updates.setter
    def check_updates(self, on: bool) -> None:
        self._s.setValue("updates/check", "on" if on else "off")

    @property
    def last_update_check(self) -> float:
        try:
            return float(self._s.value("updates/last", 0) or 0)
        except (TypeError, ValueError):
            return 0.0

    @last_update_check.setter
    def last_update_check(self, t: float) -> None:
        self._s.setValue("updates/last", float(t))

    @property
    def skipped_version(self) -> str:
        return str(self._s.value("updates/skip", "") or "")

    @skipped_version.setter
    def skipped_version(self, v: str) -> None:
        self._s.setValue("updates/skip", v)

    @property
    def superscript(self) -> bool:
        """Superscript exponent in scientific notation (set superscript) at startup."""
        return str(self._s.value("behavior/superscript", "on")) != "off"

    @superscript.setter
    def superscript(self, on: bool) -> None:
        self._s.setValue("behavior/superscript", "on" if on else "off")

    @property
    def numerics(self) -> str:
        """set numerics stata|precise at startup."""
        v = str(self._s.value("behavior/numerics", "stata"))
        return v if v in ("stata", "precise") else "stata"

    @numerics.setter
    def numerics(self, value: str) -> None:
        self._s.setValue("behavior/numerics", value if value in ("stata", "precise") else "stata")

    # interface language (phase 9)
    @property
    def language(self) -> str:
        v = str(self._s.value("behavior/language", "en"))
        return v if v in ("en", "pt") else "en"

    @language.setter
    def language(self, value: str) -> None:
        self._s.setValue("behavior/language", value if value in ("en", "pt") else "en")

    @property
    def messages(self) -> str:
        """stata (error messages in English, as in Stata) or translated."""
        v = str(self._s.value("behavior/messages", "stata"))
        return v if v in ("stata", "translated") else "stata"

    @messages.setter
    def messages(self, value: str) -> None:
        self._s.setValue("behavior/messages", value if value in ("stata", "translated") else "stata")

    @property
    def font_family(self) -> str:
        return str(self._s.value("fonts/family", default_monospace()[0]))

    @font_family.setter
    def font_family(self, value: str) -> None:
        self._s.setValue("fonts/family", value)

    @property
    def font_size(self) -> int:
        try:
            return int(self._s.value("fonts/size", default_monospace()[1]))
        except (TypeError, ValueError):
            return default_monospace()[1]

    @font_size.setter
    def font_size(self, value: int) -> None:
        self._s.setValue("fonts/size", int(value))

    # color schemes (palettes.py) of the Results window and of the Do-file Editor
    @property
    def results_palette(self) -> str:
        v = str(self._s.value("colors/results", "classic"))
        return v if v in PALETTE_KEYS else "classic"

    @results_palette.setter
    def results_palette(self, value: str) -> None:
        self._s.setValue("colors/results", value if value in PALETTE_KEYS else "classic")

    @property
    def editor_palette(self) -> str:
        v = str(self._s.value("colors/editor", "classic"))
        return v if v in PALETTE_KEYS else "classic"

    @editor_palette.setter
    def editor_palette(self, value: str) -> None:
        self._s.setValue("colors/editor", value if value in PALETTE_KEYS else "classic")

    # Do-file Editor: reopen the files of the last session (off by default)
    @property
    def reopen_files(self) -> bool:
        return str(self._s.value("editor/reopen", "off")) == "on"

    @reopen_files.setter
    def reopen_files(self, on: bool) -> None:
        self._s.setValue("editor/reopen", "on" if on else "off")

    @property
    def editor_files(self) -> list[str]:
        try:
            files = json.loads(str(self._s.value("editor/files", "[]")))
        except ValueError:
            return []
        return [f for f in files if isinstance(f, str)] if isinstance(files, list) else []

    @editor_files.setter
    def editor_files(self, files: list[str]) -> None:
        self._s.setValue("editor/files", json.dumps(list(files)))

    # session: initial folder and automatic log (both applied when OpenDTA opens)
    @property
    def start_dir(self) -> str:
        return str(self._s.value("session/start_dir", "") or "")

    @start_dir.setter
    def start_dir(self, path: str) -> None:
        self._s.setValue("session/start_dir", path)

    @property
    def auto_log(self) -> bool:
        return str(self._s.value("session/auto_log", "off")) == "on"

    @auto_log.setter
    def auto_log(self, on: bool) -> None:
        self._s.setValue("session/auto_log", "on" if on else "off")

    # confirmation before risky commands (guard.py): one switch per command, all on by default;
    # and whether the first-time notice was shown
    def confirm_enabled(self, command: str) -> bool:
        return str(self._s.value(f"safety/confirm_{command}", "on")) != "off"

    def set_confirm(self, command: str, on: bool) -> None:
        self._s.setValue(f"safety/confirm_{command}", "on" if on else "off")

    @property
    def confirm_notice_shown(self) -> bool:
        return str(self._s.value("safety/notice_shown", "off")) == "on"

    @confirm_notice_shown.setter
    def confirm_notice_shown(self, shown: bool) -> None:
        self._s.setValue("safety/notice_shown", "on" if shown else "off")

    # position of the windows and height of the Command window
    def layout(self):
        return (self._s.value("layout/geometry"), self._s.value("layout/state"),
                self._s.value("layout/splitter"))

    def save_layout(self, geometry, state, splitter) -> None:
        self._s.setValue("layout/geometry", geometry)
        self._s.setValue("layout/state", state)
        self._s.setValue("layout/splitter", splitter)

    def reset(self) -> None:
        self._s.clear()

    def sync(self) -> None:
        self._s.sync()
