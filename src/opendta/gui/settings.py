"""Preferências persistentes da interface (QSettings).

Ficam no local padrão de cada sistema: registro no Windows, .plist no macOS
e ~/.config no Linux.
"""

from __future__ import annotations

from PySide6.QtCore import QSettings

from .theme import default_monospace

ORG, APP = "OpenDTA", "OpenDTA"

APP_ICON_VARIANTS = {
    "light": "Tabela — claro (padrão)",
    "dark": "Tabela — escuro",
}


class Preferences:
    def __init__(self, settings: QSettings | None = None):
        self._s = settings or QSettings(ORG, APP)

    # ícone do aplicativo
    @property
    def app_icon(self) -> str:
        v = str(self._s.value("appearance/app_icon", "light"))
        return v if v in APP_ICON_VARIANTS else "light"

    @app_icon.setter
    def app_icon(self, value: str) -> None:
        self._s.setValue("appearance/app_icon", value if value in APP_ICON_VARIANTS else "light")

    # tema da interface (independente do sistema)
    @property
    def theme(self) -> str:
        v = str(self._s.value("appearance/theme", "light"))
        return v if v in ("light", "dark") else "light"

    @theme.setter
    def theme(self, value: str) -> None:
        self._s.setValue("appearance/theme", value if value in ("light", "dark") else "light")

    # fonte das janelas Results, Command e Review
    @property
    def hints(self) -> bool:
        """Explicação dos erros (set hints) ligada ao abrir o programa."""
        return str(self._s.value("behavior/hints", "on")) != "off"

    @hints.setter
    def hints(self, on: bool) -> None:
        self._s.setValue("behavior/hints", "on" if on else "off")

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

    # posição das janelas e altura da Command
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
