"""Ícones próprios do OpenDTA.

Os SVGs em gui/icons/ usam `currentColor`; aqui a cor é trocada pela cor do
texto do tema atual (claro ou escuro), e o estado desabilitado recebe uma
versão mais apagada. Assim um único desenho serve aos dois temas.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from PySide6.QtCore import QByteArray, QRectF, QSize, Qt
from PySide6.QtGui import QGuiApplication, QIcon, QPainter, QPalette, QPixmap
from PySide6.QtSvg import QSvgRenderer

ICON_DIR = Path(__file__).with_name("icons")
APPICON_DIR = Path(__file__).with_name("appicon")

# cores no espírito dos ícones do macOS/iOS
LIGHT_INK = "#3a3a3c"
DARK_INK = "#e5e5ea"
ACCENT_RED = "#ff453a"


def is_dark() -> bool:
    app = QGuiApplication.instance()
    if app is None:
        return False
    hints = app.styleHints()
    scheme = getattr(hints, "colorScheme", None)
    if scheme is not None:
        try:
            from PySide6.QtCore import Qt as _Qt
            value = scheme()
            if value == _Qt.ColorScheme.Dark:
                return True
            if value == _Qt.ColorScheme.Light:
                return False
        except (AttributeError, TypeError):
            pass
    return app.palette().color(QPalette.ColorRole.Window).lightness() < 128


def _render(name: str, color: str, size: int, ratio: float) -> QPixmap:
    svg = (ICON_DIR / f"{name}.svg").read_text(encoding="utf-8").replace("currentColor", color)
    renderer = QSvgRenderer(QByteArray(svg.encode("utf-8")))
    px = int(round(size * ratio))
    pm = QPixmap(px, px)
    pm.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pm)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    renderer.render(painter, QRectF(0, 0, px, px))
    painter.end()
    pm.setDevicePixelRatio(ratio)
    return pm


def _faded(pm: QPixmap, opacity: float) -> QPixmap:
    out = QPixmap(pm.size())
    out.setDevicePixelRatio(pm.devicePixelRatio())
    out.fill(Qt.GlobalColor.transparent)
    painter = QPainter(out)
    painter.setOpacity(opacity)
    painter.drawPixmap(0, 0, pm)
    painter.end()
    return out


@lru_cache(maxsize=None)
def themed_icon(name: str, dark: bool, accent: str | None = None) -> QIcon:
    base = DARK_INK if dark else LIGHT_INK
    ink = accent or base
    app = QGuiApplication.instance()
    ratio = app.devicePixelRatio() if app is not None else 1.0
    icon = QIcon()
    for size in (16, 20, 22, 24, 32):
        icon.addPixmap(_render(name, ink, size, ratio), QIcon.Mode.Normal)
        icon.addPixmap(_faded(_render(name, base, size, ratio), 0.38), QIcon.Mode.Disabled)
    return icon


def icon(name: str, accent: str | None = None) -> QIcon:
    return themed_icon(name, is_dark(), accent)


def app_icon() -> QIcon:
    """Ícone do aplicativo (variante clara ou escura, conforme o tema)."""
    variant = "opendta-dark.png" if is_dark() else "opendta.png"
    path = APPICON_DIR / variant
    if not path.exists():
        path = APPICON_DIR / "opendta.png"
    return QIcon(str(path)) if path.exists() else QIcon()


TOOLBAR_SIZE = QSize(20, 20)
