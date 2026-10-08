"""Cores e fontes da interface.

O esquema padrão segue a aparência clara do Stata 14. As cores de sintaxe
do Do-file Editor foram medidas numa captura das preferências do Stata 14
no macOS (cores do esquema padrão). Capturas de tela do macOS podem converter
o perfil de cor, então os valores podem diferir em poucos pontos dos
originais. Basta ajustar aqui.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field


@dataclass(frozen=True)
class ResultsScheme:
    background: str = "#ffffff"
    colors: dict[str, str] = field(default_factory=lambda: {
        "text": "#000000",
        "result": "#000000",
        "error": "#cc0000",
        "input": "#000000",
        "command": "#000000",
    })
    bold: frozenset[str] = frozenset({"result", "input"})


STANDARD = ResultsScheme()

DARK_RESULTS = ResultsScheme(
    background="#1e1f22",
    colors={
        "text": "#d4d4d6",
        "result": "#ffffff",
        "error": "#ff6b68",
        "input": "#ffffff",
        "command": "#d4d4d6",
    },
)

THEMES = ("light", "dark")


def results_scheme(theme: str) -> ResultsScheme:
    return DARK_RESULTS if theme == "dark" else STANDARD


# Cores da interface (janelas, painéis, menus) em cada tema. A interface
# não segue o modo claro/escuro do sistema: o tema é escolhido em
# Preferences, com o claro como padrão.
_UI = {
    "light": {
        "Window": "#ececec", "WindowText": "#1d1d1f", "Base": "#ffffff",
        "AlternateBase": "#f5f5f7", "Text": "#1d1d1f", "Button": "#f2f2f2",
        "ButtonText": "#1d1d1f", "Highlight": "#0a84ff", "HighlightedText": "#ffffff",
        "ToolTipBase": "#ffffff", "ToolTipText": "#1d1d1f", "PlaceholderText": "#8e8e93",
        "Mid": "#c8c8cc", "Dark": "#a0a0a5", "Light": "#ffffff", "Midlight": "#e5e5ea",
        "Shadow": "#8e8e93", "BrightText": "#ff3b30", "Link": "#0a84ff",
    },
    "dark": {
        "Window": "#2b2c30", "WindowText": "#e5e5ea", "Base": "#1e1f22",
        "AlternateBase": "#26272b", "Text": "#e5e5ea", "Button": "#36373b",
        "ButtonText": "#e5e5ea", "Highlight": "#0a84ff", "HighlightedText": "#ffffff",
        "ToolTipBase": "#36373b", "ToolTipText": "#e5e5ea", "PlaceholderText": "#8e8e93",
        "Mid": "#48494e", "Dark": "#1a1a1d", "Light": "#4a4b50", "Midlight": "#3a3b40",
        "Shadow": "#000000", "BrightText": "#ff453a", "Link": "#64a8ff",
    },
}
_DISABLED_TEXT = {"light": "#a1a1a6", "dark": "#6e6e73"}


def apply_theme(app, theme: str) -> None:
    """Aplica o tema claro ou escuro à aplicação inteira."""
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QPalette
    from PySide6.QtWidgets import QStyleFactory

    theme = theme if theme in THEMES else "light"
    hints = app.styleHints()
    if hasattr(hints, "setColorScheme"):           # Qt 6.8+: inclui a barra de título
        hints.setColorScheme(Qt.ColorScheme.Dark if theme == "dark" else Qt.ColorScheme.Light)
    fusion = QStyleFactory.create("Fusion")
    if fusion is not None:
        app.setStyle(fusion)
    pal = QPalette()
    for role, color in _UI[theme].items():
        pal.setColor(getattr(QPalette.ColorRole, role), QColor(color))
    for role in ("Text", "WindowText", "ButtonText"):
        pal.setColor(QPalette.ColorGroup.Disabled, getattr(QPalette.ColorRole, role),
                     QColor(_DISABLED_TEXT[theme]))
    app.setPalette(pal)


@dataclass(frozen=True)
class SyntaxScheme:
    """Realce de sintaxe do Do-file Editor (Preferences > Do-file Editor > Syntax)."""

    plain: str = "#010000"
    keywords: str = "#121779"
    comments: str = "#0e7324"
    functions: str = "#252df1"
    macros: str = "#167574"
    strings: str = "#72000b"
    compound_strings: str = "#72000b"
    numbers: str = "#263ff1"
    operators: str = "#010000"
    brace_matching: str = "#fa1b1a"     # em negrito
    # aba General
    background: str = "#ffffff"
    selection: str = "#add3fa"
    cursor: str = "#000000"
    invisibles: str = "#3178dd"
    current_line: str = "#fefdcc"
    page_guide: str = "#dcdcdc"
    matches: str = "#fcfb52"
    # números de linha
    line_numbers_background: str = "#dfdfe0"
    line_numbers: str = "#747474"


DOFILE_SYNTAX = SyntaxScheme()


def default_monospace() -> tuple[str, int]:
    """Fonte monoespaçada padrão por sistema.

    macOS: Menlo (identificada na captura das preferências do Stata);
    Windows: Courier New; Linux: DejaVu Sans Mono (mesmo desenho da Menlo).
    VERIFICAR: tamanhos padrão do Stata 14 em cada sistema.
    """
    if sys.platform == "darwin":
        return "Menlo", 12
    if sys.platform.startswith("win"):
        return "Courier New", 10
    return "DejaVu Sans Mono", 10


MONOSPACE_FALLBACKS = ["Menlo", "Courier New", "DejaVu Sans Mono", "Liberation Mono", "Consolas", "monospace"]
