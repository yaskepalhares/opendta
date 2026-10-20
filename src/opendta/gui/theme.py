"""Interface colors and fonts.

The default scheme follows Stata 14's light appearance. The syntax colors
of the Do-file Editor were measured from a screenshot of Stata 14's preferences
on macOS (colors of the default scheme). macOS screenshots may convert
the color profile, so the values may differ by a few points from the
originals. Just adjust them here.
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
        "alert": "#cc0000",     # "(N observations deleted)": the data changed (OpenDTA)
        "hint": "#6b5a00",      # error explanations (OpenDTA)
        "hint_label": "#6b5a00",  # "How to fix:", "Learn more:" (bold)
        "hint_link": "#1f5fbf",   # "help topic" links (underlined, clickable)
    })
    bold: frozenset[str] = frozenset({"result", "input", "hint_label", "alert"})


STANDARD = ResultsScheme()

DARK_RESULTS = ResultsScheme(
    background="#1e1f22",
    colors={
        "text": "#d4d4d6",
        "result": "#ffffff",
        "error": "#ff6b68",
        "alert": "#ff6b68",
        "input": "#ffffff",
        "command": "#d4d4d6",
        "hint": "#e5c07b",
        "hint_label": "#e5c07b",
        "hint_link": "#7cb7ff",
    },
)

THEMES = ("light", "dark")


def results_scheme(theme: str, palette: str = "classic") -> ResultsScheme:
    """Colors of the Results window: the chosen color scheme (palettes.py) or, for "classic", the
    Stata-like look of the theme."""
    from .palettes import resolve
    pal = resolve(palette, theme)
    if pal is None:
        return DARK_RESULTS if theme == "dark" else STANDARD
    return ResultsScheme(background=pal.bg, colors={
        "text": pal.fg, "result": pal.fg, "error": pal.error, "alert": pal.error, "input": pal.fg,
        "command": pal.fg, "hint": pal.warn, "hint_label": pal.warn, "hint_link": pal.link})


# Interface colors (windows, panels, menus) in each theme. The theme is chosen in Settings:
# light (default), dark, or system (follows the operating system).
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


def system_theme() -> str:
    """The operating system's appearance (light or dark), ignoring OpenDTA's own override."""
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QGuiApplication
    app = QGuiApplication.instance()
    if app is None:
        return "light"
    hints = app.styleHints()
    if hasattr(hints, "unsetColorScheme"):         # Qt 6.8+: drop the override to read the system
        hints.unsetColorScheme()
    scheme = hints.colorScheme() if hasattr(hints, "colorScheme") else None
    return "dark" if scheme == Qt.ColorScheme.Dark else "light"


def apply_theme(app, theme: str) -> None:
    """Applies the light or dark theme to the whole application."""
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor, QPalette
    from PySide6.QtWidgets import QStyleFactory

    if theme == "system":
        theme = system_theme()
    theme = theme if theme in THEMES else "light"
    hints = app.styleHints()
    if hasattr(hints, "setColorScheme"):           # Qt 6.8+: includes the title bar
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
    """Do-file Editor syntax highlighting (Preferences > Do-file Editor > Syntax)."""

    plain: str = "#010000"
    keywords: str = "#121779"
    comments: str = "#0e7324"
    functions: str = "#252df1"
    macros: str = "#167574"
    strings: str = "#72000b"
    compound_strings: str = "#72000b"
    numbers: str = "#263ff1"
    operators: str = "#010000"
    brace_matching: str = "#fa1b1a"     # in bold
    # aba General
    background: str = "#ffffff"
    selection: str = "#add3fa"
    cursor: str = "#000000"
    invisibles: str = "#3178dd"
    current_line: str = "#fefdcc"
    page_guide: str = "#dcdcdc"
    matches: str = "#fcfb52"
    # line numbers
    line_numbers_background: str = "#dfdfe0"
    line_numbers: str = "#747474"


DOFILE_SYNTAX = SyntaxScheme()


def default_monospace() -> tuple[str, int]:
    """Default monospaced font per system.

    macOS: Menlo (identified in the screenshot of Stata's preferences);
    Windows: Courier New; Linux: DejaVu Sans Mono (same design as Menlo).
    VERIFY: Stata 14's default sizes on each system.
    """
    if sys.platform == "darwin":
        return "Menlo", 12
    if sys.platform.startswith("win"):
        return "Courier New", 10
    return "DejaVu Sans Mono", 10


MONOSPACE_FALLBACKS = ["Menlo", "Courier New", "DejaVu Sans Mono", "Liberation Mono", "Consolas", "monospace"]
