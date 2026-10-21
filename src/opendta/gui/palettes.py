"""Color schemes for the Results window and the Do-file Editor.

"Classic Stata" is the original look and follows the light or dark theme of the interface (its colors
live in theme.py and dofile_editor.py). "High contrast" and "Soft" also follow the theme. The others
are the schemes most used in code editors; they keep their own colors in any theme. The names are
those of the public color schemes; the values are written out here.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..i18n import L


@dataclass(frozen=True)
class Palette:
    bg: str
    fg: str
    error: str
    warn: str            # error explanations
    link: str            # "help topic" links
    command: str
    string: str
    comment: str
    macro: str
    number: str
    gutter: str
    gutter_text: str
    current: str         # line of the cursor


_HC_LIGHT = Palette("#ffffff", "#000000", "#c00000", "#7a4a00", "#0033cc", "#0000cc", "#8b0000", "#006400",
                    "#7a0099", "#a04000", "#eeeeee", "#444444", "#eef3ff")
_HC_DARK = Palette("#000000", "#ffffff", "#ff6b6b", "#ffd866", "#66b3ff", "#7fc8ff", "#ffb38a", "#7ee787",
                   "#e0a8ff", "#ffd479", "#111111", "#bbbbbb", "#1c1c1c")
_SOFT_LIGHT = Palette("#f7f5f0", "#4a4a4a", "#b5524a", "#8a6d3b", "#4f7fb0", "#4f6f9f", "#a67c52", "#7d9a7d",
                      "#8f6fa8", "#b08a4a", "#efece5", "#9a968c", "#efece4")
_SOFT_DARK = Palette("#2a2b2e", "#c9c9c6", "#d98880", "#d4b97a", "#8fb4d8", "#8fb0d8", "#d4a781", "#8fb08f",
                     "#b79fcc", "#d4b97a", "#252629", "#7d7e82", "#323337")

# key -> (English label, Portuguese label, light variant, dark variant); fixed schemes have one palette
_FIXED = {
    "solarized_light": ("Solarized Light", Palette("#fdf6e3", "#657b83", "#dc322f", "#b58900", "#268bd2", "#268bd2",
                                                    "#2aa198", "#93a1a1", "#d33682", "#cb4b16", "#eee8d5",
                                                    "#93a1a1", "#eee8d5")),
    "solarized_dark": ("Solarized Dark", Palette("#002b36", "#839496", "#dc322f", "#b58900", "#268bd2", "#268bd2",
                                                  "#2aa198", "#586e75", "#d33682", "#cb4b16", "#073642",
                                                  "#586e75", "#073642")),
    "github_light": ("GitHub Light", Palette("#ffffff", "#24292f", "#cf222e", "#9a6700", "#0969da", "#cf222e",
                                              "#0a3069", "#6e7781", "#8250df", "#0550ae", "#f6f8fa",
                                              "#8c959f", "#f6f8fa")),
    "monokai": ("Monokai", Palette("#272822", "#f8f8f2", "#f92672", "#e6db74", "#66d9ef", "#f92672", "#e6db74",
                                   "#75715e", "#ae81ff", "#ae81ff", "#2d2e27", "#90908a", "#3e3d32")),
    "dracula": ("Dracula", Palette("#282a36", "#f8f8f2", "#ff5555", "#f1fa8c", "#8be9fd", "#ff79c6", "#f1fa8c",
                                   "#6272a4", "#bd93f9", "#ffb86c", "#21222c", "#6272a4", "#44475a")),
    "one_dark": ("One Dark", Palette("#282c34", "#abb2bf", "#e06c75", "#e5c07b", "#61afef", "#c678dd", "#98c379",
                                     "#5c6370", "#56b6c2", "#d19a66", "#21252b", "#636d83", "#2c313c")),
    "nord": ("Nord", Palette("#2e3440", "#d8dee9", "#bf616a", "#ebcb8b", "#88c0d0", "#81a1c1", "#a3be8c",
                             "#616e88", "#b48ead", "#d08770", "#2b303b", "#4c566a", "#3b4252")),
    "gruvbox_dark": ("Gruvbox Dark", Palette("#282828", "#ebdbb2", "#fb4934", "#fabd2f", "#83a598", "#fb4934",
                                             "#b8bb26", "#928374", "#d3869b", "#fe8019", "#32302f", "#7c6f64",
                                             "#3c3836")),
}

_FOLLOWING = {
    "high_contrast": (_HC_LIGHT, _HC_DARK),
    "soft": (_SOFT_LIGHT, _SOFT_DARK),
}

DEFAULT = "classic"
KEYS = [DEFAULT, *_FOLLOWING, *_FIXED]


def labels() -> list[tuple[str, str]]:
    """(key, name) in the language of the interface, for the selection menus."""
    names = {"classic": L("Classic Stata (follows the theme)", "Clássico Stata (segue o tema)"),
             "high_contrast": L("High contrast (follows the theme)", "Alto contraste (segue o tema)"),
             "soft": L("Soft (follows the theme)", "Suave (segue o tema)")}
    return [(k, names.get(k) or _FIXED[k][0]) for k in KEYS]


def resolve(key: str, theme: str) -> Palette | None:
    """The palette for a theme ("light" or "dark"); None for "classic" (the original colors)."""
    if key in _FOLLOWING:
        light, dark = _FOLLOWING[key]
        return dark if theme == "dark" else light
    if key in _FIXED:
        return _FIXED[key][1]
    return None
