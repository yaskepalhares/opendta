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
