"""Cores e fontes da interface.

O esquema padrão imita a aparência clara da janela Results do Stata 14:
fundo branco, texto em preto, resultados destacados, erros em vermelho.
Os valores exatos são ajustáveis aqui e serão revisados por comparação
visual (ver docs/interface.md).
"""

from __future__ import annotations

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

MONOSPACE_FAMILIES = ["Courier New", "Consolas", "Menlo", "DejaVu Sans Mono", "Liberation Mono", "monospace"]
MONOSPACE_SIZE = 10
