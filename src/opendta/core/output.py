"""Saída de texto (janela Results, console, logs).

Toda saída passa por `Output`, que mantém a coluna atual (necessária para
_column() e _continue) e repassa os trechos para os ouvintes registrados:
a janela Results da interface, o console e arquivos de log.

Estilos seguem as categorias do Stata: "text", "result", "error",
"input" e "command" (eco do comando digitado).
"""

from __future__ import annotations

from typing import Callable

Listener = Callable[[str, str], None]

STYLES = ("text", "result", "error", "input", "command")


class Output:
    def __init__(self) -> None:
        self._listeners: list[Listener] = []
        self.column = 0           # coluna atual (0 = início da linha)
        self.quiet_depth = 0      # > 0 dentro de quietly
        self.noisy_depth = 0      # noisily dentro de quietly
        self.capture_depth = 0    # > 0 dentro de capture: silencia até erros

    # -- ouvintes -------------------------------------------------------
    def add_listener(self, fn: Listener) -> None:
        self._listeners.append(fn)

    def remove_listener(self, fn: Listener) -> None:
        if fn in self._listeners:
            self._listeners.remove(fn)

    # -- estado de silêncio -----------------------------------------------
    @property
    def suppressed(self) -> bool:
        return self.quiet_depth > 0 and self.noisy_depth == 0

    # -- escrita -----------------------------------------------------------
    def write(self, text: str, style: str = "text", *, force: bool = False) -> None:
        if not text:
            return
        if force:
            # mensagens de erro aparecem sob quietly, mas não sob capture
            if self.capture_depth > 0 and self.noisy_depth == 0:
                return
        elif self.suppressed:
            return
        for fn in list(self._listeners):
            fn(text, style)
        last_nl = text.rfind("\n")
        if last_nl == -1:
            self.column += len(text)
        else:
            self.column = len(text) - last_nl - 1

    def newline(self, n: int = 1, *, force: bool = False) -> None:
        self.write("\n" * n, "text", force=force)

    def ensure_line_start(self, *, force: bool = False) -> None:
        if self.column != 0:
            self.newline(force=force)

    def error(self, message: str, rc: int, *, show_rc: bool = True,
              hints: list[str] | None = None) -> None:
        """Mensagem de erro e r(#); — mostradas mesmo sob quietly (mas não
        sob capture). `hints`: explicação do OpenDTA, entre as duas."""
        self.ensure_line_start(force=True)
        if message:
            self.write(message + "\n", "error", force=True)
        for k, line in enumerate(hints or []):
            prefix = "  \u2192 " if k == 0 else "    "
            self.write(prefix + line + "\n", "hint", force=True)
        if show_rc:
            self.write(f"r({rc});\n", "error", force=True)

    def echo_command(self, lines) -> None:
        """Eco '. comando' como a janela Results faz; linhas de continuação
        aparecem com '> '."""
        if isinstance(lines, str):
            lines = (lines,)
        self.ensure_line_start()
        first, *rest = lines
        self.write(". " + first + "\n", "command")
        for cont in rest:
            self.write("> " + cont + "\n", "command")

    def end_command(self) -> None:
        """Linha em branco depois de cada comando (como no Stata)."""
        self.ensure_line_start()
        self.write("\n")


class Capture:
    """Ouvinte que acumula o texto (usado em testes e no modo batch)."""

    def __init__(self) -> None:
        self.parts: list[tuple[str, str]] = []

    def __call__(self, text: str, style: str) -> None:
        self.parts.append((text, style))

    @property
    def text(self) -> str:
        return "".join(t for t, _ in self.parts)

    def clear(self) -> None:
        self.parts.clear()
