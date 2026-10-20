"""Text output (Results window, console, logs).

All output goes through `Output`, which keeps the current column (needed for
_column() and _continue) and passes the pieces on to the registered listeners:
the Results window of the interface, the console and log files.

Styles follow Stata's categories: "text", "result", "error",
"input" and "command" (echo of the typed command).
"""

from __future__ import annotations

import contextlib

from typing import Callable

Listener = Callable[[str, str], None]

STYLES = ("text", "result", "error", "input", "command")


class Output:
    def __init__(self) -> None:
        self._listeners: list[Listener] = []
        self.column = 0           # current column (0 = start of the line)
        self.quiet_depth = 0      # > 0 inside quietly
        self.noisy_depth = 0      # noisily inside quietly
        self.capture_depth = 0    # > 0 inside capture: silences even errors

    # -- ouvintes -------------------------------------------------------
    def add_listener(self, fn: Listener) -> None:
        self._listeners.append(fn)

    def remove_listener(self, fn: Listener) -> None:
        if fn in self._listeners:
            self._listeners.remove(fn)

    @contextlib.contextmanager
    def quietly(self):
        """Silences everything inside, even when an enclosing noisily is active (as a nested
        quietly does in Stata)."""
        self.quiet_depth += 1
        saved = self.noisy_depth
        self.noisy_depth = 0
        try:
            yield
        finally:
            self.quiet_depth -= 1
            self.noisy_depth = saved

    # -- silence state -----------------------------------------------
    @property
    def suppressed(self) -> bool:
        return self.quiet_depth > 0 and self.noisy_depth == 0

    # -- escrita -----------------------------------------------------------
    def write(self, text: str, style: str = "text", *, force: bool = False) -> None:
        if not text:
            return
        if force:
            # error messages appear under quietly, but not under capture
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
              hint: list[tuple[str, str]] = ()) -> None:
        """Error message and r(#); — shown even under quietly (but not
        under capture). `hint`: OpenDTA's explanation as (text, style) pieces, between the two."""
        self.ensure_line_start(force=True)
        if message:
            self.write(message + "\n", "error", force=True)
        for text, style in hint:
            self.write(text, style, force=True)
        if show_rc:
            self.write(f"r({rc});\n", "error", force=True)

    def echo_command(self, lines) -> None:
        """Echo '. command' as the Results window does; continuation lines
        appear with '> '."""
        if isinstance(lines, str):
            lines = (lines,)
        self.ensure_line_start()
        first, *rest = lines
        self.write(". " + first + "\n", "command")
        for cont in rest:
            self.write("> " + cont + "\n", "command")

    def end_command(self) -> None:
        """Blank line after each command (as in Stata)."""
        self.ensure_line_start()
        self.write("\n")


class Capture:
    """Listener that accumulates the text (used in tests and in batch mode)."""

    def __init__(self) -> None:
        self.parts: list[tuple[str, str]] = []

    def __call__(self, text: str, style: str) -> None:
        self.parts.append((text, style))

    @property
    def text(self) -> str:
        return "".join(t for t, _ in self.parts)

    def clear(self) -> None:
        self.parts.clear()
