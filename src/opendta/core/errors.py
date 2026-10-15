"""Errors with return codes in the r(#) style.

Codes and messages follow the public documentation of return
codes ([P] error / help r()). The exact messages should be checked
against Stata 14 in the compat/ suite.
"""

from __future__ import annotations


class StataError(Exception):
    """Error that ends a command and sets _rc."""

    def __init__(self, rc: int, message: str = "", *, hint: str | None = None):
        super().__init__(message)
        self.rc = int(rc)
        self.message = message
        # OpenDTA's own explanation (several lines); without it, core/hints.py tries one
        self.hint = hint
        # (command name, command text), attached by the interpreter
        self.context: tuple[str, str] | None = None


class ExitRequest(Exception):
    """`exit` inside a do-file or program (no error, unless rc != 0)."""

    def __init__(self, rc: int = 0, *, clear_all: bool = False):
        super().__init__(rc)
        self.rc = rc
        self.clear_all = clear_all


class LoopControl(Exception):
    """Base for `continue` and `continue, break`."""


class ContinueLoop(LoopControl):
    pass


class BreakLoop(LoopControl):
    pass


# Standard messages for some frequent codes.
STANDARD_MESSAGES: dict[int, str] = {
    100: "varlist required",
    101: "varlist not allowed",
    109: "type mismatch",
    110: "already defined",
    111: "invalid syntax",   # text shown by `error 111` in Stata (observed)
    130: "expression too long",
    132: "too many '(' or '['",
    133: "unknown function",
    198: "invalid syntax",
    199: "unrecognized command",
    459: "something that should be true of your data is not",
    498: "",
    601: "file not found",
    602: "file already exists",
    603: "file could not be opened",
    9: "assertion is false",
}


def syntax_error(message: str = "invalid syntax") -> StataError:
    return StataError(198, message)


def type_mismatch() -> StataError:
    return StataError(109, "type mismatch")
