"""Erros com códigos de retorno no padrão r(#).

Os códigos e as mensagens seguem a documentação pública de códigos de
retorno ([P] error / help r()). As mensagens exatas devem ser conferidas
contra o Stata 14 na suíte compat/.
"""

from __future__ import annotations


class StataError(Exception):
    """Erro que encerra um comando e define _rc."""

    def __init__(self, rc: int, message: str = ""):
        super().__init__(message)
        self.rc = int(rc)
        self.message = message


class ExitRequest(Exception):
    """`exit` dentro de um do-file ou programa (sem erro, a menos que rc != 0)."""

    def __init__(self, rc: int = 0, *, clear_all: bool = False):
        super().__init__(rc)
        self.rc = rc
        self.clear_all = clear_all


class LoopControl(Exception):
    """Base para `continue` e `continue, break`."""


class ContinueLoop(LoopControl):
    pass


class BreakLoop(LoopControl):
    pass


# Mensagens padrão de alguns códigos frequentes.
STANDARD_MESSAGES: dict[int, str] = {
    100: "varlist required",
    101: "varlist not allowed",
    109: "type mismatch",
    110: "already defined",
    111: "invalid syntax",   # texto exibido por `error 111` no Stata (observado)
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
