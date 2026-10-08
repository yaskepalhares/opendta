"""Macros locais e globais e sua expansão.

Expansão (manual [U] 18.3 e [P] macro):

* `nome'  -> macro local (escopo do do-file/programa em execução)
* $nome e ${nome} -> macro global
* referências aninhadas são resolvidas de dentro para fora: ``i''
* `=exp' avalia a expressão; `:função' chama uma função estendida
* `++i', `i++', `--i', `i--' incrementam/decrementam a local
* macro inexistente expande para vazio
* o texto resultante da expansão não é expandido de novo
* `"' abre aspas compostas e não é tratado como referência de macro
"""

from __future__ import annotations

import re
from typing import Callable

from ..core.errors import StataError

_GLOBAL_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]{0,31}")


class MacroStore:
    def __init__(self) -> None:
        self.globals: dict[str, str] = {}
        self._frames: list[dict[str, str]] = [{}]

    # -- escopo de locais ------------------------------------------------
    @property
    def locals(self) -> dict[str, str]:
        return self._frames[-1]

    def push_frame(self, initial: dict[str, str] | None = None) -> None:
        self._frames.append(dict(initial or {}))

    def pop_frame(self) -> None:
        if len(self._frames) > 1:
            self._frames.pop()

    @property
    def depth(self) -> int:
        return len(self._frames)

    # -- acesso ---------------------------------------------------------------
    def get_local(self, name: str) -> str:
        return self.locals.get(name, "")

    def set_local(self, name: str, value: str) -> None:
        if value == "":
            self.locals.pop(name, None)
        else:
            self.locals[name] = value

    def get_global(self, name: str) -> str:
        return self.globals.get(name, "")

    def set_global(self, name: str, value: str) -> None:
        if value == "":
            self.globals.pop(name, None)
        else:
            self.globals[name] = value

    def clear_globals(self) -> None:
        self.globals.clear()


Resolver = Callable[[str], str]


def expand(text: str, store: MacroStore, *,
           eval_inline: Resolver | None = None,
           extended: Resolver | None = None) -> str:
    """Expande todas as referências a macros em `text`.

    eval_inline recebe o texto após '=' em `=exp' e devolve a string resultante.
    extended recebe o texto após ':' em `:função' e devolve a string resultante.
    """
    s = text
    stack: list[int] = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == "`":
            if i + 1 < len(s) and s[i + 1] == '"':
                i += 2  # aspas compostas
                continue
            stack.append(i)
            i += 1
            continue
        if c == "'" and stack:
            p = stack.pop()
            inner = s[p + 1:i]
            value = _resolve_local(inner, store, eval_inline, extended)
            s = s[:p] + value + s[i + 1:]
            i = p + len(value)
            continue
        if c == "$":
            if i + 1 < len(s) and s[i + 1] == "{":
                j = s.find("}", i + 2)
                if j != -1:
                    name = s[i + 2:j]
                    value = store.get_global(name)
                    s = s[:i] + value + s[j + 1:]
                    i += len(value)
                    continue
            m = _GLOBAL_NAME.match(s, i + 1)
            if m:
                value = store.get_global(m.group(0))
                s = s[:i] + value + s[m.end():]
                i += len(value)
                continue
        i += 1
    return s


def _resolve_local(inner: str, store: MacroStore,
                   eval_inline: Resolver | None,
                   extended: Resolver | None) -> str:
    stripped = inner.strip()
    if stripped.startswith("="):
        if eval_inline is None:
            raise StataError(198, "invalid syntax")
        return eval_inline(stripped[1:])
    if stripped.startswith(":"):
        if extended is None:
            raise StataError(198, "invalid syntax")
        return extended(stripped[1:].strip())
    for op, delta in (("++", 1), ("--", -1)):
        if stripped.startswith(op) or stripped.endswith(op):
            prefix = stripped.startswith(op)
            name = stripped[2:] if prefix else stripped[:-2]
            current = store.get_local(name)
            try:
                num = float(current) if current else 0.0
            except ValueError:
                raise StataError(198, f"{name} is not a number")
            new = num + delta
            new_s = _num_str(new)
            store.set_local(name, new_s)
            return new_s if prefix else _num_str(num)
    if stripped.startswith("macval(") and stripped.endswith(")"):
        return store.get_local(stripped[7:-1].strip())
    return store.get_local(inner)


def _num_str(x: float) -> str:
    return str(int(x)) if x == int(x) else repr(x)
