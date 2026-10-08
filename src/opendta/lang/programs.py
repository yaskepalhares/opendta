"""Programas definidos pelo usuário ([P] program) e escopos de execução.

    program [define] nome [, rclass eclass sclass nclass byable(recall|onecall)
                            properties(...) sortpreserve plugin]
        ...
    end

Ao chamar `nome argumentos`, o corpo roda num novo escopo de macros locais
com `0' = tudo o que veio depois do nome e `1', `2', ... = as palavras.

Cada programa (e cada do-file) abre um Scope, que guarda o que precisa ser
desfeito na saída: variáveis, nomes e arquivos temporários, e os dados
guardados por preserve. Programas rclass/sclass acumulam os resultados de
`return`/`sreturn` e só os publicam em r()/s() quando terminam sem erro.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import ExitRequest, StataError
from .lexer import LogicalLine
from .words import split_words

if TYPE_CHECKING:
    from ..session import Session

_CLASSES = ("rclass", "eclass", "sclass", "nclass")


@dataclass
class Program:
    name: str
    lines: list[LogicalLine]
    kind: str = "nclass"                 # rclass | eclass | sclass | nclass
    byable: str = ""                     # recall | onecall | ""
    properties: str = ""
    source: str = ""                     # arquivo .ado de origem, se houver
    options: str = ""


@dataclass
class Scope:
    kind: str                            # program | do
    name: str = ""
    program: Program | None = None
    pending: dict = field(default_factory=dict)      # return/sreturn ainda não publicados
    tempvars: list[str] = field(default_factory=list)
    tempnames: list[str] = field(default_factory=list)
    tempfiles: list[Path] = field(default_factory=list)
    preserved: object = None             # ver commands/preserve.py


def parse_definition(rest: str) -> tuple[str, dict[str, str]]:
    """Texto depois de `program [define]` -> (nome, opções)."""
    head, comma, opts = rest.partition(",")
    words = head.split()
    if words and len(words[0]) >= 2 and "define".startswith(words[0]):
        words = words[1:]
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    name = words[0]
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,31}", name):
        raise StataError(198, f"{name} invalid name")
    options: dict[str, str] = {}
    for m in re.finditer(r"([A-Za-z_]+)(?:\(([^)]*)\))?", opts if comma else ""):
        options[m.group(1).lower()] = m.group(2) or ""
    return name, options


def make_program(name: str, options: dict[str, str], lines: list[LogicalLine],
                 source: str = "") -> Program:
    kind = next((c for c in _CLASSES if c in options), "nclass")
    byable = options.get("byable", "")
    if "byable" in options and byable not in ("recall", "onecall"):
        raise StataError(198, "option byable() must be recall or onecall")
    return Program(name, lines, kind, byable, options.get("properties", ""), source)


def call_program(s: "Session", prog: Program, args: str) -> None:
    """Executa o programa com os argumentos dados."""
    argv = split_words(args)
    frame = {str(k): v for k, v in enumerate(argv, start=1)}
    if args.strip():
        frame["0"] = args.strip()
    if len(s.scopes) > 60:
        raise StataError(1000, "system limit exceeded: nested program calls")   # VERIFICAR
    s.macros.push_frame(frame)
    scope = Scope("program", prog.name, prog)
    s.scopes.append(scope)
    ok = False
    try:
        try:
            s.interp.run_program(prog)
        except ExitRequest as e:
            if e.clear_all:
                raise
            if e.rc:
                raise StataError(e.rc, "")
        ok = True
    finally:
        s.scopes.pop()
        s.macros.pop_frame()
        s.close_scope(scope)
    if ok:
        if prog.kind == "rclass":
            s.r = scope.pending
        elif prog.kind == "sclass":
            s.sret = scope.pending


def current_program_scope(s: "Session") -> Scope | None:
    for sc in reversed(s.scopes):
        if sc.kind == "program":
            return sc
    return None
