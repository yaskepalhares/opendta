"""display — manual [P] display.

Diretivas suportadas: "string", `"string"', expressões, %fmt, as text|txt,
as result|res, as error|err, as input|inp, in smcl, _newline[(#)] (_n),
_column(#) (_col), _skip(#), _dup(#), _char(#), _continue (_c).
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.formats import DISPLAY_NUMERIC, format_value, parse_format
from ..lang.expr import Parser, evaluate
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_STYLES = {
    "text": "text", "txt": "text",
    "result": "result", "res": "result",
    "error": "error", "err": "error",
    "input": "input", "inp": "input",
}

_DIRECTIVE = re.compile(
    r"_(newline|n|column|col|skip|dup|char|continue|c|request)\b(\s*\(\s*([^)]*)\s*\))?")
_FORMAT = re.compile(r"%-?~?0?\d*(\.\d+)?[a-zA-Z]+c?")


@command("display", "di")
def cmd_display(s: "Session", args: str) -> None:
    out = s.output
    style = "text"
    fmt = None
    dup = 1
    cont = False
    i = 0
    text = args
    n = len(text)

    # a saída é acumulada e só escrita no fim: se houver erro, nada aparece
    pending: list[tuple[str, str]] = []
    col = out.column
    last_literal: str | None = None

    def put(piece: str, st: str) -> None:
        nonlocal col
        if not piece:
            return
        pending.append((piece, st))
        nl = piece.rfind("\n")
        col = col + len(piece) if nl == -1 else len(piece) - nl - 1

    def emit(piece: str) -> None:
        nonlocal dup
        put(piece * dup, style)
        dup = 1

    while i < n:
        while i < n and text[i] in " \t":
            i += 1
        if i >= n:
            break
        rest = text[i:]

        m = re.match(r"as\s+(\w+)", rest)
        if m and m.group(1) in _STYLES:
            style = _STYLES[m.group(1)]
            i += m.end()
            continue
        m = re.match(r"in\s+smcl\b", rest)
        if m:
            i += m.end()
            continue

        m = _DIRECTIVE.match(rest)
        if m:
            name, arg = m.group(1), m.group(3)
            num = None
            if arg is not None and arg.strip():
                num = s.eval(arg)
                if isinstance(num, str):
                    raise StataError(109, "type mismatch")
                num = int(num)
            i += m.end()
            if name in ("newline", "n"):
                put("\n" * (num if num is not None else 1), style)
            elif name in ("column", "col"):
                target = (num or 1) - 1
                if col < target:
                    put(" " * (target - col), style)
            elif name == "skip":
                put(" " * (num if num is not None else 1), style)
            elif name == "dup":
                dup = max(num or 0, 0)
            elif name == "char":
                emit(chr(num or 0))
            elif name in ("continue", "c"):
                cont = True
            elif name == "request":
                raise StataError(198, "_request() ainda não é suportado")
            continue

        m = _FORMAT.match(rest)
        if m and rest.startswith("%"):
            fmt = parse_format(m.group(0))
            i += m.end()
            continue

        # "texto" ou `"texto"' no início de uma diretiva é texto literal, não
        # o começo de uma expressão: display "a" + 1 mostra a1 (observado no
        # Stata); para concatenar strings use parênteses: display ("a" + "b")
        lit = _literal(rest)
        if lit is not None:
            text_value, consumed = lit
            last_literal = text_value
            i += consumed
            emit(format_value(text_value, fmt) if fmt is not None and fmt.kind == "s" else text_value)
            fmt = None
            continue

        # expressão
        p = Parser(rest, stop_on_unknown=True)
        if p.at_end():
            raise StataError(198, "invalid syntax")
        node = p.parse_expr()
        try:
            value = evaluate(node, s.context)
        except StataError as e:
            # display "a" + "b": o Stata responde 'a+"b" invalid name'
            # (VERIFICAR a regra geral desta mensagem)
            if e.rc == 198 and last_literal is not None and rest.lstrip().startswith("+"):
                raise StataError(198, f'{last_literal}{"".join(rest.split())} invalid name')
            raise
        consumed = p.offset
        if consumed == 0:
            raise StataError(198, "invalid syntax")
        i += consumed
        if isinstance(value, str):
            piece = format_value(value, fmt) if fmt is not None and fmt.kind == "s" else value
        elif fmt is not None:
            piece = format_value(value, fmt)
        else:
            piece = format_value(value, DISPLAY_NUMERIC, pad=False)
        fmt = None
        emit(piece)

    if not cont:
        put("\n", style)
    for piece, st in pending:
        out.write(piece, st)


def _literal(rest: str):
    """Se `rest` começa com uma string entre aspas, devolve (texto, tamanho)."""
    if rest.startswith('`"'):
        depth, j = 1, 2
        while j < len(rest) and depth:
            if rest.startswith('`"', j):
                depth += 1
                j += 2
            elif rest.startswith("\"'", j):
                depth -= 1
                j += 2
            else:
                j += 1
        if depth:
            raise StataError(198, "unmatched quote")
        return rest[2:j - 2], j
    if rest.startswith('"'):
        j = rest.find('"', 1)
        if j == -1:
            raise StataError(198, "unmatched quote")
        return rest[1:j], j + 1
    return None
