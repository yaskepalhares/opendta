"""display — manual [P] display.

Supported directives: "string", `"string"', expressions, %fmt, as text|txt,
as result|res, as error|err, as input|inp, in smcl, _newline[(#)] (_n),
_column(#) (_col), _skip(#), _dup(#), _char(#), _continue (_c), _request(macname).
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
    requested = False
    i = 0
    text = args
    n = len(text)

    # output is accumulated and written only at the end: if there is an error, nothing appears
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
        if m and m.group(1) == "request":
            # _request(macname): shows what came before as the prompt, reads a line from the
            # console and stores it in the global macname (local when it starts with "_")
            target = (m.group(3) or "").strip()
            if not re.fullmatch(r"_?[A-Za-z_]\w*", target):
                raise StataError(198, "invalid syntax")
            i += m.end()
            for piece, st in pending:
                out.write(piece, st)
            pending.clear()
            answer = s.request_line()
            out.write(answer + "\n", "input")      # VERIFY: echo of the typed text
            col = 0
            requested = True
            if target.startswith("_"):
                s.macros.set_local(target[1:], answer)
            else:
                s.macros.set_global(target, answer)
            continue
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
            continue

        m = _FORMAT.match(rest)
        if m and rest.startswith("%"):
            fmt = parse_format(m.group(0))
            i += m.end()
            continue

        # "text" or `"text"' at the start of a directive is literal text, not the start of an
        # expression: display "a" + 1 shows a1 (observed in Stata); to concatenate strings use
        # parentheses: display ("a" + "b")
        lit = _literal(rest)
        if lit is not None:
            text_value, consumed = lit
            last_literal = text_value
            i += consumed
            emit(format_value(text_value, fmt) if fmt is not None and fmt.kind == "s" else text_value)
            fmt = None
            continue

        # expression
        p = Parser(rest, stop_on_unknown=True)
        if p.at_end():
            raise StataError(198, "invalid syntax")
        node = p.parse_expr()
        try:
            value = evaluate(node, s.context)
        except StataError as e:
            # display "a" + "b": o Stata responde 'a+"b" invalid name'
            # (VERIFY a regra geral desta mensagem)
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

    if not cont and not (requested and not pending):
        put("\n", style)
    for piece, st in pending:
        out.write(piece, st)


def _literal(rest: str):
    """If `rest` starts with a quoted string, returns (text, length)."""
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
