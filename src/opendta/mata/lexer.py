"""Tokens do Mata ([M-2] syntax).

Números (1, 1.5, .5, 1e-3), valores missing (., .a–.z), strings entre
aspas ("abc" ou `"abc"'), nomes e operadores. Comentários // e /* */.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..core import missing as M


class MataSyntaxError(Exception):
    """Erro de sintaxe. `incomplete` indica que o texto acabou no meio de
    uma instrução (o bloco mata pede mais uma linha, com o prompt '>')."""

    def __init__(self, message: str, incomplete: bool = False):
        super().__init__(message)
        self.incomplete = incomplete


@dataclass
class Token:
    kind: str          # num, str, id, op, eof
    value: object
    pos: int

    def __repr__(self) -> str:
        return f"{self.kind}:{self.value!r}"


_OPS = sorted([
    "[|", "|]", ":==", ":!=", ":>=", ":<=", ":+", ":-", ":*", ":/", ":^", ":>", ":<", ":&", ":|",
    "==", "!=", ">=", "<=", "&&", "||", "++", "--", "..", "::", "->",
    "+", "-", "*", "/", "^", "'", "#", "!", "?", ":", "=", "<", ">", "&", "|", ",", "\\",
    "(", ")", "[", "]", "{", "}", ";", ".",
], key=len, reverse=True)


def tokenize(text: str) -> list[Token]:
    out: list[Token] = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in " \t\r\n":
            i += 1
            continue
        if text.startswith("//", i):
            j = text.find("\n", i)
            i = n if j == -1 else j
            continue
        if text.startswith("/*", i):
            j = text.find("*/", i + 2)
            if j == -1:
                raise MataSyntaxError("unterminated comment", incomplete=True)
            i = j + 2
            continue
        # strings
        if text.startswith('`"', i):
            j = text.find("\"'", i + 2)
            if j == -1:
                raise MataSyntaxError("unbalanced quotes")
            out.append(Token("str", text[i + 2:j], i))
            i = j + 2
            continue
        if c == '"':
            j = text.find('"', i + 1)
            if j == -1:
                raise MataSyntaxError("unbalanced quotes")
            out.append(Token("str", text[i + 1:j], i))
            i = j + 1
            continue
        # números e missing
        if c.isdigit() or (c == "." and i + 1 < n and text[i + 1].isdigit()):
            j = i
            while j < n and text[j].isdigit():
                j += 1
            if j < n and text[j] == "." and not text.startswith("..", j):
                j += 1
                while j < n and text[j].isdigit():
                    j += 1
            if j < n and text[j] in "eE":
                k = j + 1
                if k < n and text[k] in "+-":
                    k += 1
                if k < n and text[k].isdigit():
                    while k < n and text[k].isdigit():
                        k += 1
                    j = k
            out.append(Token("num", float(text[i:j]), i))
            i = j
            continue
        if c == "." and not text.startswith("..", i):
            prev = out[-1] if out else None
            member = prev is not None and (prev.kind == "id" or (prev.kind == "op" and prev.value in (")", "]")))
            if not member:
                # missing: "." ou ".a"–".z" (sem ser seguido de letra/dígito)
                if i + 1 < n and text[i + 1].islower() and not (
                        i + 2 < n and (text[i + 2].isalnum() or text[i + 2] == "_")):
                    out.append(Token("num", M.missing_code(text[i:i + 2]), i))
                    i += 2
                    continue
                out.append(Token("num", M.SYSMISS, i))
                i += 1
                continue
        if c.isalpha() or c == "_":
            j = i
            while j < n and (text[j].isalnum() or text[j] == "_"):
                j += 1
            out.append(Token("id", text[i:j], i))
            i = j
            continue
        for op in _OPS:
            if text.startswith(op, i):
                out.append(Token("op", op, i))
                i += len(op)
                break
        else:
            raise MataSyntaxError(f"invalid character {c!r}")
    out.append(Token("eof", None, n))
    return out
