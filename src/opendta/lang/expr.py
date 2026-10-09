"""Expressões: tokenização, análise sintática (Pratt) e avaliação.

Precedência, da mais forte para a mais fraca (manual [U] 13.2.5):

    !  (ou ~)          negação lógica
    ^                  potência (associativa à esquerda)
    -                  menos unário
    /  *
    -  +
    != ~= > < >= <= == relacionais
    &
    |

Regras de missing: qualquer operação aritmética com missing resulta em `.`;
nas comparações, missing é maior que qualquer número (. < .a < ... < .z);
em contexto lógico, todo valor diferente de 0 é verdadeiro, inclusive missing.

Nesta fase a avaliação é escalar. Na fase 1 o avaliador passa a operar
sobre vetores (uma posição por observação) para generate/replace/if.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any, Callable, Protocol

from ..core import missing as M
from ..core.errors import StataError, syntax_error, type_mismatch

Value = Any  # float ou str

# ---------------------------------------------------------------------------
# Tokens
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Tok:
    kind: str   # num, str, name, op, end
    value: Any
    pos: int


_NUM = re.compile(r"(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")
_MISS = re.compile(r"\.([a-z])?(?![A-Za-z0-9_])")
_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
# operadores de séries temporais numa expressão: L.y, L2.y, D.y, LD.y ([U] 13.9)
_TSNAME = re.compile(r"(?:(?:[LFDS]\d*)+\.)+[A-Za-z_][A-Za-z0-9_]*")
_OPS = ("==", "!=", "~=", ">=", "<=", "+", "-", "*", "/", "^", "!", "~",
        ">", "<", "&", "|", "(", ")", "[", "]", ",", "=")


def tokenize(text: str, *, stop_on_unknown: bool = False) -> list[Tok]:
    """Com stop_on_unknown=True, um caractere desconhecido encerra a lista de
    tokens (usado por display e pelo `if exp comando`)."""
    toks: list[Tok] = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in " \t\r\n":
            i += 1
            continue
        if c == '"':
            j = text.find('"', i + 1)
            if j == -1:
                raise StataError(198, "unmatched quote")
            toks.append(Tok("str", text[i + 1:j], i))
            i = j + 1
            continue
        if text.startswith('`"', i):
            depth, j = 1, i + 2
            while j < n and depth:
                if text.startswith('`"', j):
                    depth += 1
                    j += 2
                elif text.startswith("\"'", j):
                    depth -= 1
                    j += 2
                else:
                    j += 1
            if depth:
                raise StataError(198, "unmatched quote")
            toks.append(Tok("str", text[i + 2:j - 2], i))
            i = j
            continue
        m = _NUM.match(text, i)
        if m and (c.isdigit() or (c == "." and i + 1 < n and text[i + 1].isdigit())):
            toks.append(Tok("num", float(m.group(0)), i))
            i = m.end()
            continue
        if c == ".":
            m = _MISS.match(text, i)
            if m:
                toks.append(Tok("num", M.missing_code(m.group(0)), i))
                i = m.end()
                continue
        m = _TSNAME.match(text, i) or _NAME.match(text, i)
        if m:
            toks.append(Tok("name", m.group(0), i))
            i = m.end()
            continue
        for op in _OPS:
            if text.startswith(op, i):
                toks.append(Tok("op", op, i))
                i += len(op)
                break
        else:
            if stop_on_unknown:
                toks.append(Tok("end", None, i))
                return toks
            raise syntax_error()
    toks.append(Tok("end", None, n))
    return toks


# ---------------------------------------------------------------------------
# Árvore sintática
# ---------------------------------------------------------------------------


@dataclass
class Node:
    pass


@dataclass
class Num(Node):
    value: float


@dataclass
class Str(Node):
    value: str


@dataclass
class Name(Node):
    name: str


@dataclass
class Call(Node):
    name: str
    args: list[Node] = field(default_factory=list)
    raw: str = ""  # texto bruto do argumento (r(), e(), c())


@dataclass
class Subscript(Node):
    name: str
    index: Node | str
    index2: Node | None = None          # A[i,j]: elemento de matriz


@dataclass
class Unary(Node):
    op: str
    operand: Node


@dataclass
class Binary(Node):
    op: str
    left: Node
    right: Node


_INFIX_BP: dict[str, tuple[int, int]] = {
    "|": (10, 11),
    "&": (20, 21),
    "==": (30, 31), "!=": (30, 31), "~=": (30, 31),
    ">": (30, 31), "<": (30, 31), ">=": (30, 31), "<=": (30, 31),
    "+": (40, 41), "-": (40, 41),
    "*": (50, 51), "/": (50, 51),
    "^": (70, 71),
}
_UMINUS_BP = 60
_NOT_BP = 80

# funções cujo argumento é um nome, não uma expressão: r(mean), e(N), c(pi)
_RAW_ARG_FUNCS = {"r", "e", "c", "s"}
# _b[x], _se[x]: índice é um nome
_RAW_SUBSCRIPTS = {"_b", "_se", "_coef"}


class Parser:
    def __init__(self, text: str, *, stop_on_unknown: bool = False):
        self.text = text
        self.toks = tokenize(text, stop_on_unknown=stop_on_unknown)
        self.i = 0

    # -- utilidades ----------------------------------------------------
    @property
    def tok(self) -> Tok:
        return self.toks[self.i]

    def advance(self) -> Tok:
        t = self.toks[self.i]
        self.i += 1
        return t

    def expect_op(self, op: str) -> None:
        t = self.advance()
        if t.kind != "op" or t.value != op:
            raise syntax_error()

    def at_end(self) -> bool:
        return self.tok.kind == "end"

    @property
    def offset(self) -> int:
        """Posição no texto do próximo token não consumido."""
        return self.tok.pos

    # -- gramática -------------------------------------------------------
    def parse_full(self) -> Node:
        node = self.parse_expr()
        if not self.at_end():
            raise syntax_error()
        return node

    def parse_expr(self, min_bp: int = 0) -> Node:
        left = self.parse_prefix()
        while True:
            t = self.tok
            if t.kind != "op" or t.value not in _INFIX_BP:
                if t.kind == "op" and t.value == "=":
                    raise syntax_error()
                break
            lbp, rbp = _INFIX_BP[t.value]
            if lbp < min_bp:
                break
            self.advance()
            right = self.parse_expr(rbp)
            op = "!=" if t.value == "~=" else t.value
            left = Binary(op, left, right)
        return left

    def parse_prefix(self) -> Node:
        t = self.advance()
        if t.kind == "num":
            return Num(t.value)
        if t.kind == "str":
            return Str(t.value)
        if t.kind == "op":
            if t.value == "(":
                node = self.parse_expr()
                self.expect_op(")")
                return node
            if t.value == "-":
                return Unary("-", self.parse_expr(_UMINUS_BP))
            if t.value == "+":
                return Unary("+", self.parse_expr(_UMINUS_BP))
            if t.value in ("!", "~"):
                return Unary("!", self.parse_expr(_NOT_BP))
            raise syntax_error()
        if t.kind == "name":
            return self.parse_name(t)
        raise syntax_error()

    def parse_name(self, t: Tok) -> Node:
        name = t.value
        nxt = self.tok
        if nxt.kind == "op" and nxt.value == "(":
            if name in _RAW_ARG_FUNCS:
                self.advance()
                start = self.tok.pos
                depth = 1
                while depth:
                    tk = self.advance()
                    if tk.kind == "end":
                        raise StataError(132, "too many '(' or '['")
                    if tk.kind == "op" and tk.value == "(":
                        depth += 1
                    elif tk.kind == "op" and tk.value == ")":
                        depth -= 1
                end = self.toks[self.i - 1].pos
                return Call(name, [], raw=self.text[start:end].strip())
            self.advance()
            args: list[Node] = []
            if not (self.tok.kind == "op" and self.tok.value == ")"):
                while True:
                    args.append(self.parse_expr())
                    if self.tok.kind == "op" and self.tok.value == ",":
                        self.advance()
                        continue
                    break
            self.expect_op(")")
            return Call(name, args)
        if nxt.kind == "op" and nxt.value == "[":
            self.advance()
            if name in _RAW_SUBSCRIPTS:
                start = self.tok.pos
                while not (self.tok.kind == "op" and self.tok.value == "]"):
                    if self.tok.kind == "end":
                        raise syntax_error()
                    self.advance()
                raw = self.text[start:self.tok.pos].strip()
                self.advance()
                return Subscript(name, raw)
            idx = self.parse_expr()
            if self.tok.kind == "op" and self.tok.value == ",":
                self.advance()
                idx2 = self.parse_expr()
                self.expect_op("]")
                return Subscript(name, idx, idx2)
            self.expect_op("]")
            return Subscript(name, idx)
        return Name(name)


# funções cujo 1º argumento é o nome de uma matriz (ver commands/matrix.py)
MATRIX_FUNCS = {"rowsof", "colsof", "el", "trace", "det", "issymmetric", "rownumb",
                "colnumb", "matmissing", "mreldif"}


def parse(text: str) -> Node:
    return Parser(text).parse_full()


# ---------------------------------------------------------------------------
# Avaliação
# ---------------------------------------------------------------------------


class Context(Protocol):
    def resolve_name(self, name: str) -> Value: ...
    def resolve_subscript(self, name: str, index: Any) -> Value: ...
    def resolve_result(self, kind: str, raw: str) -> Value: ...
    def call_function(self, name: str, args: list[Value]) -> Value: ...


def _num(x: float) -> float:
    return M.normalize(x)


def _truth(v: Value) -> bool:
    if isinstance(v, str):
        raise type_mismatch()
    return v != 0


def evaluate(node: Node, ctx: Context) -> Value:
    if isinstance(node, Num):
        return node.value
    if isinstance(node, Str):
        return node.value
    if isinstance(node, Name):
        return ctx.resolve_name(node.name)
    if isinstance(node, Call):
        if node.name in _RAW_ARG_FUNCS:
            return ctx.resolve_result(node.name, node.raw)
        if node.name in MATRIX_FUNCS:
            from ..commands.matrix import matrix_scalar
            return matrix_scalar(ctx.s, node.name, node.args, lambda a: evaluate(a, ctx))
        return ctx.call_function(node.name, [evaluate(a, ctx) for a in node.args])
    if isinstance(node, Subscript):
        if node.index2 is not None:
            from ..commands.matrix import matrix_element
            return matrix_element(ctx.s, node.name, evaluate(node.index, ctx),
                                  evaluate(node.index2, ctx))
        idx = node.index if isinstance(node.index, str) else evaluate(node.index, ctx)
        return ctx.resolve_subscript(node.name, idx)
    if isinstance(node, Unary):
        v = evaluate(node.operand, ctx)
        if node.op == "!":
            return 0.0 if _truth(v) else 1.0
        if node.op == "+":
            if isinstance(v, str):
                # display "a" + "b" termina em erro r(198) no Stata
                raise StataError(198, "invalid syntax")
            return v
        if isinstance(v, str):
            raise type_mismatch()
        return M.SYSMISS if M.is_missing(v) else _num(-v)
    if isinstance(node, Binary):
        return _binary(node.op, evaluate(node.left, ctx), evaluate(node.right, ctx))
    raise TypeError(node)


def _binary(op: str, a: Value, b: Value) -> Value:
    if op in ("&", "|"):
        ta, tb = _truth(a), _truth(b)
        return float(ta and tb) if op == "&" else float(ta or tb)

    a_str, b_str = isinstance(a, str), isinstance(b, str)

    # "ab"*3 repete a string
    if op == "*" and a_str != b_str:
        s, n = (a, b) if a_str else (b, a)
        return "" if M.is_missing(n) or n < 1 else s * int(n)

    if a_str != b_str:
        raise type_mismatch()

    if op in _RELATIONAL:
        return 1.0 if _RELATIONAL[op](a, b) else 0.0

    if a_str:
        if op == "+":
            return a + b
        raise type_mismatch()

    if M.is_missing(a) or M.is_missing(b):
        return M.SYSMISS
    try:
        if op == "+":
            return _num(a + b)
        if op == "-":
            return _num(a - b)
        if op == "*":
            return _num(a * b)
        if op == "/":
            return M.SYSMISS if b == 0 else _num(a / b)
        if op == "^":
            return _power(a, b)
    except (OverflowError, ValueError, ZeroDivisionError):
        return M.SYSMISS
    raise syntax_error()


_RELATIONAL: dict[str, Callable[[Any, Any], bool]] = {
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
    ">": lambda a, b: a > b,
    "<": lambda a, b: a < b,
    ">=": lambda a, b: a >= b,
    "<=": lambda a, b: a <= b,
}


def _power(a: float, b: float) -> float:
    if a == 0 and b == 0:
        return 1.0
    if a < 0 and b != int(b):
        return M.SYSMISS
    if a == 0 and b < 0:
        return M.SYSMISS
    return _num(math.pow(a, b))


class SimpleContext:
    """Contexto mínimo para testes: só constantes e funções."""

    def __init__(self, functions: Callable[[str, list[Value]], Value] | None = None):
        self._functions = functions

    def resolve_name(self, name: str) -> Value:
        if name == "_pi":
            return math.pi
        raise StataError(111, f"{name} not found")

    def resolve_subscript(self, name: str, index: Any) -> Value:
        raise StataError(111, f"{name} not found")

    def resolve_result(self, kind: str, raw: str) -> Value:
        return M.SYSMISS

    def call_function(self, name: str, args: list[Value]) -> Value:
        if self._functions is None:
            raise StataError(133, "unknown function " + name + "()")
        return self._functions(name, args)
