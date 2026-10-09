"""Funções auxiliares dos comandos de dados: if/in, avaliação vetorizada,
contexto por observação (para replace sequencial) e mensagens."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang import functions as F
from ..lang.expr import Node, Subscript, Unary, Binary, Call, parse
from ..lang.syntax import Parsed, parse_in
from ..lang.vexpr import Groups, VectorContext, broadcast, evaluate_vec, is_str, truth

if TYPE_CHECKING:
    from ..session import Session


def groups(s: "Session") -> Groups:
    return s.by_groups if s.by_groups is not None else Groups.single(s.data.nobs)


def eval_vector(s: "Session", text: str) -> Any:
    return evaluate_vec(parse(text), VectorContext(s, groups(s)))


def eval_sample(s: "Session", text: str, mask: np.ndarray) -> Any:
    """Como eval_vector, mas as funções aleatórias só sorteiam para as
    observações da amostra (if/in), na ordem delas, como o Stata: assim
    `gen u = runiform() if x` consome a mesma sequência que no Stata."""
    s._rng_mask = mask
    try:
        return eval_vector(s, text)
    finally:
        s._rng_mask = None


def touse(s: "Session", p: Parsed) -> np.ndarray:
    n = s.data.nobs
    mask = np.ones(n, dtype=bool)
    if p.in_:
        a, b = parse_in(p.in_, n)
        inmask = np.zeros(n, dtype=bool)
        inmask[a:b] = True
        mask &= inmask
    if p.if_:
        v = eval_vector(s, p.if_)
        t = truth(v)
        mask &= np.asarray(broadcast(t, n) if not isinstance(t, np.ndarray) else t, dtype=bool)
    return mask


def plural(n: int, word: str, plural_word: str | None = None) -> str:
    return f"{n:,} {word if n == 1 else (plural_word or word + 's')}"


def references_subscript(node: Node, name: str) -> bool:
    """A expressão usa name[...]? (replace então precisa ser sequencial)."""
    if isinstance(node, Subscript):
        return node.name == name or (not isinstance(node.index, str) and references_subscript(node.index, name))
    if isinstance(node, Unary):
        return references_subscript(node.operand, name)
    if isinstance(node, Binary):
        return references_subscript(node.left, name) or references_subscript(node.right, name)
    if isinstance(node, Call):
        return any(references_subscript(a, name) for a in node.args)
    return False


class ObsContext:
    """Contexto escalar posicionado numa observação (replace sequencial)."""

    def __init__(self, s: "Session", g: Groups):
        self.s = s
        self.g = g
        self.i = 0

    def resolve_name(self, name: str) -> Any:
        ds = self.s.data
        if ds.has(name):
            return ds.get(name).value(self.i)
        if name in self.s.scalars:
            return self.s.scalars[name]
        if name == "_n":
            return float(self.g.pos[self.i] + 1)
        if name == "_N":
            return float(self.g.size[self.i])
        return self.s.context.resolve_name(name)

    def resolve_subscript(self, name: str, index: Any) -> Any:
        ds = self.s.data
        if not ds.has(name):
            return self.s.context.resolve_subscript(name, index)
        var = ds.get(name)
        if isinstance(index, str):
            raise StataError(109, "type mismatch")
        k = int(index) if not M.is_missing(index) else 0
        if 1 <= k <= self.g.size[self.i]:
            return var.value(self.g.start[self.i] + k - 1)
        return "" if var.is_string else M.SYSMISS

    def resolve_result(self, kind: str, raw: str) -> Any:
        return self.s.context.resolve_result(kind, raw)

    def call_function(self, name: str, args: list[Any]) -> Any:
        return F.call(name, args)


def check_kind(value: Any, want_string: bool) -> None:
    if is_str(value) != want_string:
        raise StataError(109, "type mismatch")
