"""Avaliação vetorizada de expressões sobre as observações do conjunto de dados.

A árvore sintática é a mesma de lang/expr.py. Aqui cada nó devolve um
escalar (float ou str) ou um vetor numpy com uma posição por observação:
float64 para números (com os códigos de missing do OpenDTA) e object para
strings.

Grupos `by`: `_n`, `_N` e subscritos como `x[_n-1]` são relativos ao grupo
da observação. Sem `by`, há um único grupo com todas as observações.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import numpy as np

from ..core import missing as M
from ..core.errors import StataError, type_mismatch
from . import functions as F
from .expr import (MATRIX_FUNCS, Binary, Call, Name, Node, Num, Str, Subscript, Unary,
                   _RAW_ARG_FUNCS, _binary)

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# Grupos
# ---------------------------------------------------------------------------

class Groups:
    def __init__(self, start: np.ndarray, size: np.ndarray):
        self.start = start            # índice da 1ª obs do grupo, por observação
        self.size = size              # tamanho do grupo, por observação
        n = len(start)
        self.pos = np.arange(n) - start   # posição 0-based no grupo

    @classmethod
    def single(cls, nobs: int) -> "Groups":
        return cls(np.zeros(nobs, dtype=np.int64), np.full(nobs, nobs, dtype=np.int64))

    @classmethod
    def from_columns(cls, columns: list[np.ndarray]) -> "Groups":
        n = len(columns[0]) if columns else 0
        if n == 0:
            return cls(np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64))
        change = np.zeros(n, dtype=bool)
        change[0] = True
        for col in columns:
            change[1:] |= col[1:] != col[:-1]
        starts = np.flatnonzero(change)
        ends = np.append(starts[1:], n)
        start = np.repeat(starts, ends - starts)
        size = np.repeat(ends - starts, ends - starts)
        return cls(start.astype(np.int64), size.astype(np.int64))

    @property
    def first(self) -> np.ndarray:
        return self.pos == 0

    @property
    def ids(self) -> np.ndarray:
        return np.cumsum(self.first) - 1


# ---------------------------------------------------------------------------
# Utilidades de tipo
# ---------------------------------------------------------------------------

def is_str(v: Any) -> bool:
    return isinstance(v, str) or (isinstance(v, np.ndarray) and v.dtype == object)


def is_vec(v: Any) -> bool:
    return isinstance(v, np.ndarray)


def broadcast(v: Any, n: int) -> np.ndarray:
    if is_vec(v):
        return v
    if isinstance(v, str):
        return np.array([v] * n, dtype=object)
    return np.full(n, float(v))


def _num_arr(v: Any) -> np.ndarray | float:
    return v if is_vec(v) else float(v)


def _clean(r: np.ndarray, miss: np.ndarray | bool) -> np.ndarray:
    """Missing nas entradas, divisões por zero e estouros viram `.`."""
    r = np.where(miss, SYS, r)
    bad = ~np.isfinite(r) | (np.abs(r) >= SYS)
    return np.where(bad, SYS, r)


def truth(v: Any) -> Any:
    if is_str(v):
        raise type_mismatch()
    return v != 0


# ---------------------------------------------------------------------------
# Contexto
# ---------------------------------------------------------------------------

class VectorContext:
    def __init__(self, session: "Session", groups: Groups | None = None):
        self.s = session
        self.ds = session.data
        self.n = self.ds.nobs
        self.groups = groups or Groups.single(self.n)

    def name(self, name: str) -> Any:
        ds = self.ds
        if ds.has(name):
            return ds.get(name).data
        if name in self.s.scalars:
            return self.s.scalars[name]
        if name == "_n":
            return (self.groups.pos + 1).astype(np.float64)
        if name == "_N":
            return self.groups.size.astype(np.float64)
        if name == "_pi":
            return math.pi
        if name == "_rc":
            return float(self.s.rc)
        if self.s.settings.get("varabbrev", "on") == "on":
            matches = [n for n in ds.names if n.startswith(name)]
            if len(matches) == 1:
                return ds.get(matches[0]).data
            if len(matches) > 1:
                raise StataError(111, f"{name} ambiguous abbreviation")
        raise StataError(111, f"{name} not found")

    def subscript(self, name: str, idx: Any) -> Any:
        if name in ("_b", "_se", "_coef"):
            return self.s.context.resolve_subscript(name, idx)
        data = self.name(name)
        if not is_vec(data):
            raise StataError(198, "invalid syntax")
        if is_str(idx):
            raise type_mismatch()
        k = np.trunc(broadcast(idx, self.n))
        g = self.groups
        valid = (k < SYS) & (k >= 1) & (k <= g.size)
        target = np.where(valid, g.start + k.clip(1, None) - 1, 0).astype(np.int64)
        if self.n == 0:
            return data[:0]
        vals = data[target]
        if data.dtype == object:
            return np.where(valid, vals, "").astype(object)
        return np.where(valid, vals, SYS)

    def call(self, name: str, args: list[Any]) -> Any:
        if name == "sum" and len(args) == 1:
            return self._running_sum(args[0])
        if name in _RANDOM:
            return self._random(name, args)
        if not any(is_vec(a) for a in args):
            return F.call(name, args)
        fast = _FAST.get(name)
        if fast is not None:
            out = fast(*args)
            if out is not NotImplemented:
                return out
        if name not in F.FUNCTIONS:
            raise StataError(133, f"unknown function {name}()")
        cols = [broadcast(a, self.n) for a in args]
        results = [F.call(name, [c[i] for c in cols]) for i in range(self.n)]
        if any(isinstance(r, str) for r in results):
            return np.array(results, dtype=object)
        return np.array(results, dtype=np.float64)

    def _random(self, name: str, args: list[Any]) -> np.ndarray:
        """Um sorteio por observação (só as da amostra em generate/replace)."""
        from ..core import rng
        lo, hi, _ = rng.DISTRIBUTIONS[name]
        if not (lo <= len(args) <= hi):
            raise StataError(198, "invalid syntax")
        if any(is_str(a) for a in args):
            raise type_mismatch()
        mask = getattr(self.s, "_rng_mask", None)
        if mask is None or len(mask) != self.n:
            return rng.draw(name, [broadcast(a, self.n) if is_vec(a) else a for a in args], self.n)
        # só as observações da amostra sorteiam (generate/replace com if/in)
        idx = np.flatnonzero(mask)
        out = np.full(self.n, SYS)
        out[idx] = rng.draw(name, [broadcast(a, self.n)[idx] if is_vec(a) else a for a in args], len(idx))
        return out

    def _running_sum(self, x: Any) -> np.ndarray:
        if is_str(x):
            raise type_mismatch()
        v = broadcast(x, self.n).astype(np.float64)
        v = np.where(v >= SYS, 0.0, v)
        total = np.cumsum(v)
        base = np.where(self.groups.start > 0, total[self.groups.start - 1], 0.0) if self.n else total
        return total - base


# ---------------------------------------------------------------------------
# Avaliação
# ---------------------------------------------------------------------------

def evaluate_vec(node: Node, ctx: VectorContext) -> Any:
    if isinstance(node, Num):
        return node.value
    if isinstance(node, Str):
        return node.value
    if isinstance(node, Name):
        return ctx.name(node.name)
    if isinstance(node, Call):
        if node.name in _RAW_ARG_FUNCS:
            return ctx.s.context.resolve_result(node.name, node.raw)
        if node.name in MATRIX_FUNCS:
            from ..commands.matrix import matrix_scalar
            return matrix_scalar(ctx.s, node.name, node.args, lambda a: evaluate_vec(a, ctx))
        return ctx.call(node.name, [evaluate_vec(a, ctx) for a in node.args])
    if isinstance(node, Subscript):
        if node.index2 is not None:
            from ..commands.matrix import matrix_element
            return matrix_element(ctx.s, node.name, evaluate_vec(node.index, ctx),
                                  evaluate_vec(node.index2, ctx))
        idx = node.index if isinstance(node.index, str) else evaluate_vec(node.index, ctx)
        return ctx.subscript(node.name, idx)
    if isinstance(node, Unary):
        v = evaluate_vec(node.operand, ctx)
        if node.op == "!":
            t = truth(v)
            return (~t).astype(np.float64) if is_vec(t) else (0.0 if t else 1.0)
        if node.op == "+":
            if is_str(v):
                raise StataError(198, "invalid syntax")
            return v
        if is_str(v):
            raise type_mismatch()
        if not is_vec(v):
            return SYS if M.is_missing(v) else M.normalize(-v)
        return np.where(v >= SYS, SYS, -v)
    if isinstance(node, Binary):
        return binary_vec(node.op, evaluate_vec(node.left, ctx), evaluate_vec(node.right, ctx))
    raise TypeError(node)


def binary_vec(op: str, a: Any, b: Any) -> Any:
    if not is_vec(a) and not is_vec(b):
        return _binary(op, a, b)

    if op in ("&", "|"):
        ta, tb = truth(a), truth(b)
        r = (ta & tb) if op == "&" else (ta | tb)
        return np.asarray(r, dtype=np.float64)

    sa, sb = is_str(a), is_str(b)
    if op == "*" and sa != sb:
        s, k = (a, b) if sa else (b, a)
        n = len(s) if is_vec(s) else len(k)
        ss, kk = broadcast(s, n), broadcast(k, n)
        return np.array(["" if (x >= SYS or x < 1) else t * int(x) for t, x in zip(ss, kk)], dtype=object)
    if sa != sb:
        raise type_mismatch()

    if sa:
        n = len(a) if is_vec(a) else len(b)
        aa, bb = broadcast(a, n), broadcast(b, n)
        if op == "+":
            return np.array([x + y for x, y in zip(aa, bb)], dtype=object)
        cmp = {"==": lambda x, y: x == y, "!=": lambda x, y: x != y, ">": lambda x, y: x > y,
               "<": lambda x, y: x < y, ">=": lambda x, y: x >= y, "<=": lambda x, y: x <= y}.get(op)
        if cmp is None:
            raise type_mismatch()
        return np.array([1.0 if cmp(x, y) else 0.0 for x, y in zip(aa, bb)])

    a, b = _num_arr(a), _num_arr(b)
    if op in ("==", "!=", ">", "<", ">=", "<="):
        r = {"==": np.equal, "!=": np.not_equal, ">": np.greater, "<": np.less,
             ">=": np.greater_equal, "<=": np.less_equal}[op](a, b)
        return np.asarray(r, dtype=np.float64)

    miss = (np.asarray(a) >= SYS) | (np.asarray(b) >= SYS)
    with np.errstate(all="ignore"):
        if op == "+":
            r = np.add(a, b)
        elif op == "-":
            r = np.subtract(a, b)
        elif op == "*":
            r = np.multiply(a, b)
        elif op == "/":
            r = np.divide(a, b)
            miss = miss | (np.asarray(b) == 0)
        elif op == "^":
            aa, bb = np.broadcast_arrays(np.asarray(a, dtype=float), np.asarray(b, dtype=float))
            r = np.power(aa, bb)
            miss = miss | ((aa < 0) & (bb != np.trunc(bb))) | ((aa == 0) & (bb < 0))
        else:
            raise StataError(198, "invalid syntax")
    return _clean(np.asarray(r, dtype=np.float64), miss)


# ---------------------------------------------------------------------------
# Funções vetorizadas mais comuns (as demais caem no caminho elemento a elemento)
# ---------------------------------------------------------------------------

def _unary_num(fn, domain=None):
    def f(x):
        if is_str(x):
            raise type_mismatch()
        x = np.asarray(x, dtype=np.float64)
        miss = x >= SYS
        if domain is not None:
            miss = miss | ~domain(np.where(miss, 1.0, x))
        with np.errstate(all="ignore"):
            r = fn(np.where(miss, 1.0, x))
        return _clean(np.asarray(r, dtype=np.float64), miss)
    return f


def _round(x, y=1.0):
    if is_str(x) or is_str(y):
        raise type_mismatch()
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    miss = (x >= SYS) | (y >= SYS)
    with np.errstate(all="ignore"):
        safe = np.where(y == 0, 1.0, y)
        r = np.where(y == 0, x, np.floor(x / safe + 0.5) * safe)
    return _clean(r, miss)


def _mod(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    miss = (x >= SYS) | (y >= SYS) | (y == 0)
    with np.errstate(all="ignore"):
        r = x - y * np.floor(x / np.where(y == 0, 1, y))
    return _clean(r, miss)


def _missing(*xs):
    out = None
    for x in xs:
        m = (np.asarray(x, dtype=object) == "") if is_str(x) else (np.asarray(x) >= SYS)
        out = m if out is None else (out | m)
    return np.asarray(out, dtype=np.float64)


def _minmax(fn):
    def f(*xs):
        if any(is_str(x) for x in xs):
            raise type_mismatch()
        n = max(len(x) for x in xs if is_vec(x))
        stack = np.vstack([broadcast(x, n) for x in xs])
        masked = np.where(stack >= SYS, np.nan, stack)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                r = fn(masked, axis=0)
        return np.where(np.isnan(r), SYS, r)
    return f


def _cond(x, a, b, c=None):
    if is_str(x):
        raise type_mismatch()
    n = max(len(v) for v in (x, a, b, c) if is_vec(v))
    xv = broadcast(x, n)
    out = np.where(xv != 0, broadcast(a, n), broadcast(b, n))
    if c is not None:
        out = np.where(xv >= SYS, broadcast(c, n), out)
    return out.astype(object) if is_str(a) else out.astype(np.float64)


def _inlist(z, *items):
    if any(is_str(i) != is_str(z) for i in items):
        raise type_mismatch()
    zz = np.asarray(z, dtype=object if is_str(z) else np.float64)
    r = np.zeros(zz.shape, dtype=bool)
    for it in items:
        r |= zz == np.asarray(it, dtype=zz.dtype)
    return r.astype(np.float64)


def _inrange(z, a, b):
    if is_str(z):
        return NotImplemented
    z = np.asarray(z, dtype=np.float64)
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    lo = np.where(a >= SYS, -np.inf, a)
    r = (z >= lo) & (z <= b) & ((z < SYS) | (b >= SYS))
    return r.astype(np.float64)


def _float(x):
    x = np.asarray(x, dtype=np.float64)
    miss = x >= SYS
    return np.where(miss, x, x.astype(np.float32).astype(np.float64))


_RANDOM = {"runiform", "runiformint", "rnormal", "rbinomial", "rpoisson", "rchi2", "rt", "rbeta",
           "rgamma", "rexponential", "rlogistic", "rweibull", "rnbinomial", "rhypergeometric"}

_FAST = {
    "abs": _unary_num(np.abs),
    "sqrt": _unary_num(np.sqrt, lambda x: x >= 0),
    "exp": _unary_num(np.exp),
    "ln": _unary_num(np.log, lambda x: x > 0),
    "log": _unary_num(np.log, lambda x: x > 0),
    "log10": _unary_num(np.log10, lambda x: x > 0),
    "int": _unary_num(np.trunc),
    "trunc": _unary_num(np.trunc),
    "floor": _unary_num(np.floor),
    "ceil": _unary_num(np.ceil),
    "round": _round,
    "mod": _mod,
    "missing": _missing,
    "mi": _missing,
    "min": _minmax(np.nanmin),
    "max": _minmax(np.nanmax),
    "cond": _cond,
    "inlist": _inlist,
    "inrange": _inrange,
    "float": _float,
}
