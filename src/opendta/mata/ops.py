"""Operadores do Mata ([M-2] op_arith, op_colon, op_join, op_range...).

Aritmética com missing: qualquer operando missing dá "."; resultados não
finitos (divisão por zero, overflow) também viram ".".
"""

from __future__ import annotations

import numpy as np

from ..core import missing as M
from .values import MV, MataError, SYS, bool_mv, clean, conformability, real, type_mismatch


def _numeric(*vals: MV) -> None:
    for v in vals:
        if v.t not in ("real", "complex"):
            raise MataError(3251, "nonnumeric found where numeric required")


def _same_type(a: MV, b: MV) -> None:
    if (a.t == "string") != (b.t == "string"):
        raise type_mismatch()


def _arith(fn, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    miss = (a >= SYS) | (b >= SYS)
    with np.errstate(all="ignore"):
        an = np.where(a >= SYS, 0.0, a)
        bn = np.where(b >= SYS, 0.0, b)
        r = fn(an, bn)
    r = clean(r)
    return np.where(miss, SYS, r)


def _c_conformable(a: MV, b: MV) -> None:
    for x, y in zip(a.a.shape, b.a.shape):
        if x != y and x != 1 and y != 1:
            raise conformability()


def _plain_conformable(a: MV, b: MV) -> None:
    """+ e - exigem o mesmo tamanho (A + 10 é erro; use :+). Observado no Stata 14."""
    if a.a.shape != b.a.shape:
        raise conformability()


def _scalar_or_same(a: MV, b: MV) -> None:
    if a.a.shape != b.a.shape and not a.is_scalar and not b.is_scalar:
        raise conformability()


_REAL_FN = {
    "+": np.add, "-": np.subtract, "*": np.multiply, "/": np.divide,
}


def _elementwise(op: str, a: MV, b: MV) -> MV:
    if a.t == "string" or b.t == "string":
        _same_type(a, b)
        if op == "+":
            r = np.frompyfunc(lambda x, y: str(x) + str(y), 2, 1)(a.a, b.a)
            return MV(np.asarray(r, dtype=object).reshape(np.broadcast(a.a, b.a).shape), "string")
        raise type_mismatch()
    _numeric(a, b)
    if a.t == "complex" or b.t == "complex":
        fn = {"+": np.add, "-": np.subtract, "*": np.multiply, "/": np.divide, "^": np.power}[op]
        with np.errstate(all="ignore"):
            return MV(fn(a.a.astype(np.complex128), b.a.astype(np.complex128)), "complex")
    if op == "^":
        def pw(x, y):
            r = np.power(x, y)
            return np.where((x < 0) & (y != np.trunc(y)), np.nan, r)
        return real(_arith(pw, a.a, b.a))
    return real(_arith(_REAL_FN[op], a.a, b.a))


def binary(op: str, a: MV, b: MV) -> MV:
    if op in (":+", ":-", ":*", ":/", ":^"):
        _c_conformable(a, b)
        return _elementwise(op[1], a, b)
    if op in ("+", "-"):
        _plain_conformable(a, b)
        return _elementwise(op, a, b)
    if op == "*":
        if a.t == "string" or b.t == "string":
            return _repeat(a, b)
        _numeric(a, b)
        if a.is_scalar or b.is_scalar:
            return _elementwise("*", a, b)
        return matmul(a, b)
    if op == "/":
        if not b.is_scalar:
            raise conformability()
        return _elementwise("/", a, b)
    if op == "^":
        _numeric(a, b)
        if a.is_scalar and b.is_scalar:
            return _elementwise("^", a, b)
        if b.is_scalar and a.rows == a.cols:
            k = b.real_scalar()
            if k == int(k) and k >= 0:
                return real(clean(np.linalg.matrix_power(np.where(a.a >= SYS, np.nan, a.a), int(k))))
        raise conformability()
    if op == "#":
        _numeric(a, b)
        return real(_arith(lambda x, y: x * y, np.kron(a.a, np.ones(b.a.shape)),
                           np.kron(np.ones(a.a.shape), b.a)))
    if op == ",":
        return hjoin(a, b)
    if op == "\\":
        return vjoin(a, b)
    if op in ("..", "::"):
        return seq(op, a, b)
    if op in ("==", "!="):
        eq = (a.t == b.t or (a.t in ("real", "complex") and b.t in ("real", "complex"))) \
            and a.a.shape == b.a.shape and bool(np.all(a.a == b.a))
        return bool_mv(eq if op == "==" else not eq)
    if op in (">", ">=", "<", "<="):
        _scalar_or_same(a, b)
        _same_type(a, b)
        r = _compare(op, a.a, b.a)
        return bool_mv(bool(np.all(r)))
    if op in (":==", ":!=", ":>", ":>=", ":<", ":<="):
        _c_conformable(a, b)
        _same_type(a, b)
        return real(_compare(op[1:], a.a, b.a).astype(np.float64))
    if op in ("&", "|"):
        x, y = _truth_scalar(a), _truth_scalar(b)
        return bool_mv((x and y) if op == "&" else (x or y))
    if op in (":&", ":|"):
        _c_conformable(a, b)
        _numeric(a, b)
        x, y = a.a != 0, b.a != 0
        return real((x & y if op == ":&" else x | y).astype(np.float64))
    raise MataError(3000, f"operator {op} not supported")


def _truth_scalar(v: MV) -> bool:
    if v.t != "real":
        raise MataError(3253, "nonreal found where real required")
    return v.scalar() != 0


def _compare(op: str, a, b) -> np.ndarray:
    import operator
    fn = {"==": operator.eq, "!=": operator.ne, ">": operator.gt, ">=": operator.ge,
          "<": operator.lt, "<=": operator.le}[op]
    return np.asarray(np.frompyfunc(fn, 2, 1)(a, b), dtype=bool) if a.dtype == object or \
        getattr(b, "dtype", None) == object else fn(a, b)


def _repeat(a: MV, b: MV) -> MV:
    s, n = (a, b) if a.t == "string" else (b, a)
    if n.t != "real":
        raise type_mismatch()
    _scalar_or_same(s, n)
    k = np.broadcast_to(n.a, np.broadcast(s.a, n.a).shape)
    ss = np.broadcast_to(s.a, k.shape)
    out = np.empty(k.shape, dtype=object)
    for idx in np.ndindex(k.shape):
        reps = k[idx]
        out[idx] = "" if reps >= SYS or reps <= 0 else str(ss[idx]) * int(reps)
    return MV(out, "string")


def matmul(a: MV, b: MV) -> MV:
    if a.cols != b.rows:
        raise conformability()
    if a.t == "complex" or b.t == "complex":
        return MV(a.a.astype(np.complex128) @ b.a.astype(np.complex128), "complex")
    x = np.where(a.a >= SYS, np.nan, a.a)
    y = np.where(b.a >= SYS, np.nan, b.a)
    with np.errstate(all="ignore"):
        return real(clean(x @ y))


def unary(op: str, v: MV) -> MV:
    if op == "-":
        _numeric(v)
        if v.t == "complex":
            return MV(-v.a, "complex")
        return real(np.where(v.a >= SYS, SYS, -v.a))
    if op == "!":
        _numeric(v)
        return real((v.a == 0).astype(np.float64))
    raise MataError(3000, f"operator {op} not supported")


def transpose(v: MV) -> MV:
    if v.t == "complex":
        return MV(np.conj(v.a.T), "complex")
    return MV(v.a.T.copy(), v.t)


def _join_types(a: MV, b: MV) -> str:
    if a.a.size == 0:
        return b.t
    if b.a.size == 0:
        return a.t
    if a.t == b.t:
        return a.t
    if {a.t, b.t} == {"real", "complex"}:
        return "complex"
    raise type_mismatch()


def hjoin(a: MV, b: MV) -> MV:
    t = _join_types(a, b)
    if a.cols == 0 and a.rows in (0, b.rows):
        return MV(b.a.copy(), t if b.a.size else b.t)
    if b.cols == 0 and b.rows in (0, a.rows):
        return MV(a.a.copy(), t if a.a.size else a.t)
    if a.rows != b.rows:
        raise conformability()
    return MV(np.hstack([a.a.astype(_dtype(t)), b.a.astype(_dtype(t))]), t)


def vjoin(a: MV, b: MV) -> MV:
    t = _join_types(a, b)
    if a.rows == 0 and a.cols in (0, b.cols):
        return MV(b.a.copy(), t if b.a.size else b.t)
    if b.rows == 0 and b.cols in (0, a.cols):
        return MV(a.a.copy(), t if a.a.size else a.t)
    if a.cols != b.cols:
        raise conformability()
    return MV(np.vstack([a.a.astype(_dtype(t)), b.a.astype(_dtype(t))]), t)


def _dtype(t: str):
    return object if t in ("string", "pointer", "struct") else (np.complex128 if t == "complex" else np.float64)


def seq(op: str, a: MV, b: MV) -> MV:
    x, y = a.real_scalar(), b.real_scalar()
    if x >= SYS or y >= SYS:
        raise MataError(3300, "argument out of range")
    step = 1.0 if y >= x else -1.0
    n = int(np.floor(abs(y - x))) + 1
    vals = x + step * np.arange(n)
    return real(vals.reshape(1, -1) if op == ".." else vals.reshape(-1, 1))


def missing_name(x: float) -> str:
    return M.missing_name(x)
