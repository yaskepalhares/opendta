"""Valores do Mata.

Toda matriz é um MV: um vetor numpy de duas dimensões e o tipo dos
elementos (eltype): real (float64, com os missing do Stata), complex,
string (objetos str), pointer (Pointer ou None) e struct (StructVal).
Um escalar é uma matriz 1×1.
"""

from __future__ import annotations

import numpy as np

from ..core import missing as M

SYS = M.SYSMISS


class MataError(Exception):
    """Erro de execução do Mata: código (3200, 3499...) e mensagem."""

    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code
        self.message = message
        self.where: list[str] = []      # pilha de funções (para a mensagem)


def conformability() -> MataError:
    return MataError(3200, "conformability error")


def type_mismatch(msg: str = "type mismatch") -> MataError:
    return MataError(3250, msg)


class Pointer:
    """Ponteiro para uma variável (frame, nome), um elemento de struct ou uma
    função (&f())."""

    def __init__(self, getter, setter=None, label: str = ""):
        self.get = getter
        self.set = setter
        self.label = label

    def __repr__(self) -> str:
        return f"<pointer {self.label}>"


class StructVal:
    def __init__(self, sname: str, fields: dict):
        self.sname = sname
        self.fields = fields


class MV:
    __slots__ = ("a", "t")

    def __init__(self, a, t: str = "real"):
        arr = np.asarray(a, dtype=object if t in ("string", "pointer", "struct") else
                         (np.complex128 if t == "complex" else np.float64))
        if arr.ndim == 0:
            arr = arr.reshape(1, 1)
        elif arr.ndim == 1:
            arr = arr.reshape(1, -1)
        self.a = arr
        self.t = t

    # tamanho e organização ----------------------------------------------
    @property
    def rows(self) -> int:
        return int(self.a.shape[0])

    @property
    def cols(self) -> int:
        return int(self.a.shape[1])

    @property
    def is_scalar(self) -> bool:
        return self.a.shape == (1, 1)

    def orgtype(self) -> str:
        r, c = self.a.shape
        if r == 1 and c == 1:
            return "scalar"
        if r == 1:
            return "rowvector"
        if c == 1:
            return "colvector"
        return "matrix"

    def copy(self) -> "MV":
        return MV(self.a.copy(), self.t)

    def __repr__(self) -> str:
        return f"MV({self.t}, {self.a.tolist()})"

    # conversões -------------------------------------------------------------
    def scalar(self):
        if not self.is_scalar:
            raise MataError(3204, f"{self.orgtype()} found where scalar required")
        return self.a[0, 0]

    def real_scalar(self) -> float:
        if self.t != "real":
            raise MataError(3253, "nonreal found where real required" if self.t != "string"
                            else "string found where real required")
        return float(self.scalar())

    def int_scalar(self) -> int:
        x = self.real_scalar()
        if x >= SYS:
            raise MataError(3300, "argument out of range")
        return int(x)

    def str_scalar(self) -> str:
        if self.t != "string":
            raise MataError(3254, f"{self.t} found where string required")
        return str(self.scalar())


def real(x) -> MV:
    return MV(np.asarray(x, dtype=np.float64), "real")


def string(x) -> MV:
    arr = np.empty((1, 1), dtype=object)
    if isinstance(x, str):
        arr[0, 0] = x
        return MV(arr, "string")
    a = np.asarray(x, dtype=object)
    return MV(a, "string")


def empty(t: str = "real", r: int = 0, c: int = 0) -> MV:
    if t in ("string",):
        a = np.empty((r, c), dtype=object)
        a[...] = ""
        return MV(a, t)
    if t in ("pointer", "struct"):
        return MV(np.empty((r, c), dtype=object), t)
    if t == "complex":
        return MV(np.zeros((r, c), dtype=np.complex128), t)
    return MV(np.full((r, c), SYS), "real")


def bool_mv(b: bool) -> MV:
    return real(1.0 if b else 0.0)


def is_missing(a: np.ndarray) -> np.ndarray:
    return a >= SYS


def clean(a: np.ndarray) -> np.ndarray:
    """nan/inf viram missing."""
    a = np.asarray(a, dtype=np.float64)
    bad = ~np.isfinite(a)
    if bad.any():
        a = np.where(bad, SYS, a)
    return a


def truth(v: MV) -> bool:
    """Condição de if/while: escalar real; diferente de zero é verdadeiro
    (missing também é verdadeiro, como no Stata)."""
    if v.t != "real":
        raise MataError(3253, "nonreal found where real required")
    return v.scalar() != 0
