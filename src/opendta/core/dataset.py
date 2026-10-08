"""Conjunto de dados em memória.

Cada variável numérica fica na memória no próprio tipo (byte = 1 byte,
int = 2, long e float = 4, double = 8), como no Stata; ver core/storage.py.
Os cálculos usam vetores float64 com os códigos de missing do OpenDTA (ver
core/missing.py), obtidos por `Variable.data`. O tipo define a faixa de
valores aceitos e a precisão: uma variável float guarda valores em
precisão simples, de modo que `gen x = 0.1` seguido de
`count if x == 0.1` dá 0.

Strings são vetores de objetos Python (str). O tipo str# acompanha o maior
comprimento em bytes (UTF-8); acima de 2045 bytes a variável vira strL.

Faixas (manual [D] data types):
    byte   -127 .. 100
    int    -32,767 .. 32,740
    long   -2,147,483,647 .. 2,147,483,620
    float  ±1.70141173319e+38 (precisão de ~7 dígitos)
    double ±8.9884656743e+307
"""

from __future__ import annotations

import re

import numpy as np

from . import missing as M
from . import storage
from .errors import StataError

NUMERIC_TYPES = ("byte", "int", "long", "float", "double")
INT_RANGES = {
    "byte": (-127, 100),
    "int": (-32767, 32740),
    "long": (-2147483647, 2147483620),
}
FLOAT_MAX = 1.7014117331926443e38   # maior float não missing (bits 0x7effffff)
STR_MAX = 2045
TYPE_ORDER = {"byte": 0, "int": 1, "long": 2, "float": 3, "double": 4}
_STR_TYPE = re.compile(r"^str(\d+|L)$")

DEFAULT_NUMERIC_FORMAT = {
    "byte": "%8.0g", "int": "%8.0g", "long": "%12.0g", "float": "%9.0g", "double": "%10.0g",
}

NAME_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]{0,31}$")
RESERVED = {
    "_all", "_b", "byte", "_coef", "_cons", "double", "float", "if", "in", "int",
    "long", "_n", "_N", "_pi", "_pred", "_rc", "_se", "_skip", "using", "with",
    "strL",
}


def is_string_type(t: str) -> bool:
    return _STR_TYPE.match(t) is not None


def str_type_for(length: int) -> str:
    return "strL" if length > STR_MAX else f"str{max(1, length)}"


def str_len(s: str) -> int:
    return len(s.encode("utf-8"))


def check_name(name: str) -> None:
    if not NAME_RE.match(name) or name in RESERVED:
        raise StataError(198, f"{name} invalid name")


def default_format(vtype: str) -> str:
    if vtype in DEFAULT_NUMERIC_FORMAT:
        return DEFAULT_NUMERIC_FORMAT[vtype]
    if vtype == "strL":
        return "%9s"
    n = int(vtype[3:])
    # VERIFICAR: formato padrão de strings criadas por generate
    return f"%{max(9, n)}s"


def fit_numeric(values: np.ndarray, vtype: str) -> tuple[np.ndarray, int]:
    """Ajusta valores ao tipo: precisão de float, truncamento em inteiros e
    faixa válida. Devolve (vetor ajustado, nº de valores que viraram missing
    por estarem fora da faixa)."""
    v = np.array(values, dtype=np.float64, copy=True)
    miss = v >= M.SYSMISS
    nonmiss = ~miss
    lost = 0
    if vtype == "double":
        return v, 0
    if vtype == "float":
        over = nonmiss & (np.abs(v) > FLOAT_MAX)
        lost = int(over.sum())
        v[over] = M.SYSMISS
        ok = nonmiss & ~over
        v[ok] = v[ok].astype(np.float32).astype(np.float64)
        return v, lost
    lo, hi = INT_RANGES[vtype]
    # VERIFICAR: não inteiros em tipos inteiros são truncados em direção a zero
    v[nonmiss] = np.trunc(v[nonmiss])
    over = nonmiss & ((v < lo) | (v > hi))
    lost = int(over.sum())
    v[over] = M.SYSMISS
    return v, lost


def smallest_type_for(values: np.ndarray) -> str:
    """Menor tipo que guarda os valores sem perda (usado por compress e promoção)."""
    v = values[values < M.SYSMISS]
    if v.size == 0:
        return "byte"
    if np.all(v == np.trunc(v)):
        lo, hi = v.min(), v.max()
        for t in ("byte", "int", "long"):
            a, b = INT_RANGES[t]
            if lo >= a and hi <= b:
                return t
    as_float = v.astype(np.float32).astype(np.float64)
    if np.all(as_float == v) and np.all(np.abs(v) <= FLOAT_MAX):
        return "float"
    return "double"


class Variable:
    """Uma variável. Numéricas ficam no tipo nativo (ver core/storage.py) em
    `raw`; `data` devolve o vetor float64 para cálculos e aceita um float64
    na atribuição. Strings ficam num vetor de objetos str (`raw` e `data`
    são o mesmo vetor).

    `data` é só para leitura: para mudar um valor use `set_value`, ou
    atribua um vetor inteiro a `data`. Para um valor só, `value(i)` evita
    decodificar a coluna inteira."""

    def __init__(self, name: str, vtype: str, data, fmt: str = "", label: str = "",
                 value_label: str = "", notes: list[str] | None = None) -> None:
        self.name = name
        self._vtype = vtype
        self.fmt = fmt or default_format(vtype)
        self.label = label
        self.value_label = value_label
        self.notes = notes if notes is not None else []
        self.raw: np.ndarray = np.empty(0)
        self.data = data

    def __repr__(self) -> str:
        return f"Variable({self.name!r}, {self._vtype!r}, n={len(self.raw)})"

    # -- tipo ----------------------------------------------------------------
    @property
    def vtype(self) -> str:
        return self._vtype

    @vtype.setter
    def vtype(self, new: str) -> None:
        if new == self._vtype:
            return
        if is_string_type(new) != is_string_type(self._vtype):
            raise ValueError("troca entre número e string exige uma variável nova")
        if not is_string_type(new):
            values = self.data
            self._vtype = new
            self.raw = storage.encode(values, new)
        else:
            self._vtype = new

    @property
    def is_string(self) -> bool:
        return is_string_type(self._vtype)

    # -- valores ---------------------------------------------------------------
    @property
    def data(self) -> np.ndarray:
        if self.is_string or self._vtype == "double":
            view = self.raw.view()
            view.flags.writeable = False
            return view
        return storage.decode(self.raw, self._vtype)

    @data.setter
    def data(self, values) -> None:
        if self.is_string:
            self.raw = np.array(values, dtype=object, copy=True)
        else:
            self.raw = storage.encode(values, self._vtype)

    def value(self, i: int):
        x = self.raw[i]
        return x if self.is_string else storage.decode_one(x, self._vtype)

    def set_value(self, i: int, x) -> None:
        self.raw[i] = x if self.is_string else storage.encode_one(float(x), self._vtype)

    def __len__(self) -> int:
        return len(self.raw)

    @property
    def nbytes(self) -> int:
        """Memória ocupada pelos valores (strings: estimativa pelo conteúdo)."""
        if self.is_string:
            return int(sum(len(s) for s in self.raw)) + 8 * len(self.raw)
        return int(self.raw.nbytes)

    def max_strlen(self) -> int:
        if not self.is_string or len(self.raw) == 0:
            return 1
        return max((str_len(s) for s in self.raw), default=1) or 1


class Dataset:
    def __init__(self) -> None:
        self.vars: list[Variable] = []
        self.nobs = 0
        self.label = ""
        self.notes: list[str] = []
        self.value_labels: dict[str, dict[int, str]] = {}
        self.sortlist: list[str] = []
        self.filename = ""             # como foi digitado (describe mostra assim)
        self.fullpath = ""             # caminho absoluto, para save sem nome
        self.timestamp = ""            # data de gravação ("13 Apr 2014 17:45")
        self.chars: dict[str, dict[str, str]] = {}   # características (_dta e variáveis)
        self.changed = False

    def copy(self) -> "Dataset":
        """Cópia independente (preserve)."""
        import copy as _copy
        new = Dataset()
        for v in self.vars:
            nv = Variable(v.name, v.vtype, np.empty(0), fmt=v.fmt, label=v.label,
                          value_label=v.value_label, notes=list(v.notes))
            nv.raw = v.raw.copy()
            new.vars.append(nv)
        new.nobs = self.nobs
        new.label = self.label
        new.notes = list(self.notes)
        new.value_labels = _copy.deepcopy(self.value_labels)
        new.sortlist = list(self.sortlist)
        new.filename, new.fullpath, new.timestamp = self.filename, self.fullpath, self.timestamp
        new.chars = _copy.deepcopy(self.chars)
        new.changed = self.changed
        return new

    # -- consulta --------------------------------------------------------
    @property
    def nvars(self) -> int:
        return len(self.vars)

    @property
    def names(self) -> list[str]:
        return [v.name for v in self.vars]

    def has(self, name: str) -> bool:
        return any(v.name == name for v in self.vars)

    def get(self, name: str) -> Variable:
        for v in self.vars:
            if v.name == name:
                return v
        raise StataError(111, f"variable {name} not found")

    def index(self, name: str) -> int:
        for i, v in enumerate(self.vars):
            if v.name == name:
                return i
        raise StataError(111, f"variable {name} not found")

    @property
    def empty(self) -> bool:
        return self.nvars == 0 and self.nobs == 0

    def width(self) -> int:
        """Bytes por observação, como em describe (VERIFICAR para strL)."""
        size = {"byte": 1, "int": 2, "long": 4, "float": 4, "double": 8}
        total = 0
        for v in self.vars:
            if v.vtype in size:
                total += size[v.vtype]
            elif v.vtype == "strL":
                total += 8
            else:
                total += int(v.vtype[3:])
        return total

    # -- alterações --------------------------------------------------------
    def _touch(self) -> None:
        self.changed = True

    def clear(self) -> None:
        self.__init__()

    def set_obs(self, n: int) -> None:
        if n < self.nobs:
            raise StataError(198, f"observations must be at least {self.nobs}")
        extra = n - self.nobs
        for v in self.vars:
            if v.is_string:
                pad = np.array([""] * extra, dtype=object)
            else:
                pad = np.full(extra, storage.missing_raw(v.vtype), dtype=v.raw.dtype)
            v.raw = np.concatenate([v.raw, pad])
        self.nobs = n
        self._touch()

    def add(self, var: Variable, *, position: int | None = None) -> None:
        if self.has(var.name):
            raise StataError(110, f"variable {var.name} already defined")
        check_name(var.name)
        if len(var) != self.nobs:
            raise ValueError("tamanho do vetor diferente do número de observações")
        if position is None:
            self.vars.append(var)
        else:
            self.vars.insert(position, var)
        self._touch()   # criar variável não desfaz a ordenação

    def drop_vars(self, names: list[str]) -> None:
        drop = set(names)
        self.vars = [v for v in self.vars if v.name not in drop]
        for n in drop:
            self.chars.pop(n, None)
        self.sortlist = [s for s in self.sortlist if s not in drop]
        if not self.vars:
            self.nobs = 0
        self._touch()

    def keep_obs(self, mask: np.ndarray) -> int:
        """Mantém as observações com mask verdadeiro; devolve quantas saíram."""
        mask = np.asarray(mask, dtype=bool)
        removed = int(self.nobs - mask.sum())
        if removed:
            for v in self.vars:
                v.raw = v.raw[mask]
            self.nobs = int(mask.sum())
            self._touch()
        return removed

    def reorder_obs(self, order: np.ndarray) -> None:
        for v in self.vars:
            v.raw = v.raw[order]
        self._touch()

    def rename(self, old: str, new: str) -> None:
        if old == new:
            return
        if self.has(new):
            raise StataError(110, f"variable {new} already defined")
        check_name(new)
        self.get(old).name = new
        self.sortlist = [new if s == old else s for s in self.sortlist]
        if old in self.chars:
            self.chars[new] = self.chars.pop(old)
        self._touch()

    def order(self, names: list[str], *, last: bool = False,
              before: str | None = None, after: str | None = None) -> None:
        moving = [self.get(n) for n in names]
        rest = [v for v in self.vars if v.name not in set(names)]
        if before or after:
            anchor = before or after
            idx = [v.name for v in rest].index(anchor) if anchor in [v.name for v in rest] else None
            if idx is None:
                raise StataError(198, f"{anchor} may not be in varlist")
            pos = idx if before else idx + 1
            self.vars = rest[:pos] + moving + rest[pos:]
        elif last:
            self.vars = rest + moving
        else:
            self.vars = moving + rest
        self._touch()

    # -- tipos -------------------------------------------------------------
    def set_numeric(self, var: Variable, values: np.ndarray, rows: np.ndarray | None = None,
                    *, promote: bool = True) -> tuple[int, str | None]:
        """Grava valores numéricos (todas as linhas ou só `rows`), promovendo o
        tipo se necessário, como faz o replace. Devolve (nº de mudanças,
        mensagem de promoção ou None)."""
        if rows is None:
            new_full = np.asarray(values, dtype=np.float64).copy()
        else:
            new_full = np.array(var.data, dtype=np.float64)
            new_full[rows] = values
        note = None
        if promote and var.vtype not in ("double",):
            target = var.vtype
            subset = new_full if rows is None else new_full[rows]
            needed = smallest_type_for(subset)
            if var.vtype == "float" and needed == "double":
                needed = "float"  # replace não promove float para double
            if TYPE_ORDER[needed] > TYPE_ORDER[var.vtype]:
                target = needed
            if target != var.vtype:
                # o formato de exibição não muda (b byte %8.0g vira float %8.0g)
                note = f"variable {var.name} was {var.vtype} now {target}"
                var.vtype = target
        fitted, _ = fit_numeric(new_full, var.vtype)
        changed = int(np.sum(fitted != var.data))
        var.data = fitted
        if changed:
            self._unsort_from(var.name)
        self._touch()
        return changed, note

    def set_string(self, var: Variable, values: np.ndarray, rows: np.ndarray | None = None
                   ) -> tuple[int, str | None]:
        new_full = np.array(var.raw, dtype=object)
        if rows is None:
            new_full = np.array(values, dtype=object)
        else:
            new_full[rows] = values
        note = None
        longest = max((str_len(s) for s in new_full), default=1) or 1
        if var.vtype != "strL":
            cur = int(var.vtype[3:])
            if longest > cur:
                target = str_type_for(longest)
                note = f"variable {var.name} was {var.vtype} now {target}"
                if var.fmt == default_format(var.vtype):
                    var.fmt = default_format(target)
                var.vtype = target
        changed = int(sum(1 for a, b in zip(new_full, var.data) if a != b))
        var.data = new_full
        if changed:
            self._unsort_from(var.name)
        self._touch()
        return changed, note

    def _unsort_from(self, name: str) -> None:
        """Alterar uma variável de ordenação desfaz a ordenação a partir dela."""
        if name in self.sortlist:
            self.sortlist = self.sortlist[:self.sortlist.index(name)]
