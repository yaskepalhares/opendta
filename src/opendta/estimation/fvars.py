"""Varlists com variáveis fatoriais e operadores de séries temporais
([U] 11.4.3 factor variables, [U] 11.4.4 time-series varlists).

    i.g          indicadores dos níveis de g (base: o menor nível)
    ib3.g        base 3;  ib(first|last|freq).g;  ibn.g sem base;  ib(#2).g
    i(1 3).g     só esses níveis;  i(1/4).g  intervalo
    c.x          contínua (necessário dentro de interações)
    a#b          interação;  a##b  fatorial completo (a b a#b)
    i.(a b)      o operador vale para cada variável do grupo
    L.x L2.x L(1/3).x F.x D.x D2.x S.x S4.x LD.x  (exigem tsset)

`expand_fv` devolve a lista de termos; `build_design` cria as colunas da
amostra, com os nomes que o Stata usa em e(b) ("2.g", "1b.g",
"2.g#c.x", "L.x") e as marcas de base (b) e de omissão (o).
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.varlist import expand as expand_names
from . import tsops

if TYPE_CHECKING:
    from ..core.dataset import Dataset

SYS = M.SYSMISS


@dataclass
class Component:
    """Um fator de um termo: variável com operadores."""
    var: str
    factor: bool                       # i. (ou nome puro dentro de #)
    ts: list = field(default_factory=list)   # [(letra, ordem)]
    base: str = ""                     # "", "n" (ibn), "first", "last", "freq", "#k", ou número
    levels: list | None = None         # i(1 3).g
    explicit_c: bool = False           # c. escrito pelo usuário

    @property
    def tsname(self) -> str:
        p = tsops.ops_name(self.ts)
        return f"{p}.{self.var}" if p else self.var

    def key(self) -> tuple:
        return (self.var, self.factor, tuple(self.ts))


@dataclass
class Term:
    comps: list[Component]

    def key(self) -> tuple:
        return tuple(sorted(c.key() for c in self.comps))

    @property
    def is_interaction(self) -> bool:
        return len(self.comps) > 1

    @property
    def has_factor(self) -> bool:
        return any(c.factor for c in self.comps)


@dataclass
class Column:
    name: str                          # nome em e(b)
    values: np.ndarray | None          # valores na amostra (None para base/omitida)
    term: int                          # índice do termo
    base: bool = False
    omitted: bool = False
    cell: tuple = ()                   # níveis dos fatores (para a tabela)
    labels: tuple = ()                 # (nome do termo para o cabeçalho, rótulo do nível)
    group: tuple = ()                  # colunas consecutivas com o mesmo grupo dividem cabeçalho
    rowlab: str = ""                   # rótulo da linha dentro do grupo ("2", "L1.")


# ---------------------------------------------------------------------------
# análise do texto
# ---------------------------------------------------------------------------

_OP = re.compile(r"""^(
      (?P<fv>i|ib(?:n|\d+|\((?:first|last|freq|\#\d+|\d+)\))|i\([^)]*\)|c|o|bn)
    | (?P<ts>(?:[LFDS](?:\d+|\(\d+/\d+\)|\(\d+(?:\s+\d+)*\))?)+)
    | (?P<mix>(?:[LFDS]\d*)*[ic](?:[LFDS]\d*)*)
)\.""", re.X | re.I)


def _split_tokens(text: str) -> list[str]:
    """Separa por espaços, respeitando parênteses: i.(a b) L(1/2).x i(1 2).g."""
    out, cur, depth = [], "", 0
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch.isspace() and depth == 0:
            if cur:
                out.append(cur)
            cur = ""
            continue
        cur += ch
    if cur:
        out.append(cur)
    return out


def _split_hash(tok: str) -> list[str]:
    """'a##b#c' -> ['a', 'b', 'c'], sem cortar dentro de parênteses (ib(#2).g)."""
    parts, cur, depth, i = [], "", 0, 0
    while i < len(tok):
        ch = tok[i]
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "#" and depth == 0:
            parts.append(cur)
            cur = ""
            i += 2 if tok.startswith("##", i) else 1
            continue
        cur += ch
        i += 1
    parts.append(cur)
    return parts


def _numlist(text: str) -> list[float]:
    vals: list[float] = []
    for part in text.replace(",", " ").split():
        if "/" in part:
            a, b = part.split("/")
            vals.extend(float(v) for v in range(int(float(a)), int(float(b)) + 1))
        else:
            vals.append(float(part))
    return vals


def _parse_ts(spec: str) -> list[list[tuple[str, int]]]:
    """'L2D' -> [[('L',2),('D',1)]]; 'L(1/3)' -> três alternativas."""
    alts: list[list[tuple[str, int]]] = [[]]
    for m in re.finditer(r"([LFDS])(\d+|\([^)]*\))?", spec, re.I):
        op = m.group(1).upper()
        arg = m.group(2)
        if arg and arg.startswith("("):
            ks = [int(v) for v in _numlist(arg[1:-1])]
        else:
            ks = [int(arg) if arg else 1]
        new = []
        for a in alts:
            for k in ks:
                new.append(a + ([(op, k)] if k != 0 else []))
        alts = new
    return alts


def _component_alts(ds: "Dataset", text: str, in_interaction: bool) -> list[Component] | None:
    """Um fator escrito (ex. 'ib2.g', 'L.x', 'c.x', 'i.(a b)') -> componentes."""
    factor = None
    explicit_c = False
    base = ""
    levels = None
    ts_alts: list[list] = [[]]
    rest = text
    while True:
        m = re.match(r"^([A-Za-z][A-Za-z0-9]*(?:\([^)]*\))?)\.", rest)
        if not m:
            break
        op = m.group(1)
        low = op.lower()
        if low == "i":
            factor = True
        elif low == "c":
            factor, explicit_c = False, True
        elif low.startswith("ib"):
            factor = True
            arg = op[2:]
            base = arg[1:-1] if arg.startswith("(") else arg
        elif low.startswith("i("):
            factor = True
            levels = _numlist(op[2:-1])
        elif low == "bn":
            factor, base = True, "n"
        elif low == "o":
            pass   # o.x: omitida pelo usuário — tratada como presente, VERIFICAR
        elif re.fullmatch(r"(?:[LFDS](?:\d+|\([^)]*\))?)+", op, re.I):
            ts_alts = [a + b for a in ts_alts for b in _parse_ts(op)]
        elif re.fullmatch(r"[A-Za-z\d]+", op) and re.search(r"[ic]", op, re.I):
            # combinações como iL.x ou cL.x
            letters = re.sub(r"[icIC]", "", op)
            factor = "i" in op.lower()
            explicit_c = "c" in op.lower()
            if letters:
                ts_alts = [a + b for a in ts_alts for b in _parse_ts(letters)]
        else:
            return None
        rest = rest[m.end():]
    if rest.startswith("(") and rest.endswith(")"):
        names = expand_names(ds, rest[1:-1], allow_empty=False)
    else:
        if not re.fullmatch(r"[A-Za-z_*?~][\w*?~-]*", rest):
            return None
        names = expand_names(ds, rest, allow_empty=False)
    if factor is None:
        factor = in_interaction
    out = []
    for nm in names:
        for ts in ts_alts:
            out.append(Component(nm, factor, list(ts), base, levels, explicit_c))
    return out


def expand_fv(ds: "Dataset", text: str) -> list[Term]:
    """Varlist com fv/ts -> termos, na ordem em que aparecem."""
    terms: list[Term] = []
    seen: set = set()

    def add(t: Term) -> None:
        k = t.key()
        if k not in seen:
            seen.add(k)
            terms.append(t)

    for tok in _split_tokens(text):
        full = "##" in re.sub(r"\([^)]*\)", "", tok)
        parts = _split_hash(tok)
        inter = len(parts) > 1
        alts_per_part = []
        for part in parts:
            comps = _component_alts(ds, part, inter)
            if comps is None:
                raise StataError(198, f"{tok} invalid name")   # VERIFICAR
            alts_per_part.append(comps)
        for combo in itertools.product(*alts_per_part):
            combo = list(combo)
            if full:
                for r_ in range(1, len(combo) + 1):
                    for sub in itertools.combinations(combo, r_):
                        comps = [Component(c.var, c.factor if r_ > 1 or c.factor else False, c.ts, c.base,
                                           c.levels, c.explicit_c) for c in sub]
                        if r_ == 1 and not sub[0].factor:
                            comps[0].explicit_c = True
                        add(Term(comps))
            else:
                add(Term(combo))
    return terms


def base_vars(terms: list[Term]) -> list[str]:
    out: list[str] = []
    for t in terms:
        for c in t.comps:
            if c.var not in out:
                out.append(c.var)
    return out


def has_fv_or_ts(text: str) -> bool:
    return bool(re.search(r"(^|[\s(#])[A-Za-z0-9()/# ]*\.", text)) or "#" in text


# ---------------------------------------------------------------------------
# colunas
# ---------------------------------------------------------------------------

def comp_values(ds: "Dataset", c: Component) -> np.ndarray:
    v = ds.get(c.var)
    if v.is_string:
        raise StataError(109, "string variables not allowed in varlist;\n"
                         f"{c.var} is a string variable")
    x = v.data.astype(np.float64)
    if c.ts:
        x = tsops.apply_ops(ds, c.ts, x)
    return x


def term_missing(ds: "Dataset", terms: list[Term]) -> np.ndarray:
    """Observações com algum valor missing nos termos (para markout)."""
    bad = np.zeros(ds.nobs, dtype=bool)
    for t in terms:
        for c in t.comps:
            bad |= comp_values(ds, c) >= SYS
    return bad


def _levels(ds: "Dataset", c: Component, x: np.ndarray) -> tuple[list[float], float | None]:
    vals = x[x < SYS]
    if np.any(vals < 0):
        raise StataError(452, f"{c.var}:  factor variables may not contain negative values")
    if np.any(vals != np.trunc(vals)):
        raise StataError(452, f"{c.var}:  factor variables may not contain noninteger values")
    levels = sorted(set(vals.tolist()))
    if c.levels is not None:
        levels = [lv for lv in levels if lv in set(c.levels)]
    if not levels:
        return [], None
    b = c.base
    if b == "n":
        return levels, None
    if b == "" or b == "first":
        base = levels[0]
    elif b == "last":
        base = levels[-1]
    elif b == "freq":
        counts = {lv: int(np.sum(vals == lv)) for lv in levels}
        base = max(levels, key=lambda lv: (counts[lv], -lv))
    elif b.startswith("#"):
        k = int(b[1:])
        base = levels[min(max(k, 1), len(levels)) - 1]
    else:
        base = float(b)
        if base not in levels:
            # VERIFICAR: o Stata aceita base ausente na amostra?
            levels = sorted(levels + [base])
    return levels, base


def _level_text(lv: float) -> str:
    return f"{int(lv)}" if lv == int(lv) else f"{lv:g}"


def build_design(ds: "Dataset", terms: list[Term], mask: np.ndarray, *,
                 constant: bool = True) -> list[Column]:
    """Colunas na amostra `mask`, na ordem do Stata, com bases marcadas."""
    present = {t.key() for t in terms}
    cols: list[Column] = []
    for ti, t in enumerate(terms):
        vals = [comp_values(ds, c)[mask] for c in t.comps]
        info = []
        for c, x in zip(t.comps, vals):
            if c.factor:
                info.append(_levels(ds, c, x))
            else:
                info.append(None)
        factor_idx = [k for k, c in enumerate(t.comps) if c.factor]
        level_lists = [info[k][0] for k in factor_idx]
        for cell in itertools.product(*level_lists) if factor_idx else [()]:
            # base: algum subconjunto de fatores no nível-base cujo termo
            # complementar está no modelo (vazio = constante)
            at_base = [factor_idx[j] for j, lv in enumerate(cell) if info[factor_idx[j]][1] == lv]
            is_base = False
            for r_ in range(1, len(at_base) + 1):
                for sub in itertools.combinations(at_base, r_):
                    rest = [t.comps[k] for k in range(len(t.comps)) if k not in sub]
                    if not rest:
                        if constant or len(t.comps) == 1:
                            is_base = True
                    elif Term(rest).key() in present:
                        is_base = True
                    if is_base:
                        break
                if is_base:
                    break
            parts = []
            y = np.ones(int(mask.sum()))
            levmap = dict(zip(factor_idx, cell))
            for k, c in enumerate(t.comps):
                if c.factor:
                    lv = levmap[k]
                    mark = ""
                    if info[k][1] == lv:
                        mark = "b"
                    elif is_base:
                        mark = "o"
                    pre = tsops.ops_name(c.ts)
                    parts.append(f"{_level_text(lv)}{mark}{pre}.{c.var}")
                    y = y * (vals[k] == lv)
                else:
                    if t.is_interaction or c.explicit_c:
                        pre = tsops.ops_name(c.ts)
                        parts.append(f"c{pre}.{c.var}" if pre else f"c.{c.var}")
                    else:
                        parts.append(c.tsname)
                    y = y * vals[k]
            name = "#".join(parts)
            if not t.is_interaction and not t.comps[0].factor:
                name = t.comps[0].tsname   # x, L.x: sem o "c."
            if t.is_interaction:
                head = "#".join(c.tsname if c.factor else f"c.{c.tsname}" for c in t.comps)
            else:
                head = t.comps[0].tsname if not t.comps[0].factor else t.comps[0].var
            lab = " ".join(_level_text(lv) for lv in cell)
            if factor_idx:
                group, rowlab = ("t", ti), lab
            elif not t.is_interaction:
                c0 = t.comps[0]
                group, rowlab = ("v", c0.var), ts_rowlabel(c0.ts)
                head = c0.tsname
            else:
                group, rowlab = (), ""
            cols.append(Column(name, None if is_base else y.astype(np.float64), ti,
                               base=is_base, cell=tuple(cell), labels=(head, lab),
                               group=group, rowlab=rowlab))
    return cols


def ts_rowlabel(ops: list) -> str:
    """Rótulo da linha de um operador na tabela: L1. L2. D1. LD. --. (VERIFICAR)."""
    if not ops:
        return "--."
    if len(ops) == 1:
        op, k = ops[0]
        return f"{op}{k}."
    return "".join(op + (str(k) if k != 1 else "") for op, k in ops) + "."


def drop_collinear(cols: list[Column], constant: bool, w: np.ndarray | None = None,
                   tol: float = 1e-13) -> list[str]:
    """Marca como omitidas as colunas colineares; devolve os nomes para as
    notas "omitted because of collinearity".

    Varre a matriz de produtos cruzados (centrada, se há constante) como o
    invsym() do Stata, que escolhe as colunas "para minimizar o erro de
    arredondamento": a cada passo, o maior elemento restante da diagonal.
    Com x3 = 2*x1, o Stata omite x1 (compat 0501). VERIFICAR a regra geral
    e a tolerância."""
    notes: list[str] = []
    live = [c for c in cols if not c.base and not c.omitted and c.values is not None]
    if not live:
        return notes
    n = len(live[0].values)
    ww = np.ones(n) if w is None else np.asarray(w, dtype=np.float64)
    X = np.column_stack([c.values for c in live])
    if constant:
        X = X - (ww @ X) / ww.sum()
    A = (X * ww[:, None]).T @ X
    k = A.shape[0]
    orig = np.diag(A).copy()
    done = np.zeros(k, dtype=bool)
    dropped = np.zeros(k, dtype=bool)
    for _ in range(k):
        cand = [j for j in range(k) if not done[j]]
        if not cand:
            break
        j = max(cand, key=lambda q: A[q, q])
        d = A[j, j]
        done[j] = True
        if d <= tol * max(orig[j], 1e-300) or orig[j] <= 0:
            dropped[j] = True
            continue
        # varredura (sweep) no pivô j
        row = A[j, :].copy()
        A = A - np.outer(row, row) / d
        A[j, :] = 0
        A[:, j] = 0
    # colunas não varridas que zeraram também são colineares
    for j in range(k):
        if not done[j] and A[j, j] <= tol * max(orig[j], 1e-300):
            dropped[j] = True
    for j, c in enumerate(live):
        if dropped[j] or orig[j] <= 0:
            c.omitted = True
            c.values = None
            notes.append(c.name)
    return notes


def omitted_name(name: str) -> str:
    """'x' -> 'o.x'; '2.g' -> '2o.g'; '2.g#c.x' -> '2o.g#co.x'."""
    parts = []
    for p in name.split("#"):
        m = re.match(r"^(\d+(?:\.\d+)?)(b?)(\w*)\.(.+)$", p)
        if m and not p.startswith("c."):
            parts.append(f"{m.group(1)}{m.group(2) or ''}o{m.group(3)}.{m.group(4)}"
                         if "o" not in m.group(2) else p)
        elif p.startswith("c."):
            parts.append("co." + p[2:])
        else:
            parts.append("o." + p)
    return "#".join(parts)
