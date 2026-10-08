"""append, merge, joinby e cross ([D] append, merge, joinby, cross).

Regras comuns ao juntar dois conjuntos de dados:

* atributos (formato, rótulos) da base na memória (master) prevalecem;
* variáveis numéricas com tipos diferentes são promovidas para um tipo que
  guarde as duas (long + float vira double); strings, para o maior str#;
* string de um lado e número do outro é erro r(106), a menos que se use
  `force` (o lado da base usada fica missing);
* rótulos de valores da base usada entram se não existirem na master; os
  que já existem são mantidos, com a nota "(label x already defined)".

merge imprime a tabela de resultados do Stata 11+ (colunas VERIFICAR).
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import sorting
from ..core.dataset import (TYPE_ORDER, Dataset, Variable, check_name, default_format,
                            is_string_type, str_len, str_type_for)
from ..core.errors import StataError
from ..core.grouping import joint_keys
from ..core.storage import missing_raw
from ..io.dta import read_dta
from ..lang.syntax import match_options, parse_standard
from ..lang.words import split_words, strip_outer_quotes
from .files import dta_path
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


# ---------------------------------------------------------------------------
# auxiliares
# ---------------------------------------------------------------------------

def load_using(text: str) -> Dataset:
    path = dta_path(text)
    if not path.exists():
        raise StataError(601, f"file {path} not found")
    return read_dta(path)


def _is_missing(v: Variable) -> np.ndarray:
    if v.is_string:
        return np.array([t == "" for t in v.raw], dtype=bool)
    return v.data >= M.SYSMISS


def common_type(a: str, b: str) -> str:
    """Tipo que guarda valores dos dois tipos (mesmo gênero: número/string)."""
    if is_string_type(a):
        if a == "strL" or b == "strL":
            return "strL"
        return str_type_for(max(int(a[3:]), int(b[3:])))
    if {a, b} == {"long", "float"}:
        return "double"
    return a if TYPE_ORDER[a] >= TYPE_ORDER[b] else b


def _mismatch(name: str, mv: Variable, uv: Variable, force: bool, what: str) -> bool:
    """Erro (ou True com force) quando um lado é string e o outro número."""
    if mv.is_string == uv.is_string:
        return False
    if force:
        return True
    # VERIFICAR texto exato
    filler = '""' if mv.is_string else "numeric missing value"
    raise StataError(106, f"variable {name} is {mv.vtype} in master but {uv.vtype} in using data\n"
                          f"    You could specify {what}'s force option to ignore this string/numeric "
                          f"mismatch.  The using\n    variable would then be treated as if it contained "
                          f"{filler}.")


def _promote(s: "Session", mv: Variable, utype: str, quiet: bool = False) -> None:
    target = common_type(mv.vtype, utype)
    if target != mv.vtype:
        old = mv.vtype
        if is_string_type(target) and mv.fmt == default_format(old):
            mv.fmt = default_format(target)
        mv.vtype = target
        if not quiet:
            # VERIFICAR: nota impressa pelo append/merge do Stata 14
            s.output.write(f"(note: variable {mv.name} was {old}, now {target} to accommodate "
                           f"using data's values)\n", "text")


def _blank(vtype: str, n: int) -> np.ndarray:
    if is_string_type(vtype):
        return np.array([""] * n, dtype=object)
    return np.full(n, missing_raw(vtype))


def _raw_as(uv: Variable, vtype: str) -> np.ndarray:
    """Valores de uv no tipo de armazenamento vtype."""
    if is_string_type(vtype):
        return np.array(uv.raw, dtype=object) if uv.is_string else _blank(vtype, len(uv))
    if uv.is_string:
        return _blank(vtype, len(uv))
    tmp = Variable("_tmp", vtype, uv.data)
    return tmp.raw


def _merge_value_labels(s: "Session", master: Dataset, using: Dataset, names_used: set[str]) -> None:
    for lab in sorted(names_used):
        if lab not in using.value_labels:
            continue
        if lab in master.value_labels:
            s.output.write(f"(label {lab} already defined)\n", "text")   # VERIFICAR
        else:
            master.value_labels[lab] = dict(using.value_labels[lab])


def _fix_strings(v: Variable) -> None:
    """Ajusta str# ao maior valor presente (depois de copiar valores)."""
    if v.is_string and v.vtype != "strL":
        longest = max((str_len(t) for t in v.raw), default=1) or 1
        if longest > int(v.vtype[3:]):
            v.vtype = str_type_for(longest)


# ---------------------------------------------------------------------------
# append
# ---------------------------------------------------------------------------

@command("append", "app")
def cmd_append(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.using is None or p.varlist.strip():
        raise StataError(100, "using required") if p.using is None else StataError(198, "invalid syntax")
    opts = match_options(p.options, {"generate": 3, "keep": 4, "nolabel": 5, "nonotes": 5, "force": 5})
    files = split_words(p.using)
    gen = str(opts["generate"]).strip() if opts.get("generate") not in (None, True) else ""
    ds = s.data
    if gen:
        if ds.has(gen):
            raise StataError(110, f"variable {gen} already defined")
        check_name(gen)
    usings = [load_using(f) for f in files]
    force = bool(opts.get("force"))
    origin = [np.zeros(ds.nobs)]
    for j, U in enumerate(usings, start=1):
        if opts.get("keep") not in (None, True):
            from ..core.varlist import expand
            keep = set(expand(U, str(opts["keep"])))
            U.drop_vars([n for n in U.names if n not in keep])
        n0, nu = ds.nobs, U.nobs
        # tipos e mismatches primeiro (nada muda se houver erro)
        bad = set()
        for uv in U.vars:
            if ds.has(uv.name) and _mismatch(uv.name, ds.get(uv.name), uv, force, "append"):
                bad.add(uv.name)
        for uv in U.vars:
            if ds.has(uv.name) and uv.name not in bad:
                _promote(s, ds.get(uv.name), uv.vtype)
        for mv in ds.vars:
            if ds.has(mv.name) and U.has(mv.name) and mv.name not in bad:
                add = _raw_as(U.get(mv.name), mv.vtype)
            else:
                add = _blank(mv.vtype, nu)
            mv.raw = np.concatenate([mv.raw, add])
            _fix_strings(mv)
        for uv in U.vars:
            if not any(v.name == uv.name for v in ds.vars):
                nv = Variable(uv.name, uv.vtype, np.empty(0), fmt=uv.fmt, label=uv.label,
                              value_label="" if opts.get("nolabel") else uv.value_label,
                              notes=[] if opts.get("nonotes") else list(uv.notes))
                nv.raw = np.concatenate([_blank(uv.vtype, n0), uv.raw])
                ds.vars.append(nv)
        ds.nobs = n0 + nu
        origin.append(np.full(nu, float(j)))
        if not opts.get("nolabel"):
            used = {v.value_label for v in U.vars if v.value_label}
            _merge_value_labels(s, ds, U, used)
        if not opts.get("nonotes") and U.notes:
            ds.notes.extend(n for n in U.notes if n not in ds.notes)
    ds.sortlist = []
    ds.changed = True
    if gen:
        ds.add(Variable(gen, "byte", np.concatenate(origin)))
    s.notify_state()


# ---------------------------------------------------------------------------
# merge
# ---------------------------------------------------------------------------

_RESULTS = {"master": [1], "masters": [1], "using": [2], "usings": [2], "match": [3], "matches": [3],
            "matched": [3], "match_update": [4], "match_updates": [4], "match_conflict": [5],
            "match_conflicts": [5]}
_MERGE_LABELS = {1: "master only (1)", 2: "using only (2)", 3: "matched (3)",
                 4: "missing updated (4)", 5: "nonmissing conflict (5)"}


def _results(text: str) -> set[int]:
    out: set[int] = set()
    for w in text.replace(",", " ").split():
        wl = w.lower()
        if wl in _RESULTS:
            out.update(_RESULTS[wl])
        elif wl in ("1", "2", "3", "4", "5"):
            out.add(int(wl))
        else:
            raise StataError(198, f"invalid results specification {w}")   # VERIFICAR
    # VERIFICAR: com update, "match" talvez inclua também 4 e 5
    return out


def _check_unique(ds: Dataset, keys: list[str], where: str) -> None:
    cols = [ds.get(k).raw.astype(str) if ds.get(k).is_string else ds.get(k).data for k in keys]
    if ds.nobs < 2:
        return
    if len(keys) == 1:
        _, counts = np.unique(cols[0], return_counts=True)
    else:
        from ..core.grouping import group_ids
        ids, _ = group_ids(ds, keys)
        counts = np.bincount(ids)
    if (counts > 1).any():
        word = "variable" if len(keys) == 1 else "variables"
        verb = "does" if len(keys) == 1 else "do"
        raise StataError(459, f"{word} {' '.join(keys)} {verb} not uniquely identify observations "
                              f"in the {where} data")


def _key_col(v: Variable) -> np.ndarray:
    return np.array(v.raw, dtype=object) if v.is_string else np.asarray(v.data, dtype=np.float64)


def _pairs(kind: str, gm: np.ndarray, gu: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Pares (linha master ou -1, linha using ou -1), ordenados pela chave."""
    order_m = np.argsort(gm, kind="stable")
    order_u = np.argsort(gu, kind="stable")
    bm = np.searchsorted(gm[order_m], np.arange(k + 1))
    bu = np.searchsorted(gu[order_u], np.arange(k + 1))
    pm: list[np.ndarray] = []
    pu: list[np.ndarray] = []
    for g in range(k):
        rm = order_m[bm[g]:bm[g + 1]]
        ru = order_u[bu[g]:bu[g + 1]]
        if len(rm) and len(ru):
            if kind == "m:m":
                L = max(len(rm), len(ru))
                pm.append(rm[np.minimum(np.arange(L), len(rm) - 1)])
                pu.append(ru[np.minimum(np.arange(L), len(ru) - 1)])
            elif len(rm) >= len(ru):          # 1:1 e m:1
                pm.append(rm)
                pu.append(np.repeat(ru[:1], len(rm)))
            else:                             # 1:m
                pm.append(np.repeat(rm[:1], len(ru)))
                pu.append(ru)
        elif len(rm):
            pm.append(rm)
            pu.append(np.full(len(rm), -1))
        else:
            pm.append(np.full(len(ru), -1))
            pu.append(ru)
    if not pm:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return np.concatenate(pm).astype(np.int64), np.concatenate(pu).astype(np.int64)


def _take(v: Variable, rows: np.ndarray, vtype: str | None = None) -> np.ndarray:
    """Valores de v nas linhas `rows` (-1 = missing), no tipo vtype."""
    vtype = vtype or v.vtype
    src = _raw_as(v, vtype) if vtype != v.vtype else v.raw
    out = _blank(vtype, len(rows))
    ok = rows >= 0
    out[ok] = src[rows[ok]]
    return out


def _report(s: "Session", codes: np.ndarray, genname: str, update: bool) -> None:
    def line(label: str, n: int, indent: int, code: int | None = None) -> None:
        left = " " * indent + label
        tail = f"  ({genname}=={code})" if code is not None else ""
        s.output.write(left, "text")
        s.output.write(f"{n:,}".rjust(45 - len(left)), "result")
        s.output.write(tail + "\n", "text")

    c = {k: int((codes == k).sum()) for k in range(1, 6)}
    w = s.output.write
    w("\n    Result                           # of obs.\n", "text")
    w("    " + "-" * 41 + "\n", "text")
    line("not matched", c[1] + c[2], 4)
    if c[1] + c[2]:
        line("from master", c[1], 8, 1)
        line("from using", c[2], 8, 2)
        w("\n", "text")
    if update:
        line("matched", c[3] + c[4] + c[5], 4)
        line("not updated", c[3], 8, 3)
        line("missing updated", c[4], 8, 4)
        line("nonmissing conflict", c[5], 8, 5)
    else:
        line("matched", c[3], 4, 3)
    w("    " + "-" * 41 + "\n", "text")


@command("merge", "mer")
def cmd_merge(s: "Session", args: str) -> None:
    m = re.match(r"^\s*(1:1|m:1|1:m|m:m)\s+(.*)$", args, re.S)
    if not m:
        # VERIFICAR: sintaxe antiga (merge varlist using) não implementada
        raise StataError(198, "invalid syntax")
    kind = m.group(1)
    p = parse_standard(m.group(2))
    if p.using is None:
        raise StataError(100, "using required")
    opts = match_options(p.options, {
        "keepusing": 5, "generate": 3, "nogenerate": 5, "nolabel": 5, "nonotes": 5,
        "update": 3, "replace": 3, "noreport": 5, "force": 5, "assert": 6, "keep": 4, "sorted": 6})
    ds = s.data
    by_n = p.varlist.strip() == "_n"
    if by_n and kind != "1:1":
        raise StataError(198, "_n may only be used with 1:1")   # VERIFICAR
    keys = [] if by_n else s.expand_varlist(p.varlist)
    if not keys and not by_n:
        raise StataError(100, "varlist required")
    U = load_using(p.using)
    for k in keys:
        if not U.has(k):
            raise StataError(111, f"variable {k} not found")    # VERIFICAR (falta na using)
    if opts.get("keepusing") not in (None, True):
        from ..core.varlist import expand
        want = set(expand(U, str(opts["keepusing"]))) | set(keys)
        U.drop_vars([n for n in U.names if n not in want])
    update = bool(opts.get("update"))
    replace = bool(opts.get("replace"))
    if replace and not update:
        raise StataError(198, "option replace requires option update")   # VERIFICAR
    genname = "_merge"
    if opts.get("generate") not in (None, True):
        genname = str(opts["generate"]).strip()
    nogen = bool(opts.get("nogenerate"))
    if not nogen and ds.has(genname):
        raise StataError(110, f"variable {genname} already defined")
    if not nogen and U.has(genname):
        raise StataError(110, f"variable {genname} already defined")   # VERIFICAR
    force = bool(opts.get("force"))
    for k in keys:
        mv, uv = ds.get(k), U.get(k)
        if mv.is_string != uv.is_string:
            raise StataError(106, f"key variable {k} is {mv.vtype} in master but {uv.vtype} "
                                  f"in using data")   # VERIFICAR
    if kind in ("1:1", "1:m") and keys:
        _check_unique(ds, keys, "master")
    if kind in ("1:1", "m:1") and keys:
        _check_unique(U, keys, "using")

    # pares de observações
    if by_n:
        n = max(ds.nobs, U.nobs)
        r = np.arange(n)
        pm = np.where(r < ds.nobs, r, -1)
        pu = np.where(r < U.nobs, r, -1)
    else:
        gm, gu, kk = joint_keys([_key_col(ds.get(k)) for k in keys],
                                [_key_col(U.get(k)) for k in keys])
        pm, pu = _pairs(kind, gm, gu, kk)
    codes = np.where(pm < 0, 2, np.where(pu < 0, 1, 3)).astype(np.int64)

    # variáveis comuns (não chave)
    common = [v.name for v in U.vars if ds.has(v.name) and v.name not in keys]
    bad = {nm for nm in common if _mismatch(nm, ds.get(nm), U.get(nm), force, "merge")}
    new_vars: list[Variable] = []
    both = (pm >= 0) & (pu >= 0)
    upd_any = np.zeros(len(pm), dtype=bool)
    conflict = np.zeros(len(pm), dtype=bool)
    for mv in ds.vars:
        name = mv.name
        in_u = U.has(name) and name not in bad
        vtype = common_type(mv.vtype, U.get(name).vtype) if in_u else mv.vtype
        vals = _take(mv, pm, vtype)
        if in_u:
            uv = U.get(name)
            uvals = _take(uv, pu, vtype)
            only_u = pm < 0
            vals[only_u] = uvals[only_u]
            if update and name not in keys:
                tmp_m = Variable("_m", vtype, np.empty(0))
                tmp_m.raw = vals
                tmp_u = Variable("_u", vtype, np.empty(0))
                tmp_u.raw = uvals
                mm, um = _is_missing(tmp_m), _is_missing(tmp_u)
                fill = both & mm & ~um
                vals[fill] = uvals[fill]
                upd_any |= fill
                if tmp_m.is_string:
                    differ = np.array([a != b for a, b in zip(tmp_m.raw, tmp_u.raw)], dtype=bool)
                else:
                    differ = tmp_m.data != tmp_u.data
                clash = both & ~mm & ~um & differ
                conflict |= clash
                if replace:
                    vals[clash] = uvals[clash]
        nv = Variable(name, vtype, np.empty(0), fmt=mv.fmt if vtype == mv.vtype else
                      (default_format(vtype) if mv.fmt == default_format(mv.vtype) else mv.fmt),
                      label=mv.label, value_label=mv.value_label, notes=list(mv.notes))
        nv.raw = vals
        _fix_strings(nv)
        new_vars.append(nv)
    for uv in U.vars:
        if ds.has(uv.name):
            continue
        nv = Variable(uv.name, uv.vtype, np.empty(0), fmt=uv.fmt, label=uv.label,
                      value_label="" if opts.get("nolabel") else uv.value_label,
                      notes=[] if opts.get("nonotes") else list(uv.notes))
        nv.raw = _take(uv, pu)
        new_vars.append(nv)
    if update:
        codes = np.where(both & conflict, 5, np.where(both & upd_any, 4, codes))

    if opts.get("assert") not in (None, True):
        allowed = _results(str(opts["assert"]))
        if not np.isin(codes, sorted(allowed)).all():
            raise StataError(9, f"merge:  after merge, not all observations from "
                                f"{strip_outer_quotes(str(opts['assert']))}")   # VERIFICAR
    keep = np.ones(len(codes), dtype=bool)
    if opts.get("keep") not in (None, True):
        keep = np.isin(codes, sorted(_results(str(opts["keep"]))))

    new = Dataset()
    new.label = ds.label
    new.notes = list(ds.notes)
    new.chars = ds.chars
    new.value_labels = ds.value_labels
    new.filename, new.fullpath, new.timestamp = ds.filename, ds.fullpath, ds.timestamp
    new.nobs = int(keep.sum())
    for nv in new_vars:
        nv.raw = nv.raw[keep]
        new.vars.append(nv)
    codes = codes[keep]
    if not opts.get("nolabel"):
        used = {v.value_label for v in U.vars if v.value_label}
        _merge_value_labels(s, new, U, used)
    if not nogen:
        new.vars.append(Variable(genname, "byte", codes.astype(np.float64), value_label="_merge"))
        new.value_labels["_merge"] = dict(_MERGE_LABELS)
    new.sortlist = list(keys) if keys else []
    new.changed = True
    s.data = new
    if not opts.get("noreport"):
        _report(s, codes, genname, update)
    s.notify_state()


# ---------------------------------------------------------------------------
# joinby
# ---------------------------------------------------------------------------

@command("joinby")
def cmd_joinby(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.using is None:
        raise StataError(100, "using required")
    opts = match_options(p.options, {"update": 2, "replace": 3, "_merge": 6, "unmatched": 3, "nolabel": 5})
    ds = s.data
    U = load_using(p.using)
    keys = s.expand_varlist(p.varlist) if p.varlist.strip() else [n for n in ds.names if U.has(n)]
    if not keys:
        raise StataError(459, "no common variables")   # VERIFICAR
    for k in keys:
        if not U.has(k):
            raise StataError(111, f"variable {k} not found")
        if ds.get(k).is_string != U.get(k).is_string:
            raise StataError(106, f"variable {k} is {ds.get(k).vtype} in master but "
                                  f"{U.get(k).vtype} in using data")   # VERIFICAR
    unmatched = str(opts.get("unmatched") or "none").strip().lower()
    if unmatched not in ("none", "both", "master", "using"):
        raise StataError(198, "invalid unmatched() option")
    gname = str(opts["_merge"]).strip() if opts.get("_merge") not in (None, True) else "_merge"
    gm, gu, k = joint_keys([_key_col(ds.get(n)) for n in keys], [_key_col(U.get(n)) for n in keys])
    order_m = np.argsort(gm, kind="stable")
    order_u = np.argsort(gu, kind="stable")
    bm = np.searchsorted(gm[order_m], np.arange(k + 1))
    bu = np.searchsorted(gu[order_u], np.arange(k + 1))
    pm_l, pu_l = [], []
    for g in range(k):
        rm = order_m[bm[g]:bm[g + 1]]
        ru = order_u[bu[g]:bu[g + 1]]
        if len(rm) and len(ru):
            pm_l.append(np.repeat(rm, len(ru)))
            pu_l.append(np.tile(ru, len(rm)))
        elif len(rm) and unmatched in ("both", "master"):
            pm_l.append(rm)
            pu_l.append(np.full(len(rm), -1))
        elif len(ru) and unmatched in ("both", "using"):
            pm_l.append(np.full(len(ru), -1))
            pu_l.append(ru)
    pm = np.concatenate(pm_l).astype(np.int64) if pm_l else np.zeros(0, dtype=np.int64)
    pu = np.concatenate(pu_l).astype(np.int64) if pu_l else np.zeros(0, dtype=np.int64)
    both = (pm >= 0) & (pu >= 0)
    new = Dataset()
    new.label = ds.label
    new.value_labels = ds.value_labels
    new.nobs = len(pm)
    update, replace = bool(opts.get("update")), bool(opts.get("replace"))
    for mv in ds.vars:
        in_u = U.has(mv.name) and U.get(mv.name).is_string == mv.is_string
        vtype = common_type(mv.vtype, U.get(mv.name).vtype) if in_u else mv.vtype
        vals = _take(mv, pm, vtype)
        if in_u:
            uvals = _take(U.get(mv.name), pu, vtype)
            vals[pm < 0] = uvals[pm < 0]
            if mv.name not in keys and update:
                tm = Variable("_m", vtype, np.empty(0))
                tm.raw = vals
                tu = Variable("_u", vtype, np.empty(0))
                tu.raw = uvals
                fill = both & _is_missing(tm) & ~_is_missing(tu)
                if replace:
                    fill = both & ~_is_missing(tu)
                vals[fill] = uvals[fill]
        nv = Variable(mv.name, vtype, np.empty(0), fmt=mv.fmt, label=mv.label, value_label=mv.value_label)
        nv.raw = vals
        _fix_strings(nv)
        new.vars.append(nv)
    for uv in U.vars:
        if ds.has(uv.name):
            continue
        nv = Variable(uv.name, uv.vtype, np.empty(0), fmt=uv.fmt, label=uv.label,
                      value_label="" if opts.get("nolabel") else uv.value_label)
        nv.raw = _take(uv, pu)
        new.vars.append(nv)
    if not opts.get("nolabel"):
        _merge_value_labels(s, new, U, {v.value_label for v in U.vars if v.value_label})
    if unmatched != "none":
        if new.has(gname):
            raise StataError(110, f"variable {gname} already defined")
        codes = np.where(pm < 0, 2, np.where(pu < 0, 1, 3)).astype(np.float64)
        new.vars.append(Variable(gname, "byte", codes))
    new.changed = True
    s.data = new
    sorting.sort(new, keys)   # VERIFICAR ordem das observações
    s.notify_state()


# ---------------------------------------------------------------------------
# cross
# ---------------------------------------------------------------------------

@command("cross")
def cmd_cross(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.using is None or p.varlist.strip() or p.options.strip():
        raise StataError(100, "using required") if p.using is None else StataError(198, "invalid syntax")
    ds = s.data
    U = load_using(p.using)
    for uv in U.vars:
        if ds.has(uv.name):
            raise StataError(110, f"variable {uv.name} already defined")   # VERIFICAR
    nm, nu = ds.nobs, U.nobs
    rm = np.repeat(np.arange(nm), nu)
    ru = np.tile(np.arange(nu), nm)
    for v in ds.vars:
        v.raw = v.raw[rm]
    for uv in U.vars:
        nv = Variable(uv.name, uv.vtype, np.empty(0), fmt=uv.fmt, label=uv.label, value_label=uv.value_label)
        nv.raw = uv.raw[ru]
        ds.vars.append(nv)
    ds.nobs = nm * nu
    _merge_value_labels(s, ds, U, {v.value_label for v in U.vars if v.value_label})
    ds.sortlist = []
    ds.changed = True
    s.notify_state()
