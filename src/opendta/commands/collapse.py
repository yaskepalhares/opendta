"""collapse, contract, expand e fillin ([D] collapse, contract, expand, fillin).

collapse troca os dados por uma linha por grupo de by() com as estatísticas
pedidas:

    collapse [(stat)] varlist [(stat) ...] [weight] [if] [in] [, by(varlist) cw fast]

e também aceita `(stat) novo=velho`. O rótulo de cada variável nova é
"(stat) rótulo-ou-nome" (VERIFICAR). collapse não imprime nada.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import stats as S
from ..core.dataset import Dataset, Variable, check_name, default_format, smallest_type_for
from ..core.errors import StataError
from ..core.grouping import group_ids, first_index
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import eval_vector, plural, touse
from .registry import command
from .summarize import weights

if TYPE_CHECKING:
    from ..session import Session

_STATS = ("mean", "median", "sd", "semean", "sebinomial", "sepoisson", "sum", "rawsum",
          "count", "percent", "max", "min", "iqr", "first", "last", "firstnm", "lastnm")
_KEEP_TYPE = ("min", "max", "first", "last", "firstnm", "lastnm")
_STRING_OK = ("first", "last", "firstnm", "lastnm", "count")


def _stat_name(token: str) -> str:
    st = token.strip("() ").lower()
    if re.fullmatch(r"p\d{1,2}(\.\d+)?", st) and 0 < float(st[1:]) < 100:
        return st
    if st not in _STATS:
        raise StataError(198, f"{st} invalid statistic")   # VERIFICAR mensagem
    return st


def _parse_clist(s: "Session", text: str) -> list[tuple[str, str, str]]:
    """Lista de (estatística, variável nova, variável de origem)."""
    toks = re.findall(r"\([^)]*\)|=|[^\s=()]+", text)
    out: list[tuple[str, str, str]] = []
    stat = "mean"
    i = 0
    while i < len(toks):
        t = toks[i]
        if t.startswith("("):
            stat = _stat_name(t)
            i += 1
            continue
        if t == "=":
            raise StataError(198, "invalid syntax")
        if i + 1 < len(toks) and toks[i + 1] == "=":
            if i + 2 >= len(toks):
                raise StataError(198, "invalid syntax")
            new, src = t, s.expand_varlist(toks[i + 2])
            if len(src) != 1:
                raise StataError(103, "too many variables specified")
            check_name(new)
            out.append((stat, new, src[0]))
            i += 3
            continue
        for name in s.expand_varlist(t):
            out.append((stat, name, name))
        i += 1
    if not out:
        raise StataError(100, "varlist required")
    seen: set[str] = set()
    for _, new, _ in out:
        if new in seen:
            # VERIFICAR mensagem
            raise StataError(198, f"{new} specified more than once")
        seen.add(new)
    return out


class _Groups:
    """Observações usadas, agrupadas: `g` é o id do grupo (0..k-1) de cada
    observação usada, na ordem dos dados."""

    def __init__(self, ids: np.ndarray, k: int, w: np.ndarray | None, wtype: str):
        self.rows = np.flatnonzero(ids >= 0)
        self.g = ids[self.rows]
        self.k = k
        self.w = None if w is None else w[self.rows]
        self.wtype = wtype


def _compute(stat: str, x: np.ndarray, G: _Groups) -> np.ndarray:
    """Estatística `stat` de x (alinhado com G.rows) para cada grupo."""
    k = G.k
    out = np.full(k, M.SYSMISS)
    ok = x < M.SYSMISS
    g = G.g[ok]
    xv = x[ok]
    w = np.ones(len(xv)) if G.w is None else G.w[ok]
    n_obs = np.bincount(g, minlength=k).astype(np.float64)
    sum_w = np.bincount(g, weights=w, minlength=k)
    if G.wtype in ("aweight", "pweight"):
        # pesos normalizados para somar o nº de observações do grupo
        with np.errstate(all="ignore"):
            wn = w * (n_obs / sum_w)[g]
        N = n_obs
    else:
        wn = w
        N = sum_w
    has = n_obs > 0
    if stat in ("first", "last", "firstnm", "lastnm"):
        rows_g = G.g if stat in ("first", "last") else g
        vals = x if stat in ("first", "last") else xv
        idx = np.arange(len(rows_g))
        pick = np.full(k, -1, dtype=np.int64)
        order = idx[::-1] if stat.startswith("first") else idx
        pick[rows_g[order]] = order
        got = pick >= 0
        out[got] = vals[pick[got]]
        return out
    if stat == "count":
        return N.copy() if G.w is not None else n_obs
    if stat == "percent":
        total = N.sum()
        return N / total * 100 if total else np.zeros(k)
    if stat == "rawsum":
        return np.bincount(g, weights=xv, minlength=k)
    if stat == "sum":
        return np.bincount(g, weights=wn * xv, minlength=k)
    with np.errstate(all="ignore"):
        mean = np.bincount(g, weights=wn * xv, minlength=k) / np.bincount(g, weights=wn, minlength=k)
    if stat == "mean":
        out[has] = mean[has]
        return out
    if stat in ("sd", "semean"):
        dev2 = wn * (xv - mean[g]) ** 2
        m2 = np.bincount(g, weights=dev2, minlength=k)
        ok2 = has & (N > 1)
        sd = np.full(k, M.SYSMISS)
        sd[ok2] = np.sqrt(m2[ok2] / (N[ok2] - 1))
        if stat == "sd":
            return sd
        out[ok2] = sd[ok2] / np.sqrt(N[ok2])
        return out
    if stat == "sebinomial":
        ok2 = has & (mean >= 0) & (mean <= 1)
        out[ok2] = np.sqrt(mean[ok2] * (1 - mean[ok2]) / N[ok2])
        return out
    if stat == "sepoisson":
        ok2 = has & (mean >= 0)
        out[ok2] = np.sqrt(mean[ok2] / N[ok2])
        return out
    if stat in ("min", "max"):
        if stat == "min":
            acc = np.full(k, np.inf)
            np.minimum.at(acc, g, xv)
        else:
            acc = np.full(k, -np.inf)
            np.maximum.at(acc, g, xv)
        out[has] = acc[has]
        return out
    # percentis (median, iqr, p#): ordena por (grupo, valor)
    order = np.lexsort((xv, g))
    gs, xs, ws = g[order], xv[order], w[order]
    bounds = np.searchsorted(gs, np.arange(k + 1))
    for j in range(k):
        a, b = bounds[j], bounds[j + 1]
        if a == b:
            continue
        wj = None if G.w is None else ws[a:b]
        if stat == "median":
            out[j] = S.percentile(xs[a:b], 50, wj)
        elif stat == "iqr":
            out[j] = S.percentile(xs[a:b], 75, wj) - S.percentile(xs[a:b], 25, wj)
        else:
            out[j] = S.percentile(xs[a:b], float(stat[1:]), wj)
    return out


def _string_stat(stat: str, x: np.ndarray, G: _Groups) -> np.ndarray | None:
    k = G.k
    if stat == "count":
        return None
    out = np.array([""] * k, dtype=object)
    nonempty = np.array([t != "" for t in x], dtype=bool)
    sel = np.ones(len(x), dtype=bool) if stat in ("first", "last") else nonempty
    idx = np.flatnonzero(sel)
    order = idx[::-1] if stat.startswith("first") else idx
    pick = np.full(k, -1, dtype=np.int64)
    pick[G.g[order]] = order
    got = pick >= 0
    out[got] = x[pick[got]]
    return out


@command("collapse")
def cmd_collapse(s: "Session", args: str) -> None:
    p = parse_standard(args)
    clist_text = p.varlist if p.exp is None else f"{p.varlist} = {p.exp}"
    opts = match_options(p.options, {"by": 2, "cw": 2, "fast": 4})
    ds = s.data
    items = _parse_clist(s, clist_text)
    by = s.expand_varlist(str(opts["by"])) if opts.get("by") not in (None, True) else []
    for _, new, _ in items:
        if new in by:
            raise StataError(198, f"{new} may not be both a by() variable and a target")  # VERIFICAR
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight", "iweight", "pweight"))
    if opts.get("cw"):
        for _, _, src in items:
            v = ds.get(src)
            mask &= (v.data < M.SYSMISS) if not v.is_string else np.array([t != "" for t in v.raw], dtype=bool)
    if not mask.any():
        raise StataError(2000, "no observations")
    ids, k = group_ids(ds, by, mask)
    G = _Groups(ids, k, w, wtype)
    first = first_index(ids, k)
    new = Dataset()
    new.nobs = k
    new.label = ds.label
    used_labels: set[str] = set()
    for name in by:
        v = ds.get(name)
        nv = Variable(name, v.vtype, np.empty(0), fmt=v.fmt, label=v.label, value_label=v.value_label)
        nv.raw = v.raw[first].copy()
        new.vars.append(nv)
        if v.value_label:
            used_labels.add(v.value_label)
    for stat, name, src in items:
        v = ds.get(src)
        lab = f"({stat}) {v.label or src}"     # VERIFICAR
        if v.is_string:
            if stat not in _STRING_OK:
                raise StataError(109, "type mismatch")
            res = _string_stat(stat, v.raw[G.rows], G)
            if res is not None:
                nv = Variable(name, v.vtype, res, fmt=v.fmt, label=lab)
                new.vars.append(nv)
                continue
            x = np.array([M.SYSMISS if t == "" else 0.0 for t in v.raw[G.rows]])
        else:
            x = np.asarray(v.data, dtype=np.float64)[G.rows]
        res = _compute(stat, x, G)
        fmt = v.fmt if not v.is_string else ""
        vlab = ""
        if stat in _KEEP_TYPE:
            vtype = v.vtype
            vlab = v.value_label
        elif stat == "count":
            vtype = "long" if np.all(res == np.trunc(res)) else "double"   # VERIFICAR
            fmt = ""
        elif stat in ("sum", "rawsum"):
            vtype, fmt = "double", ""                                       # VERIFICAR
        else:
            vtype = "double" if v.vtype == "double" else "float"           # VERIFICAR
            if v.vtype in ("byte", "int", "long"):
                fmt = ""
        nv = Variable(name, vtype, res, fmt=fmt or default_format(vtype), label=lab, value_label=vlab)
        if vlab:
            used_labels.add(vlab)
        new.vars.append(nv)
    new.value_labels = {n: dict(t) for n, t in ds.value_labels.items() if n in used_labels}
    new.sortlist = list(by)
    new.changed = True
    s.data = new
    s.notify_state()


# ---------------------------------------------------------------------------
# contract
# ---------------------------------------------------------------------------

@command("contract")
def cmd_contract(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"freq": 1, "cfreq": 2, "percent": 1, "cpercent": 2,
                                     "float": 2, "format": 3, "zero": 1, "nomiss": 5})
    ds = s.data
    names = s.expand_varlist(p.varlist)
    if not names:
        raise StataError(100, "varlist required")
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    w, wtype, mask = weights(s, p, mask, ("fweight",))
    if opts.get("nomiss"):
        for nm in names:
            v = ds.get(nm)
            mask &= (v.data < M.SYSMISS) if not v.is_string else np.array([t != "" for t in v.raw], dtype=bool)
    out_names = {}
    for key, default in (("freq", "_freq"), ("cfreq", None), ("percent", None), ("cpercent", None)):
        val = opts.get(key)
        if val is True:
            raise StataError(198, "invalid syntax")
        if val or default:
            nm = str(val).strip() if val else default
            check_name(nm)
            if nm in names:
                raise StataError(110, f"variable {nm} already defined")
            out_names[key] = nm
    ids, k = group_ids(ds, names, mask)
    rows = np.flatnonzero(ids >= 0)
    ww = np.ones(len(rows)) if w is None else w[rows]
    freq = np.bincount(ids[rows], weights=ww, minlength=k)
    first = first_index(ids, k)
    keys = {nm: ds.get(nm).raw[first] for nm in names}
    if opts.get("zero") and k:
        levels = [np.unique(ds.get(nm).raw[rows].astype(str) if ds.get(nm).is_string
                            else ds.get(nm).raw[rows]) for nm in names]
        grids = np.meshgrid(*levels, indexing="ij")
        combos = [g.ravel() for g in grids]
        have = {tuple(keys[nm][j] for nm in names): freq[j] for j in range(k)}
        freq = np.array([have.get(tuple(c[j] for c in combos), 0.0) for j in range(len(combos[0]))])
        keys = {}
        for nm, c in zip(names, combos):
            v = ds.get(nm)
            keys[nm] = np.array(c, dtype=object) if v.is_string else np.asarray(c, dtype=v.raw.dtype)
        k = len(freq)
    new = Dataset()
    new.nobs = k
    new.label = ds.label
    used = set()
    for nm in names:
        v = ds.get(nm)
        nv = Variable(nm, v.vtype, np.empty(0), fmt=v.fmt, label=v.label, value_label=v.value_label)
        nv.raw = keys[nm].copy()
        new.vars.append(nv)
        if v.value_label:
            used.add(v.value_label)
    total = freq.sum()
    ptype = "float" if opts.get("float") else "double"
    pfmt = strip_outer_quotes(str(opts["format"])) if opts.get("format") not in (None, True) else "%8.2f"
    ftype = smallest_type_for(freq) if k else "byte"   # VERIFICAR tipo de _freq
    if ftype in ("float", "double"):
        ftype = "long"
    columns = {
        "freq": (freq, ftype, "", "Frequency"),
        "cfreq": (np.cumsum(freq), ftype, "", "Cumulative frequency"),
        "percent": (freq / total * 100 if total else freq * 0, ptype, pfmt, "Percent"),
        "cpercent": (np.cumsum(freq) / total * 100 if total else freq * 0, ptype, pfmt,
                     "Cumulative percent"),
    }
    for key in ("freq", "cfreq", "percent", "cpercent"):    # VERIFICAR ordem e rótulos
        if key in out_names:
            vals, vt, fmt, lab = columns[key]
            new.vars.append(Variable(out_names[key], vt, vals, fmt=fmt, label=lab))
    new.value_labels = {n: dict(t) for n, t in ds.value_labels.items() if n in used}
    new.sortlist = list(names)
    new.changed = True
    s.data = new
    s.notify_state()


# ---------------------------------------------------------------------------
# expand
# ---------------------------------------------------------------------------

@command("expand")
def cmd_expand(s: "Session", args: str) -> None:
    text = args.strip()
    if text.startswith("="):
        text = text[1:]
    p = parse_standard(text)
    opts = match_options(p.options, {"generate": 1})
    ds = s.data
    if not p.varlist.strip():
        raise StataError(198, "invalid syntax")
    gen = str(opts["generate"]).strip() if opts.get("generate") not in (None, True) else ""
    if gen:
        if ds.has(gen):
            raise StataError(110, f"variable {gen} already defined")
        check_name(gen)
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    k = np.broadcast_to(np.asarray(eval_vector(s, p.varlist), dtype=np.float64), (ds.nobs,))
    # VERIFICAR: valores não inteiros são truncados; < 1 ou missing mantêm a obs.
    copies = np.where(mask & (k < M.SYSMISS) & (k >= 1), np.trunc(np.minimum(k, 2**31)), 1).astype(np.int64)
    extra = np.repeat(np.arange(ds.nobs), copies - 1)
    n0 = ds.nobs
    if len(extra):
        for v in ds.vars:
            v.raw = np.concatenate([v.raw, v.raw[extra]])
        ds.nobs = n0 + len(extra)
        ds.sortlist = []
        ds.changed = True
    if gen:
        flag = np.zeros(ds.nobs)
        flag[n0:] = 1
        ds.add(Variable(gen, "byte", flag))
    s.output.write(f"({plural(len(extra), 'observation')} created)\n", "text")
    s.notify_state()


# ---------------------------------------------------------------------------
# fillin
# ---------------------------------------------------------------------------

@command("fillin")
def cmd_fillin(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.options.strip() or p.if_ or p.in_:
        raise StataError(198, "invalid syntax")
    ds = s.data
    names = s.expand_varlist(p.varlist)
    if len(names) < 2:
        raise StataError(102, "too few variables specified")   # VERIFICAR
    if ds.has("_fillin"):
        raise StataError(110, "variable _fillin already defined")
    levels = []
    for nm in names:
        v = ds.get(nm)
        if v.is_string:
            levels.append(np.unique(v.raw.astype(str)).astype(object))
        else:
            levels.append(np.unique(v.raw))
    present = set(zip(*[ds.get(nm).raw.astype(str) if ds.get(nm).is_string else ds.get(nm).raw
                        for nm in names])) if ds.nobs else set()
    grids = np.meshgrid(*[np.arange(len(lv)) for lv in levels], indexing="ij")
    combos = [lv[g.ravel()] for lv, g in zip(levels, grids)]
    missing_rows = [j for j in range(len(combos[0]) if combos else 0)
                    if tuple(c[j] for c in combos) not in present]
    n0 = ds.nobs
    add = len(missing_rows)
    if add:
        ds.set_obs(n0 + add)
        for nm, c in zip(names, combos):
            v = ds.get(nm)
            v.raw[n0:] = c[missing_rows]
    flag = np.zeros(ds.nobs)
    flag[n0:] = 1
    ds.add(Variable("_fillin", "byte", flag))
    from ..core import sorting
    sorting.sort(ds, names)
    s.notify_state()
