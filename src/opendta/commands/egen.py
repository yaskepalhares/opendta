"""egen ([D] egen): funções de grupo e de linha.

    egen [type] newvar = fcn(arguments) [if] [in] [, by(varlist) options]

Funções de grupo (calculadas dentro de by): count, mean, sd, total/sum,
min, max, median, mode, pctile, iqr, kurt, skew, mad, mdev, std, rank.
Outras: group, tag, seq, fill, cut, concat, ends, diff.
Funções de linha: rowmean, rowtotal, rowmin, rowmax, rowsd, rownonmiss,
rowmiss, rowfirst, rowlast, rowmedian, rowpctile, anycount, anymatch,
anyvalue.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import stats as S
from ..core.grouping import group_ids
from ..core.dataset import Variable, check_name, smallest_type_for, str_type_for
from ..core.errors import StataError
from ..lang.functions import _string
from ..lang.syntax import Parsed, parse_options, parse_standard
from ..lang.words import parse_numlist, strip_outer_quotes
from ._util import eval_vector, plural, touse
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_TYPES = ("byte", "int", "long", "float", "double")


def _group_ids(s: "Session", by: list[str], mask: np.ndarray) -> np.ndarray:
    """Id de grupo (0..k-1) por observação, em ordem crescente das chaves;
    -1 fora da máscara."""
    return group_ids(s.data, by, mask)[0]


def _by_reduce(x: np.ndarray, ids: np.ndarray, fn) -> np.ndarray:
    """Aplica fn(valores do grupo) e espalha o resultado pelo grupo."""
    out = np.full(len(x), M.SYSMISS)
    for g in np.unique(ids[ids >= 0]):
        sel = ids == g
        out[sel] = fn(x[sel])
    return out


def _vals(s: "Session", expr: str) -> np.ndarray:
    v = eval_vector(s, expr)
    if isinstance(v, str) or (isinstance(v, np.ndarray) and v.dtype == object):
        raise StataError(109, "type mismatch")
    return np.broadcast_to(np.asarray(v, dtype=np.float64), (s.data.nobs,)).copy()


def _nonmiss(f):
    def g(x):
        x = x[x < M.SYSMISS]
        return f(x) if len(x) else M.SYSMISS
    return g


def _mode(x: np.ndarray, how: str) -> float:
    x = x[x < M.SYSMISS]
    if len(x) == 0:
        return M.SYSMISS
    vals, counts = np.unique(x, return_counts=True)
    top = vals[counts == counts.max()]
    if len(top) > 1 and how == "":
        return M.SYSMISS
    return float(top.min() if how in ("minmode", "") else top.max())


def _rank(x: np.ndarray, how: str) -> np.ndarray:
    ok = x < M.SYSMISS
    out = np.full(len(x), M.SYSMISS)
    v = x[ok]
    if len(v) == 0:
        return out
    from scipy.stats import rankdata
    if how == "field":
        r = rankdata(-v, method="min")
    elif how == "track":
        r = rankdata(v, method="min")
    elif how == "unique":
        r = rankdata(v, method="ordinal")
    else:
        r = rankdata(v, method="average")
    out[ok] = r
    return out


def _fill(nums: list[float], n: int) -> np.ndarray:
    """fill(numlist): progressão ou padrão que se repete. Procura o menor
    período p tal que nums[i+p] - nums[i] seja constante (d); o valor da
    obs. k é nums[k % p] + (k // p)·d. Ex.: 1 2 3 → 1 2 3 4...;
    10 20 10 20 → repete; 1 1 2 2 → 1 1 2 2 3 3..."""
    if not nums:
        raise StataError(198, "invalid syntax")
    if len(nums) == 1:
        return np.full(n, float(nums[0]))
    a = np.asarray(nums, dtype=np.float64)
    for p in range(1, len(a)):
        d = a[p:] - a[:-p]
        if np.allclose(d, d[0]):
            k = np.arange(n)
            return a[k % p] + (k // p) * d[0]
    k = np.arange(n)
    return a[k % len(a)]


def _text_opt(o: dict, key: str, default: str) -> str:
    v = o.get(key)
    return default if v in (None, True) else strip_outer_quotes(str(v))


def _row_vars(s: "Session", arg: str) -> list[Variable]:
    names = s.expand_varlist(arg)
    return [s.data.get(n) for n in names]


def _matrix_rows(vars_: list[Variable]) -> np.ndarray:
    return np.column_stack([v.data for v in vars_]) if vars_ else np.zeros((0, 0))


@command("egen", byable=True)
def cmd_egen(s: "Session", args: str) -> None:
    m = re.match(r"^\s*(?:(byte|int|long|float|double|str\d+|strL)\s+)?([A-Za-z_]\w*)\s*=\s*"
                 r"([A-Za-z_]\w*)\s*\((.*)$", args, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    vtype, name, fn = m.group(1), m.group(2), m.group(3)
    rest = m.group(4)
    # argumento: até o parêntese que fecha
    depth, k = 1, 0
    while k < len(rest) and depth:
        if rest[k] == "(":
            depth += 1
        elif rest[k] == ")":
            depth -= 1
        elif rest[k] == '"':
            q = rest.find('"', k + 1)
            k = len(rest) - 1 if q == -1 else q
        k += 1
    if depth:
        raise StataError(198, "invalid syntax")
    arg, tail = rest[:k - 1], rest[k:]
    p = parse_standard(tail)
    o = {n.lower(): (a if a is not None else True) for n, a in parse_options(p.options)} \
        if p.options.strip() else {}
    ds = s.data
    if ds.has(name):
        raise StataError(110, f"variable {name} already defined")
    check_name(name)
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    by = s.expand_varlist(str(o["by"])) if o.get("by") else list(getattr(s, "_by_keys", []))
    ids = _group_ids(s, by, mask)
    n = ds.nobs
    is_string = False
    default_type = s.settings.get("type", "float")
    fn_l = fn.lower()

    if fn_l in ("count", "mean", "sd", "total", "sum", "min", "max", "median", "mode", "pctile",
                "iqr", "kurt", "skew", "mad", "mdev", "std", "rank", "diff"):
        if fn_l == "diff":
            vs = _row_vars(s, arg)
            X = _matrix_rows(vs)
            res = np.where(np.all(X == X[:, :1], axis=1), 0.0, 1.0)
            res[~mask] = M.SYSMISS
            vtype = vtype or "byte"
        else:
            x = _vals(s, arg)
            x = np.where(mask, x, M.SYSMISS)
            if fn_l == "count":
                res = _by_reduce(x, ids, lambda a: float((a < M.SYSMISS).sum()))
                vtype = vtype or "long"   # VERIFICAR
            elif fn_l == "mean":
                res = _by_reduce(x, ids, _nonmiss(lambda a: float(a.mean())))
            elif fn_l == "sd":
                res = _by_reduce(x, ids, _nonmiss(lambda a: float(a.std(ddof=1)) if len(a) > 1 else M.SYSMISS))
            elif fn_l in ("total", "sum"):
                keep_missing = bool(o.get("missing"))
                res = _by_reduce(x, ids, lambda a: (M.SYSMISS if keep_missing and not (a < M.SYSMISS).any()
                                                    else float(a[a < M.SYSMISS].sum())))
            elif fn_l == "min":
                res = _by_reduce(x, ids, _nonmiss(lambda a: float(a.min())))
            elif fn_l == "max":
                res = _by_reduce(x, ids, _nonmiss(lambda a: float(a.max())))
            elif fn_l == "median":
                res = _by_reduce(x, ids, _nonmiss(lambda a: S.percentile(np.sort(a), 50)))
            elif fn_l == "pctile":
                q = float(o.get("p", 50))
                res = _by_reduce(x, ids, _nonmiss(lambda a: S.percentile(np.sort(a), q)))
            elif fn_l == "iqr":
                res = _by_reduce(x, ids, _nonmiss(lambda a: S.percentile(np.sort(a), 75)
                                                  - S.percentile(np.sort(a), 25)))
            elif fn_l == "mode":
                how = "minmode" if o.get("minmode") else "maxmode" if o.get("maxmode") else ""
                res = _by_reduce(x, ids, lambda a: _mode(a, how))
            elif fn_l == "kurt":
                res = _by_reduce(x, ids, _nonmiss(lambda a: S.moments(a).kurtosis))
            elif fn_l == "skew":
                res = _by_reduce(x, ids, _nonmiss(lambda a: S.moments(a).skewness))
            elif fn_l == "mad":
                res = _by_reduce(x, ids, _nonmiss(
                    lambda a: S.percentile(np.sort(np.abs(a - S.percentile(np.sort(a), 50))), 50)))
            elif fn_l == "mdev":
                res = _by_reduce(x, ids, _nonmiss(lambda a: float(np.abs(a - a.mean()).mean())))
            elif fn_l == "std":
                mean = _by_reduce(x, ids, _nonmiss(lambda a: float(a.mean())))
                sd = _by_reduce(x, ids, _nonmiss(lambda a: float(a.std(ddof=1)) if len(a) > 1 else M.SYSMISS))
                target_mean = float(o.get("mean", 0))
                target_sd = float(o.get("sd", o.get("std", 1)))
                ok = (x < M.SYSMISS) & (sd < M.SYSMISS) & (sd > 0)
                res = np.full(n, M.SYSMISS)
                res[ok] = target_mean + target_sd * (x[ok] - mean[ok]) / sd[ok]
            else:   # rank
                how = "field" if o.get("field") else "track" if o.get("track") else \
                    "unique" if o.get("unique") else ""
                res = np.full(n, M.SYSMISS)
                for g in np.unique(ids[ids >= 0]):
                    sel = ids == g
                    res[sel] = _rank(x[sel], how)
            res = np.where(ids >= 0, res, M.SYSMISS)
    elif fn_l == "group":
        vs = s.expand_varlist(arg)
        sel = mask.copy()
        if not o.get("missing"):
            for nm in vs:
                v = ds.get(nm)
                sel &= S.valid(v.data) if not v.is_string else np.array([x != "" for x in v.raw])
        gids = _group_ids(s, vs, sel)
        res = np.where(gids >= 0, gids + 1.0, M.SYSMISS)
        vtype = vtype or ("int" if gids.max(initial=0) < 32000 else "long")   # VERIFICAR
        if o.get("label"):
            _group_labels(s, name, vs, gids)
    elif fn_l == "tag":
        vs = s.expand_varlist(arg)
        sel = mask.copy()
        if not o.get("missing"):
            for nm in vs:
                v = ds.get(nm)
                sel &= S.valid(v.data) if not v.is_string else np.array([x != "" for x in v.raw])
        gids = _group_ids(s, vs, sel)
        res = np.zeros(n)
        seen: set[int] = set()
        for i in range(n):
            g = int(gids[i])
            if g >= 0 and g not in seen:
                seen.add(g)
                res[i] = 1.0
        vtype = vtype or "byte"
    elif fn_l == "seq":
        start = int(o.get("from", 1))
        to = o.get("to")
        block = int(o.get("block", 1))
        res = np.full(n, M.SYSMISS)
        for g in np.unique(ids[ids >= 0]):
            idx = np.flatnonzero(ids == g)
            k = np.arange(len(idx)) // block
            if to is not None:
                span = int(to) - start
                vals = start + (k % (abs(span) + 1)) * (1 if span >= 0 else -1)
            else:
                vals = start + k
            res[idx] = vals
    elif fn_l == "fill":
        res = _fill(parse_numlist(arg), n)
    elif fn_l == "cut":
        x = _vals(s, arg)
        res = np.full(n, M.SYSMISS)
        if o.get("at"):
            at = parse_numlist(str(o["at"]))
            k = np.searchsorted(at, x, side="right") - 1
            ok = mask & (x < M.SYSMISS) & (k >= 0) & (k < len(at) - 1)
            res[ok] = np.array(at)[k[ok]] if not o.get("icodes") else k[ok]
        elif o.get("group"):
            g = int(o["group"])
            ok = mask & (x < M.SYSMISS)
            xs = np.sort(x[ok])
            cuts = [S.percentile(xs, 100 * j / g) for j in range(1, g)]
            k = np.searchsorted(cuts, x[ok], side="right")
            res[ok] = k
        else:
            raise StataError(198, "at() or group() required")
    elif fn_l == "concat":
        vs = _row_vars(s, arg)
        punct = _text_opt(o, "punct", "")
        fmt = _text_opt(o, "format", "")
        maxlen = int(o["maxlength"]) if o.get("maxlength") not in (None, True) else None
        from .inspect import cell_text
        texts = []
        for i in range(n):
            if not mask[i]:
                texts.append("")
                continue
            parts = []
            for v in vs:
                if v.is_string:
                    parts.append(v.raw[i])
                elif o.get("decode") and v.value_label:
                    parts.append(cell_text(ds, v, i, use_labels=True).strip())
                else:
                    # VERIFICAR: formato padrão de concat para variáveis numéricas
                    parts.append(_string(v.value(i), fmt) if fmt else _string(v.value(i)))
            t = punct.join(parts)
            texts.append(t[:maxlen] if maxlen else t)
        res = texts
        is_string = True
    elif fn_l == "ends":
        v = _row_vars(s, arg)[0]
        punct = _text_opt(o, "punct", " ")
        texts = []
        for i in range(n):
            t = v.raw[i] if mask[i] else ""
            t = t.strip() if punct == " " else t
            if o.get("tail"):
                texts.append(t.split(punct, 1)[1] if punct in t else "")
            elif o.get("last"):
                texts.append(t.rsplit(punct, 1)[-1])
            else:
                texts.append(t.split(punct, 1)[0])
        res = texts
        is_string = True
    elif fn_l.startswith("row") or fn_l.startswith("any"):
        vs = _row_vars(s, arg)
        X = _matrix_rows(vs).astype(np.float64)
        miss = X >= M.SYSMISS
        Xn = np.where(miss, np.nan, X)
        with np.errstate(all="ignore"):
            if fn_l == "rowmean":
                res = np.nanmean(Xn, axis=1)
            elif fn_l == "rowtotal":
                res = np.nansum(Xn, axis=1)
                if o.get("missing"):
                    res[miss.all(1)] = np.nan
            elif fn_l == "rowmin":
                res = np.nanmin(Xn, axis=1)
            elif fn_l == "rowmax":
                res = np.nanmax(Xn, axis=1)
            elif fn_l == "rowsd":
                res = np.nanstd(Xn, axis=1, ddof=1)
            elif fn_l == "rownonmiss":
                res = (~miss).sum(1).astype(float)
            elif fn_l == "rowmiss":
                res = miss.sum(1).astype(float)
            elif fn_l == "rowmedian":
                res = np.array([S.percentile(np.sort(r[~np.isnan(r)]), 50) if (~np.isnan(r)).any()
                                else np.nan for r in Xn])
            elif fn_l == "rowpctile":
                q = float(o.get("p", 50))
                res = np.array([S.percentile(np.sort(r[~np.isnan(r)]), q) if (~np.isnan(r)).any()
                                else np.nan for r in Xn])
            elif fn_l == "rowfirst":
                res = np.array([r[~np.isnan(r)][0] if (~np.isnan(r)).any() else np.nan for r in Xn])
            elif fn_l == "rowlast":
                res = np.array([r[~np.isnan(r)][-1] if (~np.isnan(r)).any() else np.nan for r in Xn])
            elif fn_l in ("anycount", "anymatch", "anyvalue"):
                vals = parse_numlist(str(o.get("values", "")))
                hit = np.isin(X, vals)
                if fn_l == "anycount":
                    res = hit.sum(1).astype(float)
                elif fn_l == "anymatch":
                    res = hit.any(1).astype(float)
                else:
                    if len(vs) != 1:
                        raise StataError(103, "too many variables specified")
                    res = np.where(hit[:, 0], X[:, 0], np.nan)
            else:
                raise StataError(133, f"unknown egen function {fn}()")
        res = np.where(np.isnan(res) | ~mask, M.SYSMISS, res)
        if fn_l in ("anycount", "anymatch", "rownonmiss", "rowmiss"):
            vtype = vtype or "byte"   # VERIFICAR
    else:
        raise StataError(133, f"unknown egen function {fn}()")

    if is_string:
        longest = max((len(t.encode("utf-8")) for t in res), default=1) or 1
        nv = Variable(name, vtype if vtype and vtype.startswith("str") else str_type_for(longest),
                      np.array(res, dtype=object))
        ds.add(nv)
        nmiss = sum(1 for t in res if t == "")
    else:
        res = np.asarray(res, dtype=np.float64)
        t = vtype or default_type
        if t not in _TYPES:
            raise StataError(109, "type mismatch")
        if t in ("byte", "int", "long") and not vtype:
            t = max(t, smallest_type_for(res), key=lambda z: _TYPES.index(z))
        nv = Variable(name, t, res)
        ds.add(nv)
        nmiss = int((nv.data >= M.SYSMISS).sum())
    if nmiss:
        s.output.write(f"({plural(nmiss, 'missing value')} generated)\n", "text")
    s.notify_state()


def _group_labels(s: "Session", name: str, vs: list[str], gids: np.ndarray) -> None:
    from .inspect import cell_text
    ds = s.data
    lab: dict[int, str] = {}
    for i in range(ds.nobs):
        g = int(gids[i])
        if g >= 0 and g + 1 not in lab:
            lab[g + 1] = " ".join(cell_text(ds, ds.get(v), i) for v in vs)
    ds.value_labels[name] = lab
