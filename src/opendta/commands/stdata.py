"""Survival-time data manipulation: stsplit, stjoin, stfill, stgen, stvary, stbase and snapspan
([ST] stsplit, [ST] stjoin, [ST] stfill, [ST] stgen, [ST] stvary, [ST] stbase, [ST] snapspan).

The array work is in core/survival.py. Details that the manuals do not show carry the VERIFY mark
(compat/do/1901_survival.do)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import survival as SV
from ..core.dataset import Variable
from ..core.errors import StataError
from ..core.grouping import group_ids
from ..estimation.st import STInfo, header, require_st
from ..lang.syntax import match_options, parse_standard
from ..lang.words import parse_numlist
from ._util import eval_vector, touse

if TYPE_CHECKING:
    from ..core.dataset import Dataset
    from ..session import Session

SYS = M.SYSMISS
_ST_VARS = ("_st", "_d", "_t", "_t0")


def _need_id(info: STInfo, cmd: str) -> None:
    if not info.idvar:
        raise StataError(459, f"{cmd} requires that the data be stset with id()")   # VERIFY


def _take(ds: "Dataset", src: np.ndarray) -> None:
    """Rebuilds the observations as ds[src] (repeating or reordering rows)."""
    for v in ds.vars:
        v.raw = v.raw[src]
    ds.nobs = len(src)
    ds.esample = None
    ds.sortlist = []
    ds._touch()


def _put(ds: "Dataset", name: str, values: np.ndarray, rows: np.ndarray) -> None:
    var = ds.get(name)
    full = np.array(var.data, dtype=np.float64)
    full[rows] = values
    ds.set_numeric(var, full)


def _origin(s: "Session", info: STInfo, n: int) -> np.ndarray:
    """The origin of analysis time of each observation, in the units of the time variable."""
    if info.origin and info.origin != "min":
        return np.broadcast_to(np.asarray(eval_vector(s, info.origin), dtype=np.float64), (n,)).copy()
    return np.zeros(n)


def _st_arrays(s: "Session"):
    ds = s.data
    return (ds.get("_st").data == 1, ds.get("_t0").data.astype(np.float64),
            ds.get("_t").data.astype(np.float64), ds.get("_d").data.astype(np.float64))


def _subjects(s: "Session", info: STInfo, rows: np.ndarray):
    """Rows (sorted by subject, entry and exit time), their subject numbers and the bounds of
    each subject's run."""
    ds = s.data
    ok = np.zeros(ds.nobs, dtype=bool)
    ok[rows] = True
    if info.idvar:
        gid_all = group_ids(ds, [info.idvar], ok)[0]
    else:
        gid_all = np.where(ok, np.arange(ds.nobs), -1)
    t0 = ds.get("_t0").data.astype(np.float64)
    t = ds.get("_t").data.astype(np.float64)
    gid = gid_all[rows]
    o = SV.time_order(gid, t0[rows], t[rows])
    rows_s = rows[o]
    gid_s = gid[o]
    return rows_s, gid_s, SV.group_bounds(gid_s)


# ---------------------------------------------------------------------------
# stsplit
# ---------------------------------------------------------------------------

def _after_ref(s: "Session", info: STInfo, spec: str, rows: np.ndarray, t: np.ndarray, org: np.ndarray):
    """Reference point of each record for after(): at() values are measured from it, in analysis
    time units."""
    m = re.match(r"^\s*(time|t)\s*=\s*(.*)$", spec)
    unit, text = (m.group(1), m.group(2)) if m else ("time", spec)
    kind = "exp"
    m = re.match(r"^\s*(min|asis)\s*\((.*)\)\s*$", text)
    if m:
        kind, text = m.group(1), m.group(2)
    n = s.data.nobs
    v = np.broadcast_to(np.asarray(eval_vector(s, text), dtype=np.float64), (n,)).copy()
    if kind == "min":
        gid = group_ids(s.data, [info.idvar], None)[0]
        for g in np.unique(gid[rows]):
            sel = rows[gid[rows] == g]
            good = v[sel][v[sel] < SYS]
            v[sel] = good.min() if len(good) else SYS
    if unit == "time":
        v = (v - org) / info.scale
    return v


def cmd_stsplit(s: "Session", args: str) -> None:
    ds = s.data
    info = require_st(s)
    _need_id(info, "stsplit")
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) > 1:
        raise StataError(103, "too many variables specified")
    o = match_options(p.options, {"at": 2, "every": 2, "after": 2, "trim": 4, "strata": 3,
                                  "riskset": 3, "nopreserve": 3})
    newvar = names[0] if names else ""
    at_txt = str(o["at"]).strip() if isinstance(o.get("at"), str) else ""
    failures = at_txt == "failures"
    if failures:
        if newvar:
            raise StataError(101, "varlist not allowed")
        if o.get("every") or o.get("after") or o.get("trim"):
            raise StataError(198, "at(failures) may not be combined with every(), after() or "
                             "trim")   # VERIFY
    else:
        if o.get("strata") or o.get("riskset"):
            raise StataError(198, "options strata() and riskset() require at(failures)")   # VERIFY
        if o.get("every") and not newvar:
            raise StataError(100, "varlist required")
        if bool(at_txt) == bool(o.get("every")):
            raise StataError(198, "must specify either at() or every()")   # VERIFY
    n = ds.nobs
    sel = (ds.get("_st").data == 1) & touse(s, p)
    rows = np.flatnonzero(sel)
    _st, t0_all, t_all, d_all = _st_arrays(s)
    org = _origin(s, info, n)
    # -- cut points of each selected record
    at = np.array([])
    rs_info = None
    if failures:
        svars = s.expand_varlist(str(o["strata"])) if isinstance(o.get("strata"), str) else []
        sg = group_ids(ds, svars, sel)[0] if svars else np.where(sel, 0, -1)
        ft = {int(h): np.unique(t_all[sel & (sg == h) & (d_all == 1)]) for h in np.unique(sg[sel])}
        ordered = sorted(ft)
        offs, tot = {}, 0
        for h in ordered:
            offs[h] = tot
            tot += len(ft[h])
        s.output.write(f"({tot} failure times)\n", "text")
        cuts = [ft[int(sg[i])] for i in rows]
        rs_info = (sg, ft, offs)
    else:
        if o.get("every"):
            step = float(s.eval(str(o["every"])))
            if not step > 0:
                raise StataError(198, "every() must be positive")   # VERIFY
            hi = float(t_all[rows].max()) if len(rows) else 0.0
            lo = float((_after_ref(s, info, str(o["after"]), rows, t_all, org)[rows]).min()) \
                if o.get("after") and len(rows) else 0.0
            at = np.arange(1, int(np.ceil((hi - lo) / step)) + 2) * step
        else:
            big = float(t_all[rows].max()) * 2 + 1e6 if len(rows) else 1e6
            at = np.array(parse_numlist(re.sub(r"\bmax\b", repr(big), at_txt)), dtype=np.float64)
            at = np.unique(at)
        ref = np.zeros(n)
        if o.get("after"):
            ref = _after_ref(s, info, str(o["after"]), rows, t_all, org)
        cuts = [ref[i] + at for i in rows]
    parent, a, b, last = SV.split_pieces(t0_all[rows], t_all[rows], cuts)
    parent_row = rows[parent]
    reps = np.ones(n, dtype=np.int64)
    reps[rows] = np.bincount(parent, minlength=len(rows))
    created = int(len(parent) - len(rows))
    # -- the new values of the pieces
    lab = None if failures or not newvar else SV.category_labels(b, cuts, parent, at)
    excl = np.zeros(len(parent), dtype=bool)
    if o.get("trim") and not failures:
        lo_c = np.array([c.min() for c in cuts])[parent]
        hi_c = np.array([c.max() for c in cuts])[parent]
        excl = (b <= lo_c) | (a >= hi_c)
    rsid = None
    if failures and o.get("riskset"):
        sg, ft, offs = rs_info
        rsid = np.full(len(parent), SYS)
        for j in range(len(parent)):
            h = int(sg[parent_row[j]])
            f = ft[h]
            k = np.searchsorted(f, b[j])
            if k < len(f) and f[k] == b[j]:
                rsid[j] = offs[h] + k + 1
    tv = ds.get(info.timevar).data.astype(np.float64)
    org_rec = tv - t_all * info.scale       # origin of each record, from the stored times
    src = np.repeat(np.arange(n), reps)
    pos = np.flatnonzero(np.repeat(sel, reps))
    fail_new = np.where(last, d_all[parent_row], 0.0)
    _take(ds, src)
    _put(ds, "_t", b, pos)
    _put(ds, "_t0", a, pos)
    _put(ds, "_d", fail_new, pos)
    nonlast = pos[~last]
    ot = org_rec[parent_row]
    _put(ds, info.timevar, (ot + b * info.scale)[~last], nonlast)
    if info.time0:
        first = a == t0_all[parent_row]
        _put(ds, info.time0, (ot + a * info.scale)[~first], pos[~first])
    if info.failvar and ds.has(info.failvar) and not last.all():
        _put(ds, info.failvar, np.zeros(int((~last).sum())), nonlast)   # VERIFY: censored value 0
    if excl.any():
        for nm, val in (("_st", 0.0), ("_d", SYS), ("_t", SYS), ("_t0", SYS)):
            _put(ds, nm, np.full(int(excl.sum()), val), pos[excl])
    if newvar:
        full = np.full(ds.nobs, SYS)
        full[pos] = lab
        ds.add(Variable(newvar, "float", full.astype(np.float32).astype(np.float64)))
    if rsid is not None:
        full = np.full(ds.nobs, SYS)
        full[pos] = rsid
        ds.add(Variable(str(o["riskset"]).strip(), "float", full))
    s.output.write(f"({created} observations (episodes) created)\n", "text")
    s.notify_state()


# ---------------------------------------------------------------------------
# stjoin
# ---------------------------------------------------------------------------

def _same_as_previous(ds: "Dataset", rows: np.ndarray, skip: set[str]) -> np.ndarray:
    """For rows in order: do all the variables not in `skip` have the values of the previous row?"""
    same = np.ones(max(len(rows) - 1, 0), dtype=bool)
    for v in ds.vars:
        if v.name in skip:
            continue
        x = v.raw[rows]
        same &= (x[1:] == x[:-1]) if not v.is_string else np.array([a == b for a, b in zip(x[1:], x[:-1])],
                                                                    dtype=bool)
    return same


def cmd_stjoin(s: "Session", args: str) -> None:
    ds = s.data
    info = require_st(s)
    _need_id(info, "stjoin")
    o = match_options(args.strip().lstrip(","), {"censored": 4}) if args.strip() else {}
    cens = None
    if isinstance(o.get("censored"), str):
        cens = np.array(parse_numlist(str(o["censored"])), dtype=np.float64)
    elif info.failvals is not None and info.failvar:
        cens = np.array([0.0])
    st_, t0, t, d = _st_arrays(s)
    rows, gid, _bounds = _subjects(s, info, np.flatnonzero(st_))
    censored = d == 0
    if cens is not None and info.failvar and ds.has(info.failvar):
        censored = censored & np.isin(ds.get(info.failvar).data.astype(np.float64), cens)
    skip = set(_ST_VARS) | {info.timevar, info.time0, info.failvar}
    skip.discard("")
    gid_all = np.zeros(ds.nobs, dtype=np.int64)
    gid_all[rows] = gid
    same = _same_as_previous(ds, rows, skip)
    join = SV.join_chains(rows, gid_all, t0, t, censored, same)
    drop = np.zeros(ds.nobs, dtype=bool)
    if join.any():
        # each chain keeps its first record, which takes the end of the last one
        head = np.arange(len(rows))
        for k in range(1, len(rows)):
            if join[k]:
                head[k] = head[k - 1]
        for k in np.flatnonzero(join):
            drop[rows[k]] = True
        for k in range(len(rows)):
            h = head[k]
            if h != k and (k + 1 == len(rows) or head[k + 1] != h):
                dst, srcr = rows[h], rows[k]
                for nm in ("_t", "_d", info.timevar, info.failvar):
                    if nm and ds.has(nm):
                        ds.get(nm).raw[dst] = ds.get(nm).raw[srcr]
    removed = ds.keep_obs(~drop)
    s.output.write(f"({removed} obs. eliminated)\n", "text")
    if not removed and not isinstance(o.get("censored"), str):
        s.output.write("(perhaps you need to specify stjoin's censored() option)\n", "text")
    s.notify_state()


# ---------------------------------------------------------------------------
# stfill and stvary
# ---------------------------------------------------------------------------

def _value_vectors(ds: "Dataset", name: str):
    v = ds.get(name)
    if v.is_string:
        x = np.array(v.raw, dtype=object)
        return x, np.array([a == "" for a in x], dtype=bool), True
    x = v.data.astype(np.float64)
    return x, x >= SYS, False


def cmd_stfill(s: "Session", args: str) -> None:
    ds = s.data
    info = require_st(s)
    p = parse_standard(args)
    o = match_options(p.options, {"baseline": 1, "forward": 1, "noshow": 6})
    if bool(o.get("baseline")) == bool(o.get("forward")):
        raise StataError(198, "must specify either baseline or forward")   # VERIFY
    names = s.expand_varlist(p.varlist)
    if not names:
        raise StataError(100, "varlist required")
    allowed = touse(s, p)
    rows, _gid, bounds = _subjects(s, info, np.flatnonzero(ds.get("_st").data == 1))
    if not o.get("noshow"):
        header(s, info)
    out = s.output
    out.write("replace all values with value at earliest observed:\n" if o.get("baseline")
              else "replace missing values with previously observed values:\n", "text")
    for nm in names:
        x, miss, is_str = _value_vectors(ds, nm)
        xs, ms = x[rows], miss[rows]
        new = xs.copy()
        if o.get("baseline"):
            for a, b in zip(bounds[:-1], bounds[1:]):
                new[a:b] = xs[a]
        else:
            src = SV.fill_forward(xs, ms, bounds)
            fill = ms & (src >= 0)
            new[fill] = xs[src[fill]]
        full = x.copy()
        ok = allowed[rows]
        full[rows[ok]] = new[ok]
        changed = int(np.sum([a != b for a, b in zip(full, x)])) if is_str else int(np.sum(full != x))
        if changed:
            (ds.set_string if is_str else ds.set_numeric)(ds.get(nm), full)
        out.write(f"{nm:>16}:  {changed:,} real changes made\n", "text")   # VERIFY: layout
    s.notify_state()


def cmd_stvary(s: "Session", args: str) -> None:
    ds = s.data
    info = require_st(s)
    p = parse_standard(args)
    o = match_options(p.options, {"noshow": 6})
    mask = touse(s, p)
    if p.varlist.strip():
        names = s.expand_varlist(p.varlist)
    else:
        skip = set(_ST_VARS) | {info.idvar, info.failvar}   # VERIFY: variables left out without a varlist
        names = [v.name for v in ds.vars if v.name not in skip]
    rows, _gid, bounds = _subjects(s, info, np.flatnonzero((ds.get("_st").data == 1) & mask))
    if not o.get("noshow"):
        header(s, info)
    out = s.output
    out.write(f"{'':15}subjects for whom the variable is\n", "text")
    out.write(f"{'':49}never    always   sometimes\n", "text")
    out.write(f"{'variable':>12} |  constant    varying{'':13}missing   missing   missing\n", "text")
    out.write("-" * 13 + "+" + "-" * 62 + "\n", "text")
    last = (0, 0, 0, 0, 0)
    for nm in names:
        x, miss, _ = _value_vectors(ds, nm)
        last = SV.subject_summary(x[rows], miss[rows], bounds)
        out.write(f"{nm[:12]:>12} |", "text")
        out.write("".join(f"{v:>{w}}" for v, w in zip(last, (10, 11, 20, 10, 10))) + "\n", "result")
    s.r = dict(zip(("cons", "varies", "never", "always", "miss"), (float(v) for v in last)))


# ---------------------------------------------------------------------------
# stgen
# ---------------------------------------------------------------------------

_TYPES = ("byte", "int", "long", "float", "double")
_NOARG = ("nfailures", "ngaps", "gaplen", "hasgap")
_FUNCS = ("ever", "never", "always", "min", "max", "when", "when0", "count", "count0", "minage", "maxage",
          "avgage") + _NOARG


def _time_bounds(s: "Session", info: STInfo, rows: np.ndarray):
    """Start and end of each record, in the units of the time variable."""
    ds = s.data
    t = ds.get("_t").data.astype(np.float64)
    t0 = ds.get("_t0").data.astype(np.float64)
    end = ds.get(info.timevar).data.astype(np.float64)
    org = end - t * info.scale
    start = ds.get(info.time0).data.astype(np.float64) if info.time0 else org + t0 * info.scale
    return start[rows], end[rows]


def cmd_stgen(s: "Session", args: str) -> None:
    ds = s.data
    info = require_st(s)
    left, eq, right = args.partition("=")
    toks = left.split()
    if not eq or not toks:
        raise StataError(198, "invalid syntax")
    vtype = "float"
    if toks[0] in _TYPES:
        vtype = toks.pop(0)
    if len(toks) != 1:
        raise StataError(198, "invalid syntax")
    newvar = toks[0]
    m = re.match(r"^\s*([A-Za-z0-9_]+)\s*\((.*)\)\s*$", right.strip(), re.S)
    if not m or m.group(1) not in _FUNCS:
        raise StataError(198, f"{right.strip()} invalid function")   # VERIFY
    fn, arg = m.group(1), m.group(2).strip()
    if fn in _NOARG and arg:
        raise StataError(198, f"{fn}() takes no arguments")   # VERIFY
    if fn not in _NOARG and not arg:
        raise StataError(198, f"{fn}() requires an expression")   # VERIFY
    if ds.has(newvar):
        raise StataError(110, f"variable {newvar} already defined")
    rows, _gid, bounds = _subjects(s, info, np.flatnonzero(ds.get("_st").data == 1))
    _st, t0_all, t_all, d_all = _st_arrays(s)
    n = ds.nobs
    res_s = np.full(len(rows), SYS)
    if fn not in _NOARG:
        v = np.broadcast_to(np.asarray(eval_vector(s, arg), dtype=np.float64), (n,))[rows]
        true = (v != 0) & (v < SYS)
    ng, gl = SV.gaps(t0_all[rows], t_all[rows], bounds)
    start, end = _time_bounds(s, info, rows)
    for a, b in zip(bounds[:-1], bounds[1:]):
        sl = slice(a, b)
        if fn == "ever":
            res_s[sl] = float(true[sl].any())
        elif fn == "never":
            res_s[sl] = float(not true[sl].any())
        elif fn == "always":
            res_s[sl] = float(true[sl].all())
        elif fn in ("min", "max"):
            good = v[sl][v[sl] < SYS]
            res_s[sl] = (good.min() if fn == "min" else good.max()) if len(good) else SYS
        elif fn in ("when", "when0"):
            hit = np.flatnonzero(true[sl])
            res_s[sl] = (end if fn == "when" else start)[sl][hit[0]] if len(hit) else SYS
        elif fn in ("count", "count0"):
            res_s[sl] = SV.cumulative_before(true[sl], np.array([0, b - a]), fn == "count0")
        elif fn in ("minage", "maxage", "avgage"):
            ref = {"minage": start, "maxage": end, "avgage": (start + end) / 2}[fn][sl]
            res_s[sl] = np.where(v[sl] < SYS, ref - v[sl], SYS)
        elif fn == "nfailures":
            f = (d_all[rows][sl] == 1).astype(float)
            res_s[sl] = np.cumsum(f) - f
        elif fn == "ngaps":
            res_s[sl] = ng[sl]
        elif fn == "gaplen":
            res_s[sl] = gl[sl]
        else:   # hasgap
            res_s[sl] = float(ng[b - 1] > 0)
    full = np.full(n, SYS)
    full[rows] = res_s
    ds.add(Variable(newvar, vtype, full))
    s.notify_state()


# ---------------------------------------------------------------------------
# snapspan
# ---------------------------------------------------------------------------

def cmd_snapspan(s: "Session", args: str) -> None:
    ds = s.data
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) < 3:
        raise StataError(100 if not names else 102, "varlist required" if not names
                         else "too few variables specified")
    o = match_options(p.options, {"generate": 1, "replace": 1})
    if not isinstance(o.get("generate"), str):
        raise StataError(198, "option generate() required")
    idv, tv = names[0], names[1]
    idv, tv = s.expand_varlist(idv)[0], s.expand_varlist(tv)[0]
    events = s.expand_varlist(" ".join(names[2:]))
    if ds.get(tv).is_string:
        raise StataError(109, "type mismatch")
    if ds.changed and not o.get("replace"):
        raise StataError(4, "no; data in memory would be lost")
    new0 = str(o["generate"]).strip()
    if ds.has(new0):
        raise StataError(110, f"variable {new0} already defined")
    gid = group_ids(ds, [idv], None)[0]
    single = np.bincount(gid)[gid] == 1
    if single.any():
        s.output.write(f"note:  {int(single.sum())} obs. ({100 * single.mean():.1f}%) have only a single "
                       "record; they will be ignored\n", "text")
        if single.all():
            raise StataError(2000, "no observations")
        ds.keep_obs(~single)
        gid = group_ids(ds, [idv], None)[0]
    tt = ds.get(tv).data.astype(np.float64)
    order = np.lexsort((tt, gid))
    _take(ds, order)
    gid, tt = gid[order], tt[order]
    prev = np.zeros(len(order), dtype=bool)
    prev[1:] = gid[1:] == gid[:-1]
    shift_src = np.where(prev, np.arange(len(order)) - 1, -1)
    keep = {idv, tv, *events}
    for v in ds.vars:
        if v.name in keep:
            continue
        old = v.raw.copy()
        if v.is_string:
            v.raw = np.array([old[j] if j >= 0 else "" for j in shift_src], dtype=object)
        else:
            from ..core import storage
            v.raw = np.where(shift_src >= 0, old[np.maximum(shift_src, 0)], storage.missing_raw(v.vtype)) \
                .astype(old.dtype)
    t0 = np.where(prev, tt[np.maximum(shift_src, 0)], SYS)
    tvar = ds.get(tv)
    ds.add(Variable(new0, tvar.vtype, t0, fmt=tvar.fmt), position=ds.index(tv))
    ds.sortlist = [idv, tv]
    s.notify_state()


# ---------------------------------------------------------------------------
# stbase
# ---------------------------------------------------------------------------

def cmd_stbase(s: "Session", args: str) -> None:
    ds = s.data
    info = require_st(s)
    p = parse_standard(args)
    o = match_options(p.options, {"at": 2, "gap": 3, "replace": 1, "noshow": 6, "nopreserve": 3})
    if ds.changed and not o.get("replace"):
        raise StataError(4, "no; data in memory would be lost")
    if o.get("gap") and not isinstance(o.get("at"), str):
        raise StataError(198, "option gap() requires at()")   # VERIFY
    mask = touse(s, p)
    rows, gid, bounds = _subjects(s, info, np.flatnonzero((ds.get("_st").data == 1) & mask))
    if not o.get("noshow"):
        header(s, info)
    st_, t0_all, t_all, d_all = _st_arrays(s)
    if isinstance(o.get("at"), str):
        _stbase_at(s, info, o, rows, bounds, float(s.eval(str(o["at"]))))
        return
    if np.all(np.diff(bounds) == 1):
        s.output.write(f"nothing to do; data has one record per id({info.idvar}) and so is already "
                       "baseline\n", "text")
        return
    # -- baseline values, then the records that can be joined
    skip = set(_ST_VARS) | {info.timevar, info.time0, info.failvar, info.idvar}
    skip.discard("")
    for v in ds.vars:
        if v.name in skip:
            continue
        x = v.raw[rows].copy()
        for a, b in zip(bounds[:-1], bounds[1:]):
            x[a:b] = x[a]
        v.raw[rows] = x
    ng, _gl = SV.gaps(t0_all[rows], t_all[rows], bounds)
    after_fail = False
    for a, b in zip(bounds[:-1], bounds[1:]):
        d = d_all[rows][a:b]
        after_fail |= bool(d[:-1].any()) if b - a > 1 else False
    gid_all = np.zeros(ds.nobs, dtype=np.int64)
    gid_all[rows] = gid
    same = np.ones(max(len(rows) - 1, 0), dtype=bool)
    join = SV.join_chains(rows, gid_all, t0_all, t_all, d_all == 0, same)
    drop = np.zeros(ds.nobs, dtype=bool)
    head = np.arange(len(rows))
    for k in range(1, len(rows)):
        if join[k]:
            head[k] = head[k - 1]
            drop[rows[k]] = True
    for k in range(len(rows)):
        h = head[k]
        if h != k and (k + 1 == len(rows) or head[k + 1] != h):
            for nm in ("_t", "_d", info.timevar, info.failvar):
                if nm and ds.has(nm):
                    ds.get(nm).raw[rows[h]] = ds.get(nm).raw[rows[k]]
    ds.keep_obs(~drop)
    left = SV.group_bounds(gid[~drop[rows]])
    multi = bool(np.any(np.diff(left) > 1))
    out = s.output
    out.write("notes:\n", "text")
    k = 1
    out.write(f"{k}. {'there were gaps' if ng.any() else 'no gaps'}\n", "text")
    if after_fail:
        k += 1
        out.write(f"{k}. there were multiple failures or reentries after failures\n", "text")   # VERIFY
    if multi:
        k += 1
        out.write(f"{k}. baseline data has multiple records per id({info.idvar})\n", "text")   # VERIFY
    k += 1
    out.write(f"{k}. all records have covariate values at baseline\n", "text")
    s.notify_state()


def _stbase_at(s: "Session", info: STInfo, o: dict, rows: np.ndarray, bounds: np.ndarray, at: float) -> None:
    """stbase, at(#): one record per subject, as a cross-section at analysis time #."""
    ds = s.data
    _st, t0_all, t_all, d_all = _st_arrays(s)
    keep_rows, time0, gapv, tend, nfail = [], [], [], [], []
    ng, gl = SV.gaps(t0_all[rows], t_all[rows], bounds)
    for a, b in zip(bounds[:-1], bounds[1:]):
        entry = t0_all[rows][a:b].min()
        exits = t_all[rows][a:b]
        d = d_all[rows][a:b]
        fails = exits[(d == 1) & (exits <= at)]
        if entry >= at:
            continue
        if not len(fails) and exits.max() < at:
            continue
        keep_rows.append(rows[a])
        time0.append(entry)
        tend.append(min(exits.max(), at))
        nfail.append(float(len(fails)))
        gl_at = 0.0
        for k in range(a + 1, b):
            if t0_all[rows][k] > t_all[rows][k - 1] and t0_all[rows][k] < at:
                gl_at += min(t0_all[rows][k], at) - t_all[rows][k - 1]
        gapv.append(gl_at)
    keep = np.zeros(ds.nobs, dtype=bool)
    keep[np.array(keep_rows, dtype=int)] = True
    # names of the new variables
    t0name = info.time0 or next(nm for nm in ("t0", "time0") if not ds.has(nm) or nm == info.time0)
    gapname = str(o["gap"]).strip() if isinstance(o.get("gap"), str) else \
        next(nm for nm in ("gap", "gaptime") if not ds.has(nm))
    failname = info.failvar or next(nm for nm in ("failure", "d") if not ds.has(nm))
    tname = info.timevar
    vals = {t0name: np.array(time0), gapname: np.array(gapv), tname: np.array(tend), failname: np.array(nfail)}
    order = np.array(keep_rows, dtype=int)
    drop_names = [v.name for v in ds.vars if v.name in _ST_VARS]
    _take(ds, order)
    ds.drop_vars(drop_names)
    for nm, arr in vals.items():
        if ds.has(nm):
            ds.drop_vars([nm])
        ds.add(Variable(nm, "double" if nm == tname else "float", arr))
    ch = ds.chars.setdefault("_dta", {})
    for k in [k for k in ch if k.startswith("st_")] + ["_dta"]:
        ch.pop(k, None)
    s.output.write(f"data now cross-section at time {at:g}\n", "text")
    s.run_command(f"summarize {t0name} {gapname} {tname} {failname}")
    s.notify_state()
