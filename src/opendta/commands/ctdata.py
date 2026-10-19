"""Count-time data and case-control sampling of survival data: ct/ctset, cttost, sttoct and sttocc
([ST] ctset, [ST] cttost, [ST] sttoct, [ST] sttocc).

Details that the manuals do not show carry the VERIFY mark (compat/do/1901_survival.do)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import rng
from ..core.dataset import Dataset, Variable, smallest_type_for
from ..core.errors import StataError
from ..core.grouping import group_ids
from ..estimation.st import _stset, header, require_st, st_data
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .strates import _result_dataset

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_KEYS = ("ct_t", "ct_nf", "ct_nc", "ct_ne", "ct_by", "ct_show")


def _ct(s: "Session") -> dict:
    ch = s.data.chars.get("_dta", {})
    if "ct_t" not in ch:
        raise StataError(119, "data not ct")
    return {k: ch.get(k, "") for k in _KEYS}


def _show(s: "Session", c: dict) -> None:
    ds = s.data
    out = s.output
    name = ds.fullpath or ds.filename
    if name:
        out.write(f"dataset name:  {name}\n", "text")
    out.write(f"{'time:':>17}  {c['ct_t']}\n", "text")
    out.write(f"{'no. fail:':>17}  {c['ct_nf']}\n", "text")
    out.write(f"{'no. lost:':>17}  {c['ct_nc'] or '-- (meaning 0 lost)'}\n", "text")
    out.write(f"{'no. enter:':>17}  {c['ct_ne'] or '-- (meaning all enter at time 0)'}\n", "text")
    if c["ct_by"]:
        out.write(f"{'by:':>17}  {c['ct_by']}\n", "text")   # VERIFY: alignment of the labels


def _check(s: "Session", c: dict) -> None:
    ds = s.data
    t = ds.get(c["ct_t"]).data
    if np.any(t[t < SYS] < 0):
        raise StataError(459, f"{c['ct_t']} has negative values")   # VERIFY
    for key in ("ct_nf", "ct_nc", "ct_ne"):
        if c[key]:
            x = ds.get(c[key]).data
            if np.any(x >= SYS) or np.any(x < 0):
                raise StataError(459, f"{c[key]} has missing or negative values")   # VERIFY
    if c["ct_ne"]:
        by = c["ct_by"].split()
        gid, k = group_ids(ds, by)
        out = ds.get(c["ct_nf"]).data + (ds.get(c["ct_nc"]).data if c["ct_nc"] else 0)
        ent = ds.get(c["ct_ne"]).data
        for g in range(k):
            if ent[gid == g].sum() != out[gid == g].sum():
                raise StataError(459, "entering and exiting subjects do not balance")   # VERIFY


def cmd_ctset(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "show": 4, "noshow": 6, "clear": 5})
    ds = s.data
    ch = ds.chars.setdefault("_dta", {})
    if o.get("clear"):
        for k in _KEYS:
            ch.pop(k, None)
        return
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if not names:
        if o.get("show") or o.get("noshow"):
            _ct(s)
            ch["ct_show"] = "1" if o.get("show") else "0"
            return
        c = _ct(s)
        _show(s, c)
        _check(s, c)
        return
    if len(names) > 4:
        raise StataError(103, "too many variables specified")
    if len(names) < 2:
        raise StataError(102, "too few variables specified")
    for n in names:
        if ds.get(n).is_string:
            raise StataError(109, "type mismatch")
    by = s.expand_varlist(str(o["by"])) if o.get("by") else []
    c = {"ct_t": names[0], "ct_nf": names[1], "ct_nc": names[2] if len(names) > 2 else "",
         "ct_ne": names[3] if len(names) > 3 else "", "ct_by": " ".join(by),
         "ct_show": "0" if o.get("noshow") else "1"}
    _check(s, c)
    for k in _KEYS:
        ch.pop(k, None)
    ch.update(c)
    _show(s, c)


def cmd_ct(s: "Session", args: str) -> None:
    if args.strip():
        raise StataError(198, "invalid syntax")
    cmd_ctset(s, "")


def _free_name(ds: Dataset, given, candidates) -> str:
    if given:
        if ds.has(str(given)):
            raise StataError(110, f"{given} already defined")
        return str(given)
    for nm in candidates:
        if not ds.has(nm):
            return nm
    raise StataError(110, "variable already defined; specify wvar() and t0()")   # VERIFY


# ---------------------------------------------------------------------------
# cttost
# ---------------------------------------------------------------------------

def cmd_cttost(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"t0": 2, "wvar": 1, "clear": 5, "nopreserve": 5})
    c = _ct(s)
    ds = s.data
    if ds.changed and not (o.get("clear") or o.get("nopreserve")):
        raise StataError(4, "no; data in memory would be lost")
    _check(s, c)
    wname = _free_name(ds, o.get("wvar"), ("w", "pop", "weight", "wgt"))
    t0name = _free_name(ds, o.get("t0"), ("t0", "time0", "etime")) if c["ct_ne"] else ""
    tv, nf = ds.get(c["ct_t"]), ds.get(c["ct_nf"])
    n = ds.nobs
    nc = ds.get(c["ct_nc"]).data if c["ct_nc"] else np.zeros(n)
    t, f = tv.data, nf.data
    by = c["ct_by"].split()
    src, tt, dd, ww, ee = [], [], [], [], []
    if not c["ct_ne"]:
        # failures first and then censorings at each time, everyone entering at 0
        for i in range(n):
            for d, cnt in ((1, f[i]), (0, nc[i])):
                if cnt > 0:
                    src.append(i)
                    tt.append(t[i])
                    dd.append(d)
                    ww.append(cnt)
    else:
        ne = ds.get(c["ct_ne"]).data
        gid, k = group_ids(ds, by)
        # each exit draws its entry time at random among the entries before it (VERIFY: Stata
        # does not document the draw)
        for g in range(k):
            idx = np.flatnonzero(gid == g)
            et = np.unique(t[idx][ne[idx] > 0])
            remain = np.array([ne[idx][t[idx] == x].sum() for x in et], dtype=float)
            units = []
            for i in idx[np.argsort(t[idx], kind="stable")]:
                units += [(i, 1)] * int(f[i]) + [(i, 0)] * int(nc[i])
            u = rng.RNG.uniform(len(units)) if units else []
            tally: dict = {}
            for (i, d), ui in zip(units, u):
                avail = remain[:np.searchsorted(et, t[i], "left")]
                if avail.sum() <= 0:
                    raise StataError(459, "entering and exiting subjects cannot be matched")   # VERIFY
                j = int(np.searchsorted(np.cumsum(avail), np.floor(ui * avail.sum()) + 1, "left"))
                remain[j] -= 1
                tally[(i, d, j)] = tally.get((i, d, j), 0) + 1
            for (i, d, j), cnt in tally.items():
                src.append(i)
                tt.append(t[i])
                dd.append(d)
                ww.append(cnt)
                ee.append(et[j])
    cols = [(c["ct_nf"], dd, "%9.0g"), (wname, ww, "%9.0g")]
    if t0name:
        cols.append((t0name, ee, tv.fmt))
    new = _result_dataset(s, by + [c["ct_t"]], src, [])
    for nm, vals, f_ in cols:
        arr = np.array(vals, dtype=np.float64)
        vt = "byte" if nm == c["ct_nf"] else ("double" if nm == t0name else smallest_type_for(arr))
        new.vars.append(Variable(nm, vt, arr, fmt=f_ if vt in ("double", "float") or nm != c["ct_nf"] else "%8.0g"))
    new.nobs = len(src)
    s.data = new
    s.notify_state()
    text = f"{c['ct_t']} [fweight={wname}], failure({c['ct_nf']})"
    if t0name:
        text += f" time0({t0name})"
    _stset(s, text, report=True)


# ---------------------------------------------------------------------------
# sttoct
# ---------------------------------------------------------------------------

def cmd_sttoct(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "replace": 1, "noshow": 6})
    info = require_st(s)
    names = p.varlist.split()
    if len(names) not in (2, 3):
        raise StataError(198, "invalid syntax")
    ds = s.data
    for nm in names:
        if ds.has(nm):
            raise StataError(110, f"{nm} already defined")
    by = s.expand_varlist(str(o["by"])) if o.get("by") else []
    sd = st_data(s, touse(s, p))
    rows = np.flatnonzero(sd.mask)
    t, t0, d, w = sd.t, sd.t0, sd.d, sd.w
    gid, k = group_ids(ds, by, sd.mask)
    gid = gid[rows]
    if len(names) == 2:
        keep = _subject_ends(sd, gid, by, ds, rows, info)
    else:
        keep = np.ones(len(t), dtype=bool)
    times = np.unique(np.concatenate([t[keep], t0[keep]])) if len(names) == 3 else np.unique(t[keep])
    out_rows = []
    for g in range(k):
        sel = gid == g
        first = int(rows[np.flatnonzero(sel)[0]])
        for x in (times if len(names) == 3 else np.unique(t[keep & sel])):
            nfail = float(w[keep & sel & (t == x) & (d == 1)].sum())
            ncens = float(w[keep & sel & (t == x) & (d == 0)].sum())
            nent = float(w[sel & (t0 == x)].sum()) if len(names) == 3 else 0.0
            if nfail or ncens or nent:
                out_rows.append((first, x, nfail, ncens, nent))
    firsts = [r[0] for r in out_rows]
    cols = [("t", [r[1] for r in out_rows], "%9.0g")]
    for j, nm in enumerate(names):
        cols.append((nm, [r[2 + j] for r in out_rows], "%9.0g"))
    new = _result_dataset(s, by, firsts, [])
    for nm, vals, f_ in cols:
        arr = np.array(vals, dtype=np.float64)
        new.vars.append(Variable(nm, "double" if nm == "t" else smallest_type_for(arr), arr, fmt=f_))
    new.nobs = len(firsts)
    s.data = new
    s.notify_state()
    ch = new.chars.setdefault("_dta", {})
    ch.update({"ct_t": "t", "ct_nf": names[0], "ct_nc": names[1], "ct_ne": names[2] if len(names) == 3 else "",
               "ct_by": " ".join(by), "ct_show": "0" if o.get("noshow") else "1"})


def _subject_ends(sd, gid, by, ds, rows, info) -> np.ndarray:
    """The last record of each subject (all records if there is no id()); checks the conditions
    for the form without newentvar: everyone enters at 0, no gaps, no recurring risk, and by()
    constant within subject (VERIFY the messages)."""
    t, t0, d = sd.t, sd.t0, sd.d
    if sd.ids is None:
        if np.any(t0 != 0):
            raise StataError(459, "sttoct requires newentvar when subjects do not all enter at time 0")   # VERIFY
        return np.ones(len(t), dtype=bool)
    keep = np.zeros(len(t), dtype=bool)
    for sid in np.unique(sd.ids):
        j = np.flatnonzero(sd.ids == sid)
        j = j[np.argsort(t[j], kind="stable")]
        if t0[j[0]] != 0 or np.any(t0[j[1:]] != t[j[:-1]]) or np.any(d[j[:-1]] != 0) \
                or len(np.unique(gid[j])) > 1:
            raise StataError(459, "sttoct requires newentvar for these data")   # VERIFY
        keep[j[-1]] = True
    return keep


# ---------------------------------------------------------------------------
# sttocc
# ---------------------------------------------------------------------------

def cmd_sttocc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"match": 1, "number": 1, "nodots": 3, "generate": 1})
    info = require_st(s)
    ds = s.data
    number = int(float(str(o["number"]))) if o.get("number") else 1
    if number < 1:
        raise StataError(198, "number() must be a positive integer")   # VERIFY
    gen = str(o["generate"]).split() if o.get("generate") else ["_case", "_set", "_time"]
    if len(gen) != 3:
        raise StataError(198, "generate() requires three names")
    for nm in gen:
        if ds.has(nm):
            raise StataError(110, f"{nm} already defined")
    mvars = s.expand_varlist(str(o["match"])) if o.get("match") else []
    keep_names = s.expand_varlist(p.varlist) if p.varlist.strip() else \
        [n for n in ds.names if n not in ("_st", "_d", "_t", "_t0")]
    sd = st_data(s, touse(s, p))
    rows = np.flatnonzero(sd.mask)
    t, t0, d = sd.t, sd.t0, sd.d
    ids = sd.ids if sd.ids is not None else np.arange(len(t))
    mg = group_ids(ds, mvars, sd.mask)[0][rows] if mvars else np.zeros(len(t), dtype=int)
    cases = np.flatnonzero(d == 1)
    out = s.output
    header(s, info)
    if mvars:
        out.write(f"{'matching for':>19}:  {' '.join(mvars)}\n", "text")
    ft, cnt = np.unique(t[cases], return_counts=True)
    ties = int((cnt > 1).sum())
    if ties:
        out.write(f"There were {ties} tied times involving failure(s)\n", "text")
        out.write("- tied failure times split at random\n", "text")
    out.write(f"There are {len(cases)} cases\n", "text")
    out.write(f"Sampling {number} controls for each case \n", "text")
    order = np.full(len(t), 2.0)
    order[cases] = rng.RNG.uniform(len(cases)) if len(cases) else []
    picked, setid, ctime, is_case = [], [], [], []
    for k, c in enumerate(cases, start=1):
        tf = t[c]
        elig = (t0 <= tf) & ((t > tf) | ((t == tf) & (order > order[c]))) & (ids != ids[c]) & (mg == mg[c])
        cand = np.flatnonzero(elig)
        if len(cand):
            _, first = np.unique(ids[cand], return_index=True)
            cand = cand[np.sort(first)]
            if len(cand) < number:
                out.write(f" Warning: sample requested greater than population. Only {len(cand)} "
                          "controls selected\n", "text")
            u = rng.RNG.uniform(len(cand))
            cand = cand[np.argsort(u, kind="stable")[:number]]
        for j, flag in [(c, 1)] + [(j, 0) for j in cand]:
            picked.append(rows[j])
            setid.append(k)
            ctime.append(tf)
            is_case.append(flag)
    ordr = np.lexsort((np.array(is_case), np.array(setid)))
    picked = [picked[i] for i in ordr]
    setid, ctime, is_case = ([x[i] for i in ordr] for x in (setid, ctime, is_case))
    new = _result_dataset(s, keep_names, picked,
                          [(gen[0], is_case, "%8.0g"), (gen[1], setid, "%12.0g"), (gen[2], ctime, "%8.0g")])
    for nm, vt, lab in ((gen[0], "byte", "0 for controls; 1 for cases"), (gen[1], "long", "case-control ID"),
                        (gen[2], smallest_type_for(new.get(gen[2]).data), "analysis time of the case's failure")):
        v = new.get(nm)
        v.vtype = vt
        v.label = lab
    new.sortlist = [gen[1], gen[0]]
    new.changed = True
    s.data = new
    s.notify_state()
