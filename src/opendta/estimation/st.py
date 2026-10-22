"""Survival-time data: stset, st, stdescribe and stsum ([ST] stset,
[ST] stdescribe, [ST] stsum).

stset keeps its settings in the _dta characteristics (st_bt, st_bd,
st_id...) and creates _st, _d, _t and _t0, as described in [ST] stset.
The layout of the report follows the manual's examples; details that the
manual does not show carry the VERIFY mark (compat/do/0901_st.do).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_options

if TYPE_CHECKING:
    from ..core.dataset import Dataset
    from ..session import Session

SYS = M.SYSMISS
LINE = "-" * 78


def _g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _count(x: float) -> str:
    """Counts and totals in the stset report: integers without commas, otherwise about six
    significant digits (192.286). VERIFY: the exact format Stata uses for non-integer totals."""
    if float(x).is_integer():
        return f"{int(x)}"
    return format_value(float(x), "%9.0g", pad=False).strip()


def _plural(k: int, one: str, many: str) -> str:
    return one if k == 1 else many


# exclusion reasons: (singular, plural, probable error, continuation line)
_EXCL = {
    "begins_after_fail": ("observation begins on or after (first) failure",
                          "observations begin on or after (first) failure", False, ""),
    "ends_before_enter": ("observation ends on or before enter()",
                          "observations end on or before enter()", False, ""),
    "begins_after_exit": ("observation begins on or after exit",
                          "observations begin on or after exit", False, ""),   # VERIFY
}


# ---------------------------------------------------------------------------
# settings
# ---------------------------------------------------------------------------

@dataclass
class STInfo:
    timevar: str
    time0: str
    failvar: str
    failvals: list[float] | None          # None: any nonzero nonmissing value
    idvar: str
    origin: str                            # "time exp" text (without the word time)
    enter: str
    exit: str                              # "failure", "time ." or "time exp"
    scale: float
    ifexp: str
    wtype: str
    wvar: str
    show: bool
    args: str


def st_info(ds: "Dataset") -> STInfo | None:
    ch = ds.chars.get("_dta", {})
    if ch.get("_dta") != "st" or not all(ds.has(v) for v in ("_st", "_d", "_t", "_t0")):
        return None
    ev = ch.get("st_ev", "")
    return STInfo(ch.get("st_bt", ""), ch.get("st_bt0", ""), ch.get("st_bd", ""),
                  [float(v) for v in ev.split()] if ev else None, ch.get("st_id", ""),
                  ch.get("st_orig", ""), ch.get("st_enter", ""), ch.get("st_exit", ""),
                  float(ch.get("st_bs", "1") or 1), ch.get("st_ifexp", ""), ch.get("st_wt", ""),
                  ch.get("st_wv", ""), ch.get("st_show", "") != "noshow", ch.get("st_args", ""))


def require_st(s: "Session") -> STInfo:
    info = st_info(s.data)
    if info is None:
        if s.data.chars.get("_dta", {}).get("ct_t"):
            raise StataError(119, "data not st (data are marked as being ct data)")
        raise StataError(119, "data not st")
    return info


def _fail_text(info: STInfo) -> str:
    if not info.failvar:
        return ""
    if info.failvals is None:
        return info.failvar
    return f"{info.failvar} == {' '.join(_g(v) for v in info.failvals)}"


def _analysis_text(info: STInfo, word: str = "") -> str:
    """The analysis-time expression; the stset report writes the literal word "time" and the st
    header the variable name (tt/7 against time/7)."""
    t = word or info.timevar
    if info.origin:
        t = f"({t}-origin)"
    if info.scale != 1:
        t = f"{t}/{_g(info.scale)}"
    return t


def header(s: "Session", info: STInfo | None = None, *, force: bool = False) -> None:
    """The lines that st commands print first (unless stset ..., noshow)."""
    info = info or require_st(s)
    if not info.show and not force:
        return
    out = s.output
    out.write("\n", "text")
    ft = _fail_text(info)
    out.write(f"{'failure _d':>19}:  {ft if ft else '1 (meaning all fail)'}\n", "text")
    out.write(f"{'analysis time _t':>19}:  {_analysis_text(info)}\n", "text")
    if info.origin:
        out.write(f"{'origin':>19}:  time {info.origin}\n", "text")
    if info.enter:
        out.write(f"{'enter on or after':>19}:  time {info.enter}\n", "text")
    if info.exit and info.exit != "failure":
        out.write(f"{'exit on or before':>19}:  {info.exit}\n", "text")
    if info.idvar:
        out.write(f"{'id':>19}:  {info.idvar}\n", "text")
    if info.ifexp:
        out.write(f"{'keep if exp':>19}:  {info.ifexp}\n", "text")   # VERIFY
    if info.wtype:
        out.write(f"{'weight':>19}:  [{info.wtype}={info.wvar}]\n", "text")   # VERIFY
    out.write("\n", "text")


# ---------------------------------------------------------------------------
# stset
# ---------------------------------------------------------------------------

_ST_VARS = ("_st", "_d", "_t", "_t0")
_ST_CHARS = ("_dta", "st_ver", "st_id", "st_bt", "st_bt0", "st_bd", "st_ev", "st_enter", "st_exit",
             "st_orig", "st_bs", "st_ifexp", "st_wt", "st_wv", "st_show", "st_args", "st_s")


def _clear(ds: "Dataset") -> None:
    ch = ds.chars.setdefault("_dta", {})
    if ch.get("_dta") == "st":
        for k in _ST_CHARS:
            ch.pop(k, None)
    have = [v for v in _ST_VARS if ds.has(v)]
    if have:
        ds.drop_vars(have)


def _time_exp(text: str, what: str) -> str:
    t = text.strip()
    if t.startswith("time "):
        return t[5:].strip()
    if t == "min" and what == "origin":
        return "min"
    if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", t):
        return t
    raise StataError(198, f"{what}() invalid")   # VERIFY


def cmd_stset(s: "Session", args: str) -> None:
    ds = s.data
    text = args.strip()
    if text.startswith(","):
        opts = match_options(text[1:], {"clear": 5, "noshow": 6, "show": 4})
        if opts.get("clear"):
            _clear(ds)
            s.notify_state()
            return
        info = require_st(s)
        ch = ds.chars["_dta"]
        if opts.get("noshow"):
            ch["st_show"] = "noshow"
        elif opts.get("show"):
            ch["st_show"] = ""
        return
    if not text:
        info = require_st(s)
        _stset(s, info.args, report=True)
        return
    _stset(s, text, report=True)


def _stset(s: "Session", text: str, *, report: bool) -> None:
    from ..commands._util import eval_vector
    from ..core.varlist import resolve_name
    from ..lang.syntax import parse_standard
    ds = s.data
    p = parse_standard(text)
    spec = {"failure": 1, "id": 2, "origin": 2, "enter": 2, "exit": 2, "scale": 2, "time0": 5,
            "if": 2, "noshow": 6, "show": 4, "ever": 4, "never": 5, "after": 5, "before": 6}
    o = {}
    for name, arg in (parse_options(p.options) if p.options.strip() else []):
        nm = "if" if name == "if" else name
        for full, k in spec.items():
            if len(nm) >= k and full.startswith(nm):
                o[full] = arg if arg is not None else True
                break
        else:
            raise StataError(198, f"option {name} not allowed")
    words = p.varlist.split()
    if len(words) != 1:
        raise StataError(102 if not words else 103,
                         "too few variables specified" if not words else "too many variables specified")
    if p.if_ or p.in_:
        raise StataError(101, "if not allowed" if p.if_ else "in range not allowed")   # VERIFY
    tvar = resolve_name(ds, words[0])
    if ds.get(tvar).is_string:
        raise StataError(109, "type mismatch")
    n = ds.nobs
    # remove previous settings first (their variables would be confused with the data)
    old_vars = [v for v in _ST_VARS if ds.has(v)]
    if old_vars:
        ds.drop_vars(old_vars)
    t_raw = ds.get(tvar).data.astype(np.float64)
    failvar, failvals = "", None
    if isinstance(o.get("failure"), str):
        ftxt = o["failure"].strip()
        if "==" in ftxt:
            a, b = ftxt.split("==", 1)
            failvar = resolve_name(ds, a.strip())
            failvals = [float(s.eval(v)) for v in b.split()]
        else:
            failvar = resolve_name(ds, ftxt)
        if ds.get(failvar).is_string:
            raise StataError(109, "type mismatch")
    idvar = resolve_name(ds, str(o["id"]).strip()) if isinstance(o.get("id"), str) else ""
    time0 = resolve_name(ds, str(o["time0"]).strip()) if isinstance(o.get("time0"), str) else ""
    scale = float(s.eval(str(o["scale"]))) if isinstance(o.get("scale"), str) else 1.0
    if not scale > 0:
        raise StataError(198, "scale() must be positive")   # VERIFY
    origin = _time_exp(o["origin"], "origin") if isinstance(o.get("origin"), str) else ""
    enter = _time_exp(o["enter"], "enter") if isinstance(o.get("enter"), str) else ""
    exit_ = ""
    if isinstance(o.get("exit"), str):
        e = o["exit"].strip()
        exit_ = "failure" if "failure".startswith(e) and len(e) >= 1 else \
            ("time " + _time_exp(e, "exit"))
    else:
        exit_ = "failure" if failvar else "failure"
    ifexp = o["if"].strip() if isinstance(o.get("if"), str) else ""
    wtype, wvar = "", ""
    if p.weight is not None:
        wtype = {"fweight": "fweight", "pweight": "pweight", "iweight": "iweight"}.get(p.weight[0], "")
        if not wtype:
            raise StataError(101, f"{p.weight[0]}s not allowed")
        wvar = resolve_name(ds, p.weight[1].strip())

    def num(expr: str) -> np.ndarray:
        v = eval_vector(s, expr)
        return np.broadcast_to(np.asarray(v, dtype=np.float64), (n,)).copy()

    # failure indicator
    if failvar:
        fv = ds.get(failvar).data.astype(np.float64)
        if failvals is None:
            dfail = (fv != 0) & (fv < SYS)
        else:
            dfail = np.isin(fv, failvals)
    else:
        dfail = np.ones(n, dtype=bool)
    # exclusions, in the order they are reported: [count, singular, plural, probable error,
    # continuation line]
    excl: list[list] = []
    keep = np.ones(n, dtype=bool)

    def drop(mask: np.ndarray, why, error: bool = False, cont: str = "") -> None:
        nonlocal keep
        one, many = (why, why) if isinstance(why, str) else why
        m = keep & mask
        if m.any():
            for e in excl:
                if e[1] == one:
                    e[0] += int(m.sum())
                    break
            else:
                excl.append([int(m.sum()), one, many, error, cont])
        keep &= ~mask

    def reason(key: str):
        one, many, _err, _c = _EXCL[key]
        return (one, many)

    drop(t_raw >= SYS, f"event time missing ({tvar}>=.)", True)
    if ifexp:
        iv = num(ifexp)
        drop(~((iv != 0) & (iv < SYS)), "ignored per if()")   # VERIFY
    w = None
    if wvar:
        w = ds.get(wvar).data.astype(np.float64)
        drop(~((w >= 0) & (w < SYS)), "weights invalid")   # VERIFY
    org = np.zeros(n)
    if origin:
        if origin == "min":
            org = np.zeros(n)
        else:
            org = num(origin)
            drop(org >= SYS, "origin missing")   # VERIFY
    t0_raw = np.full(n, np.nan)
    if time0:
        t0_raw = ds.get(time0).data.astype(np.float64)
        drop(t0_raw >= SYS, f"entry time missing ({time0}>=.)", True)   # VERIFY
    ids = ds.get(idvar).data if idvar else None
    if idvar:
        miss_id = (ids == "") if ds.get(idvar).is_string else (ids.astype(np.float64) >= SYS)
        drop(np.asarray(miss_id, dtype=bool), f"obs. with missing {idvar}")   # VERIFY
    # begin of each record, in raw time units
    begin = np.where(np.isnan(t0_raw), org, t0_raw) if time0 else org.copy()
    order = np.arange(n)
    if idvar and not time0:
        # records of a subject are consecutive spells: each begins where the previous ended;
        # Stata excludes both records that end at the same instant (compat 0901)
        key_id = np.array([str(v) for v in ids]) if ds.get(idvar).is_string else ids.astype(np.float64)
        order = np.lexsort((t_raw, key_id))
        prev_i, prev_id, prev_t = None, None, None
        same = np.zeros(n, dtype=bool)
        for i in order:
            if not keep[i]:
                continue
            if prev_id is not None and key_id[i] == prev_id:
                if t_raw[i] == prev_t:
                    same[i] = True
                    same[prev_i] = True
                begin[i] = prev_t
            prev_i, prev_id, prev_t = i, key_id[i], t_raw[i]
        drop(same, "multiple records at same instant", True, f"({tvar}[_n-1]=={tvar})")
    elif idvar and time0:
        key_id = np.array([str(v) for v in ids]) if ds.get(idvar).is_string else ids.astype(np.float64)
        order = np.lexsort((t_raw, key_id))
    else:
        key_id = None
    if time0:
        drop(t0_raw >= t_raw, f"entry on or after exit ({time0}>{tvar})", True)
    if idvar and time0:
        overl = np.zeros(n, dtype=bool)
        prev_id, prev_t = None, None
        for i in order:
            if not keep[i]:
                continue
            if prev_id is not None and key_id[i] == prev_id and begin[i] < prev_t:
                overl[i] = True
            prev_id, prev_t = key_id[i], t_raw[i]
        drop(overl, "overlapping records", True, f"({tvar}[_n-1]>{time0})")   # VERIFY
    # origin and enter() exclusions share one line (compat 0901: 4 + 6 reported as 10)
    if origin:
        drop(t_raw <= org, reason("ends_before_enter"))
    drop(t_raw <= begin, reason("ends_before_enter"))
    if enter:
        ent = num(enter)
        drop(keep & (t_raw <= ent), reason("ends_before_enter"))
        begin = np.where(keep & (ent < SYS) & (ent > begin), ent, begin)
    tend = t_raw.copy()
    d = dfail.copy()
    # exit
    if exit_.startswith("time "):
        etxt = exit_[5:].strip()
        if etxt != ".":
            ex = num(etxt)
            drop(keep & (begin >= ex), reason("begins_after_exit"))
            cut = keep & (tend > ex)
            tend = np.where(cut, ex, tend)
            d = np.where(cut, False, d)
    elif exit_ == "failure" and idvar:
        # after the first failure of a subject, its later records are not used
        fail_t: dict = {}
        for i in order:
            if keep[i] and d[i]:
                k = key_id[i]
                if k not in fail_t:
                    fail_t[k] = t_raw[i]
        after = np.array([keep[i] and key_id[i] in fail_t and begin[i] >= fail_t[key_id[i]]
                          for i in range(n)], dtype=bool)
        drop(after, reason("begins_after_fail"))
    _t = (tend - org) / scale
    _t0 = np.maximum((begin - org) / scale, 0.0)
    # store
    from ..core.dataset import Variable
    is_double = ds.get(tvar).vtype == "double"
    ttype = "double" if is_double else "float"
    st_ = keep.astype(np.float64)
    dd = np.where(keep, d.astype(np.float64), SYS)
    tt = np.where(keep, _t, SYS)
    tt0 = np.where(keep, _t0, SYS)
    if ttype == "float":
        tt[keep] = tt[keep].astype(np.float32)
        tt0[keep] = tt0[keep].astype(np.float32)
    ds.add(Variable("_st", "byte", st_))
    ds.add(Variable("_d", "byte", dd))
    ds.add(Variable("_t", ttype, tt))
    ds.add(Variable("_t0", ttype, tt0))
    ch = ds.chars.setdefault("_dta", {})
    for k in _ST_CHARS:
        ch.pop(k, None)
    shown = "noshow" if o.get("noshow") else ""
    ch.update({"_dta": "st", "st_ver": "2", "st_bt": tvar, "st_bd": failvar,
               "st_ev": " ".join(_g(v) for v in failvals) if failvals else "", "st_id": idvar,
               "st_bt0": time0, "st_orig": origin, "st_enter": enter, "st_exit": exit_,
               "st_bs": _g(scale) if scale != 1 else "", "st_ifexp": ifexp, "st_wt": wtype,
               "st_wv": wvar, "st_show": shown, "st_args": text})
    for k in [k for k, v in ch.items() if k.startswith("st_") and v == ""]:
        ch.pop(k)
    s.notify_state()
    if report:
        _report(s, st_info(ds), keep, excl, n, key_id, dd, tt, tt0, w)


def _report(s: "Session", info: STInfo, keep, excl, n, key_id, dd, tt, tt0, w) -> None:
    out = s.output
    _report_head(out, info)
    out.write("\n" + LINE + "\n", "text")
    _report_body(out, info, keep, excl, n, key_id, dd, tt, tt0, w)


def _report_head(out, info: STInfo) -> None:
    out.write("\n", "text")
    if info.idvar:
        out.write(f"{'id':>18}:  {info.idvar}\n", "text")
    ft = _fail_text(info)
    if ft:
        out.write(f"{'failure event':>18}:  {ft if info.failvals else ft + ' != 0 & ' + ft + ' < .'}\n", "text")
    else:
        out.write(f"{'failure event':>18}:  (assumed to fail at time={info.timevar})\n", "text")
    if info.time0:
        interval = f"({info.time0}, {info.timevar}]"
    elif info.idvar:
        interval = f"({info.timevar}[_n-1], {info.timevar}]"
    else:
        interval = f"(0, {info.timevar}]" if not info.origin else f"(origin, {info.timevar}]"
    out.write(f"{'obs. time interval':>18}:  {interval}\n", "text")
    if info.enter:
        out.write(f"{'enter on or after':>18}:  time {info.enter}\n", "text")
    out.write(f"{'exit on or before':>18}:  {info.exit if info.exit != 'failure' else 'failure'}\n", "text")
    if info.origin or info.scale != 1:
        out.write(f"{'t for analysis':>18}:  {_analysis_text(info, 'time')}\n", "text")
    if info.origin:
        out.write(f"{'origin':>18}:  time {info.origin}\n", "text")
    if info.ifexp:
        out.write(f"{'keep if exp':>18}:  {info.ifexp}\n", "text")   # VERIFY
    if info.wtype:
        out.write(f"{'weight':>18}:  [{info.wtype}={info.wvar}]\n", "text")   # VERIFY


def _report_body(out, info: STInfo, keep, excl, n, key_id, dd, tt, tt0, w) -> None:
    def line(v, txt, style="result"):
        out.write(f"{_count(v):>11}", style)
        out.write(f"  {txt}\n", "text")

    line(n, "total observations")
    if excl:
        for k, one, many, err, cont in excl:
            txt = _plural(k, one, many)
            if err:
                txt = f"{txt:<51}PROBABLE ERROR"
            line(k, txt)
            if cont:
                out.write(f"{'':13}{cont}\n", "text")
    else:
        line(0, "exclusions")
    out.write(LINE + "\n", "text")
    nk = int(keep.sum())
    wk = w[keep] if w is not None and info.wtype == "fweight" else None
    if wk is not None:
        # frequency weights: physical observations first, then their weighted count (cttost manual)
        line(nk, "physical observations remaining, equal to")
        line(float(wk.sum()), "weighted observations, representing")
    else:
        line(nk, "observations remaining, representing")
    if info.idvar:
        line(len(set(key_id[keep].tolist())), "subjects")
    nfail = float(np.sum(dd[keep] * (wk if wk is not None else 1)))
    if info.idvar:
        kind = "single-failure-per-subject" if info.exit == "failure" else "multiple-failure-per-subject"
    else:
        kind = "single-record/single-failure"
    line(nfail, f"failures in {kind} data")
    risk = float(np.sum((tt[keep] - tt0[keep]) * (wk if wk is not None else 1)))
    line(risk, "total analysis time at risk and under observation")
    if nk:
        for lab, val in (("at risk from t =", 0.0), ("earliest observed entry t =", float(tt0[keep].min())),
                         ("last observed exit t =", float(tt[keep].max()))):
            out.write(f"{lab:>64}", "text")
            out.write(f"{_g(val):>10}\n", "result")


def cmd_st(s: "Session", args: str) -> None:
    """st shows the stset command and the first lines of its report ([ST] st_is)."""
    info = require_st(s)
    out = s.output
    out.write(f"-> stset {info.args}\n", "text")
    _report_head(out, info)


def cmd_streset(s: "Session", args: str) -> None:
    info = require_st(s)
    extra = args.strip()
    base = info.args
    if extra:
        if extra.startswith(","):
            base = base + (" " if "," in base else ", ") + extra[1:].strip()
        else:
            raise StataError(198, "invalid syntax")   # VERIFY
    _stset(s, base, report=True)


# ---------------------------------------------------------------------------
# shared data access
# ---------------------------------------------------------------------------

@dataclass
class STData:
    mask: np.ndarray        # observations used (over the whole dataset)
    t: np.ndarray
    t0: np.ndarray
    d: np.ndarray
    w: np.ndarray           # frequency weights (ones without weights)
    ids: np.ndarray | None  # subject identifiers (over mask), or None
    info: STInfo


def st_data(s: "Session", mask: np.ndarray | None = None) -> STData:
    info = require_st(s)
    ds = s.data
    st_ = ds.get("_st").data
    m = (st_ == 1)
    if mask is not None:
        m &= mask
    t = ds.get("_t").data[m].astype(np.float64)
    t0 = ds.get("_t0").data[m].astype(np.float64)
    d = (ds.get("_d").data[m] == 1).astype(np.float64)
    w = np.ones(int(m.sum()))
    if info.wvar and ds.has(info.wvar):
        w = ds.get(info.wvar).data[m].astype(np.float64)
    ids = None
    if info.idvar and ds.has(info.idvar):
        raw = ds.get(info.idvar).data[m]
        ids = np.array([str(v) for v in raw]) if ds.get(info.idvar).is_string else raw.astype(np.float64)
    return STData(m, t, t0, d, w, ids, info)


def n_subjects(sd: STData) -> float:
    if sd.ids is None:
        return float(sd.w.sum()) if sd.info.wtype == "fweight" else float(len(sd.t))
    return float(len(np.unique(sd.ids)))


def _qualifiers(s: "Session", args: str, spec: dict):
    """[varlist] [if] [in] [, options] for the st commands; returns (varlist, mask, options)."""
    from ..commands._util import touse
    from ..lang.syntax import parse_standard
    p = parse_standard(args)
    o = match_options(p.options, spec) if p.options.strip() else {}
    return p.varlist.strip(), touse(s, p), o


# ---------------------------------------------------------------------------
# stdescribe
# ---------------------------------------------------------------------------

def cmd_stdescribe(s: "Session", args: str) -> None:
    _vl, mask, o = _qualifiers(s, args, {"weight": 6, "noshow": 6})
    info = require_st(s)
    if not o.get("noshow"):
        header(s, info)
    sd = st_data(s, mask)
    out = s.output
    if sd.ids is not None:
        keys, inv = np.unique(sd.ids, return_inverse=True)
    else:
        keys, inv = np.arange(len(sd.t)), np.arange(len(sd.t))
    G = len(keys)
    nrec = np.bincount(inv, minlength=G).astype(float)
    entry = np.full(G, np.inf)
    np.minimum.at(entry, inv, sd.t0)
    exitt = np.full(G, -np.inf)
    np.maximum.at(exitt, inv, sd.t)
    risk = np.bincount(inv, weights=sd.t - sd.t0, minlength=G)
    fails = np.bincount(inv, weights=sd.d, minlength=G)
    span = exitt - entry
    gap = span - risk
    has_gap = gap > 1e-12 * np.maximum(span, 1)
    out.write(f"{'':35}|-------------- per subject --------------|\n", "text")
    out.write("Category                   total        mean         min     median        max\n", "text")
    out.write(LINE + "\n", "text")

    def row(lab, total, vals=None):
        out.write(f"{lab:<22}", "text")
        out.write(f"{_g(total, '%10.0g') if total is not None else '':>10}", "result")
        if vals is not None and not len(vals):
            out.write(f" {'.':>11} {'.':>11} {'.':>10} {'.':>10}", "result")
        elif vals is not None:
            out.write(f" {_g(float(np.mean(vals))):>11} {_g(float(np.min(vals))):>11}"
                      f" {_g(float(np.median(vals))):>10} {_g(float(np.max(vals))):>10}", "result")
        out.write("\n", "text")

    row("no. of subjects", float(G))
    row("no. of records", float(len(sd.t)), nrec)
    out.write("\n", "text")
    row("(first) entry time", None, entry)
    row("(final) exit time", None, exitt)
    out.write("\n", "text")
    row("subjects with gap", float(has_gap.sum()))
    # with multiple-record data the row shows "." statistics when no subject has a gap
    # (compat 0901); with one record per subject it shows only the total
    row("time on gap if gap", float(gap[has_gap].sum()),
        gap[has_gap] if has_gap.any() or sd.ids is not None else None)
    row("time at risk", float(risk.sum()), risk)
    out.write("\n", "text")
    row("failures", float(fails.sum()), fails)
    out.write(LINE + "\n", "text")
    s.r = {"N_sub": float(G), "N_total": float(len(sd.t)), "N_fail": float(fails.sum()),
           "risk": float(risk.sum()), "gap": float(gap[has_gap].sum()), "N_gap": float(has_gap.sum()),
           "t0_min": float(entry.min()) if G else SYS, "t_max": float(exitt.max()) if G else SYS}


# ---------------------------------------------------------------------------
# stsum
# ---------------------------------------------------------------------------

def cmd_stsum(s: "Session", args: str) -> None:
    from .km import km, quantile_time
    _vl, mask, o = _qualifiers(s, args, {"by": 2, "noshow": 6})
    info = require_st(s)
    if not o.get("noshow"):
        header(s, info)
    out = s.output
    byv = str(o["by"]).split() if isinstance(o.get("by"), str) else []
    if len(byv) > 1:
        raise StataError(103, "too many variables specified")   # VERIFY: stsum accepts by(varlist)
    sd_all = st_data(s, mask)
    rows = []
    if byv:
        from ..core.varlist import resolve_name
        bv = resolve_name(s.data, byv[0])
        vals = s.data.get(bv).data
        sub = vals[sd_all.mask]
        levels = sorted(set(sub.tolist()) - ({SYS} if not s.data.get(bv).is_string else set()),
                        key=lambda x: (isinstance(x, str), x))
        for lv in levels:
            m = sd_all.mask & (vals == lv)
            rows.append((_level_label(s, bv, lv), st_data(s, m)))
    rows.append(("total", sd_all))
    head = (byv[0][:8] if byv else "")
    out.write(f"{'':9}|               incidence       no. of    |------ Survival time -----|\n", "text")
    out.write(f"{head:<9}| time at risk     rate        subjects        25%       50%       75%\n", "text")
    out.write("-" * 9 + "+" + "-" * 69 + "\n", "text")
    results = {}
    for i, (lab, sd) in enumerate(rows):
        if lab == "total" and byv:
            out.write("-" * 9 + "+" + "-" * 69 + "\n", "text")
        risk = float(np.sum((sd.t - sd.t0) * sd.w))
        nf = float(np.sum(sd.d * sd.w))
        rate = nf / risk if risk > 0 else SYS
        nsub = n_subjects(sd)
        k = km(sd.t, sd.t0, sd.d, sd.w)
        qs = [quantile_time(k, p) for p in (0.25, 0.5, 0.75)]
        out.write(f"{lab:>8} |", "text")
        out.write(f" {_g(risk, '%12.0g'):>12} {_g(rate):>10} {_g(nsub, '%12.0g'):>13}"
                  f" {_g(qs[0]):>10} {_g(qs[1]):>9} {_g(qs[2]):>9}\n", "result")
        results = {"risk": risk, "ir": rate, "N_sub": nsub, "p25": qs[0], "p50": qs[1], "p75": qs[2]}
    s.r = results


def _level_label(s: "Session", var: str, lv) -> str:
    v = s.data.get(var)
    if isinstance(lv, str):
        return lv[:8]
    if v.value_label:
        lab = s.data.value_labels.get(v.value_label, {}).get(int(lv)) if float(lv).is_integer() else None
        if lab:
            return lab[:8]
    return _g(lv)
