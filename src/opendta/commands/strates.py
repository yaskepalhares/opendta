"""Rates and life tables for survival data: strate, stptime, stci, stir, stmh, stmc and ltable
([ST] strate, [ST] stptime, [ST] stci, [ST] stir, [ST] stmh, [ST] ltable).

The statistics are in core/survival.py and core/epi.py; this module reads the data, writes the
tables and stores r(). Layouts the manuals do not show carry the VERIFY mark
(compat/do/1901_survival.do)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import epi as E
from ..core import missing as M
from ..core import survival as SV
from ..core.dataset import Dataset, Variable
from ..core.errors import StataError
from ..estimation.km import km, km_ci, quantile_time, _by_groups
from ..estimation.st import header, require_st, st_data
from ..estimation.stfit import fmt
from ..lang.syntax import match_options, parse_standard
from ..lang.words import parse_numlist
from ._util import touse
from .summarize import _level, weights

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
NAN = float("nan")


def _not_available(o, names) -> None:
    for k in names:
        if o.get(k):
            raise StataError(198, f"option {k} is not available yet in OpenDTA")


def _sysmis(x: float) -> float:
    return SYS if x is None or not math.isfinite(x) else float(x)


def _label_text(s: "Session", names, key_rows) -> str:
    from .inspect import cell_text
    ds = s.data
    return ", ".join(cell_text(ds, ds.get(n), key_rows) for n in names)


def _groups(s: "Session", names: list[str], rows: np.ndarray):
    """[(positions in rows, first dataset row)] for each combination of the variables in `names`
    (missing values excluded; VERIFY), in sorted order. Without variables, a single group."""
    if not names:
        return [(np.arange(len(rows)), int(rows[0]) if len(rows) else 0)]
    ds = s.data
    cols = [ds.get(n) for n in names]
    keys: dict = {}
    for j, i in enumerate(rows):
        key, ok = [], True
        for c in cols:
            v = c.raw[i] if c.is_string else float(c.data[i])
            ok &= (v != "") if c.is_string else (v < SYS)
            key.append(v)
        if ok:
            keys.setdefault(tuple(key), []).append(j)
    return [(np.array(v), int(rows[v[0]])) for _, v in sorted(keys.items())]


def _result_dataset(s: "Session", names: list[str], firsts: list[int], cols) -> Dataset:
    """The by variables at the first observation of each group followed by double columns
    [(name, values, format)]."""
    ds = s.data
    new = Dataset()
    new.nobs = len(firsts)
    idx = np.array(firsts, dtype=int)
    used = set()
    for nm in names:
        v = ds.get(nm)
        nv = Variable(nm, v.vtype, np.empty(0), fmt=v.fmt, label=v.label, value_label=v.value_label)
        nv.raw = v.raw[idx].copy()
        new.vars.append(nv)
        if v.value_label:
            used.add(v.value_label)
    new.value_labels = {k: dict(ds.value_labels[k]) for k in used if k in ds.value_labels}
    for nm, vals, f in cols:
        new.vars.append(Variable(nm, "double", np.array([_sysmis(x) for x in vals], dtype=np.float64), fmt=f))
    return new


def _show_or_save(s: "Session", new: Dataset, o) -> None:
    """Lists the result dataset (noobs) or, with output(), saves it and leaves the data alone."""
    if o.get("output"):
        from ..io.dta import write_dta
        from .files import dta_path
        head, _, rest = str(o["output"]).partition(",")
        path = dta_path(head)
        if path.exists() and "replace" not in rest:
            raise StataError(602, f"file {path} already exists")
        write_dta(new, path)
        return
    from .inspect import cmd_list
    old = s.data
    s.data = new
    try:
        cmd_list(s, ", noobs")
    finally:
        s.data = old
        s.notify_state()


def _cluster_of(s: "Session", sd, o) -> np.ndarray:
    if o.get("cluster"):
        name = s.expand_varlist(str(o["cluster"]))[0]
        v = s.data.get(name)
        return np.array([str(x) for x in v.raw[sd.mask]]) if v.is_string else v.data[sd.mask]
    if sd.ids is not None:
        return sd.ids
    raise StataError(198, "jackknife requires cluster() or data stset with id()")   # VERIFY


def _rate_ci(D, T, scale, level, jack, num, den, cluster):
    """Rate = failures per `T`/scale of time and its limits: Poisson (ln rate with se 1/sqrt(D))
    or the jackknife over clusters (num, den and cluster are per observation)."""
    if jack:
        r, lo, hi = SV.jackknife_log_ci(num, den, cluster, level)
        return r * scale, lo * scale, hi * scale
    rate = D / T * scale if T > 0 else NAN
    lo, hi = SV.poisson_rate_ci(D, rate, level)
    return rate, lo, hi


# ---------------------------------------------------------------------------
# strate
# ---------------------------------------------------------------------------

_STRATE_OPTS = {"per": 3, "smr": 3, "cluster": 2, "jackknife": 4, "output": 3, "level": 3,
                "graph": 5, "noshow": 6}


def cmd_strate(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _STRATE_OPTS)
    _not_available(o, ("graph",))
    info = require_st(s)
    level = _level(s, o)
    per = float(o["per"]) if o.get("per") else 1.0
    if o.get("smr") and not o.get("per"):
        raise StataError(198, "option smr() requires per()")   # VERIFY
    if o.get("cluster") and not o.get("jackknife"):
        raise StataError(198, "option cluster() requires jackknife")   # VERIFY
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    sd = st_data(s, touse(s, p))
    rows = np.flatnonzero(sd.mask)
    ref = None
    if o.get("smr"):
        ref = s.data.get(s.expand_varlist(str(o["smr"]))[0]).data[sd.mask]
    T = (sd.t - sd.t0) * sd.w
    Dn = sd.d * sd.w
    cluster = _cluster_of(s, sd, o) if o.get("jackknife") else None
    z = SV._z(level)
    firsts, vals = [], []
    for pos, first in _groups(s, names, rows):
        D, Tg = float(Dn[pos].sum()), float(T[pos].sum())
        if ref is not None:
            Eg = float((ref[pos] * T[pos]).sum()) / per
            smr = D / Eg if Eg > 0 else NAN
            f = math.exp(z / math.sqrt(D)) if D > 0 else NAN
            vals.append((D, Eg, smr, smr / f, smr * f))
        else:
            rate, lo, hi = _rate_ci(D, Tg / per, 1.0, level, bool(o.get("jackknife")), Dn[pos],
                                    T[pos] / per, cluster[pos] if cluster is not None else None)
            vals.append((D, Tg / per, rate, lo, hi))
        firsts.append(first)
    if ref is not None:
        cn = ("_D", "_E", "_SMR", "_Lower", "_Upper")
    else:
        cn = ("_D", "_Y", "_Rate", "_Lower", "_Upper")
    cols = [(c, [v[k] for v in vals], "%9.0g") for k, c in enumerate(cn)]
    new = _result_dataset(s, names, firsts, cols)
    if not o.get("noshow"):
        header(s, info)
    out = s.output
    kind = "SMRs" if ref is not None else "rates"
    out.write(f"Estimated {kind} and lower/upper bounds of {level:g}% confidence intervals\n", "text")   # VERIFY
    out.write(f"({int(sd.w.sum()):,} records included in the analysis)\n\n", "text")   # VERIFY
    _show_or_save(s, new, o)


# ---------------------------------------------------------------------------
# stptime
# ---------------------------------------------------------------------------

_STPTIME_OPTS = {"per": 3, "by": 2, "title": 3, "at": 2, "cluster": 2, "jackknife": 4,
                 "level": 3, "noshow": 6, "smr": 3, "dd": 2, "output": 3}


def cmd_stptime(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _STPTIME_OPTS)
    _not_available(o, ("smr", "dd", "output"))
    info = require_st(s)
    level = _level(s, o)
    per = float(o["per"]) if o.get("per") else 1.0
    if o.get("cluster") and not o.get("jackknife"):
        raise StataError(198, "option cluster() requires jackknife")   # VERIFY
    names = s.expand_varlist(str(o["by"])) if o.get("by") else []
    sd = st_data(s, touse(s, p))
    rows = np.flatnonzero(sd.mask)
    cuts = sorted(parse_numlist(str(o["at"]))) if o.get("at") else []
    cluster = _cluster_of(s, sd, o) if o.get("jackknife") else None
    title = str(o["title"]).strip('"') if isinstance(o.get("title"), str) else "person-time"
    if not o.get("noshow"):
        header(s, info)
    out = s.output
    lw = 12
    out.write(f"{'Cohort':>{lw}} |{title.capitalize():>13}{'Failures':>11}{'Rate':>12}"
              f"   [{level:g}% Conf. Interval]\n", "text")   # VERIFY
    out.write("-" * (lw + 1) + "+" + "-" * (58 + len(f"{level:g}") - 2) + "\n", "text")

    def piece(label: str, pos: np.ndarray, lo: float | None = None, hi: float | None = None):
        t, t0, d, w = sd.t[pos], sd.t0[pos], sd.d[pos], sd.w[pos]
        if lo is None:
            tt, dd = (t - t0) * w, d * w
        else:
            top = np.minimum(t, hi if hi is not None else np.inf)
            tt = np.maximum(top - np.maximum(t0, lo), 0) * w
            dd = d * w * (t > lo) * ((t <= hi) if hi is not None else True)
        D, T = float(dd.sum()), float(tt.sum())
        rate, l_, u_ = _rate_ci(D, T, per, level, bool(o.get("jackknife")), dd, tt,
                                cluster[pos] if cluster is not None else None)
        out.write(f"{label:>{lw}} |", "text")
        out.write(f"{fmt(T, '%12.0g'):>13}{fmt(D, '%9.0g'):>11}{fmt(rate, '%10.0g'):>12}"
                  f"{fmt(l_, '%10.0g'):>13}{fmt(u_, '%10.0g'):>11}\n", "result")
        return T, D, rate, l_, u_

    bounds = [(None, None)] if not cuts else (
        [(a, b) for a, b in zip(cuts[:-1], cuts[1:])] + [(cuts[-1], None)])
    for pos, first in (_groups(s, names, rows) if names else []):
        base = _label_text(s, names, first)
        for a, b in bounds:
            lab = base if a is None else f"{base} {_interval_label(a, b)}"
            piece(lab, pos, a, b)
    if names:
        out.write("-" * (lw + 1) + "+" + "-" * (58 + len(f"{level:g}") - 2) + "\n", "text")
    everyone = np.arange(len(rows))
    if cuts:
        for a, b in bounds:
            res = piece(_interval_label(a, b), everyone, a, b)
    res = piece("total", everyone)
    s.r = {"ptime": res[0], "failures": res[1], "rate": _sysmis(res[2]),
           "lb": _sysmis(res[3]), "ub": _sysmis(res[4])}


def _interval_label(a: float, b: float | None) -> str:
    return f"({a:g} - {b:g}]" if b is not None else f"> {a:g}"   # VERIFY


# ---------------------------------------------------------------------------
# stci
# ---------------------------------------------------------------------------

_STCI_OPTS = {"by": 2, "median": 3, "p": 1, "rmean": 2, "emean": 2, "ccorr": 2, "level": 3,
              "graph": 5, "noshow": 6}


def cmd_stci(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _STCI_OPTS)
    _not_available(o, ("graph",))
    info = require_st(s)
    level = _level(s, o)
    if sum(bool(o.get(k)) for k in ("median", "p", "rmean", "emean")) > 1:
        raise StataError(198, "options median, p(), rmean, and emean are mutually exclusive")   # VERIFY
    pct = float(o["p"]) if isinstance(o.get("p"), str) else 50.0
    if not 0 < pct < 100:
        raise StataError(198, "p() invalid")   # VERIFY
    mode = "rmean" if o.get("rmean") else "emean" if o.get("emean") else "p"
    byvars = s.expand_varlist(str(o["by"])) if o.get("by") else []
    mask = touse(s, p)
    groups = _by_groups(s, mask, byvars)
    if byvars:
        groups.append(("total", st_data(s, mask)))
    z = SV._z(level)
    if not o.get("noshow"):
        header(s, info)
    out = s.output
    lw = 12
    head = {"p": f"{pct:g}%", "rmean": "restricted mean", "emean": "extended mean"}[mode]
    out.write(f"{'':>{lw}} |       no. of\n", "text")
    extra = "" if mode == "emean" else f"   Std. Err.   [{level:g}% Conf. Interval]"
    out.write(f"{'':>{lw}} |      subjects{head:>17}{extra}\n", "text")   # VERIFY
    out.write("-" * (lw + 1) + "+" + "-" * 60 + "\n", "text")
    res = {}
    underestimated = False
    for lab, sd in groups:
        label = ", ".join(x.split(" = ")[-1] for x in lab.split(", ")) if lab else "total"
        k = km(sd.t, sd.t0, sd.d, sd.w)
        nsub = float(sd.w.sum()) if sd.ids is None else float(len(np.unique(sd.ids)))
        star = ""
        if mode == "p":
            q = 1 - pct / 100
            tp = quantile_time(k, pct / 100)
            _se, slo, shi = km_ci(k.S, k.gw, level)
            slo, shi = np.where(k.S <= 0, 0.0, slo), np.where(k.S <= 0, 0.0, shi)

            def first(arr):
                hit = np.flatnonzero((k.d > 0) & (arr <= q + 1e-12))
                return float(k.times[hit[0]]) if hit.size else SYS
            lb, ub = first(slo), first(shi)   # S limits reach 1 - p at the earlier / later time
            se = SV.km_percentile_se(k.times, k.S, k.gw, tp, pct / 100) if tp < SYS else NAN
            est = tp
            res = {"N_sub": nsub, f"p{pct:g}": tp, "lb": lb, "ub": ub, "se": _sysmis(se)}
        else:
            tmax = float(sd.t.max())
            mean, se = SV.km_rmean(k.times, k.S, k.d, k.n, tmax, bool(o.get("ccorr")))
            if mode == "emean":
                mean, se = SV.km_emean(k.times, k.S, mean, tmax), NAN
            est, lb, ub = mean, mean - z * se, mean + z * se
            if sd.d[sd.t == tmax].sum() == 0 and mode == "rmean":
                star, underestimated = "*", True
            res = {"N_sub": nsub, mode: mean, "se": _sysmis(se)}
        out.write(f"{label:>{lw}} |", "text")
        out.write(f"{fmt(nsub, '%12.0g'):>13}{fmt(est, '%10.0g') + star:>{17 if mode != 'p' else 10}}", "result")
        if mode != "emean":
            out.write(f"{fmt(se, '%10.0g'):>12}{fmt(lb, '%10.0g'):>13}{fmt(ub, '%10.0g'):>11}", "result")
        out.write("\n", "result")
    if underestimated:
        out.write("   (*) largest observed analysis time is censored, mean is underestimated.\n", "text")   # VERIFY
    s.r = {k: _sysmis(v) for k, v in res.items()}   # of the last row (total with by())


# ---------------------------------------------------------------------------
# stir, stmh, stmc
# ---------------------------------------------------------------------------

def _rate_options(p, extra=None):
    from .epitab import _IR_OPTS, _check_by, _not_yet
    spec = dict(_IR_OPTS)
    spec["noshow"] = 6
    spec.update(extra or {})
    o = match_options(p.options, spec)
    _not_yet(o)
    _check_by(o)
    return o


def _exposure(s: "Session", p, sd):
    """The two-level exposure variable restricted to the sample: (name, values, exposed value,
    unexposed value). Exposed is the larger value (VERIFY)."""
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(198, "invalid syntax" if not names else "too many variables specified")
    v = s.data.get(names[0])
    if v.is_string:
        raise StataError(109, "type mismatch")
    x = v.data[sd.mask]
    levels = sorted(set(x[x < SYS].tolist()))
    if len(levels) != 2:
        raise StataError(198, f"{names[0]} must have exactly two nonmissing values")   # VERIFY
    return names[0], x, levels[1], levels[0]


def _cells(sd, x, hi, lo, sel):
    a = float((sd.d * sd.w)[sel & (x == hi)].sum())
    b = float((sd.d * sd.w)[sel & (x == lo)].sum())
    t1 = float(((sd.t - sd.t0) * sd.w)[sel & (x == hi)].sum())
    t0 = float(((sd.t - sd.t0) * sd.w)[sel & (x == lo)].sum())
    return a, b, t1, t0


def _strata_of(s: "Session", sd, o):
    """[(label, selection over the sample)] for by(); one stratum with everyone otherwise."""
    from .epitab import _strata
    if not o.get("by"):
        return None, [("", np.ones(len(sd.t), dtype=bool))]
    name, st = _strata(s, str(o["by"]), sd.mask)
    return name, [(lab, sel[sd.mask]) for lab, sel in st]


def cmd_stir(s: "Session", args: str) -> None:
    from .epitab import _ir_by, _ir_table, _std_spec
    p = parse_standard(args)
    o = _rate_options(p)
    info = require_st(s)
    sd = st_data(s, touse(s, p))
    name, x, hi, lo = _exposure(s, p, sd)
    if not o.get("noshow"):
        header(s, info)
    s.output.write(f"note: Exposed <-> {name}=={hi:g} and Unexposed <-> {name}=={lo:g}\n", "text")   # VERIFY
    byname, strata = _strata_of(s, sd, o)
    if byname is None:
        _ir_table(s, _cells(sd, x, hi, lo, strata[0][1]), o, name, ("Failure", "Time"))
        return
    cells = [(lab, _cells(sd, x, hi, lo, sel)) for lab, sel in strata]
    full = [(lab, sel) for lab, sel in strata]
    _ir_by(s, byname, cells, o, _std_spec(s, o, full))


def _mh_report(s: "Session", title: str, cells, level: float) -> None:
    est, lo, hi = E.mh_irr(cells, level)
    num = sum(a for a, _b, _t1, _t0 in cells)
    ex = sum((a + b) * t1 / (t1 + t0) for a, b, t1, t0 in cells if t1 + t0)
    var = sum((a + b) * t1 * t0 / (t1 + t0) ** 2 for a, b, t1, t0 in cells if t1 + t0)
    chi2 = (num - ex) ** 2 / var if var > 0 else NAN
    pv = E.st.chi2.sf(chi2, 1) if chi2 == chi2 else NAN
    out = s.output
    out.write(f"\n{title}\n\n", "text")
    out.write(f"{'Rate ratio':>12}{'chi2(1)':>12}{'Pr>chi2':>10}   [{level:g}% Conf. Interval]\n", "text")   # VERIFY
    out.write(f"{fmt(est, '%10.0g'):>12}{fmt(chi2, '%9.2f'):>12}{fmt(pv, '%7.4f'):>10}"
              f"{fmt(lo, '%10.0g'):>13}{fmt(hi, '%10.0g'):>11}\n", "result")
    s.r = {"rr": _sysmis(est), "lb": _sysmis(lo), "ub": _sysmis(hi), "chi2": _sysmis(chi2), "p": _sysmis(pv)}


def cmd_stmh(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = _rate_options(p)
    info = require_st(s)
    sd = st_data(s, touse(s, p))
    name, x, hi, lo = _exposure(s, p, sd)
    if not o.get("noshow"):
        header(s, info)
    _, strata = _strata_of(s, sd, o)
    cells = [_cells(sd, x, hi, lo, sel) for _lab, sel in strata]
    _mh_report(s, "Mantel-Haenszel estimate of the rate ratio", cells, _level(s, o))   # VERIFY


def cmd_stmc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = _rate_options(p)
    info = require_st(s)
    sd = st_data(s, touse(s, p))
    name, x, hi, lo = _exposure(s, p, sd)
    if not o.get("noshow"):
        header(s, info)
    _, strata = _strata_of(s, sd, o)
    cells = []
    for _lab, sel in strata:
        for tf in np.unique(sd.t[sel & (sd.d > 0)]):
            risk = sel & (sd.t0 < tf) & (sd.t >= tf)
            ev = sel & (sd.t == tf) & (sd.d > 0)
            # the risk set counts play the role of person-time (Mantel-Cox; VERIFY the variance)
            cells.append((float(sd.w[ev & (x == hi)].sum()), float(sd.w[ev & (x == lo)].sum()),
                          float(sd.w[risk & (x == hi)].sum()), float(sd.w[risk & (x == lo)].sum())))
    _mh_report(s, "Mantel-Cox estimate of the rate ratio", cells, _level(s, o))   # VERIFY


# ---------------------------------------------------------------------------
# ltable
# ---------------------------------------------------------------------------

_LTABLE_OPTS = {"by": 2, "intervals": 1, "level": 3, "survival": 2, "failure": 1, "hazard": 1,
                "noadjust": 3, "notable": 3, "test": 2, "graph": 5, "saving": 3, "tvid": 2}


def _four(x: float) -> str:
    return "." if x is None or not math.isfinite(x) else f"{x:.4f}"


def cmd_ltable(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _LTABLE_OPTS)
    _not_available(o, ("graph", "saving", "tvid"))
    level = _level(s, o)
    names = s.expand_varlist(p.varlist)
    if not 1 <= len(names) <= 2:
        raise StataError(198, "invalid syntax" if not names else "too many variables specified")
    ds = s.data
    mask = touse(s, p)
    w, _wt, mask = weights(s, p, mask, ("fweight",))
    tv = ds.get(names[0])
    mask = mask & (tv.data < SYS)
    dv = ds.get(names[1]) if len(names) == 2 else None
    if dv is not None:
        mask = mask & (dv.data < SYS)
    rows = np.flatnonzero(mask)
    if not len(rows):
        raise StataError(2000, "no observations")
    t = tv.data[rows]
    if np.any(t < 0):
        raise StataError(198, f"{tv.name} has negative values")   # VERIFY
    dead = dv.data[rows] if dv is not None else np.ones(len(rows))
    wt = w[rows] if w is not None else np.ones(len(rows))
    spec: float | list = 1.0
    if o.get("intervals"):
        nums = [float(x) for x in parse_numlist(str(o["intervals"]))]
        spec = nums[0] if len(nums) == 1 else nums
        if np.isscalar(spec) and spec <= 0:
            raise StataError(198, "intervals() invalid")   # VERIFY
    kind = "failure" if o.get("failure") else "hazard" if o.get("hazard") else "survival"
    adjust = not o.get("noadjust")
    byname = s.expand_varlist(str(o["by"]))[0] if o.get("by") else ""
    groups = _groups(s, [byname] if byname else [], rows)
    out = s.output
    for pos, first in groups:
        if byname:
            out.write(f"\n{byname} = {_label_text(s, [byname], first)}\n", "text")   # VERIFY
        tab = SV.life_table(t[pos], dead[pos], wt[pos], spec, adjust, level)
        if o.get("notable"):
            continue
        _ltable_rows(s, tab, kind, level)
    if o.get("test"):
        if not byname:
            raise StataError(198, "option test requires by()")   # VERIFY
        g = np.zeros(len(rows), dtype=int)
        for k, (pos, _f) in enumerate(groups):
            g[pos] = k
        G = len(groups)
        d_g = np.array([float((wt * (dead != 0))[g == k].sum()) for k in range(G)])
        T_g = np.array([float((wt * t)[g == k].sum()) for k in range(G)])
        for title, (chi2, df) in (("Log-rank test for equality of survivor functions",
                                   SV.logrank_chi2(t, dead, wt, g, G)),
                                  ("Likelihood-ratio test of equal exponential hazards",
                                   SV.lr_exponential(d_g, T_g))):   # VERIFY: the tests of ltable
            pv = float(E.st.chi2.sf(chi2, df))
            out.write(f"\n{title}\n", "text")
            out.write(f"{'chi2(' + str(df) + ') =':>16}", "text")
            out.write(f"{chi2:>9.2f}", "result")
            out.write("  Pr>chi2 = ", "text")
            out.write(f"{pv:.4f}\n", "result")
        s.r = {"chi2": _sysmis(chi2), "df": float(df), "p": _sysmis(pv)}


def _ltable_rows(s: "Session", tab, kind: str, level: float) -> None:
    out = s.output
    if kind == "hazard":
        cols = ("Hazard", "Std. Err.")
    else:
        cols = ("Survival" if kind == "survival" else "Cum. Fail.", "Std. Err.")
    out.write(f"{'Beg.':>15}\n", "text")
    out.write(f"{'Interval':>12}{'Total':>8}{'Deaths':>9}{'Lost':>7}{cols[0]:>10}{'Error':>10}"
              f"{'[' + format(level, 'g') + '% Conf. Int.]':>22}\n", "text")   # VERIFY
    out.write("-" * 79 + "\n", "text")
    for r in tab:
        hi = f"{r['hi']:g}" if r["hi"] is not None else "+"
        if kind == "survival":
            est, se, lb, ub = r["S"], r["se"], r["lb"], r["ub"]
        elif kind == "failure":
            est, se = 1 - r["S"], r["se"]
            lb, ub = 1 - r["ub"], 1 - r["lb"]
        else:
            est, se, lb, ub = r["h"], r["hse"], r["hlb"], r["hub"]
        out.write(f"{r['lo']:>6g}{hi:>6}{r['N']:>8g}{r['d']:>9g}{r['m']:>7g}", "result")
        out.write(f"{_four(est):>10}{_four(se):>10}{_four(lb):>11}{_four(ub):>11}\n", "result")
    out.write("-" * 79 + "\n", "text")
