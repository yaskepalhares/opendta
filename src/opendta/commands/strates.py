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


def _sig(x: float, digits: int) -> str:
    """%.{digits}g without the leading zero of numbers below one; . for missing."""
    if x is None or not math.isfinite(x):
        return "."
    t = f"{x:.{digits}g}"
    return t[1:] if t.startswith("0.") else "-" + t[2:] if t.startswith("-0.") else t


def _int_type(vals) -> str:
    """Smallest integer type that holds the values (double when any is not an integer)."""
    a = np.array([v for v in vals if v is not None and math.isfinite(v)], dtype=float)
    if len(a) and np.all(a == np.floor(a)):
        m = float(np.max(np.abs(a)))
        return "byte" if m <= 100 else "int" if m <= 32740 else "long" if m <= 2147483620 else "double"
    return "double"


def _digits(x: float) -> int:
    return len(str(int(abs(x)))) if math.isfinite(x) else 1


def _rate_formats(per: float, y, rates) -> tuple[str, str]:
    """Display formats of the person-time and of the rates in the strate table (VERIFY for per()
    other than 1 and 100: the decimals of the rates fall by one for each power of ten)."""
    d = max(0, 6 - int(round(math.log10(per)))) if per > 0 else 6
    wy = _digits(max([v for v in y if math.isfinite(v)], default=0)) + 5
    wr = _digits(max([v for v in rates if math.isfinite(v)], default=0)) + 1 + d
    return f"%{wy}.4f", f"%{wr}.{d}f"


def _show_or_save(s: "Session", shown: Dataset, o, saved: Dataset | None = None) -> None:
    """Lists the result dataset (noobs); with output(), also saves `saved` and leaves the data alone."""
    from .inspect import cmd_list
    old = s.data
    s.data = shown
    try:
        cmd_list(s, ", noobs")
    finally:
        s.data = old
        s.notify_state()
    if o.get("output"):
        from ..io.dta import write_dta
        from .files import dta_path
        head, _, rest = str(o["output"]).partition(",")
        path = dta_path(head)
        if path.exists() and "replace" not in rest:
            raise StataError(602, f"file {path} already exists")
        write_dta(saved, path)


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
        cn = ("D", "E", "SMR", "Lower", "Upper")
    else:
        cn = ("D", "Y", "Rate", "Lower", "Upper")
    fy, fr = _rate_formats(per, [v[1] for v in vals], [x for v in vals for x in v[2:]])
    shown = _result_dataset(s, names, firsts, [
        (cn[0], [v[0] for v in vals], "%4.0f"), (cn[1], [v[1] for v in vals], fy),
        *[(c, [v[k] for v in vals], fr) for k, c in enumerate(cn) if k >= 2]])
    saved = None
    if o.get("output"):
        lv = f"{level:g}"
        labels = ("# Failures", "Person-time observation" if ref is None else "Expected failures",
                  "Rate estimate" if ref is None else "SMR estimate",
                  f"Lower {lv}% confidence limit", f"Upper {lv}% confidence limit")
        saved = _result_dataset(s, names, firsts, [
            ("_" + c, [v[k] for v in vals], f) for k, (c, f) in enumerate(
                zip(cn, ("%4.0f", "%8.2f", "%8.6f", "%8.6f", "%8.6f")))])
        for var, lab in zip(saved.vars[len(names):], labels):
            var.label = lab
        dv = saved.vars[len(names)]
        dv.vtype = _int_type([v[0] for v in vals])
    if not o.get("noshow"):
        header(s, info)
    out = s.output
    kind = "SMRs" if ref is not None else "rates"
    per_txt = f" (per {per:g})" if o.get("per") and per != 1 else ""
    out.write(f"Estimated {kind}{per_txt} and lower/upper bounds of {level:g}% confidence intervals\n", "text")
    out.write(f"({int(sd.w.sum()):,} records included in the analysis)\n", "text")
    if o.get("jackknife"):
        keys = {tuple(r) for r in zip(cluster.tolist(), *[
            [str(x) for x in s.data.get(n).raw[rows]] for n in names])} if names else set(cluster.tolist())
        out.write("The jackknife was used to calculate confidence intervals\n", "text")
        out.write(f"{len(keys):<7}total clusters identified by "
                  f"{'cluster()' if o.get('cluster') else 'id'}{' and the strata variables' if names else ''}\n", "text")
    out.write("\n", "text")
    _show_or_save(s, shown, o, saved)
    out.write("\n\n", "text")


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
    lv = f"{level:g}"
    bylab = names[0] if len(names) == 1 else "Cohort" if not names else ", ".join(names)
    left = bool(cuts) and bool(names)
    rule = "-" * 11 + "+" + "-" * (59 + len(lv) - 2)
    out.write("\n", "text")
    out.write(f"{bylab:<11}|" if (left or not names) else f"{bylab:>10} |", "text")
    out.write(f"{title:>13}{'failures':>11}{'rate':>12}   [{lv}% Conf. Interval]\n", "text")
    out.write(rule + "\n", "text")

    def piece(label: str, pos: np.ndarray, lo: float | None = None, hi: float | None = None, lab_left=False):
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
        out.write(f"{label:<11}|" if lab_left else f"{label:>10} |", "text")
        out.write(f"{fmt(T, '%12.0g'):>13}{fmt(D, '%9.0g'):>11}{fmt(rate, '%10.0g'):>12}"
                  f"{fmt(l_, '%9.0g'):>11}{fmt(u_, '%9.0g'):>12}\n", "result")
        return T, D, rate, l_, u_

    bounds = [(None, None)] if not cuts else (
        [(a, b) for a, b in zip([0.0] + cuts[:-1], cuts)] + [(cuts[-1], None)])
    groups = _groups(s, names, rows) if names else []
    for k, (pos, first) in enumerate(groups):
        base = _label_text(s, names, first)
        if cuts:
            out.write(f"{base:<11}|\n", "text")
            for a, b in bounds:
                piece(_interval_label(a, b), pos, a, b, lab_left=True)
        else:
            piece(base, pos)
        if k < len(groups) - 1 and cuts:
            out.write(rule + "\n", "text")
    if names:
        out.write(rule + "\n", "text")
    everyone = np.arange(len(rows))
    if cuts and not names:
        for a, b in bounds:
            piece(_interval_label(a, b), everyone, a, b, lab_left=True)
        out.write(rule + "\n", "text")
    res = piece("total", everyone)
    s.r = {"ptime": res[0], "failures": res[1], "rate": _sysmis(res[2]),
           "lb": _sysmis(res[3]), "ub": _sysmis(res[4])}


def _interval_label(a: float, b: float | None) -> str:
    if b is None:
        return f"     >{a:>4g} "
    return f"({a:<3g} - {b:>3g}]"   # VERIFY


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
    lv = f"{level:g}"
    bylab = ", ".join(byvars)
    out.write("\n", "text")
    if mode == "emean":
        out.write(f"{'':>12} |    no. of    extended\n{bylab:<12} |  subjects        mean\n", "text")
        rule = "-" * 13 + "+" + "-" * 22
    else:
        head = {"p": f"{pct:g}%", "rmean": "restricted", "emean": "extended"}[mode]
        if mode == "p":
            out.write(f"{'':>12} |    no. of \n", "text")
            out.write(f"{bylab:<12} |  subjects{head:>12}{'Std. Err.':>14}     [{lv}% Conf. Interval]\n", "text")
        else:
            out.write(f"{'':>12} |    no. of  restricted\n", "text")
            out.write(f"{bylab:<12} |  subjects        mean      Std. Err.    [{lv}% Conf. Interval]\n", "text")
        rule = "-" * 13 + "+" + "-" * (61 + len(lv) - 2)
    out.write(rule + "\n", "text")
    res = {}
    underestimated = False
    for lab, sd in groups:
        if lab == "total" and byvars:
            out.write(rule + "\n", "text")
        label = ", ".join(x.split(" = ")[-1] for x in lab.split(", ")) if lab else "total"
        k = km(sd.t, sd.t0, sd.d, sd.w)
        nsub = float(sd.w.sum()) if sd.ids is None else float(len(np.unique(sd.ids)))
        star = "   "
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
                star, underestimated = "(*)", True
            res = {"N_sub": nsub, mode: mean, "se": _sysmis(se)}
        out.write(f"{label:>12} |", "text")
        if mode == "p":
            out.write(f"{fmt(nsub, '%10.0g'):>10}{_sig(est, 7):>12}{_sig(se, 7):>14}"
                      f"{_sig(lb, 7):>13}{_sig(ub, 7):>11}", "result")
        elif mode == "rmean":
            out.write(f"{fmt(nsub, '%10.0g'):>10}{_sig(est, 7):>12}{star}{_sig(se, 7):>11}"
                      f"{_sig(lb, 6):>13}{_sig(ub, 6):>11}", "result")
        else:
            out.write(f"{fmt(nsub, '%10.0g'):>10}{_sig(est, 7):>12}", "result")
        out.write("\n", "result")
    if underestimated:
        out.write("\n(*) largest observed analysis time is censored, mean is underestimated\n", "text")   # VERIFY
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


def _exposure(s: "Session", p, sd, allow_strata: bool = False):
    """The two-level exposure variable restricted to the sample: (name, values, exposed value,
    unexposed value). Exposed is the larger value (VERIFY)."""
    names = s.expand_varlist(p.varlist)
    if not names:
        raise StataError(198, "invalid syntax")
    if len(names) > 1 and not allow_strata:
        raise StataError(198, "too many variables specified")
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
    s.output.write(f"note:  Exposed <-> {name}=={hi:g} and Unexposed <-> {name}=={lo:g}\n", "text")   # VERIFY
    byname, strata = _strata_of(s, sd, o)
    if byname is None:
        _ir_table(s, _cells(sd, x, hi, lo, strata[0][1]), o, name, ("Failure", "Time"))
        return
    cells = [(lab, _cells(sd, x, hi, lo, sel)) for lab, sel in strata]
    full = [(lab, sel) for lab, sel in strata]
    _ir_by(s, byname, cells, o, _std_spec(s, o, full))


def _mh_report(s: "Session", lines: list[str], cells, level: float, rr_line: bool = True) -> None:
    est, lo, hi = E.mh_irr(cells, level)
    num = sum(a for a, _b, _t1, _t0 in cells)
    ex = sum((a + b) * t1 / (t1 + t0) for a, b, t1, t0 in cells if t1 + t0)
    var = sum((a + b) * t1 * t0 / (t1 + t0) ** 2 for a, b, t1, t0 in cells if t1 + t0)
    chi2 = (num - ex) ** 2 / var if var > 0 else NAN
    pv = E.st.chi2.sf(chi2, 1) if chi2 == chi2 else NAN
    out = s.output
    lv = f"{level:g}"
    out.write("\n" + "\n".join(lines) + "\n\n", "text")
    rule = "        " + "-" * (58 + len(lv) - 2)
    if rr_line:
        out.write(f"RR estimate, and lower and upper {lv}% confidence limits\n\n", "text")
    out.write(f"{rule}\n", "text")
    out.write(f"             RR       chi2        P>chi2      [{lv}% Conf. Interval]\n{rule}\n", "text")
    out.write(f"{est:>15.3f}{chi2:>11.2f}{pv:>14.4f}{lo:>14.3f}{hi:>11.3f}\n", "result")
    out.write(f"{rule}\n", "text")
    s.r = {"rr": _sysmis(est), "lb": _sysmis(lo), "ub": _sysmis(hi), "chi2": _sysmis(chi2), "p": _sysmis(pv)}


def _strata_cells(s: "Session", p, o, sd, names):
    """Strata selections: the combinations of the variables after the exposure, or by()."""
    if len(names) > 1:
        rows = np.flatnonzero(sd.mask)
        sels = []
        for pos, _first in _groups(s, names[1:], rows):
            sel = np.zeros(len(rows), dtype=bool)
            sel[pos] = True
            sels.append(sel)
        return sels
    return [sel for _lab, sel in _strata_of(s, sd, o)[1]]


def cmd_stmh(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = _rate_options(p)
    info = require_st(s)
    sd = st_data(s, touse(s, p))
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    name, x, hi, lo = _exposure(s, p, sd, True)
    if not o.get("noshow"):
        header(s, info)
    cells = [_cells(sd, x, hi, lo, sel) for sel in _strata_cells(s, p, o, sd, names)]
    lines = ["Mantel-Haenszel estimate of the rate ratio", f"  comparing {name}=={hi:g} vs. {name}=={lo:g}"]
    if len(names) > 1:
        lines.append(f"  controlling for {' '.join(names[1:])}")
    _mh_report(s, lines, cells, _level(s, o))


def cmd_stmc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = _rate_options(p)
    info = require_st(s)
    sd = st_data(s, touse(s, p))
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    name, x, hi, lo = _exposure(s, p, sd, True)
    if not o.get("noshow"):
        header(s, info)
    cells = []
    for sel in _strata_cells(s, p, o, sd, names):
        for tf in np.unique(sd.t[sel & (sd.d > 0)]):
            risk = sel & (sd.t0 < tf) & (sd.t >= tf)
            ev = sel & (sd.t == tf) & (sd.d > 0)
            # the risk set counts play the role of person-time (Mantel-Cox; VERIFY the variance)
            cells.append((float(sd.w[ev & (x == hi)].sum()), float(sd.w[ev & (x == lo)].sum()),
                          float(sd.w[risk & (x == hi)].sum()), float(sd.w[risk & (x == lo)].sum())))
    strata = f" and {' '.join(names[1:])}" if len(names) > 1 else ""
    lines = ["Mantel-Cox comparisons", "", "Mantel-Haenszel estimates of the rate ratio",
             f"   comparing {name}=={hi:g} vs. {name}=={lo:g}",
             f"   controlling for time (by clicks){strata}", "",
             f"Overall Mantel-Haenszel estimate, controlling for time{' and' if strata else ''}"]
    if strata:
        lines.append(" ".join(names[1:]))
    _mh_report(s, lines, cells, _level(s, o), rr_line=False)


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
    out.write("\n", "text")
    if not o.get("notable"):
        _ltable_head(s, kind, level)
    for pos, first in groups:
        if byname and not o.get("notable"):
            out.write(f"{byname} = {_label_text(s, [byname], first)}\n", "text")
        tab = SV.life_table(t[pos], dead[pos], wt[pos], spec, adjust, level)
        if not o.get("notable"):
            _ltable_rows(s, tab, kind)
    if not o.get("notable"):
        out.write("-" * 79 + "\n", "text")
    if o.get("test"):
        if not byname:
            raise StataError(198, "option test requires by()")   # VERIFY
        g = np.zeros(len(rows), dtype=int)
        for k, (pos, _f) in enumerate(groups):
            g[pos] = k
        G = len(groups)
        d_g = np.array([float((wt * (dead != 0))[g == k].sum()) for k in range(G)])
        T_g = np.array([float((wt * t)[g == k].sum()) for k in range(G)])
        lr, df_lr = SV.lr_exponential(d_g, T_g)
        chi2, df, Ob, Ex = SV.logrank_chi2(t, dead, wt, g, G)
        out.write("\n", "text")
        out.write(f"Likelihood-ratio test statistic of homogeneity (group={byname}):\n", "text")
        out.write(f"chi2( {df_lr} ) = {_sig(lr, 7)},   P = {_sig(float(E.st.chi2.sf(lr, df_lr)), 7)}\n \n", "result")
        out.write(f"Logrank test of homogeneity (group={byname}):\n\n\n", "text")
        out.write("Log-rank test for equality of survivor functions\n\n", "text")
        labs = [_label_text(s, [byname], f) for _p, f in groups]
        L = max(len(byname), len("Total"), *(len(x) for x in labs)) + 1
        out.write(f"{'':<{L}}|   Events         Events\n{byname:<{L}}|  observed       expected\n", "text")
        out.write("-" * L + "+" + "-" * 25 + "\n", "text")
        for lab, o_, e_ in zip(labs, Ob, Ex):
            out.write(f"{lab:<{L}}|{o_:>10g}{e_:>15.2f}\n", "result")
        out.write("-" * L + "+" + "-" * 25 + "\n", "text")
        out.write(f"{'Total':<{L}}|{Ob.sum():>10g}{Ex.sum():>15.2f}\n\n", "result")
        pv = float(E.st.chi2.sf(chi2, df))
        out.write(f"            chi2({df}) = {chi2:>10.2f}\n            Pr>chi2 = {pv:>10.4f}\n", "result")
        s.r = {"chi2": _sysmis(chi2), "df": float(df), "p": _sysmis(pv)}


def _ltable_head(s: "Session", kind: str, level: float) -> None:
    out = s.output
    lv = f"[{level:g}% Conf. Int.]"
    if kind == "hazard":
        out.write("                 Beg.     Cum.     Std.                Std.\n", "text")
        out.write(f"   Interval     Total   Failure   Error    Hazard     Error    {lv}\n", "text")
    else:
        c1 = "Survival" if kind == "survival" else "Failure"
        top = ("                 Beg.                                 Std." if kind == "survival" else
               "                 Beg.                      Cum.       Std.")
        out.write(top + "\n", "text")
        out.write(f"   Interval     Total   Deaths   Lost    {c1:<8}    Error     {lv}\n", "text")
    out.write("-" * 79 + "\n", "text")


def _ltable_rows(s: "Session", tab, kind: str) -> None:
    out = s.output
    for r in tab:
        hi = f"{r['hi']:g}" if r["hi"] is not None else "+"
        out.write(f"{r['lo']:>5g}{hi:>6}{r['N']:>10g}", "result")
        if kind == "hazard":
            cum, cse = 1 - r["S"], r["se"]
            out.write(f"{_four(cum):>10}{_four(cse):>8}{_four(r['h']):>10}{_four(r['hse']):>10}"
                      f"{_four(max(r['hlb'], 0.0) if math.isfinite(r['hlb']) else r['hlb']):>10}{_four(r['hub']):>10}\n", "result")
            continue
        out.write(f"{r['d']:>9g}{r['m']:>7g}", "result")
        if kind == "survival":
            est, se, lb, ub = r["S"], r["se"], r["lb"], r["ub"]
        else:
            est, se, lb, ub = 1 - r["S"], r["se"], 1 - r["ub"], 1 - r["lb"]
        out.write(f"{_four(est):>11}{_four(se):>10}{_four(lb):>11}{_four(ub):>10}\n", "result")
