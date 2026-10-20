"""Agreement, effect-size and transformation commands of the base-statistics block: ameans
([R] ameans), icc ([R] icc), loneway ([R] loneway), esize/esizei ([R] esize), lnskew0
([R] lnskew0), ladder ([R] ladder) and boxcox ([R] boxcox).

The computations are in opendta.core.basestats; layouts the manual does not show are marked
VERIFY and checked in compat/do/1903_base_stats.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import basestats as B
from ..core import missing as M
from ..core import stats as S
from ..core.dataset import Variable
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from ..lazy import lazy
from ._util import eval_vector, touse
from .nptests import _royston_adjust, _set_r, _sk_z
from .registry import command
from .summarize import _level, g9, name12, numeric_vars, weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _lvl_text(level: float) -> str:
    return f"{level:g}"


# ---------------------------------------------------------------------------
# ameans
# ---------------------------------------------------------------------------

@command("ameans", "ameans", byable=True)
def cmd_ameans(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"level": 1, "cluster": 7}) if p.options.strip() else {}
    if p.weight is not None:
        raise StataError(101, "weights not allowed")       # VERIFY: weights are not supported yet
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    mask = touse(s, p)
    level = _level(s, o)
    out = s.output
    out.write(f"\n{'Variable':>12} |{'Type':>10}{'Obs':>13}{'Mean':>14}"
              f"{'[' + _lvl_text(level) + '% Conf. Interval]':>27}\n", "text")
    out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
    last = None
    for v in vars_:
        x = v.data[mask & S.valid(v.data)]
        if len(x) == 0:
            raise StataError(2000, "no observations")
        res = B.three_means(x, level)
        last = res
        first = True
        for key, lab in (("a", "Arithmetic"), ("g", "Geometric"), ("h", "Harmonic")):
            if key not in res:
                continue
            n, m, se, lo, hi = res[key]
            out.write(f"{name12(v.name) if first else '':>12} |", "text")
            out.write(f"{lab:>11}{n:>8}{g9(m, '%9.0g'):>12}{g9(lo, '%9.0g'):>13}"
                      f"{g9(hi, '%9.0g'):>13}\n", "result")
            first = False
    out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
    if last is not None:
        pairs = []
        for key, sfx in (("a", ""), ("g", "_g"), ("h", "_h")):
            if key in last:
                n, m, se, lo, hi = last[key]
                pairs += [(f"N{sfx}", float(n)), (f"mean{sfx}", m), (f"lb{sfx}", lo),
                          (f"ub{sfx}", hi), (f"se{sfx}", se)]
        _set_r(s, pairs)


# ---------------------------------------------------------------------------
# icc
# ---------------------------------------------------------------------------

@command("icc", "icc")
def cmd_icc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"mixed": 1, "oneway": 1, "random": 1, "absolute": 1,
                                  "consistency": 1, "level": 1, "testvalue": 1, "format": 1}) \
        if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if len(names) not in (2, 3):
        raise StataError(102, "too few variables specified" if len(names) < 2 else "too many variables specified")
    if o.get("testvalue"):
        raise StataError(198, "option testvalue() is not available yet in OpenDTA")
    model = "oneway" if o.get("oneway") or len(names) == 2 else "mixed" if o.get("mixed") else "random"
    if len(names) == 2 and (o.get("mixed") or o.get("random") or o.get("absolute")
                            or o.get("consistency")):
        model = "oneway"
    kind = "absolute" if o.get("absolute") else "consistency"
    if model == "oneway":
        kind = "absolute"
    mask = touse(s, p)
    cols = [s.data.get(n) for n in names]
    for c in cols:
        mask &= S.valid(c.data)
    y = cols[0].data[mask]
    tg = cols[1].data[mask]
    tl, ti = np.unique(tg, return_inverse=True)
    n = len(tl)
    if len(names) == 3:
        rl, ri = np.unique(cols[2].data[mask], return_inverse=True)
        k = len(rl)
    else:
        cnt = np.bincount(ti)
        k = int(cnt.max())
        ri = np.zeros(len(y), dtype=int)
        for t in range(n):
            idx = np.flatnonzero(ti == t)
            ri[idx] = np.arange(len(idx))
        rl = np.arange(k)
    Y = np.full((n, k), np.nan)
    Y[ti, ri] = y
    if np.isnan(Y).any():
        good = ~np.isnan(Y).any(1)
        Y = Y[good]
        if Y.shape[0] < 2:
            raise StataError(2001, "insufficient observations")
        n = Y.shape[0]
    if n < 2 or k < 2:
        raise StataError(2001, "insufficient observations")
    level = _level(s, o)
    r = B.icc(Y, model, kind, level)
    out = s.output
    title = {"oneway": "One-way random-effects model",
             "random": "Two-way random-effects model",
             "mixed": "Two-way mixed-effects model"}[model]
    out.write("\nIntraclass correlations\n", "text")
    out.write(title + "\n", "text")
    out.write(("Absolute agreement" if kind == "absolute" else "Consistency of agreement") + "\n\n",
              "text")
    if model == "oneway":
        out.write(f"{'Random effects: ' + names[1]:<34}{'Number of targets':<19}= ", "text")
    elif model == "random":
        out.write(f"{'Random effects: ' + names[1]:<34}{'Number of targets':<19}= ", "text")
    else:
        out.write(f"{'Random effects: ' + names[1]:<34}{'Number of targets':<19}= ", "text")
    out.write(f"{n:>9}\n", "result")
    if len(names) == 3:
        lab = ("Random" if model == "random" else "Fixed") + " effects: " + names[2]
        out.write(f"{lab:<34}{'Number of raters':<19}= ", "text")
        out.write(f"{k:>9}\n", "result")
    else:
        out.write(f"{'':<34}{'Number of raters':<19}= ", "text")
        out.write(f"{k:>9}\n", "result")
    out.write("\n" + "-" * 62 + "\n", "text")
    out.write(f"{'':>23}|{'ICC':>12}{'[' + _lvl_text(level) + '% Conf. Interval]':>26}\n", "text")
    out.write("-" * 23 + "+" + "-" * 38 + "\n", "text")
    for lab, key in (("Individual", "ind"), ("Average", "avg")):
        out.write(f"{lab:>22} |", "text")
        out.write(f"{r[key]:>12.7f}{r[key + '_lb']:>13.7f}{r[key + '_ub']:>13.7f}\n", "result")
    out.write("-" * 62 + "\n", "text")
    out.write("F test that\n", "text")
    out.write(f"  ICC=0.00: F({r['df1']:.1f},{r['df2']:.1f}) = {r['F']:>6.2f}"
              f"{'Prob > F =':>16} {r['p']:.3f}\n\n", "text")
    out.write("Note: ICCs estimate correlations between individual measurements and between\n"
              "      average measurements made on the same target.\n", "text")
    _set_r(s, [("level", level), ("N_target", float(n)), ("N_rater", float(k)),
               ("icc_i", r["ind"]), ("icc_i_F", r["F"]), ("icc_i_df1", float(r["df1"])),
               ("icc_i_df2", float(r["df2"])), ("icc_i_p", r["p"]), ("icc_i_lb", r["ind_lb"]),
               ("icc_i_ub", r["ind_ub"]), ("icc_avg", r["avg"]), ("icc_avg_F", r["F"]),
               ("icc_avg_df1", float(r["df1"])), ("icc_avg_df2", float(r["df2"])),
               ("icc_avg_p", r["p"]), ("icc_avg_lb", r["avg_lb"]), ("icc_avg_ub", r["avg_ub"])])


# ---------------------------------------------------------------------------
# loneway
# ---------------------------------------------------------------------------

@command("loneway", "loneway")
def cmd_loneway(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"median": 1, "exact": 1, "level": 1}) if p.options.strip() else {}
    if o.get("median") or o.get("exact"):
        raise StataError(198, f"option {'median' if o.get('median') else 'exact'} is not available yet in OpenDTA")
    names = s.expand_varlist(p.varlist)
    if len(names) != 2:
        raise StataError(102, "too few variables specified" if len(names) < 2 else "too many variables specified")
    yv, gv = s.data.get(names[0]), s.data.get(names[1])
    mask = touse(s, p) & S.valid(yv.data)
    if gv.is_string:
        gcol = np.asarray(gv.raw, dtype=object)
        mask &= np.array([x != "" for x in gcol])
        codes = np.unique(gcol[mask], return_inverse=True)[1]
    else:
        mask &= S.valid(gv.data)
        codes = np.unique(gv.data[mask], return_inverse=True)[1]
    if mask.sum() == 0:
        raise StataError(2000, "no observations")
    level = _level(s, o)
    r = B.loneway(yv.data[mask], codes, level)
    if r["G"] < 2:
        raise StataError(198, "need at least 2 groups")        # VERIFY
    out = s.output
    head = f"One-way Analysis of Variance for {yv.name}:"
    if yv.label:
        head += f" {yv.label}"
    out.write("\n" + head + "\n", "text")
    out.write(f"{'':>48}{'Number of obs':<14}= ", "text")
    out.write(f"{r['N']:>9}\n", "result")
    out.write(f"{'':>48}{'R-squared':<14}= ", "text")
    out.write(f"{r['r2']:>9.4f}\n\n", "result")
    out.write(f"{'Source':<20}{'SS':>12}{'df':>8}{'MS':>14}{'F':>10}{'Prob > F':>10}\n", "text")
    out.write("-" * 74 + "\n", "text")
    out.write(f"{'Between ' + gv.name:<20}", "text")
    out.write(f"{r['ssb']:>12.7g}{r['dfb']:>8}{r['msb']:>14.7g}{r['F']:>10.2f}{r['p']:>10.4f}\n", "result")
    out.write(f"{'Within ' + gv.name:<20}", "text")
    out.write(f"{r['ssw']:>12.7g}{r['dfw']:>8}{r['msw']:>14.7g}\n", "result")
    out.write("-" * 74 + "\n", "text")
    out.write(f"{'Total':<20}", "text")
    out.write(f"{r['sst']:>12.7g}{r['N'] - 1:>8}{r['sst'] / (r['N'] - 1):>14.7g}\n\n", "result")
    out.write(f"{'':>9}{'Intraclass':>12}{'Asy.':>12}\n", "text")
    out.write(f"{'':>9}{'correlation':>12}{'S.E.':>12}{'[' + _lvl_text(level) + '% Conf. Interval]':>28}\n", "text")
    out.write(f"{'':>9}" + "-" * 48 + "\n", "text")
    out.write(f"{'':>9}{r['rho']:>12.5f}{r['se']:>12.5f}{r['lb']:>14.5f}{r['ub']:>14.5f}\n\n", "result")
    out.write(f"{'':>9}{'Estimated SD of ' + gv.name + ' effect':<39}", "text")
    out.write(f"{r['sd_b']:>12.7g}\n", "result")
    out.write(f"{'':>9}{'Estimated SD within ' + gv.name:<39}", "text")
    out.write(f"{r['sd_w']:>12.7g}\n", "result")
    out.write(f"{'':>9}{'Estimated SD of ' + yv.name:<39}", "text")
    out.write(f"{r['sd']:>12.7g}\n", "result")
    out.write(f"{'':>9}{'Est. reliability of a ' + gv.name + ' mean':<39}", "text")
    out.write(f"{r['rel']:>12.5f}\n", "result")
    out.write(f"{'':>9}{'(evaluated at n=' + f'{r['n0']:.2f}' + ')':>40}\n", "text")
    _set_r(s, [("N", float(r["N"])), ("rho", r["rho"]), ("lb", r["lb"]), ("ub", r["ub"]),
               ("F", r["F"]), ("df_m", float(r["dfb"])), ("df_r", float(r["dfw"])),
               ("sd_w", r["sd_w"]), ("sd_b", r["sd_b"]), ("rho_se", r["se"]), ("N_g", float(r["G"])),
               ("Var_w", r["msw"]), ("Var_b", r["sd_b"] ** 2)])


# ---------------------------------------------------------------------------
# esize, esizei
# ---------------------------------------------------------------------------

_ES_ROWS = (("cohen", "Cohen's d"), ("hedges", "Hedges's g"), ("glass1", "Glass's Delta 1"),
            ("glass2", "Glass's Delta 2"), ("pbcorr", "Point-Biserial r"))


def _es_options(o: dict) -> list[str]:
    keys = [k for k in ("cohensd", "hedgesg", "glassdelta", "pbcorr") if o.get(k)]
    if o.get("all"):
        keys = ["cohensd", "hedgesg", "glassdelta", "pbcorr"]
    if not keys:
        keys = ["cohensd", "hedgesg"]
    want: list[str] = []
    for k in keys:
        want += {"cohensd": ["cohen"], "hedgesg": ["hedges"], "glassdelta": ["glass1", "glass2"],
                 "pbcorr": ["pbcorr"]}[k]
    return want


def _es_print(s: "Session", n1: int, m1: float, s1: float, n2: int, m2: float, s2: float,
              o: dict, lab1: str, lab2: str) -> None:
    level = _level(s, o)
    if o.get("unequal") or o.get("welch"):
        raise StataError(198, "options unequal and welch are not available yet in OpenDTA")
    r = B.effect_sizes(n1, m1, s1, n2, m2, s2, level)
    want = _es_options(o)
    out = s.output
    out.write("\nEffect size based on mean comparison\n\n", "text")
    out.write(f"{'Obs per group:':>31}\n", "text")
    out.write(f"{lab1:>31} = {n1:>9}\n{lab2:>31} = {n2:>9}\n", "text")
    out.write("-" * 17 + "+" + "-" * 50 + "\n", "text")
    out.write(f"{'Effect Size':>16} |{'Estimate':>11}{'[' + _lvl_text(level) + '% Conf. Interval]':>34}\n", "text")
    out.write("-" * 17 + "+" + "-" * 50 + "\n", "text")
    for key, lab in _ES_ROWS:
        if key in want:
            est, lo, hi = r[key]
            out.write(f"{lab:>16} |", "text")
            out.write(f"{est:>11.7f}{lo:>14.7f}{hi:>14.7f}\n", "result")
    out.write("-" * 17 + "+" + "-" * 50 + "\n", "text")
    pairs = [("level", level), ("df", float(r["df"]))]
    names = {"cohen": "d", "hedges": "g", "glass1": "delta1", "glass2": "delta2", "pbcorr": "r_pb"}
    for key in want:
        est, lo, hi = r[key]
        pairs += [(names[key], est), (f"lb_{names[key]}", lo), (f"ub_{names[key]}", hi)]
    _set_r(s, pairs)


@command("esize", "esize")
def cmd_esize(s: "Session", args: str) -> None:
    a = args.lstrip()
    if not a.startswith("twosample"):
        raise StataError(199, "unrecognized command: esize")
    p = parse_standard(a[len("twosample"):])
    o = match_options(p.options, {"by": 2, "cohensd": 2, "hedgesg": 2, "glassdelta": 2,
                                  "pbcorr": 2, "all": 3, "level": 1, "unequal": 3, "welch": 1}) \
        if p.options.strip() else {}
    if not o.get("by"):
        raise StataError(198, "option by() required")
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    from .nptests import _groups
    _, grp = _groups(s, o["by"], mask)
    if len(grp) != 2:
        raise StataError(420, "group variable must have exactly 2 distinct values")    # VERIFY
    (g1, t1, m1), (g2, t2, m2) = grp
    x1, x2 = v.data[m1], v.data[m2]
    if len(x1) < 2 or len(x2) < 2:
        raise StataError(2001, "insufficient observations")
    _es_print(s, len(x1), float(x1.mean()), float(x1.std(ddof=1)), len(x2), float(x2.mean()),
              float(x2.std(ddof=1)), o, t1, t2)


@command("esizei", "esizei")
def cmd_esizei(s: "Session", args: str) -> None:
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    a = body.split()
    if len(a) != 6:
        raise StataError(198, "invalid syntax")
    o = match_options(opts, {"cohensd": 2, "hedgesg": 2, "glassdelta": 2, "pbcorr": 2, "all": 3,
                             "level": 1, "unequal": 3, "welch": 1}) if opts.strip() else {}
    n1, m1, s1, n2, m2, s2 = (float(t) for t in a)
    _es_print(s, int(n1), m1, s1, int(n2), m2, s2, o, "Group 1", "Group 2")


# ---------------------------------------------------------------------------
# lnskew0
# ---------------------------------------------------------------------------

@command("lnskew0", "lnskew0")
def cmd_lnskew0(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"level": 1, "zero": 1}) if p.options.strip() else {}
    new = p.varlist.split()
    if len(new) != 1 or p.exp is None:
        raise StataError(198, "invalid syntax")
    name = new[0]
    if s.data.has(name):
        raise StataError(110, f"{name} already defined")
    mask = touse(s, p)
    val = np.broadcast_to(np.asarray(eval_vector(s, p.exp), dtype=float), (s.data.nobs,))
    mask &= S.valid(val)
    x = val[mask]
    if len(x) < 4:
        raise StataError(2001, "insufficient observations")
    level = _level(s, o)
    k, lo, hi, sign = B.lnskew(x, None, level)
    if math.isnan(k):
        raise StataError(198, "could not find a solution")      # VERIFY
    col = np.full(s.data.nobs, SYS)
    col[mask] = np.log(x - k) if sign == 1 else np.log(k - x)
    s.data.add(Variable(name, "float", col))
    skew = B.skewness(col[mask])
    skew = 0.0 if abs(skew) < 5e-7 else skew
    out = s.output
    lab = f"ln({p.exp.strip()}-k)" if sign == 1 else f"ln(k-{p.exp.strip()})"
    out.write(f"\n{'Transform':>16} |{'k':>10}{'[' + _lvl_text(level) + '% Conf. Interval]':>30}{'Skewness':>12}\n", "text")
    out.write("-" * 17 + "+" + "-" * 51 + "\n", "text")
    out.write(f"{lab:>16} |", "text")
    out.write(f"{k:>10.5g}{lo:>14.5g}{hi:>14.5g}{skew:>12.6f}\n", "result")
    _set_r(s, [("gamma", k), ("lb", lo), ("ub", hi), ("skewness", skew), ("level", level)])
    s.r["N"] = float(len(x))


# ---------------------------------------------------------------------------
# ladder
# ---------------------------------------------------------------------------

_LADDER = (("cubic", "x^3", lambda x: x ** 3), ("square", "x^2", lambda x: x ** 2),
           ("identity", "x", lambda x: x), ("square root", "sqrt(x)", np.sqrt),
           ("log", "log(x)", np.log), ("1/(square root)", "1/sqrt(x)", lambda x: 1 / np.sqrt(x)),
           ("inverse", "1/x", lambda x: 1 / x), ("1/square", "1/(x^2)", lambda x: 1 / x ** 2),
           ("1/cubic", "1/(x^3)", lambda x: 1 / x ** 3))


def ladder_table(x: np.ndarray, adjust: bool = True):
    """(label, formula, values or None, chi2, p) for each rung."""
    rows = []
    n = len(x)
    for lab, form, f in _LADDER:
        with np.errstate(all="ignore"):
            t = f(x)
        if not np.isfinite(t).all():
            rows.append((lab, form, None, math.nan, math.nan))
            continue
        z1, z2 = _sk_z(t)
        k2 = z1 * z1 + z2 * z2
        if adjust:
            chi2, pv = _royston_adjust(k2, n)
        else:
            chi2, pv = k2, float(st.chi2.sf(k2, 2))
        rows.append((lab, form, t, chi2, pv))
    return rows


@command("ladder", "ladder", byable=True)
def cmd_ladder(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"generate": 1, "noadjust": 3, "histogram": 4}) \
        if p.options.strip() else {}
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    x = v.data[mask]
    if len(x) < 8:
        raise StataError(2001, "insufficient observations")
    rows = ladder_table(x, not o.get("noadjust"))
    out = s.output
    out.write(f"\n{'Transformation':<22}{'formula':<14}{'chi2(2)':>12}{'P(chi2)':>12}\n", "text")
    out.write("-" * 60 + "\n", "text")
    for lab, form, t, chi2, pv in rows:
        out.write(f"{lab:<22}{form:<14}", "text")
        if t is None or not math.isfinite(chi2):
            out.write(f"{'.':>12}{'.' if t is None else f'{pv:.3f}':>12}\n", "result")
        else:
            out.write(f"{chi2:>12.2f}{pv:>12.3f}\n", "result")
    if o.get("generate"):
        new = o["generate"]
        if s.data.has(new):
            raise StataError(110, f"{new} already defined")
        ok = [(chi2 if math.isfinite(chi2) else math.inf, i) for i, (_l, _f, t, chi2, _p) in enumerate(rows)
              if t is not None]
        best = min(ok)[1]
        col = np.full(s.data.nobs, SYS)
        col[mask] = rows[best][2]
        s.data.add(Variable(new, "float", col))


# ---------------------------------------------------------------------------
# boxcox
# ---------------------------------------------------------------------------

@command("boxcox", "boxcox")
def cmd_boxcox(s: "Session", args: str) -> None:
    from ..estimation.results import CoefRow, Estimates, coef_table, post
    p = parse_standard(args)
    o = match_options(p.options, {"model": 1, "notrans": 3, "lrtest": 2, "level": 1, "nolog": 3,
                                  "iterate": 2}) if p.options.strip() else {}
    mdl = str(o.get("model", "lhsonly"))
    if mdl != "lhsonly" or o.get("notrans"):
        raise StataError(198, "only model(lhsonly) without notrans() is available in OpenDTA")  # VERIFY
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    mask = touse(s, p)
    for v in vars_:
        mask &= S.valid(v.data)
    y = vars_[0].data[mask]
    if (y <= 0).any():
        raise StataError(198, "dependent variable must be positive")
    n = len(y)
    X = np.column_stack([v.data[mask] for v in vars_[1:]] + [np.ones(n)])
    fit = B.boxcox_fit(y, X)
    th, se, ll = fit["theta"], fit["se"], fit["ll"]
    level = _level(s, o)
    out = s.output
    ll1 = fit["ll_at"](1.0)
    lr = 2 * (ll - ll1)
    out.write("\n", "text")
    out.write(f"{'':>49}{'Number of obs':<16}= ", "text")
    out.write(f"{n:>9}\n", "result")
    out.write(f"{'':>49}{'LR chi2(1)':<16}= ", "text")
    out.write(f"{lr:>9.2f}\n", "result")
    out.write(f"{'Log likelihood = ' + f'{ll:.7g}':<49}{'Prob > chi2':<16}= ", "text")
    out.write(f"{float(st.chi2.sf(lr, 1)):>9.3f}\n\n", "result")
    rows = [CoefRow("coef", "/theta", th, se, name="theta")]
    coef_table(s, vars_[0].name, rows, stat="z", df=None, level=level)
    if len(vars_) > 1:
        out.write("\nEstimates of scale-variant parameters\n", "text")
        out.write("-" * 13 + "+" + "-" * 12 + "\n", "text")
        out.write(f"{'':>12} |{'Coef.':>11}\n", "text")
        out.write("-" * 13 + "+" + "-" * 12 + "\n", "text")
        out.write(f"{'Notrans':<12} |\n", "text")
        names = [v.name for v in vars_[1:]] + ["_cons"]
        for nm, b in zip(names, fit["b"]):
            out.write(f"{name12(nm):>12} |", "text")
            out.write(f"{b:>11.7g}\n", "result")
        out.write("-" * 13 + "+" + "-" * 12 + "\n", "text")
        out.write(f"{'/sigma':>12} |", "text")
        out.write(f"{fit['sigma']:>11.7g}\n", "result")
        out.write("-" * 13 + "+" + "-" * 12 + "\n", "text")
    out.write("\n" + "-" * 58 + "\n", "text")
    out.write(f"{'Test':<12}{'Restricted':>16}{'LR statistic':>16}{'P-value':>14}\n", "text")
    out.write(f"{'H0:':<12}{'log likelihood':>16}{'chi2':>16}{'Prob > chi2':>14}\n", "text")
    out.write("-" * 58 + "\n", "text")
    tests = []
    for t0 in (-1, 0, 1):
        llr = fit["ll_at"](float(t0))
        chi = 2 * (ll - llr)
        pv = float(st.chi2.sf(chi, 1))
        out.write(f"{'theta = ' + str(t0):<12}", "text")
        out.write(f"{llr:>16.7g}{chi:>16.2f}{pv:>14.3f}\n", "result")
        tests.append((t0, llr, chi, pv))
    out.write("-" * 58 + "\n", "text")
    sample = mask.copy()
    est = Estimates(cmd="boxcox", depvar=vars_[0].name, names=["/theta"], b=np.array([th]),
                    V=np.array([[se * se]]), N=n, stat="z", rows=rows)
    est.level = level
    est.title = "Box-Cox regression"
    post(s, est, [("N", float(n)), ("ll", ll), ("chi2", lr), ("theta", th)],
         [("cmd", "boxcox"), ("depvar", vars_[0].name), ("model", "lhsonly")], sample)
