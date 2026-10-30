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
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    mask = touse(s, p)
    from .summarize import weights
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight"))
    level = _level(s, o)
    out = s.output
    out.write(f"\n{'Variable':>12} |    Type             Obs        Mean"
              f"       [{_lvl_text(level)}% Conf. Interval]\n", "text")
    sep = "-" * 13 + "+" + "-" * 63 + "\n"
    out.write(sep, "text")
    last = None
    for i, v in enumerate(vars_):
        ok = mask & S.valid(v.data)
        x = v.data[ok]
        if len(x) == 0:
            raise StataError(2000, "no observations")
        res = B.three_means(x, level, w[ok] if w is not None else None, wtype or "")
        last = res
        first = True
        if i > 0:
            out.write(sep, "text")
        for key, lab in (("a", "Arithmetic"), ("g", "Geometric"), ("h", "Harmonic")):
            if key not in res:
                continue
            n, m, se, lo, hi = res[key]
            out.write(f"{name12(v.name) if first else '':>12} |", "text")
            out.write(f"{lab:>11}{n:>13}{g9(m, '%9.0g'):>12}{g9(lo, '%9.0g'):>16}"
                      f"{g9(hi, '%9.0g'):>11} \n", "result")
            first = False
    out.write(sep if len(vars_) > 1 else "-" * 77 + "\n", "text")      # VERIFY: closing rule
    if last is not None:
        # Stata stores the variance of x, of ln x and of 1/x (N times the squared standard error
        # of their means) and the number of positive values, listed from ub_h back to N (compat
        # 1913)
        pairs = []
        for key, sfx in (("a", ""), ("g", "_g"), ("h", "_h")):
            if key in last:
                n, m, se, lo, hi = last[key]
                scale = 1.0 if key == "a" else m if key == "g" else m * m
                var = (se / scale) ** 2 * n
                head = [("N", float(n))] if key == "a" else [("N_pos", float(n))] if key == "g" else []
                pairs = [(f"ub{sfx}", hi), (f"lb{sfx}", lo), (f"Var{sfx}", var), (f"mean{sfx}", m)] + \
                    list(reversed(head)) + pairs
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
    r0 = 0.0
    if o.get("testvalue") not in (None, True):
        try:
            r0 = float(s.eval(str(o["testvalue"])))
        except (StataError, TypeError, ValueError):
            raise StataError(198, "option testvalue() incorrectly specified")
        if not 0 <= r0 < 1:
            raise StataError(198, "testvalue() must be in [0, 1)")                  # VERIFY
    model = "oneway" if o.get("oneway") or len(names) == 2 else "mixed" if o.get("mixed") else "random"
    if len(names) == 2 and (o.get("mixed") or o.get("random") or o.get("absolute")
                            or o.get("consistency")):
        model = "oneway"
    # VERIFY: the default of the two-way mixed model is consistency, the other models absolute
    kind = "consistency" if (model == "mixed" and not o.get("absolute")) or o.get("consistency") \
        else "absolute"
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
        # targets with fewer ratings are left out, with a note (compat 1913)
        good = ~np.isnan(Y).any(1)
        dropped = int((~good).sum())
        why = "of unbalanced data" if len(names) == 2 else "not rated by all raters"
        s.output.write(f"({dropped} target{'s' if dropped > 1 else ''} omitted from computation because "
                       f"{why})\n", "text")
        Y = Y[good]
        if Y.shape[0] < 2:
            raise StataError(2001, "insufficient observations")
        n = Y.shape[0]
    if n < 2 or k < 2:
        raise StataError(2001, "insufficient observations")
    level = _level(s, o)
    r = B.icc(Y, model, kind, level)
    tests = B.icc_test(r, model, kind, r0)
    out = s.output
    title = {"oneway": "One-way random-effects model",
             "random": "Two-way random-effects model",
             "mixed": "Two-way mixed-effects model"}[model]
    out.write("\nIntraclass correlations\n", "text")
    out.write(title + "\n", "text")
    out.write(("Absolute agreement" if kind == "absolute" else "Consistency of agreement") + "\n\n",
              "text")
    out.write(f"{'Random effects: ' + names[1]:<33}{'Number of targets':<18}= ", "text")
    out.write(f"{n:>9}\n", "result")
    if len(names) == 3:
        lab = f"{'Random' if model == 'random' else 'Fixed':>6} effects: " + names[2]
        out.write(f"{lab:<33}{'Number of raters':<18}= ", "text")
    else:
        out.write(f"{'':<33}{'Number of raters':<18}= ", "text")
    out.write(f"{k:>9}\n", "result")
    out.write("\n" + "-" * 62 + "\n", "text")
    out.write(f"{names[0]:>22} |{'ICC':>11}{'[' + _lvl_text(level) + '% Conf. Interval]':>27}\n", "text")
    out.write("-" * 23 + "+" + "-" * 38 + "\n", "text")
    for lab, key in (("Individual", "ind"), ("Average", "avg")):
        out.write(f"{lab:>22} |", "text")
        out.write(f"{g9(r[key], '%9.0g'):>11}{g9(r[key + '_lb'], '%9.0g'):>15}"
                  f"{g9(r[key + '_ub'], '%9.0g'):>12}\n", "result")
    out.write("-" * 62 + "\n", "text")
    out.write("F test that\n", "text")
    # one line per ICC when the two tests differ: ICC(1) and ICC(k) (compat 1913)
    for lab, key in ((("ICC", "i"),) if r0 == 0 else (("ICC(1)", "i"), ("ICC(k)", "avg"))):
        F, df1, df2, pv = tests[key]
        out.write(f"  {lab}={r0:.2f}: F({df1:.1f}, {df2:.1f}) = {F:.2f}".ljust(46)
                  + f"Prob > F = {pv:.3f}\n", "text")
    out.write("\n", "text")
    if not (model == "mixed" and kind == "absolute"):      # no note there (compat 1913)
        out.write("Note: ICCs estimate correlations between individual measurements\n"
                  "      and between average measurements made on the same target.\n", "text")
    macros = [("model", {"oneway": "one-way random effects", "random": "two-way random effects",
                         "mixed": "two-way mixed effects"}[model]),
              ("depvar", names[0]), ("target", names[1])]
    if len(names) == 3:
        macros.append(("rater", names[2]))
    macros.append(("type", kind))                       # VERIFY: macros of the two-way models
    ti_, ta_ = tests["i"], tests["avg"]
    _set_r(s, [("N_target", float(n)), ("N_rater", float(k)),
               ("icc_i", r["ind"]), ("icc_i_F", ti_[0]), ("icc_i_df1", ti_[1]),
               ("icc_i_df2", ti_[2]), ("icc_i_p", ti_[3]), ("icc_i_lb", r["ind_lb"]),
               ("icc_i_ub", r["ind_ub"]), ("icc_avg", r["avg"]), ("icc_avg_F", ta_[0]),
               ("icc_avg_df1", ta_[1]), ("icc_avg_df2", ta_[2]),
               ("icc_avg_p", ta_[3]), ("icc_avg_lb", r["avg_lb"]), ("icc_avg_ub", r["avg_ub"]),
               ("testvalue", r0), ("level", level)] + macros)


# ---------------------------------------------------------------------------
# loneway
# ---------------------------------------------------------------------------

@command("loneway", "loneway")
def cmd_loneway(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"mean": 1, "median": 1, "exact": 1, "level": 1}) if p.options.strip() else {}
    if o.get("mean") and o.get("median"):
        raise StataError(198, "options mean and median may not be combined")     # VERIFY
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
    if o.get("exact") and len(set(np.bincount(codes).tolist())) > 1:
        raise StataError(198, "Balanced data (equal groups) required for exact confidence interval.")  # VERIFY: rc
    r = B.loneway(yv.data[mask], codes, level, center="mean" if o.get("mean") else "median" if o.get("median") else "",
                  exact=bool(o.get("exact")))
    if r["G"] < 2:
        raise StataError(198, "need at least 2 groups")        # VERIFY
    out = s.output
    head = f"One-way Analysis of Variance for {yv.name}: " + (yv.label or "")
    out.write("\n" + " " * ((80 - len(head)) // 2) + head + "\n\n", "text")
    out.write(f"{'':>45}{'Number of obs':<14}= ", "text")
    out.write(f"{r['N']:>12}\n", "result")
    out.write(f"{'':>45}{'R-squared':<14}= ".replace("R-squared    ", "    R-squared"), "text")
    out.write(f"{r['r2']:>12.4f}\n\n", "result")
    out.write(f"{'Source':>10}{'SS':>18}{'df':>11}{'MS':>8}{'F':>13}{'Prob > F':>13}\n", "text")
    out.write("-" * 73 + "\n", "text")
    out.write(f"{'Between ' + gv.name:<14}", "text")
    out.write(f"{g9(r['ssb'], '%10.0g'):>18}{r['dfb']:>7}{g9(r['msb'], '%10.0g'):>13}"
              f"{r['F']:>10.2f}{r['p']:>11.4f}\n", "result")
    out.write(f"{'Within ' + gv.name:<14}", "text")
    out.write(f"{g9(r['ssw'], '%10.0g'):>18}{r['dfw']:>7}{g9(r['msw'], '%10.0g'):>13}\n", "result")
    out.write("-" * 73 + "\n", "text")
    out.write(f"{'Total':<14}", "text")
    out.write(f"{g9(r['sst'], '%10.0g'):>18}{r['N'] - 1:>7}"
              f"{g9(r['sst'] / (r['N'] - 1), '%10.0g'):>13}\n\n", "result")
    trunc = r["rho_raw"] < 0
    if o.get("mean") or o.get("median"):
        # with mean or median there is no standard error nor interval (compat 1913)
        r["se"] = r["lb"] = r["ub"] = math.nan
        # the reliability follows the correlation shown (compat 1913)
        n0 = r["n0"]
        r["rel"] = n0 * r["rho"] / (1 + (n0 - 1) * r["rho"])
    exact = "        ------ Exact -----" if o.get("exact") else ""
    out.write(f"{'':>9}{'Intraclass':>10}{'Asy.':>11}{exact}\n", "text")
    out.write(f"{'':>9}{'correlation':>11}{'S.E.':>10}{'':>7}[{_lvl_text(level)}% Conf. Interval]\n", "text")
    out.write(f"{'':>9}" + "-" * 48 + "\n", "text")
    def f5(x: float, w: int) -> str:
        return f"{'.':>{w}}" if math.isnan(x) else f"{x:>{w}.5f}"
    out.write(f"{'':>9}{r['rho']:>10.5f}" + ("*" if trunc else " ")
              + f"{f5(r['se'], 11)}{f5(r['lb'], 14)}{f5(r['ub'], 12)}\n\n", "result")
    rel_trunc = r["rel_raw"] < 0
    for lab, val, star in (("Estimated SD of " + gv.name + " effect", r["sd_b"], False),
                           ("Estimated SD within " + gv.name, r["sd_w"], False),
                           ("Est. reliability of a " + gv.name + " mean", r["rel"], rel_trunc)):
        txt = "." if math.isnan(val) else (f"{val:.5f}" if lab.startswith("Est. rel") else g9(val, "%9.0g"))
        out.write(f"{'':>9}{lab:<39}", "text")
        out.write(f"{txt:>9}" + ("*" if star else "") + "\n", "result")
    out.write(f"{'':>14}(evaluated at n={r['n0']:.2f})\n", "text")
    if trunc or rel_trunc:
        out.write("\n(*) Truncated at zero.\n", "text")
    # in the order of return list (compat 1913); missing for mean and median (VERIFY)
    nan = lambda x: SYS if math.isnan(x) else x                    # noqa: E731
    _set_r(s, [("lb", nan(r["lb"])), ("ub", nan(r["ub"])), ("se", nan(r["se"])), ("sd_w", r["sd_w"]),
               ("sd_b", nan(r["sd_b"])), ("rho_t", r["rel"]), ("rho", r["rho"]), ("N", float(r["N"]))])


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
              o: dict, lab1: str, lab2: str, empty: tuple = ()) -> None:
    """empty: labels of other values of the by() variable, listed with 0 observations."""
    level = _level(s, o)
    r = B.effect_sizes(n1, m1, s1, n2, m2, s2, level, unequal=bool(o.get("unequal")), welch=bool(o.get("welch")))
    want = _es_options(o)
    out = s.output
    uneq = bool(o.get("unequal") or o.get("welch"))
    out.write("\nEffect size based on mean comparison" + (", unequal variances" if uneq else "") + "\n\n", "text")
    out.write(f"{'Obs per group:':>45}\n", "text")
    out.write(f"{lab1:>44} = {n1:>10}\n{lab2:>44} = {n2:>10}\n", "text")
    for lab in empty:
        out.write(f"{lab:>44} = {0:>10}\n", "text")
    out.write("-" * 57 + "\n", "text")
    out.write(f"{'Effect Size':>19} |{'Estimate':>11}{'[' + _lvl_text(level) + '% Conf. Interval]':>25}\n", "text")
    out.write("-" * 20 + "+" + "-" * 36 + "\n", "text")
    for key, lab in _ES_ROWS:
        if key in want:
            est, lo, hi = r[key]
            out.write(f"{lab:>19} |", "text")
            out.write(f"{g9(est):>11}{g9(lo):>13}{g9(hi):>12}\n", "result")
    out.write("-" * 57 + "\n", "text")
    if uneq:
        lab = ("Welch's" if o.get("welch") else "Satterthwaite's") + " degrees of freedom"
        out.write(f"{lab:>46} = ", "text")
        out.write(f"{r['df']:8.4f}\n", "result")
    pairs = []
    names = {"cohen": "d", "hedges": "g", "glass1": "delta1", "glass2": "delta2", "pbcorr": "r_pb"}
    for key in want:
        est, lo, hi = r[key]
        pairs += [(names[key], est), (f"lb_{names[key]}", lo), (f"ub_{names[key]}", hi)]
    if uneq:
        pairs.append(("df_t", r["df"]))
    pairs += [("N_1", float(n1)), ("N_2", float(n2)), ("level", level)]
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
    by = s.expand_varlist(o["by"])[0]
    lab = lambda g: f"{by}=={g:g}" if not isinstance(g, str) else f"{by}=={g}"      # noqa: E731
    # the other values of the by() variable in the data are listed with 0 observations (compat
    # 1913: esize ... if j < 3, by(j) lists j==3 and j==4; VERIFY)
    bv = s.data.get(by)
    if bv.is_string:
        others = sorted({x for x in bv.raw if x != ""} - {g1, g2})
    else:
        others = sorted({float(x) for x in bv.data[S.valid(bv.data)]} - {float(g1), float(g2)})
    _es_print(s, len(x1), float(x1.mean()), float(x1.std(ddof=1)), len(x2), float(x2.mean()),
              float(x2.std(ddof=1)), o, lab(g1), lab(g2), tuple(lab(g) for g in others))


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
    out = s.output
    ex = p.exp.strip()
    # VERIFY: for a negatively skewed variable Stata shows ln(-x-k) with k and its limits negated
    if sign == 1:
        lab, kd, lod, hid = f"ln({ex}-k)", k, lo, hi
        lot = "-infinity" if math.isnan(lo) else g9(lo)
        hit = "infinity" if math.isnan(hi) else g9(hi)
    else:
        lab, kd, lod, hid = f"ln(-{ex}-k)", -k, -hi, -lo
        lot = "-infinity" if math.isnan(hi) else g9(-hi)
        hit = "infinity" if math.isnan(lo) else g9(-lo)
    out.write(f"\n{'Transform':>16} |{'k':>10}{'[' + _lvl_text(level) + '% Conf. Interval]':>25}{'Skewness':>15}\n", "text")
    out.write("-" * 17 + "+" + "-" * 50 + "\n", "text")
    out.write(f"{lab:>16} |", "text")
    if math.isnan(lo) and math.isnan(hi):
        out.write(f"{g9(kd):>11}{'(not calculated)':>22}{g9(skew):>17}\n", "result")
    else:
        out.write(f"{g9(kd):>11}{lot:>13}{hit:>11}{g9(skew):>15}\n", "result")
    _set_r(s, [("gamma", kd), ("lb", lod), ("ub", hid), ("skewness", skew), ("level", level)])
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
    out.write(f"\n{'Transformation':<23}{'formula':<7}{'chi2(2)':>22}{'P(chi2)':>14}\n", "text")
    out.write("-" * 66 + "\n", "text")
    for lab, form, t, chi2, pv in rows:
        out.write(f"{lab:<23}{form.replace('x', v.name):<13}", "text")
        pt = "." if math.isnan(pv) else f"{pv:.3f}"
        if t is None or not math.isfinite(chi2):
            out.write(f"{'.':>15}{pt:>13}\n", "result")
        else:
            out.write(f"{chi2:>15.2f}{pt:>13}\n", "result")
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
        out.write(f"\n({new} = {rows[best][1].replace('x', v.name)}   generated)\n", "text")


# ---------------------------------------------------------------------------
# boxcox
# ---------------------------------------------------------------------------

@command("boxcox", "boxcox")
def cmd_boxcox(s: "Session", args: str) -> None:
    from ..estimation.results import CoefRow, Estimates, coef_table, post
    p = parse_standard(args)
    o = match_options(p.options, {"model": 1, "notrans": 3, "lrtest": 2, "level": 1, "nolog": 3,
                                  "iterate": 2}) if p.options.strip() else {}
    mdl = str(o.get("model", "lhsonly")).strip()
    if mdl not in ("lhsonly", "rhsonly", "lambda", "theta"):
        raise StataError(198, f"model({mdl}) invalid")                              # VERIFY
    if mdl != "lhsonly":
        _boxcox_general(s, p, o, mdl, args)
        return
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    if o.get("notrans") not in (None, True):
        # the notrans() variables come first among the coefficients (compat 1913)
        extra = [v for v in numeric_vars(s, str(o["notrans"]), allow_string=False) if v not in vars_]
        vars_ = vars_[:1] + extra + vars_[1:]
    mask = touse(s, p)
    for v in vars_:
        mask &= S.valid(v.data)
    y = vars_[0].data[mask]
    if (y <= 0).any():
        raise StataError(198, "dependent variable must be positive")
    n = len(y)
    X = np.column_stack([v.data[mask] for v in vars_[1:]] + [np.ones(n)])
    fit = B.boxcox_fit(y, X)
    fit0 = B.boxcox_fit(y, np.ones((n, 1)))            # comparison model: constant only
    th, se, ll = fit["theta"], fit["se"], fit["ll"]
    level = _level(s, o)
    out = s.output
    ll0 = fit0["ll"]
    lr = 2 * (ll - ll0)
    ncov = len(vars_) - 1
    out.write("Fitting comparison model\n\nFitting full model\n\n", "text")
    out.write(f"{'':>50}{'Number of obs':<16}= ", "text")
    out.write(f"{n:>10}\n", "result")
    out.write(f"{'':>50}{'LR chi2(' + str(ncov) + ')':<16}= ", "text")
    out.write(f"{lr:>10.2f}\n", "result")
    out.write(f"{'Log likelihood = ' + f'{ll:.8g}':<50}{'Prob > chi2':<16}= ", "text")
    out.write((f"{float(st.chi2.sf(lr, ncov)):>10.3f}" if ncov > 0 else f"{'.':>10}") + "\n", "result")
    out.write(" \n", "text")
    rows = [CoefRow("coef", "/theta", th, se, name="theta")]
    coef_table(s, vars_[0].name, rows, stat="z", df=None, level=level)
    out.write(" \nEstimates of scale-variant parameters\n", "text")
    out.write("-" * 28 + "\n", "text")
    out.write(f"{'':>13}|{'Coef.':>11}\n", "text")
    sep = "-" * 13 + "+" + "-" * 14 + "\n"
    out.write(sep, "text")
    out.write(f"{'Notrans':<13}|\n", "text")
    names = [v.name for v in vars_[1:]] + ["_cons"]
    for nm, b in zip(names, fit["b"]):
        out.write(f"{name12(nm):>12} |", "text")
        out.write(f"{g9(b):>11}\n", "result")
    out.write(sep, "text")
    out.write(f"{'/sigma':>12} |", "text")
    out.write(f"{g9(fit['sigma']):>11}\n", "result")
    out.write("-" * 28 + "\n\n", "text")
    out.write("-" * 57 + "\n", "text")
    out.write("   Test         Restricted     LR statistic      P-value\n", "text")
    out.write("    H0:       log likelihood       chi2       Prob > chi2\n", "text")
    out.write("-" * 57 + "\n", "text")
    tests = {}
    for t0, key in ((-1, "tm1"), (0, "t0"), (1, "t1")):
        llr = fit["ll_at"](float(t0))
        chi = 2 * (ll - llr)
        pv = float(st.chi2.sf(chi, 1))
        out.write(f"{'theta = ' + format(t0, '2d'):<10}", "text")
        out.write(f"{g9(llr, '%10.0g'):>16}{chi:>13.2f}{pv:>16.3f}\n", "result")
        tests[key] = (llr, chi, pv)
    out.write("-" * 57 + "\n", "text")
    sample = mask.copy()
    k = len(names) + 2
    enames = ["/theta"] + names + ["/sigma"]
    est = Estimates(cmd="boxcox", depvar=vars_[0].name, names=enames,
                    b=np.concatenate([[th], fit["b"], [fit["sigma"]]]), V=fit["V"], N=n,
                    stat="z", rows=rows)
    est.eqnames = ["theta"] + ["Notrans"] * len(names) + ["sigma"]
    est.level = level
    est.title = "Box-Cox regression"
    scal = [("N", float(n)), ("df_m", float(k)), ("ll0", ll0), ("df_r", 3.0), ("chi2", lr), ("ll", ll)]
    for key in ("t1", "tm1", "t0"):
        llr, chi, pv = tests[key]
        scal += [(f"chi2_{key}", chi), (f"p_{key}", pv), (f"ll_{key}", llr)]
    scal += [("ic", float(fit["nit"])), ("rc", 0.0), ("rank", 1.0)]      # VERIFY: e(ic), e(rank)
    cl = f"boxcox {args.strip()}"
    post(s, est, scal,
         [("cmd", "boxcox"), ("marginsnotok", "_ALL"), ("cmdline", cl), ("ntrans", "yes"),
          ("model", "lhsonly"), ("predict", "boxcox_p"), ("chi2type", "LR"),
          ("depvar", vars_[0].name), ("properties", "b V")], sample)


def _boxcox_general(s: "Session", p, o: dict, mdl: str, args: str) -> None:
    """boxcox with model(rhsonly), model(lambda) or model(theta): the regressors not in
    notrans() are transformed by lambda. VERIFY: the layout, the tests of the theta model and
    the names in e() (compat 1913)."""
    from ..estimation.results import CoefRow, Estimates, coef_table, post
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    notrans = numeric_vars(s, str(o["notrans"]), allow_string=False) if o.get("notrans") not in (None, True) else []
    # every regressor of the varlist is transformed; one also in notrans() is dropped from there
    # as collinear (compat 1913)
    trans = vars_[1:]
    for v in [v for v in notrans if v in trans]:
        s.output.write(f"note: {v.name} dropped because of collinearity\n", "text")
        notrans.remove(v)
    if not trans:
        raise StataError(198, "no variables to transform")                         # VERIFY
    mask = touse(s, p)
    for v in vars_ + notrans:
        mask &= S.valid(v.data)
    y = vars_[0].data[mask]
    n = len(y)
    if mdl != "rhsonly" and (y <= 0).any():
        raise StataError(411, "dependent variable must be positive")              # VERIFY
    for v in trans:
        if (v.data[mask] <= 0).any():
            raise StataError(411, f"{v.name} contains observations that are not strictly positive")  # VERIFY: rc
    Xt = np.column_stack([v.data[mask] for v in trans])
    Xn = np.column_stack([v.data[mask] for v in notrans] + [np.ones(n)])
    fit = B.boxcox_general(y, Xt, Xn, mdl)
    log0 = []
    if mdl == "rhsonly":
        # no iterations for the comparison model (VERIFY)
        r0 = y - y.mean()
        ll0 = -n / 2 * (math.log(2 * math.pi) + 1) - n / 2 * math.log(float(r0 @ r0) / n)
    else:
        # the constant-only model with y transformed, by Newton-Raphson from 1 (compat 1913)
        nr0 = B.profile_newton(lambda t: B.boxcox_ll(float(t[0]), y, np.ones((n, 1))), [1.0])
        ll0, log0 = (nr0.ll, nr0.log) if nr0.converged else (B.boxcox_fit(y, np.ones((n, 1)))["ll"], [])
    ll = fit["ll"]
    lr = 2 * (ll - ll0)
    ncov = len(trans) + len(notrans)
    level = _level(s, o)
    out = s.output
    from ..estimation import ml as ML
    for head, log in (("Fitting comparison model", log0), ("Fitting full model", fit["log"])):
        out.write(head + "\n\n", "text")
        if log and not o.get("nolog"):
            ML.print_log(out, ML.MLResult(None, 0.0, None, None, 0, True, log))
            out.write("\n", "text")
    out.write(f"{'':>50}{'Number of obs':<16}= ", "text")
    out.write(f"{n:>10}\n", "result")
    out.write(f"{'':>50}{'LR chi2(' + str(ncov) + ')':<16}= ", "text")
    out.write(f"{lr:>10.2f}\n", "result")
    out.write(f"{'Log likelihood = ' + f'{ll:.8g}':<50}{'Prob > chi2':<16}= ", "text")
    out.write(f"{float(st.chi2.sf(lr, ncov)):>10.3f}\n", "result")
    out.write(" \n", "text")
    npar = 2 if mdl == "theta" else 1
    V = fit["V"]
    se = np.sqrt(np.maximum(np.diag(V), 0))
    rows = [CoefRow("coef", "/lambda", fit["lambda"], se[1 if mdl == "theta" else 0], name="lambda")]
    if mdl == "theta":
        rows.append(CoefRow("coef", "/theta", fit["theta"], se[0], name="theta"))
    coef_table(s, vars_[0].name, rows, stat="z", df=None, level=level)
    out.write(" \nEstimates of scale-variant parameters\n", "text")
    out.write("-" * 28 + "\n", "text")
    out.write(f"{'':>13}|{'Coef.':>11}\n", "text")
    sep = "-" * 13 + "+" + "-" * 14 + "\n"
    out.write(sep, "text")
    bt, bn = fit["b"][:len(trans)], fit["b"][len(trans):]
    nnames = [v.name for v in notrans] + ["_cons"]
    out.write(f"{'Notrans':<13}|\n", "text")
    for nm, b in zip(nnames, bn):
        out.write(f"{name12(nm):>12} |", "text")
        out.write(f"{g9(b):>11}\n", "result")
    out.write(sep, "text")
    out.write(f"{'Trans':<13}|\n", "text")
    for v, b in zip(trans, bt):
        out.write(f"{name12(v.name):>12} |", "text")
        out.write(f"{g9(b):>11}\n", "result")
    out.write(sep, "text")
    out.write(f"{'/sigma':>12} |", "text")
    out.write(f"{g9(fit['sigma']):>11}\n", "result")
    out.write("-" * 28 + "\n\n", "text")
    out.write("-" * 57 + "\n", "text")
    out.write("   Test         Restricted     LR statistic      P-value\n", "text")
    out.write("    H0:       log likelihood       chi2       Prob > chi2\n", "text")
    out.write("-" * 57 + "\n", "text")
    lab = "theta=lambda" if mdl == "theta" else "lambda"
    tests = {}
    for t0, key in ((-1, "tm1"), (0, "t0"), (1, "t1")):
        llr = fit["ll_at"](float(t0))
        chi = 2 * (ll - llr)
        pv = float(st.chi2.sf(chi, npar))
        # the restricted log likelihood ends in column 26 (compat 1913)
        tl = f"{lab} = {t0:2d}"
        out.write(tl, "text")
        out.write(f"{g9(llr, '%10.0g'):>{26 - len(tl)}}{chi:>13.2f}{pv:>16.3f}\n", "result")
        tests[key] = (llr, chi, pv)
    out.write("-" * 57 + "\n", "text")
    # e(b): the transformation parameters, the coefficients and sigma
    names = (["/theta"] if mdl == "theta" else []) + ["/lambda"] + nnames + [v.name for v in trans] + ["/sigma"]
    eqs = (["theta"] if mdl == "theta" else []) + ["lambda"] + ["Notrans"] * len(nnames) + ["Trans"] * len(trans) + ["sigma"]
    order = list(range(npar)) + [npar + len(trans) + j for j in range(len(nnames))] + \
        [npar + j for j in range(len(trans))] + [V.shape[0] - 1]
    bvec = np.concatenate([fit["params"], bn, bt, [fit["sigma"]]])
    est = Estimates(cmd="boxcox", depvar=vars_[0].name, names=names, b=bvec,
                    V=V[np.ix_(order, order)], N=n, stat="z", rows=rows)
    est.eqnames = eqs
    est.level = level
    est.title = "Box-Cox regression"
    # the order of model(lambda) in compat 1913: df_m is the number of parameters and df_r 3, as
    # with the default model (VERIFY: rhsonly and theta)
    scal = [("N", float(n)), ("df_m", float(len(bvec))), ("df_r", 3.0), ("chi2", lr), ("ll0", ll0), ("ll", ll)]
    for key in ("t1", "tm1", "t0"):
        llr, chi, pv = tests[key]
        scal += [(f"chi2_{key}", chi), (f"p_{key}", pv), (f"ll_{key}", llr)]
    scal += [("ic", float(fit["nit"])), ("rc", 0.0), ("rank", float(npar))]
    post(s, est, scal,
         [("cmd", "boxcox"), ("marginsnotok", "_ALL"), ("model", mdl), ("cmdline", f"boxcox {args.strip()}"),
          ("ntrans", "yes"), ("predict", "boxcox_p"),
          ("chi2type", "LR"), ("depvar", vars_[0].name), ("properties", "b V")], mask)
