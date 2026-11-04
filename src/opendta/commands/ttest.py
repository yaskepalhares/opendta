"""ci, ttest and prtest ([R] ci, ttest, prtest): confidence intervals for means, proportions and
variances, and the one-sample, paired and two-sample tests of means and proportions."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M, stats as S
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.functions import _abbrev
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .data import _by_header
from .registry import command
from .summarize import _level, by_groups, g9, name12, numeric_vars, weights

if TYPE_CHECKING:
    from ..session import Session


@command("ci", "ci", byable=True)
def cmd_ci(s: "Session", args: str) -> None:
    """ci [means|proportions|variances] varlist ([R] ci, syntax of Stata 14.1; ci varlist alone
    means ci means)."""
    from .immediate import _prop_ci, variance_row, variance_table_head
    t = args.strip()
    kind = "means"
    head = t.split(" ", 1)[0].split(",", 1)[0]
    for k in ("means", "proportions", "variances"):
        if len(head) >= 4 and k.startswith(head) or head == k[:-1]:
            kind = k
            t = t[len(head):].strip()
            break
    p = parse_standard(t)
    o = match_options(p.options, {"level": 1, "binomial": 1, "poisson": 1, "exact": 2,
                                  "wald": 2, "wilson": 2, "agresti": 2, "jeffreys": 2,
                                  "exposure": 3, "separator": 3, "total": 3, "sd": 2,
                                  "bonett": 3}) if p.options.strip() else {}
    if o.get("binomial"):
        kind = "proportions"
    level = _level(s, o)
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("aweight", "fweight") if kind == "means" else ("fweight",))
    out = s.output
    method = next((m for m in ("wald", "wilson", "agresti", "jeffreys") if o.get(m)), "exact")
    titles = {"exact": "-- Binomial Exact --", "wald": "-- Binomial Wald ---",
              "wilson": "------ Wilson ------", "agresti": "-- Agresti-Coull ---",
              "jeffreys": "----- Jeffreys -----"}
    poisson = kind == "means" and bool(o.get("poisson"))
    expo = s.data.get(s.expand_varlist(o["exposure"])[0]) if isinstance(o.get("exposure"), str) else None
    for group, first in by_groups(s):
        if first is not None:
            _by_header(s, first)
        mask = base & group
        if kind == "variances":
            variance_table_head(s, level, bool(o.get("sd")), bool(o.get("bonett")))
            for v in vars_:
                sel = mask & S.valid(v.data)
                xs = v.data[sel]
                m = S.moments(xs, w[sel] if w is not None else None, wtype)
                n = int(m.N)
                if n < 2:
                    raise StataError(2001, "insufficient observations")
                xb = None
                if o.get("bonett"):
                    xb = np.repeat(xs, w[sel].astype(int)) if w is not None else xs
                variance_row(s, v.name, n, m.var, level, bool(o.get("sd")), xb)
            continue
        lv = f"[{level:g}% Conf. Interval]"
        if kind == "proportions":
            out.write("\n" + " " * 48 + titles[method].rjust(30) + "\n", "text")
        elif poisson:
            out.write("\n" + " " * 48 + "-- Poisson  Exact --".rjust(30) + "\n", "text")   # VERIFY
        else:
            out.write("\n", "text")
        obs = "Exposure" if poisson else "Obs"
        out.write(f"    Variable |{obs:>11}        Mean    Std. Err.       {lv}\n", "text")
        out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
        for v in vars_:
            x = v.data
            sel = mask & S.valid(x)
            if expo is not None:
                sel &= S.valid(expo.data)
            xs = x[sel]
            ws = w[sel] if w is not None else None
            m = S.moments(xs, ws, wtype)
            n = m.N
            if kind == "proportions":
                if not np.all((xs == 0) | (xs == 1)):
                    raise StataError(2000, f"{v.name} is not a 0/1 variable")   # VERIFY
                k = float((xs * (ws if ws is not None else 1)).sum())
                mean = k / n if n else M.SYSMISS
                se = math.sqrt(mean * (1 - mean) / n) if n else M.SYSMISS
                lo, hi = _prop_ci(k, n, level, method) if n else (M.SYSMISS, M.SYSMISS)
            elif poisson:
                # exact Poisson interval for the count per unit of exposure ([R] ci)
                k = float((xs * (ws if ws is not None else 1)).sum())
                E = float((expo.data[sel] * (ws if ws is not None else 1)).sum()) if expo is not None else n
                n = E
                mean = k / E if E else M.SYSMISS
                se = math.sqrt(k) / E if E else M.SYSMISS
                a = (100 - level) / 200
                from scipy import stats as st
                lo = float(st.chi2.ppf(a, 2 * k)) / 2 / E if k > 0 else 0.0
                hi = float(st.chi2.ppf(1 - a, 2 * (k + 1))) / 2 / E
            else:
                mean = m.mean
                se = m.sd / math.sqrt(n) if n > 0 and m.sd < M.SYSMISS else M.SYSMISS
                lo, hi = S.t_ci(mean, se, n - 1, level) if n > 1 else (M.SYSMISS, M.SYSMISS)
            out.write(f"{name12(v.name):>12} |", "text")
            out.write(f" {n:>10,.0f}" + "   " + g9(mean) + "   " + g9(se) + "    " + g9(lo).rjust(12)
                      + "   " + g9(hi).rjust(9) + "\n", "result")
            s.r = {"N": float(n), "mean": mean, "se": se, "lb": lo, "ub": hi, "level": level}
            if kind == "proportions":
                s.r["proportion"] = s.r.pop("mean")


def _clopper(k: float, n: float, level: float) -> tuple[float, float]:
    from scipy import stats as st
    a = 1 - level / 100
    lo = 0.0 if k == 0 else float(st.beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(st.beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def _ttest_table(s: "Session", rows: list[tuple[str, float, float, float, float, float, float]],
                 level: float, *, title: str, foot_rows: list | None = None,
                 head: str = "Variable", close: bool = True) -> None:
    out = s.output
    out.write(f"\n{title}\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(f"{head:>8} |     Obs        Mean    Std. Err.   Std. Dev.   [{level:g}% Conf. Interval]\n",
              "text")
    out.write("-" * 9 + "+" + "-" * 68 + "\n", "text")

    def line(name, n, mean, se, sd, lo, hi):
        out.write(f"{_abbrev(name, 8):>8} |", "text")
        out.write(f"{n:>8,.0f}" + "   " + g9(mean) + "   " + g9(se) + "   " + g9(sd) + "   " + g9(lo)
                  + "   " + g9(hi) + "\n", "result")
    for r in rows:
        line(*r)
    if foot_rows:
        out.write("-" * 9 + "+" + "-" * 68 + "\n", "text")
        for r in foot_rows:
            line(*r)
    if close:
        out.write("-" * 78 + "\n", "text")


def _pvals(t: float, df: float | None) -> tuple[float, float, float]:
    from scipy import stats as st
    dist = st.t(df) if df is not None else st.norm()
    lower = float(dist.cdf(t))
    upper = float(dist.sf(t))
    two = float(2 * dist.sf(abs(t)))
    return lower, two, upper


def _ha_lines(s: "Session", diff_name: str, hyp: str, stat: str, t: float, df: float | None,
              two_op: str = ">") -> None:
    """The three alternative hypotheses: the probabilities start at columns
    1, 28 and 60; each "Ha:" is centered over its probability."""
    # VERIFY: positions checked against a single Stata 14 example
    out = s.output
    lo, two, up = _pvals(t, df)
    T = stat.upper()
    probs = [f"Pr({T} < {stat}) = {lo:.4f}", f"Pr(|{T}| {two_op} |{stat}|) = {two:.4f}",
             f"Pr({T} > {stat}) = {up:.4f}"]
    has = [f"Ha: {diff_name} < {hyp}", f"Ha: {diff_name} != {hyp}", f"Ha: {diff_name} > {hyp}"]
    starts = [1, 28, 60]
    line_p = [" "] * 80
    line_h = [" "] * 80
    for st_, pr, ha in zip(starts, probs, has):
        for k, ch in enumerate(pr):
            line_p[st_ + k] = ch
        hs = max(0, st_ + (len(pr) - len(ha) + 1) // 2)
        for k, ch in enumerate(ha):
            if hs + k < len(line_h):
                line_h[hs + k] = ch
    out.write("\n" + "".join(line_h).rstrip() + "\n", "text")
    out.write("".join(line_p).rstrip() + "\n", "result")


@command("ttest", byable=True)
def cmd_ttest(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "unequal": 3, "welch": 1, "level": 1,
                                  "unpaired": 4}) if p.options.strip() else {}
    level = _level(s, o)
    mask = touse(s, p)
    if p.exp is None and "==" in p.varlist:          # ttest x == 15 / ttest x == y
        left, _, right = p.varlist.partition("==")
        p.varlist, p.exp = left.strip(), right.strip()
    if p.exp is None:
        if not o.get("by"):
            raise StataError(198, "option by() required")   # VERIFY
        _ttest_by(s, p.varlist, str(o["by"]).strip(), mask, level, bool(o.get("unequal")),
                  bool(o.get("welch")))
        return
    v = numeric_vars(s, p.varlist, allow_string=False)
    if len(v) != 1:
        raise StataError(103, "too many variables specified")
    v = v[0]
    rhs = p.exp.strip()
    ds = s.data
    if ds.has(rhs) or (rhs and not _is_number(rhs) and _try_var(s, rhs)):
        other = ds.get(_try_var(s, rhs))
        if o.get("unpaired"):
            _ttest_two(s, v.name, other.name, v.data[mask & S.valid(v.data)],
                       other.data[mask & S.valid(other.data)], level, bool(o.get("unequal")),
                       bool(o.get("welch")), names=(v.name, other.name))
            return
        _ttest_paired(s, v, other, mask, level)
        return
    mu = float(s.eval(rhs))
    x = v.data[mask & S.valid(v.data)]
    m = S.moments(x)
    _ttest_one(s, v.name, len(x), m.mean, m.sd, mu, level)


def _ttest_one(s: "Session", name: str, n: int, mean: float, sd: float, mu: float, level: float) -> None:
    """One-sample t test from the statistics (ttest and ttesti)."""
    se = sd / math.sqrt(n) if n > 1 else M.SYSMISS
    lo, hi = S.t_ci(mean, se, n - 1, level)
    _ttest_table(s, [(name, n, mean, se, sd, lo, hi)], level, title="One-sample t test")
    t = (mean - mu) / se if se and se < M.SYSMISS else M.SYSMISS
    out = s.output
    out.write(f"    mean = mean({name})".ljust(66) + f"t = {t:>8.4f}\n", "text")
    out.write(f"Ho: mean = {_num(mu)}".ljust(49) + f"degrees of freedom = {n - 1:>8}\n", "text")
    _ha_lines(s, "mean", _num(mu), "t", t, n - 1)
    lo_p, two, up = _pvals(t, n - 1)
    s.r = {"N_1": float(n), "mu_1": mean, "sd_1": sd, "se": se, "t": t, "df_t": float(n - 1),
           "p_l": lo_p, "p": two, "p_u": up, "level": level}


def _num(x: float) -> str:
    return format_value(x, "%9.0g", pad=False).strip()


def _is_number(t: str) -> bool:
    try:
        float(t)
        return True
    except ValueError:
        return False


def _try_var(s: "Session", name: str) -> str | None:
    try:
        names = s.expand_varlist(name)
    except StataError:
        return None
    return names[0] if len(names) == 1 else None


def _ttest_paired(s, v, other, mask, level):
    sel = mask & S.valid(v.data) & S.valid(other.data)
    x, y = v.data[sel], other.data[sel]
    d = x - y
    n = len(d)
    rows = []
    for name, arr in ((v.name, x), (other.name, y)):
        m = S.moments(arr)
        se = m.sd / math.sqrt(n)
        rows.append((name, n, m.mean, se, m.sd, *S.t_ci(m.mean, se, n - 1, level)))
    md = S.moments(d)
    sed = md.sd / math.sqrt(n)
    foot = [("diff", n, md.mean, sed, md.sd, *S.t_ci(md.mean, sed, n - 1, level))]
    _ttest_table(s, rows, level, title="Paired t test", foot_rows=foot)
    t = md.mean / sed if sed else M.SYSMISS
    out = s.output
    out.write(f"     mean(diff) = mean({v.name} - {other.name})".ljust(66) + f"t = {t:>8.4f}\n", "text")
    out.write(" Ho: mean(diff) = 0".ljust(49) + f"degrees of freedom = {n - 1:>8}\n", "text")
    _ha_lines(s, "mean(diff)", "0", "t", t, n - 1)
    lo_p, two, up = _pvals(t, n - 1)
    s.r = {"N_1": float(n), "N_2": float(n), "mu_1": rows[0][2], "mu_2": rows[1][2], "t": t,
           "df_t": float(n - 1), "p": two, "p_l": lo_p, "p_u": up, "se": sed,
           "sd_1": rows[0][4], "sd_2": rows[1][4], "level": level}


def _ttest_by(s, varlist, byvar, mask, level, unequal, welch):
    v = numeric_vars(s, varlist, allow_string=False)[0]
    g = s.data.get(s.expand_varlist(byvar)[0])
    sel = mask & S.valid(v.data) & (S.valid(g.data) if not g.is_string else np.array([x != "" for x in g.raw]))
    groups = sorted(set(g.data[sel].tolist()))
    if len(groups) != 2:
        raise StataError(420, "more than 2 groups found, only 2 allowed")   # VERIFY
    from .inspect import cell_text
    names = []
    arrays = []
    for gv in groups:
        idx = np.flatnonzero(sel & (g.data == gv))
        names.append(cell_text(s.data, g, int(idx[0])))
        arrays.append(v.data[idx])
    _ttest_two(s, v.name, g.name, arrays[0], arrays[1], level, unequal, welch, names=names)


def _ttest_two(s, vname, gname, x1, x2, level, unequal, welch, *, names):
    n1, n2 = len(x1), len(x2)
    m1, m2 = S.moments(x1), S.moments(x2)
    se1, se2 = m1.sd / math.sqrt(n1), m2.sd / math.sqrt(n2)
    allx = np.concatenate([x1, x2])
    ma = S.moments(allx)
    sea = ma.sd / math.sqrt(len(allx))
    diff = m1.mean - m2.mean
    if unequal or welch:
        v1, v2 = m1.var / n1, m2.var / n2
        sed = math.sqrt(v1 + v2)
        if welch:
            df = -2 + (v1 + v2) ** 2 / (v1 ** 2 / (n1 + 1) + v2 ** 2 / (n2 + 1))
        else:
            df = (v1 + v2) ** 2 / (v1 ** 2 / (n1 - 1) + v2 ** 2 / (n2 - 1))
    else:
        sp2 = ((n1 - 1) * m1.var + (n2 - 1) * m2.var) / (n1 + n2 - 2)
        sed = math.sqrt(sp2 * (1 / n1 + 1 / n2))
        df = n1 + n2 - 2
    rows = [(str(names[0]), n1, m1.mean, se1, m1.sd, *S.t_ci(m1.mean, se1, n1 - 1, level)),
            (str(names[1]), n2, m2.mean, se2, m2.sd, *S.t_ci(m2.mean, se2, n2 - 1, level))]
    foot = [("combined", len(allx), ma.mean, sea, ma.sd, *S.t_ci(ma.mean, sea, len(allx) - 1, level))]
    title = "Two-sample t test with unequal variances" if (unequal or welch) else \
        "Two-sample t test with equal variances"
    # the "diff" line comes right after "combined", without a full line between them
    _ttest_table(s, rows, level, title=title, foot_rows=foot, head="Group", close=False)
    lo, hi = S.t_ci(diff, sed, df, level)
    out = s.output
    out.write("-" * 9 + "+" + "-" * 68 + "\n", "text")
    out.write(f"{'diff':>8} |" + " " * 8 + "   " + g9(diff) + "   " + g9(sed) + " " * 15 + g9(lo)
              + "   " + g9(hi) + "\n", "result")
    out.write("-" * 78 + "\n", "text")
    t = diff / sed
    out.write(f"    diff = mean({names[0]}) - mean({names[1]})".ljust(66) + f"t = {t:>8.4f}\n", "text")
    dflabel = "Satterthwaite's degrees of freedom" if unequal else (
        "Welch's degrees of freedom" if welch else "degrees of freedom")
    if unequal or welch:
        # non-integer degrees of freedom in %8.0g ("5.58439"; Stata 14)
        out.write("Ho: diff = 0".ljust(78 - len(dflabel) - 3 - 8) + f"{dflabel} = "
                  + format_value(df, "%8.0g") + "\n", "text")
    else:
        out.write("Ho: diff = 0".ljust(49) + f"degrees of freedom = {df:>8.0f}\n", "text")
    _ha_lines(s, "diff", "0", "t", t, df)
    lo_p, two, up = _pvals(t, df)
    s.r = {"N_1": float(n1), "N_2": float(n2), "mu_1": m1.mean, "mu_2": m2.mean, "sd_1": m1.sd,
           "sd_2": m2.sd, "se": sed, "t": t, "df_t": float(df), "p": two, "p_l": lo_p, "p_u": up,
           "level": level}


@command("prtest", byable=True)
def cmd_prtest(s: "Session", args: str) -> None:
    """prtest varname == #p, prtest varname1 == varname2 and prtest varname, by(groupvar)."""
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "level": 1}) if p.options.strip() else {}
    level = _level(s, o)
    mask = touse(s, p)
    if p.exp is None and "==" in p.varlist:          # prtest b == 0.5
        left, _, right = p.varlist.partition("==")
        p.varlist, p.exp = left.strip(), right.strip()
    v = numeric_vars(s, p.varlist, allow_string=False)[0]

    def zero_one(x, name):
        if not np.all((x == 0) | (x == 1)):
            raise StataError(2000, f"{name} is not a 0/1 variable")   # VERIFY

    if isinstance(o.get("by"), str):
        if p.exp is not None:
            raise StataError(198, "invalid syntax")
        g = s.data.get(s.expand_varlist(o["by"])[0])
        sel = mask & S.valid(v.data) & (S.valid(g.data) if not g.is_string
                                         else np.array([x != "" for x in g.raw]))
        groups = sorted(set(g.data[sel].tolist()) if not g.is_string else set(np.asarray(g.raw)[sel]))
        if len(groups) != 2:
            raise StataError(420, "more than 2 groups found, only 2 allowed")   # VERIFY
        from .inspect import cell_text
        col = g.data if not g.is_string else np.asarray(g.raw)
        xs, names = [], []
        for gv in groups:
            idx = np.flatnonzero(sel & (col == gv))
            names.append(cell_text(s.data, g, int(idx[0])))
            xs.append(v.data[idx])
            zero_one(v.data[idx], v.name)
        _prtest_two(s, names, xs[0], xs[1], level)
        return
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    other = _try_var(s, p.exp.strip())
    if other is not None:                       # prtest x == y: two variables
        w2 = s.data.get(other)
        x1 = v.data[mask & S.valid(v.data)]
        x2 = w2.data[mask & S.valid(w2.data)]
        zero_one(x1, v.name)
        zero_one(x2, w2.name)
        _prtest_two(s, [v.name, w2.name], x1, x2, level)
        return
    x = v.data[mask & S.valid(v.data)]
    zero_one(x, v.name)
    p0 = float(s.eval(p.exp))
    n = len(x)
    ph = float(x.mean()) if n else M.SYSMISS
    _prtest_one(s, v.name, n, ph, p0, level)


def _prtest_two(s: "Session", names, x1, x2, level: float, *, stats: tuple | None = None) -> None:
    """Two-sample test of proportions ([R] prtest; layout VERIFY, case compat 0805)."""
    if stats is None:
        n1, n2 = len(x1), len(x2)
        p1, p2 = float(np.mean(x1)), float(np.mean(x2))
    else:
        n1, p1, n2, p2 = stats
    if n1 == 0 or n2 == 0:
        raise StataError(2000, "no observations")
    zc = S.z_crit(level)
    se1, se2 = math.sqrt(p1 * (1 - p1) / n1), math.sqrt(p2 * (1 - p2) / n2)
    diff = p1 - p2
    se = math.sqrt(se1 ** 2 + se2 ** 2)
    pp = (p1 * n1 + p2 * n2) / (n1 + n2)
    se0 = math.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
    z = diff / se0 if se0 > 0 else M.SYSMISS
    out = s.output
    a, b = str(names[0]), str(names[1])
    title = "Two-sample test of proportions"
    out.write("\n" + title + f"{name12(a)}: Number of obs = {n1:>8,}".rjust(78 - len(title)) + "\n", "text")
    out.write(f"{name12(b)}: Number of obs = {n2:>8,}".rjust(78) + "\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(f"{'Variable':>12} |       Mean   Std. Err.      z    P>|z|     [{level:g}% Conf. Interval]\n",
              "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    for nm, ph, se_ in ((a, p1, se1), (b, p2, se2)):
        out.write(f"{name12(nm):>12} |", "text")
        out.write("  " + g9(ph) + "  " + g9(se_) + " " * 21 + g9(ph - zc * se_) + "   "
                  + g9(ph + zc * se_) + "\n", "result")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    out.write(f"{'diff':>12} |", "text")
    out.write("  " + g9(diff) + "  " + g9(se) + " " * 21 + g9(diff - zc * se) + "   "
              + g9(diff + zc * se) + "\n", "result")
    lo_p, two, up = _pvals(z, None)
    out.write(" " * 12 + " |  under Ho:  ", "text")
    out.write(g9(se0) + f"{z:>9.2f}" + f"{two:>8.3f}" + "\n", "result")
    out.write("-" * 78 + "\n", "text")
    lab = f"        diff = prop({a}) - prop({b})"
    out.write(lab.ljust(66) + f"z = {z:>8.4f}\n", "text")
    out.write("    Ho: diff = 0\n", "text")
    # Stata 14 labels the two-sided p of the two-sample test "Pr(|Z| < |z|)" (compat 0805)
    _ha_lines(s, "diff", "0", "z", z, None, two_op="<")
    s.r = {"N_1": float(n1), "N_2": float(n2), "P_1": p1, "P_2": p2, "z": z, "p": two,
           "p_l": lo_p, "p_u": up}


def _prtest_one(s: "Session", name: str, n: int, ph: float, p0: float, level: float) -> None:
    """One-proportion test from the statistics (prtest and prtesti)."""
    se = math.sqrt(ph * (1 - ph) / n)
    z = (ph - p0) / math.sqrt(p0 * (1 - p0) / n)
    zc = S.z_crit(level)
    out = s.output
    title = "One-sample test of proportion"
    out.write("\n" + title + f"{name}: Number of obs = {n:>8,}".rjust(78 - len(title)) + "\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(f"    Variable |       Mean   Std. Err.                     [{level:g}% Conf. Interval]\n",
              "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    out.write(f"{name12(name):>12} |", "text")
    out.write("  " + g9(ph) + "  " + g9(se) + " " * 21 + g9(ph - zc * se) + "   " + g9(ph + zc * se) + "\n",
              "result")
    out.write("-" * 78 + "\n", "text")
    out.write(f"    p = proportion({name})".ljust(66) + f"z = {z:>8.4f}\n", "text")
    # the hypothesis value shows a leading zero ("Ho: p = 0.5"; Stata 14)
    out.write(f"Ho: p = {p0:g}\n", "text")
    _ha_lines(s, "p", f"{p0:g}", "z", z, None)
    lo_p, two, up = _pvals(z, None)
    s.r = {"N_1": float(n), "P_1": ph, "z": z, "p": two, "p_l": lo_p, "p_u": up}
