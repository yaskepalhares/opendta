"""summarize, tabstat, ci, ttest, prtest, correlate/pwcorr, centile, pctile,
xtile e _pctile ([R] summarize e afins).

Layouts reproduzidos a partir do Stata 14; detalhes não documentados têm
a marca VERIFICAR e casos em compat/do/0301_summarize.do.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import stats as S
from ..core.dataset import Variable
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.functions import _abbrev
from ..lang.syntax import Parsed, match_options, parse_standard
from ._util import eval_vector, touse
from .data import _by_header
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def g9(x: float, fmt: str = "%9.0g") -> str:
    if x >= M.SYSMISS:
        return M.missing_name(x)
    return format_value(float(x), fmt, pad=True)


def name12(name: str) -> str:
    return _abbrev(name, 12)


def weights(s: "Session", p: Parsed, mask: np.ndarray, allowed: tuple[str, ...]):
    """(vetor de pesos ou None, tipo) e máscara sem pesos missing/negativos."""
    if p.weight is None:
        return None, "", mask
    wtype, wexp = p.weight
    if wtype == "weight":
        wtype = allowed[0]
    if wtype not in allowed:
        raise StataError(101, f"{wtype}s not allowed")
    w = np.broadcast_to(np.asarray(eval_vector(s, wexp), dtype=np.float64), (s.data.nobs,)).copy()
    ok = (w < M.SYSMISS)
    if wtype == "fweight":
        if np.any((w[mask & ok] != np.trunc(w[mask & ok]))):
            raise StataError(401, "may not use noninteger frequency weights")
        if np.any(w[mask & ok] < 0):
            raise StataError(402, "negative weights encountered")
    mask = mask & ok & (w > 0)
    return w, wtype, mask


def by_groups(s: "Session"):
    """Pares (máscara do grupo, índice da 1ª obs.) para o prefixo by; um só
    par com tudo verdadeiro quando não há by."""
    n = s.data.nobs
    if s.by_groups is None or n == 0:
        yield np.ones(n, dtype=bool), None
        return
    ids = s.by_groups.ids
    for gid in range(int(ids.max()) + 1):
        sel = ids == gid
        yield sel, int(np.flatnonzero(sel)[0])


def numeric_vars(s: "Session", varlist: str, *, allow_string: bool = True) -> list[Variable]:
    ds = s.data
    names = s.expand_varlist(varlist) if varlist.strip() else list(ds.names)
    out = [ds.get(n) for n in names]
    if not allow_string:
        for v in out:
            if v.is_string:
                raise StataError(109, "type mismatch")
    return out


# ---------------------------------------------------------------------------
# summarize
# ---------------------------------------------------------------------------

_SUM_HEADER = "    Variable |        Obs        Mean    Std. Dev.       Min        Max"
_SUM_LINE = "-" * 13 + "+" + "-" * 57


@command("summarize", "su", byable=True)
def cmd_summarize(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"detail": 1, "meanonly": 4, "format": 3, "separator": 3,
                                  "noformat": 5}) if p.options.strip() else {}
    vars_ = numeric_vars(s, p.varlist)
    sep = int(str(o.get("separator", 5))) if o.get("separator") is not None else 5
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("aweight", "fweight", "iweight"))
    if o.get("detail") and wtype == "iweight":
        raise StataError(101, "iweights not allowed with detail")   # VERIFICAR
    out = s.output
    for group, first in by_groups(s):
        if first is not None:
            _by_header(s, first)
        mask = base & group
        if o.get("meanonly"):
            for v in vars_:
                _summ_results(s, v, mask, w, wtype, detail=False)
            continue
        if o.get("detail"):
            for v in vars_:
                _detail(s, v, mask, w, wtype)
            continue
        out.write("\n" + _SUM_HEADER + "\n" + _SUM_LINE + "\n", "text")
        for k, v in enumerate(vars_):
            if sep and k and k % sep == 0:
                out.write(_SUM_LINE + "\n", "text")
            m = _summ_results(s, v, mask, w, wtype, detail=False)
            out.write(f"{name12(v.name):>12} |", "text")
            fmt = v.fmt if o.get("format") and not v.is_string else "%9.0g"
            if m is None or m.N == 0:
                out.write(f" {0:>10,}\n", "result")
                continue
            line = f" {m.N:>10,.0f}" + "   " + g9(m.mean, fmt).rjust(9) + "   " + g9(m.sd, "%9.0g").rjust(9) \
                + "  " + g9(m.min, fmt).rjust(9) + "  " + g9(m.max, fmt).rjust(9)
            out.write(line + "\n", "result")


def _summ_results(s: "Session", v: Variable, mask, w, wtype, *, detail: bool):
    """Calcula e publica r() para a variável; devolve os momentos."""
    if v.is_string:
        s.r = {"N": 0.0, "sum_w": 0.0, "sum": 0.0}
        return None
    x = v.data
    sel = mask & S.valid(x)
    xs = x[sel]
    ws = w[sel] if w is not None else None
    m = S.moments(xs, ws, wtype)
    r = {"N": float(m.N), "sum_w": float(m.sum_w), "mean": m.mean, "Var": m.var, "sd": m.sd,
         "min": m.min, "max": m.max, "sum": float(m.sum)}
    if m.N == 0:
        r = {"N": 0.0, "sum_w": 0.0, "sum": 0.0}
    s.r = r
    return m


def _detail(s: "Session", v: Variable, mask, w, wtype) -> None:
    out = s.output
    title = v.label or v.name
    out.write("\n" + title.center(61).rstrip() + "\n" + "-" * 61 + "\n", "text")
    m = _summ_results(s, v, mask, w, wtype, detail=True)
    if m is None or m.N == 0:
        out.write("no observations\n", "text")   # VERIFICAR
        return
    x = v.data
    sel = mask & S.valid(x)
    order = np.argsort(x[sel], kind="stable")
    xs = x[sel][order]
    ws = (w[sel][order] if w is not None else None)
    pw = ws if wtype in ("fweight", "aweight", "pweight") else None
    pct = {q: S.percentile(xs, q, pw) for q in (1, 5, 10, 25, 50, 75, 90, 95, 99)}
    smallest = list(xs[:4])
    largest = list(xs[-4:]) if len(xs) >= 4 else list(xs)
    largest = [None] * (4 - len(largest)) + largest

    def val(x) -> str:
        return "" if x is None else g9(float(x)).strip()

    stats_right = {10: ("Obs", f"{m.N:,.0f}"), 25: ("Sum of Wgt.", f"{m.sum_w:,.0f}" if m.sum_w == int(m.sum_w) else g9(m.sum_w).strip()),
                   50: ("Mean", val(m.mean)), "L": ("Std. Dev.", val(m.sd)), 90: ("Variance", val(m.var)),
                   95: ("Skewness", val(m.skewness)), 99: ("Kurtosis", val(m.kurtosis))}

    def row(label: str, q, third: str) -> str:
        line = f"{label:>3}" + f"{val(pct[q]) if q is not None else '':>13}" + f"{third:>15}"
        key = q if q is not None else "L"
        if key in stats_right:
            name, value = stats_right[key]
            line = line.ljust(31) + "       " + f"{name:<12}" + f"{value:>11}"
        return line.rstrip()

    out.write("      Percentiles      Smallest\n", "text")
    for k, q in enumerate((1, 5, 10, 25)):
        third = val(smallest[k]) if k < len(smallest) else ""
        out.write(row(f"{q}%", q, third) + "\n", "result")
    out.write("\n" + row("50%", 50, "") + "\n", "result")
    out.write(row("", None, "Largest") + "\n", "result")
    for k, q in enumerate((75, 90, 95, 99)):
        out.write(row(f"{q}%", q, val(largest[k])) + "\n", "result")
    r = s.r
    r.update({f"p{q}": pct[q] for q in pct})
    r["skewness"], r["kurtosis"] = m.skewness, m.kurtosis


# ---------------------------------------------------------------------------
# ci / ttest / prtest
# ---------------------------------------------------------------------------

def _level(s: "Session", o: dict) -> float:
    return float(o.get("level", s.settings.get("level", "95")))


@command("ci", "ci", byable=True)
def cmd_ci(s: "Session", args: str) -> None:
    t = args.strip()
    kind = "means"
    for k in ("means", "proportions", "variances"):
        if t.split(" ", 1)[0] in (k, k[:4], k[:-1]):
            kind = k
            t = t.split(" ", 1)[1] if " " in t else ""
            break
    p = parse_standard(t)
    o = match_options(p.options, {"level": 1, "binomial": 1, "poisson": 1, "exact": 2,
                                  "wald": 2, "wilson": 2, "agresti": 2, "jeffreys": 2,
                                  "exposure": 3, "separator": 3}) if p.options.strip() else {}
    if o.get("binomial"):
        kind = "proportions"
    level = _level(s, o)
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("aweight", "fweight"))
    out = s.output
    head_se = "Std. Err." if kind == "means" else "Std. Err."
    for group, first in by_groups(s):
        if first is not None:
            _by_header(s, first)
        mask = base & group
        lv = f"[{level:g}% Conf. Interval]"
        extra = "" if kind == "means" else ("          -- Binomial Exact --" if kind == "proportions" else "")
        if extra:
            out.write("\n" + " " * 48 + extra.strip().rjust(30) + "\n", "text")
        else:
            out.write("\n", "text")
        out.write(f"    Variable |        Obs        Mean    {head_se}       {lv}\n", "text")
        out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
        for v in vars_:
            x = v.data
            sel = mask & S.valid(x)
            xs = x[sel]
            ws = w[sel] if w is not None else None
            m = S.moments(xs, ws, wtype)
            n = m.N
            if kind == "proportions":
                if not np.all((xs == 0) | (xs == 1)):
                    raise StataError(2000, f"{v.name} is not a 0/1 variable")   # VERIFICAR
                k = float((xs * (ws if ws is not None else 1)).sum())
                phat = k / n if n else M.SYSMISS
                se = math.sqrt(phat * (1 - phat) / n) if n else M.SYSMISS
                lo, hi = _clopper(k, n, level)
                mean = phat
            else:
                mean = m.mean
                se = m.sd / math.sqrt(n) if n > 0 and m.sd < M.SYSMISS else M.SYSMISS
                lo, hi = S.t_ci(mean, se, n - 1, level) if n > 1 else (M.SYSMISS, M.SYSMISS)
            out.write(f"{name12(v.name):>12} |", "text")
            out.write(f" {n:>10,.0f}" + "   " + g9(mean) + "   " + g9(se) + "    " + g9(lo).rjust(11)
                      + "   " + g9(hi).rjust(9) + "\n", "result")
            s.r = {"N": float(n), "mean": mean, "se": se, "lb": lo, "ub": hi, "level": level}


def _clopper(k: float, n: float, level: float) -> tuple[float, float]:
    from scipy import stats as st
    a = 1 - level / 100
    lo = 0.0 if k == 0 else float(st.beta.ppf(a / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(st.beta.ppf(1 - a / 2, k + 1, n - k))
    return lo, hi


def _ttest_table(s: "Session", rows: list[tuple[str, float, float, float, float, float, float]],
                 level: float, *, title: str, foot_rows: list | None = None) -> None:
    out = s.output
    out.write(f"\n{title}\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(f"{'Variable':>8} |     Obs        Mean    Std. Err.   Std. Dev.   [{level:g}% Conf. Interval]\n",
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
    out.write("-" * 78 + "\n", "text")


def _pvals(t: float, df: float | None) -> tuple[float, float, float]:
    from scipy import stats as st
    dist = st.t(df) if df is not None else st.norm()
    lower = float(dist.cdf(t))
    upper = float(dist.sf(t))
    two = float(2 * dist.sf(abs(t)))
    return lower, two, upper


def _ha_lines(s: "Session", diff_name: str, hyp: str, stat: str, t: float, df: float | None) -> None:
    """As três hipóteses alternativas: as probabilidades começam nas colunas
    1, 28 e 60; cada "Ha:" fica centrado sobre a sua probabilidade."""
    # VERIFICAR: posições conferidas com um único exemplo do Stata 14
    out = s.output
    lo, two, up = _pvals(t, df)
    T = stat.upper()
    probs = [f"Pr({T} < {stat}) = {lo:.4f}", f"Pr(|{T}| > |{stat}|) = {two:.4f}",
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
            raise StataError(198, "option by() required")   # VERIFICAR
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
    n = len(x)
    m = S.moments(x)
    se = m.sd / math.sqrt(n) if n > 1 else M.SYSMISS
    lo, hi = S.t_ci(m.mean, se, n - 1, level)
    _ttest_table(s, [(v.name, n, m.mean, se, m.sd, lo, hi)], level, title="One-sample t test")
    t = (m.mean - mu) / se if se and se < M.SYSMISS else M.SYSMISS
    out = s.output
    out.write(f"    mean = mean({v.name})".ljust(66) + f"t = {t:>8.4f}\n", "text")
    out.write(f"Ho: mean = {_num(mu)}".ljust(49) + f"degrees of freedom = {n - 1:>8}\n", "text")
    _ha_lines(s, "mean", _num(mu), "t", t, n - 1)
    lo_p, two, up = _pvals(t, n - 1)
    s.r = {"N_1": float(n), "mu_1": m.mean, "sd_1": m.sd, "se": se, "t": t, "df_t": float(n - 1),
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
    out.write(f" mean(diff) = mean({v.name} - {other.name})".ljust(66) + f"t = {t:>8.4f}\n", "text")
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
        raise StataError(420, "more than 2 groups found, only 2 allowed")   # VERIFICAR
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
    _ttest_table(s, rows, level, title=title, foot_rows=foot)
    lo, hi = S.t_ci(diff, sed, df, level)
    out = s.output
    out.write(f"{'diff':>8} |" + " " * 8 + "   " + g9(diff) + "   " + g9(sed) + " " * 12 + "    " + g9(lo)
              + "   " + g9(hi) + "\n", "result")
    out.write("-" * 78 + "\n", "text")
    t = diff / sed
    out.write(f"    diff = mean({names[0]}) - mean({names[1]})".ljust(66) + f"t = {t:>8.4f}\n", "text")
    dflabel = "Satterthwaite's degrees of freedom" if unequal else (
        "Welch's degrees of freedom" if welch else "degrees of freedom")
    if unequal or welch:
        out.write("Ho: diff = 0".ljust(78 - len(dflabel) - 3 - 8) + f"{dflabel} = {df:>8.4f}\n", "text")
    else:
        out.write("Ho: diff = 0".ljust(49) + f"degrees of freedom = {df:>8.0f}\n", "text")
    _ha_lines(s, "diff", "0", "t", t, df)
    lo_p, two, up = _pvals(t, df)
    s.r = {"N_1": float(n1), "N_2": float(n2), "mu_1": m1.mean, "mu_2": m2.mean, "sd_1": m1.sd,
           "sd_2": m2.sd, "se": sed, "t": t, "df_t": float(df), "p": two, "p_l": lo_p, "p_u": up,
           "level": level}


@command("prtest", byable=True)
def cmd_prtest(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "level": 1}) if p.options.strip() else {}
    level = _level(s, o)
    mask = touse(s, p)
    if p.exp is None and "==" in p.varlist:          # prtest b == 0.5
        left, _, right = p.varlist.partition("==")
        p.varlist, p.exp = left.strip(), right.strip()
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    x = v.data[mask & S.valid(v.data)]
    if not np.all((x == 0) | (x == 1)):
        raise StataError(2000, f"{v.name} is not a 0/1 variable")   # VERIFICAR
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    p0 = float(s.eval(p.exp))
    n = len(x)
    ph = float(x.mean()) if n else M.SYSMISS
    se = math.sqrt(ph * (1 - ph) / n)
    z = (ph - p0) / math.sqrt(p0 * (1 - p0) / n)
    zc = S.z_crit(level)
    out = s.output
    out.write("\nOne-sample test of proportion" + f"{v.name}: Number of obs = {n:>10,}".rjust(49) + "\n",
              "text")
    out.write("-" * 78 + "\n", "text")
    out.write(f"{'Variable':>8} |       Mean   Std. Err.                     [{level:g}% Conf. Interval]\n",
              "text")
    out.write("-" * 9 + "+" + "-" * 68 + "\n", "text")
    out.write(f"{_abbrev(v.name, 8):>8} |", "text")
    out.write("   " + g9(ph) + "   " + g9(se) + " " * 25 + g9(ph - zc * se) + "   " + g9(ph + zc * se) + "\n",
              "result")
    out.write("-" * 78 + "\n", "text")
    out.write(f"    p = proportion({v.name})".ljust(66) + f"z = {z:>8.4f}\n", "text")
    out.write(f"Ho: p = {_num(p0)}\n", "text")
    _ha_lines(s, "p", _num(p0), "z", z, None)
    lo_p, two, up = _pvals(z, None)
    s.r = {"N_1": float(n), "P_1": ph, "z": z, "p": two, "p_l": lo_p, "p_u": up}


# ---------------------------------------------------------------------------
# correlate / pwcorr
# ---------------------------------------------------------------------------

def _corr_table(s: "Session", names: list[str], cells, *, extra=None, star_p=None) -> None:
    """cells[i][j] (j <= i): texto da correlação; extra: linhas por célula (sig., obs)."""
    out = s.output
    k = len(names)
    for start in range(0, k, 7):          # VERIFICAR: até 7 colunas por bloco
        cols = list(range(start, min(k, start + 7)))
        out.write("\n" + " " * 13 + "|" + "".join(f"{_abbrev(names[j], 8):>9}" for j in cols) + "\n",
                  "text")
        out.write("-" * 13 + "+" + "-" * (9 * len(cols)) + "\n", "text")
        for i in range(start, k):
            out.write(f"{name12(names[i]):>12} |", "text")
            out.write("".join(f"{cells[i][j]:>9}" for j in cols if j <= i) + "\n", "result")
            if extra:
                for ex in extra:
                    out.write(" " * 12 + " |", "text")
                    out.write("".join(f"{ex[i][j]:>9}" for j in cols if j <= i) + "\n", "result")
                out.write(" " * 12 + " |\n", "text")


@command("correlate", "cor", byable=True)
def cmd_correlate(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"means": 1, "covariance": 3, "wrap": 1,
                                  "noformat": 5}) if p.options.strip() else {}
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    mask = touse(s, p)
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight"))
    for v in vars_:
        mask &= S.valid(v.data)
    X = np.column_stack([v.data[mask] for v in vars_]) if vars_ else np.zeros((0, 0))
    n = int(mask.sum())
    if n == 0:
        raise StataError(2000, "no observations")
    ww = w[mask] if w is not None else np.ones(n)
    if wtype == "aweight":
        ww = ww * n / ww.sum()
    N = float(ww.sum())
    mean = (ww[:, None] * X).sum(0) / N
    D = X - mean
    cov = (ww[:, None] * D).T @ D / (N - 1)
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)
    s.output.write(f"(obs={n:,})\n", "text")
    names = [v.name for v in vars_]
    if o.get("means"):
        out = s.output
        out.write("\n    Variable |         Mean    Std. Dev.          Min          Max\n", "text")
        out.write("-" * 13 + "+" + "-" * 52 + "\n", "text")
        for j, v in enumerate(vars_):
            out.write(f"{name12(v.name):>12} |", "text")
            out.write(f"   {g9(mean[j])}    {g9(sd[j])}    {g9(X[:, j].min())}    {g9(X[:, j].max())}\n",
                      "result")
    mat = cov if o.get("covariance") else corr
    cells = [[(format_value(float(mat[i, j]), "%9.4f" if not o.get("covariance") else "%9.0g",
                            pad=False).strip()) for j in range(len(names))] for i in range(len(names))]
    _corr_table(s, names, cells)
    s.r = {"N": float(n)}
    if len(names) == 2:
        s.r["rho"] = float(corr[0, 1])
        if o.get("covariance"):
            s.r["cov_12"] = float(cov[0, 1])
        s.r["Var_1"], s.r["Var_2"] = float(cov[0, 0]), float(cov[1, 1])
    from .matrix import Matrix, _store
    _store(s)  # garante s.matrices
    s.r["C"] = Matrix(mat.copy(), names, names)


@command("pwcorr", "pwcorr", byable=True)
def cmd_pwcorr(s: "Session", args: str) -> None:
    from scipy import stats as st
    p = parse_standard(args)
    o = match_options(p.options, {"obs": 3, "sig": 3, "star": 2, "bonferroni": 1,
                                  "sidak": 2, "print": 2, "listwise": 4,
                                  "casewise": 4}) if p.options.strip() else {}
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    mask = touse(s, p)
    if o.get("listwise") or o.get("casewise"):
        for v in vars_:
            mask &= S.valid(v.data)
    k = len(vars_)
    R = np.full((k, k), np.nan)
    P = np.full((k, k), np.nan)
    Nn = np.zeros((k, k), dtype=int)
    for i in range(k):
        for j in range(i + 1):
            sel = mask & S.valid(vars_[i].data) & S.valid(vars_[j].data)
            x, y = vars_[i].data[sel], vars_[j].data[sel]
            n = len(x)
            Nn[i, j] = n
            if n > 2 and x.std() > 0 and y.std() > 0:
                r = float(np.corrcoef(x, y)[0, 1])
                R[i, j] = r
                t = r * math.sqrt((n - 2) / max(1e-300, 1 - r * r))
                P[i, j] = float(2 * st.t.sf(abs(t), n - 2)) if i != j else M.SYSMISS
            elif n > 0 and i == j:
                R[i, j] = 1.0
    m = k * (k - 1) / 2
    if o.get("bonferroni") and m:
        P = np.minimum(1.0, P * m)
    elif o.get("sidak") and m:
        P = np.minimum(1.0, 1 - (1 - P) ** m)
    star = float(str(o["star"])) if o.get("star") not in (None, True) else None
    pr = float(str(o["print"])) if o.get("print") not in (None, True) else None

    def cell(i, j):
        if np.isnan(R[i, j]):
            return "."
        if pr is not None and i != j and P[i, j] > pr:
            return ""
        txt = format_value(float(R[i, j]), "%9.4f", pad=False).strip()
        if star is not None and i != j and P[i, j] <= star:
            txt += "*"
        return txt
    names = [v.name for v in vars_]
    cells = [[cell(i, j) for j in range(k)] for i in range(k)]
    extra = []
    if o.get("sig"):
        extra.append([["" if (i == j or np.isnan(P[i, j])) else f"{P[i, j]:.4f}" for j in range(k)]
                      for i in range(k)])
    if o.get("obs"):
        extra.append([[f"{Nn[i, j]:,}" for j in range(k)] for i in range(k)])
    _corr_table(s, names, cells, extra=extra or None)
    if k == 2:
        s.r = {"N": float(Nn[1, 0]), "rho": float(R[1, 0])}
    from .matrix import Matrix
    C = np.where(np.isnan(R), M.SYSMISS, R)
    C = np.tril(C) + np.tril(C, -1).T
    s.r = {**(s.r if k == 2 else {}), "C": Matrix(C, names, names)}


# ---------------------------------------------------------------------------
# centile / pctile / xtile / _pctile
# ---------------------------------------------------------------------------

@command("centile", byable=True)
def cmd_centile(s: "Session", args: str) -> None:
    from scipy import stats as st
    p = parse_standard(args)
    o = match_options(p.options, {"centile": 1, "cci": 3, "normal": 1, "meansd": 1,
                                  "level": 1}) if p.options.strip() else {}
    level = _level(s, o)
    from ..lang.words import parse_numlist
    cents = parse_numlist(str(o["centile"])) if o.get("centile") else [50.0]
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    mask = touse(s, p)
    out = s.output
    head = "-- Binom. Interp. --" if not (o.get("normal") or o.get("meansd")) else "-- Normal, based on --"
    out.write("\n" + " " * 54 + head + "\n", "text")
    out.write(f"    Variable |       Obs  Percentile    Centile        [{level:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 61 + "\n", "text")
    alpha = 1 - level / 100
    for v in vars_:
        xs = np.sort(v.data[mask & S.valid(v.data)])
        n = len(xs)
        for k, c in enumerate(cents):
            q = c / 100
            val = S.centile_value(xs, c)
            # VERIFICAR: limites binomiais interpolados
            j = int(st.binom.ppf(alpha / 2, n, q))
            kk = int(st.binom.ppf(1 - alpha / 2, n, q)) + 1
            lo = float(xs[max(0, min(n - 1, j - 1))]) if n else M.SYSMISS
            hi = float(xs[max(0, min(n - 1, kk - 1))]) if n else M.SYSMISS
            name = f"{name12(v.name):>12}" if k == 0 else " " * 12
            out.write(name + " |", "text")
            obs = f"{n:>10,}" if k == 0 else " " * 10
            out.write(f"{obs}{c:>12g}    {g9(val)}       {g9(lo)}   {g9(hi)}\n", "result")
            s.r = {"N": float(n), f"c_{k + 1}": val, f"lb_{k + 1}": lo, f"ub_{k + 1}": hi}


def _pct_values(xs: np.ndarray, w, qs: list[float]) -> list[float]:
    return [S.percentile(xs, q, w) for q in qs]


@command("_pctile")
def cmd__pctile(s: "Session", args: str) -> None:
    from ..lang.words import parse_numlist
    p = parse_standard(args)
    o = match_options(p.options, {"nquantiles": 2, "percentiles": 1, "altdef": 3}) if p.options.strip() else {}
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p)
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight", "pweight"))
    sel = mask & S.valid(v.data)
    order = np.argsort(v.data[sel], kind="stable")
    xs = v.data[sel][order]
    ws = w[sel][order] if w is not None else None
    if o.get("percentiles"):
        qs = parse_numlist(str(o["percentiles"]))
    else:
        nq = int(str(o.get("nquantiles", 2)))
        qs = [100 * k / nq for k in range(1, nq)]
    s.r = {f"r{k}": val for k, val in enumerate(_pct_values(xs, ws, qs), start=1)}


@command("pctile")
def cmd_pctile(s: "Session", args: str) -> None:
    from ..lang.words import parse_numlist
    p = parse_standard(args)
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"nquantiles": 2, "percentiles": 1, "genp": 4,
                                  "altdef": 3}) if p.options.strip() else {}
    newname = p.varlist.split()[-1]
    vtype = p.varlist.split()[0] if len(p.varlist.split()) == 2 else "float"
    x = np.broadcast_to(np.asarray(eval_vector(s, p.exp), dtype=np.float64), (s.data.nobs,))
    mask = touse(s, p) & S.valid(x)
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight", "pweight"))
    order = np.argsort(x[mask], kind="stable")
    xs = x[mask][order]
    ws = w[mask][order] if w is not None else None
    if o.get("percentiles"):
        qs = parse_numlist(str(o["percentiles"]))
    else:
        nq = int(str(o.get("nquantiles", 2)))
        qs = [100 * k / nq for k in range(1, nq)]
    vals = _pct_values(xs, ws, qs)
    n = s.data.nobs
    if len(vals) > n:
        raise StataError(198, "nquantiles() too large for the number of observations")   # VERIFICAR
    col = np.full(n, M.SYSMISS)
    col[: len(vals)] = vals
    s.data.add(Variable(newname, vtype, col))
    if o.get("genp"):
        pcol = np.full(n, M.SYSMISS)
        pcol[: len(qs)] = qs
        s.data.add(Variable(str(o["genp"]).strip(), "float", pcol))
    s.notify_state()


@command("xtile")
def cmd_xtile(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"nquantiles": 2, "cutpoints": 2, "altdef": 3}) if p.options.strip() else {}
    newname = p.varlist.split()[-1]
    x = np.broadcast_to(np.asarray(eval_vector(s, p.exp), dtype=np.float64), (s.data.nobs,))
    mask = touse(s, p) & S.valid(x)
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight", "pweight"))
    if o.get("cutpoints"):
        cv = s.data.get(s.expand_varlist(str(o["cutpoints"]))[0]).data
        cuts = np.sort(cv[S.valid(cv)])
    else:
        nq = int(str(o.get("nquantiles", 2)))
        order = np.argsort(x[mask], kind="stable")
        xs = x[mask][order]
        ws = w[mask][order] if w is not None else None
        cuts = np.array(_pct_values(xs, ws, [100 * k / nq for k in range(1, nq)]))
    col = np.full(s.data.nobs, M.SYSMISS)
    # categoria = 1 + nº de pontos de corte estritamente menores que x
    col[mask] = 1 + np.searchsorted(cuts, x[mask], side="left")
    s.data.add(Variable(newname, "byte" if len(cuts) < 100 else "int", col))
    s.notify_state()


# ---------------------------------------------------------------------------
# tabstat
# ---------------------------------------------------------------------------

_TS_LABEL = {"mean": "mean", "count": "N", "n": "N", "sum": "sum", "max": "max", "min": "min",
             "range": "range", "sd": "sd", "variance": "variance", "var": "variance", "cv": "cv",
             "semean": "se(mean)", "skewness": "skewness", "kurtosis": "kurtosis", "median": "p50",
             "iqr": "iqr"}


def _ts_stats(spec: str) -> list[str]:
    out = []
    for w in spec.split():
        w = w.lower()
        if w == "q":
            out += ["p25", "p50", "p75"]
        elif w in _TS_LABEL or (w.startswith("p") and w[1:].isdigit()):
            out.append(w)
        else:
            raise StataError(198, f"{w} invalid statistic")   # VERIFICAR
    return out


def _ts_value(stat: str, x: np.ndarray, w, wtype) -> float:
    m = S.moments(x, w, wtype)
    if stat in ("count", "n"):
        return float(m.N)
    if m.N == 0:
        return M.SYSMISS
    simple = {"mean": m.mean, "sum": float(m.sum), "max": m.max, "min": m.min,
              "range": m.max - m.min, "sd": m.sd, "variance": m.var, "var": m.var,
              "skewness": m.skewness, "kurtosis": m.kurtosis}
    if stat in simple:
        return simple[stat]
    if stat == "cv":
        return m.sd / m.mean if m.mean else M.SYSMISS
    if stat == "semean":
        return m.sd / math.sqrt(m.N) if m.N > 0 and m.sd < M.SYSMISS else M.SYSMISS
    order = np.argsort(x, kind="stable")
    xs = x[order]
    ws = w[order] if w is not None and wtype in ("fweight", "aweight") else None
    if stat == "median":
        return S.percentile(xs, 50, ws)
    if stat == "iqr":
        return S.percentile(xs, 75, ws) - S.percentile(xs, 25, ws)
    return S.percentile(xs, float(stat[1:]), ws)


@command("tabstat")
def cmd_tabstat(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "statistics": 1, "stats": 5, "columns": 1, "format": 1,
                                  "nototal": 3, "missing": 1, "longstub": 1, "save": 2,
                                  "casewise": 1, "labelwidth": 3, "varwidth": 4,
                                  "noseparator": 5}) if p.options.strip() else {}
    stats = _ts_stats(str(o.get("statistics", o.get("stats", "mean"))))
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("aweight", "fweight"))
    if o.get("casewise"):
        for v in vars_:
            base &= S.valid(v.data)
    fmt = str(o.get("format", "%9.0g")).strip() if o.get("format") not in (None, True) else "%9.0g"
    cols_stats = str(o.get("columns", "variables")).strip().lower().startswith("s")
    out = s.output

    def value(stat, v, mask):
        sel = mask & S.valid(v.data)
        return _ts_value(stat, v.data[sel], w[sel] if w is not None else None, wtype)

    def fmtv(x):
        return "." if x >= M.SYSMISS else format_value(x, fmt, pad=False).strip()

    if not o.get("by"):
        if cols_stats:
            heads = [_TS_LABEL.get(st_, st_) for st_ in stats]
            out.write("\n    variable |" + "".join(f"{h:>10}" for h in heads) + "\n", "text")
            out.write("-" * 13 + "+" + "-" * (10 * len(heads)) + "\n", "text")
            for v in vars_:
                out.write(f"{name12(v.name):>12} |", "text")
                out.write("".join(f"{fmtv(value(st_, v, base)):>10}" for st_ in stats) + "\n", "result")
            out.write("-" * (14 + 10 * len(heads)) + "\n", "text")
        else:
            out.write("\n   stats |" + "".join(f"{_abbrev(v.name, 8):>10}" for v in vars_) + "\n", "text")
            out.write("-" * 9 + "+" + "-" * (10 * len(vars_)) + "\n", "text")
            for st_ in stats:
                out.write(f"{_TS_LABEL.get(st_, st_):>8} |", "text")
                out.write("".join(f"{fmtv(value(st_, v, base)):>10}" for v in vars_) + "\n", "result")
            out.write("-" * (10 + 10 * len(vars_)) + "\n", "text")
        return
    # by(): uma linha (ou bloco) por categoria
    from .tabulate import _categories, _codes
    bv = s.data.get(s.expand_varlist(str(o["by"]).strip())[0])
    cats, labels, ok = _categories(s, bv, base, missing=bool(o.get("missing")), nolabel=False)
    codes = _codes(bv, cats)
    W = max([8, len(bv.name)] + [len(x) for x in labels])
    names = ", ".join(_TS_LABEL.get(st_, st_) for st_ in stats)
    if len(stats) > 1 or len(vars_) > 1:
        lead = f"Summary statistics: {names}" if len(stats) > 1 else \
            f"Summary for variables: {' '.join(v.name for v in vars_)}"
        out.write(f"\n{lead}\n", "text")
        label = f" ({bv.label})" if bv.label else ""
        out.write(f"  by categories of: {bv.name}{label}\n", "text")
    out.write("\n" + f"{_abbrev(bv.name, W):>{W}} |" + "".join(
        f"{_abbrev(v.name, 8):>10}" for v in vars_) + "\n", "text")
    sep = "-" * (W + 1) + "+" + "-" * (10 * len(vars_)) + "\n"
    out.write(sep, "text")
    groups = [(labels[k], ok & (codes == k)) for k in range(len(cats))]
    if not o.get("nototal"):
        groups.append(("Total", ok))
    for gi, (label, mask) in enumerate(groups):
        if label == "Total" and gi:
            out.write(sep, "text")
        for k, st_ in enumerate(stats):
            out.write(f"{(label if k == 0 else ''):>{W}} |", "text")
            out.write("".join(f"{fmtv(value(st_, v, mask)):>10}" for v in vars_) + "\n", "result")
        if len(stats) > 1 and gi < len(groups) - 1 and groups[gi + 1][0] != "Total":
            out.write(sep, "text")
    out.write("-" * (W + 2 + 10 * len(vars_)) + "\n", "text")
