"""summarize, tabstat, ci, ttest, prtest, correlate/pwcorr, centile, pctile,
xtile and _pctile ([R] summarize and related).

Layouts reproduced from Stata 14; undocumented details have
the VERIFY mark and cases in compat/do/0301_summarize.do.
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
    """(weight vector or None, type) and a mask without missing/negative weights."""
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
    """Pairs (group mask, index of the 1st obs.) for the by prefix; a single
    pair with everything true when there is no by."""
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
        raise StataError(101, "iweights not allowed with detail")   # VERIFY
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
    """Computes and publishes r() for the variable; returns the moments."""
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
        out.write("no observations\n", "text")   # VERIFY
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


# ---------------------------------------------------------------------------
# correlate / pwcorr
# ---------------------------------------------------------------------------

def _corr_table(s: "Session", names: list[str], cells, *, extra=None, star_p=None) -> None:
    """cells[i][j] (j <= i): correlation text; extra: lines per cell (sig., obs)."""
    out = s.output
    k = len(names)
    for start in range(0, k, 7):          # VERIFY: up to 7 columns per block
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
    cells = [[(format_value(float(mat[i, j]), "%9.4f" if not o.get("covariance") else "%8.0g",
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
    out.write("\n" + " " * 55 + head + "\n", "text")
    out.write(f"    Variable |       Obs  Percentile    Centile        [{level:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 61 + "\n", "text")
    alpha = 1 - level / 100
    held_any = False
    s.r = {}
    for v in vars_:
        xs = np.sort(v.data[mask & S.valid(v.data)])
        n = len(xs)
        for k, c in enumerate(cents):
            q = c / 100
            val = S.centile_value(xs, c)
            (lo, lo_held), (hi, hi_held) = (_binom_limit(xs, q, alpha / 2),
                                            _binom_limit(xs, q, 1 - alpha / 2))
            held = lo_held or hi_held
            held_any |= held
            name = f"{name12(v.name):>12}" if k == 0 else " " * 12
            out.write(name + " |", "text")
            obs = f"{n:>10,}" if k == 0 else " " * 10
            out.write(f"{obs}{c:>11g}   {g9(val)}       {g9(lo)}   {g9(hi)}" + ("*" if held else "") + "\n",
                      "result")
            s.r.update({"N": float(n), f"c_{k + 1}": val, f"lb_{k + 1}": lo, f"ub_{k + 1}": hi})
    if held_any:
        out.write("\n* Lower (upper) confidence limit held at minimum (maximum) of sample\n", "text")


def _binom_limit(xs: np.ndarray, q: float, target: float) -> tuple[float, bool]:
    """Interpolated binomial confidence limit of centile q ([R] centile).
    Finds k with F(k-1) < target <= F(k), F the distribution function of the
    Binomial(n, q), and interpolates between x_(k) and x_(k+1) by the fraction
    (target - F(k-1)) / (F(k) - F(k-1)). Outside the sample, the limit stays
    at the minimum or maximum (marked with * in the output). Checked against
    compat/expected/0301_summarize.log."""
    from scipy import stats as st
    n = len(xs)
    if n == 0:
        return M.SYSMISS, False
    F = st.binom.cdf(np.arange(n + 1), n, q)
    k = int(np.searchsorted(F, target - 1e-12, side="left"))
    if k <= 0:
        return float(xs[0]), True
    if k >= n:
        return float(xs[-1]), True
    f = (target - F[k - 1]) / (F[k] - F[k - 1])
    return float(xs[k - 1] + f * (xs[k] - xs[k - 1])), False


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
        raise StataError(198, "nquantiles() too large for the number of observations")   # VERIFY
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
    # category = 1 + number of cut points strictly less than x
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
            raise StataError(198, f"{w} invalid statistic")   # VERIFY
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
    # by(): one line (or block) per category
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
