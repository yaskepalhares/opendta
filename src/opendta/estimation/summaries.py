"""mean, proportion, ratio and total ([R] mean, [R] proportion, [R] ratio, [R] total), with or
without the svy prefix.

Each estimate is a weighted statistic of one variable in one over() group (a "cell") with its
linearized score. What differs between `mean x` and `svy: mean x` is only who supplies the
weights, turns the scores into a variance and writes the header and e(): an engine. SRS (here)
serves the plain commands; svy.SurveyEngine serves the svy prefix.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Protocol

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lazy import lazy
from .results import Estimates, g9, post

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


@dataclass(frozen=True)
class StatInfo:
    title: str      # header title ("Mean estimation")
    column: str     # column of the table and e(depvar) under svy ("Mean")


STATS = {"mean": StatInfo("Mean estimation", "Mean"),
         "proportion": StatInfo("Proportion estimation", "Proportion"),
         "total": StatInfo("Total estimation", "Total"),
         "ratio": StatInfo("Ratio estimation", "Ratio")}


# ---------------------------------------------------------------------------
# the request: variables, over() and the estimation sample
# ---------------------------------------------------------------------------

@dataclass
class Request:
    stat: str
    args: str
    p: Parsed
    o: dict
    level: float
    names: list[str]       # one per estimated variable (ratio: the names of the ratios)
    num: list[str]         # the variables (ratio: the numerators)
    den: list[str]         # ratio: the denominators
    over: str
    mask: np.ndarray       # if/in, without missing values in the variables and in over()

    @property
    def info(self) -> StatInfo:
        return STATS[self.stat]

    @property
    def allv(self) -> list[str]:
        return self.num + self.den


def parse_request(s: "Session", stat: str, args: str) -> Request:
    from ..commands._util import touse
    from ..core.varlist import expand
    p = parse_standard(args)
    o = match_options(p.options, {"over": 2, "level": 1, "vce": 3, "nolegend": 4, "noheader": 6,
                                  "missing": 4, "percent": 3, "citype": 3}) if p.options.strip() else {}
    ds = s.data
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    mask = touse(s, p)
    if stat == "ratio":
        pairs = re.findall(r"(?:\(?\s*(\w+)\s*:)?\s*(\w+)\s*/\s*(\w+)\s*\)?", p.varlist)
        if not pairs:
            raise StataError(198, "invalid syntax")
        names = [nm or f"_ratio_{k + 1}" for k, (nm, _a, _b) in enumerate(pairs)]
        num = [expand(ds, a)[0] for _, a, _b in pairs]
        den = [expand(ds, b)[0] for _, _a, b in pairs]
    else:
        names = expand(ds, p.varlist, allow_empty=False)
        num, den = list(names), []
    over = str(o.get("over") or "").split(",")[0].strip()
    for n in num + den + (expand(ds, over) if over else []):
        mask &= ds.get(n).data < SYS
    return Request(stat, args, p, o, level, names, num, den, over, mask)


# ---------------------------------------------------------------------------
# cells and their estimates
# ---------------------------------------------------------------------------

def _value_label(ds, var, v: float) -> str:
    lab = ds.value_labels.get(var.value_label, {}).get(int(v)) if var.value_label else None
    return lab or (f"{int(v)}" if v == int(v) else f"{v:g}")


def _over_groups(s: "Session", over: str, rows: np.ndarray) -> list[tuple[str, np.ndarray, tuple]]:
    """(label, selection within rows, key) of each over() group; one group without over()."""
    n = int(rows.sum())
    if not over:
        return [("", np.ones(n, bool), ())]
    from ..core.varlist import expand
    ds = s.data
    names = expand(ds, over)
    cols = [ds.get(nm).data[rows] for nm in names]
    groups = []
    for key in sorted(set(zip(*[c.tolist() for c in cols]))):
        sel = np.ones(n, bool)
        for c, k_ in zip(cols, key):
            sel &= c == k_
        labs = [_value_label(ds, ds.get(nm), k_) for nm, k_ in zip(names, key)]
        groups.append((" ".join(labs), sel, key))
    return groups


def _legend(over: str, groups) -> list[str]:
    """over() legend: "         masc: sexo = masc" (compat 0601)."""
    if not over:
        return []
    ovars = over.split()
    legend, seen = [], set()
    for gname, _sel, _key in groups:
        if gname in seen:
            continue
        seen.add(gname)
        desc = ", ".join(f"{v} = {lab}" for v, lab in zip(ovars, gname.split(" "))) \
            if len(ovars) > 1 else f"{ovars[0]} = {gname}"
        legend.append(f"{gname[:13]:>13}: {desc}")
    return legend


@dataclass
class Cell:
    label: tuple[str, str]       # (variable, over() group or category)
    y: np.ndarray
    sel: np.ndarray              # the over() group
    x: np.ndarray | None = None  # ratio: the denominator


def _cells(s: "Session", req: Request, rows: np.ndarray, groups) -> list[Cell]:
    """proportion: one indicator per category of each variable; the others: the variables.
    Missing values (outside the estimation sample) count as 0: their weight is 0."""
    ds = s.data
    if req.stat == "proportion":
        out = []
        for n in req.names:
            v = ds.get(n)
            x = v.data[rows]
            for lv in sorted(set(x.tolist())):
                y, cat = (x == lv).astype(float), _value_label(ds, v, lv)
                out += [Cell((n, cat + (f" {g}" if req.over else "")), y, sel) for g, sel, _ in groups]
        return out

    def values(name):
        x = ds.get(name).data[rows].astype(float)
        return np.where(x < SYS, x, 0.0)
    xs = [values(n) for n in req.den] or [None] * len(req.num)
    return [Cell((name, g), values(n), sel, x)
            for name, n, x in zip(req.names, req.num, xs) for g, sel, _ in groups]


def estimate(stat: str, cells: list[Cell], w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Value and linearized score (n x k) of each cell; proportion is the mean of an indicator."""
    vals, scores = [], []
    n = len(w)
    for c in cells:
        ws = w * c.sel
        if stat == "total":
            m, z = float(ws @ c.y), ws * c.y
        else:
            W = float(ws @ c.x) if stat == "ratio" else ws.sum()
            x = c.x if stat == "ratio" else 1.0
            ok = W != 0 if stat == "ratio" else W > 0
            m = float(ws @ c.y) / W if ok else SYS
            z = ws * (c.y - m * x) / W if ok else np.zeros(n)
        vals.append(m)
        scores.append(z)
    return np.array(vals), np.column_stack(scores) if scores else np.zeros((n, 0))


@dataclass
class Fit:
    cells: list[Cell]
    groups: list
    legend: list[str]
    vals: np.ndarray
    Z: np.ndarray

    @property
    def labels(self) -> list[tuple[str, str]]:
        return [c.label for c in self.cells]

    @property
    def colnames(self) -> list[str]:
        return [f"{a}:{b}" if b else a for a, b in self.labels]


# ---------------------------------------------------------------------------
# engines: weights, variance, header and e()
# ---------------------------------------------------------------------------

class Engine(Protocol):
    rows: np.ndarray     # observations whose scores enter the variance
    w: np.ndarray        # their weights (0 outside the estimation sample)
    n_obs: int

    def variance(self, req: Request, fit: Fit, refit: Callable) -> tuple[np.ndarray, float | None, str]: ...
    def header(self, s: "Session", req: Request, df) -> None: ...
    def post(self, s: "Session", est: Estimates, req: Request, fit: Fit, df) -> None: ...


class SRS:
    """The plain commands: simple random sample, optionally weighted."""

    def __init__(self, s: "Session", req: Request):
        rows, w = req.mask, None
        if req.p.weight:
            from ..commands.summarize import weights
            w_all, _wt, rows = weights(s, req.p, rows, ("pweight", "aweight", "fweight", "iweight"))
            w = w_all[rows]
        self.wtype = req.p.weight[0] if req.p.weight else ""
        self.rows = rows
        self.w = (w if w is not None else np.ones(int(rows.sum()))).astype(float)
        self.n_obs = int(self.w.sum()) if self.wtype == "fweight" else int(rows.sum())

    def variance(self, req, fit, refit):
        # V = N/(N-1) Σ z z' (with fweight, each score counts w times)
        Z, w = fit.Z, self.w
        robust = str(req.o.get("vce", "")).startswith("r")
        fw = self.wtype == "fweight"
        cnt = w if fw else np.ones(len(w))
        nn = float(cnt.sum())
        Zr = Z / np.where(w > 0, w, 1)[:, None] if fw else Z
        if req.stat in ("total", "ratio") or robust or self.wtype == "pweight":
            # linearization: centered scores, V = n/(n-1) Σ (z - z̄)(z - z̄)'; with pweight the
            # standard error is linearized, but Stata shows no label
            zbar = (cnt @ Zr) / nn if req.stat == "total" else np.zeros(Z.shape[1])
            Zc = Zr - zbar
            V = nn / (nn - 1) * (Zc * cnt[:, None]).T @ Zc
        else:
            # "analytic" mean/proportion: s²/n within each over() group, with no covariance between
            # groups (compat 0601: .1209974 = sd/sqrt(n_g))
            V = (Zr * cnt[:, None]).T @ Zr
            Wg = np.array([float(cnt @ c.sel) for c in fit.cells])
            c = np.sqrt(np.where(Wg > 1, Wg / (Wg - 1), np.nan))
            V = V * np.outer(c, c)
        vcetype = "Linearized" if req.stat == "ratio" else "Robust" if robust else ""
        return V, int(nn) - 1, vcetype

    def header(self, s, req, df):
        out = s.output
        out.write("\n", "text")
        out.write(f"{req.info.title:<34}{'Number of obs':<16}= ", "text")
        out.write(f"{self.n_obs:>10,}\n", "result")

    def post(self, s, est, req, fit, df):
        scalars = [("N", float(self.n_obs)), ("df_r", float(df) if df is not None else SYS)]
        macros = [("cmd", req.stat), ("cmdline", f"{req.stat} {req.args.strip()}".strip()),
                  ("properties", "b V"), ("varlist", " ".join(req.allv))]
        post(s, est, scalars, macros, self.rows.copy())


# ---------------------------------------------------------------------------
# the command
# ---------------------------------------------------------------------------

def summary_stat(s: "Session", stat: str, args: str,
                 engine: Callable[["Session", Request], Engine] = SRS) -> None:
    req = parse_request(s, stat, args)
    eng = engine(s, req)
    groups = _over_groups(s, req.over, eng.rows)
    cells = _cells(s, req, eng.rows, groups)
    vals, Z = estimate(stat, cells, eng.w)
    fit = Fit(cells, groups, _legend(req.over, groups), vals, Z)
    V, df, vcetype = eng.variance(req, fit, lambda w: estimate(stat, cells, w)[0])
    se = np.sqrt(np.maximum(np.diag(V), 0)) if V.size else np.zeros(0)
    if np.isnan(V).any():
        se = np.full(len(vals), np.nan)
    out = s.output
    if not req.o.get("noheader"):
        eng.header(s, req, df)
        if req.over:
            out.write("\n" + "".join(ln + "\n" for ln in fit.legend), "text")
        if stat == "ratio":
            out.write("\n" + "".join(f"{rn[:13]:>13}: {a}/{b}\n" for rn, a, b in zip(req.names, req.num, req.den)),
                      "text")
        out.write("\n", "text")

    def table(ss, *_a, **_kw):
        summary_table(ss, req.info.column, fit.labels, vals, se, df, req.level, vcetype, req.over)
    table(s)
    est = Estimates(stat, "", fit.colnames, vals, V, eng.n_obs, "t", df)
    est.rows = []
    est.level = req.level
    est.display = table
    eng.post(s, est, req, fit, df)


def summary_table(s: "Session", column: str, labels, vals, se, df, level, vcetype: str, over: str) -> None:
    out = s.output
    lv = f"{level:g}"
    width = 62
    q = 1 - (1 - level / 100) / 2
    crit = st.t.ppf(q, df) if df and df > 0 else st.norm.ppf(q)
    out.write("-" * width + "\n", "text")
    if vcetype == "Bootstrap" and not (df and df > 0):
        # bootstrap replicate weights: observed value and normal-based interval (compat 0603)
        out.write(f"{'':>12} |   Observed   Bootstrap         Normal-based\n", "text")
    elif vcetype:
        start = 28 + (9 - len(vcetype)) // 2
        out.write(f"{'':>12} |" + " " * (start - 14) + vcetype + "\n", "text")
    head = "Over" if over else ""
    out.write(f"{head:>12} | {column:>10}   Std. Err.     [{lv}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 48 + "\n", "text")
    last = None
    for (var, grp), b, e in zip(labels, vals, se):
        if grp:
            if var != last:
                # variable header on the left; between over() blocks, a separator line (compat
                # 0601)
                if last is not None and over:
                    out.write("-" * 13 + "+" + "-" * 48 + "\n", "text")
                out.write(f"{var[:12]:<12} |" + "\n", "text")
                last = var
            lab = grp
        else:
            lab = var
        out.write(f"{lab[:12]:>12} |", "text")
        if e is None or not np.isfinite(e):
            out.write(f"  {g9(b):>9}  {'.':>9}     {'.':>9}   {'.':>9}\n", "result")
            continue
        if column == "Proportion" and 0 < b < 1 and vcetype in ("Linearized", ""):
            # CI on the logit scale (proportion's default in Stata 14, VERIFY)
            half = crit * e / (b * (1 - b))
            lo, hi = special.expit(special.logit(b) - half), special.expit(special.logit(b) + half)
        else:
            lo, hi = b - crit * e, b + crit * e
        out.write(f"  {g9(b):>9}  {g9(e):>9}     {g9(lo):>9}   {g9(hi):>9}\n", "result")
    out.write("-" * width + "\n", "text")


def cmd_mean(s: "Session", args: str) -> None:
    summary_stat(s, "mean", args)


def cmd_proportion(s: "Session", args: str) -> None:
    summary_stat(s, "proportion", args)


def cmd_total(s: "Session", args: str) -> None:
    summary_stat(s, "total", args)


def cmd_ratio(s: "Session", args: str) -> None:
    summary_stat(s, "ratio", args)
