"""oneway, anova and tabulate, summarize() ([R] oneway, [R] anova,
[R] tabulate, summarize()).

anova uses partial sums of squares (type III): each term is tested with
the others in the model, in a deviation (sum-to-zero) coding, which is what
Stata reports as Partial SS. The e(b) and e(V) coefficients come from a
regress with the same terms (factors with i.), so that predict, test and
margins work after anova."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail as _chi2_tail
from ..lang.functions import _abbrev
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x: float, fmt: str = "%10.0g") -> str:
    if x is None or x >= SYS or not math.isfinite(x):
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _group_values(s: "Session", name: str, mask: np.ndarray):
    from ..commands.inspect import cell_text
    v = s.data.get(name)
    col = np.asarray(v.raw if v.is_string else v.data)
    vals = sorted(set(col[mask].tolist()), key=lambda x: (isinstance(x, str), x))
    labels = [cell_text(s.data, v, int(np.flatnonzero(mask & (col == x))[0])) for x in vals]
    return col, vals, labels


# ---------------------------------------------------------------------------
# oneway
# ---------------------------------------------------------------------------

def cmd_oneway(s: "Session", args: str) -> None:
    from ..commands._util import touse
    p = parse_standard(args)
    o = match_options(p.options, {"tabulate": 3, "means": 1, "standard": 2, "freq": 1, "obs": 2,
                                  "noanova": 4, "nolabel": 5, "wrap": 4, "missing": 1,
                                  "bonferroni": 1, "scheffe": 1, "sidak": 2})
    names = s.expand_varlist(p.varlist)
    if len(names) != 2:
        raise StataError(102 if len(names) < 2 else 103,
                         "too few variables specified" if len(names) < 2 else "too many variables specified")
    yv, gv = s.data.get(names[0]), s.data.get(names[1])
    if yv.is_string:
        raise StataError(109, "type mismatch")
    mask = touse(s, p) & (yv.data < SYS)
    if not gv.is_string:
        mask &= gv.data < SYS
    from ..commands.summarize import weights
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight"))
    y = yv.data.astype(float)
    col, vals, labels = _group_values(s, gv.name, mask)
    out = s.output
    if o.get("tabulate"):
        _summary_table(s, gv, names[0], [(lab, mask & (col == v)) for v, lab in zip(vals, labels)], mask,
                       y, w, o)
    ww = np.ones(len(y)) if w is None else w
    N = float(ww[mask].sum())
    ybar = float((ww[mask] * y[mask]).sum() / N)
    ssb = ssw = 0.0
    sds = []
    ns = []
    for v in vals:
        sel = mask & (col == v)
        n_g = float(ww[sel].sum())
        m_g = float((ww[sel] * y[sel]).sum() / n_g)
        ssb += n_g * (m_g - ybar) ** 2
        ssw += float((ww[sel] * (y[sel] - m_g) ** 2).sum())
        ns.append(n_g)
        sds.append(float((ww[sel] * (y[sel] - m_g) ** 2).sum() / (n_g - 1)) if n_g > 1 else SYS)
    k = len(vals)
    dfb, dfw = k - 1, N - k
    msb = ssb / dfb if dfb > 0 else SYS
    msw = ssw / dfw if dfw > 0 else SYS
    F = msb / msw if dfb > 0 and msw > 0 and msw < SYS else SYS
    pF = float(st.f.sf(F, dfb, dfw)) if F < SYS else SYS
    tss = ssb + ssw
    if not o.get("noanova"):
        out.write("\n" + " " * 24 + "Analysis of Variance\n", "text")
        out.write("    Source              SS         df      MS            F     Prob > F\n", "text")
        out.write("-" * 72 + "\n", "text")
        out.write(f"{'Between groups':<14}", "text")
        out.write(f"{_g(ssb, '%11.0g'):>16}{dfb:>7.0f}{_g(msb, '%11.0g'):>13}{F:>10.2f}{pF:>11.4f}\n", "result")
        out.write(f"{' Within groups':<14}", "text")
        out.write(f"{_g(ssw, '%11.0g'):>16}{dfw:>7.0f}{_g(msw, '%11.0g'):>13}\n", "result")
        out.write("-" * 72 + "\n", "text")
        out.write(f"{'    Total':<14}", "text")
        out.write(f"{_g(tss, '%11.0g'):>16}{N - 1:>7.0f}{_g(tss / (N - 1), '%11.0g'):>13}\n", "result")
    # Bartlett: only when all groups have 2 or more observations
    chi_b = SYS
    if all(n_g > 1 for n_g in ns) and all(v > 0 for v in sds):
        sp2 = sum((n_g - 1) * v for n_g, v in zip(ns, sds)) / dfw
        num = dfw * math.log(sp2) - sum((n_g - 1) * math.log(v) for n_g, v in zip(ns, sds))
        den = 1 + (sum(1 / (n_g - 1) for n_g in ns) - 1 / dfw) / (3 * (k - 1))
        chi_b = num / den
        pb = float(_chi2_tail(chi_b, k - 1))
        if not o.get("noanova"):
            out.write(f"\nBartlett's test for equal variances:  chi2({k - 1}) = {chi_b:>8.4f}  Prob>chi2 = {pb:.3f}\n",
                      "text")
    # return list order of Stata 14 (compat 0802)
    s.r = dict(reversed([("N", N), ("mss", ssb), ("df_m", float(dfb)), ("rss", ssw), ("df_r", dfw),
                         ("F", F), ("chi2bart", chi_b), ("df_bart", float(k - 1))]))


def _summary_table(s: "Session", gv, yname: str, groups, mask, y, w, o) -> None:
    """Table of means, standard deviations and frequencies by group (oneway, tabulate and
    tabulate, summarize())."""
    show = [k for k in ("means", "standard", "freq", "obs") if o.get(k)]
    if not show:
        show = ["means", "standard", "freq"]
    for k in ("means", "standard", "freq", "obs"):
        if o.get("no" + k) and k in show:
            show.remove(k)
    heads = {"means": "Mean", "standard": "Std. Dev.", "freq": "Freq.", "obs": "Obs."}
    yvar = s.data.get(yname)
    title = f"Summary of {yvar.label or yname}"
    if yvar.label:
        title = f"Summary of {yvar.label}"
    width = 12 * len(show)
    out = s.output
    rowlab = gv.label or gv.name
    words = rowlab.split()
    lines: list[str] = []
    cur = ""
    for wd in words:
        if cur and len(cur) + 1 + len(wd) > 11:
            lines.append(cur)
            cur = wd
        else:
            cur = (cur + " " + wd).strip()
    if cur:
        lines.append(cur)
    lines = [ln[:11] for ln in lines] or [""]
    head_lines = [""] * max(0, 2 - len(lines)) + lines
    out.write("\n", "text")
    tlines = _wrap(title, width)
    nhead = max(len(tlines) + 1, 2, len(lines))
    head_lines = [""] * (nhead - len(lines)) + lines
    for k, ln in enumerate(head_lines):
        right = ""
        ti = k - (nhead - 1 - len(tlines))
        if 0 <= ti < len(tlines):
            pad = width - len(tlines[ti])
            right = (" " * ((pad + 1) // 2) + tlines[ti]).rstrip()
        elif k == nhead - 1:
            right = "".join(f"{heads[c]:>12}" for c in show)
        out.write(f"{ln:>11} |{right}\n", "text")
    out.write("-" * 12 + "+" + "-" * width + "\n", "text")
    ww = np.ones(len(y)) if w is None else w

    def cells(sel):
        x, wx = y[sel], ww[sel]
        n = float(wx.sum())
        m = float((wx * x).sum() / n) if n else SYS
        sd = float(math.sqrt((wx * (x - m) ** 2).sum() / (n - 1))) if n > 1 else SYS
        vals = {"means": _g(m), "standard": _g(sd), "freq": f"{n:,.0f}", "obs": f"{int(sel.sum()):,}"}
        return "".join(f"{vals[c]:>12}" for c in show)

    for lab, sel in groups:
        out.write(f"{_abbrev(str(lab), 11):>11} |", "text")
        out.write(cells(sel) + "\n", "result")
    out.write("-" * 12 + "+" + "-" * width + "\n", "text")
    out.write(f"{'Total':>11} |", "text")
    out.write(cells(mask) + "\n", "result")


def _wrap(text: str, width: int) -> list[str]:
    """Wraps the title in lines of up to `width` characters (compat 0802:
    "Summary of" / "price" over a single column)."""
    if len(text) <= width:
        return [text]
    lines, cur = [], ""
    for w in text.split():
        if cur and len(cur) + 1 + len(w) > width:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    if cur:
        lines.append(cur)
    return lines


def tab_summarize(s: "Session", gv, yname: str, mask: np.ndarray, w, o: dict) -> None:
    """tabulate var1, summarize(var2): one line per value of var1."""
    yv = s.data.get(yname)
    if yv.is_string:
        raise StataError(109, "type mismatch")
    mask = mask & (yv.data < SYS)
    if not gv.is_string and not o.get("missing"):
        mask &= gv.data < SYS
    col, vals, labels = _group_values(s, gv.name, mask)
    if o.get("nolabel"):
        labels = [format_value(float(v), gv.fmt, pad=False).strip() if not isinstance(v, str) else v for v in vals]
    _summary_table(s, gv, yv.name, [(lab, mask & (col == v)) for v, lab in zip(vals, labels)], mask,
                   yv.data.astype(float), w, o)


# ---------------------------------------------------------------------------
# anova
# ---------------------------------------------------------------------------

_TERM = re.compile(r"^(?:(c|i)\.)?([A-Za-z_]\w*)$")


class _Term:
    """A model term: crossed components (a#b) optionally nested in other variables (a|b)."""

    def __init__(self, comps: list[tuple[str, bool]], nest: list[str]):
        self.comps, self.nest = comps, nest

    @property
    def name(self) -> str:
        left = "#".join(n for n, _ in self.comps)
        return left + ("|" + "#".join(self.nest) if self.nest else "")

    @property
    def variables(self) -> list[str]:
        return [n for n, _ in self.comps] + list(self.nest)

    def key(self):
        return (tuple(sorted(self.comps)), tuple(sorted(self.nest)))


def _parse_terms(s: "Session", text: str, continuous: set[str]):
    """Segments of terms separated by "/" (error terms). Each term: _Term."""
    segments: list[list[_Term]] = [[]]
    seen = set()
    for tok in text.replace("/", " / ").split():
        if tok == "/":
            segments.append([])
            continue
        left, bar, right = tok.partition("|")
        nest = [s.expand_varlist(n)[0] for n in right.split("#")] if bar else []
        if "##" in left:
            from itertools import combinations
            comps = [_comp(s, p_, continuous) for p_ in left.split("##")]
            new = [list(c) for size in range(1, len(comps) + 1) for c in combinations(comps, size)]
        else:
            new = [[_comp(s, p_, continuous) for p_ in left.split("#")]]
        for c in new:
            term = _Term(c, nest)
            if term.key() not in seen:
                seen.add(term.key())
                segments[-1].append(term)
    return segments


def _comp(s: "Session", tok: str, continuous: set[str]) -> tuple[str, bool]:
    m = _TERM.match(tok.strip())
    if not m:
        raise StataError(198, f"{tok} invalid term")   # VERIFY
    prefix, name = m.group(1), m.group(2)
    name = s.expand_varlist(name)[0]
    return name, (prefix == "c" or name in continuous)


def _effect_columns(x: np.ndarray, levels: list) -> np.ndarray:
    """Deviation coding: k-1 columns, the last category = -1."""
    k = len(levels)
    cols = np.zeros((len(x), max(k - 1, 0)))
    for j, lv in enumerate(levels[:-1]):
        cols[:, j] = (x == lv).astype(float)
    if k > 1:
        cols[x == levels[-1], :] = -1.0
    return cols


def _crossed(comps, data: dict, rows: np.ndarray) -> np.ndarray:
    mats = []
    for name, cont in comps:
        x, levels = data[name]
        x = x[rows]
        if cont:
            mats.append((x - x.mean())[:, None])
        else:
            mats.append(_effect_columns(x, sorted(set(x.tolist()))))
    out = mats[0]
    for m_ in mats[1:]:
        out = np.einsum("ij,ik->ijk", out, m_).reshape(len(out), -1)
    return out


def _term_matrix(t: _Term, data: dict) -> np.ndarray:
    n = len(next(iter(data.values()))[0])
    if not t.nest:
        return _crossed(t.comps, data, np.arange(n))
    # nested: the coding of the components is built within each cell of the nesting
    # variables, so the term does not repeat their effects
    cells = np.zeros(n, dtype=np.int64)
    for v in t.nest:
        x = data[v][0]
        cells = cells * (len(data[v][1]) + 1) + np.searchsorted(data[v][1], x)
    blocks = []
    for c in np.unique(cells):
        rows = np.flatnonzero(cells == c)
        sub = _crossed(t.comps, data, rows)
        full = np.zeros((n, sub.shape[1]))
        full[rows] = sub
        blocks.append(full)
    return np.column_stack(blocks) if blocks else np.zeros((n, 0))


def _rss(X: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[float, int]:
    sw = np.sqrt(w)
    Xw, yw = X * sw[:, None], y * sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    r = yw - Xw @ beta
    return float(r @ r), int(np.linalg.matrix_rank(Xw))


def cmd_anova(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..commands.summarize import weights
    p = parse_standard(args)
    o = match_options(p.options, {"partial": 4, "sequential": 3, "continuous": 4, "repeated": 6,
                                  "noconstant": 4, "dropemptycells": 4, "category": 3,
                                  "bse": 3, "bseunit": 4, "grouping": 3})
    sequential = bool(o.get("sequential"))
    body = p.varlist.replace("/", " / ")
    words = body.split()
    if len(words) < 2:
        raise StataError(102, "too few variables specified")
    dep = s.expand_varlist(words[0])[0]
    continuous = set(s.expand_varlist(o["continuous"])) if isinstance(o.get("continuous"), str) else set()
    segments = _parse_terms(s, " ".join(words[1:]), continuous)
    terms = [t_ for seg in segments for t_ in seg]
    if not terms:
        raise StataError(102, "too few variables specified")
    allvars = sorted({n for t_ in terms for n in t_.variables})
    mask = touse(s, p)
    yv = s.data.get(dep)
    mask &= yv.data < SYS
    for n in allvars:
        v = s.data.get(n)
        if v.is_string:
            raise StataError(109, "type mismatch")
        mask &= v.data < SYS
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight"))
    y = yv.data[mask].astype(float)
    ww = np.ones(len(y)) if w is None else w[mask].astype(float)
    if wtype == "aweight":
        ww = ww * len(ww) / ww.sum()
    data = {}
    for n in allvars:
        x = s.data.get(n).data[mask].astype(float)
        data[n] = (x, sorted(set(x.tolist())))
    blocks = [_term_matrix(t_, data) for t_ in terms]
    one = np.ones((len(y), 1))
    Xfull = np.column_stack([one] + blocks)
    rss, rank = _rss(Xfull, y, ww)
    N = float(ww.sum())
    ybar = float((ww * y).sum() / N)
    tss = float((ww * (y - ybar) ** 2).sum())
    mss = tss - rss
    df_m = rank - 1
    df_r = N - rank
    mse = rss / df_r if df_r > 0 else SYS
    ss_df = []
    prev_rss, prev_rank = _rss(one, y, ww)
    for k, t_ in enumerate(terms):
        if sequential:
            r_k, rk_k = _rss(np.column_stack([one] + blocks[:k + 1]), y, ww)
            ss_df.append((prev_rss - r_k, rk_k - prev_rank))
            prev_rss, prev_rank = r_k, rk_k
        else:
            Xr = np.column_stack([one] + [b for j, b in enumerate(blocks) if j != k])
            rss_k, rank_k = _rss(Xr, y, ww)
            ss_df.append((rss_k - rss, rank - rank_k))
    index = {id(t_): k for k, t_ in enumerate(terms)}

    def ms_of(k):
        ss, df = ss_df[k]
        return ss / df if df > 0 else SYS

    def test(k, err_ms, err_df):
        ss, df = ss_df[k]
        ms = ms_of(k)
        F = ms / err_ms if df > 0 and err_ms < SYS and err_ms > 0 else SYS
        pv = float(st.f.sf(F, df, err_df)) if F < SYS else SYS
        return F, pv

    # groups: terms of segment i tested against the first term of segment i+1 (or the residual
    # when that segment is empty: a trailing "/"); the last segment against the residual
    groups = []
    tested_against = {}
    for i, seg in enumerate(segments):
        rest = seg if i == 0 else seg[1:]
        if i + 1 < len(segments):
            nxt = segments[i + 1]
            err = nxt[0] if nxt else None
            members = seg
            groups.append((members, err))
            for t_ in members:
                tested_against[index[id(t_)]] = None if err is None else index[id(err)]
        else:
            groups.append((rest, None))
            for t_ in rest:
                tested_against[index[id(t_)]] = None
    rows = []
    for k, t_ in enumerate(terms):
        e = tested_against.get(k, "none")
        if e == "none":
            F, pv = SYS, SYS
        elif e is None:
            F, pv = test(k, mse, df_r)
        else:
            F, pv = test(k, ms_of(e), ss_df[e][1])
        ss, df = ss_df[k]
        rows.append((t_.name, ss, df, ms_of(k), F, pv))
    F_m = (mss / df_m) / mse if df_m > 0 and mse < SYS and mse > 0 else SYS
    p_m = float(st.f.sf(F_m, df_m, df_r)) if F_m < SYS else SYS
    r2 = mss / tss if tss > 0 else SYS
    r2a = 1 - (1 - r2) * (N - 1) / df_r if df_r > 0 and r2 < SYS else SYS
    rmse = math.sqrt(mse) if mse < SYS else SYS

    # e(b) and e(V): regress with the same terms (factors with i.; a|b as i.a#i.b)
    def fv(t_):
        parts = [("c." if c else "i.") + n for n, c in t_.comps] + ["i." + n for n in t_.nest]
        return "#".join(parts)
    rhs = " ".join(fv(t_) for t_ in terms)
    wtxt = f" [{wtype}={p.weight[1]}]" if p.weight else ""
    cond = f" if {p.if_}" if p.if_ else ""
    rng = f" in {p.in_}" if p.in_ else ""
    s.output.quiet_depth += 1
    try:
        s.interp._execute_expanded(f"regress {dep} {rhs}{wtxt}{cond}{rng}")
    finally:
        s.output.quiet_depth -= 1

    out = s.output
    out.write("\n", "text")
    out.write(" " * 25 + "Number of obs = ", "text")
    out.write(f"{N:>10,.0f}", "result")
    out.write("    R-squared     = ", "text")
    out.write(f"{r2:>7.4f}\n", "result")
    out.write(" " * 25 + "Root MSE      = ", "text")
    out.write(f"{_g(rmse, '%8.0g'):>10}", "result")
    out.write("    Adj R-squared = ", "text")
    out.write(f"{r2a:>7.4f}\n", "result")
    head = "   Seq. SS" if sequential else "Partial SS"
    out.write("\n" + f"{'Source':>24} | {head}         df         MS        F    Prob>F\n", "text")
    longest = max([len(r[0]) for r in rows] + [10])
    dash = max(11, longest + 1)
    sep = " " * (25 - dash) + "-" * dash + "+" + "-" * 52 + "\n"
    out.write(sep, "text")

    def line(label, ss, df, ms, F=None, pv=None):
        out.write(f"{_abbrev(label, 24) if len(label) > 24 else label:>24} |", "text")
        tail = f" {_g(ss):>10}{df:>11.0f}{_g(ms):>12}"
        if F is not None:
            tail += f"{F:>10.2f}{pv:>8.4f}" if F < SYS else f"{'.':>10}{'.':>8}"
        out.write(tail + "\n", "result")

    line("Model", mss, df_m, mss / df_m if df_m else SYS, F_m, p_m)
    out.write(f"{'':>24} |\n", "text")
    for gi, (members, err) in enumerate(groups):
        last = gi == len(groups) - 1
        for t_ in members:
            name, ss, df, ms, F, pv = rows[index[id(t_)]]
            line(name, ss, df, ms, F, pv)
        if err is not None:
            k = index[id(err)]
            line(err.name, ss_df[k][0], ss_df[k][1], ms_of(k))
            out.write(sep, "text")
        elif not last:
            out.write(f"{'':>24} |\n", "text")
            line("Residual", rss, df_r, mse)
            out.write(sep, "text")
    if groups[-1][0] or len(groups) == 1:
        out.write(f"{'':>24} |\n", "text")
    line("Residual", rss, df_r, mse)
    out.write(sep, "text")
    line("Total", tss, N - 1, tss / (N - 1))

    repeated = []
    if isinstance(o.get("repeated"), str):
        repeated = s.expand_varlist(o["repeated"])
        _repeated(s, o, terms, ss_df, data, y, ww, repeated, tested_against, index, rss, df_r)

    # e() of anova in Stata 14, in ereturn list order (compat 0802)
    e0 = dict(s.e)
    seen: list[str] = []
    for t_ in terms:
        for n_ in t_.variables:
            if n_ not in seen:
                seen.append(n_)
    scalars = [("N", N), ("df_m", float(df_m)), ("df_r", float(df_r)), ("F", F_m), ("r2", r2),
               ("rmse", rmse), ("mss", mss), ("rss", rss), ("r2_a", r2a), ("ll", e0.get("ll", SYS)),
               ("ll_0", e0.get("ll_0", SYS)), ("rank", e0.get("rank", SYS))]
    for k, (name, ss, df, ms, F, pv) in enumerate(rows, start=1):
        scalars += [(f"ss_{k}", ss), (f"df_{k}", float(df)), (f"F_{k}", F)]
        e = tested_against.get(k - 1)
        if isinstance(e, int):     # VERIFY: names of the error-term results
            scalars += [(f"ssdenom_{k}", ss_df[e][0]), (f"dfdenom_{k}", float(ss_df[e][1]))]
    macros = [("cmdline", f"anova {args.strip()}"), ("depvar", dep), ("cmd", "anova"),
              ("varnames", " ".join(seen)), ("estat_cmd", "anova_estat")]
    macros += [(f"term_{k}", rows[k - 1][0]) for k in range(len(rows), 0, -1)]
    if repeated:
        macros.append(("repvars", " ".join(repeated)))
    macros += [("sstype", "sequential" if sequential else "partial"), ("properties", "b V"),
               ("predict", "regres_p"), ("model", "ols")]
    e: dict = {}
    for k_, v in reversed(macros):
        e[k_] = v
    for k_, v in reversed(scalars):
        e[k_] = float(v)
    for k_ in ("V", "b"):
        if k_ in e0:
            e[k_] = e0[k_]
    for k_, v in getattr(s, "_anova_extra", {}).items():
        e[k_] = v
    s._anova_extra = {}
    s.e = e
    est = getattr(s, "_estimates", None)
    if est is not None:
        est.cmd = "anova"


def _cells(data: dict, names: list[str]) -> np.ndarray:
    n = len(next(iter(data.values()))[0])
    ids = np.zeros(n, dtype=np.int64)
    for v in names:
        ids = ids * (len(data[v][1]) + 1) + np.searchsorted(data[v][1], data[v][0])
    return np.unique(ids, return_inverse=True)[1]


def _repeated(s, o, terms, ss_df, data, y, w, repeated, tested_against, index, rss, df_r) -> None:
    """repeated(): Huynh-Feldt, Greenhouse-Geisser and Box's conservative corrections of the
    tests of the repeated variables ([R] anova, "Repeated-measures ANOVA"). The covariance of
    the repeated measures is computed on the subject (between-subjects error unit) by level
    matrix, pooled within the between-subjects groups (VERIFY: details for designs with more
    than one repeated variable)."""
    out = s.output
    # between-subjects error term: given by bse(), or the term without repeated variables that,
    # together with them, identifies the observations with the fewest levels
    bse = None
    if isinstance(o.get("bse"), str):
        want = o["bse"].strip().replace(" ", "")
        bse = next((t_ for t_ in terms if t_.name == want), None)
        if bse is None:
            raise StataError(198, f"bse({o['bse']}) not a term in the model")   # VERIFY
    else:
        best = None
        for t_ in terms:
            if set(t_.variables) & set(repeated) or any(c for _, c in t_.comps):
                continue
            ids = _cells(data, t_.variables + repeated)
            if len(np.unique(ids)) != len(ids):
                continue
            nlev = len(np.unique(_cells(data, t_.variables)))
            if best is None or nlev < best[0]:
                best = (nlev, t_)
        if best is None:
            raise StataError(421, "could not determine between-subjects error term; use bse() "
                                  "option")   # VERIFY
        bse = best[1]
    unit_var = (s.expand_varlist(o["bseunit"])[0] if isinstance(o.get("bseunit"), str)
                else bse.comps[0][0])
    unit = _cells(data, bse.variables)
    nunits = int(unit.max()) + 1
    k_bse = index[id(bse)]
    out.write(f"\nBetween-subjects error term:  {bse.name}\n", "text")
    out.write(f"{'Levels:':>28}  ", "text")
    out.write(f"{nunits:<10}", "result")
    out.write(f"({ss_df[k_bse][1]:.0f} df)\n", "text")
    out.write(f"{'Lowest b.s.e. variable:':>28}  {unit_var}\n", "text")
    # between-subjects groups: model variables constant within each subject
    if isinstance(o.get("grouping"), str):
        gvars = s.expand_varlist(o["grouping"])
    else:
        gvars = []
        for v in sorted({n for t_ in terms for n in t_.variables}):
            if v in repeated or v == unit_var:
                continue
            x = data[v][0]
            if all(len(set(x[unit == u].tolist())) == 1 for u in range(nunits)):
                gvars.append(v)
    if gvars:
        out.write(f"{'Covariance pooled over:':>28}  {' '.join(gvars)}   (for repeated variable)\n",
                  "text")
    group = _cells(data, gvars) if gvars else np.zeros(len(y), dtype=np.int64)
    extra = {}
    for r_i, rv in enumerate(repeated, start=1):
        levels = data[rv][1]
        kl = len(levels)
        lev = np.searchsorted(levels, data[rv][0])
        Y = np.full((nunits, kl), np.nan)
        sums = np.zeros((nunits, kl))
        cnt = np.zeros((nunits, kl))
        np.add.at(sums, (unit, lev), w * y)
        np.add.at(cnt, (unit, lev), w)
        with np.errstate(invalid="ignore", divide="ignore"):
            Y = sums / cnt
        ok = ~np.isnan(Y).any(axis=1)
        ug = np.zeros(nunits, dtype=np.int64)
        ug[unit] = group
        Yc = Y[ok]
        gg = ug[ok]
        D = np.zeros_like(Yc)
        for gcode in np.unique(gg):
            m_ = gg == gcode
            D[m_] = Yc[m_] - Yc[m_].mean(axis=0)
        n_sub, n_grp = int(ok.sum()), len(np.unique(gg))
        S_ = D.T @ D / max(n_sub - n_grp, 1)
        # orthonormal contrasts orthogonal to the constant
        C = np.linalg.qr(np.column_stack([np.ones(kl), np.eye(kl)[:, :kl - 1]]))[0][:, 1:].T
        Mx = C @ S_ @ C.T
        km1 = kl - 1
        gg_eps = float(np.trace(Mx) ** 2 / (km1 * np.trace(Mx @ Mx))) if km1 > 0 else 1.0
        hf = (n_sub * km1 * gg_eps - 2) / (km1 * (n_sub - n_grp - km1 * gg_eps))
        box = 1 / km1 if km1 > 0 else 1.0
        out.write(f"\nRepeated variable: {rv}\n", "text")
        pad = " " * 42
        out.write(f"{pad}Huynh-Feldt epsilon        = ", "text")
        out.write(f"{hf:7.4f}\n", "result")
        if hf > 1:
            out.write(f"{pad}*Huynh-Feldt epsilon reset to 1.0000\n", "text")
        hf_used = min(hf, 1.0)
        out.write(f"{pad}Greenhouse-Geisser epsilon = ", "text")
        out.write(f"{gg_eps:7.4f}\n", "result")
        out.write(f"{pad}Box's conservative epsilon = ", "text")
        out.write(f"{box:7.4f}\n", "result")
        out.write("\n" + " " * 44 + "------------ Prob > F ------------\n", "text")
        out.write(f"{'Source':>24} |     df      F    Regular    H-F      G-G      Box\n", "text")
        sep = " " * 14 + "-" * 11 + "+" + "-" * 52 + "\n"
        out.write(sep, "text")
        shown_err = set()
        for k, t_ in enumerate(terms):
            if rv not in t_.variables:
                continue
            ss, df = ss_df[k]
            e = tested_against.get(k)
            if isinstance(e, int):
                ems, edf, ename = ss_df[e][0] / ss_df[e][1], ss_df[e][1], terms[e].name
            else:
                ems, edf, ename = rss / df_r, df_r, "Residual"
            F = (ss / df) / ems if df > 0 and ems > 0 else SYS
            ps = [float(st.f.sf(F, df * eps, edf * eps)) if F < SYS else SYS
                  for eps in (1.0, hf_used, gg_eps, box)]
            out.write(f"{t_.name:>24} |", "text")
            out.write(f"{df:>7.0f}{F:>9.2f}" + "".join(f"{x:>9.4f}" for x in ps) + "\n", "result")
            shown_err.add((ename, edf))
        for ename, edf in shown_err:
            out.write(f"{ename:>24} |", "text")
            out.write(f"{edf:>7.0f}\n", "result")
        out.write(sep, "text")
        extra.update({f"gg_{r_i}": gg_eps, f"hf_{r_i}": hf_used,
                      f"box_{r_i}": box})   # VERIFY: names of these e() results
    s._anova_extra = extra
