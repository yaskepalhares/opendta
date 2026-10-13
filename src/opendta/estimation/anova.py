"""oneway, anova e tabulate, summarize() ([R] oneway, [R] anova,
[R] tabulate, summarize()).

anova usa somas de quadrados parciais (tipo III): cada termo é testado com
os demais no modelo, numa codificação por desvios (soma zero), que é o que
o Stata reporta como Partial SS. Os coeficientes e(b) e e(V) vêm de um
regress com os mesmos termos (fatores com i.), para que predict, test e
margins funcionem depois de anova."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.functions import _abbrev
from ..lang.syntax import match_options, parse_standard

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
        out.write(f"{_g(ssb):>16}{dfb:>7.0f}{_g(msb):>13}{F:>10.2f}{pF:>11.4f}\n", "result")
        out.write(f"{' Within groups':<14}", "text")
        out.write(f"{_g(ssw):>16}{dfw:>7.0f}{_g(msw):>13}\n", "result")
        out.write("-" * 72 + "\n", "text")
        out.write(f"{'    Total':<14}", "text")
        out.write(f"{_g(tss):>16}{N - 1:>7.0f}{_g(tss / (N - 1)):>13}\n", "result")
    # Bartlett: só com todos os grupos com 2 ou mais observações
    chi_b = SYS
    if all(n_g > 1 for n_g in ns) and all(v > 0 for v in sds):
        sp2 = sum((n_g - 1) * v for n_g, v in zip(ns, sds)) / dfw
        num = dfw * math.log(sp2) - sum((n_g - 1) * math.log(v) for n_g, v in zip(ns, sds))
        den = 1 + (sum(1 / (n_g - 1) for n_g in ns) - 1 / dfw) / (3 * (k - 1))
        chi_b = num / den
        pb = float(st.chi2.sf(chi_b, k - 1))
        if not o.get("noanova"):
            out.write(f"\nBartlett's test for equal variances:  chi2({k - 1}) = {chi_b:>8.4f}  Prob>chi2 = {pb:.3f}\n",
                      "text")
    s.r = {"N": N, "df_r": dfw, "F": F, "df_m": float(dfb), "rss": ssw, "mss": ssb,
           "chi2bart": chi_b, "df_bart": float(k - 1)}   # VERIFICAR ordem


def _summary_table(s: "Session", gv, yname: str, groups, mask, y, w, o) -> None:
    """Tabela de médias, desvios e frequências por grupo (oneway, tabulate e
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
    for k, ln in enumerate(head_lines):
        right = ""
        if k == len(head_lines) - 2:
            right = title.center(width).rstrip()
        elif k == len(head_lines) - 1:
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


def tab_summarize(s: "Session", gv, yname: str, mask: np.ndarray, w, o: dict) -> None:
    """tabulate var1, summarize(var2): uma linha por valor de var1."""
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


def _parse_terms(s: "Session", text: str, continuous: set[str]):
    """Lista de termos: cada termo é uma lista de (variável, contínua?)."""
    terms: list[list[tuple[str, bool]]] = []
    for tok in text.split():
        if "##" in tok:
            parts = tok.split("##")
            comps = [_comp(s, p_, continuous) for p_ in parts]
            k = len(comps)
            for size in range(1, k + 1):
                from itertools import combinations
                for combo in combinations(comps, size):
                    terms.append(list(combo))
        elif "#" in tok:
            terms.append([_comp(s, p_, continuous) for p_ in tok.split("#")])
        else:
            terms.append([_comp(s, tok, continuous)])
    seen, out = set(), []
    for t in terms:
        key = tuple(sorted(t))
        if key not in seen:
            seen.add(key)
            out.append(t)
    return out


def _comp(s: "Session", tok: str, continuous: set[str]) -> tuple[str, bool]:
    m = _TERM.match(tok.strip())
    if not m:
        raise StataError(198, f"{tok} invalid term")   # VERIFICAR
    prefix, name = m.group(1), m.group(2)
    name = s.expand_varlist(name)[0]
    return name, (prefix == "c" or name in continuous)


def _effect_columns(x: np.ndarray, levels: list) -> np.ndarray:
    """Codificação por desvios: k-1 colunas, a última categoria = -1."""
    k = len(levels)
    cols = np.zeros((len(x), max(k - 1, 0)))
    for j, lv in enumerate(levels[:-1]):
        cols[:, j] = (x == lv).astype(float)
    if k > 1:
        cols[x == levels[-1], :] = -1.0
    return cols


def _term_matrix(t, data: dict) -> np.ndarray:
    mats = []
    for name, cont in t:
        x, levels = data[name]
        if cont:
            mats.append((x - x.mean())[:, None])
        else:
            mats.append(_effect_columns(x, levels))
    out = mats[0]
    for m_ in mats[1:]:
        out = np.einsum("ij,ik->ijk", out, m_).reshape(len(out), -1)
    return out


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
                                  "noconstant": 4, "dropemptycells": 4, "category": 3})
    if o.get("sequential"):
        raise StataError(198, "option sequential not implemented")   # VERIFICAR
    words = p.varlist.split()
    if len(words) < 2:
        raise StataError(102, "too few variables specified")
    dep = s.expand_varlist(words[0])[0]
    continuous = set(s.expand_varlist(o["continuous"])) if isinstance(o.get("continuous"), str) else set()
    terms = _parse_terms(s, " ".join(words[1:]), continuous)
    allvars = sorted({n for t in terms for n, _ in t})
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
    blocks = [_term_matrix(t, data) for t in terms]
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
    rows = []
    for k, t in enumerate(terms):
        Xr = np.column_stack([one] + [b for j, b in enumerate(blocks) if j != k])
        rss_k, rank_k = _rss(Xr, y, ww)
        ss = rss_k - rss
        df = rank - rank_k
        ms = ss / df if df > 0 else SYS
        F = ms / mse if df > 0 and mse < SYS and mse > 0 else SYS
        pv = float(st.f.sf(F, df, df_r)) if F < SYS else SYS
        name = "#".join(n for n, _ in t)
        rows.append((name, ss, df, ms, F, pv))
    F_m = (mss / df_m) / mse if df_m > 0 and mse < SYS and mse > 0 else SYS
    p_m = float(st.f.sf(F_m, df_m, df_r)) if F_m < SYS else SYS
    r2 = mss / tss if tss > 0 else SYS
    r2a = 1 - (1 - r2) * (N - 1) / df_r if df_r > 0 and r2 < SYS else SYS
    rmse = math.sqrt(mse) if mse < SYS else SYS

    # e(b) e e(V): regress com os mesmos termos (fatores com i.)
    rhs = " ".join("#".join(("c." if c else "i.") + n for n, c in t) for t in terms)
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
    out.write(f"{_g(rmse, '%10.0g'):>10}", "result")
    out.write("    Adj R-squared = ", "text")
    out.write(f"{r2a:>7.4f}\n", "result")
    out.write("\n" + f"{'Source':>24} | Partial SS         df         MS        F    Prob>F\n", "text")
    sep = " " * 14 + "-" * 11 + "+" + "-" * 52 + "\n"
    out.write(sep, "text")

    def line(label, ss, df, ms, F=None, pv=None):
        out.write(f"{_abbrev(label, 24) if len(label) > 24 else label:>24} |", "text")
        tail = f" {_g(ss):>10}{df:>11.0f}{_g(ms):>12}"
        if F is not None:
            tail += f"{F:>10.2f}{pv:>8.4f}" if F < SYS else f"{'.':>10}{'.':>8}"
        out.write(tail + "\n", "result")

    line("Model", mss, df_m, mss / df_m if df_m else SYS, F_m, p_m)
    out.write(f"{'':>24} |\n", "text")
    for name, ss, df, ms, F, pv in rows:
        line(name, ss, df, ms, F, pv)
    out.write(f"{'':>24} |\n", "text")
    line("Residual", rss, df_r, mse)
    out.write(sep, "text")
    line("Total", tss, N - 1, tss / (N - 1))

    e = dict(s.e)
    e["cmd"] = "anova"
    e["varnames"] = " ".join(allvars)
    e["depvar"] = dep
    for k, (name, ss, df, ms, F, pv) in enumerate(rows, start=1):
        e[f"term_{k}"] = name
        e[f"ss_{k}"] = ss
        e[f"df_{k}"] = float(df)
        e[f"F_{k}"] = F
    e.update({"N": N, "mss": mss, "rss": rss, "df_m": float(df_m), "df_r": float(df_r), "F": F_m,
              "r2": r2, "r2_a": r2a, "rmse": rmse})
    s.e = e
