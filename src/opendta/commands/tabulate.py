"""tabulate (uma e duas vias), tab1, tab2 e table ([R] tabulate, [R] table).

Layouts reproduzidos a partir do Stata 14; detalhes não documentados têm a
marca VERIFICAR e casos em compat/do/0302_tabulate.do.
"""

from __future__ import annotations

import math
import re
import textwrap
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import stats as S
from ..core.dataset import Variable
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .data import _by_header
from .registry import command
from .summarize import by_groups, g9, weights

if TYPE_CHECKING:
    from ..session import Session


# ---------------------------------------------------------------------------
# categorias
# ---------------------------------------------------------------------------

def _categories(s: "Session", v: Variable, mask: np.ndarray, *, missing: bool, nolabel: bool):
    """Valores distintos (ordenados, missing no fim) e os rótulos exibidos."""
    if v.is_string:
        vals = v.raw
        ok = mask if missing else mask & np.array([x != "" for x in vals], dtype=bool)
        cats = sorted(set(vals[ok].tolist()))
        if "" in cats:
            cats.remove("")
            cats.append("")
        labels = [c for c in cats]
        return cats, labels, ok
    x = v.data
    ok = mask if missing else mask & S.valid(x)
    cats = sorted(set(x[ok].tolist()))
    lab = s.data.value_labels.get(v.value_label, {}) if v.value_label and not nolabel else {}
    labels = []
    for c in cats:
        if c >= M.SYSMISS:
            labels.append(M.missing_name(c))
        elif c == int(c) and int(c) in lab:
            labels.append(lab[int(c)])
        else:
            labels.append(format_value(c, v.fmt, pad=False).strip())
    return cats, labels, ok


def _codes(v: Variable, cats: list) -> np.ndarray:
    """Índice da categoria de cada observação (-1 fora)."""
    lookup = {c: k for k, c in enumerate(cats)}
    src = v.raw.tolist() if v.is_string else v.data.tolist()
    return np.array([lookup.get(x, -1) for x in src], dtype=np.int64)


def _header_lines(text: str, width: int) -> list[str]:
    """Rótulo da variável quebrado em linhas da largura da 1ª coluna."""
    lines = textwrap.wrap(text, width) or [""]
    return lines[-2:] if len(lines) > 2 else lines   # VERIFICAR: no máximo duas linhas


def _fmt_freq(x: float) -> str:
    return f"{x:,.0f}" if x == int(x) else g9(x).strip()


# ---------------------------------------------------------------------------
# uma via
# ---------------------------------------------------------------------------

def _oneway(s: "Session", v: Variable, mask, w, wtype, o: dict) -> None:
    out = s.output
    cats, labels, ok = _categories(s, v, mask, missing=bool(o.get("missing")),
                                   nolabel=bool(o.get("nolabel")))
    codes = _codes(v, cats)
    ww = w if w is not None else np.ones(s.data.nobs)
    freq = np.zeros(len(cats))
    sel = ok & (codes >= 0)
    np.add.at(freq, codes[sel], ww[sel])
    total = float(freq.sum())
    if total == 0:
        out.write("no observations\n", "text")
        s.r = {"N": 0.0, "r": 0.0}
        return
    order = list(range(len(cats)))
    if o.get("sort"):
        order.sort(key=lambda k: -freq[k])
    title = v.label or v.name
    W = max([10] + [len(lb) for lb in labels])
    W = min(W, 30)
    out.write("\n", "text")
    head = _header_lines(title, W)
    for k, line in enumerate(head):
        last = k == len(head) - 1
        out.write(f"{line:>{W}} |" + ("      Freq.     Percent        Cum." if last else "") + "\n", "text")
    sep = "-" * (W + 1) + "+" + "-" * 35 + "\n"
    out.write(sep, "text")
    cum = 0.0
    for k in order:
        cum += freq[k]
        out.write(f"{labels[k][:W]:>{W}} |", "text")
        out.write(f"{_fmt_freq(freq[k]):>11}{100 * freq[k] / total:>12.2f}{100 * cum / total:>12.2f}\n",
                  "result")
    out.write(sep, "text")
    out.write(f"{'Total':>{W}} |", "text")
    out.write(f"{_fmt_freq(total):>11}{100.0:>12.2f}\n", "result")
    s.r = {"N": total, "r": float(len(cats))}
    _store_mats(s, o, cell=freq[:, None], rows=cats, v=v)
    if o.get("generate"):
        stub = str(o["generate"]).strip()
        for k, c in enumerate(cats, start=1):
            col = np.where(codes == k - 1, 1.0, 0.0)
            col[~(mask & (S.valid(v.data) if not v.is_string else True))] = M.SYSMISS
            nv = Variable(f"{stub}{k}", "byte", col)
            nv.label = f"{v.name}=={labels[k - 1]}"
            s.data.add(nv)
        s.notify_state()


def _store_mats(s: "Session", o: dict, *, cell, rows, v, cols=None) -> None:
    from .matrix import Matrix, _store
    if o.get("matcell"):
        _store(s)[str(o["matcell"]).strip()] = Matrix(np.asarray(cell, dtype=np.float64))
    if o.get("matrow") and not v.is_string:
        _store(s)[str(o["matrow"]).strip()] = Matrix(np.asarray(rows, dtype=np.float64)[:, None])
    if o.get("matcol") and cols is not None:
        _store(s)[str(o["matcol"]).strip()] = Matrix(np.asarray(cols, dtype=np.float64)[None, :])


# ---------------------------------------------------------------------------
# duas vias
# ---------------------------------------------------------------------------

def _twoway(s: "Session", rv: Variable, cv: Variable, mask, w, wtype, o: dict) -> None:
    out = s.output
    miss, nolab = bool(o.get("missing")), bool(o.get("nolabel"))
    rcats, rlabels, rok = _categories(s, rv, mask, missing=miss, nolabel=nolab)
    ccats, clabels, cok = _categories(s, cv, mask, missing=miss, nolabel=nolab)
    ok = rok & cok
    rc, cc = _codes(rv, rcats), _codes(cv, ccats)
    sel = ok & (rc >= 0) & (cc >= 0)
    ww = w if w is not None else np.ones(s.data.nobs)
    T = np.zeros((len(rcats), len(ccats)))
    np.add.at(T, (rc[sel], cc[sel]), ww[sel])
    N = float(T.sum())
    if N == 0:
        out.write("no observations\n", "text")
        s.r = {"N": 0.0}
        return
    rowt, colt = T.sum(1), T.sum(0)
    E = np.outer(rowt, colt) / N
    shows = []
    if not o.get("nofreq"):
        shows.append("frequency")
    if o.get("expected"):
        shows.append("expected frequency")
    if o.get("row"):
        shows.append("row percentage")
    if o.get("column"):
        shows.append("column percentage")
    if o.get("cell"):
        shows.append("cell percentage")
    if len(shows) > 1 and not o.get("nokey"):
        wkey = max(len(x) for x in shows)
        out.write("\n+" + "-" * (wkey + 2) + "+\n", "text")
        out.write("| " + "Key".ljust(wkey) + " |\n", "text")
        out.write("|" + "-" * (wkey + 2) + "|\n", "text")
        for x in shows:
            out.write("| " + x.center(wkey) + " |\n", "text")
        out.write("+" + "-" * (wkey + 2) + "+\n", "text")

    def values(i: int | None, j: int | None) -> list[str]:
        f = (T[i, j] if i is not None and j is not None else
             rowt[i] if i is not None else colt[j] if j is not None else N)
        rtot = rowt[i] if i is not None else N
        ctot = colt[j] if j is not None else N
        res = []
        for x in shows:
            if x == "frequency":
                res.append(_fmt_freq(f))
            elif x == "expected frequency":
                e = E[i, j] if i is not None and j is not None else f
                res.append(f"{e:.1f}")
            elif x == "row percentage":
                res.append(f"{100 * f / rtot:.2f}" if rtot else ".")
            elif x == "column percentage":
                res.append(f"{100 * f / ctot:.2f}" if ctot else ".")
            else:
                res.append(f"{100 * f / N:.2f}")
        return res

    W = min(30, max([10] + [len(x) for x in rlabels]))
    rtitle = rv.label or rv.name
    ctitle = cv.label or cv.name
    per_panel = max(1, (79 - (W + 2) - 11) // 11)   # VERIFICAR: quebra em painéis
    ncol = len(ccats)
    for start in range(0, ncol, per_panel):
        cols = list(range(start, min(ncol, start + per_panel)))
        last_panel = cols[-1] == ncol - 1
        block = 11 * len(cols)
        out.write("\n", "text")
        rhead = _header_lines(rtitle, W)
        above = rhead[:-1] if len(rhead) > 1 else [""]
        for k, left in enumerate(above):
            mid = ctitle.center(block).rstrip() if k == len(above) - 1 else ""
            out.write(f"{left:>{W}} |{mid}\n", "text")
        lab = "".join(f"{clabels[j][:10]:>10} " for j in cols)
        tot = "|     Total" if last_panel else ""
        out.write(f"{rhead[-1]:>{W}} |" + lab + tot + "\n", "text")
        sep = "-" * (W + 1) + "+" + "-" * block + ("+" + "-" * 10 if last_panel else "") + "\n"
        out.write(sep, "text")
        for i in range(len(rcats)):
            cells = [values(i, j) for j in cols]
            tot = values(i, None)
            for line in range(len(shows)):
                label = rlabels[i][:W] if line == 0 else ""
                out.write(f"{label:>{W}} |", "text")
                out.write("".join(f"{c[line]:>10} " for c in cells), "result")
                if last_panel:
                    out.write("|", "text")
                    out.write(f"{tot[line]:>10} ", "result")
                out.write("\n", "text")
            if len(shows) > 1:
                out.write(f"{'':>{W}} |" + " " * block + ("|" if last_panel else "") + "\n", "text")
        out.write(sep, "text")
        cells = [values(None, j) for j in cols]
        tot = values(None, None)
        for line in range(len(shows)):
            label = "Total" if line == 0 else ""
            out.write(f"{label:>{W}} |", "text")
            out.write("".join(f"{c[line]:>10} " for c in cells), "result")
            if last_panel:
                out.write("|", "text")
                out.write(f"{tot[line]:>10} ", "result")
            out.write("\n", "text")
    s.r = {"N": N, "r": float(len(rcats)), "c": float(len(ccats))}
    _tests(s, T, E, N, o)
    _store_mats(s, o, cell=T, rows=rcats, v=rv, cols=ccats if not cv.is_string else None)


def _tests(s: "Session", T: np.ndarray, E: np.ndarray, N: float, o: dict) -> None:
    from scipy import stats as st
    out = s.output
    r, c = T.shape
    df = (r - 1) * (c - 1)
    lines = []
    if o.get("chi2") or o.get("all"):
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2 = float(np.nansum((T - E) ** 2 / E))
        p = float(st.chi2.sf(chi2, df)) if df > 0 else M.SYSMISS
        lines.append((f"Pearson chi2({df})", f"{chi2:8.4f}", f"Pr = {p:5.3f}"))
        s.r.update({"chi2": chi2, "p": p})
    if o.get("lrchi2") or o.get("all"):
        with np.errstate(divide="ignore", invalid="ignore"):
            lr = float(2 * np.nansum(np.where(T > 0, T * np.log(T / E), 0.0)))
        p = float(st.chi2.sf(lr, df)) if df > 0 else M.SYSMISS
        lines.append((f"likelihood-ratio chi2({df})", f"{lr:8.4f}", f"Pr = {p:5.3f}"))
        s.r.update({"chi2_lr": lr, "p_lr": p})
    if o.get("v") or o.get("all"):
        with np.errstate(divide="ignore", invalid="ignore"):
            chi2 = float(np.nansum((T - E) ** 2 / E))
        V = math.sqrt(chi2 / (N * (min(r, c) - 1))) if min(r, c) > 1 else M.SYSMISS
        if r == 2 and c == 2:
            V = float((T[0, 0] * T[1, 1] - T[0, 1] * T[1, 0])
                      / math.sqrt(np.prod(T.sum(0)) * np.prod(T.sum(1))))
        lines.append(("Cramér's V", f"{V:8.4f}", ""))
        s.r["CramersV"] = V
    if o.get("gamma") or o.get("all"):
        g, ase = _gamma(T)
        lines.append(("gamma", f"{g:8.4f}", f"ASE = {ase:.3f}"))
        s.r.update({"gamma": g, "ase_gam": ase})
    if o.get("taub") or o.get("all"):
        tb, ase = _taub(T)
        lines.append(("Kendall's tau-b", f"{tb:8.4f}", f"ASE = {ase:.3f}"))
        s.r.update({"taub": tb, "ase_taub": ase})
    if lines:
        out.write("\n", "text")
        for name, value, extra in lines:
            out.write(f"{name:>26} = ", "text")
            out.write(value + ("   " + extra if extra else "") + "\n", "result")
    if o.get("exact"):
        if np.any(T != np.trunc(T)):
            raise StataError(198, "exact requires frequency counts")   # VERIFICAR
        if not lines:
            out.write("\n", "text")
        if (r, c) == (2, 2):
            _, p2 = st.fisher_exact(T.astype(int))
            # 1-sided: na direção observada (VERIFICAR)
            _, pl = st.fisher_exact(T.astype(int), alternative="less")
            _, pg = st.fisher_exact(T.astype(int), alternative="greater")
            p1 = min(pl, pg)
            out.write(f"{'Fisher' + chr(39) + 's exact':>26} = ", "text")
            out.write(f"{p2:>21.3f}\n", "result")
            out.write(f"{'1-sided Fisher' + chr(39) + 's exact':>26} = ", "text")
            out.write(f"{p1:>21.3f}\n", "result")
            s.r.update({"p_exact": float(p2), "p1_exact": float(p1)})
        else:
            p2 = fisher_rxc(T.astype(int))
            out.write(f"{'Fisher' + chr(39) + 's exact':>26} = ", "text")
            out.write(f"{p2:>21.3f}\n", "result")
            s.r["p_exact"] = float(p2)


def fisher_rxc(T: np.ndarray, limit: int = 5_000_000) -> float:
    """Teste exato de Fisher para tabelas r x c: soma as probabilidades
    (hipergeométrica multivariada) das tabelas com as mesmas margens e
    probabilidade menor ou igual à observada, enumerando-as uma a uma."""
    from math import lgamma
    rows = [int(x) for x in T.sum(1)]
    cols = [int(x) for x in T.sum(0)]
    n = sum(rows)
    const = sum(lgamma(x + 1) for x in rows) + sum(lgamma(x + 1) for x in cols) - lgamma(n + 1)

    def logp(table_cells) -> float:
        return const - sum(lgamma(x + 1) for x in table_cells)

    obs = logp(T.ravel().tolist())
    tol = 1e-7
    total = 0.0
    count = [0]
    r, c = len(rows), len(cols)

    def fill_row(i: int, colrem: list[int], acc: float) -> None:
        nonlocal total
        if i == r - 1:
            cells = colrem
            if sum(cells) != rows[i]:
                return
            lp = acc - sum(lgamma(x + 1) for x in cells)
            count[0] += 1
            if count[0] > limit:
                raise StataError(908, "table too large for the exact test; "
                                      "OpenDTA enumerates every table with the same margins")
            if lp <= obs + tol:
                total += math.exp(lp)
            return

        def place(j: int, left: int, cur: list[int], colrem: list[int]) -> None:
            if j == c - 1:
                if left <= colrem[j]:
                    row_cells = cur + [left]
                    rest = [colrem[k] - row_cells[k] for k in range(c)]
                    fill_row(i + 1, rest, acc - sum(lgamma(x + 1) for x in row_cells))
                return
            hi = min(left, colrem[j])
            for x in range(hi + 1):
                place(j + 1, left - x, cur + [x], colrem)
        place(0, rows[i], [], colrem)

    fill_row(0, cols, const)
    return min(1.0, total)


def _concordance(T: np.ndarray):
    r, c = T.shape
    P = Q = 0.0
    A = np.zeros_like(T)    # concordantes para cada célula
    D = np.zeros_like(T)
    for i in range(r):
        for j in range(c):
            A[i, j] = T[:i, :j].sum() + T[i + 1:, j + 1:].sum()
            D[i, j] = T[:i, j + 1:].sum() + T[i + 1:, :j].sum()
    P = float((T * A).sum())
    Q = float((T * D).sum())
    return P, Q, A, D


def _gamma(T: np.ndarray) -> tuple[float, float]:
    P, Q, A, D = _concordance(T)
    if P + Q == 0:
        return M.SYSMISS, M.SYSMISS
    g = (P - Q) / (P + Q)
    ase = 4 / (P + Q) ** 2 * math.sqrt(float((T * (Q * A - P * D) ** 2).sum()))
    return g, ase


def _taub(T: np.ndarray) -> tuple[float, float]:
    P, Q, A, D = _concordance(T)
    n = T.sum()
    wr = n ** 2 - float((T.sum(1) ** 2).sum())
    wc = n ** 2 - float((T.sum(0) ** 2).sum())
    if wr <= 0 or wc <= 0:
        return M.SYSMISS, M.SYSMISS
    w = math.sqrt(wr * wc)
    tb = (P - Q) / w
    # ASE (Agresti); VERIFICAR contra o Stata
    rows, cols = T.sum(1), T.sum(0)
    vij = np.zeros_like(T)
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            vij[i, j] = 2 * w * (A[i, j] - D[i, j]) + tb * (rows[i] * wc + cols[j] * wr)
    ase = math.sqrt(max(0.0, float((T * vij ** 2).sum()) - n ** 3 * tb ** 2 * (wr + wc) ** 2)) / w ** 2
    return tb, ase


# ---------------------------------------------------------------------------
# comandos
# ---------------------------------------------------------------------------

def _opts(text: str) -> dict:
    # a opção de Cramér é V maiúsculo
    return match_options(re.sub(r"\bV\b", "v", text), _OPTS) if text.strip() else {}


_OPTS = {"missing": 1, "nolabel": 5, "sort": 4, "generate": 3, "matcell": 7, "matrow": 7,
         "matcol": 7, "chi2": 2, "lrchi2": 3, "v": 1, "gamma": 1, "taub": 2, "exact": 2,
         "row": 1, "column": 3, "cell": 4, "expected": 5, "nofreq": 5, "nokey": 5, "all": 3,
         "plot": 4, "subpop": 3, "wrap": 4, "nolog": 5, "firstonly": 5}


@command("tabulate", "ta", byable=True)
def cmd_tabulate(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = _opts(p.options)
    names = s.expand_varlist(p.varlist)
    if not 1 <= len(names) <= 2:
        raise StataError(103 if len(names) > 2 else 100,
                         "too many variables specified" if len(names) > 2 else "varlist required")
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("fweight", "aweight", "iweight"))
    vs = [s.data.get(n) for n in names]
    for group, first in by_groups(s):
        if first is not None:
            _by_header(s, first)
        mask = base & group
        if len(vs) == 1:
            _oneway(s, vs[0], mask, w, wtype, o)
        else:
            _twoway(s, vs[0], vs[1], mask, w, wtype, o)


@command("tab1")
def cmd_tab1(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = _opts(p.options)
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("fweight", "aweight", "iweight"))
    for n in s.expand_varlist(p.varlist):
        s.output.write(f"\n-> tabulation of {n}  \n", "text")
        _oneway(s, s.data.get(n), base, w, wtype, o)


@command("tab2")
def cmd_tab2(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = _opts(p.options)
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("fweight", "aweight", "iweight"))
    names = s.expand_varlist(p.varlist)
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            s.output.write(f"\n-> tabulation of {names[i]} by {names[j]}  \n", "text")
            _twoway(s, s.data.get(names[i]), s.data.get(names[j]), base, w, wtype, o)


# ---------------------------------------------------------------------------
# table (sintaxe do Stata 14)
# ---------------------------------------------------------------------------

_STATS = {"freq", "mean", "sd", "sum", "rawsum", "count", "n", "max", "min", "median", "iqr",
          "p1", "p2", "p5", "p10", "p25", "p50", "p75", "p90", "p95", "p98", "p99"}


def _stat_value(stat: str, x: np.ndarray, w: np.ndarray | None) -> float:
    if stat in ("count", "n"):
        return float(len(x))
    if len(x) == 0:
        return M.SYSMISS
    m = S.moments(x, w, "fweight" if w is not None else "")
    if stat == "mean":
        return m.mean
    if stat == "sd":
        return m.sd
    if stat in ("sum", "rawsum"):
        return float(m.sum)
    if stat == "max":
        return m.max
    if stat == "min":
        return m.min
    xs = np.sort(x)
    if stat == "median":
        return S.percentile(xs, 50)
    if stat == "iqr":
        return S.percentile(xs, 75) - S.percentile(xs, 25)
    if stat.startswith("p"):
        return S.percentile(xs, float(stat[1:]))
    return M.SYSMISS


@command("table", byable=True)
def cmd_table(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"contents": 1, "format": 1, "row": 3, "col": 3, "missing": 1,
                                  "by": 2, "center": 3, "left": 4, "cellwidth": 5,
                                  "stubwidth": 5}) if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if not 1 <= len(names) <= 2:
        raise StataError(198, "table with more than two variables is not yet supported by OpenDTA")
    contents = str(o.get("contents", "freq")).split()
    specs: list[tuple[str, str | None]] = []
    k = 0
    while k < len(contents):
        st = contents[k]
        if st not in _STATS:
            raise StataError(198, f"{st} invalid statistic")   # VERIFICAR
        if st == "freq":
            specs.append(("freq", None))
            k += 1
        else:
            specs.append((st, s.expand_varlist(contents[k + 1])[0]))
            k += 2
    if len(specs) > 5:
        raise StataError(198, "too many statistics")
    fmt = str(o.get("format", "%9.0g")).strip()
    base = touse(s, p)
    w, wtype, base = weights(s, p, base, ("fweight", "aweight", "iweight", "pweight"))
    rv = s.data.get(names[0])
    cv = s.data.get(names[1]) if len(names) == 2 else None
    out = s.output
    rcats, rlabels, rok = _categories(s, rv, base, missing=bool(o.get("missing")), nolabel=False)
    rc = _codes(rv, rcats)
    if cv is not None:
        ccats, clabels, cok = _categories(s, cv, base, missing=bool(o.get("missing")), nolabel=False)
        cc = _codes(cv, ccats)
    else:
        ccats, clabels, cok, cc = [None], [], base, np.zeros(s.data.nobs, dtype=np.int64)
    ok = rok & cok

    def cell(sel: np.ndarray) -> list[str]:
        res = []
        for st, var in specs:
            if st == "freq":
                n = float((w[sel] if w is not None else np.ones(int(sel.sum()))).sum())
                res.append(_fmt_freq(n) if n else "")
            else:
                x = s.data.get(var).data
                ss = sel & S.valid(x)
                val = _stat_value(st, x[ss], w[ss] if w is not None else None)
                res.append("" if val >= M.SYSMISS else format_value(val, fmt, pad=False).strip())
        return res

    rows = list(range(len(rcats))) + (["Total"] if o.get("row") else [])
    cols = list(range(len(ccats))) + (["Total"] if o.get("col") and cv is not None else [])
    W = max([9] + [len(x) for x in rlabels] + [len(rv.name)])
    CW = 11
    title_lines = textwrap.wrap(rv.label or rv.name, W) or [rv.name]
    if cv is None:
        heads = ["Freq." if st == "freq" else f"{st}({var})" for st, var in specs]
        width = W + 2 + CW * len(heads)
        out.write("\n" + "-" * width + "\n", "text")
        for k, line in enumerate(title_lines):
            last = k == len(title_lines) - 1
            out.write(f"{line:<{W}} |" + ("".join(f"{h:>{CW}}" for h in heads) if last else "") + "\n", "text")
        out.write("-" * (W + 1) + "+" + "-" * (CW * len(heads)) + "\n", "text")
        for r in rows:
            sel = ok & (rc == r) if r != "Total" else ok
            label = rlabels[r] if r != "Total" else "Total"
            out.write(f"{label:>{W}} |", "text")
            out.write("".join(f"{v:>{CW}}" for v in cell(sel)) + "\n", "result")
        out.write("-" * width + "\n", "text")
        return
    ctitle = cv.label or cv.name
    block = CW * len(cols)
    width = W + 2 + block
    out.write("\n" + "-" * width + "\n", "text")
    out.write(f"{'':<{W}} |{ctitle.center(block).rstrip()}\n", "text")
    out.write(f"{title_lines[-1]:<{W}} |" + "".join(
        f"{(clabels[c] if c != 'Total' else 'Total'):>{CW}}" for c in cols) + "\n", "text")
    out.write("-" * (W + 1) + "+" + "-" * block + "\n", "text")
    for r in rows:
        rsel = ok & (rc == r) if r != "Total" else ok
        vals = [cell(rsel & (cc == c) if c != "Total" else rsel) for c in cols]
        for line in range(len(specs)):
            label = (rlabels[r] if r != "Total" else "Total") if line == 0 else ""
            out.write(f"{label:>{W}} |", "text")
            out.write("".join(f"{v[line]:>{CW}}" for v in vals) + "\n", "result")
    out.write("-" * width + "\n", "text")
