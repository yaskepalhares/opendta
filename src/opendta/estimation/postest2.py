"""More postestimation: testnl, predictnl, linktest, contrast, pwcompare, pwmean and suest
([R] testnl, [R] predictnl, [R] linktest, [R] contrast, [R] pwcompare, [R] pwmean, [R] suest).

Derivatives are numerical (as in margins.py). The margins behind contrast and pwcompare are
linear predictions averaged `asbalanced` (the default of contrast) or `asobserved`. Layouts follow
the manuals; the VERIFY marks are the details they do not show (compat/do/1903)."""

from __future__ import annotations

import itertools
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from ..lazy import lazy
from . import fvars
from .results import CoefRow, Estimates, coef_table, crit, current, g9, post, pvalue

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# expressions of the coefficients
# ---------------------------------------------------------------------------

def _coef_fn(s: "Session", est, expr: str):
    """fn(b) -> value of the expression with _b[name] (and _se[name]) taken from b."""
    from .postest import _coef_index
    base_se = np.sqrt(np.maximum(np.diag(est.V), 0))

    def fn(bb):
        def rep(m):
            kind, name = m.group(1), m.group(2)
            idx = _coef_index(est, s, name)
            if kind == "_se":
                return repr(float(base_se[idx]))
            return f"({float(bb[idx])!r})"
        t = re.sub(r"(_b|_coef|_se)\[([^\]]+)\]", rep, expr)
        return float(s.eval(t))
    return fn


def _jac(f, b: np.ndarray) -> np.ndarray:
    from .margins import _jacobian
    return _jacobian(f, b)


# ---------------------------------------------------------------------------
# testnl
# ---------------------------------------------------------------------------

def cmd_testnl(s: "Session", args: str) -> None:
    est = current(s)
    head, comma, opts = args.partition(",")
    if comma and opts.strip():
        match_options(opts, {"iterate": 4, "df": 2})
        raise StataError(198, "options of testnl are not available yet in OpenDTA")
    text = head.strip()
    if text.startswith("("):
        from .postest import _split_parens
        hyps = _split_parens(text)
    else:
        hyps = [text]
    fns, shown = [], []
    for h in hyps:
        lhs, eq_sign, rhs = h.partition("=") if "==" not in h else h.partition("==")
        if eq_sign:
            rhs = rhs.lstrip("=")
            e = f"({lhs}) - ({rhs})"
        else:
            e = f"({lhs})"
        fns.append(_coef_fn(s, est, e))
        shown.append(h.strip())
    b0 = np.asarray(est.b, dtype=float)
    g = np.array([f(b0) for f in fns])
    J = np.vstack([_jac(f, b0) for f in fns])
    GVG = J @ est.V @ J.T
    out = s.output
    out.write("\n", "text")
    for i, h in enumerate(shown, start=1):
        out.write(f"  ({i})  {h}\n", "text")
    out.write("\n", "text")
    q = len(fns)
    W = float(g @ np.linalg.solve(GVG, g))
    if est.stat == "t" and est.df_r:
        F = W / q
        dfr = float(est.df_r)
        p = float(st.f.sf(F, q, dfr))
        lab = f"F({q},{int(dfr)}) ="
        out.write(f"{lab:>24}", "text")
        out.write(f"{F:>12.2f}\n", "result")
        out.write(f"{'Prob > F =':>24}", "text")
        out.write(f"{p:>12.4f}\n", "result")
        s.r = {"p": p, "df_r": dfr, "df": float(q), "F": F}
    else:
        p = float(chi2_tail(W, q))
        out.write(f"{'chi2(' + str(q) + ') =':>24}", "text")
        out.write(f"{W:>12.2f}\n", "result")
        out.write(f"{'Prob > chi2 =':>24}", "text")
        out.write(f"{p:>12.4f}\n", "result")
        s.r = {"p": p, "df": float(q), "chi2": W}


# ---------------------------------------------------------------------------
# predictnl
# ---------------------------------------------------------------------------

_PNL = {"se": 2, "variance": 3, "ci": 2, "force": 5, "g": 1, "wald": 4, "p": 1, "level": 1, "iterate": 4}


def cmd_predictnl(s: "Session", args: str) -> None:
    est = current(s)
    from ..commands._util import eval_vector, touse
    from ..core.dataset import Variable
    from .margins import _transform, _xb
    from .postest import esample
    p = parse_standard(args)
    m = re.match(r"\s*(?:(byte|int|long|float|double)\s+)?(\w+)\s*=\s*(.+)$", p.varlist + (" = " + p.exp if p.exp else "")
                 if "=" not in p.varlist else p.varlist, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    vtype, newvar, expr = m.group(1) or "float", m.group(2), m.group(3).strip()
    o = match_options(p.options, _PNL) if p.options.strip() else {}
    ds = s.data
    if ds.has(newvar):
        raise StataError(110, f"variable {newvar} already defined")
    sample = esample(s) & touse(s, p)
    rows = np.flatnonzero(sample)
    if not len(rows):
        raise StataError(2000, "no observations")
    f_default, _ = _transform(est, "")
    names = {"__pnl_xb": "xb", "__pnl_pr": "pr"}

    def prep(text: str) -> str:
        text = re.sub(r"predict\(\s*xb\s*\)|\bxb\(\s*\)", "__pnl_xb", text)
        text = re.sub(r"predict\(\s*\)", "__pnl_pr", text)
        return text
    expr_t = prep(expr)

    def evaluate(bb: np.ndarray) -> np.ndarray:
        xb = _xb(s, est, bb, rows, {})
        for nm, vals in (("__pnl_xb", xb), ("__pnl_pr", f_default(xb, bb))):
            full = np.full(ds.nobs, SYS)
            full[rows] = vals
            if ds.has(nm):
                ds.get(nm).data = full
            else:
                ds.add(Variable(nm, "double", full))

        def rep(mt):
            return f"({float(bb[_cidx(s, est, mt.group(2))])!r})"
        t = re.sub(r"(_b|_coef)\[([^\]]+)\]", rep, expr_t)
        v = np.broadcast_to(np.asarray(eval_vector(s, t), dtype=float), (ds.nobs,))
        return v[rows]
    b0 = np.asarray(est.b, dtype=float)
    try:
        val = evaluate(b0)
        need_se = any(o.get(k) for k in ("se", "variance", "ci", "wald", "p"))
        J = None
        if need_se:
            J = np.zeros((len(rows), b0.size))
            for j in range(b0.size):
                h = 1e-4 * (abs(b0[j]) + 1e-3)

                def D(hh, j=j):
                    a, c = b0.copy(), b0.copy()
                    a[j] += hh
                    c[j] -= hh
                    return (evaluate(a) - evaluate(c)) / (2 * hh)
                J[:, j] = (4 * D(h / 2) - D(h)) / 3
    finally:
        for nm in names:
            if ds.has(nm):
                ds.drop_vars([nm])

    def put(name: str, vals: np.ndarray, vt: str = vtype) -> None:
        full = np.full(ds.nobs, SYS)
        full[rows] = vals
        ds.add(Variable(name, vt, full))
    put(newvar, val)
    if need_se:
        var = np.einsum("ij,jk,ik->i", J, np.asarray(est.V, dtype=float), J)
        lev = float(o["level"]) if isinstance(o.get("level"), str) else getattr(est, "level", 95.0)
        if isinstance(o.get("se"), str):
            put(o["se"], np.sqrt(np.maximum(var, 0)))
        if isinstance(o.get("variance"), str):
            put(o["variance"], var)
        if isinstance(o.get("ci"), str):
            lo, hi = o["ci"].split()
            z = crit("z", lev, None)
            se_ = np.sqrt(np.maximum(var, 0))
            put(lo, val - z * se_)
            put(hi, val + z * se_)
        if isinstance(o.get("wald"), str):
            put(o["wald"], val ** 2 / var)
        if isinstance(o.get("p"), str):
            put(o["p"], 2 * st.norm.sf(np.abs(val) / np.sqrt(var)))
    n_missing = 0
    if isinstance(val, np.ndarray):
        n_missing = int(np.sum(~np.isfinite(val)))
    if n_missing:
        s.output.write(f"({n_missing} missing value{'s' if n_missing != 1 else ''} generated)\n", "text")


def _cidx(s: "Session", est, name: str) -> int:
    from .postest import _coef_index
    return _coef_index(est, s, name)


# ---------------------------------------------------------------------------
# linktest
# ---------------------------------------------------------------------------

def cmd_linktest(s: "Session", args: str) -> None:
    est = current(s)
    from ..core.dataset import Variable
    from .postest import _restore, _snapshot, xb_all
    cmdline = s.e.get("cmdline")
    cmdname = s.e.get("cmd")
    dep = s.e.get("depvar") or est.depvar
    if not isinstance(cmdline, str) or not cmdname or not dep:
        raise StataError(321, "linktest requires a single-equation estimation command")    # VERIFY
    c = find_top(cmdline, ",")
    opts = "" if c == -1 else cmdline[c + 1:].strip()
    body = cmdline if c == -1 else cmdline[:c]
    word, _, rest = body.strip().partition(" ")
    p = parse_standard(rest)
    ds = s.data
    for nm in ("_hat", "_hatsq"):
        if ds.has(nm):
            raise StataError(110, f"variable {nm} already defined")      # VERIFY
    snap = _snapshot(s)
    xb, ok = xb_all(s, est)
    smp = np.asarray(ds.esample, dtype=bool) & ok
    hat = np.where(smp, xb, SYS)
    ds.add(Variable("_hat", "float", hat))
    ds.add(Variable("_hatsq", "float", np.where(smp, xb * xb, SYS)))
    parts = [word, dep, "_hat", "_hatsq"]
    if p.weight is not None:
        parts.append(f"[{p.weight[0]}={p.weight[1]}]")
    if p.if_:
        parts.append(f"if {p.if_}")
    if p.in_:
        parts.append(f"in {p.in_}")
    text = " ".join(parts) + (f", {opts}" if opts else "")
    try:
        s.interp._execute_expanded(text)
    finally:
        changed = ds.changed
        ds.drop_vars([n for n in ("_hat", "_hatsq") if ds.has(n)])
        ds.changed = changed
        _restore(s, snap)


# ---------------------------------------------------------------------------
# margins of factor terms for contrast and pwcompare
# ---------------------------------------------------------------------------

def _term_cells(s: "Session", est, tokens: list[str], asbalanced: bool, rows: np.ndarray):
    """For a factor term (a list of variables): the cell levels, the matrix A (cells × k) with
    margin = A b, and the number of levels of each factor."""
    from .margins import _jacobian, _xb
    ds = s.data
    factor_vars = {c.var for t in est.terms for c in t.comps if c.factor}
    for v in tokens:
        if v not in factor_vars:
            raise StataError(198, f"{v} is not a factor variable in the model")      # VERIFY
    levels = [sorted(set(ds.get(v).data[rows].tolist())) for v in tokens]
    base = {}
    if asbalanced:
        for t in est.terms:
            for comp in t.comps:
                if comp.var in tokens or comp.var in base:
                    continue
                vals = fvars.comp_values(ds, comp)[rows]
                if comp.factor:
                    lv = sorted(set(vals.tolist()))
                    base[comp.var] = {v: 1.0 / len(lv) for v in lv}
                else:
                    base[comp.var] = np.full(len(rows), float(np.mean(vals)))
    b = np.asarray(est.b, dtype=float)
    cells = list(itertools.product(*levels))

    def margin_all(bb):
        out = []
        for cell in cells:
            ov = {k: (dict(v) if isinstance(v, dict) else v) for k, v in base.items()}
            for v, lv in zip(tokens, cell):
                ov[v] = np.full(len(rows), lv)
            out.append(float(np.mean(_xb(s, est, bb, rows, ov))))
        return np.array(out)
    A = _jacobian(margin_all, b)
    return cells, A, [len(lv) for lv in levels]


def _label_level(s: "Session", var: str, lv: float) -> str:
    from .. import core
    return f"{int(lv)}" if float(lv).is_integer() else f"{lv:g}"


_OPS = ("ar", "r", "a", "g", "h", "j")


def _op_matrix(op: str, n: int):
    """Contrast matrix and labels for n levels (labels use the level positions 1..n)."""
    rows, labs = [], []
    if op in ("", "r"):
        for j in range(1, n):
            v = np.zeros(n)
            v[j], v[0] = 1, -1
            rows.append(v)
            labs.append((j, 0, "vs"))
    elif op == "a":
        for j in range(1, n):
            v = np.zeros(n)
            v[j], v[j - 1] = 1, -1
            rows.append(v)
            labs.append((j, j - 1, "vs"))
    elif op == "ar":
        for j in range(n - 1):
            v = np.zeros(n)
            v[j], v[j + 1] = 1, -1
            rows.append(v)
            labs.append((j, j + 1, "vs"))
    elif op == "g":
        for j in range(n):
            v = np.full(n, -1.0 / n)
            v[j] += 1
            rows.append(v)
            labs.append((j, None, "mean"))
        rows, labs = rows[:-1], labs[:-1]
    elif op == "h":
        for j in range(n - 1):
            v = np.zeros(n)
            v[j] = 1
            v[j + 1:] = -1.0 / (n - 1 - j)
            rows.append(v)
            labs.append((j, None, "next"))
    elif op == "j":
        for j in range(1, n):
            v = np.zeros(n)
            v[j] = 1
            v[:j] = -1.0 / j
            rows.append(v)
            labs.append((j, None, "prev"))
    else:
        raise StataError(198, f"operator {op}. is not available yet in OpenDTA")
    return np.array(rows), labs


def _split_terms(text: str):
    """`r.g a.h g#h` -> [(operator, [vars])]."""
    out = []
    for tok in text.split():
        m = re.match(r"^(?:(" + "|".join(_OPS) + r")\.)?(\w+(?:#\w+)*)$", tok)
        if not m:
            raise StataError(198, f"{tok} is not allowed here")
        op = m.group(1) or ""
        vars_ = [v.split(".")[-1] for v in m.group(2).split("#")]
        out.append((op, vars_, tok))
    return out


def _wald_L(L: np.ndarray, m: np.ndarray, Vm: np.ndarray):
    d = L @ m
    LVL = L @ Vm @ L.T
    rank = int(np.linalg.matrix_rank(LVL))
    return float(d @ np.linalg.pinv(LVL) @ d), rank


def cmd_contrast(s: "Session", args: str) -> None:
    est = current(s)
    p = parse_standard(args)
    o = match_options(p.options, {"asbalanced": 3, "asobserved": 3, "effects": 2, "level": 1, "nowald": 3,
                                  "noeffects": 3, "pveffects": 2, "cieffects": 2, "mcompare": 2}) \
        if p.options.strip() else {}
    if o.get("mcompare"):
        raise StataError(198, "option mcompare() is not available yet in OpenDTA")
    from .postest import esample
    rows = np.flatnonzero(esample(s))
    asbal = not o.get("asobserved")
    terms = _split_terms(p.varlist)
    if not terms:
        raise StataError(198, "no terms")
    out = s.output
    df_r = est.df_r if est.stat == "t" and est.df_r else None
    stat_name = "F" if df_r else "chi2"
    lev = float(o["level"]) if isinstance(o.get("level"), str) else getattr(est, "level", 95.0)
    results = []
    eff_rows = []
    for op, vars_, tok in terms:
        cells, A, nl = _term_cells(s, est, vars_, asbal, rows)
        Vm = A @ np.asarray(est.V, dtype=float) @ A.T
        m = A @ np.asarray(est.b, dtype=float)
        mats = [_op_matrix(op, n) for n in nl]
        if len(vars_) == 1:
            L, labs = mats[0]
            labels = []
            levs = sorted(set(c[0] for c in cells))
            for (a_, b_, kind) in labs:
                la = _label_level(s, vars_[0], levs[a_])
                if kind == "vs":
                    labels.append(f"({la} vs {_label_level(s, vars_[0], levs[b_])})")
                elif kind == "mean":
                    labels.append(f"({la} vs mean)")
                elif kind == "next":
                    labels.append(f"({la} vs next)")      # VERIFY
                else:
                    labels.append(f"({la} vs previous)")  # VERIFY
        else:
            L = np.array([np.kron(*[mm[0][i] for mm in mats[:1]], mats[1][0][j])
                          for i in range(len(mats[0][0])) for j in range(len(mats[1][0]))]) \
                if len(vars_) == 2 else None
            if L is None:
                raise StataError(198, "interactions of more than two factors are not available yet")
            labels = []
        joint, rank = _wald_L(L, m, Vm)
        results.append((tok, labels, L, m, Vm, joint, rank))
    out.write("\nContrasts of marginal linear predictions\n\n", "text")
    out.write(f"Margins      : {'asbalanced' if asbal else 'asobserved'}\n\n", "text")
    hdr = f"{'df':>11}{stat_name:>12}{'P>' + stat_name:>11}"
    width = 13 + 1 + 34
    out.write("-" * width + "\n", "text")
    out.write(f"{'':>13}|{hdr}\n", "text")
    out.write("-" * 13 + "+" + "-" * 34 + "\n", "text")

    def line(label: str, df: int, w: float) -> None:
        if df_r:
            f_ = w / df
            pv = float(st.f.sf(f_, df, df_r))
            val = f_
        else:
            pv = float(chi2_tail(w, df))
            val = w
        out.write(f"{label:>12} |", "text")
        out.write(f"{df:>11}{val:>12.2f}{pv:>11.4f}\n", "result")
    for tok, labels, L, m, Vm, joint, rank in results:
        op_given = "." in tok
        if not op_given or not labels:
            line(tok.split(".")[-1], rank, joint)
            continue
        out.write(f"{tok.split('.')[-1]:>12} |\n", "text")
        for lab, row in zip(labels, L):
            w, rk = _wald_L(row.reshape(1, -1), m, Vm)
            line(lab, rk, w)
            eff_rows.append((lab, float(row @ m), float(np.sqrt(max(row @ Vm @ row, 0)))))
        line("Joint", rank, joint)
    out.write(f"{'':>13}|\n", "text")
    if df_r:
        out.write(f"{'Denominator':>12} |{int(df_r):>11}\n", "text")
    out.write("-" * width + "\n", "text")
    if o.get("effects") and eff_rows:
        out.write("\n", "text")
        rws = [CoefRow("coef", lab, b_, se_) for lab, b_, se_ in eff_rows]
        coef_table(s, "", rws, stat="t" if df_r else "z", df=df_r, level=lev, coef_title="Contrast")
    s.r = {}


# ---------------------------------------------------------------------------
# pairwise comparisons
# ---------------------------------------------------------------------------

def _adjust(method: str, tval: float, df, m: int, k: int, level: float):
    """(adjusted p, critical value) for a comparison under the multiple-comparison method."""
    a = 1 - level / 100
    p = 2 * (st.t.sf(abs(tval), df) if df else st.norm.sf(abs(tval)))
    if method in ("", "noadjust"):
        return p, (st.t.ppf(1 - a / 2, df) if df else st.norm.ppf(1 - a / 2))
    if method == "bonferroni":
        return min(1.0, p * m), (st.t.ppf(1 - a / (2 * m), df) if df else st.norm.ppf(1 - a / (2 * m)))
    if method == "sidak":
        a2 = 1 - (1 - a) ** (1.0 / m)
        return 1 - (1 - p) ** m, (st.t.ppf(1 - a2 / 2, df) if df else st.norm.ppf(1 - a2 / 2))
    if method == "scheffe":
        f_ = tval ** 2 / (k - 1)
        crit_ = np.sqrt((k - 1) * st.f.ppf(1 - a, k - 1, df))
        return float(st.f.sf(f_, k - 1, df)), float(crit_)
    if method == "tukey":
        q = abs(tval) * np.sqrt(2)
        return float(st.studentized_range.sf(q, k, df)), float(st.studentized_range.ppf(1 - a, k, df) / np.sqrt(2))
    raise StataError(198, f"mcompare({method}) is not available yet in OpenDTA")


def _pairwise_table(s: "Session", term: str, labels: list[str], means: np.ndarray, Vm: np.ndarray,
                    df, o: dict, lev: float, title_lines: list[str]) -> None:
    out = s.output
    method = str(o.get("mcompare") or "noadjust").lower()
    k = len(labels)
    pairs = [(j, i) for i in range(k) for j in range(i + 1, k)]
    m = len(pairs)
    eff = o.get("effects")
    pv_eff = o.get("pveffects")
    adj = method not in ("noadjust", "")
    rows, tvals = [], []
    for (j, i) in pairs:
        d = means[j] - means[i]
        se = float(np.sqrt(max(Vm[j, j] + Vm[i, i] - 2 * Vm[i, j], 0)))
        tv = d / se if se > 0 else np.nan
        p, c = _adjust(method, tv, df, m, k, lev) if se > 0 else (np.nan, np.nan)
        rows.append((f"{labels[j]} vs {labels[i]}", d, se, tv, p, c))
    for ln in title_lines:
        out.write(ln, "text")
    cw = 13 + 1 + 64
    lv = f"{lev:g}"
    tag = {"noadjust": "Unadjusted", "bonferroni": "Bonferroni", "sidak": "Sidak", "scheffe": "Scheffe",
           "tukey": "Tukey"}.get(method, method.title())
    if adj:
        out.write("-" * cw + "\n", "text")
        out.write(f"{'':>13}|   Number of\n{'':>13}| Comparisons\n", "text")
        out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
        out.write(f"{term:>12} |{m:>12}\n", "result")
        out.write("-" * cw + "\n", "text")
        out.write("\n", "text")
    out.write("-" * cw + "\n", "text")
    if eff or pv_eff:
        cols = ("   Contrast   Std. Err.      " + ("t" if df else "z") + "    P>|" + ("t" if df else "z") + "|"
                if not adj else "")
        out.write(f"{'':>13}|{' ' * 38}{tag}\n" if not adj else f"{'':>13}|{' ' * 40}{tag} \n", "text")
        out.write(f"{'':>13}|   Contrast   Std. Err.      {'t' if df else 'z'}    P>|{'t' if df else 'z'}|     "
                  f"[{lv}% Conf. Interval]\n", "text")
    else:
        out.write(f"{'':>13}|{' ' * 28}{tag}\n", "text")
        out.write(f"{'':>13}|   Contrast   Std. Err.     [{lv}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    out.write(f"{term:>12} |\n", "text")
    for lab, d, se, tv, p, c in rows:
        out.write(f"{lab:>12}  |" if False else f"{lab:>12} |", "text")
        lo, hi = d - c * se, d + c * se
        if eff or pv_eff:
            out.write(f" {g9(d):>10} {g9(se):>10}  {tv:>7.2f}  {p:>6.3f}    {g9(lo):>9}   {g9(hi):>9}\n", "result")
        else:
            out.write(f" {g9(d):>10} {g9(se):>10}     {g9(lo):>9}   {g9(hi):>9}\n", "result")
    out.write("-" * cw + "\n", "text")


_PWOPT = {"mcompare": 2, "effects": 2, "cieffects": 2, "pveffects": 2, "level": 1, "asbalanced": 3,
          "asobserved": 3, "groups": 2, "sort": 2, "unequal": 3, "over": 2}


def cmd_pwcompare(s: "Session", args: str) -> None:
    est = current(s)
    p = parse_standard(args)
    o = match_options(p.options, _PWOPT) if p.options.strip() else {}
    if o.get("groups") or o.get("sort"):
        raise StataError(198, "options groups and sort are not available yet in OpenDTA")
    from .postest import esample
    rows = np.flatnonzero(esample(s))
    asbal = not o.get("asobserved")
    lev = float(o["level"]) if isinstance(o.get("level"), str) else getattr(est, "level", 95.0)
    df = est.df_r if est.stat == "t" and est.df_r else None
    for _, vars_, tok in _split_terms(p.varlist):
        cells, A, nl = _term_cells(s, est, vars_, asbal, rows)
        Vm = A @ np.asarray(est.V, dtype=float) @ A.T
        m = A @ np.asarray(est.b, dtype=float)
        labels = ["#".join(_label_level(s, v, lv) for v, lv in zip(vars_, c)) for c in cells]
        head = ["\nPairwise comparisons of marginal linear predictions\n\n",
                f"Margins      : {'asbalanced' if asbal else 'asobserved'}\n\n"]
        _pairwise_table(s, "#".join(vars_), labels, m, Vm, df, o, lev, head)


def cells_of(s: "Session", overs: list[str], mask: np.ndarray):
    """Cells of the over() variables (sorted) with their labels, and the cell number of each row."""
    from ..commands.inspect import cell_text
    ds = s.data
    cols = []
    for v in overs:
        var = ds.get(v)
        cols.append(np.asarray(var.raw, dtype=object) if var.is_string else var.data)
    keys = {}
    for i in np.flatnonzero(mask):
        keys.setdefault(tuple(c[i] for c in cols), i)
    order = sorted(keys)
    labels = ["#".join(cell_text(ds, ds.get(v), keys[k]).strip() for v, _ in zip(overs, k)) for k in order]
    ids = np.full(ds.nobs, -1)
    pos = {k: j for j, k in enumerate(order)}
    for i in np.flatnonzero(mask):
        ids[i] = pos[tuple(c[i] for c in cols)]
    return labels, ids


def cmd_pwmean(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _PWOPT) if p.options.strip() else {}
    if not isinstance(o.get("over"), str):
        raise StataError(198, "option over() required")
    from ..commands.summarize import numeric_vars
    v = numeric_vars(s, p.varlist, allow_string=False)
    if len(v) != 1:
        raise StataError(103, "too many variables specified")
    from ..commands._util import touse
    mask = touse(s, p) & (v[0].data < SYS)
    overs = s.expand_varlist(o["over"])
    labels, ids = cells_of(s, overs, mask)
    k = len(labels)
    y = v[0].data
    n_i = np.array([np.sum(mask & (ids == g)) for g in range(k)], dtype=float)
    m_i = np.array([y[mask & (ids == g)].mean() for g in range(k)])
    ss = sum(((y[mask & (ids == g)] - m_i[g]) ** 2).sum() for g in range(k))
    N = int(n_i.sum())
    lev = float(o["level"]) if isinstance(o.get("level"), str) else float(s.settings.get("level", 95))
    if o.get("unequal"):
        s2 = np.array([y[mask & (ids == g)].var(ddof=1) for g in range(k)])
        Vm = np.diag(s2 / n_i)
        df = None       # Satterthwaite df per pair is not implemented
        raise StataError(198, "option unequal is not available yet in OpenDTA")
    df = N - k
    Vm = np.diag(np.full(k, ss / df) / n_i)
    head = ["\nPairwise comparisons of means with equal variances\n\n",
            f"{'Over':<13}: {' '.join(overs)}\n"]
    _pairwise_table(s, "#".join(overs), labels, m_i, Vm, df, o, lev, head)
    s.r = {}


# ---------------------------------------------------------------------------
# suest
# ---------------------------------------------------------------------------

def _scores_hessian(s: "Session", est):
    """Per-observation scores (n × p) over all rows, Hessian (p × p), parameter names and the
    equation names, for regress, logit, probit and poisson."""
    from .postest import design_all, xb_all
    cmd = getattr(est, "ml_cmd", est.cmd)
    if cmd == "logistic":
        cmd = "logit"
    ds = s.data
    X, ok = design_all(s, est)
    sample = np.asarray(ds.esample, dtype=bool) & ok if ds.esample is not None else ok
    y = ds.get(est.depvar).data
    n, k = X.shape
    b = np.asarray(est.b, dtype=float)
    xb = X @ b[:k] if len(b) >= k else X @ b
    if getattr(est, "offset_name", None):
        nm = est.offset_name
        xb = xb + (np.log(ds.get(nm[3:-1]).data) if nm.startswith("ln(") else ds.get(nm).data)
    smp = sample.copy()
    if cmd == "regress":
        e = np.where(smp, y - xb, 0.0)
        nn = smp.sum()
        sig2 = float((e[smp] ** 2).sum() / nn)
        S = np.column_stack([X * (e / sig2)[:, None], np.where(smp, 0.5 * (e ** 2 / sig2 - 1), 0.0)])
        H = np.zeros((k + 1, k + 1))
        Xs = X[smp]
        H[:k, :k] = -Xs.T @ Xs / sig2
        H[:k, k] = H[k, :k] = -(Xs.T @ e[smp]) / sig2
        H[k, k] = -0.5 * (e[smp] ** 2).sum() / sig2
        bnew = np.concatenate([b[:k], [np.log(sig2)]])
        names = list(est.names[:k]) + ["_cons"]
        eqs = ["mean"] * k + ["lnvar"]
        return S * smp[:, None], H, bnew, names, eqs, smp
    if cmd == "logit":
        pr = 1 / (1 + np.exp(-xb))
        S = X * np.where(smp, y - pr, 0.0)[:, None]
        w = np.where(smp, pr * (1 - pr), 0.0)
    elif cmd == "probit":
        eta = xb
        Phi, phi = st.norm.cdf(eta), st.norm.pdf(eta)
        lam1, lam0 = phi / Phi, phi / (1 - Phi)
        sc = np.where(y == 1, lam1, -lam0)
        S = X * np.where(smp, sc, 0.0)[:, None]
        w = np.where(smp, np.where(y == 1, lam1 * (eta + lam1), lam0 * (lam0 - eta)), 0.0)
    elif cmd == "poisson":
        mu = np.exp(xb)
        S = X * np.where(smp, y - mu, 0.0)[:, None]
        w = np.where(smp, mu, 0.0)
    else:
        raise StataError(322, f"suest is not available yet in OpenDTA after {cmd}")
    H = -(X * w[:, None]).T @ X
    eqs = [est.depvar] * k
    return S, H, b[:k], list(est.names[:k]), eqs, smp


def cmd_suest(s: "Session", args: str) -> None:
    from .postest import _restore, _snapshot, _stored
    c = find_top(args, ",")
    names = (args if c == -1 else args[:c]).split()
    o = match_options("" if c == -1 else args[c + 1:], {"vce": 3, "svy": 3, "cluster": 2}) if c != -1 else {}
    if o.get("svy"):
        raise StataError(198, "option svy is not available yet in OpenDTA")
    if len(names) < 2:
        raise StataError(198, "at least two models are required")      # VERIFY
    store = _stored(s)
    snap0 = _snapshot(s)
    parts = []
    try:
        for nm in names:
            if nm not in store:
                raise StataError(111, f"estimation result {nm} not found")
            _restore(s, store[nm])
            parts.append((nm,) + _scores_hessian(s, s._estimates))
    finally:
        _restore(s, snap0)
    sizes = [len(p[3]) for p in parts]
    total = sum(sizes)
    n = s.data.nobs
    S = np.hstack([p[1] for p in parts])
    H = np.zeros((total, total))
    pos = 0
    for p, k in zip(parts, sizes):
        H[pos:pos + k, pos:pos + k] = p[2]
        pos += k
    b = np.concatenate([p[3] for p in parts])
    cnames = [x for p in parts for x in p[4]]
    ceqs = [f"{p[0]}_{q}" for p in parts for q in p[5]]
    union = np.zeros(n, dtype=bool)
    for p in parts:
        union |= p[6]
    clustvar = ""
    if isinstance(o.get("cluster"), str):
        clustvar = o["cluster"]
    elif isinstance(o.get("vce"), str) and o["vce"].split()[0].lower().startswith("cl"):
        clustvar = o["vce"].split()[1]
    Su = S[union]
    if clustvar:
        cv = s.data.get(clustvar)
        key = np.asarray(cv.raw if cv.is_string else cv.data)[union]
        ids = np.unique(key, return_inverse=True)[1]
        Sg = np.zeros((ids.max() + 1, total))
        np.add.at(Sg, ids, Su)
        G = Sg.shape[0]
        adj = G / (G - 1)
        meat = adj * Sg.T @ Sg
    else:
        N = int(union.sum())
        meat = N / (N - 1) * Su.T @ Su
    Hinv = np.linalg.inv(-H)
    V = Hinv @ meat @ Hinv
    Nobs = int(union.sum())
    est = Estimates("suest", "", cnames, b, V, Nobs, "z", None)
    est.eqnames = ceqs
    est.vcetype = "Robust"
    est.title = "Simultaneous results for " + ", ".join(names)
    est.level = float(s.settings.get("level", 95))
    rows = []
    last = None
    for j, (nm_, eq_) in enumerate(zip(cnames, ceqs)):
        if eq_ != last:
            if last is not None:
                rows.append(CoefRow("sep"))
            rows.append(CoefRow("eq", eq_))
            last = eq_
        rows.append(CoefRow("coef", nm_, float(b[j]), float(np.sqrt(max(V[j, j], 0))), name=nm_))
    est.rows = rows
    scalars = [("N", float(Nobs))]
    if clustvar:
        scalars.append(("N_clust", float(G)))
    macros = [("cmd", "suest"), ("names", " ".join(names)), ("vce", "robust"), ("vcetype", "Robust"),
              ("title", est.title), ("properties", "b V")]
    if clustvar:
        macros.append(("clustvar", clustvar))
    post(s, est, list(reversed(scalars)), list(reversed(macros)), union)
    out = s.output
    out.write(f"\n{est.title}\n\n", "text")
    out.write(f"{'':<48}{'Number of obs':<18}= ", "text")
    out.write(f"{Nobs:>10,}\n", "result")
    out.write("\n", "text")
    if clustvar:
        out.write(f"{'':>48}(Std. Err. adjusted for {G} clusters in {clustvar})\n", "text")
    coef_table(s, "", rows, stat="z", df=None, level=est.level, vcetype="Robust")
