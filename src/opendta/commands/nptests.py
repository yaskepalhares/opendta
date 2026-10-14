"""Testes não paramétricos e de normalidade: ranksum, signrank, kwallis,
spearman, swilk, sktest e sdtest ([R] ranksum, signrank, kwallis,
spearman, swilk, sktest, sdtest).

Fórmulas dos manuais (seção Methods and formulas). Os layouts vêm de
saídas do Stata 14 de memória e ainda não foram conferidos com logs: as
linhas têm a marca VERIFICAR e o caso compat/do/0801_testes.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core import stats as S
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .registry import command
from .summarize import _corr_table, _level, _ttest_table, g9, name12, numeric_vars

if TYPE_CHECKING:
    from ..session import Session


def _groups(s: "Session", gname: str, mask: np.ndarray):
    """Valores distintos de by() na amostra, com o texto (rótulo) de cada um."""
    from .inspect import cell_text
    g = s.data.get(s.expand_varlist(gname)[0])
    if g.is_string:
        col = np.asarray(g.raw, dtype=object)
        ok = np.array([x != "" for x in col])
    else:
        col = g.data
        ok = col < M.SYSMISS
    mask = mask & ok
    vals = sorted(set(col[mask].tolist()))
    out = []
    for v in vals:
        sel = mask & (col == v)
        out.append((v, cell_text(s.data, g, int(np.flatnonzero(sel)[0])), sel))
    return g, out


def _tie_sum(x: np.ndarray) -> float:
    _, counts = np.unique(x, return_counts=True)
    return float(((counts ** 3) - counts).sum())


def _rank_table(s: "Session", head: str, cols: tuple[str, str, str], rows, total) -> None:
    """Tabela de ranksum/signrank: grupo | obs, soma dos postos, esperado."""
    out = s.output
    # VERIFICAR: larguras (obs %9.0g, somas %12.0g)
    out.write(f"{name12(head):>12} |{cols[0]:>9}{cols[1]:>12}{cols[2]:>12}\n", "text")
    out.write("-" * 13 + "+" + "-" * 33 + "\n", "text")
    for name, n, sobs, sexp in rows:
        out.write(f"{name12(str(name)):>12} |", "text")
        out.write(f"{n:>9.0f}{g9(sobs, '%12.0g'):>12}{g9(sexp, '%12.0g'):>12}\n", "result")
    out.write("-" * 13 + "+" + "-" * 33 + "\n", "text")
    name, n, sobs, sexp = total
    out.write(f"{name:>12} |", "text")
    out.write(f"{n:>9.0f}{g9(sobs, '%12.0g'):>12}{g9(sexp, '%12.0g'):>12}\n", "result")


def _variance_lines(s: "Session", lines: list[tuple[str, float]], total: float) -> None:
    out = s.output
    out.write("\n", "text")
    for label, v in lines:
        if abs(v) < 0.005:
            v = 0.0     # evita "-0.00" (compat 0801)
        out.write(f"{label:<21}", "text")
        out.write(f"{v:>10.2f}\n", "result")
    out.write(" " * 21 + "-" * 10 + "\n", "text")
    out.write(f"{'adjusted variance':<21}", "text")
    out.write(f"{total:>10.2f}\n", "result")


def _z_lines(s: "Session", ho: str, z: float, p: float) -> None:
    out = s.output
    out.write(f"\n{ho}\n", "text")
    out.write(f"{'z':>14} = ", "text")
    out.write(f"{z:>7.3f}\n", "result")
    out.write(f"{'Prob > |z|':>14} = ", "text")
    out.write(f"{p:>8.4f}\n", "result")


# ---------------------------------------------------------------------------
# ranksum
# ---------------------------------------------------------------------------

@command("ranksum", "ranksum")
def cmd_ranksum(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "porder": 2, "exact": 2})
    if not isinstance(o.get("by"), str):
        raise StataError(198, "option by() required")
    v = numeric_vars(s, p.varlist, allow_string=False)
    if len(v) != 1:
        raise StataError(103, "too many variables specified")
    v = v[0]
    mask = touse(s, p) & S.valid(v.data)
    g, groups = _groups(s, o["by"], mask)
    if len(groups) != 2:
        raise StataError(420, f"{len(groups)} groups found, 2 required")   # VERIFICAR
    sel = groups[0][2] | groups[1][2]
    x = v.data[sel]
    r = st.rankdata(x)
    in1 = groups[0][2][sel]
    n1, n2 = int(in1.sum()), int((~in1).sum())
    N = n1 + n2
    r1, r2 = float(r[in1].sum()), float(r[~in1].sum())
    e1, e2 = n1 * (N + 1) / 2, n2 * (N + 1) / 2
    var_u = n1 * n2 * (N + 1) / 12
    adj = -n1 * n2 * _tie_sum(x) / (12 * N * (N - 1)) if N > 1 else 0.0
    var_a = var_u + adj
    z = (r1 - e1) / math.sqrt(var_a) if var_a > 0 else M.SYSMISS
    pv = 2 * float(st.norm.sf(abs(z))) if z < M.SYSMISS else M.SYSMISS
    out = s.output
    out.write("\nTwo-sample Wilcoxon rank-sum (Mann-Whitney) test\n\n", "text")
    _rank_table(s, g.name, ("obs", "rank sum", "expected"),
                [(groups[0][1], n1, r1, e1), (groups[1][1], n2, r2, e2)],
                ("combined", N, r1 + r2, e1 + e2))
    _variance_lines(s, [("unadjusted variance", var_u), ("adjustment for ties", adj)], var_a)
    ho = f"Ho: {v.name}({g.name}=={groups[0][1]}) = {v.name}({g.name}=={groups[1][1]})"
    _z_lines(s, ho, z, pv)
    porder = (r1 - n1 * (n1 + 1) / 2) / (n1 * n2)
    if o.get("porder"):
        out.write(f"\nP{{{v.name}({g.name}=={groups[0][1]}) > {v.name}({g.name}=={groups[1][1]})}} = ", "text")
        out.write(f"{porder:.3f}\n", "result")
    g1 = float(groups[0][0]) if not isinstance(groups[0][0], str) else groups[0][0]
    # ordem do return list do Stata 14 (compat 0801); porder só com a opção
    _set_r(s, ([("porder", porder)] if o.get("porder") else [])
           + [("N_2", float(n2)), ("N_1", float(n1)), ("Var_a", var_a), ("z", z), ("sum_exp", e1),
              ("sum_obs", r1), ("group1", g1)])


# ---------------------------------------------------------------------------
# signrank
# ---------------------------------------------------------------------------

@command("signrank", "signrank")
def cmd_signrank(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.exp is None and "==" in p.varlist:
        a, _, b = p.varlist.partition("==")
        p.varlist, p.exp = a.strip(), b.strip()
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    rhs = p.exp.strip()
    mask = touse(s, p) & S.valid(v.data)
    try:
        other = s.data.get(s.expand_varlist(rhs)[0]) if not _is_num(rhs) else None
    except StataError:
        other = None
    if other is not None:
        mask &= S.valid(other.data)
        d = v.data[mask] - other.data[mask]
        rhs_name = other.name
    else:
        val = float(s.eval(rhs))
        d = v.data[mask] - val
        rhs_name = rhs
    n = len(d)
    if n == 0:
        raise StataError(2000, "no observations")
    ad = np.abs(d)
    r = st.rankdata(ad)
    pos, neg, zer = d > 0, d < 0, d == 0
    sp, sn, sz = float(r[pos].sum()), float(r[neg].sum()), float(r[zer].sum())
    tot = n * (n + 1) / 2
    e = (tot - sz) / 2
    n0 = int(zer.sum())
    var_u = n * (n + 1) * (2 * n + 1) / 24
    ties = -_tie_sum(ad[~zer]) / 48 if (~zer).any() else 0.0
    zeros = -n0 * (n0 + 1) * (2 * n0 + 1) / 24
    var_a = var_u + ties + zeros
    z = (sp - e) / math.sqrt(var_a) if var_a > 0 else M.SYSMISS
    pv = 2 * float(st.norm.sf(abs(z))) if z < M.SYSMISS else M.SYSMISS
    out = s.output
    out.write("\nWilcoxon signed-rank test\n\n", "text")
    _rank_table(s, "sign", ("obs", "sum ranks", "expected"),
                [("positive", int(pos.sum()), sp, e), ("negative", int(neg.sum()), sn, e),
                 ("zero", n0, sz, sz)], ("all", n, tot, tot))
    # VERIFICAR: ordem das linhas de ajuste (empates, zeros)
    _variance_lines(s, [("unadjusted variance", var_u), ("adjustment for ties", ties),
                        ("adjustment for zeros", zeros)], var_a)
    _z_lines(s, f"Ho: {v.name} = {rhs_name}", z, pv)
    _set_r(s, [("N_tie", float(n0)), ("N_neg", float(neg.sum())), ("N_pos", float(pos.sum())),
               ("Var_a", var_a), ("z", z), ("sum_neg", sn), ("sum_pos", sp)])


def _is_num(t: str) -> bool:
    try:
        float(t)
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------------------
# kwallis
# ---------------------------------------------------------------------------

@command("kwallis", "kwallis")
def cmd_kwallis(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2})
    if not isinstance(o.get("by"), str):
        raise StataError(198, "option by() required")
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    g, groups = _groups(s, o["by"], mask)
    sel = np.zeros_like(mask)
    for _, _, m in groups:
        sel |= m
    x = v.data[sel]
    N = len(x)
    r = st.rankdata(x)
    rows = []
    H = 0.0
    for val, lab, m in groups:
        rk = r[m[sel]]
        rows.append((lab, int(m.sum()), float(rk.sum())))
        H += rk.sum() ** 2 / len(rk)
    H = 12 / (N * (N + 1)) * H - 3 * (N + 1)
    df = len(groups) - 1
    corr = 1 - _tie_sum(x) / (N ** 3 - N)
    Ht = H / corr if corr > 0 else M.SYSMISS
    out = s.output
    out.write("\nKruskal-Wallis equality-of-populations rank test\n\n", "text")
    w0 = max(len(g.name), max(len(str(lab)) for lab, _, _ in rows))
    w1 = max(3, max(len(str(n)) for _, n, _ in rows))
    w2 = max(8, max(len(f"{rs:.2f}") for _, _, rs in rows))
    border = "  +" + "-" * (w0 + w1 + w2 + 8) + "+\n"
    out.write(border, "text")
    out.write(f"  | {g.name:>{w0}} | {'Obs':>{w1}} | {'Rank Sum':>{w2}} |\n", "text")
    out.write(f"  |{'-' * (w0 + 2)}+{'-' * (w1 + 2)}+{'-' * (w2 + 2)}|\n", "text")
    for lab, n, rs in rows:
        out.write(f"  | {str(lab):>{w0}} | {n:>{w1}} | {rs:>{w2}.2f} |\n", "text")
    out.write(border, "text")
    out.write(f"\nchi-squared = {H:>9.3f} with {df} d.f.\n", "text")
    out.write(f"probability = {st.chi2.sf(H, df):>10.4f}\n", "text")
    out.write(f"\nchi-squared with ties = {Ht:>9.3f} with {df} d.f.\n", "text")
    out.write(f"probability = {st.chi2.sf(Ht, df):>10.4f}\n", "text")
    _set_r(s, [("chi2_adj", Ht), ("df", float(df)), ("chi2", H)])


# ---------------------------------------------------------------------------
# spearman
# ---------------------------------------------------------------------------

@command("spearman", "spearman")
def cmd_spearman(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"stats": 2, "print": 2, "star": 2, "bonferroni": 1, "sidak": 2,
                                  "pw": 2, "matrix": 3})
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    if len(vars_) < 2:
        raise StataError(102, "too few variables specified")
    mask = touse(s, p)
    for v in vars_:
        mask &= S.valid(v.data)
    n = int(mask.sum())
    if n < 2:
        raise StataError(2001, "insufficient observations")
    R = np.column_stack([st.rankdata(v.data[mask]) for v in vars_])
    rho = np.corrcoef(R, rowvar=False)
    out = s.output

    def pval(r: float) -> float:
        if abs(r) >= 1:
            return 0.0
        t = r * math.sqrt((n - 2) / (1 - r * r))
        return float(2 * st.t.sf(abs(t), n - 2))

    if len(vars_) == 2 and not o.get("matrix"):
        r = float(rho[0, 1])
        pv = pval(r)
        out.write("\n", "text")
        # compat 0801: rótulos em 14 colunas, valores em %8.0f e %12.4f
        out.write(f"{'Number of obs':>14} = ", "text")
        out.write(f"{n:>7}\n", "result")
        out.write("Spearman's rho = ", "text")
        out.write(f"{r:>12.4f}\n", "result")
        out.write(f"\nTest of Ho: {vars_[0].name} and {vars_[1].name} are independent\n", "text")
        out.write(f"{'Prob > |t|':>14} = ", "text")
        out.write(f"{pv:>12.4f}\n", "result")
        _set_r(s, [("N", float(n)), ("rho", r), ("p", pv)])
        return
    out.write(f"(obs={n:,})\n", "text")
    names = [v.name for v in vars_]
    k = len(names)
    cells = [[f"{rho[i, j]:.4f}" for j in range(k)] for i in range(k)]
    extra = None
    if o.get("stats") or o.get("print"):
        extra = [[[f"{pval(rho[i, j]):.4f}" if i != j else "" for j in range(k)] for i in range(k)]]
    _corr_table(s, names, cells, extra=extra)
    from .matrix import Matrix
    s.r = {"N": float(n), "Rho": Matrix(rho.copy(), names, names)}


# ---------------------------------------------------------------------------
# swilk e sktest
# ---------------------------------------------------------------------------

def _swilk_stats(x: np.ndarray) -> tuple[float, float, float, float]:
    """W (algoritmo AS R94, o mesmo do Stata), V, z e Prob>z (Royston 1992)."""
    return _royston(float(st.shapiro(x).statistic), len(x))


def _royston(W: float, n: int) -> tuple[float, float, float, float]:
    ln = math.log(n)
    if n >= 12:
        mu = 0.0038915 * ln ** 3 - 0.083751 * ln ** 2 - 0.31082 * ln - 1.5861
        sigma = math.exp(0.0030302 * ln ** 2 - 0.082676 * ln - 0.4803)
        y = math.log(1 - W) if W < 1 else -math.inf
        V = (1 - W) / math.exp(mu)
    else:
        gamma = 0.459 * n - 2.273
        mu = 0.544 - 0.39978 * n + 0.025054 * n ** 2 - 0.0006714 * n ** 3
        sigma = math.exp(1.3822 - 0.77857 * n + 0.062767 * n ** 2 - 0.0020322 * n ** 3)
        y = -math.log(gamma - math.log(1 - W)) if W < 1 else -math.inf
        V = (1 - W) / math.exp(mu)   # VERIFICAR: V para n < 12
    z = (y - mu) / sigma
    return W, V, z, float(st.norm.sf(z))


@command("swilk", "swilk")
def cmd_swilk(s: "Session", args: str) -> None:
    p = parse_standard(args)
    match_options(p.options, {"generate": 1, "lnnormal": 2, "noties": 5})
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    base = touse(s, p)
    out = s.output
    out.write("\n" + " " * 19 + "Shapiro-Wilk W test for normal data\n\n", "text")
    out.write(f"    Variable |{'Obs':>11}{'W':>8}{'V':>12}{'z':>10}{'Prob>z':>13}\n", "text")
    out.write("-" * 13 + "+" + "-" * 54 + "\n", "text")
    for v in vars_:
        x = v.data[base & S.valid(v.data)]
        n = len(x)
        if n < 4:
            raise StataError(2001, "insufficient observations")
        W, V, z, pz = _swilk_stats(x)
        out.write(f"{name12(v.name):>12} |", "text")
        out.write(f"{n:>11}{W:>11.5f}{V:>11.3f}{z:>10.3f}{pz:>11.5f}\n", "result")
        _set_r(s, [("p", pz), ("z", z), ("V", V), ("W", W), ("N", float(n))])


def _sk_z(x: np.ndarray) -> tuple[float, float]:
    """Estatísticas z de assimetria (D'Agostino 1970) e curtose (Anscombe e
    Glynn 1983), como em [R] sktest."""
    n = len(x)
    d = x - x.mean()
    m2, m3, m4 = (d ** 2).mean(), (d ** 3).mean(), (d ** 4).mean()
    b1 = m3 / m2 ** 1.5
    b2 = m4 / m2 ** 2
    Y = b1 * math.sqrt((n + 1) * (n + 3) / (6 * (n - 2)))
    beta2 = 3 * (n * n + 27 * n - 70) * (n + 1) * (n + 3) / ((n - 2) * (n + 5) * (n + 7) * (n + 9))
    W2 = -1 + math.sqrt(2 * (beta2 - 1))
    delta = 1 / math.sqrt(math.log(math.sqrt(W2)))
    alpha = math.sqrt(2 / (W2 - 1))
    z1 = delta * math.log(Y / alpha + math.sqrt((Y / alpha) ** 2 + 1))
    E = 3 * (n - 1) / (n + 1)
    var = 24 * n * (n - 2) * (n - 3) / ((n + 1) ** 2 * (n + 3) * (n + 5))
    xk = (b2 - E) / math.sqrt(var)
    sb1 = 6 * (n * n - 5 * n + 2) / ((n + 7) * (n + 9)) * math.sqrt(6 * (n + 3) * (n + 5) / (n * (n - 2) * (n - 3)))
    A = 6 + 8 / sb1 * (2 / sb1 + math.sqrt(1 + 4 / sb1 ** 2))
    t = (1 - 2 / A) / (1 + xk * math.sqrt(2 / (A - 4)))
    z2 = ((1 - 2 / (9 * A)) - math.copysign(abs(t) ** (1 / 3), t)) / math.sqrt(2 / (9 * A))
    return z1, z2


@command("sktest", "sktest")
def cmd_sktest(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"noadjust": 5})
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    base = touse(s, p)
    out = s.output
    out.write("\n" + " " * 20 + "Skewness/Kurtosis tests for Normality\n", "text")
    out.write(" " * 58 + "------ joint ------\n", "text")
    lab = "chi2(2)" if o.get("noadjust") else "adj chi2(2)"
    out.write(f"    Variable |{'Obs':>11}  Pr(Skewness)  Pr(Kurtosis){lab:>12}   Prob>chi2\n", "text")
    out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
    for v in vars_:
        x = v.data[base & S.valid(v.data)]
        n = len(x)
        if n < 8:
            raise StataError(2001, "insufficient observations")
        z1, z2 = _sk_z(x)
        ps, pk = 2 * float(st.norm.sf(abs(z1))), 2 * float(st.norm.sf(abs(z2)))
        K2 = z1 * z1 + z2 * z2
        if o.get("noadjust"):
            pj = float(st.chi2.sf(K2, 2))
        else:
            K2, pj = _royston_adjust(K2, n)
        out.write(f"{name12(v.name):>12} |", "text")
        k2t = "." if not math.isfinite(K2) else f"{K2:.2f}"
        out.write(f"{n:>11}{ps:>11.4f}{pk:>14.4f}{k2t:>12}{pj:>15.4f}\n", "result")
        s.r = {"P_chi2": pj, "chi2": K2, "P_kurt": pk, "P_skew": ps, "N": float(n)}


# ---------------------------------------------------------------------------
# sdtest
# ---------------------------------------------------------------------------

def _set_r(s: "Session", pairs: list) -> None:
    """r() na ordem em que o Stata lista (return list mostra na ordem
    inversa da gravação)."""
    s.r = dict(reversed(pairs))


def _royston_adjust(K2: float, n: int) -> tuple[float, float]:
    """Correção de Royston (1991c) do K² de D'Agostino ([R] sktest, Methods
    and formulas). Devolve (adj chi2, P)."""
    Zc = -float(st.norm.ppf(math.exp(-0.5 * K2)))
    ln = math.log(n)
    Zt = 0.55 * n ** 0.2 - 0.21
    a1 = (-5 + 3.46 * ln) * math.exp(-1.37 * ln)
    b1 = 1 + (0.854 - 0.148 * ln) * math.exp(-0.55 * ln)
    a2 = a1 - (2.13 / (1 - 2.37 * ln)) * Zt
    b2 = 2.13 / (1 - 2.37 * ln) + b1
    Z = Zc if Zc < -1 else (a1 + b1 * Zc if Zc < Zt else a2 + b2 * Zc)
    P = float(st.norm.sf(Z))
    return (-2 * math.log(P) if P > 0 else math.inf), P


def _ha3(s: "Session", has: list[str], probs: list[str]) -> None:
    """Linhas Ha/Pr do sdtest (compat 0801): Pr nas colunas 2, 29 e 60; cada
    Ha centrado em 10,5, 39,5 e 69,5 (arredondado para baixo)."""
    starts = [2, 29, 60]
    centers = [10.5, 39.5, 69.5]
    line_p = [" "] * 90
    line_h = [" "] * 90
    for st_, c, pr, ha in zip(starts, centers, probs, has):
        for k, ch in enumerate(pr):
            line_p[st_ + k] = ch
        hs = int(math.floor(c - len(ha) / 2))
        for k, ch in enumerate(ha):
            line_h[hs + k] = ch
    s.output.write("\n" + "".join(line_h).rstrip() + "\n", "text")
    s.output.write("".join(line_p).rstrip() + "\n", "result")


def _row(name, x):
    n = len(x)
    m = S.moments(x)
    se = m.sd / math.sqrt(n) if n > 1 else M.SYSMISS
    lo, hi = S.t_ci(m.mean, se, n - 1, 95.0) if n > 1 else (M.SYSMISS, M.SYSMISS)
    return [name, n, m.mean, se, m.sd, lo, hi]


@command("sdtest", "sdtest")
def cmd_sdtest(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "level": 1})
    level = _level(s, o)
    if p.exp is None and "==" in p.varlist:
        a, _, b = p.varlist.partition("==")
        p.varlist, p.exp = a.strip(), b.strip()
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)

    def row(name, x):
        r = _row(name, x)
        if r[1] > 1:
            r[5], r[6] = S.t_ci(r[2], r[3], r[1] - 1, level)
        return tuple(r)

    out = s.output
    if isinstance(o.get("by"), str):
        g, groups = _groups(s, o["by"], mask)
        if len(groups) != 2:
            raise StataError(420, f"{len(groups)} groups found, 2 required")   # VERIFICAR
        x1, x2 = v.data[groups[0][2]], v.data[groups[1][2]]
        allx = np.concatenate([x1, x2])
        _ttest_table(s, [row(groups[0][1], x1), row(groups[1][1], x2)], level, title="Variance ratio test",
                     foot_rows=[row("combined", allx)], head="Group")
        sd1, sd2 = S.moments(x1).sd, S.moments(x2).sd
        f = sd1 ** 2 / sd2 ** 2
        df1, df2 = len(x1) - 1, len(x2) - 1
        out.write(f"    ratio = sd({groups[0][1]}) / sd({groups[1][1]})".ljust(66) + f"f = {f:>8.4f}\n", "text")
        out.write("Ho: ratio = 1".ljust(49) + "degrees of freedom = " + f"{df1}, {df2}".rjust(8) + "\n", "text")
        lo = float(st.f.cdf(f, df1, df2))
        up = float(st.f.sf(f, df1, df2))
        two = 2 * min(lo, up)
        mid = f"2*Pr(F < f) = {two:.4f}" if lo < up else f"2*Pr(F > f) = {two:.4f}"
        _ha3(s, ["Ha: ratio < 1", "Ha: ratio != 1", "Ha: ratio > 1"],
             [f"Pr(F < f) = {lo:.4f}", mid, f"Pr(F > f) = {up:.4f}"])
        _set_r(s, [("sd", S.moments(allx).sd), ("p", two), ("p_l", lo), ("p_u", up), ("sd_2", sd2),
                   ("sd_1", sd1), ("F", f), ("df_2", float(df2)), ("df_1", float(df1)),
                   ("N", float(len(allx)))])
        return
    if p.exp is None:
        raise StataError(198, "option by() required")   # VERIFICAR
    rhs = p.exp.strip()
    x = v.data[mask]
    if not _is_num(rhs):
        other = s.data.get(s.expand_varlist(rhs)[0])
        m2 = mask & S.valid(other.data)
        x1, x2 = v.data[m2], other.data[m2]
        allx = np.concatenate([x1, x2])
        _ttest_table(s, [row(v.name, x1), row(other.name, x2)], level, title="Variance ratio test",
                     foot_rows=[row("combined", allx)])
        sd1, sd2 = S.moments(x1).sd, S.moments(x2).sd
        f = sd1 ** 2 / sd2 ** 2
        df1 = df2 = len(x1) - 1
        out.write(f"    ratio = sd({v.name}) / sd({other.name})".ljust(66) + f"f = {f:>8.4f}\n", "text")
        out.write("Ho: ratio = 1".ljust(49) + "degrees of freedom = " + f"{df1}, {df2}".rjust(8) + "\n", "text")
        lo, up = float(st.f.cdf(f, df1, df2)), float(st.f.sf(f, df1, df2))
        two = 2 * min(lo, up)
        mid = f"2*Pr(F < f) = {two:.4f}" if lo < up else f"2*Pr(F > f) = {two:.4f}"
        _ha3(s, ["Ha: ratio < 1", "Ha: ratio != 1", "Ha: ratio > 1"],
             [f"Pr(F < f) = {lo:.4f}", mid, f"Pr(F > f) = {up:.4f}"])
        s.r = {"p": two, "p_u": up, "p_l": lo, "F": f, "df_2": float(df2), "df_1": float(df1),
               "sd_2": sd2, "sd_1": sd1, "N_2": float(len(x2)), "N_1": float(len(x1))}
        return
    sd0 = float(s.eval(rhs))
    n = len(x)
    _ttest_table(s, [row(v.name, x)], level, title="One-sample test of variance")
    sd = S.moments(x).sd
    c = (n - 1) * sd ** 2 / sd0 ** 2
    num = format(sd0, "g")
    out.write(f"    sd = sd({v.name})".ljust(59) + f"c = chi2 = {c:>8.4f}\n", "text")
    out.write(f"Ho: sd = {num}".ljust(49) + f"degrees of freedom = {n - 1:>8}\n", "text")
    lo, up = float(st.chi2.cdf(c, n - 1)), float(st.chi2.sf(c, n - 1))
    two = 2 * min(lo, up)
    mid = f"2*Pr(C < c) = {two:.4f}" if lo < up else f"2*Pr(C > c) = {two:.4f}"
    _ha3(s, [f"Ha: sd < {num}", f"Ha: sd != {num}", f"Ha: sd > {num}"],
         [f"Pr(C < c) = {lo:.4f}", mid, f"Pr(C > c) = {up:.4f}"])
    s.r = {"p": two, "p_u": up, "p_l": lo, "chi2": c, "df": float(n - 1), "sd": sd, "N": float(n)}
