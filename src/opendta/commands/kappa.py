"""Interrater agreement ([R] kappa): kap (two raters, one variable each; or
several raters, one variable each, rating the same outcomes) and kappa
(frequency form: one variable per outcome holding how many raters chose it).

Two raters: Cohen's kappa, optionally weighted (wgt(w) linear, wgt(w2)
quadratic), with the standard error under the null of Fleiss, Cohen and
Everitt (1969). More raters: Fleiss's kappa per outcome and combined, with the
null standard errors of Fleiss, Nee and Landis (1979). Layout details are
marked VERIFY and checked in compat/do/1402_kappa_roc.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from ._util import touse
from .epitab import _set_r
from .registry import command
from .summarize import name12, weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

NAN = float("nan")
SYS = M.SYSMISS


def _num(v: float) -> str:
    return f"{v:g}"


def cohen(table: np.ndarray, W: np.ndarray) -> tuple[float, float, float, float]:
    """(observed agreement, expected agreement, kappa, standard error under H0) for a KxK table
    of counts and a KxK agreement-weight matrix (identity for unweighted kappa)."""
    n = table.sum()
    p = table / n
    r, c = p.sum(1), p.sum(0)
    po = float((W * p).sum())
    pe = float((W * np.outer(r, c)).sum())
    kappa = (po - pe) / (1 - pe) if pe < 1 else NAN
    wr = W @ c          # mean weight of row category i
    wc = W.T @ r        # mean weight of column category j
    term = float((np.outer(r, c) * (W - (wr[:, None] + wc[None, :])) ** 2).sum())
    var0 = (term - pe * pe) / (n * (1 - pe) ** 2) if pe < 1 else NAN
    return po, pe, kappa, math.sqrt(var0) if var0 > 0 else NAN


def fleiss(counts: np.ndarray) -> tuple[list[float], list[float], float, float, bool]:
    """Fleiss's kappa from an N x K matrix of rater counts. Returns per-outcome kappas and null
    standard errors, the combined kappa and its null standard error, and whether the number
    of raters is constant (the standard errors need it)."""
    N, K = counts.shape
    m_i = counts.sum(1)
    const = bool(np.all(m_i == m_i[0]))
    m = float(m_i.mean())
    pj = counts.sum(0) / m_i.sum()
    qj = 1 - pj
    kj = []
    for j in range(K):
        num = float((counts[:, j] * (m_i - counts[:, j])).sum())
        den = float((m_i * (m_i - 1)).sum()) * pj[j] * qj[j]
        kj.append(1 - num / den if den else NAN)
    pq = float((pj * qj).sum())
    num = float((counts * (m_i[:, None] - counts)).sum())
    den = float((m_i * (m_i - 1)).sum()) * pq
    kappa = 1 - num / den if den else NAN
    if const and m > 1 and pq > 0:
        se_j = [math.sqrt(2 / (N * m * (m - 1)))] * K
        inner = pq * pq - float((pj * qj * (qj - pj)).sum())
        se = math.sqrt(2) / (pq * math.sqrt(N * m * (m - 1))) * math.sqrt(max(inner, 0.0))
    else:
        se_j, se = [NAN] * K, NAN
    return kj, se_j, kappa, se, const


def _zp(k: float, se: float) -> tuple[float, float]:
    if not se > 0 or math.isnan(k):
        return NAN, NAN
    z = k / se
    return z, float(st.norm.sf(z))


def _fmt(x: float, f: str) -> str:
    return "." if x is None or math.isnan(x) else format(x, f)


@command("kap", "kap")
def cmd_kap(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) < 2:
        raise StataError(198, "too few variables specified")
    o = match_options(p.options, {"tab": 3, "wgt": 1, "absolute": 3})
    ds = s.data
    cols = [ds.get(n) for n in names]
    if any(c.is_string for c in cols):
        raise StataError(109, "type mismatch")
    mask = touse(s, p)
    for c in cols:
        mask &= c.data < SYS
    w, _, mask = weights(s, p, mask, ("fweight",))
    w_given = w
    if w is None:
        w = np.ones(ds.nobs)
    levels = sorted(set(np.concatenate([c.data[mask] for c in cols]).tolist()))
    K = len(levels)
    if K < 2:
        raise StataError(459, "ratings do not vary")   # VERIFY
    if len(names) > 2:
        pos = {v: i for i, v in enumerate(levels)}
        idx = np.flatnonzero(mask)
        counts = np.zeros((len(idx), K))
        for r, i in enumerate(idx):
            for c in cols:
                counts[r, pos[float(c.data[i])]] += 1
        counts = np.repeat(counts, w[idx].astype(int), axis=0)
        s.output.write(f"\nThere are {len(names)} raters per subject:\n", "text")
        _fleiss_out(s, counts, [_num(v) for v in levels])
        return
    T = np.zeros((K, K))
    pos = {v: i for i, v in enumerate(levels)}
    for i in np.flatnonzero(mask):
        T[pos[float(cols[0].data[i])], pos[float(cols[1].data[i])]] += w[i]
    W = np.eye(K)
    wname = o.get("wgt")
    if isinstance(wname, str):
        x = np.asarray(levels, float) if o.get("absolute") else np.arange(K, dtype=float)
        span = x.max() - x.min()
        if wname == "w":
            W = 1 - np.abs(x[:, None] - x[None, :]) / span
        elif wname == "w2":
            W = 1 - ((x[:, None] - x[None, :]) / span) ** 2
        elif wname in getattr(s, "matrices", {}):
            W = np.asarray(s.matrices[wname].data, float)
        else:
            raise StataError(111, f"weight {wname} not found")   # VERIFY (kapwgt)
    out = s.output
    if o.get("tab"):
        # the two-way table of tabulate (compat 1402)
        from .tabulate import _twoway
        _twoway(s, cols[0], cols[1], mask, w_given, "fweight" if w_given is not None else "", {})
    if isinstance(wname, str):
        out.write("\nRatings weighted by:\n", "text")
        for i in range(K):
            out.write("".join(f"{W[i, j]:>9.4f}" for j in range(K)) + "\n", "result")
    po, pe, k, se = cohen(T, W)
    z, pv = _zp(k, se)
    out.write("\n" + " " * 13 + "Expected\n", "text")
    out.write("Agreement   Agreement     Kappa   Std. Err.         Z      Prob>Z\n", "text")
    out.write("-" * 65 + "\n", "text")
    out.write(f"{100 * po:>7.2f}%{100 * pe:>11.2f}%{_fmt(k, '>11.4f')}{_fmt(se, '>11.4f')}"
              f"{_fmt(z, '>11.2f')}{_fmt(pv, '>12.4f')}\n", "result")
    # return list order of Stata 14 (compat 1402)
    _set_r(s, [("N", T.sum()), ("prop_o", po), ("prop_e", pe), ("kappa", k), ("z", z), ("se", se)])


def _fleiss_out(s: "Session", counts: np.ndarray, labels: list[str]) -> None:
    kj, sej, k, se, const = fleiss(counts)
    out = s.output
    K = counts.shape[1]
    if K == 2:
        z, pv = _zp(k, se)
        out.write("\n" + " " * 9 + "Kappa        Z     Prob>Z\n", "text")
        out.write(" " * 7 + "-" * 27 + "\n", "text")
        out.write(f"{_fmt(k, '>14.4f')}{_fmt(z, '>9.2f')}{_fmt(pv, '>11.4f')}\n", "result")
        _set_r(s, [("prop_e", NAN), ("prop_o", NAN), ("se", se), ("kappa", k), ("z", z), ("N", counts.shape[0])])
        return
    out.write("\n" + f"{'Outcome':>16} |    Kappa          Z     Prob>Z\n", "text")
    out.write("-" * 17 + "+" + "-" * 31 + "\n", "text")
    for lab, kk, ss in zip(labels, kj, sej):
        z, pv = _zp(kk, ss)
        out.write(f"{name12(lab):>16} |", "text")
        out.write(f"{_fmt(kk, '>10.4f')}{_fmt(z, '>11.2f')}{_fmt(pv, '>10.4f')}\n", "result")
    out.write("-" * 17 + "+" + "-" * 31 + "\n", "text")
    z, pv = _zp(k, se)
    out.write(f"{'combined':>16} |", "text")
    out.write(f"{_fmt(k, '>10.4f')}{_fmt(z, '>11.2f')}{_fmt(pv, '>10.4f')}\n", "result")
    if not const:
        out.write("\nNote: Number of ratings per subject vary; cannot calculate test statistics.\n", "text")
    # Stata 14 returns only these two (compat 1402)
    _set_r(s, [("kappa", k), ("z", z)])


@command("kappa", "kappa")
def cmd_kappa(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) < 2:
        raise StataError(198, "too few variables specified")
    match_options(p.options, {})
    ds = s.data
    cols = [ds.get(n) for n in names]
    mask = touse(s, p)
    for c in cols:
        mask &= c.data < SYS
    counts = np.column_stack([c.data[mask] for c in cols]).astype(float)
    if np.any(counts < 0):
        raise StataError(459, "negative counts not allowed")   # VERIFY
    m_i = counts.sum(1)
    if np.any(m_i < 2):
        raise StataError(459, "each subject needs at least two ratings")   # VERIFY
    _fleiss_out(s, counts, names)
