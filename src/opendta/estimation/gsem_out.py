"""Output and prediction after gsem: the header, the table with the equations, the variances and
covariances of the latent variables and the constraints, and predict."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lazy import lazy
from .gsem import GModel, _FAM_TITLE, _LINK_TITLE
from .results import Estimates, crit, g9

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

SYS = M.SYSMISS


def _auxlabel(m: GModel, i: int) -> str:
    """/cut1 when one response has it, cut1:y when several do (VERIFY)."""
    P = m.P
    same = [j for j in range(len(P.name)) if P.name[j] == P.name[i] and P.kind[j] == P.kind[i]]
    return f"/{P.name[i]}" if len(same) == 1 else f"{P.name[i]}:{P.eq[i]}"


def _constraint_lines(m: GModel):
    lines = []
    for r in m.resps:
        nm = r.eq.name
        for t in r.eq.terms:
            pass
        for j, (i, v) in enumerate(r.cref):
            if i == -1:
                L = m.lat[j]
                lab = f"{L}[{m.cvar}]" if m.cvar else L
                lines.append(f"[{nm}]{lab} = {v:g}")
        for k, (i, v) in enumerate(r.xref):
            if i < 0:
                xs = [t for t in r.eq.terms if t[0] == "x"]
                lines.append(f"[{nm}]{xs[k][1]} = {v:g}")
    for j, (i, v) in enumerate(m.vref):
        if i < 0:
            lines.append(f"var({m.lat[j]}) = {v:g}")
    # shared labels
    by_param: dict[int, list[str]] = {}
    for r in m.resps:
        for j, (i, _v) in enumerate(r.cref):
            if i >= 0:
                by_param.setdefault(i, []).append(f"[{r.eq.name}]{m.lat[j]}")
        xs = [t for t in r.eq.terms if t[0] == "x"]
        for k, (i, _v) in enumerate(r.xref):
            if i >= 0 and k < len(xs):
                by_param.setdefault(i, []).append(f"[{r.eq.name}]{xs[k][1]}")
    for i, items in by_param.items():
        for other in items[1:]:
            lines.append(f"{items[0]} - {other} = 0")
    return lines


def display_gsem(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                 level: float | None = None, **_kw) -> None:
    out = s.output
    m: GModel = est.model
    x = est.extra
    lv = level or est.level
    if header:
        out.write("\n" + f"{'Generalized structural equation model':<48}{'Number of obs':<18}= ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        many = len(m.resps) > 1
        for r in m.resps:
            e = r.eq
            right = [""] * 4
            if e.family == "gamma":
                # VERIFY: the count and the time at risk with several responses
                right[0] = f"{'No. of failures':<18}= {len(r.y):>10,}"
                right[1] = f"{'Time at risk':<18}= {g9(float(np.sum(r.y)), '%10.0g'):>10}"
            lines_ = [f"{'Response':<15}: {e.name}", f"{'Family':<15}: {_FAM_TITLE[e.family]}"]
            if e.family == "nbinomial":
                lines_.append(f"{'Dispersion':<15}: mean")
            if e.family == "gamma":
                lines_.append(f"{'Form':<15}: accelerated failure-time")
            lines_.append(f"{'Link':<15}: {_LINK_TITLE[e.link]}")
            if many:
                out.write("\n", "text")
            for k, ln in enumerate(lines_):
                out.write((f"{ln:<48}{right[k]}" if k < len(right) and right[k] else ln) + "\n", "text")
        if many:
            out.write("\n", "text")
        out.write("Log likelihood = ", "text")
        out.write(f"{g9(x['ll'], '%10.0g'):>10}\n", "result")
        out.write("\n", "text")
        for k, ln in enumerate(_constraint_lines(m), 1):
            out.write(f" ({k:2d})  {ln}\n", "text")
    if not table:
        return
    rep, Vr = x["rep"], x["Vr"]
    se = np.sqrt(np.maximum(np.diag(Vr), 0))
    c = crit("z", lv, None)
    P = m.P
    labels = []
    for i, kd in enumerate(P.kind):
        if kd in ("lvar",):
            labels.append(f"var({P.eq[i][4:-1]})")
        elif kd == "gvar":
            labels.append(f"var(e.{P.eq[i]})")
        elif kd == "cov":
            labels.append(P.eq[i])
        elif kd in ("cut", "aux"):
            labels.append(_auxlabel(m, i))
    W = max([12] + [len(lb) - 1 for lb in labels])
    dash = "-" * (W + 1) + "+" + "-" * 64 + "\n"
    out.write("-" * (W + 66) + "\n", "text")
    out.write(f"{'':>{W}} |      Coef.   Std. Err.      z    P>|z|     [{lv:g}% Conf. Interval]\n", "text")
    out.write(dash, "text")

    def row(label, b, s_, kind="coef"):
        if kind in ("var", "cov"):
            out.write(f"{label:>{W + 1}}|", "text")
        else:
            out.write(f"{label:>{W}} |", "text")
        if kind == "constrained":
            out.write(f"  {g9(b):>9}  (constrained)\n", "result")
            return
        if s_ > 0 and np.isfinite(s_):
            if kind == "var":
                try:
                    lo, hi = b * math.exp(-c * s_ / b), b * math.exp(c * s_ / b)
                except OverflowError:
                    lo, hi = float("nan"), float("nan")
                out.write(f"  {g9(b):>9}  {g9(s_):>9}" + " " * 17 + f"    {g9(lo):>9}   {g9(hi):>9}\n",
                          "result")
                return
            z = b / s_
            pv = 2 * float(st.norm.sf(abs(z)))
            out.write(f"  {g9(b):>9}  {g9(s_):>9}  {z:>7.2f}  {pv:>6.3f}    {g9(b - c * s_):>9}"
                      f"   {g9(b + c * s_):>9}\n", "result")
        elif kind == "var":
            out.write(f"  {g9(b):>9}  {'.':>9}" + " " * 17 + f"    {'.':>9}   {'.':>9}\n", "result")
        else:
            out.write(f"  {g9(b):>9}  {'.':>9}\n", "result")

    first = True
    for e in m.eqs:
        if not first:
            out.write(dash, "text")
        first = False
        out.write(f"{e.name + ' <-':<{W}} |\n", "text")
        if e.kind == "resp":
            r = next(rr for rr in m.resps if rr.eq is e)
            xi = 0
            for t in e.terms:
                if t[0] == "x":
                    i, v = r.xref[xi]
                    xi += 1
                    row(t[1], rep[i], se[i]) if i >= 0 else row(t[1], v, 0, "constrained")
                else:
                    j = m.lat.index(t[1])
                    i, v = r.cref[j]
                    lab = f"{t[1]}[{m.cvar}]" if t[0] == "clat" else t[1]
                    # a cluster-level latent variable sits between blank lines (case 1909)
                    if t[0] == "clat":
                        out.write(f"{'':>{W}} |\n", "text")
                    row(lab, rep[i], se[i]) if i >= 0 else row(lab, v, 0, "constrained")
                    if t[0] == "clat":
                        out.write(f"{'':>{W}} |\n", "text")
            if len(r.xref) > xi:
                i = r.xref[xi][0]
                row("_cons", rep[i], se[i])
            own = [i for i in r.aux if P.kind[i] in ("cut", "aux")]
            if own:
                # the ancillary parameters of the response, under its name (case 1909)
                out.write(dash, "text")
                out.write(f"{e.name:<{W}} |\n", "text")
                for i in own:
                    row(f"/{P.name[i]}", rep[i], se[i])
        else:
            j = m.lat.index(e.name)
            for k, nm in enumerate(m.gnames[j]):
                i, v = m.gref[j][k]
                row(nm, rep[i], se[i]) if i >= 0 else row(nm, v, 0, "constrained")
    rows_after = [i for i, kd in enumerate(P.kind) if kd == "lvar"]
    for i in rows_after:
        out.write(dash, "text") if i == rows_after[0] else None
        row(f"var({P.eq[i][4:-1]})", rep[i], se[i], "var")
    for j, (i, v) in enumerate(m.vref):
        if i < 0:
            row(f"var({m.lat[j]})", v, 0, "constrained")
    covs = [i for i, kd in enumerate(P.kind) if kd == "cov"]
    for i in covs:
        out.write(dash, "text") if i == covs[0] else None
        row(P.eq[i], rep[i], se[i], "cov")
    gvars = [i for i, kd in enumerate(P.kind) if kd == "gvar"]
    for n_, i in enumerate(gvars):
        out.write(dash, "text") if n_ == 0 else None
        row(f"var(e.{P.eq[i]})", rep[i], se[i], "var")
    out.write("-" * (W + 66) + "\n", "text")


def _predict(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    """predict newvar, mu(y) | xb(y) | latent(L): the mean and the linear predictor of a
    response (latent variables at their empirical Bayes means) and the latent means."""
    m: GModel = est.model
    th = est.extra["th"]
    opts = getattr(est, "_pred_opts", None) or {}
    if stat in ("mu", "pr", "eta", "xb") and isinstance(opts.get(stat), str):
        # Stata 14: predict newvar, mu outcome(y) (case 1909)
        raise StataError(198, f"option {stat}() incorrectly specified")
    if stat == "eta":
        stat = "xb"
    target = opts.get(stat) if stat == "latent" else opts.get("outcome")
    target = target.strip() if isinstance(target, str) and target.strip() else None
    n = s.data.nobs
    rows = np.nonzero(m.sample)[0]
    out = np.full(n, SYS)
    eb = m.posterior_means(th) if m.q else np.zeros((m.U, 0))
    if stat == "latent":
        L = target or (m.lat[0] if m.lat else None)
        if L not in m.lat:
            raise StataError(198, f"{L} is not a latent variable")   # VERIFY
        out[rows] = eb[m.unit, m.lat.index(L)]
        note = "" if opts.get("ebmeans") else "(option ebmeans assumed)\n"
        return out, note + f"(using {m.nint} quadrature points)"
    r = next((rr for rr in m.resps if rr.eq.name == target), None) if target else m.resps[0]
    if r is None:
        raise StataError(198, f"{target} is not an outcome of the model")   # VERIFY
    ds = s.data
    e = r.eq
    cols = [ds.get(t[1]).data[rows].astype(float) for t in e.terms if t[0] == "x"]
    beta = np.array([m._v(th, x) for x in r.xref])
    Xf = np.column_stack(cols + ([np.ones(m.N)] if len(r.xref) > len(cols) else [])) \
        if (cols or len(r.xref) > len(cols)) else np.zeros((m.N, 0))
    eta = (Xf @ beta if len(beta) else np.zeros(m.N))
    for v in (e.offset,):
        if v:
            eta = eta + ds.get(v).data[rows]
    if e.exposure:
        eta = eta + np.log(ds.get(e.exposure).data[rows])
    if stat in ("mu", "") and m.q:
        c = np.array([m._v(th, x) for x in r.cref])
        eta = eta + eb[m.unit] @ c
    if stat == "xb" or e.link == "identity":
        val = eta
    elif e.link == "log":
        val = np.exp(eta)
    elif e.link == "logit":
        val = special.expit(eta)
    elif e.link == "probit":
        val = special.ndtr(eta)
    else:
        val = -np.expm1(-np.exp(eta))
    if e.family == "ordinal" and stat != "xb":
        raise StataError(198, "predict mu is not available for ordinal outcomes")   # OpenDTA
    out[rows] = val
    return out, ""
