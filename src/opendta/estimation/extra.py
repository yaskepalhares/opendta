"""lrtest, hausman, xtsum, xtdescribe e qreg ([R] lrtest, [R] hausman,
[XT] xtsum, [XT] xtdescribe, [R] qreg).

Fórmulas dos manuais; layouts de memória com a marca VERIFICAR até a
conferência com o caso compat/do/0804_comandos.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from .results import CoefRow, Estimates, coef_table, post

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or x >= SYS or not np.isfinite(x):
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _opts(text: str, spec: dict) -> tuple[str, dict]:
    c = find_top(text, ",")
    body, opt = (text, "") if c == -1 else (text[:c], text[c + 1:])
    return body.strip(), (match_options(opt, spec) if opt.strip() else {})


def _get_stored(s: "Session", name: str) -> dict:
    """Resultados guardados (estimates store) ou os atuais ('.')."""
    from .postest import _snapshot, _stored
    if name in (".", ""):
        if not s.e or "cmd" not in s.e:
            raise StataError(301, "last estimates not found")
        return _snapshot(s)
    st_ = _stored(s)
    if name not in st_:
        raise StataError(111, f"estimation result {name} not found")
    return st_[name]


# ---------------------------------------------------------------------------
# lrtest
# ---------------------------------------------------------------------------

def cmd_lrtest(s: "Session", args: str) -> None:
    body, o = _opts(args, {"stats": 3, "dir": 3, "df": 2, "force": 5})
    names = body.split()
    if not names:
        raise StataError(198, "invalid syntax")
    if len(names) == 1:
        names = [names[0], "."]
    a, b = names[:2]
    ea, eb = _get_stored(s, a)["e"], _get_stored(s, b)["e"]
    for nm, e in ((a, ea), (b, eb)):
        if not isinstance(e.get("ll"), float):
            raise StataError(498, f"{nm} does not contain scalar e(ll)")   # VERIFICAR
    ka = float(ea.get("rank", ea.get("df_m", 0) + 1))
    kb = float(eb.get("rank", eb.get("df_m", 0) + 1))
    if ka == kb and not isinstance(o.get("df"), str):
        raise StataError(498, "df(unrestricted) = df(restricted) = %d" % ka)   # VERIFICAR
    if float(ea.get("N", 0)) != float(eb.get("N", 0)) and not o.get("force"):
        raise StataError(498, "observations differ: %d vs. %d" % (ea.get("N", 0), eb.get("N", 0)))
    # o modelo com mais parâmetros é o irrestrito
    (u, eu, ku), (r, er, kr) = ((a, ea, ka), (b, eb, kb)) if ka > kb else ((b, eb, kb), (a, ea, ka))
    df = int(o["df"]) if isinstance(o.get("df"), str) else int(round(ku - kr))
    chi2 = 2 * (eu["ll"] - er["ll"])
    p = chi2_tail(chi2, df) if chi2 >= 0 else 1.0
    out = s.output
    out.write("\n", "text")
    out.write(f"{'Likelihood-ratio test':<54}LR chi2({df}){'':<{max(1, 3 - len(str(df)))}}= ", "text")
    out.write(f"{chi2:>9.2f}\n", "result")
    out.write(f"{'(Assumption: ' + r + ' nested in ' + u + ')':<54}Prob > chi2 = ", "text")
    out.write(f"{p:>9.4f}\n", "result")
    s.r = dict(reversed([("p", p), ("chi2", chi2), ("df", float(df))]))   # ordem do compat 0804


# ---------------------------------------------------------------------------
# hausman
# ---------------------------------------------------------------------------

def cmd_hausman(s: "Session", args: str) -> None:
    body, o = _opts(args, {"constant": 4, "alleqs": 3, "sigmamore": 5, "sigmaless": 5, "df": 2,
                           "skipeqs": 4, "equations": 2, "force": 5})
    names = body.split()
    if len(names) == 1:
        names.append(".")
    if len(names) != 2:
        raise StataError(198, "invalid syntax")
    sa, sb = _get_stored(s, names[0]), _get_stored(s, names[1])
    ea, eb = sa["e"], sb["e"]
    ba, bb = ea["b"], eb["b"]
    Va, Vb = np.asarray(ea["V"].data, float), np.asarray(eb["V"].data, float)
    ca, cb = list(ba.colnames), list(bb.colnames)
    common = [c for c in ca if c in cb and (o.get("constant") or c != "_cons")
              and Va[ca.index(c), ca.index(c)] > 0 and Vb[cb.index(c), cb.index(c)] > 0]
    if not common:
        raise StataError(498, "no coefficients in common")   # VERIFICAR
    ia = [ca.index(c) for c in common]
    ib = [cb.index(c) for c in common]
    b1 = np.asarray(ba.data, float).ravel()[ia]
    b2 = np.asarray(bb.data, float).ravel()[ib]
    V1 = Va[np.ix_(ia, ia)]
    V2 = Vb[np.ix_(ib, ib)]
    if o.get("sigmamore") or o.get("sigmaless"):
        # VERIFICAR: escala de V pelo sigma do modelo eficiente (ou do consistente)
        sa2 = float(eb.get("sigma_e", eb.get("rmse", 1.0))) ** 2
        sc2 = float(ea.get("sigma_e", ea.get("rmse", 1.0))) ** 2
        f = sa2 / sc2 if o.get("sigmamore") else sc2 / sa2
        if o.get("sigmamore"):
            V1 = V1 * f
        else:
            V2 = V2 * f
    d = b1 - b2
    Vd = V1 - V2
    w, U = np.linalg.eigh((Vd + Vd.T) / 2)
    tol = max(abs(w).max(), 1e-300) * len(w) * 2.2e-16
    pos = w > tol
    rank = int(pos.sum())
    Vinv = (U[:, pos] / w[pos]) @ U[:, pos].T
    chi2 = float(d @ Vinv @ d)
    df = int(o["df"]) if isinstance(o.get("df"), str) else rank
    p = chi2_tail(chi2, df)
    out = s.output
    na = names[0] if names[0] != "." else "."
    nb = names[1] if names[1] != "." else "."
    out.write("\n" + " " * 17 + "---- Coefficients ----\n", "text")
    # larguras conferidas no compat 0804
    out.write(f"{'':>12} |{'(b)':>9}{'(B)':>13}{'(b-B)':>17}     sqrt(diag(V_b-V_B))\n", "text")
    out.write(f"{'':>12} |{na[:10]:>9}{nb[:10]:>13}{'Difference':>19}{'S.E.':>14}\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    diag = np.diag(Vd)
    for j, c in enumerate(common):
        se = math.sqrt(diag[j]) if diag[j] > 0 else SYS
        out.write(f"{c[:12]:>12} |", "text")
        out.write(f"{_g(b1[j]):>12}{_g(b2[j]):>13}{_g(d[j]):>16}{_g(se):>16}\n", "result")
    out.write("-" * 78 + "\n", "text")
    cmd_a = ea.get("cmd", "")
    cmd_b = eb.get("cmd", "")
    out.write(f"{'b = consistent under Ho and Ha; obtained from ' + cmd_a:>78}\n", "text")
    out.write(f"{'B = inconsistent under Ha, efficient under Ho; obtained from ' + cmd_b:>78}\n", "text")
    out.write("\n    Test:  Ho:  difference in coefficients not systematic\n\n", "text")
    out.write(f"{'chi2(' + str(df) + ')':>25} = (b-B)'[(V_b-V_B)^(-1)](b-B)\n", "text")
    out.write(f"{'':>25} = ", "text")
    out.write(f"{chi2:>11.2f}\n", "result")
    out.write(f"{'Prob>chi2':>25} = ", "text")
    out.write(f"{p:>11.4f}\n", "result")
    if rank < len(common):
        out.write(f"{'(V_b-V_B is not positive definite)':>78}\n", "text")
    s.r = dict(reversed([("rank", float(rank)), ("chi2", chi2), ("df", float(df)), ("p", p)]))   # compat 0804


# ---------------------------------------------------------------------------
# xtsum, xtdescribe
# ---------------------------------------------------------------------------

def _panel(s: "Session", opt_i: str | None = None):
    from .tsops import panel_info, ts_info
    ds = s.data
    pvar = opt_i or panel_info(ds)
    if not pvar:
        raise StataError(459, "panel variable not set; use xtset varname ...")
    ti = ts_info(ds)
    return pvar, (ti.tvar if ti and ti.tvar != pvar else "")


def cmd_xtsum(s: "Session", args: str) -> None:
    from ..commands._util import touse
    p = parse_standard(args)
    o = match_options(p.options, {"i": 1}) if p.options.strip() else {}
    pvar, _t = _panel(s, str(o["i"]).strip() if isinstance(o.get("i"), str) else None)
    ds = s.data
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else [n for n in ds.names if n != pvar]
    base = touse(s, p) & (ds.get(pvar).data < SYS)
    ids = ds.get(pvar).data
    out = s.output
    out.write("\n" + f"{'Variable':<17}|{'Mean':>10}{'Std. Dev.':>12}{'Min':>10}{'Max':>11} |{'Observations':>16}\n",
              "text")
    out.write("-" * 17 + "+" + "-" * 44 + "+" + "-" * 16 + "\n", "text")
    for k, n in enumerate(names):
        v = ds.get(n)
        if v.is_string:
            continue
        x = v.data.astype(float)
        m = base & (x < SYS)
        xs, gi = x[m], ids[m]
        N = len(xs)
        if N == 0:
            continue
        uniq, inv, cnt = np.unique(gi, return_inverse=True, return_counts=True)
        G = len(uniq)
        means = np.bincount(inv, weights=xs) / cnt
        grand = xs.mean()
        within = xs - means[inv] + grand
        sd_o = xs.std(ddof=1) if N > 1 else SYS
        sd_b = means.std(ddof=1) if G > 1 else SYS
        sd_w = within.std(ddof=1) if N > 1 else SYS
        balanced = bool(np.all(cnt == cnt[0]))
        tlab = "T" if balanced else "T-bar"
        tval = f"{int(cnt[0]):>7}" if balanced else f"{N / G:>7.4g}"
        rows = [(f"{n[:8]:<8} overall", _g(grand), _g(sd_o), _g(xs.min()), _g(xs.max()), "N", f"{N:>7}"),
                (f"{'':<8} between", "", _g(sd_b), _g(means.min()), _g(means.max()), "n", f"{G:>7}"),
                (f"{'':<8} within ", "", _g(sd_w), _g(within.min()), _g(within.max()), tlab, tval)]
        if k:
            out.write(f"{'':<17}|{'':<44}|\n", "text")
        for lab, mn, sd, lo, hi, ol, ov in rows:
            out.write(f"{lab:<16} |", "text")
            out.write(f"{mn:>10}{sd:>11}{lo:>11}{hi:>11}", "result")
            out.write(f" | {ol:>5} = ", "text")
            out.write(f"{ov}\n", "result")   # layout conferido no compat 0804


def cmd_xtdescribe(s: "Session", args: str) -> None:
    from ..commands._util import touse
    p = parse_standard(args)
    o = match_options(p.options, {"i": 1, "t": 1, "patterns": 1}) if p.options.strip() else {}
    pvar, tvar = _panel(s)
    if not tvar:
        raise StataError(459, "time variable not set; use xtset varname timevar")   # VERIFICAR
    ds = s.data
    m = touse(s, p) & (ds.get(pvar).data < SYS) & (ds.get(tvar).data < SYS)
    ids, tt = ds.get(pvar).data[m], ds.get(tvar).data[m]
    uniq_i = np.unique(ids)
    uniq_t = np.unique(tt)
    n, T = len(uniq_i), len(uniq_t)
    from .tsops import ts_info
    delta = ts_info(ds).delta if ts_info(ds) else 1.0
    span = int(round((uniq_t.max() - uniq_t.min()) / delta)) + 1
    out = s.output

    def lst(vals):
        v = [format_value(float(x), "%9.0g", pad=False).strip() for x in vals]
        return ", ".join(v) if len(v) <= 3 else f"{v[0]}, {v[1]}, ..., {v[-1]}"
    out.write("\n", "text")
    out.write(f"{pvar:>8}:  {lst(uniq_i):<50}n = ", "text")
    out.write(f"{n:>10}\n", "result")
    out.write(f"{tvar:>8}:  {lst(uniq_t):<50}T = ", "text")
    out.write(f"{T:>10}\n", "result")
    out.write(f"{'':>11}Delta({tvar}) = {delta:g} unit\n", "text")
    out.write(f"{'':>11}Span({tvar})  = {span} periods\n", "text")
    uniqpairs = len(set(zip(ids.tolist(), tt.tolist()))) == len(ids)
    if uniqpairs:
        out.write(f"{'':>11}({pvar}*{tvar} uniquely identifies each observation)\n", "text")
    else:
        out.write(f"{'':>11}({pvar}*{tvar} does not uniquely identify observations)\n", "text")
    counts = np.array([np.sum(ids == i) for i in uniq_i])
    qs = [np.min(counts)] + [S_pct(counts, q) for q in (5, 25, 50, 75, 95)] + [np.max(counts)]
    out.write("\nDistribution of T_i:   min      5%     25%       50%       75%     95%     max\n", "text")
    widths = (26, 8, 8, 10, 10, 8, 8)     # compat 0804
    out.write("".join(f"{int(q):>{wd}}" for wd, q in zip(widths, qs)) + "\n", "result")
    # padrões de presença (1 = observado, . = ausente) por período
    tindex = {t: k for k, t in enumerate(uniq_t)}
    pats: dict[str, int] = {}
    for i in uniq_i:
        row = ["."] * T
        for t in tt[ids == i]:
            row[tindex[t]] = "1"
        key = "".join(row)
        pats[key] = pats.get(key, 0) + 1
    order = sorted(pats.items(), key=lambda kv: (-kv[1], kv[0]))
    npat = int(o["patterns"]) if isinstance(o.get("patterns"), str) else 9
    out.write("\n     Freq.  Percent    Cum. |  Pattern\n", "text")
    out.write(" " + "-" * 27 + "+" + "-" * (T + 2) + "\n", "text")
    cum = 0.0
    shown = order[:npat]
    for key, f in shown:
        pct = 100 * f / n
        cum += pct
        out.write(f"{f:>9}{pct:>10.2f}{cum:>8.2f} |  {key}\n", "result")
    rest = sum(f for _, f in order[npat:])
    if rest:
        pct = 100 * rest / n
        cum += pct
        out.write(f"{rest:>9}{pct:>10.2f}{100.0:>8.2f} | (other patterns)\n", "result")
    out.write(" " + "-" * 27 + "+" + "-" * (T + 2) + "\n", "text")
    out.write(f"{n:>9}{100.0:>10.2f}{'':>8} |  {'X' * T}\n", "result")
    s.r = dict(reversed([("n", float(n)), ("N", float(len(ids))), ("max", float(np.max(counts))),
                         ("min", float(np.min(counts)))]))   # VERIFICAR


def S_pct(x: np.ndarray, q: float) -> float:
    from ..core.stats import percentile
    return percentile(np.sort(x.astype(float)), q)


# ---------------------------------------------------------------------------
# qreg
# ---------------------------------------------------------------------------

def _qfit(X: np.ndarray, y: np.ndarray, tau: float, w: np.ndarray | None = None) -> np.ndarray:
    """min Σ w ρ_τ(y - Xb) por programação linear (HiGHS)."""
    from scipy.optimize import linprog
    n, k = X.shape
    ww = np.ones(n) if w is None else w
    c = np.concatenate([np.zeros(2 * k), tau * ww, (1 - tau) * ww])
    A = np.hstack([X, -X, np.eye(n), -np.eye(n)])
    res = linprog(c, A_eq=A, b_eq=y, bounds=[(0, None)] * (2 * k + 2 * n), method="highs")
    if not res.success:
        raise StataError(430, "convergence not achieved")
    z = res.x
    return z[:k] - z[k:2 * k]


def cmd_qreg(s: "Session", args: str) -> None:
    from .regress import Sample
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "qreg", args)
        return
    smp = Sample(s, args, opts_spec={"quantile": 1, "vce": 3, "level": 1, "wlsiter": 3, "nolog": 4},
                 weights=("aweight", "fweight", "iweight", "pweight"))
    o = smp.o
    q = float(o["quantile"]) if isinstance(o.get("quantile"), str) else 0.5
    if q >= 1:
        q = q / 100
    if not 0 < q < 1:
        raise StataError(198, "quantile() must be between 0 and 1")   # VERIFICAR
    from . import fvars
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [c for c in smp.cols if c.values is not None]
    y = smp.y.astype(float)
    n = len(y)
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(n)] if smp.constant else []))
    k = X.shape[1]
    w = smp.w[smp.mask].astype(float) if smp.w is not None else None
    b = _qfit(X, y, q, w)
    r = y - X @ b
    ww = np.ones(n) if w is None else w
    rho = lambda u, t: u * (t - (u < 0))   # noqa: E731
    min_dev = float(ww @ rho(r, q))
    # desvio do modelo só com a constante: quantil de y
    # o Stata usa a observação de ordem ceil(n·q) de y (com pesos, a primeira
    # cujo peso acumulado chega a q·W): 5.36 e 3.44 no compat 0804
    o_ = np.argsort(y, kind="stable")
    cw = np.cumsum(ww[o_])
    yq = float(y[o_][min(int(np.searchsorted(cw, q * cw[-1] - 1e-9 * cw[-1])), n - 1)])
    raw_dev = float(ww @ rho(y - yq, q))
    pr2 = 1 - min_dev / raw_dev if raw_dev > 0 else SYS
    # erros-padrão "iid" com densidade pelo método fitted e banda de
    # Hall-Sheather ([R] qreg, Methods and formulas; VERIFICAR)
    alpha = 0.05
    zq = st.norm.ppf(q)
    h = n ** (-1 / 3) * st.norm.ppf(1 - alpha / 2) ** (2 / 3) * \
        ((1.5 * st.norm.pdf(zq) ** 2) / (2 * zq * zq + 1)) ** (1 / 3)
    lo_q, hi_q = max(q - h, 1e-4), min(q + h, 1 - 1e-4)
    xbar = (ww @ X) / ww.sum()
    s_hat = float(xbar @ (_qfit(X, y, hi_q, w) - _qfit(X, y, lo_q, w))) / (hi_q - lo_q)
    XtX = (X * ww[:, None]).T @ X
    V = q * (1 - q) * s_hat ** 2 * np.linalg.inv(XtX)
    df_r = n - k
    names = [c.name for c in smp.cols] + (["_cons"] if smp.constant else [])
    from .results import rows_from_columns
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([len(smp.cols)] if smp.constant else [])
    kk = len(names)
    bf, Vf = np.zeros(kk), np.zeros((kk, kk))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = V
    bmap = {c.name: bf[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(Vf[i, i], 0))) for i, c in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    est = Estimates("qreg", smp.depname, names, bf, Vf, n, "t", df_r)
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = ""
    est.extra = {"q": q, "raw": raw_dev, "min": min_dev, "yq": yq, "r2_p": pr2, "N": n}
    est.display = _display_qreg
    from .postest import xb_all

    def _pred(ss, ee, stat, mask):
        xb, ok = xb_all(ss, ee)
        if stat in ("", "xb"):
            return np.where(ok, xb, SYS), ("(option xb assumed; fitted values)" if not stat else "")
        if stat in ("residuals", "r"):
            yy = ss.data.get(ee.depvar).data
            return np.where(ok & (yy < SYS), yy - xb, SYS), ""
        raise StataError(198, f"option {stat} not allowed")
    est.predict = _pred
    # e() na ordem do ereturn list do Stata 14 (compat 0804)
    scalars = [("rank", float(k)), ("sparsity", s_hat), ("bwidth", h),
               ("df_m", float(k - (1 if smp.constant else 0))), ("df_r", float(df_r)),
               ("f_r", 1 / s_hat if s_hat else SYS), ("N", float(n)), ("sum_w", float(ww.sum())),
               ("q_v", yq), ("q", q), ("sum_rdev", raw_dev), ("sum_adev", min_dev), ("convcode", 0.0)]
    macros = [("cmd", "qreg"), ("cmdline", f"qreg {args.strip()}"), ("predict", "qreg_p"),
              ("properties", "b V"), ("marginsnotok", "stdp stddp residuals"), ("vce", "iid"),
              ("denmethod", "fitted"), ("bwmethod", "hsheather"), ("depvar", smp.depname)]
    post(s, est, scalars, macros, smp.mask)
    _display_qreg(s, est)


def _display_qreg(s: "Session", est, *, header: bool = True, table: bool = True, level=None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        q = x["q"]
        title = ("Median regression" if abs(q - 0.5) < 1e-12
                 else f"{format_value(q, '%9.0g', pad=False).strip()} Quantile regression")
        out.write("\n", "text")
        out.write(f"{title:<52}Number of obs = ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        out.write(f"  Raw sum of deviations {_g(x['raw'], '%9.0g'):>8} (about {_g(x['yq'], '%10.0g').strip()})\n", "text")
        out.write(f"  Min sum of deviations {_g(x['min'], '%9.0g'):>8}{'':<20}Pseudo R2     = ", "text")
        out.write(f"{x['r2_p']:>10.4f}\n", "result")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level)
