"""sem with group(): results, display and e() of multiple-group models ([SEM] sem group options).

Each group has its own copy of the model; the parameters of the classes in ginvariant() are
constrained equal across groups (default: measurement coefficients and intercepts). The expanded
e(b) has one column per group and parameter, named after the factor-variable notation
(1bn.group#c.F); the covariance matrix of the expanded vector is T V T' with T the map from the
global to the expanded parameters. The layout of the table and of the header is VERIFY (compat
case 1909)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from .results import Estimates, crit, post
from .sem import Model, _g, _param_names, _standardized, st

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _gname(lev: float) -> str:
    return f"{lev:g}"


def _col_name(model: Model, item, lev: float, gvar: str, first: bool) -> tuple[str, str]:
    """(equation, name) of an item of group `lev` in the expanded e(b)."""
    kind, a, b = item
    prefix = f"{_gname(lev)}{'bn' if first else ''}.{gvar}"
    eq, nm = _param_names(model)[model.free.index(item)]
    if kind == "path":
        return eq, f"{prefix}#c.{nm}"
    return eq, prefix


def finish_groups(s: "Session", args: str, prob, spec, th, obj, V, conv, o, info: dict) -> None:
    models, datas = prob.models, prob.data
    G = len(models)
    gvar, levels = info["gname"], info["levels"]
    meth = info["meth"]
    N = info["N"]
    model = models[0]
    # observed exogenous variables are conditioned on unless noxconditional (VERIFY)
    oxset = {n for n in model.exog if n in spec.observed}
    xcond = bool(oxset) and not o.get("noxconditional") and meth == "ml"
    hidden = set()
    if xcond:
        hidden = {it for it in model.free if it[0] in ("mean", "var", "cov")
                  and it[1] in oxset and it[2] in oxset}
    T = prob.expand_matrix()
    th_exp = T @ th
    V_exp = T @ V @ T.T
    cols = []          # (group, item index, eq, name)
    keep = []
    r = 0
    for g, m in enumerate(models):
        for i, it in enumerate(m.free):
            if it not in hidden:
                eq, nm = _col_name(m, it, levels[g], gvar, g == 0)
                cols.append((g, i, eq, nm))
                keep.append(r)
            r += 1
    est = Estimates("sem", " ".join(model.endog), [c[3] for c in cols], th_exp[keep].copy(),
                    V_exp[np.ix_(keep, keep)], N, "z", None)
    est.eqnames = [c[2] for c in cols]
    est.model = model
    est.level = info["level"]
    ll = info["ll"]
    est.extra = dict(N=N, ll=ll, df=info["df"], chi2_ms=info["chi2_ms"], p_ms=info["p_ms"],
                     chi2_bs=info["chi2_bs"], df_bs=info["df_bs"], p_bs=info["p_bs"],
                     mean=datas[0].mean, cov=datas[0].cov, ll_s=info["ll_s"],
                     standardized=bool(o.get("standardized")), th=th.copy(), V=V, xcond=xcond,
                     method=meth, prob=prob, mask=info["mask"], gvar=gvar, levels=levels,
                     hidden=hidden, T=T, cols=cols, groups=True)
    est.display = display_groups
    est.predict = None
    from .semestat import ESTAT
    est.estat = dict(ESTAT)
    ilog = list(prob.ilog) if prob.ilog else [ll]
    pobs = model.p
    df_s = pobs * (pobs + 3) // 2
    scalars = [("converged", 1.0 if conv else 0.0), ("rc", 0.0), ("ic", float(len(ilog) - 1)),
               ("N_groups", float(G)), ("modelmeans", 1.0), ("N", float(N)), ("N_missing", 0.0),
               ("ll", ll), ("critvalue", ll if meth != "adf" else -obj / max(N, 1)),
               ("df_m", float(prob.nglob)), ("rank", float(np.linalg.matrix_rank(est.V))),
               ("estimates", 1.0), ("critvalue_s", info["ll_s"]), ("critvalue_b", info["ll_b"]),
               ("df_s", float(G * df_s)), ("p_bs", info["p_bs"]), ("df_bs", float(info["df_bs"])),
               ("chi2_bs", info["chi2_bs"]), ("p_ms", info["p_ms"]), ("df_ms", float(info["df"])),
               ("chi2_ms", info["chi2_ms"])]
    if meth == "adf":
        scalars = [(k, v) for k, v in scalars if k not in ("ll", "critvalue_s", "critvalue_b")]
    macros = [("cmdline", f"sem {args.strip()}"), ("groupvar", gvar),
              ("cmd", "sem"), ("predict", "sem_p"), ("estat_cmd", "sem_estat"),
              ("vcetype", "OIM"), ("vce", "oim"), ("title", "Structural equation model"),
              ("method", meth), ("data", "raw"), ("technique", "nr"), ("properties", "b V")]
    mats = []
    post(s, est, scalars, macros, info["mask"], extra_mats=mats)
    display_groups(s, est)


def constraint_lines(est: Estimates) -> list[str]:
    """The linear constraints of the model as Stata lists them: equalities across groups, then
    the loadings fixed at 1 and the other fixed paths (VERIFY order and labels)."""
    x = est.extra
    prob = x["prob"]
    gvar, levels = x["gvar"], x["levels"]
    models = prob.models
    eqs = []
    # equalities: items sharing a global entry
    first: dict[int, tuple[int, int]] = {}
    for g, m in enumerate(models):
        for i, gid in enumerate(prob.gmap[g]):
            if x["hidden"] and m.free[i] in x["hidden"]:
                continue
            if gid in first:
                g0, i0 = first[gid]
                eq0, nm0 = _col_name(models[g0], models[g0].free[i0], levels[g0], gvar, g0 == 0)
                eq1, nm1 = _col_name(m, m.free[i], levels[g], gvar, g == 0)
                eqs.append(f"[{eq0}]{nm0} - [{eq1}]{nm1} = 0")
            else:
                first[gid] = (g, i)
    # Stata's order (case 1909): the fixed paths of the first group, the equalities, the fixed
    # paths of the other groups, then the latent means fixed at 0
    fixed = []
    for g, m in enumerate(models):
        spec = m.spec
        pre = f"{_gname(levels[g])}{'bn' if g == 0 else ''}.{gvar}"
        lines = []
        anchors = sorted(m.anchor.items(), key=lambda la: la[0] not in m.endog)
        for L, a in anchors:
            lines.append(f"[{a}]{pre}#c.{L} = 1")
        for (a, b), v in spec.fixed.items():
            lines.append(f"[{b}]{pre}#c.{a} = {v:g}")
        fixed.append(lines)
    means = []
    for g, m in enumerate(models):
        pre = f"{_gname(levels[g])}{'bn' if g == 0 else ''}.{gvar}"
        for L in m.spec.latent:
            if L in m.exog and ("mean", L, L) not in m.free:
                means.append(f"[mean({L})]{pre} = 0")
    return fixed[0] + eqs + [x for lines in fixed[1:] for x in lines] + means


def display_groups(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                   level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    prob = x["prob"]
    models = prob.models
    spec = models[0].spec
    gvar, levels = x["gvar"], x["levels"]
    lv = level or est.level
    G = len(models)
    if header:
        out.write("\n" + f"{'Structural equation model':<48}{'Number of obs':<18}= ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        out.write(f"{'Grouping variable  = ' + gvar:<48}{'Number of groups':<18}= ", "text")
        out.write(f"{G:>10,}\n", "result")
        out.write(f"Estimation method  = {x['method']}\n", "text")
        if x["method"] != "adf":
            out.write("Log likelihood     = ", "text")
            out.write(f"{_g(x['ll'], '%10.0g'):>10}\n", "result")
        out.write("\n", "text")
        for k, line in enumerate(constraint_lines(est), 1):
            out.write(f" ({k:2d})  {line}\n", "text")
        if not constraint_lines(est):
            pass
    if not table:
        return
    th = x["th"]
    V = x["V"]
    c = crit("z", lv, None)
    std = x.get("standardized")
    T = x["T"]
    th_exp = T @ th
    se_exp = np.sqrt(np.maximum(np.diag(T @ V @ T.T), 0))
    if std:
        def f(t):
            return np.concatenate([_standardized(m, t[prob.gmap[g]]) for g, m in enumerate(models)])
        full = f(th)
        J = np.zeros((len(full), len(th)))
        for i in range(len(th)):
            e = np.zeros(len(th))
            h = 1e-6 * (abs(th[i]) + 1)
            e[i] = h
            J[:, i] = (f(th + e) - f(th - e)) / (2 * h)
        th_exp = full
        se_exp = np.sqrt(np.maximum(np.diag(J @ V @ J.T), 0))
    offs = np.cumsum([0] + [len(m.free) for m in models])

    def val(g, item):
        i = models[g].free.index(item)
        return th_exp[offs[g] + i], se_exp[offs[g] + i]
    W = 12
    dash = "-" * (W + 1) + "+" + "-" * 64 + "\n"
    starred = [False]
    out.write("-" * (W + 66) + "\n", "text")
    out.write(f"{'':>{W}} |                 OIM\n", "text")
    out.write(f"{'Standardized' if std else '':>{W}} |      Coef.   Std. Err.      z    P>|z|     "
              f"[{lv:g}% Conf. Interval]\n", "text")
    out.write(dash, "text")

    def row(label, b, s_, kind="coef"):
        out.write(f"{label:>{W}} |" if label == "[*]" else f"{label:>{W - 1}}  |", "text")
        if kind == "constrained":
            out.write(f"  {_g(b):>9}  (constrained)\n", "result")
            return
        if s_ > 0 and np.isfinite(s_):
            if kind == "var":
                lo, hi = b * math.exp(-c * s_ / b), b * math.exp(c * s_ / b)
                out.write(f"  {_g(b):>9}  {_g(s_):>9}" + " " * 17 + f"    {_g(lo):>9}   {_g(hi):>9}\n",
                          "result")
                return
            z = b / s_
            pv = 2 * float(st.norm.sf(abs(z)))
            out.write(f"  {_g(b):>9}  {_g(s_):>9}  {z:>7.2f}  {pv:>6.3f}    {_g(b - c * s_):>9}"
                      f"   {_g(b + c * s_):>9}\n", "result")
        else:
            out.write(f"  {_g(b):>9}  {'.':>9}\n", "result")

    def same(item):
        """Constrained equal across groups: one global entry, or fixed at one value everywhere."""
        if all(item in m.free for m in models):
            ids = {int(prob.gmap[g][models[g].free.index(item)]) for g in range(G)}
            return len(ids) == 1
        return False

    def group_rows(fn, item=None, fixed_same=False):
        if fixed_same or (item is not None and same(item)):
            starred[0] = True
            fn(0, "[*]")
            return
        for g in range(G):
            fn(g, f"{_gname(levels[g])}")

    m0 = models[0]
    eqs = [n for n in m0.names if n in m0.endog]
    meas = [n for n in eqs if n in spec.observed and all(a in spec.latent for a, b in spec.paths if b == n)]
    struct = [n for n in eqs if n not in meas]
    for section, members in (("Structural", struct), ("Measurement", meas)):
        if not members:
            continue
        out.write(f"{section:<{W}} |\n", "text")
        for j, y in enumerate(members):
            if j:
                out.write("  " + "-" * (W - 1) + "+" + "-" * 64 + "\n", "text")
            out.write(f"  {y + ' <-':<{W - 2}} |\n", "text")
            for a, b in spec.paths:
                if b != y:
                    continue
                out.write(f"{a:>{W}} |\n", "text")

                def one(g, lab, a=a, b=b):
                    m = models[g]
                    it = ("path", a, b)
                    if it in m.free:
                        row(lab, *val(g, it))
                    else:
                        row(lab, spec.fixed.get((a, b), 1.0), 0, "constrained")
                it = ("path", a, b)
                fixed_all = all(it not in m.free for m in models)
                group_rows(one, it, fixed_same=fixed_all)
            if y in spec.observed:
                out.write(f"{'_cons':>{W}} |\n", "text")
                group_rows(lambda g, lab, y=y: row(lab, *val(g, ("mean", y, y))), ("mean", y, y))
        out.write(dash, "text")
    means = [n for n in m0.exog if n in spec.observed and ("mean", n, n) not in x["hidden"]]
    lat_means = [n for n in spec.latent if any(("mean", n, n) in m.free for m in models)]
    for n in means + lat_means:
        out.write(f"{'mean(' + n + ')':>{W + 1}}|\n", "text")

        def one(g, lab, n=n):
            it = ("mean", n, n)
            if it in models[g].free:
                row(lab, *val(g, it))
            else:
                row(lab, 0.0, 0, "constrained")
        group_rows(one, ("mean", n, n))
    if means or lat_means:
        out.write(dash, "text")
    var_order = ([n for n in m0.names if n in m0.endog and n in spec.observed]
                 + [n for n in m0.names if n in m0.endog and n not in spec.observed]
                 + [n for n in m0.names if n not in m0.endog and n not in spec.observed]
                 + [n for n in m0.names if n not in m0.endog and n in spec.observed])
    for n in var_order:
        it = ("var", n, n)
        if it in x["hidden"]:
            continue
        lab = f"var(e.{n})" if n in m0.endog else f"var({n})"
        out.write(f"{lab:>{W + 1}}|\n", "text")

        def one(g, lg, it=it, n=n):
            m = models[g]
            if it in m.free:
                row(lg, *val(g, it), "var")
            elif n in m.var_fixed:
                row(lg, m.var_fixed[n], 0, "constrained")
        group_rows(one, it)
    covs = [it for it in m0.free if it[0] == "cov" and it not in x["hidden"]]
    if covs:
        out.write(dash, "text")
    for it in covs:
        _, a, b = it
        ea = f"e.{a}" if a in m0.endog else a
        eb = f"e.{b}" if b in m0.endog else b
        out.write(f"{'cov(' + ea + ',' + eb + ')':>{W + 1}}|\n", "text")
        group_rows(lambda g, lab, it=it: row(lab, *val(g, it)), it)
    out.write("-" * (W + 66) + "\n", "text")
    if starred[0]:
        out.write("Note: [*] identifies parameter estimates constrained to be equal across\n      groups.\n",
                  "text")
    df = x["df"]
    pv = "     ." if df == 0 else f"{x['p_ms']:>6.4f}"
    lab = f"chi2({df})"
    out.write(f"LR test of model vs. saturated: {lab:<10}= {x['chi2_ms']:>9.2f}, Prob > chi2 = {pv}\n",
              "text")
