"""estat scoretests, estat mindices and estat ginvariant after sem ([SEM] estat scoretests,
estat mindices, estat ginvariant): Lagrange multiplier tests of omitted paths and of the
constraints across groups."""

from __future__ import annotations

import copy
import math

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lazy import lazy
from .sem import Model, _dmoments
from .semestat import MINCHI2, _fv, _nobs, _opts

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

SYS = M.SYSMISS


def _extended(prob, split=frozenset(), split_all=frozenset()):
    from .semfit import Problem
    return Problem(prob.models, prob.data, prob.method, prob.inv, split=set(split), split_all=set(split_all))


def _theta_in(prob, th_glob, ext):
    t = np.zeros(ext.nglob)
    for g in range(len(prob.models)):
        t[ext.gmap[g]] = th_glob[prob.gmap[g]]
    return t


def _info(ext, t):
    """Information of the extended problem: expected information for ml (as Stata's score tests
    and modification indices, case 1909), minus the Hessian otherwise."""
    if ext.method != "ml":
        return -ext.hess(t)
    I = np.zeros((ext.nglob, ext.nglob))
    for g, (m, d) in enumerate(zip(ext.models, ext.data)):
        ids = ext.gmap[g]
        Ig = _eim(m, t[ids], d.N)
        np.add.at(I, (ids[:, None], ids[None, :]), Ig)
    return I


def _lm(prob, th_glob, ext, t=None) -> tuple[float, int]:
    t = _theta_in(prob, th_glob, ext) if t is None else t
    sc = ext.grad(t)
    info = _info(ext, t)
    stat = float(sc @ np.linalg.pinv(info) @ sc)
    return max(stat, 0.0), ext.nglob - prob.nglob


def _th_glob(est):
    x = est.extra
    return x.get("th_glob", x["th"])


def _label(model: Model, it) -> str:
    kind, a, b = it
    if kind == "path":
        return f"{b} <- {a}"
    if kind == "mean":
        return f"mean({a})" if a in model.exog else f"{a} <- _cons"
    if kind == "var":
        return f"var(e.{a})" if a in model.endog else f"var({a})"
    ea = f"e.{a}" if a in model.endog else a
    eb = f"e.{b}" if b in model.endog else b
    return f"cov({ea},{eb})"


def estat_scoretests(s, est, args):
    from ..commands.matrix import Matrix
    _, o = _opts(args, {"minchi2": 3})
    minchi2 = float(o["minchi2"]) if o.get("minchi2") else MINCHI2
    prob = est.extra["prob"]
    th = _th_glob(est)
    out = s.output
    cons = []      # (text, ext problem)
    for g, m in enumerate(prob.models):
        groups: dict[str, list] = {}
        for it, lab in m.ties.items():
            groups.setdefault(lab, []).append(it)
        for lab, items in groups.items():
            for it in items[1:]:
                text = f"[{items[0][2]}]{items[0][1]} - [{it[2]}]{it[1]} = 0"
                cons.append((text, _extended(prob, split={(g, it)})))
    rows = []
    for k, (text, ext) in enumerate(cons, 1):
        stat, df = _lm(prob, th, ext)
        rows.append((k, stat, max(df, 1)))
    shown = [r for r in rows if r[1] >= minchi2]
    if not shown:
        from ..core.formats import format_value
        # Stata's own value of invchi2tail(1, .05) (case 1909)
        lim = "3.841458820694123" if not o.get("minchi2") else format_value(minchi2, "%18.0g").strip()
        out.write(f"(no score tests to report; all chi2 values less than {lim})\n", "text")
    else:
        # VERIFY: the layout when some test is shown
        out.write("\nScore tests for linear constraints\n\n", "text")
        for k, (text, _ext) in enumerate(cons, 1):
            out.write(f" ({k:2d}) {text}\n", "text")
        out.write(f"\n{'':>8}{'chi2':>10}{'df':>6}{'P>chi2':>9}\n", "text")
        for k, stat, df in shown:
            out.write(f" ({k:2d}) ", "text")
            out.write(f"{stat:>9.3f}{df:>6}{float(st.chi2.sf(stat, df)):>9.2f}\n", "result")
    s.r = {"N_groups": float(len(prob.data)), "nobs": _nobs(s, est)}
    if rows:
        s.r["Cns_sctest"] = Matrix(np.array([[r[1], r[2], float(st.chi2.sf(r[1], r[2]))] for r in rows]),
                                   [f"c{r[0]}" for r in rows], ["chi2", "df", "p"])


def _augment(model: Model, item):
    kind, a, b = item
    spec = copy.deepcopy(model.spec)
    covs = list(model.covs_in)
    if kind == "path":
        spec.paths.append((a, b))
    else:
        covs.append((a, b))
    return Model(spec, covs, dict(model.var_fixed))


def _eim(model: Model, th: np.ndarray, N: int) -> np.ndarray:
    """Expected information of the normal likelihood: N (dmu' Si dmu + tr(Si dS Si dS) / 2)."""
    mu, Sig, dm = _dmoments(model, th)
    Si = np.linalg.inv(Sig)
    k = len(dm)
    I = np.zeros((k, k))
    for a in range(k):
        ma, Sa = dm[a]
        for b in range(a + 1):
            mb, Sb = dm[b]
            v = 0.0
            if ma is not None and mb is not None:
                v += float(ma @ Si @ mb)
            if Sa is not None and Sb is not None:
                v += 0.5 * float(np.sum((Si @ Sa) * (Si @ Sb).T))
            I[a, b] = I[b, a] = N * v
    return I


def estat_mindices(s, est, args):
    from ..commands.matrix import Matrix
    from .semfit import Problem, item_class, parse_ginvariant
    _, o = _opts(args, {"showpclass": 5, "minchi2": 3})
    minchi2 = float(o["minchi2"]) if o.get("minchi2") else MINCHI2
    want = parse_ginvariant("all" if str(o.get("showpclass") or "all").strip() == "all"
                            else str(o["showpclass"]))
    if est.extra.get("groups"):
        raise StataError(198, "estat mindices after sem with group() is not available yet in OpenDTA")
    prob = est.extra["prob"]
    model = est.model
    spec = model.spec
    th = est.extra["th"]
    data = prob.data[0]
    mu, Sig, (A, S, m, B, Sv, mv) = model.implied(th)
    sd = np.sqrt(np.maximum(np.diag(Sv), 1e-300))
    cands = []
    exo_obs = {n for n in model.exog if n in spec.observed}
    for to in model.names:
        if to not in model.endog:
            continue
        for fr in model.names:
            # VERIFY (case 1909): Stata does not consider paths out of observed endogenous variables
            if fr == to or (fr, to) in spec.paths or (fr in spec.observed and fr in model.endog):
                continue
            cands.append(("path", fr, to))
    err = [n for n in model.names if n in model.endog]
    have = {tuple(sorted(c, key=lambda n: model.idx[n])) for c in model.covs_in}
    for i, a in enumerate(err):
        for b in err[i + 1:]:
            key = tuple(sorted((a, b), key=lambda n: model.idx[n]))
            if key in have:
                continue
            cands.append(("cov", key[0], key[1]))
    rows = []
    for it in cands:
        try:
            aug = _augment(model, it)
        except StataError:
            continue
        cls = item_class(aug, it)
        if cls not in want:
            continue
        th_aug = np.array([th[model.free.index(x)] if x in model.free else 0.0 for x in aug.free])
        sub = Problem([aug], [data], prob.method, set())
        # the observed exogenous variables stay unconstrained as in the fitted model
        gr = sub.grad(th_aug)
        info = _eim(aug, th_aug, data.N) if prob.method == "ml" else -sub.hess(th_aug)
        p = aug.free.index(it)
        c = [i for i in range(len(th_aug)) if i != p]
        ieff = info[p, p] - info[p, c] @ np.linalg.pinv(info[np.ix_(c, c)]) @ info[c, p]
        if ieff <= 1e-12:
            continue
        mi = gr[p] ** 2 / ieff
        epc = gr[p] / ieff
        if it[0] == "path":
            sepc = epc * sd[model.idx[it[1]]] / sd[model.idx[it[2]]]
        else:
            va, vb = S[model.idx[it[1]], model.idx[it[1]]], S[model.idx[it[2]], model.idx[it[2]]]
            sepc = epc / math.sqrt(va * vb) if va > 0 and vb > 0 else SYS
        rows.append((it, cls, mi, epc, sepc))
    order = {n: i for i, n in enumerate(model.names)}
    shown = sorted([r for r in rows if r[2] >= minchi2],
                   key=lambda r: (r[0][0] != "path", order[r[0][2]], order[r[0][1]]))
    out = s.output
    out.write("\nModification indices\n\n", "text")
    if not shown:
        out.write(f"(no modification indices to report; all chi2 values less than {minchi2})\n", "text")
    else:
        line = "-" * 14 + "+" + "-" * 46
        out.write("-" * 61 + "\n", "text")
        out.write(f"{'':14}|{'Standard':>46}\n", "text")
        out.write(f"{'':14}|{'MI':>10}{'df':>7}{'P>MI':>7}{'EPC':>11}{'EPC':>11}\n", "text")
        out.write(line + "\n", "text")

        def vals(mi, epc, sepc):
            return (f"{mi:>10.3f}{1:>7}{float(st.chi2.sf(mi, 1)):>7.2f}{_fv(epc, '%9.0g'):>11}"
                    f"{_fv(sepc, '%9.0g'):>11}\n")
        paths = [x for x in shown if x[0][0] == "path"]
        covs = [x for x in shown if x[0][0] != "path"]
        # VERIFY: the Structural section (paths into latent variables) is not in case 1909
        for sec, members in (("Structural", [x for x in paths if x[0][2] in spec.latent]),
                             ("Measurement", [x for x in paths if x[0][2] in spec.observed])):
            if not members:
                continue
            out.write(f"{sec:<14}|\n", "text")
            cur = None
            for it, cls, mi, epc, sepc in members:
                if it[2] != cur:
                    cur = it[2]
                    out.write(f"  {cur + ' <-':<12}|\n", "text")
                out.write(f"{it[1][-12:]:>13} |", "text")
                out.write(vals(mi, epc, sepc), "result")
            out.write(line + "\n", "text")
        for it, cls, mi, epc, sepc in covs:
            out.write(f"{_label(model, it):<14}|", "text")
            out.write(vals(mi, epc, sepc), "result")
        out.write("-" * 61 + "\n", "text")
        out.write("EPC = expected parameter change\n", "text")
    classes = ["scoef", "scons", "mcoef", "mcons", "serrvar", "merrvar", "smerrcov", "meanex", "covex"]
    r = {"N_groups": 1.0, "nobs": _nobs(s, est)}
    if shown:
        r["mindices_pclass"] = Matrix(np.array([[classes.index(x[1]) + 1] for x in shown], dtype=float),
                                      [_label(model, x[0]) for x in shown], ["pclass"])
        r["mindices"] = Matrix(np.array([[x[2], 1, float(st.chi2.sf(x[2], 1)), x[3], x[4]] for x in shown]),
                               [_label(model, x[0]) for x in shown], ["MI", "df", "P>MI", "EPC", "stdEPC"])
    s.r = r


def estat_ginvariant(s, est, args):
    from ..commands.matrix import Matrix
    from .semfit import PCLASSES, item_class, parse_ginvariant
    if not est.extra.get("groups"):
        raise StataError(322, "estat ginvariant is only available after sem with group()")   # VERIFY
    _, o = _opts(args, {"showpclass": 5, "class": 3, "legend": 3})
    show = parse_ginvariant("all" if str(o.get("showpclass") or "all").strip() == "all"
                            else str(o["showpclass"]))
    prob = est.extra["prob"]
    th = _th_glob(est)
    Vg = est.extra["V"]
    G = len(prob.models)
    m0 = prob.models[0]
    out = s.output
    # items present in every group, by class
    items = [it for it in m0.free if all(it in m.free for m in prob.models)
             and not (it in est.extra["hidden"])]
    # loadings fixed for identification (same value in every group) are tested by freeing them
    # in the groups after the first (Stata reports a score test for them, case 1909)
    anchors = [("path", m0.names[j], m0.names[i]) for (i, j), v in m0.fixed_A.items()
               if all(m.fixed_A.get((i, j)) == v for m in prob.models)]
    res = []
    for it in anchors:
        stat, df = _lm_anchor(prob, th, it)
        cls = item_class(m0, it)
        if cls in show:
            res.append((it, cls, (np.nan, 0, np.nan), (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)))
    for it in items:
        i0 = m0.free.index(it)
        cls = prob.classes[0][i0]
        if cls not in show:
            continue
        ids = [prob.gmap[g][prob.models[g].free.index(it)] for g in range(G)]
        wald = (np.nan, 0, np.nan)
        score = (np.nan, 0, np.nan)
        if len(set(ids)) == 1:
            ext = _extended(prob, split_all={it})
            stat, df = _lm(prob, th, ext)
            score = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        else:
            R = np.zeros((G - 1, prob.nglob))
            for g in range(1, G):
                R[g - 1, ids[g]] += 1
                R[g - 1, ids[0]] -= 1
            d = R @ th
            Vd = R @ Vg @ R.T
            stat = float(d @ np.linalg.pinv(Vd) @ d)
            df = int(np.linalg.matrix_rank(Vd))
            wald = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        res.append((it, cls, wald, score))
    out.write("\nTests for group invariance of parameters\n\n", "text")
    if o.get("class"):
        _class_tests(s, est, prob, th, Vg, show)
    else:
        res = _ginv_order(m0, res)
        out.write("-" * 78 + "\n", "text")
        out.write(f"{'':12} |{'':12}Wald Test{'':23}Score Test\n", "text")
        out.write(f"{'':12} |{'chi2':>10}{'df':>11}{'p>chi2':>10}{'chi2':>11}{'df':>12}{'p>chi2':>10}\n", "text")
        out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
        eqs = [x for x in res if x[0][0] in ("path", "mean")]
        rest = [x for x in res if x[0][0] not in ("path", "mean")]
        spec = m0.spec
        for sec, pool in (("Structural", spec.latent), ("Measurement", spec.observed)):
            # VERIFY: the Structural section is not in case 1909
            members = [x for x in eqs if x[0][2 if x[0][0] == "path" else 1] in pool]
            if not members:
                continue
            out.write(f"{sec:<12} |\n", "text")
            cur = None
            for it, cls, wald, score in members:
                y = it[2] if it[0] == "path" else it[1]
                if y != cur:
                    if cur is not None:
                        out.write("  " + "-" * 11 + "+" + "-" * 64 + "\n", "text")
                    cur = y
                    out.write(f"  {y + ' <-':<10} |\n", "text")
                lab = it[1] if it[0] == "path" else "_cons"
                out.write(f"{lab[-12:]:>12} |", "text")
                out.write(_tcell(wald) + _tcell(score, wide=True) + "\n", "result")
        if rest:
            out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
            for it, cls, wald, score in rest:
                out.write(f"{_label(m0, it)[-13:]:>13}|", "text")
                out.write(_tcell(wald) + _tcell(score, wide=True) + "\n", "result")
        out.write("-" * 78 + "\n", "text")
    r = {"N_groups": float(G), "nobs": _nobs(s, est)}
    if res:
        r["test"] = Matrix(np.array([[x[2][0], x[2][1], x[2][2], x[3][0], x[3][1], x[3][2]] for x in res],
                                    dtype=float), [_label(m0, x[0]) for x in res],
                           ["Wald_chi2", "Wald_df", "Wald_p", "score_chi2", "score_df", "score_p"])
        r["test_pclass"] = Matrix(np.array([[PCLASSES.index(x[1]) + 1] for x in res], dtype=float),
                                  [_label(m0, x[0]) for x in res], ["pclass"])
    s.r = r


def _tcell(t, wide=False) -> str:
    """Wald cells are 10, 11 and 10 wide; score cells 11, 12 and 10 (case 1909)."""
    stat, df, p = t
    w1, w2 = (11, 12) if wide else (10, 11)
    if df == 0 or not np.isfinite(stat):
        return f"{'.':>{w1}}{'.':>{w2}}{'.':>10}"
    return f"{stat:>{w1}.3f}{df:>{w2}}{p:>10.4f}"


def _ginv_order(m0, res):
    """Equations in the order of the variables, each with its paths and then its intercept;
    variances and covariances after them."""
    pos = {n: i for i, n in enumerate(m0.names)}

    def key(x):
        kind, a, b = x[0]
        if kind == "path":
            return (0, pos[b], 0, pos[a])
        if kind == "mean":
            return (0, pos[a], 1, 0)
        return (1, 0 if kind == "var" else 1, pos[a], pos[b])
    return sorted(res, key=key)


def _lm_anchor(prob, th_glob, it):
    """Score test of a path fixed in every group: the path becomes free (at its fixed value) in
    the groups after the first, which leaves the model identified."""
    from .semfit import Problem
    kind, a, b = it
    models = []
    for g, m in enumerate(prob.models):
        if g == 0:
            models.append(m)
            continue
        mm = copy.copy(m)
        mm.fixed_A = dict(m.fixed_A)
        key = (m.idx[b], m.idx[a])
        v = mm.fixed_A.pop(key)
        mm.free = list(m.free) + [it]
        mm.nfree = len(mm.free)
        mm._anchor_value = v
        models.append(mm)
    ext = Problem(models, prob.data, prob.method, prob.inv,
                  split={(g, it) for g in range(1, len(models))})
    t = np.zeros(ext.nglob)
    for g, m in enumerate(models):
        local = prob.group_theta(th_glob, g)
        vals = np.concatenate([local, [m._anchor_value]]) if g else local
        t[ext.gmap[g]] = vals
    return _lm(prob, th_glob, ext, t=t)


def _class_tests(s, est, prob, th, Vg, show):
    from .semfit import PCLASSES
    out = s.output
    G = len(prob.models)
    m0 = prob.models[0]
    out.write("-" * 76 + "\n", "text")
    out.write(f"{'Class':>14} |{'Wald Test':^30}|{'Score Test':^30}\n", "text")
    out.write(f"{'':>14} |{'chi2':>10}{'df':>6}{'p':>9}   |{'chi2':>10}{'df':>6}{'p':>9}\n", "text")
    out.write("-" * 15 + "+" + "-" * 60 + "\n", "text")
    for cls in PCLASSES:
        if cls not in show:
            continue
        its = [it for i, it in enumerate(m0.free) if prob.classes[0][i] == cls
               and all(it in m.free for m in prob.models) and it not in est.extra["hidden"]]
        if not its:
            continue
        constrained = cls in prob.inv
        wald = score = (np.nan, 0, np.nan)
        if constrained:
            ext = _extended(prob, split_all=set(its))
            stat, df = _lm(prob, th, ext)
            score = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        else:
            rows = []
            for it in its:
                ids = [prob.gmap[g][prob.models[g].free.index(it)] for g in range(G)]
                for g in range(1, G):
                    r = np.zeros(prob.nglob)
                    r[ids[g]] += 1
                    r[ids[0]] -= 1
                    rows.append(r)
            R = np.array(rows)
            d = R @ th
            Vd = R @ Vg @ R.T
            stat = float(d @ np.linalg.pinv(Vd) @ d)
            df = int(np.linalg.matrix_rank(Vd))
            wald = (stat, df, float(st.chi2.sf(stat, df)) if df > 0 else np.nan)
        out.write(f"{cls:>14} |", "text")
        out.write(_tcell(wald).replace("", "") + "   |" + _tcell(score) + "\n", "result")
    out.write("-" * 76 + "\n", "text")
