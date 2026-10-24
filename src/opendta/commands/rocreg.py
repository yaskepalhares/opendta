"""ROC regression ([R] rocreg): nonparametric ROC estimation with covariate-adjusted
control distributions and summaries (AUC, ROC(t), inverse ROC, partial AUC), and the
parametric (probit GLM) ROC model with ROC covariates; standard errors and intervals by a
bootstrap stratified on the case/control status. The statistics are in
opendta.estimation.rocreg; layout details are marked VERIFY and checked in
compat/do/1601_epidemiology.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..estimation import rocreg as RR
from ..lang.syntax import match_options, parse_standard
from ..lang.words import parse_numlist
from ..lazy import lazy
from ._util import touse
from .registry import command
from .summarize import _level

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_OPTS = {"probit": 2, "ml": 2, "roccov": 4, "fprpts": 4, "ctrlfprall": 5, "ctrlcov": 5,
         "ctrlmodel": 5, "pvc": 3, "tiecorrected": 3, "auc": 3, "roc": 3, "invroc": 4, "pauc": 4,
         "nobootstrap": 6, "bseed": 5, "breps": 5, "nobstrata": 6, "nodots": 4, "level": 1,
         "cluster": 2}


@command("rocreg", "rocreg")
def cmd_rocreg(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.weight is not None:
        raise StataError(101, "weights not allowed")
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) < 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, _OPTS)
    level = _level(s, o)
    probit = bool(o.get("probit"))
    refname, classes = names[0], names[1:]
    ctrlcov = s.expand_varlist(str(o["ctrlcov"])) if isinstance(o.get("ctrlcov"), str) else []
    roccov = s.expand_varlist(str(o["roccov"])) if isinstance(o.get("roccov"), str) else []
    if roccov and not probit:
        raise StataError(198, "option roccov() requires probit")   # VERIFY
    if probit and len(classes) > 1:
        raise StataError(103, "too many variables specified")   # VERIFY: one classifier per model
    if probit and any(o.get(k) for k in ("auc", "roc", "invroc", "pauc")):
        raise StataError(198, "summary options auc, roc(), invroc() and pauc() are not allowed with probit")   # VERIFY
    model = str(o.get("ctrlmodel") or "strata").strip()
    if model not in ("strata", "linear"):
        raise StataError(198, "ctrlmodel() must be strata or linear")
    pvc = str(o.get("pvc") or "empirical").strip()
    if pvc not in ("empirical", "normal"):
        raise StataError(198, "pvc() must be empirical or normal")
    ds = s.data
    clustvar = s.expand_varlist(str(o["cluster"]))[0] if isinstance(o.get("cluster"), str) else ""
    cols = {n: ds.get(n) for n in [refname] + classes + ctrlcov + roccov + ([clustvar] if clustvar else [])}
    mask = touse(s, p)
    for n, v in cols.items():
        if v.is_string:
            raise StataError(109, "type mismatch")
        mask &= v.data < SYS
    ref = cols[refname].data[mask] != 0
    if ref.all() or not ref.any():
        raise StataError(2000, "no observations")   # VERIFY: needs cases and controls
    clust = cols[clustvar].data[mask] if clustvar else None
    if o.get("ml"):
        if not probit:
            raise StataError(198, "option ml requires probit")                  # VERIFY
        _rocreg_ml(s, o, args, refname, classes, ctrlcov, roccov, cols, mask, ref, clust, clustvar, level)
        return
    spec = _summary_spec(o) if not probit else []

    def statistic(idx: np.ndarray) -> np.ndarray:
        """The estimates for the observations idx (positions within the sample)."""
        d = ref[idx]
        out = []
        for c in classes:
            y = cols[c].data[mask][idx]
            pv = _pv(y[d], y[~d], [cols[v].data[mask][idx] for v in ctrlcov], d, model, pvc,
                     bool(o.get("tiecorrected")))
            ok = ~np.isnan(pv)
            if probit:
                Z = np.column_stack([cols[v].data[mask][idx][d][ok] for v in roccov]) if roccov \
                    else np.zeros((int(ok.sum()), 0))
                t = np.unique(1 - pv[ok]) if o.get("ctrlfprall") else RR.fpr_points(_fprpts(o))
                t = t[(t > 0) & (t < 1)]
                g = RR.probit_glm(pv[ok], Z, t)
                # i_cons, the ROC covariates, s_cons; and the AUC without covariates
                kz = Z.shape[1]
                b = np.concatenate([[g[kz]], g[:kz], [g[kz + 1]]])
                if not kz:
                    b = np.append(b, st.norm.cdf(b[0] / math.sqrt(1 + b[-1] ** 2)))
                out.append(b)
            else:
                out.append(RR.summaries(pv[ok], spec))
        return np.concatenate(out)

    n = int(mask.sum())
    observed = statistic(np.arange(n))
    colnames, eqs = [], []
    for c in classes:
        if probit:
            # VERIFY: names of e(b) after probit
            colnames += ["_cons"] + roccov + ["_cons"] + ([] if roccov else ["auc"])
            eqs += [c] * (1 + len(roccov)) + ["probit"] + ([] if roccov else [c])
        else:
            labels = [_colname(k, j, spec) for j, (k, _) in enumerate(spec)]
            colnames += labels
            eqs += [c] * len(labels)
    k = len(observed)
    boot = None
    if not o.get("nobootstrap"):
        boot = _bootstrap(s, o, statistic, ref, observed, level, clust)
        V = np.cov(boot["reps"], rowvar=False, ddof=1).reshape(k, k)
    else:
        V = np.zeros((k, k))
    _display(s, o, refname, classes, spec, colnames, observed, boot, level, n, probit, ctrlcov, roccov,
             model, pvc, mask)
    _post(s, o, args, refname, classes, ctrlcov, roccov, colnames, eqs, observed, V, boot, n, mask, probit,
          model, pvc, spec)


def _fprpts(o) -> int:
    k = int(o["fprpts"]) if isinstance(o.get("fprpts"), str) else 10
    if k < 2:
        raise StataError(198, "fprpts() must be at least 2")   # VERIFY
    return k


def _summary_spec(o) -> list[tuple[str, float | None]]:
    """The summaries in Stata's order: ROC(t), inverse ROC, partial AUC, then the AUC
    (compat 1601)."""
    spec: list[tuple[str, float | None]] = []
    for kind in ("roc", "invroc", "pauc"):
        if isinstance(o.get(kind), str):
            for x in parse_numlist(str(o[kind])):
                if not 0 < x < 1 and not (kind == "pauc" and x == 1):
                    raise StataError(125, f"{kind}() values must be between 0 and 1")   # VERIFY
                spec.append((kind, float(x)))
    if o.get("auc") or not spec:
        spec.append(("auc", None))
    return spec


def _colname(kind: str, j: int, spec) -> str:
    if kind == "auc":
        return "auc"
    k = sum(1 for q, _ in spec[:j + 1] if q == kind)
    return f"{kind}_{k}"


def _pv(case_y, ctrl_y, ctrl_cov, d, model, pvc, tiecorrected) -> np.ndarray:
    """Percentile values of the cases, adjusted for the control covariates."""
    if not ctrl_cov:
        return RR.pv_normal(case_y, ctrl_y) if pvc == "normal" else RR.pv_empirical(case_y, ctrl_y, tiecorrected)
    if model == "linear":
        X = np.column_stack(ctrl_cov)
        return RR.pv_linear(case_y, X[d], ctrl_y, X[~d], pvc, tiecorrected)
    # strata of the control covariates: each case against the controls of its stratum;
    # cases without controls in their stratum are left out (VERIFY)
    key = np.column_stack(ctrl_cov)
    kc, kd = key[~d], key[d]
    pv = np.full(len(case_y), np.nan)
    for row in np.unique(kd, axis=0):
        sd, sc = (kd == row).all(1), (kc == row).all(1)
        if sc.any():
            pv[sd] = RR.pv_normal(case_y[sd], ctrl_y[sc]) if pvc == "normal" else \
                RR.pv_empirical(case_y[sd], ctrl_y[sc], tiecorrected)
    return pv


def _bootstrap(s: "Session", o, statistic, ref, observed, level, clust=None) -> dict:
    from ..estimation.bootstrap import _draw
    reps = int(o["breps"]) if isinstance(o.get("breps"), str) else 1000
    if reps < 2:
        raise StataError(198, "breps() must be an integer greater than 1")   # VERIFY
    if isinstance(o.get("bseed"), str):
        s.interp._execute_expanded(f"set seed {str(o['bseed']).strip()}")
    out = s.output
    groups = [np.arange(len(ref))] if o.get("nobstrata") else [np.flatnonzero(~ref), np.flatnonzero(ref)]
    if clust is not None:
        # cluster(): whole clusters are drawn, within the case/control strata when no cluster has
        # both (VERIFY: Stata's rule for mixed clusters)
        ids = np.unique(clust)
        members = {c: np.flatnonzero(clust == c) for c in ids}
        mixed = any(ref[m].any() and not ref[m].all() for m in members.values())
        if o.get("nobstrata") or mixed:
            cgroups = [ids]
        else:
            cgroups = [np.array([c for c in ids if not ref[members[c]].any()]),
                       np.array([c for c in ids if ref[members[c]].all()])]
    dots = not o.get("nodots")
    if dots:
        out.write(f"\nBootstrap replications ({reps})\n", "text")
        out.write("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 \n", "text")
    R = np.full((reps, len(observed)), np.nan)
    for r in range(reps):
        if clust is None:
            idx = np.concatenate([g[_draw(len(g))] for g in groups])
        else:
            idx = np.concatenate([members[c] for g in cgroups for c in g[_draw(len(g))]])
        ok = True
        try:
            if ref[idx].all() or not ref[idx].any():
                raise StataError(2000, "no observations")
            R[r] = statistic(idx)
        except (StataError, np.linalg.LinAlgError, ValueError, ZeroDivisionError):
            ok = False
        if np.isnan(R[r]).any():
            ok = False
        if dots:
            out.write("." if ok else "x", "text")
            if (r + 1) % 50 == 0:
                out.write(f"{r + 1:>6}\n", "text")
    if dots and reps % 50:
        out.write("\n", "text")
    good = R[~np.isnan(R).any(axis=1)]
    if len(good) < 2:
        raise StataError(2000, "insufficient observations to compute bootstrap standard errors")   # VERIFY
    se = good.std(axis=0, ddof=1)
    z = float(st.norm.ppf(0.5 + level / 200))
    a = (100 - level) / 2
    prop = np.clip((good < observed).mean(axis=0), 1e-12, 1 - 1e-12)
    z0 = st.norm.ppf(prop)
    bc = [(np.percentile(good[:, j], st.norm.cdf(2 * z0[j] - z) * 100),
           np.percentile(good[:, j], st.norm.cdf(2 * z0[j] + z) * 100)) for j in range(len(observed))]
    return {"reps": good, "n": len(good), "se": se, "bias": good.mean(axis=0) - observed, "z0": z0,
            "b_bs": good.mean(axis=0),
            "normal": (observed - z * se, observed + z * se),
            "percentile": (np.percentile(good, a, axis=0), np.percentile(good, 100 - a, axis=0)),
            "bc": (np.array([x[0] for x in bc]), np.array([x[1] for x in bc]))}


# ---------------------------------------------------------------------------
# display and e() (layout and names from compat 1601)
# ---------------------------------------------------------------------------

_TITLES = {"auc": "Area under the ROC curve", "roc": "ROC curve", "invroc": "False-positive rate",
           "pauc": "Partial area under the ROC curve"}
_HEAD = {"auc": "AUC", "roc": "ROC", "invroc": "invROC", "pauc": "pAUC"}


def _f(x) -> str:
    return "." if x is None or not math.isfinite(x) else _g9(x)


def _g9(x: float) -> str:
    from .summarize import g9
    return g9(float(x)).strip()


def _display(s: "Session", o, refname, classes, spec, colnames, observed, boot, level, n, probit,
             ctrlcov, roccov, model, pvc, mask) -> None:
    out = s.output
    title = "Parametric ROC estimation" if probit else "Nonparametric ROC estimation"
    out.write("\n", "text")
    out.write(f"{(title if boot is None else 'Bootstrap results'):<48}Number of obs     = ", "text")
    out.write(f"{n:>10,}\n", "result")
    if boot is not None:
        out.write(f"{'':<48}Replications      = ", "text")
        out.write(f"{boot['n']:>10,}\n", "result")
        out.write(f"\n{title}\n", "text")
    out.write("\n", "text")
    if ctrlcov:
        out.write(f"Covariate control      : {'stratification' if model == 'strata' else 'linear regression'}\n",
                  "text")
        out.write(f"Control variables      : {' '.join(ctrlcov)}\n", "text")
    tie = ", corrected for ties" if o.get("tiecorrected") else ""
    out.write(f"Control standardization: {pvc}{tie}\n", "text")
    if probit:
        out.write(f"ROC method             : {'parametric':<23}Link: probit\n", "text")
    else:
        out.write("ROC method             : empirical\n", "text")
    if ctrlcov and model == "linear":
        # the regression of each classifier on the control covariates among the controls
        from ..core.dataset import Variable
        ds = s.data
        tmp = "__rocreg_ctrl"
        flag = np.zeros(ds.nobs)
        flag[mask & (ds.get(refname).data == 0)] = 1
        ds.add(Variable(tmp, "byte", flag))
        try:
            for c in classes:
                out.write(f"\nStatus    : {refname}\nClassifier: {c}\nCovariate control adjustment model:\n",
                          "text")
                s.interp._execute_expanded(f"regress {c} {' '.join(ctrlcov)} if {tmp}")
        finally:
            ds.drop_vars([tmp])
    per = len(colnames) // len(classes)
    if probit:
        for ci, c in enumerate(classes):
            base = ci * per
            out.write(f"\n   Status    : {refname}\n   Classifier: {c}\n   ROC Model : \n", "text")
            _table_head(out, c, level)
            rows = [("_cons", base)] + [(v, base + 1 + j) for j, v in enumerate(roccov)]
            for lab, j in rows:
                _row(out, lab, observed[j], boot, j, level, tags=False)
            out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
            out.write(f"{'probit':>12} |\n", "text")
            out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
            j = base + 1 + len(roccov)
            _row(out, "_cons", observed[j], boot, j, level, tags=False)
            out.write("-" * 78 + "\n", "text")
            if not roccov:
                out.write("\n", "text")
                _table_head(out, "AUC", level)
                _row(out, "", observed[j + 1], boot, j + 1, level)
                out.write("-" * 78 + "\n", "text")
        return
    kinds = list(dict.fromkeys(k for k, _ in spec))
    for kind in kinds:
        out.write(f"\n{_TITLES[kind]}\n", "text")
        for ci, c in enumerate(classes):
            out.write(f"\n   Status    : {refname}\n   Classifier: {c}\n", "text")
            _table_head(out, _HEAD[kind], level)
            for j, (k, x) in enumerate(spec):
                if k == kind:
                    _row(out, "" if kind == "auc" else _g9(x), observed[ci * per + j], boot, ci * per + j, level)
            out.write("-" * 78 + "\n", "text")


def _table_head(out, head: str, level: float) -> None:
    out.write("-" * 78 + "\n", "text")
    out.write(f"{'':>12} |    Observed               Bootstrap\n", "text")
    out.write(f"{head[:12]:>12} |       Coef.       Bias    Std. Err.     [{level:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")


def _row(out, label, obs, boot, j, level, tags: bool = True) -> None:
    """Observed value, bootstrap bias and standard error, and the normal (N), percentile (P)
    and bias-corrected (BC) intervals; missing without the bootstrap."""
    if boot is None:
        bias = se = None
        cis = [(None, None)] * 3
    else:
        bias, se = boot["bias"][j], boot["se"][j]
        cis = [(boot[k][0][j], boot[k][1][j]) for k in ("normal", "percentile", "bc")]
    out.write(f"{label:>12} |", "text")
    out.write(f"{_f(obs):>12}{_f(bias):>11}{_f(se):>12}{_f(cis[0][0]):>13}{_f(cis[0][1]):>11}", "result")
    if not tags:
        # VERIFY: the bootstrap rows of the probit model
        out.write("\n", "text")
        return
    out.write("  (N)\n", "text")
    for tag, (lo, hi) in zip(("  (P)", " (BC)"), cis[1:]):
        out.write(f"{'':>12} |{'':>35}", "text")
        out.write(f"{_f(lo):>13}{_f(hi):>11}", "result")
        out.write(tag + "\n", "text")


def _post(s: "Session", o, args, refname, classes, ctrlcov, roccov, colnames, eqs, observed, V, boot, n, mask,
          probit, model, pvc, spec) -> None:
    """e() in Stata 14's ereturn list order (compat 1601; VERIFY the names with ctrlcov(),
    roc(), invroc() and pauc())."""
    from ..estimation.results import Estimates, post
    from .matrix import Matrix
    est = Estimates("rocreg", refname, colnames, np.asarray(observed, float), V, n, stat="z", eqnames=eqs)
    est.title = "Parametric ROC estimation" if probit else "Nonparametric ROC estimation"
    level = _level(s, o)

    def mat(row, rown="y1"):
        m = Matrix(np.asarray(row, float).reshape(1, -1), [rown], list(colnames))
        m.coleq = list(eqs)
        return m

    def ci(pair):
        m = Matrix(np.vstack(pair), ["ll", "ul"], list(colnames))
        m.coleq = list(eqs)
        return m
    stats = {"auc": "auc" if any(k == "auc" for k, _ in spec) else ""}
    for kind in ("roc", "invroc", "pauc"):
        vals = [x for k, x in spec if k == kind]
        stats[kind] = " ".join(_g9(x) for x in vals)
    common = [("pvc", pvc), ("ctrlcov", " ".join(ctrlcov)), ("ctrlmodel", model if ctrlcov else "")]
    if boot is None:
        scalars = [("N", float(n))] if probit else [("level", level), ("N", float(n))]
        macros = [("title", est.title), ("nobootstrap", "nobootstrap"), ("cmdline", f"rocreg {args.strip()}"),
                  ("refvar", refname), ("classvars", " ".join(classes)), ("cmd", "rocreg")] + common
        if probit:
            macros += [("fprpts", str(_fprpts(o))), ("predict", "rocreg_p"), ("probit", "probit"),
                       ("roccov", " ".join(roccov))]
        else:
            macros += [(k, stats[k]) for k in ("auc", "roc", "invroc", "pauc")]
        macros += [("np", str(len(colnames))), ("properties", "b")]
        post(s, est, scalars, macros, mask)
        s.e.pop("V", None)
        return
    scalars = [("N", float(n)), ("rank", float(np.linalg.matrix_rank(V))), ("level", level)]
    exps = " ".join(f"[{e}]_b[{c}]" for e, c in zip(eqs, colnames))
    macros = [("properties", "b V"), ("exp1", exps), ("vce", "bootstrap")]
    macros += [(k, stats[k]) for k in ("auc", "roc", "invroc", "pauc")] if not probit else \
        [("probit", "probit"), ("roccov", " ".join(roccov))]
    if isinstance(o.get("cluster"), str):
        macros.append(("clustvar", str(o["cluster"]).strip()))            # VERIFY: name and place
    macros += [("pvc", pvc), ("ctrlcov", " ".join(ctrlcov)), ("ctrlmodel", model if ctrlcov else ""),
               ("cmdline", f"rocreg {args.strip()}"), ("classvars", " ".join(classes)), ("refvar", refname),
               ("breps", str(boot["n"])), ("cmd", "rocreg"), ("title", est.title)]
    extra = [("ci_bc", ci(boot["bc"])), ("ci_percentile", ci(boot["percentile"])),
             ("ci_normal", ci(boot["normal"])), ("se", mat(boot["se"])), ("z0", mat(boot["z0"])),
             ("bias", mat(boot["bias"])), ("reps", mat(np.full(len(colnames), float(boot["n"])))),
             ("b_bs", mat(boot["b_bs"]))]
    post(s, est, scalars, macros, mask, extra_mats=extra)


# ---------------------------------------------------------------------------
# probit ml: the binormal model by maximum likelihood (VERIFY: layout and e() names, compat 1912)
# ---------------------------------------------------------------------------

def _rocreg_ml(s: "Session", o, args, refname, classes, ctrlcov, roccov, cols, mask, ref, clust, clustvar,
               level) -> None:
    from ..estimation import ml
    from ..estimation.results import CoefRow, Estimates, coef_table, post
    model = str(o.get("ctrlmodel") or "linear").strip()
    if model != "linear":
        raise StataError(198, "ctrlmodel(strata) not allowed with ml")             # VERIFY
    out = s.output
    n = int(mask.sum())
    X = np.column_stack([cols[v].data[mask] for v in roccov]) if roccov else np.zeros((n, 0))
    Z = np.column_stack([cols[v].data[mask] for v in ctrlcov] + [np.ones(n)])
    kx = X.shape[1]
    fits = []
    for c in classes:
        y = cols[c].data[mask].astype(float)
        out.write("\n", "text")
        res, llobs = RR.binormal_ml(y, ref, X, Z)
        ml.print_log(out, res)
        if not res.converged:
            raise StataError(430, "convergence not achieved")
        fits.append((res, llobs))
    # parameters of the ROC model of each classifier: intercept, ROC covariates and slope; the
    # control model (mean and log standard deviation) stays out of e(b) (VERIFY)
    keep = list(range(kx + 2))
    blocks_V = []
    scores = []
    for res, llobs in fits:
        Vf = np.linalg.inv(-res.H)
        blocks_V.append((Vf + Vf.T) / 2)
        if clust is not None:
            scores.append(RR.obs_scores(llobs, res.b))
    k_all = sum(len(r.b) for r, _ in fits)
    if clust is None:
        Vfull = np.zeros((k_all, k_all))
        pos = 0
        for B in blocks_V:
            Vfull[pos:pos + len(B), pos:pos + len(B)] = B
            pos += len(B)
    else:
        # cluster-robust sandwich over the stacked scores of every classifier
        from scipy.linalg import block_diag
        A = block_diag(*blocks_V)
        S = np.hstack(scores)
        ids = np.unique(clust)
        G = len(ids)
        U = np.vstack([S[clust == g].sum(axis=0) for g in ids])
        Vfull = A @ (U.T @ U) @ A * G / (G - 1)
    idx = []
    pos = 0
    for res, _ in fits:
        idx += [pos + j for j in keep]
        pos += len(res.b)
    b = np.concatenate([r.b[keep] for r, _ in fits])
    V = Vfull[np.ix_(idx, idx)]
    names, eqs = [], []
    for c in classes:
        names += ["_cons"] + roccov + ["_cons"]
        eqs += [c] * (1 + kx) + ["probit"]
    ll = float(sum(r.ll for r, _ in fits))
    out.write("\n", "text")
    out.write(f"{'Parametric ROC estimation':<48}Number of obs     = ", "text")
    out.write(f"{n:>10,}\n", "result")
    out.write("\n", "text")
    if ctrlcov:
        out.write("Covariate control      : linear regression\n", "text")
        out.write(f"Control variables      : {' '.join(ctrlcov)}\n", "text")
    out.write("Control standardization: normal\n", "text")
    out.write(f"{'ROC method             : ML':<47}Link: probit\n", "text")
    if clust is not None:
        out.write(f"{'':<47}(Std. Err. adjusted for {len(np.unique(clust)):,} clusters in {clustvar})\n",
                  "text")
    per = kx + 2
    aucs = []
    for ci, c in enumerate(classes):
        base = ci * per
        out.write(f"\n   Status    : {refname}\n   Classifier: {c}\n   ROC Model : \n", "text")
        se = np.sqrt(np.maximum(np.diag(V), 0))
        rows = [CoefRow("eq", c)]
        rows += [CoefRow("coef", lab, b[base + j], se[base + j]) for j, lab in enumerate(["_cons"] + roccov)]
        rows += [CoefRow("sep"), CoefRow("eq", "probit"),
                 CoefRow("coef", "_cons", b[base + kx + 1], se[base + kx + 1])]
        coef_table(s, "", rows, stat="z", df=None, level=level,
                   vcetype="Robust" if clust is not None else "")
        if not roccov:
            a_, b_ = b[base], b[base + 1]
            Vab = V[np.ix_([base, base + 1], [base, base + 1])]
            val = float(st.norm.cdf(a_ / math.sqrt(1 + b_ * b_)))
            gr = np.array([st.norm.pdf(a_ / math.sqrt(1 + b_ * b_)) / math.sqrt(1 + b_ * b_),
                           -st.norm.pdf(a_ / math.sqrt(1 + b_ * b_)) * a_ * b_ / (1 + b_ * b_) ** 1.5])
            aucs.append((val, float(math.sqrt(max(gr @ Vab @ gr, 0)))))
            out.write("\n", "text")
            coef_table(s, "AUC", [CoefRow("coef", "", val, aucs[-1][1])], stat="z", df=None, level=level)
    est = Estimates("rocreg", refname, names, b, V, n, stat="z", eqnames=eqs)
    est.title = "Parametric ROC estimation"
    scalars = [("N", float(n)), ("ll", ll), ("rank", float(np.linalg.matrix_rank(V)))]
    if clust is not None:
        scalars.append(("N_clust", float(len(np.unique(clust)))))
    macros = [("title", est.title), ("cmdline", f"rocreg {args.strip()}"), ("refvar", refname),
              ("classvars", " ".join(classes)), ("cmd", "rocreg"), ("ml", "ml"), ("probit", "probit"),
              ("roccov", " ".join(roccov)), ("ctrlcov", " ".join(ctrlcov)),
              ("ctrlmodel", "linear" if ctrlcov else ""), ("pvc", "normal"),
              ("vce", "cluster" if clust is not None else "oim"), ("clustvar", clustvar),
              ("predict", "rocreg_p"), ("properties", "b V")]
    post(s, est, scalars, macros, mask)
