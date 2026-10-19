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
    for k in ("ml", "cluster"):
        if o.get(k):
            raise StataError(198, f"option {k} of rocreg is not available yet in OpenDTA")
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
    cols = {n: ds.get(n) for n in [refname] + classes + ctrlcov + roccov}
    mask = touse(s, p)
    for n, v in cols.items():
        if v.is_string:
            raise StataError(109, "type mismatch")
        mask &= v.data < SYS
    ref = cols[refname].data[mask] != 0
    if ref.all() or not ref.any():
        raise StataError(2000, "no observations")   # VERIFY: needs cases and controls
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
                out.append(RR.probit_glm(pv[ok], Z, t))
            else:
                out.append(RR.summaries(pv[ok], spec))
        return np.concatenate(out)

    n = int(mask.sum())
    observed = statistic(np.arange(n))
    colnames, eqs = [], []
    for c in classes:
        labels = roccov + ["i_cons", "s_cons"] if probit else [_colname(k, j, spec) for j, k in enumerate(spec)]
        colnames += labels
        eqs += [c] * len(labels)
    k = len(observed)
    boot = None
    if not o.get("nobootstrap"):
        boot = _bootstrap(s, o, statistic, ref, observed, level)
        V = np.cov(boot["reps"], rowvar=False, ddof=1).reshape(k, k)
    else:
        V = np.zeros((k, k))
    _display(s, o, refname, classes, spec, colnames, observed, boot, level, n, probit, ctrlcov, roccov,
             model, pvc)
    _post(s, o, refname, classes, ctrlcov, roccov, colnames, eqs, observed, V, boot, n, mask, probit,
          model, pvc)


def _fprpts(o) -> int:
    k = int(o["fprpts"]) if isinstance(o.get("fprpts"), str) else 10
    if k < 2:
        raise StataError(198, "fprpts() must be at least 2")   # VERIFY
    return k


def _summary_spec(o) -> list[tuple[str, float | None]]:
    spec: list[tuple[str, float | None]] = []
    if o.get("auc"):
        spec.append(("auc", None))
    for kind in ("roc", "invroc", "pauc"):
        if isinstance(o.get(kind), str):
            for x in parse_numlist(str(o[kind])):
                if not 0 < x < 1 and not (kind == "pauc" and x == 1):
                    raise StataError(125, f"{kind}() values must be between 0 and 1")   # VERIFY
                spec.append((kind, float(x)))
    return spec or [("auc", None)]


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


def _bootstrap(s: "Session", o, statistic, ref, observed, level) -> dict:
    from ..estimation.bootstrap import _draw
    reps = int(o["breps"]) if isinstance(o.get("breps"), str) else 1000
    if reps < 2:
        raise StataError(198, "breps() must be an integer greater than 1")   # VERIFY
    if isinstance(o.get("bseed"), str):
        s.interp._execute_expanded(f"set seed {str(o['bseed']).strip()}")
    out = s.output
    groups = [np.arange(len(ref))] if o.get("nobstrata") else [np.flatnonzero(~ref), np.flatnonzero(ref)]
    dots = not o.get("nodots")
    if dots:
        out.write(f"\nBootstrap replications ({reps})\n", "text")
        out.write("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 \n", "text")
    R = np.full((reps, len(observed)), np.nan)
    for r in range(reps):
        idx = np.concatenate([g[_draw(len(g))] for g in groups])
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
    return {"reps": good, "n": len(good), "se": se, "bias": good.mean(axis=0) - observed,
            "normal": (observed - z * se, observed + z * se),
            "percentile": (np.percentile(good, a, axis=0), np.percentile(good, 100 - a, axis=0)),
            "bc": (np.array([x[0] for x in bc]), np.array([x[1] for x in bc]))}


# ---------------------------------------------------------------------------
# display and e()
# ---------------------------------------------------------------------------

_TITLES = {"auc": "Area under the ROC curve", "roc": "ROC curve", "invroc": "False-positive rate",
           "pauc": "Partial area under the ROC curve"}
_HEAD = {"auc": "AUC", "roc": "ROC", "invroc": "invROC", "pauc": "pAUC"}


def _f(x: float) -> str:
    return "." if x is None or not math.isfinite(x) else _g9(x)


def _g9(x: float) -> str:
    from .summarize import g9
    return g9(float(x))


def _header(s: "Session", title: str, n: int, boot) -> None:
    out = s.output
    out.write("\n", "text")
    items = [("Number of obs", f"{n:,}")]
    if boot is not None:
        items.append(("Replications", f"{boot['n']:,}"))
    for j, (lab, val) in enumerate(items):
        out.write(f"{(title if j == 0 else ''):<48}{lab:<18}= ", "text")
        out.write(f"{val:>10}\n", "result")


def _display(s: "Session", o, refname, classes, spec, colnames, observed, boot, level, n, probit,
             ctrlcov, roccov, model, pvc) -> None:
    out = s.output
    _header(s, "Bootstrap results" if boot is not None else "", n, boot)
    out.write("\n" + ("Parametric ROC estimation" if probit else "Nonparametric ROC estimation") + "\n\n", "text")
    if ctrlcov:
        out.write(f"Covariate control       : {'stratification' if model == 'strata' else 'linear regression'}\n", "text")
        out.write(f"Control variables       : {' '.join(ctrlcov)}\n", "text")
    out.write(f"Control standardization : {pvc}\n", "text")
    out.write(f"ROC method              : {'GLM' if probit else 'empirical'}\n", "text")
    if probit:
        out.write("ROC link                : probit\n", "text")
        if roccov:
            out.write(f"ROC covariates          : {' '.join(roccov)}\n", "text")
    per = len(colnames) // len(classes)
    for ci, c in enumerate(classes):
        sl = slice(ci * per, (ci + 1) * per)
        groups = [("coef", list(range(per)))] if probit else _groups(spec)
        for kind, js in groups:
            out.write("\n" + ("" if probit else _TITLES[kind] + "\n\n"), "text")
            out.write(f"   Status    : {refname}\n   Classifier: {c}\n\n", "text")
            _table(s, kind, [(colnames[sl][j], spec[j][1] if not probit else None, observed[sl][j],
                              None if boot is None else ci * per + j) for j in js], boot, level, probit)


def _groups(spec):
    order, out = [], {}
    for j, (kind, _) in enumerate(spec):
        if kind not in out:
            order.append(kind)
            out[kind] = []
        out[kind].append(j)
    return [(k, out[k]) for k in order]


def _table(s: "Session", kind: str, items, boot, level: float, probit: bool) -> None:
    """Observed values with the bootstrap bias, standard error and the normal (N), percentile
    (P) and bias-corrected (BC) intervals (VERIFY: layout)."""
    out = s.output
    head = "" if probit else _HEAD[kind]
    out.write("-" * 13 + "-" * 66 + "\n", "text")
    if boot is None:
        out.write(f"{head:>12} |      Coef.\n", "text")
    else:
        out.write(f"{'':>12} |   Observed{'Bootstrap':>23}\n", "text")
        out.write(f"{head:>12} |      Coef.       Bias    Std. Err.     [{level:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")
    for k, (name, x, obs, j) in enumerate(items):
        if probit:
            label = name
        elif kind == "auc":
            label = ""
        else:
            label = _g9(x).strip()
        out.write(f"{label:>12} |", "text")
        if j is None:
            out.write(f" {_f(obs):>10}\n", "result")
            continue
        bias, se = boot["bias"][j], boot["se"][j]
        rows = [("N", boot["normal"]), ("P", boot["percentile"]), ("BC", boot["bc"])]
        for r, (tag, (lo, hi)) in enumerate(rows):
            if r == 0:
                out.write(f" {_f(obs):>10} {_f(bias):>10}   {_f(se):>10}   ", "result")
            else:
                out.write(f"{'':>12} |{'':>38}", "text")
            out.write(f"{_f(lo[j]):>10}  {_f(hi[j]):>10}", "result")
            out.write(f"  ({tag})\n" if tag != "BC" else f" ({tag})\n", "text")
    out.write("-" * 79 + "\n", "text")


def _post(s: "Session", o, refname, classes, ctrlcov, roccov, colnames, eqs, observed, V, boot, n, mask,
          probit, model, pvc) -> None:
    from ..estimation.results import Estimates, post
    from .matrix import Matrix
    est = Estimates("rocreg", refname, colnames, np.asarray(observed, float), V, n, stat="z", eqnames=eqs)
    est.title = "Parametric ROC estimation" if probit else "Nonparametric ROC estimation"

    def mat(row, rown="y1"):
        m = Matrix(np.asarray(row, float).reshape(1, -1), [rown], list(colnames))
        m.coleq = list(eqs)
        return m

    def ci(pair):
        m = Matrix(np.vstack(pair), ["ll", "ul"], list(colnames))
        m.coleq = list(eqs)
        return m

    extra = []
    scalars = [("N", float(n))]
    if boot is not None:
        extra = [("bias", mat(boot["bias"])), ("se", mat(boot["se"])), ("ci_normal", ci(boot["normal"])),
                 ("ci_percentile", ci(boot["percentile"])), ("ci_bc", ci(boot["bc"]))]
        scalars.append(("breps", float(boot["n"])))
    # VERIFY: names and order of the e() results
    macros = [("cmd", "rocreg"), ("refvar", refname), ("classvars", " ".join(classes)),
              ("ctrlcov", " ".join(ctrlcov)), ("ctrlmodel", model if ctrlcov else ""),
              ("pvc", pvc), ("roccov", " ".join(roccov)), ("roccmd", "glm" if probit else ""),
              ("vce", "bootstrap" if boot is not None else ""), ("title", est.title), ("properties", "b V")]
    post(s, est, scalars, macros, mask, extra_mats=extra)
