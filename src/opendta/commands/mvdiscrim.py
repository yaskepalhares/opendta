"""Discriminant analysis ([MV] discrim lda, qda, knn, logistic and [MV] candisc) with its estat
subcommands and predict.

The statistics are in core/discrim.py. Layouts follow the manuals' examples; what they do not
show carries a VERIFY mark and is settled by the reference log of compat case 1906."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as sps

from ..core import discrim as DC
from ..core import missing as M
from ..core import mvtests as T
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lang.words import parse_numlist
from .mvcluster import _measure
from .mvcommon import g as gnum, matrix, numeric_names, post_e, sample
from .mvfactor import _split
from .registry import command
from .tabulate import _categories, _codes

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_TITLE = {"lda": "Linear discriminant analysis", "qda": "Quadratic discriminant analysis",
          "knn": "Kth-nearest-neighbor discriminant analysis", "logistic": "Logistic discriminant analysis"}


def _groups(s: "Session", p, names: list[str], group) -> tuple:
    """(X, w, weight type, mask, group name, categories, labels, codes) of the estimation sample."""
    if not isinstance(group, str):
        raise StataError(198, "option group() required")
    gname = s.expand_varlist(group)[0]
    gv = s.data.get(gname)
    X, w, wtype, mask = sample(s, p, names, extra=[] if gv.is_string else [gname])
    if gv.is_string:
        keep = gv.raw[mask] != ""
        X = X[keep]
        w = None if w is None else w[keep]
        mask[np.flatnonzero(mask)[~keep]] = False
    cats, labels, _ = _categories(s, gv, mask, missing=False, nolabel=False)
    if len(cats) < 2:
        raise StataError(2001, "insufficient groups")                                     # VERIFY
    return X, w, wtype, mask, gname, cats, labels, _codes(gv, cats)[mask]


def _priors(text, cnt: np.ndarray) -> np.ndarray:
    """priors(): equal (default), proportional or one value per group."""
    if text is None:
        return DC.priors_from("equal", cnt)
    t = str(text).strip()
    if t in ("equal", "proportional"):
        return DC.priors_from(t, cnt)
    v = parse_numlist(t.replace(",", " "))
    if len(v) != len(cnt) or min(v) < 0 or sum(v) <= 0:
        raise StataError(198, "invalid priors()")                                         # VERIFY
    return DC.priors_from(v, cnt)


@command("discrim", "discrim")
def cmd_discrim(s: "Session", args: str) -> None:
    words = args.strip().split(None, 1)
    if not words:
        raise StataError(198, "invalid syntax")
    sub = words[0].rstrip(",")
    if sub not in _TITLE:
        raise StataError(199, f"unrecognized discrim subcommand {sub}")                   # VERIFY
    p = parse_standard(args.strip()[len(words[0]):])
    spec = {"group": 1, "priors": 2, "notable": 5, "lootable": 3, "nopriors": 4}
    if sub == "knn":
        spec.update({"k": 1, "measure": 3})
    o = match_options(p.options, spec) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 1)
    X, w, wtype, mask, gname, cats, labels, codes = _groups(s, p, names, o.get("group"))
    G = len(cats)
    extra: dict = {}
    if sub in ("lda", "qda"):
        st = DC.group_stats(X, codes, G, w)
        if np.any(st["n"] < 2):
            raise StataError(2001, "insufficient observations in a group")                # VERIFY
        pri = _priors(o.get("priors"), st["n"])
        sc, _ = (DC.lda_scores if sub == "lda" else DC.qda_scores)(X, st, pri)
        P = DC.log_posteriors(sc)
        if sub == "lda":
            extra["cd"] = _canonical(st, names)
    else:
        cnt = np.bincount(codes, minlength=G).astype(float)
        st = {"n": cnt, "means": np.vstack([X[codes == g].mean(0) for g in range(G)])}
        pri = _priors(o.get("priors"), cnt)
        if sub == "knn":
            if not isinstance(o.get("k"), str):
                raise StataError(198, "option k() required")
            meas, power = _measure(o.get("measure"), "L2")
            extra = {"k": int(float(o["k"])), "measure": meas, "power": power}
            P = DC.knn_posteriors(X, codes, G, extra["k"], pri, X, meas, power)
        else:
            P = _logistic_probs(s, p, w, gname, names, mask, extra)
            # the logit model is fitted with the sample proportions as priors; Stata reweights by priors() (VERIFY)
            P = P * (pri * cnt.sum() / cnt)[None, :]
            P = P / P.sum(1, keepdims=True)
    fit = {"sub": sub, "names": names, "gname": gname, "cats": cats, "labels": labels, "codes": codes, "X": X, "w": w,
           "G": G, "priors": pri, "st": st, "P": P, "loo": None, "N": X.shape[0], "mask": mask, "extra": extra,
           "cd": extra.get("cd")}
    if not o.get("notable"):
        _class_table(s, fit, P, "Resubstitution classification summary", o, title=True)
        if o.get("lootable"):
            _class_table(s, fit, _loo(fit), "Leave-one-out classification summary", o)
    _post_discrim(s, fit, args, wtype)


def _logistic_probs(s: "Session", p, w, gname: str, names: list[str], mask: np.ndarray, extra_e: dict) -> np.ndarray:
    """Posterior probabilities of a multinomial logit fitted to the same observations."""
    if w is not None:
        raise StataError(101, "weights not allowed")                                      # VERIFY
    if s.data.get(gname).is_string:
        raise StataError(109, "type mismatch")
    cond = (f" if {p.if_}" if p.if_ else "") + (f" in {p.in_}" if p.in_ else "")
    from ..core.output import Capture
    out = s.output
    saved, cap = out._listeners, Capture()
    out._listeners = [cap]
    try:
        rc = s.run(f"mlogit {gname} {' '.join(names)}{cond}")
    finally:
        out._listeners = saved
    if rc:
        raise StataError(rc, "mlogit failed")
    for ln in cap.text.splitlines():
        if ln.startswith("Iteration"):
            out.write(ln + "\n", "text")
    est = s._estimates
    extra_e["ibaseout"] = float(s.e.get("ibaseout", 0))
    keep = [i for i, q in enumerate(est.eqnames) if np.any(est.b[[j for j, r in enumerate(est.eqnames) if r == q]] != 0)]
    extra_e["b"] = est.b[keep].copy()
    extra_e["bnames"] = [f"{est.eqnames[i]}:{est.names[i]}" for i in keep]
    Xd = np.column_stack([s.data.get(n).data[mask] for n in names] + [np.ones(int(mask.sum()))])
    eta = np.column_stack([Xd @ est.b[[i for i, q in enumerate(est.eqnames) if q == e]]
                           for e in dict.fromkeys(est.eqnames)])
    eta -= eta.max(1, keepdims=True)
    return np.exp(eta) / np.exp(eta).sum(1, keepdims=True)


def _loo(fit: dict) -> np.ndarray:
    """Leave-one-out posterior probabilities, computed when first asked for."""
    if fit["loo"] is None:
        sub, G = fit["sub"], fit["G"]
        if sub in ("lda", "qda"):
            sc = DC.loo_scores(sub, fit["X"], fit["codes"], G, fit["w"], fit["priors"])
            fit["loo"] = DC.log_posteriors(sc)
        elif sub == "knn":
            e = fit["extra"]
            fit["loo"] = DC.knn_posteriors(fit["X"], fit["codes"], G, e["k"], fit["priors"], fit["X"], e["measure"],
                                           e["power"], loo=True)
        else:
            raise StataError(198, "leave-one-out classification is not available after discrim logistic")   # VERIFY
    return fit["loo"]


_PAD = " " * 4


def _class_table(s: "Session", fit: dict, P: np.ndarray, caption: str, o: dict, title: bool = False,
                 loo: bool = False) -> None:
    out = s.output
    labs = [str(x) for x in fit["labels"]]
    G = fit["G"]
    Tr = DC.class_table(fit["codes"], DC.classify(P), G, fit["w"])
    if title:
        out.write(f"\n{_TITLE[fit['sub']]}\n", "text")
    out.write(f"{caption}\n\n{_PAD}+---------+\n{_PAD}| Key     |\n{_PAD}|---------|\n{_PAD}| Number  |\n"
              f"{_PAD}| Percent |\n{_PAD}+---------+\n", "text")
    W = max(8, max(len(lb) for lb in labs) + 1)
    rule = f"{_PAD}{'-' * 13}+{'-' * (W * G)}+{'-' * 7}"
    bar = f"{'':>17}|{'':>{W * G}}|"
    out.write(f"{'':>17}| {'LOO ' if loo else ''}Classified\n", "text")
    out.write(f"{_PAD}{'True ' + fit['gname']:<12} |" + "".join(f"{lb:>{W - 1}} " for lb in labs) + f"|{'Total':>7}\n{rule}\n", "text")
    rows = list(Tr) + [Tr.sum(0)]
    for i_, (name, row) in enumerate(zip(labs + ["Total"], rows)):
        tot = row.sum()
        out.write(f"{_PAD}{name[:12]:>12} |" + "".join(f"{gnum(x):>{W - 1}} " for x in row) + f"|{gnum(tot):>7}\n", "text")
        out.write(f"{'':>17}|" + "".join(f"{100 * x / tot:>{W - 1}.2f} " for x in row) + f"|{100.0:>7.2f}\n", "text")
        out.write((rule if i_ == G - 1 else bar) + "\n", "text")
    if not o.get("nopriors"):
        out.write(f"{_PAD}{'Priors':>12} |" + "".join(f" {x:.4f} " for x in fit["priors"]) + "|\n", "text")
    out.write("\n", "text")


def _canonical(st: dict, names: list[str]) -> dict:
    """Canonical discriminant functions of the groups in `st`, with the signs of Stata (the first
    coefficient of every function is positive, VERIFY) and the matrices candisc stores."""
    cd = DC.candisc(st)
    coef = cd["coef"].copy()
    for j in range(coef.shape[1]):
        if coef[0, j] < 0:
            coef[:, j] *= -1
    n, means = st["n"], st["means"]
    N = n.sum()
    grand = (n[:, None] * means).sum(0) / N
    S = st["pooled"]
    sd_w = np.sqrt(np.diag(S))
    sd_t = np.sqrt(np.diag(cd["W"] + cd["B"]) / (N - 1.0))
    cd.update({"coef": coef, "std": coef * sd_w[:, None], "totalstd": coef * sd_t[:, None],
               "unstd": np.vstack([coef, -grand @ coef]), "structure": (S @ coef) / sd_w[:, None],
               "cmeans": (means - grand) @ coef})
    return cd


def _post_discrim(s: "Session", fit: dict, args: str, wtype, cmd: str = "discrim") -> None:
    st = fit["st"]
    sub = fit["sub"]
    cn = [str(c) for c in fit["labels"]]
    names = fit["names"]
    G, p = fit["G"], len(names)
    gvals = [c if isinstance(c, (int, float)) else i + 1.0 for i, c in enumerate(fit["cats"])]
    scalars = []
    mats = []
    if sub == "lda":
        cd = fit["cd"]
        f = len(cd["ev"])
        ev = cd["ev"]
        N, nvar = float(st["n"].sum()), p
        rows = []
        cum = 0.0
        for j in range(f):
            r = T.mv_statistics(ev[j:], nvar - j, G - 1 - j, N - G)[0]
            cum += ev[j] / ev.sum()
            rows.append([cd["canon_corr"][j], ev[j], ev[j] / ev.sum(), cum, r["stat"], r["F"], r["df1"], r["df2"], r["p"],
                         {"e": 1.0, "a": 2.0, "u": 3.0}[r["flag"]]])
        fc = [f"function{j + 1}" for j in range(f)]
        Sinv = np.linalg.inv(st["pooled"])
        w_, V_ = np.linalg.eigh(cd["W"])
        C = np.vstack([Sinv @ st["means"].T, (-0.5 * np.einsum("gi,ij,gj->g", st["means"], Sinv, st["means"]) + np.log(fit["priors"]))[None, :]])
        scalars.append(("f", f))
    scalars += [("N", fit["N"]), ("N_groups", G), ("k", p)]
    if sub == "knn":
        scalars.append(("k_nn", fit["extra"]["k"]))
    if sub == "logistic":
        scalars.append(("ibaseout", fit["extra"]["ibaseout"]))
    lab = " ".join(f'`"{c}"\'' for c in cn)
    pred = f"discrim_{sub}_p"
    est_cmd = f"discrim_{sub}_estat"
    title = "Canonical linear discriminant analysis" if cmd == "candisc" else _TITLE[sub]
    macros = [("cmd", cmd), ("marginsnotok", "_ALL"), ("cmdline", f"{cmd} " + args.strip())]
    if cmd == "discrim":
        macros.append(("subcmd", sub))
    macros += [("title", title), ("predict", pred), ("estat_cmd", est_cmd)]
    if sub == "logistic":
        macros += [("ties", "missing"), ("grouplabels", lab)]
    else:
        macros += [("grouplabels", lab), ("ties", "missing")]
    macros += [("varlist", " ".join(names)), ("groupvar", fit["gname"]),
               ("properties", "b noV" if sub == "logistic" else "nob noV eigen" if sub == "lda" else "nob noV")]
    if sub == "logistic":
        b = fit["extra"]["b"]
        mats.append(("b", matrix(b.reshape(1, -1), ["y1"], fit["extra"]["bnames"])))
    mats += [("grouppriors", matrix(fit["priors"].reshape(1, -1), ["grouppriors"], cn)),
             ("groupcounts", matrix(st["n"].reshape(1, -1), ["groupcounts"], cn)),
             ("groupvalues", matrix(np.array(gvals, dtype=float).reshape(1, -1), ["groupvalues"], cn))]
    if sub == "lda":
        n = st["n"]
        sscp = [st["covs"][g] * (n[g] - 1.0) for g in range(G)]
        Wm = cd["W"]
        lamS, VS = np.linalg.eigh(Sinv)
        mats += [("C", matrix(C, names + ["_cons"], cn)),
                 ("W_eigvecs", matrix(V_[:, ::-1], names, [f"c{i + 1}" for i in range(p)])),
                 ("W_eigvals", matrix(w_[::-1].reshape(1, -1), ["W_eigvals"], [f"c{i + 1}" for i in range(p)])),
                 ("candisc_stat", matrix(np.array(rows), fc, [f"c{i + 1}" for i in range(10)])),     # VERIFY: columns
                 ("cmeans", matrix(cd["cmeans"], cn, fc)), ("canstruct", matrix(cd["structure"], names, fc)),
                 ("L_totalstd", matrix(cd["totalstd"], names, fc)), ("L_std", matrix(cd["std"], names, fc)),
                 ("L_unstd", matrix(cd["unstd"], names + ["_cons"], fc)), ("L_raw", matrix(cd["coef"], names, fc)),
                 ("Ev", matrix(ev.reshape(1, -1), ["Ev"], fc)),
                 ("sqrtSinv", matrix(VS @ np.diag(np.sqrt(lamS)) @ VS.T, names, names)),
                 ("Sinv", matrix(Sinv, names, names)), ("S", matrix(st["pooled"], names, names)),
                 ("means", matrix(st["means"], cn, names)), ("SSCP_B", matrix(cd["B"], names, names)),
                 ("SSCP_W", matrix(Wm, names, names))]
        mats += [(f"SSCP_W{g + 1}", matrix(sscp[g], names, names)) for g in reversed(range(G))]
        mats += [("SSCP_T", matrix(cd["B"] + Wm, names, names))]
    elif sub != "logistic":
        mats.append(("means", matrix(st["means"], cn, names)))
    post_e(s, "discrim", scalars, macros, mats, fit["mask"], fit["N"],
           estat={"classtable": (5, _estat_classtable), "errorrate": (3, _estat_errorrate), "list": (1, _estat_list),
                  "grdistances": (3, _estat_grdist), "covariance": (3, _estat_cov)},
           predict=_predict_one,
           predict_multi={"pr": lambda ss, ee: _predict_multi(ss, ee, "pr"),
                          "mahalanobis": lambda ss, ee: _predict_multi(ss, ee, "mahalanobis")},
           predict_opts={"classification": 5, "mahalanobis": 3}, extra={"discrim": fit})


# ---------------------------------------------------------------------------
# estat
# ---------------------------------------------------------------------------

def _fit(est) -> dict:
    return est.extra["discrim"]


def _estat_classtable(s, est, rest):
    fit = _fit(est)
    opt = _split(rest)[1]
    o = match_options(opt, {"loo": 3, "nopriors": 4}) if opt.strip() else {}
    if o.get("loo"):
        _class_table(s, fit, _loo(fit), "Leave-one-out classification table", o, loo=True)
    else:
        _class_table(s, fit, fit["P"], "Resubstitution classification table", o)


def _rate_row(label: str, vals, tot) -> str:
    return f"{_PAD}{label:>12} |" + "".join(f"{gnum(x):>10} " for x in vals) + f"|{gnum(tot) if tot is not None else '':>10}\n"


def _estat_errorrate(s, est, rest):
    fit = _fit(est)
    opt = _split(rest)[1]
    o = match_options(opt, {"pp": 2, "nopriors": 4}) if opt.strip() else {}
    G = fit["G"]
    labs = [str(x) for x in fit["labels"]]
    out = s.output
    rule = f"{_PAD}{'-' * 13}+{'-' * (11 * G)}+{'-' * 10}\n"
    head = "".join(f"{lb:>10} " for lb in labs) + f"|{'Total':>10}\n"
    if o.get("pp"):
        rates = DC.smoothed_error_rates(fit["P"], fit["codes"], G, fit["priors"], fit["w"])
        out.write(f"\nError rate estimated from posterior probabilities\n\n{'':>17}| {fit['gname']}\n", "text")
        out.write(f"{_PAD}{'Error Rate':<12} |{head}{rule}", "text")
        out.write(_rate_row("Stratified", *rates["stratified"]), "text")
        out.write(f"{'':>17}|{'':>{11 * G}}|{'':>10}\n", "text")
        out.write(_rate_row("Unstratified", *rates["unstratified"]), "text")
        err, tot = rates["stratified"]
    else:
        Tr = DC.class_table(fit["codes"], DC.classify(fit["P"]), G, fit["w"])
        err, tot = DC.error_rates(Tr, fit["priors"])
        out.write(f"\nError rate estimated by error count\n\n{'':>17}| {fit['gname']}\n{'':>17}|{head}{rule}", "text")
        out.write(_rate_row("Error rate", err, tot), "text")
    if not o.get("nopriors"):
        out.write(rule + _rate_row("Priors", fit["priors"], None), "text")
    out.write("\n", "text")
    s.r = {"Rate": matrix(np.append(err, tot).reshape(1, -1), ["r1"], labs + ["Total"]),
           "Priors": matrix(fit["priors"].reshape(1, -1), ["r1"], labs)}


def _estat_list(s, est, rest):
    fit = _fit(est)
    opt = _split(rest)[1]
    o = match_options(opt, {"misclassified": 3, "probabilities": 4}) if opt.strip() else {}
    P = fit["P"]
    pred = DC.classify(P)
    sel = np.flatnonzero(pred != fit["codes"]) if o.get("misclassified") else np.arange(fit["N"])
    obs = np.flatnonzero(fit["mask"]) + 1
    labs = [str(x) for x in fit["labels"]]
    G = fit["G"]
    ow = max(4, len(str(int(obs.max()))))
    inner = 8 * G
    line = lambda ch: f"|{ch * (ow + 1)}+{ch * 18}+{ch * inner}|\n"
    out = s.output
    border = f"+{'-' * (ow + 1 + 1 + 18 + 1 + inner)}+\n"
    out.write(border, "text")
    out.write(f"|{'':>{ow + 1}}| {'Classification':<17}| {'Probabilities':<{inner - 1}}|\n", "text")
    out.write(f"|{'':>{ow + 1}}|{'':>18}|{'':>{inner}}|\n", "text")
    out.write(f"|{'Obs.':>{ow + 1}}|{'True  Class.':^18}|" + "".join(f"{lb:>7} " for lb in labs) + "|\n", "text")
    out.write(line("-"), "text")
    for r_, i in enumerate(sel):
        if r_ and r_ % 5 == 0:
            out.write(line("-"), "text")
        star = "*" if pred[i] != fit["codes"][i] else " "
        out.write(f"|{obs[i]:>{ow}} |{labs[fit['codes'][i]][:7]:>7}{labs[pred[i]][:7]:>8} {star} |"
                  + "".join(f" {P[i, j]:.4f} " for j in range(G)) + "|\n", "text")
    out.write(border, "text")
    out.write("* indicates misclassified observations\n", "text")


def _mat_table(out, title: str, corner: str, cols: list[str], rows: list[str], M_, *, lower: bool = False,
               top: str | None = None, trim: int = 0, right: bool = False, fmt=None, blank: bool = True) -> None:
    """Matrix as estat shows it: a stub of 12 characters and cells of 10 followed by a space."""
    fmt = fmt or gnum
    if title:
        out.write(f"\n{title}\n" + ("\n" if blank else ""), "text")
    if top is not None:
        out.write(f"{'':>17}| {top}\n", "text")
    out.write(f"{_PAD}{corner:>12} |" if right else f"{_PAD}{corner:<12} |", "text")
    out.write("".join(f"{c:>10} " for c in cols) + "\n", "text")
    out.write(f"{_PAD}{'-' * 13}+{'-' * (11 * len(cols) - trim)}\n", "text")
    for i, r_ in enumerate(rows):
        m = i + 1 if lower else len(cols)
        out.write(f"{_PAD}{r_[:12]:>12} |" + "".join(f"{fmt(M_[i, j]):>10} " for j in range(m)) + "\n", "text")


def _estat_grdist(s, est, rest):
    fit = _fit(est)
    if fit["sub"] != "lda":
        raise StataError(321, "estat grdistances is only available after discrim lda")   # VERIFY
    st = fit["st"]
    G, p, n = fit["G"], len(fit["names"]), st["n"]
    labs = [str(x) for x in fit["labels"]]
    D2 = DC.mahalanobis_sq(st["means"], st["means"], st["pooled"])
    df2 = n.sum() - G - p + 1
    F = np.zeros((G, G))
    for i in range(G):
        for j in range(G):
            F[i, j] = df2 / ((n.sum() - G) * p) * n[i] * n[j] / (n[i] + n[j]) * D2[i, j]
    pv = sps.f.sf(F, p, df2)
    pv[np.arange(G), np.arange(G)] = 1.0
    _mat_table(s.output, "Mahalanobis squared distances between groups", fit["gname"], labs, labs, D2, lower=True,
               top=fit["gname"], trim=1)
    s.output.write("\n", "text")
    s.r = {"dist": matrix(D2, labs, labs), "F": matrix(F, labs, labs), "p": matrix(pv, labs, labs)}


def _estat_cov(s, est, rest):
    fit = _fit(est)
    names = fit["names"]
    if fit["sub"] == "qda":
        s.output.write("\nGroup covariance matrices\n", "text")
        for lb, S in zip(fit["labels"], fit["st"]["covs"]):
            s.output.write(f"\n  {fit['gname']}: {lb}\n", "text")
            _mat_table(s.output, "", "", names, names, S, lower=True)
    else:
        _mat_table(s.output, "Pooled within-group covariance matrix", "", names, names, fit["st"]["pooled"], lower=True)
        s.r = {"Sigma": matrix(fit["st"]["pooled"], names, names)}
    s.output.write("\n", "text")


# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

def _scores(s: "Session", fit: dict):
    """(posterior probabilities, squared Mahalanobis distances or None, usable rows) for every observation."""
    ds = s.data
    G = fit["G"]
    if fit["sub"] == "logistic":
        if len(fit["mask"]) != ds.nobs:
            raise StataError(198, "the data have changed since discrim logistic")          # VERIFY
        full = np.zeros((ds.nobs, G))
        full[fit["mask"]] = fit["P"]
        return full, None, fit["mask"].copy()
    X = np.column_stack([ds.get(n).data.astype(float) for n in fit["names"]])
    ok = np.all(X < SYS, axis=1)
    X = np.where(ok[:, None], X, 0.0)
    if fit["sub"] == "knn":
        e = fit["extra"]
        return DC.knn_posteriors(fit["X"], fit["codes"], G, e["k"], fit["priors"], X, e["measure"], e["power"]), None, ok
    sc, d2 = (DC.lda_scores if fit["sub"] == "lda" else DC.qda_scores)(X, fit["st"], fit["priors"])
    return DC.log_posteriors(sc), d2, ok


def _predict_one(s, est, stat, mask):
    fit = _fit(est)
    if stat not in ("", "classification"):
        raise StataError(198, f"option {stat} not allowed")
    P, _, ok = _scores(s, fit)
    gvals = np.array([c if isinstance(c, (int, float)) else i + 1.0 for i, c in enumerate(fit["cats"])], dtype=float)
    note = "" if stat else "(option classification assumed; group classification)"            # VERIFY
    return np.where(ok, gvals[DC.classify(P)], SYS), note


def _predict_multi(s, est, kind):
    fit = _fit(est)
    if kind == "mahalanobis" and fit["sub"] not in ("lda", "qda"):
        raise StataError(198, "option mahalanobis is only available after discrim lda and qda")   # VERIFY
    P, d2, ok = _scores(s, fit)
    A = P if kind == "pr" else d2
    G = fit["G"]
    return G, lambda mask: [np.where(ok, A[:, j], SYS) for j in range(G)]


# ---------------------------------------------------------------------------
# candisc
# ---------------------------------------------------------------------------

def _gw(x: float, width: int) -> str:
    """%g with as many significant digits as fit in `width` characters, without a leading zero."""
    for d in range(9, 0, -1):
        t_ = f"{x:.{d}g}"
        t_ = t_.replace("0.", ".", 1) if t_.startswith(("0.", "-0.")) else t_
        if len(t_) <= width:
            return t_
    return t_


@command("candisc", "candisc")
def cmd_candisc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"group": 1, "priors": 2, "notable": 5, "lootable": 3, "nopriors": 4}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 1)
    X, w, wtype, mask, gname, cats, labels, codes = _groups(s, p, names, o.get("group"))
    G = len(cats)
    st = DC.group_stats(X, codes, G, w)
    if np.any(st["n"] < 2):
        raise StataError(2001, "insufficient observations in a group")                    # VERIFY
    pri = _priors(o.get("priors"), st["n"])
    sc, _ = DC.lda_scores(X, st, pri)
    P = DC.log_posteriors(sc)
    cd = _canonical(st, names)
    fit = {"sub": "lda", "names": names, "gname": gname, "cats": cats, "labels": labels, "codes": codes, "X": X, "w": w,
           "G": G, "priors": pri, "st": st, "P": P, "loo": None, "N": X.shape[0], "mask": mask, "extra": {}, "cd": cd}
    ev = cd["ev"]
    N, m, nvar = float(st["n"].sum()), len(ev), len(names)
    out = s.output
    out.write("\nCanonical linear discriminant analysis\n\n", "text")
    out.write(f"{'':>6}|{'':>33}| Like-\n{'':>6}| Canon.   Eigen-     Variance    | lihood\n", "text")
    out.write("  Fcn | Corr.    value   Prop.   Cumul. | Ratio     F      df1    df2  Prob>F\n", "text")
    out.write(f"  ----+{'-' * 33}+{'-' * 36}\n", "text")
    cum = 0.0
    flags = []
    for j in range(m):
        r = T.mv_statistics(ev[j:], nvar - j, G - 1 - j, N - G)[0]
        cum += ev[j] / ev.sum()
        flags.append(r["flag"])
        out.write(f"{j + 1:>5} | {cd['canon_corr'][j]:.4f}  {_gw(ev[j], 7):>7}  {ev[j] / ev.sum():.4f}  {cum:.4f} | "
                  f"{r['stat']:.4f}{_gw(r['F'], 6):>8}{gnum(r['df1']):>6}{gnum(r['df2']):>7}  {r['p']:.4f} {r['flag']}\n", "text")
    out.write(f"  {'-' * 75}\n", "text")
    keys = {"e": "e = exact F", "a": "a = approximate F", "u": "u = upper bound on F"}      # VERIFY: footer with several flags
    note = ", ".join(keys[k_] for k_ in "eau" if k_ in flags)
    left = "  Ho: this and smaller canon. corr. are zero;"
    out.write(f"{left}{note:>{77 - len(left)}}\n", "text")
    fcn = [f"function{j + 1}" for j in range(m)]
    labs = [str(x) for x in labels]
    _mat_table(out, "Standardized canonical discriminant function coefficients", "", fcn, names, cd["std"])
    _mat_table(out, "Canonical structure", "", fcn, names, cd["structure"])
    _mat_table(out, "Group means on canonical variables", gname, fcn, labs, cd["cmeans"], right=True)
    _class_table(s, fit, P, "Resubstitution classification summary", o, title=False) if not o.get("notable") else None
    _post_discrim(s, fit, args, wtype, cmd="candisc")
