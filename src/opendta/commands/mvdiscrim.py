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
from .mvfactor import _print_matrix, _split
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
            P = _logistic_probs(s, p, w, gname, names, mask)
    fit = {"sub": sub, "names": names, "gname": gname, "cats": cats, "labels": labels, "codes": codes, "X": X, "w": w,
           "G": G, "priors": pri, "st": st, "P": P, "loo": None, "N": X.shape[0], "mask": mask, "extra": extra}
    if not o.get("notable"):
        _class_table(s, fit, P, "Resubstitution classification summary", o, title=True)
        if o.get("lootable"):
            _class_table(s, fit, _loo(fit), "Leave-one-out classification summary", o)
    _post_discrim(s, fit, args, wtype)


def _logistic_probs(s: "Session", p, w, gname: str, names: list[str], mask: np.ndarray) -> np.ndarray:
    """Posterior probabilities of a multinomial logit fitted to the same observations."""
    if w is not None:
        raise StataError(101, "weights not allowed")                                      # VERIFY
    if s.data.get(gname).is_string:
        raise StataError(109, "type mismatch")
    cond = (f" if {p.if_}" if p.if_ else "") + (f" in {p.in_}" if p.in_ else "")
    rc = s.run(f"quietly mlogit {gname} {' '.join(names)}{cond}")
    if rc:
        raise StataError(rc, "mlogit failed")
    est = s._estimates
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


def _class_table(s: "Session", fit: dict, P: np.ndarray, caption: str, o: dict, title: bool = False) -> None:
    """Classification summary (VERIFY: layout)."""
    out = s.output
    labs = [str(x) for x in fit["labels"]]
    G = fit["G"]
    Tr = DC.class_table(fit["codes"], DC.classify(P), G, fit["w"])
    if title:
        out.write(f"\n{_TITLE[fit['sub']]}\n", "text")
    out.write(f"{caption}\n\n  +---------+\n  | Key     |\n  |---------|\n  | Number  |\n  | Percent |\n  +---------+\n\n", "text")
    W = max(8, max(len(lb) for lb in labs)) + 1
    rule = "-" * 13 + "+" + "-" * (W * G + 1) + "+" + "-" * 10
    out.write(f"{'':>12} |{'Classified':>{W * G // 2 + 5}}\n", "text")
    out.write(f"{'True group':>12} |" + "".join(f"{lb:>{W}}" for lb in labs) + f" |{'Total':>9}\n{rule}\n", "text")
    rows = list(Tr) + [Tr.sum(0)]
    for name, row in zip(labs + ["Total"], rows):
        tot = row.sum()
        out.write(f"{name[:12]:>12} |" + "".join(f"{gnum(x):>{W}}" for x in row) + f" |{gnum(tot):>9}\n", "text")
        out.write(f"{'':>12} |" + "".join(f"{100 * x / tot:>{W}.2f}" for x in row) + f" |{100.0:>9.2f}\n{rule}\n", "text")
    if not o.get("nopriors"):
        out.write(f"{'Priors':>12} |" + "".join(f"{x:>{W}.4f}" for x in fit["priors"]) + " |\n", "text")
    out.write("\n", "text")


def _post_discrim(s: "Session", fit: dict, args: str, wtype) -> None:
    st = fit["st"]
    sub = fit["sub"]
    cn = [str(c) for c in fit["labels"]]
    scalars = [("N", fit["N"]), ("N_groups", fit["G"]), ("k", len(fit["names"]))]
    if sub == "knn":
        scalars.append(("k_nn", fit["extra"]["k"]))
    macros = [("cmdline", "discrim " + args.strip()), ("cmd", "discrim"), ("subcmd", sub), ("groupvar", fit["gname"]),
              ("varlist", " ".join(fit["names"])), ("title", _TITLE[sub]), ("estat_cmd", "discrim_estat"),
              ("predict", "discrim_p"), ("properties", "nob noV"), ("wtype", wtype)]
    mats = [("groupcounts", matrix(st["n"].reshape(-1, 1), cn, ["count"])),
            ("grpmeans", matrix(st["means"], cn, fit["names"])),
            ("priors", matrix(fit["priors"].reshape(1, -1), ["priors"], cn))]
    if sub == "lda":
        mats.append(("Sigma", matrix(st["pooled"], fit["names"], fit["names"])))
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
        _class_table(s, fit, _loo(fit), "Leave-one-out classification summary", o)
    else:
        _class_table(s, fit, fit["P"], "Resubstitution classification summary", o)


def _estat_errorrate(s, est, rest):
    fit = _fit(est)
    opt = _split(rest)[1]
    o = match_options(opt, {"pp": 2, "nopriors": 4}) if opt.strip() else {}
    G = fit["G"]
    labs = [str(x) for x in fit["labels"]]
    if o.get("pp"):
        err, tot = DC.smoothed_error_rates(fit["P"], fit["codes"], G, fit["priors"])
        title = "Error rate estimated by smoothing"
    else:
        Tr = DC.class_table(fit["codes"], DC.classify(fit["P"]), G, fit["w"])
        err, tot = DC.error_rates(Tr, fit["priors"])
        title = "Error rate estimated by error count"
    out = s.output
    W = max(9, max(len(lb) for lb in labs) + 1)
    out.write(f"\n{title}\n\n{'':>12} |" + "".join(f"{lb:>{W}}" for lb in labs) + f" |{'Total':>9}\n", "text")
    out.write("-" * 13 + "+" + "-" * (W * G + 1) + "+" + "-" * 10 + "\n", "text")
    out.write(f"{'Error rate':>12} |" + "".join(f"{x:>{W}.4f}" for x in err) + f" |{tot:>9.4f}\n", "text")
    if not o.get("nopriors"):
        out.write(f"{'Priors':>12} |" + "".join(f"{x:>{W}.4f}" for x in fit["priors"]) + " |\n", "text")
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
    probs = bool(o.get("probabilities"))
    head = f"{'Obs':>6} | {'True':>10} {'Classified':>11}" + ("".join(f"{lb[:9]:>10}" for lb in labs) if probs else "")
    out = s.output
    out.write(f"\n{head}\n{'-' * 7}+{'-' * (len(head) - 7)}\n", "text")
    for i in sel:
        star = "*" if pred[i] != fit["codes"][i] else " "
        line = f"{obs[i]:>6} | {labs[fit['codes'][i]][:10]:>10} {star}{labs[pred[i]][:10]:>10}"
        out.write(line + ("".join(f"{P[i, j]:>10.4f}" for j in range(fit["G"])) if probs else "") + "\n", "text")
    out.write("\n", "text")


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
    _print_matrix(s.output, "Mahalanobis squared distances between group means", labs, D2, labs)
    _print_matrix(s.output, "F statistics", labs, F, labs)
    _print_matrix(s.output, "Prob > F", labs, pv, labs)
    s.output.write("\n", "text")
    s.r = {"dist": matrix(D2, labs, labs), "F": matrix(F, labs, labs), "p": matrix(pv, labs, labs)}


def _estat_cov(s, est, rest):
    fit = _fit(est)
    names = fit["names"]
    if fit["sub"] == "qda":
        for lb, S in zip(fit["labels"], fit["st"]["covs"]):
            _print_matrix(s.output, f"Covariance matrix of group {lb}", names, S)
    else:
        _print_matrix(s.output, "Pooled within-group covariance matrix", names, fit["st"]["pooled"])
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

@command("candisc", "candisc")
def cmd_candisc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"group": 1}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 1)
    X, w, wtype, mask, gname, cats, labels, codes = _groups(s, p, names, o.get("group"))
    G = len(cats)
    st = DC.group_stats(X, codes, G, w)
    cd = DC.candisc(st)
    ev = cd["ev"]
    N, m, nvar = float(st["n"].sum()), len(ev), len(names)
    out = s.output
    cols = [11, 9, 9, 9, 9, 10, 6, 9, 9]                                                  # VERIFY: layout
    out.write("\nCanonical linear discriminant analysis\n\n", "text")
    out.write(f"{'':>4} |" + "".join(f"{h:>{c}}" for h, c in zip(["", "Canon.", "", "", "Like-", "Approx", "", "", ""], cols)) + "\n", "text")
    out.write(f"{'Fcn':>4} |" + "".join(f"{h:>{c}}" for h, c in zip(
        ["Eigenvalue", "Corr.", "Prop.", "Cumul.", "lihood", "F", "df1", "df2", "Prob>F"], cols)) + "\n", "text")
    rule = "-" * 5 + "+" + "-" * sum(cols)
    out.write(rule + "\n", "text")
    cum = 0.0
    for j in range(m):
        r = T.mv_statistics(ev[j:], nvar - j, G - 1 - j, N - G)[0]
        cum += ev[j] / ev.sum()
        out.write(f"{j + 1:>4} |{ev[j]:>11.4f}{cd['canon_corr'][j]:>9.4f}{ev[j] / ev.sum():>9.4f}{cum:>9.4f}"
                  f"{r['stat']:>9.4f}{r['F']:>10.4f}{r['df1']:>6.0f}{r['df2']:>9.1f}{r['p']:>9.4f} {r['flag']}\n", "text")
    out.write(f"{rule}\n{'':>5}H0: this and smaller canonical correlations are zero; e = exact F, a = approximate F\n", "text")
    fcn = [f"function{j + 1}" for j in range(m)]
    sd = np.sqrt(np.diag(st["pooled"]))
    _print_matrix(out, "Standardized canonical discriminant function coefficients", names, cd["coef"] * sd[:, None], fcn)
    _print_matrix(out, "Canonical discriminant function coefficients", names, cd["coef"], fcn)
    out.write("\n", "text")
    macros = [("cmdline", "candisc " + args.strip()), ("cmd", "candisc"), ("groupvar", gname), ("varlist", " ".join(names)),
              ("title", "Canonical linear discriminant analysis"), ("properties", "nob noV"), ("wtype", wtype)]
    mats = [("canon_corr", matrix(cd["canon_corr"].reshape(1, -1), ["canon_corr"], fcn)),
            ("eigenvalues", matrix(ev.reshape(1, -1), ["eigenvalues"], fcn)),
            ("coef", matrix(cd["coef"], names, fcn))]
    post_e(s, "candisc", [("N", N), ("N_groups", G)], macros, mats, mask, int(N))
