"""Correspondence analysis (ca, mca), multidimensional scaling (mds, mdsmat), procrustes and biplot.

The statistics are in core/ordination.py. Layouts follow the manuals' examples; what they do not
show carries a VERIFY mark and is settled by the reference log of compat case 1906."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import cluster as K
from ..core import missing as M
from ..core import ordination as O
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from .mvcluster import _measure
from .mvcommon import chi2_p, g as gnum, matrix, numeric_names, post_e, sample
from .mvfactor import _header, _right
from .registry import command
from .summarize import weights
from .tabulate import _categories, _codes
from ._util import touse

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _cat_sample(s: "Session", p, names: list[str]) -> tuple:
    """(codes per variable, labels per variable, weights, weight type, mask) of the observations
    with every categorical variable present."""
    ds = s.data
    mask = touse(s, p)
    for n in names:
        v = ds.get(n)
        mask &= (v.raw != "") if v.is_string else (v.data < SYS)
    w, wtype, mask = weights(s, p, mask, ("fweight",))
    if not mask.any():
        raise StataError(2000, "no observations")
    codes, labels = [], []
    for n in names:
        v = ds.get(n)
        cats, labs, _ = _categories(s, v, mask, missing=False, nolabel=False)
        codes.append(_codes(v, cats)[mask])
        labels.append([str(x) for x in labs])
    return codes, labels, (None if w is None else w[mask]), wtype, mask


def _dims(o: dict, available: int, default: int = 2) -> int:
    f = int(float(o["dimensions"])) if isinstance(o.get("dimensions"), str) else default
    if f < 1:
        raise StataError(198, "dimensions() must be positive")                            # VERIFY
    return min(f, available)


# ---------------------------------------------------------------------------
# points tables shared by ca and mca
# ---------------------------------------------------------------------------

def _point_table(out, title: str, labs: list[str], st: dict, C: np.ndarray, f: int) -> None:
    """Mass, quality and inertia share of every category and, for each retained dimension, its
    coordinate, squared correlation and contribution (VERIFY: layout)."""
    out.write(f"\n{title}\n\n", "text")
    out.write(f"{'':>15} |{'overall':^30}|" + "|".join(f"{'dimension_' + str(j + 1):^30}" for j in range(f)) + "|\n", "text")
    out.write(f"{'Categories':>15} |{'mass':>9}{'quality':>10}{'%inert':>10} |"
              + " |".join(f"{'coord':>9}{'sqcorr':>10}{'contrib':>10}" for _ in range(f)) + " |\n", "text")
    out.write("-" * 16 + "+" + "-" * 30 + ("+" + "-" * 30) * f + "\n", "text")
    for i, lab in enumerate(labs):
        line = f"{lab[:15]:>15} |{st['mass'][i]:>9.3f}{st['quality'][i]:>10.3f}{st['inertia'][i]:>10.3f} |"
        line += " |".join(f"{C[i, j]:>9.3f}{st['sqcorr'][i, j]:>10.3f}{st['contrib'][i, j]:>10.3f}" for j in range(f)) + " |"
        out.write(line + "\n", "text")


def _coord_table(out, title: str, labs: list[str], C: np.ndarray, f: int) -> None:
    out.write(f"\n{title}\n\n{'Categories':>15} |" + "".join(f"{'dim ' + str(j + 1):>10}" for j in range(f)) + "\n", "text")
    out.write("-" * 16 + "+" + "-" * (10 * f + 1) + "\n", "text")
    for i, lab in enumerate(labs):
        out.write(f"{lab[:15]:>15} |" + "".join(f"{C[i, j]:>10.4f}" for j in range(f)) + "\n", "text")


def _estat_table(s, est, rest):
    ctx = est.extra["ord"]
    for title, labs, st, C in ctx["tables"]:
        _point_table(s.output, title, labs, st, C, ctx["f"])
    s.output.write("\n", "text")


def _estat_coordinates(s, est, rest):
    ctx = est.extra["ord"]
    for title, labs, st, C in ctx["tables"]:
        _coord_table(s.output, "Coordinates in " + title.split(" in ", 1)[-1], labs, C, ctx["f"])
    s.output.write("\n", "text")
    s.r = {f"{k}": matrix(C, labs, [f"dim_{j + 1}" for j in range(ctx["f"])])
           for k, (_, labs, _, C) in zip(ctx["keys"], ctx["tables"])}


# ---------------------------------------------------------------------------
# ca
# ---------------------------------------------------------------------------

@command("ca", "ca")
def cmd_ca(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"dimensions": 3, "normalize": 4}) if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if len(names) != 2:
        raise StataError(102 if len(names) < 2 else 103, "too few variables specified" if len(names) < 2 else "too many variables specified")
    codes, labels, w, wtype, mask = _cat_sample(s, p, names)
    R, Cn = len(labels[0]), len(labels[1])
    if min(R, Cn) < 2:
        raise StataError(2001, "insufficient observations")                               # VERIFY
    N = np.zeros((R, Cn))
    np.add.at(N, (codes[0], codes[1]), 1.0 if w is None else w)
    norm = str(o.get("normalize", "symmetric")).strip()
    if norm not in ("symmetric", "standard", "row", "column", "principal"):
        raise StataError(198, "invalid normalize()")                                      # VERIFY
    fit = O.ca_fit(N)
    f = _dims(o, len(fit["sv"]))
    Rc, Cc = O.ca_coordinates(fit, norm)
    rs = O.ca_point_stats(fit["r"], fit["A"], fit["sv"], f)
    cs = O.ca_point_stats(fit["c"], fit["B"], fit["sv"], f)
    n, total = fit["n"], fit["total"]
    df = (R - 1) * (Cn - 1)
    chi2 = n * total
    out = s.output
    _header(out, "Correspondence analysis", [
        ("", _right("Number of obs", gnum(n, "%10.0g"))), ("", _right(f"Pearson chi2({df})", f"{chi2:.2f}")),
        ("", _right("Prob > chi2", f"{chi2_p(chi2, df):.4f}")), ("", _right("Total inertia", f"{total:.4f}"))])
    out.write(f"\n{R:>12} active rows{'':>25}{Cn:>6} active columns\n{0:>12} supplementary rows{'':>18}{0:>6} supplementary columns\n"
              f"{f:>12} dimensions\n", "text")                                            # VERIFY
    out.write(f"\n{'':>14}{'singular':>10}{'principal':>12}{'':>9}{'':>9}{'cumulative':>12}\n"
              f"{'Dimension':>14}{'value':>10}{'inertia':>12}{'chi2':>9}{'percent':>9}{'percent':>12}\n", "text")
    rule = "  " + "-" * 66
    out.write(rule + "\n", "text")
    cum = 0.0
    for j in range(len(fit["sv"])):
        pct = fit["inertia"][j] / total
        cum += pct
        out.write(f"{'dim ' + str(j + 1):>14}{fit['sv'][j]:>10.4f}{fit['inertia'][j]:>12.4f}{n * fit['inertia'][j]:>9.2f}"
                  f"{100 * pct:>9.2f}{100 * cum:>12.2f}\n", "text")
    out.write(f"{rule}\n{'total':>14}{'':>10}{total:>12.4f}{chi2:>9.2f}{100.0:>9.2f}\n", "text")
    tables = [(f"Statistics for row categories in {norm} normalization", [str(x) for x in labels[0]], rs, Rc[:, :f]),
              (f"Statistics for column categories in {norm} normalization", [str(x) for x in labels[1]], cs, Cc[:, :f])]
    for title, labs, st, C in tables:
        _point_table(out, title, labs, st, C, f)
    out.write("\n", "text")
    ctx = {"tables": tables, "f": f, "keys": ["row", "col"]}
    scalars = [("N", n), ("f", f), ("inertia", total), ("pchi2", chi2), ("df", df)]
    macros = [("cmdline", "ca " + args.strip()), ("cmd", "ca"), ("rowvar", names[0]), ("colvar", names[1]),
              ("normalize", norm), ("title", "Correspondence analysis"), ("estat_cmd", "ca_estat"), ("properties", "nob noV"),
              ("wtype", wtype)]
    mats = [("sv", matrix(fit["sv"].reshape(1, -1), ["sv"], [f"dim{j + 1}" for j in range(len(fit["sv"]))])),
            ("Ev", matrix(fit["inertia"].reshape(1, -1), ["Ev"], [f"dim{j + 1}" for j in range(len(fit["sv"]))])),
            ("rcoord", matrix(Rc[:, :f], labels[0], [f"dim{j + 1}" for j in range(f)])),
            ("ccoord", matrix(Cc[:, :f], labels[1], [f"dim{j + 1}" for j in range(f)]))]
    post_e(s, "ca", scalars, macros, mats, mask, int(round(n)),
           estat={"table": (3, _estat_table), "coordinates": (5, _estat_coordinates)}, extra={"ord": ctx})


# ---------------------------------------------------------------------------
# mca
# ---------------------------------------------------------------------------

@command("mca", "mca")
def cmd_mca(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"dimensions": 3, "normalize": 4}) if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if len(names) < 2:
        raise StataError(102, "too few variables specified")
    codes, labels, w, wtype, mask = _cat_sample(s, p, names)
    if w is not None:
        raise StataError(101, "weights not allowed")                                      # VERIFY
    norm = str(o.get("normalize", "adjusted")).strip()
    if norm not in ("adjusted", "principal"):
        raise StataError(198, "invalid normalize()")                                      # VERIFY
    fit = O.mca_fit(codes, [len(x) for x in labels])
    f = _dims(o, len(fit["sv"]))
    n = len(codes[0])
    if norm == "adjusted":
        pin, total = fit["adjusted"], fit["adj_total"]
        shown = max(int(np.sum(pin > 0)), 1)
        coord = fit["B"] * np.sqrt(pin)
    else:
        pin, total = fit["inertia"], fit["total"]
        shown = len(pin)
        coord = fit["B"] * fit["sv"]
    cats = [f"{nm}:{lb}" for nm, labs in zip(names, labels) for lb in labs]
    st = O.ca_point_stats(fit["c"], fit["B"], fit["sv"], f)
    out = s.output
    _header(out, "Multiple correspondence analysis", [
        ("", _right("Number of obs", gnum(n, "%10.0g"))), ("", _right("Total inertia", f"{total:.4f}")),
        ("", _right("Number of axes", str(f)))])                                            # VERIFY
    out.write(f"\nMethod: Burt/{norm} inertias\n", "text")
    out.write(f"\n{'':>14}{'principal':>12}{'cumul':>10}\n{'Dimension':>14}{'inertia':>12}{'percent':>10}\n", "text")
    rule = "  " + "-" * 34
    out.write(rule + "\n", "text")
    cum = 0.0
    for j in range(shown):
        cum += pin[j] / total
        out.write(f"{'dim ' + str(j + 1):>14}{pin[j]:>12.4f}{100 * cum:>10.2f}\n", "text")
    out.write(f"{rule}\n{'total':>14}{total:>12.4f}\n", "text")
    title = f"Statistics for column categories in {norm} normalization"
    _point_table(out, title, cats, st, coord[:, :f], f)
    out.write("\n", "text")
    ctx = {"tables": [(title, cats, st, coord[:, :f])], "f": f, "keys": ["coord"]}
    scalars = [("N", n), ("f", f), ("inertia", total), ("k_active", len(names))]
    macros = [("cmdline", "mca " + args.strip()), ("cmd", "mca"), ("varlist", " ".join(names)), ("normalize", norm),
              ("title", "Multiple correspondence analysis"), ("estat_cmd", "ca_estat"), ("properties", "nob noV")]
    mats = [("Ev", matrix(pin[:shown].reshape(1, -1), ["Ev"], [f"dim{j + 1}" for j in range(shown)])),
            ("coord", matrix(coord[:, :f], cats, [f"dim{j + 1}" for j in range(f)]))]
    post_e(s, "mca", scalars, macros, mats, mask, n,
           estat={"table": (3, _estat_table), "coordinates": (5, _estat_coordinates)}, extra={"ord": ctx})


# ---------------------------------------------------------------------------
# mds, mdsmat
# ---------------------------------------------------------------------------

def _mds_run(s: "Session", D: np.ndarray, ids: list[str], o: dict, desc: str, cmdline: str, cmd: str, mask: np.ndarray) -> None:
    n = D.shape[0]
    method = str(o.get("method", "classical")).strip()
    if method not in ("classical", "modern"):
        raise StataError(198, f"method {method} is not available yet in OpenDTA")        # VERIFY
    dims = int(float(o["dimension"])) if isinstance(o.get("dimension"), str) else 2
    cl = O.mds_classical(D, dims)
    dims = cl["dims"]
    out = s.output
    ev = cl["ev"]
    Y = cl["Y"]
    scalars = [("N", n), ("p", dims)]
    mats = []
    if method == "classical":
        npos = int(np.sum(ev > 0))
        m1, m2 = O.mardia(ev, dims)
        a, sq = np.abs(ev).sum(), (ev ** 2).sum()
        out.write(f"\nClassical metric multidimensional scaling\n    dissimilarity: {desc}\n\n", "text")
        out.write(f"{'Eigenvalues > 0':<45}{_right('Number of obs', gnum(n, '%10.0g'))}\n", "text")
        out.write(f"{'Retained dimensions = ' + str(dims):<45}\n", "text")
        out.write(f"\n   {'-' * 70}\n{'':>24}|{'abs(eigenvalue)':>26}{'(eigenvalue)^2':>20}\n"
                  f"{'Dimension':>24}|{'Eigenvalue':>12}{'Percent':>9}{'Cumul.':>9}{'Percent':>9}{'Cumul.':>9}\n"
                  f"   {'-' * 21}+{'-' * 48}\n", "text")                                  # VERIFY
        c1 = c2 = 0.0
        for j in range(dims):
            c1 += 100 * abs(ev[j]) / a
            c2 += 100 * ev[j] ** 2 / sq
            out.write(f"{j + 1:>24}|{ev[j]:>12.5g}{100 * abs(ev[j]) / a:>9.2f}{c1:>9.2f}{100 * ev[j] ** 2 / sq:>9.2f}{c2:>9.2f}\n", "text")
        out.write(f"   {'-' * 70}\n   Mardia fit measure 1 = {m1:.4f}\n   Mardia fit measure 2 = {m2:.4f}\n\n", "text")
        scalars += [("np", npos), ("trace", float(ev[ev > 0].sum())), ("mardia1", m1), ("mardia2", m2)]
        mats.append(("Ev", matrix(ev.reshape(1, -1), ["Ev"], [f"c{j + 1}" for j in range(n)])))
    else:
        Y, iters = O.mds_modern(D, Y)
        U, sv, Vt = np.linalg.svd(Y, full_matrices=False)
        Y = Y @ Vt.T
        st = O.stress1(D, Y)
        out.write(f"\nModern multidimensional scaling\n    dissimilarity: {desc}\n    Loss criterion: stress = raw_stress/norm(D)\n"
                  f"    Transformation: identity (no transformation)\n\n", "text")
        out.write(f"{'':>45}{_right('Number of obs', gnum(n, '%10.0g'))}\n{'':>45}{_right('Dimensions', str(dims))}\n"
                  f"{'':>45}{_right('Loss criterion', f'{st:.4f}')}\n\n", "text")           # VERIFY
        scalars += [("stress", st), ("iter", iters)]
    cn = [f"dim{j + 1}" for j in range(dims)]
    mats += [("Y", matrix(Y, ids, cn)), ("D", matrix(D, ids, ids))]
    macros = [("cmdline", f"{cmd} " + cmdline.strip()), ("cmd", cmd), ("method", method), ("title", f"{method.capitalize()} multidimensional scaling"),
              ("estat_cmd", "mds_estat"), ("predict", "mds_p"), ("properties", "nob noV")]
    post_e(s, cmd, scalars, macros, mats, mask, n, estat={"config": (3, _estat_config)},
           predict_multi={"": lambda ss, ee: _mds_predict(ss, ee), "config": lambda ss, ee: _mds_predict(ss, ee)},
           predict_opts={"config": 3}, extra={"mds": {"Y": Y, "ids": ids, "mask": mask}})


@command("mds", "mds")
def cmd_mds(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"id": 2, "dimension": 3, "measure": 3, "method": 3, "noplot": 3}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 1)
    X, _, _, mask = sample(s, p, names, allowed=())
    meas, power = _measure(o.get("measure"), "L2")
    D = K.pairwise(X, None, meas, power)
    if isinstance(o.get("id"), str):
        iv = s.data.get(s.expand_varlist(o["id"])[0])
        ids = [str(x) if iv.is_string else gnum(x) for x in (iv.raw if iv.is_string else iv.data)[mask]]
    else:
        ids = [str(i + 1) for i in np.flatnonzero(mask)]
    desc = f"{meas}, computed on {len(names)} variable{'s' if len(names) != 1 else ''}"
    _mds_run(s, D, ids, o, desc, args, "mds", mask)


@command("mdsmat", "mdsmat")
def cmd_mdsmat(s: "Session", args: str) -> None:
    from .matrix import get_matrix
    p = parse_standard(args)
    o = match_options(p.options, {"dimension": 3, "method": 3, "noplot": 3, "names": 5}) if p.options.strip() else {}
    mat = get_matrix(s, p.varlist.strip())
    if mat.rows != mat.cols or not np.allclose(mat.data, mat.data.T):
        raise StataError(505, "matrix not symmetric")
    _mds_run(s, mat.data, list(mat.colnames), o, "matrix " + p.varlist.strip(), args, "mdsmat", np.zeros(s.data.nobs, dtype=bool))


def _estat_config(s, est, rest):
    d = est.extra["mds"]
    Y = d["Y"]
    out = s.output
    out.write(f"\nConfiguration\n\n{'':>15} |" + "".join(f"{'dim' + str(j + 1):>10}" for j in range(Y.shape[1])) + "\n", "text")
    out.write("-" * 16 + "+" + "-" * (10 * Y.shape[1] + 1) + "\n", "text")
    for i, idn in enumerate(d["ids"]):
        out.write(f"{idn[:15]:>15} |" + "".join(f"{Y[i, j]:>10.4f}" for j in range(Y.shape[1])) + "\n", "text")
    out.write("\n", "text")


def _mds_predict(s, est):
    d = est.extra["mds"]
    Y, mask = d["Y"], d["mask"]
    if len(mask) != s.data.nobs or not mask.any():
        raise StataError(198, "the data have changed since mds")                           # VERIFY

    def fn(m):
        cols = []
        for j in range(Y.shape[1]):
            v = np.full(len(mask), SYS)
            v[mask] = Y[:, j]
            cols.append(v)
        return cols
    return Y.shape[1], fn


# ---------------------------------------------------------------------------
# procrustes
# ---------------------------------------------------------------------------

@command("procrustes", "procrustes")
def cmd_procrustes(s: "Session", args: str) -> None:
    m = re.match(r"\s*\(([^)]*)\)\s*\(([^)]*)\)(.*)$", args, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    tgt, src = s.expand_varlist(m.group(1)), s.expand_varlist(m.group(2))
    p = parse_standard(m.group(3))
    o = match_options(p.options, {"transform": 5, "noscale": 3}) if p.options.strip() else {}
    if str(o.get("transform", "orthogonal")).strip() != "orthogonal":
        raise StataError(198, "only transform(orthogonal) is available in OpenDTA")      # VERIFY
    if len(tgt) != len(src):
        raise StataError(198, "the two variable lists must have the same number of variables")   # VERIFY
    for nm in tgt + src:
        if s.data.get(nm).is_string:
            raise StataError(109, "type mismatch")
    X, _, _, mask = sample(s, p, tgt + src, allowed=())
    k = len(tgt)
    Yt, Xs = X[:, :k], X[:, k:]
    r = O.procrustes_orthogonal(Yt, Xs, scale=not o.get("noscale"))
    n = X.shape[0]
    out = s.output
    _header(out, "Procrustes analysis: orthogonal transformation", [
        ("", _right("Number of obs", gnum(n, "%10.0g"))), ("", _right("Dilation (c)", f"{r['c']:.4f}"))])   # VERIFY
    out.write(f"\n    Target: {' '.join(tgt)}\n    Source: {' '.join(src)}\n", "text")
    out.write(f"\n    Sum of squared differences (rss) = {r['rss']:.6g}\n    Total sum of squares (ss)        = {r['ss']:.6g}\n"
              f"    Procrustes statistic (rss/ss)    = {r['rss'] / r['ss'] if r['ss'] else float('nan'):.6g}\n", "text")
    from .mvfactor import _print_matrix
    _print_matrix(out, "Orthogonal rotation matrix", src, r["T"], tgt)
    out.write("\n", "text")
    fitted = r["fitted"]
    scalars = [("N", n), ("ss", r["ss"]), ("rss", r["rss"]), ("P", r["rss"] / r["ss"] if r["ss"] else SYS),
               ("rmse", float(np.sqrt(r["rss"] / (n * k)))), ("c", r["c"])]
    macros = [("cmdline", "procrustes " + args.strip()), ("cmd", "procrustes"), ("transform", "orthogonal"),
              ("title", "Procrustes analysis"), ("predict", "procrustes_p"), ("properties", "nob noV")]
    mats = [("A", matrix(r["T"], src, tgt)), ("trans", matrix(r["trans"].reshape(1, -1), ["trans"], tgt))]

    def pred(ss, ee, resid=False):
        ds = ss.data
        X2 = np.column_stack([ds.get(n_).data.astype(float) for n_ in src])
        ok = np.all(X2 < SYS, axis=1)
        F = r["c"] * X2 @ r["T"] + r["trans"]
        if resid:
            F = np.column_stack([ds.get(n_).data.astype(float) for n_ in tgt]) - F
        return k, (lambda m: [np.where(ok, F[:, j], SYS) for j in range(k)])
    post_e(s, "procrustes", scalars, macros, mats, mask, n,
           predict_multi={"": lambda ss, ee: pred(ss, ee), "fitted": lambda ss, ee: pred(ss, ee),
                          "residuals": lambda ss, ee: pred(ss, ee, True)}, predict_opts={"fitted": 3})


# ---------------------------------------------------------------------------
# biplot
# ---------------------------------------------------------------------------

@command("biplot", "biplot")
def cmd_biplot(s: "Session", args: str) -> None:
    from ..graphics.model import Graph, Layer, Panel
    from .graphcmd import finish
    p = parse_standard(args)
    o = match_options(p.options, {"std": 3, "dim": 3, "rowlabel": 4, "norow": 5, "nocolumn": 5, "name": 3, "saving": 3, "nodraw": 3}) \
        if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    X, _, _, mask = sample(s, p, names, allowed=())
    pair = [int(float(x)) - 1 for x in str(o["dim"]).replace(",", " ").split()] if isinstance(o.get("dim"), str) else [0, 1]
    if len(pair) != 2 or min(pair) < 0 or max(pair) >= len(names):
        raise StataError(198, "invalid dim()")                                            # VERIFY
    b = O.biplot_coords(X, std=bool(o.get("std")))
    a, c = b["rows"][:, pair], b["cols"][:, pair]
    layers = []
    if not o.get("norow"):
        lab = None
        if isinstance(o.get("rowlabel"), str):
            lv = s.data.get(s.expand_varlist(o["rowlabel"])[0])
            lab = [str(x) for x in (lv.raw if lv.is_string else lv.data)[mask]]
        layers.append(Layer("scatter", a[:, 0].tolist(), a[:, 1].tolist(), {"mlabel": lab} if lab else {}, {}, ""))
    if not o.get("nocolumn"):
        xs, ys = [], []
        for x_, y_ in c:
            xs += [0.0, float(x_), float("nan")]
            ys += [0.0, float(y_), float("nan")]
        layers.append(Layer("line", xs, ys, {}, {}, ""))
        layers.append(Layer("scatter", c[:, 0].tolist(), c[:, 1].tolist(), {"mlabel": list(names)}, {}, ""))
    tot = (b["sv"] ** 2).sum()
    opts = {"_xtitle": f"Dimension {pair[0] + 1} ({100 * b['sv'][pair[0]] ** 2 / tot:.1f}%)",
            "_ytitle": f"Dimension {pair[1] + 1} ({100 * b['sv'][pair[1]] ** 2 / tot:.1f}%)",
            "title": {"lines": ["Biplot"]}, "legend": {"off": True}}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), o)
