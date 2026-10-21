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

_PAD = " " * 4


def _center(text: str, width: int) -> str:
    left = (width - len(text)) // 2
    return " " * left + text + " " * (width - len(text) - left)


def _point_table(out, title: str, sections: list, labs: list[str], st: dict, C: np.ndarray, f: int) -> None:
    """Mass, quality and inertia share of every category and, for each retained dimension, its
    coordinate, squared correlation and contribution. `sections` lists (variable, first, last + 1)
    of the categories that belong together."""
    out.write(f"\n{title}\n\n", "text")
    out.write(f"{'':>17}|{_center('overall', 27)}|" + "|".join(_center(f"dimension_{j + 1}", 27) for j in range(f)) + "\n", "text")
    cells = lambda *h: "".join(f"{x:>8} " for x in h)
    out.write(f"{'Categories':>16} |{cells('mass', 'quality', '%inert')}|"
              + "|".join(cells("coord", "sqcorr", "contrib") for _ in range(f)) + "\n", "text")
    rule = f"{_PAD}{'-' * 13}+{'-' * 27}" + ("+" + "-" * 27) * f
    out.write(rule + "\n", "text")
    num = lambda *v: "".join(f"{x:>8.3f} " for x in v)
    for k_, (var, lo, hi) in enumerate(sections):
        if k_:
            out.write(rule + "\n", "text")
        out.write(f"{_PAD}{var[:12]:<12} |{'':>27}" + f"|{'':>27}" * f + "\n", "text")
        for i in range(lo, hi):
            out.write(f"{labs[i][:15]:>16} |{num(st['mass'][i], st['quality'][i], st['inertia'][i])}|"
                      + "|".join(num(C[i, j], st['sqcorr'][i, j], st['contrib'][i, j]) for j in range(f)) + "\n", "text")
    out.write(f"{_PAD}{'-' * (13 + 1 + 27 * (f + 1) + f)}\n", "text")


def _matlist(out, title: str, corner: str, cols: list[str], rows: list[str], M_, fmt, *, lower: bool = False,
             bottom: bool = True) -> None:
    """Matrix with a stub of 12 characters and cells of 10 followed by a space (estat tables)."""
    out.write(f"\n{title}\n\n{_PAD}{corner:>12} |" + "".join(f"{c:>10} " for c in cols) + "\n", "text")
    out.write(f"{_PAD}{'-' * 13}+{'-' * (11 * len(cols))}\n", "text")
    for i, r_ in enumerate(rows):
        m = i + 1 if lower else len(cols)
        out.write(f"{_PAD}{r_[:12]:>12} |" + "".join(f"{fmt(M_[i, j]):>10} " for j in range(m)) + "\n", "text")
    if bottom:
        out.write(f"{_PAD}{'-' * (13 + 1 + 11 * len(cols))}\n", "text")


def _estat_table(s, est, rest):
    ctx = est.extra["ord"]
    fit, f = ctx["fit"], ctx["f"]
    P = np.outer(fit["r"], fit["c"]) * (1.0 + (fit["A"][:, :f] * fit["sv"][:f]) @ fit["B"][:, :f].T)
    _matlist(s.output, f"Approximation for dim = {f} (normalized to overall sum = 1)", "",
             ctx["labels"][1], ctx["labels"][0], P, lambda x: f"{x:.4f}")
    s.output.write("\n", "text")


def _estat_coordinates(s, est, rest):
    ctx = est.extra["ord"]
    out = s.output
    f = ctx["f"]
    cols = [f"dim{j + 1}" for j in range(f)]
    if ctx["cmd"] == "ca":
        for side, var, labs, C in zip(("row", "column"), ctx["vars"], ctx["labels"], ctx["coords"]):
            _matlist(out, f"Coordinates in {side} space in {ctx['norm']} normalization", var, cols, labs, C, lambda x: f"{x:.4f}")
        s.r = {"coord_row": matrix(ctx["coords"][0], ctx["labels"][0], cols), "coord_col": matrix(ctx["coords"][1], ctx["labels"][1], cols)}
    else:
        C = ctx["coords"]
        w = 14
        out.write(f"\nColumn {ctx['norm']} coordinates\n\n", "text")
        out.write(f"{_PAD}+{'-' * (w + 1 + 11 * f)}+\n", "text")
        out.write(f"{_PAD}|{'Categories':>13} |" + "".join(f"{c:>10} " for c in cols) + "|\n", "text")
        out.write(f"{_PAD}|{'-' * w}+{'-' * (11 * f)}|\n", "text")
        for k_, (var, lo, hi) in enumerate(ctx["sections"]):
            if k_:
                out.write(f"{_PAD}|{'-' * w}+{'-' * (11 * f)}|\n", "text")
            out.write(f"{_PAD}|{var[:w]:<{w}}|{'':>{11 * f}}|\n", "text")
            for i in range(lo, hi):
                out.write(f"{_PAD}|{ctx['labels'][i][:13]:>13} |" + "".join(f"{gnum(C[i, j]):>10} " for j in range(f)) + "|\n", "text")
        out.write(f"{_PAD}+{'-' * (w + 1 + 11 * f)}+\n", "text")
        s.r = {"coord": matrix(C, ctx["labels"], cols)}
    out.write("\n", "text")


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
    nd = len(fit["sv"])
    expl = 100.0 * fit["inertia"][:f].sum() / total
    out = s.output
    blk = lambda lab, val: f"{lab:<19}={val:>11}"
    out.write(f"\n{'Correspondence analysis':<46}{blk('Number of obs', gnum(n, '%10.0g'))}\n", "text")
    out.write(f"{'':<46}{blk(f'Pearson chi2({df})', f'{chi2:.2f}')}\n{'':<46}{blk('Prob > chi2', f'{chi2_p(chi2, df):.4f}')}\n", "text")
    out.write(f"{'':<46}{blk('Total inertia', f'{total:.4f}')}\n", "text")
    out.write(f"{f'{R} active rows':<46}{blk('Number of dim.', str(f))}\n{f'{Cn} active columns':<46}{blk('Expl. inertia (%)', f'{expl:.2f}')}\n".replace("\n    ", "\n    ").replace(f"\n{R} active", f"\n{R} active"), "text") if False else None
    out.write(f"{f'    {R} active rows':<46}{blk('Number of dim.', str(f))}\n", "text")
    out.write(f"{f'    {Cn} active columns':<46}{blk('Expl. inertia (%)', f'{expl:.2f}')}\n", "text")
    out.write(f"\n{'':>16}|   singular    principal                             cumul\n", "text")
    out.write(f"{'Dimension':>15} |    value       inertia           chi2    percent   percent\n", "text")
    rule = f"{_PAD}{'-' * 12}+{'-' * 60}"
    out.write(rule + "\n", "text")
    cum = 0.0
    for j in range(nd):
        pct = 100.0 * fit["inertia"][j] / total
        cum += pct
        out.write(f"{'dim ' + str(j + 1):>15} |{gnum(fit['sv'][j]):>11}{gnum(fit['inertia'][j]):>12}{n * fit['inertia'][j]:>15.2f}"
                  f"{pct:>11.2f}{cum:>10.2f}\n", "text")
    out.write(f"{rule}\n{'total':>15} |{'':>11}{gnum(total):>12}{chi2:>15.2f}{'100':>8}\n", "text")
    labs = [str(x) for x in labels[0]] + [str(x) for x in labels[1]]
    st = {k: np.concatenate([rs[k], cs[k]]) for k in ("mass", "inertia", "quality", "sqcorr", "contrib")}
    C = np.vstack([Rc[:, :f], Cc[:, :f]])
    _point_table(out, f"Statistics for row and column categories in {norm} normalization",
                 [(names[0], 0, R), (names[1], R, R + Cn)], labs, st, C, f)
    out.write("\n", "text")
    ctx = {"cmd": "ca", "fit": fit, "f": f, "norm": norm, "vars": names, "labels": [[str(x) for x in labels[0]], [str(x) for x in labels[1]]],
           "coords": [Rc[:, :f], Cc[:, :f]]}
    scalars = [("N", n), ("inertia", total), ("f", f), ("pinertia", expl / 100.0), ("X2", chi2), ("X2_df", df), ("X2_p", chi2_p(chi2, df))]
    macros = [("cmdline", "ca " + args.strip()), ("ca_data", "variables"), ("varlist", " ".join(names)), ("cmd", "ca"),
              ("marginsnotok", "_ALL"), ("estat_cmd", "ca_estat"), ("predict", "ca_p"), ("title", "correspondence analysis"),
              ("Cname", names[1]), ("Rname", names[0]), ("norm", norm), ("sv_unique", "1"), ("properties", "nob noV eigen")]
    dn = [f"dim{j + 1}" for j in range(nd)]
    gs = lambda stt, Cx: np.column_stack([stt["mass"], stt["quality"], stt["inertia"]]
                                         + [np.column_stack([Cx[:, j], stt["sqcorr"][:, j], stt["contrib"][:, j]]) for j in range(f)])
    gcols = ["mass", "quality", "inertia"] + [f"{k_}{j + 1}" for j in range(f) for k_ in ("coord", "sqcorr", "contrib")]
    P = N / n
    mats = [("Ccoding", matrix(np.arange(1, Cn + 1, dtype=float).reshape(1, -1), ["Ccoding"], [str(x) for x in labels[1]])),
            ("Rcoding", matrix(np.arange(1, R + 1, dtype=float).reshape(1, -1), ["Rcoding"], [str(x) for x in labels[0]])),
            ("GSC", matrix(gs(cs, Cc[:, :f]), [str(x) for x in labels[1]], gcols)),
            ("GSR", matrix(gs(rs, Rc[:, :f]), [str(x) for x in labels[0]], gcols)),
            ("TC", matrix(fit["B"], [str(x) for x in labels[1]], dn)), ("TR", matrix(fit["A"], [str(x) for x in labels[0]], dn)),
            ("C", matrix(Cc, [str(x) for x in labels[1]], dn)), ("R", matrix(Rc, [str(x) for x in labels[0]], dn)),
            ("Sv", matrix(fit["sv"].reshape(1, -1), ["Sv"], dn)),
            ("c", matrix(fit["c"].reshape(-1, 1), [str(x) for x in labels[1]], ["c"])),
            ("r", matrix(fit["r"].reshape(-1, 1), [str(x) for x in labels[0]], ["r"])),
            ("P", matrix(P, [str(x) for x in labels[0]], [str(x) for x in labels[1]]))]       # VERIFY: contents
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
    norm = str(o.get("normalize", "standard")).strip()
    if norm not in ("standard", "principal"):
        raise StataError(198, "invalid normalize()")                                      # VERIFY
    fit = O.mca_fit(codes, [len(x) for x in labels])
    f = _dims(o, len(fit["sv"]))
    n = len(codes[0])
    pin, total = fit["adjusted"], fit["adj_total"]
    shown = max(int(np.sum(pin > 0)), 1)
    coord = fit["B"] if norm == "standard" else fit["B"] * np.sqrt(pin)          # VERIFY: principal coordinates
    cats = [lb for labs in labels for lb in labs]
    bounds = np.cumsum([0] + [len(x) for x in labels])
    sections = [(nm, int(bounds[i]), int(bounds[i + 1])) for i, nm in enumerate(names)]
    st = O.mca_point_stats(fit, f)
    out = s.output
    blk = lambda lab, val: f"{lab:<18}={val:>11}"
    out.write(f"\n{'Multiple/Joint correspondence analysis':<47}{blk('Number of obs', gnum(n, '%10.0g'))}\n", "text")
    out.write(f"{'':<47}{blk('Total inertia', gnum(total, '%10.0g'))}\n", "text")
    out.write(f"{'    Method: Burt/adjusted inertias':<47}{blk('Number of axes', str(f))}\n", "text")
    out.write(f"\n{'':>16}|   principal               cumul\n{'Dimension':>15} |    inertia     percent   percent\n", "text")
    rule = f"{_PAD}{'-' * 12}+{'-' * 34}"
    out.write(rule + "\n", "text")
    cum = 0.0
    for j in range(shown):
        cum += 100.0 * pin[j] / total
        out.write(f"{'dim ' + str(j + 1):>15} |{gnum(pin[j]):>12}{100.0 * pin[j] / total:>10.2f}{cum:>11.2f}\n", "text")
    out.write(f"{rule}\n{'Total':>15} |{gnum(total):>12}{100.0:>10.2f}\n", "text")
    title = f"Statistics for column categories in {norm} normalization"
    _point_table(out, title, sections, cats, st, coord[:, :f], f)
    out.write("\n", "text")
    ctx = {"cmd": "mca", "f": f, "norm": norm, "labels": cats, "coords": coord[:, :f], "sections": sections}
    scalars = [("f", f), ("ev_unique", 1), ("inertia", total), ("adjust", 1), ("N", n)]
    macros = [("cmdline", "mca " + args.strip()), ("cmd", "mca"), ("marginsnotok", "_ALL"),
              ("title", "Multiple/Joint correspondence analysis"), ("method", "burt"), ("estat_cmd", "mca_estat"),
              ("predict", "mca_p"), ("crossed", "0"), ("defs", " / ".join(names)), ("names", " ".join(names)),
              ("norm", norm), ("properties", "nob noV eigen")]
    dn = [f"dim{j + 1}" for j in range(f)]
    J = len(cats)
    gen = np.column_stack([st["mass"], st["quality"], st["inertia"]]
                          + [np.column_stack([coord[:, j], st["sqcorr"][:, j], st["contrib"][:, j]]) for j in range(f)])
    gcols = ["mass", "quality", "inertia"] + [f"{k_}{j + 1}" for j in range(f) for k_ in ("coord", "sqcorr", "contrib")]
    mats = [(f"Coding{i + 1}", matrix(np.arange(1, len(labels[i]) + 1, dtype=float).reshape(-1, 1), labels[i], ["Coding"]))
            for i in reversed(range(len(names)))]
    mats += [("inertia_e", matrix(pin[:shown].reshape(1, -1), ["inertia_e"], [f"dim{j + 1}" for j in range(shown)])),
             ("Ev", matrix(pin[:shown].reshape(1, -1), ["Ev"], [f"dim{j + 1}" for j in range(shown)])),
             ("rSCW", matrix(fit["B"][:, :f] * np.sqrt(fit["c"])[:, None], cats, dn)),
             ("cGS", matrix(gen, cats, gcols)), ("F", matrix(fit["B"][:, :f] * fit["sv"][:f], cats, dn)),
             ("A", matrix(fit["B"][:, :f], cats, dn)),
             ("cInertia", matrix((fit["c"] * (fit["B"] * fit["sv"]) .__pow__(2).sum(1)).reshape(-1, 1), cats, ["cInertia"])),
             ("cDist", matrix(((fit["B"] * fit["sv"]) ** 2).sum(1).reshape(-1, 1), cats, ["cDist"])),
             ("cMass", matrix(fit["c"].reshape(-1, 1), cats, ["cMass"]))]               # VERIFY: contents
    post_e(s, "mca", scalars, macros, mats, mask, n,
           estat={"coordinates": (5, _estat_coordinates)}, extra={"ord": ctx})


# ---------------------------------------------------------------------------
# mds, mdsmat
# ---------------------------------------------------------------------------

def _mds_run(s: "Session", D: np.ndarray, ids: list[str], o: dict, desc: str, cmdline: str, cmd: str, mask: np.ndarray,
             meta: dict | None = None) -> None:
    meta = meta or {}
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
    if method == "classical":
        npos = int(np.sum(ev > 1e-12 * ev[0]))                           # VERIFY: tolerance
        m1, m2 = O.mardia(ev, dims)
        a, sq = np.abs(ev).sum(), (ev ** 2).sum()
        out.write(f"\nClassical metric multidimensional scaling\n    {desc}\n\n", "text")
        out.write(f"{'':>45}{'Number of obs':<21}={gnum(n, '%10.0g'):>11}\n", "text")
        out.write(f"    {'Eigenvalues > 0':<21}={npos:>10}         {'Mardia fit measure 1':<21}={m1:>11.4f}\n", "text")
        out.write(f"    {'Retained dimensions':<21}={dims:>10}         {'Mardia fit measure 2':<21}={m2:>11.4f}\n", "text")
        rule = "-" * 74
        out.write(f"\n    {rule}\n{'':>17}|{'abs(eigenvalue)':>34}{'(eigenvalue)^2':>24}\n"
                  f"{'Dimension':>15}  |{'Eigenvalue':>12}{'Percent':>13}{'Cumul.':>10}{'Percent':>14}{'Cumul.':>10}\n"
                  f"    {'-' * 13}+{'-' * 60}\n", "text")                               # VERIFY
        c1 = c2 = 0.0
        for j in range(npos):
            if j == dims and j > 0:
                out.write(f"    {'-' * 13}+{'-' * 60}\n", "text")
            c1 += 100 * abs(ev[j]) / a
            c2 += 100 * ev[j] ** 2 / sq
            out.write(f"{j + 1:>15}  |{format(ev[j], '.8g'):>12}{100 * abs(ev[j]) / a:>13.2f}{c1:>10.2f}"
                      f"{100 * ev[j] ** 2 / sq:>14.2f}{c2:>10.2f}\n", "text")
        out.write(f"    {rule}\n", "text")
        scalars = [("mardia2", m2), ("mardia1", m1), ("np", npos), ("p", dims), ("addcons", 0), ("rc", 0), ("N", n)]
        mats = []
    else:
        out.write("(loss(stress) assumed)\n(transform(identity) assumed)\n\n", "text")
        trace: list[float] = []
        Y, iters = O.mds_modern(D, Y, trace=trace)
        U, sv, Vt = np.linalg.svd(Y, full_matrices=False)
        Y = Y @ Vt.T
        st = O.stress1(D, Y)
        for k_, v_ in enumerate(trace, 1):
            out.write(f"Iteration {k_}:   stress = {gnum(v_, '%10.0g'):>10}\n", "text")
        out.write(f"\nModern multidimensional scaling\n    {desc}\n\n"
                  f"    Loss criterion: stress = raw_stress/norm(distances)\n"
                  f"    Transformation: identity (no transformation)\n\n", "text")
        out.write(f"{'':>49}{'Number of obs':<17}={gnum(n, '%10.0g'):>11}\n{'':>49}{'Dimensions':<17}={dims:>11}\n"
                  f"{'    Normalization: principal':<49}{'Loss criterion':<17}={st:>11.4f}\n", "text")      # VERIFY
        scalars = [("stress", st), ("iter", iters), ("p", dims), ("N", n)]
        mats = []
    cn = [f"dim{j + 1}" for j in range(dims)]
    macros = [("cmdline", f"{cmd} " + cmdline.strip()), ("estat_cmd", "mds_estat"), ("cmd", cmd), ("marginsnotok", "_ALL"),
              ("predict", "mds_p"), ("dtype", meta.get("dtype", "dissimilarity")), ("dname", meta.get("dname", "")),
              ("duplicates", "0")]
    if cmd == "mds":
        macros += [("idtype", meta.get("idtype", "int")), ("id", meta.get("idname", "")), ("varlist", meta.get("varlist", ""))]
    macros += [("method", method), ("unique", "1"), ("norm", "principal"),
               ("properties", "nob noV eigen" if method == "classical" else "nob noV")]
    if cmd == "mds":
        mats += [("idcoding", matrix(np.arange(1.0, n + 1).reshape(-1, 1), ids, ["idcoding"])),
                 ("coding", matrix(np.tile([0.0, 1.0], (len(meta.get("names", [])), 1)), meta.get("names", []), ["shift", "scale"]))]
    mats += [("D", matrix(D, ids, ids)), ("Y", matrix(Y, ids, cn)),
             ("linearf", matrix(np.array([[0.0, 1.0]]), ["linearf"], ["c1", "c2"]))]                  # VERIFY
    if method == "classical":
        mats.append(("Ev", matrix(ev.reshape(1, -1), ["Ev"], [f"c{j + 1}" for j in range(n)])))
    post_e(s, cmd, scalars, macros, mats, mask, n, estat={"config": (3, _estat_config)},
           predict_multi={"": lambda ss, ee: _mds_predict(ss, ee), "config": lambda ss, ee: _mds_predict(ss, ee)},
           predict_opts={"config": 3}, extra={"mds": {"Y": Y, "ids": ids, "mask": mask, "idname": meta.get("idname", "")}})


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
    desc = f"dissimilarity: {meas}, computed on {len(names)} variable{'s' if len(names) != 1 else ''}"
    meta = {"dname": meas, "varlist": " ".join(names), "names": names, "idname": s.expand_varlist(o["id"])[0] if isinstance(o.get("id"), str) else ""}
    if isinstance(o.get("id"), str) and iv.is_string:
        meta["idtype"] = "str"
    _mds_run(s, D, ids, o, desc, args, "mds", mask, meta)


@command("mdsmat", "mdsmat")
def cmd_mdsmat(s: "Session", args: str) -> None:
    from .matrix import get_matrix
    p = parse_standard(args)
    o = match_options(p.options, {"dimension": 3, "method": 3, "noplot": 3, "names": 5}) if p.options.strip() else {}
    mat = get_matrix(s, p.varlist.strip())
    if mat.rows != mat.cols or not np.allclose(mat.data, mat.data.T):
        raise StataError(505, "matrix not symmetric")
    _mds_run(s, mat.data, list(mat.colnames), o, "dissimilarity matrix: " + p.varlist.strip(), args, "mdsmat", np.zeros(s.data.nobs, dtype=bool),
             {"dname": p.varlist.strip()})


def _estat_config(s, est, rest):
    d = est.extra["mds"]
    Y = d["Y"]
    k = Y.shape[1]
    out = s.output
    out.write(f"\nConfiguration in {k}-dimensional Euclidean space (principal normalization)\n\n"
              f"{d['idname'] or 'obs':>9} |" + "".join(f"{'dim' + str(j + 1):>13} " for j in range(k)) + "\n", "text")     # VERIFY
    out.write(f"    {'-' * 6}+{'-' * (14 * k)}\n", "text")
    for i, idn in enumerate(d["ids"]):
        out.write(f"{idn[:9]:>9} |" + "".join(f"{Y[i, j]:>13.4f} " for j in range(k)) + "\n", "text")
    out.write(f"    {'-' * (7 + 14 * k)}\n", "text")


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

def _g7(x: float) -> str:
    """Seven significant digits, no leading zero: how Stata's matrix tables show a number."""
    if not np.isfinite(x):
        return "."
    r = f"{x:.7g}"
    return r.replace("0.", ".", 1) if r.startswith("0.") else r.replace("-0.", "-.", 1)


def _proc_matrix(out, title: str, corner: str, rows: list[str], cols: list[str], M_) -> None:
    w = 13 + 1 + 11 * len(cols)
    out.write(f"{title}\n\n    {'-' * w}\n{corner:>16} |" + "".join(f"{c_:>10} " for c_ in cols) + "\n"
              f"    {'-' * 13}+{'-' * (11 * len(cols))}\n", "text")
    for i, r_ in enumerate(rows):
        out.write(f"{r_:>16} |" + "".join(f"{_g7(M_[i, j]):>10} " for j in range(len(cols))) + "\n", "text")
    out.write(f"    {'-' * w}\n\n", "text")


@command("procrustes", "procrustes")
def cmd_procrustes(s: "Session", args: str) -> None:
    m = re.match(r"\s*\(([^)]*)\)\s*\(([^)]*)\)(.*)$", args, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    tgt, src = s.expand_varlist(m.group(1)), s.expand_varlist(m.group(2))
    p = parse_standard(m.group(3))
    o = match_options(p.options, {"transform": 5, "norho": 3, "noscale": 3}) if p.options.strip() else {}
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
    r = O.procrustes_orthogonal(Yt, Xs, scale=not o.get("norho"))          # VERIFY: noscale is accepted and ignored
    n = X.shape[0]
    out = s.output
    df_m = k * k - k * (k + 1) // 2 + k + (0 if o.get("norho") else 1)
    df_r = n * k - df_m
    P = r["rss"] / r["ss"] if r["ss"] else float("nan")
    rmse = float(np.sqrt(r["rss"] / df_r)) if df_r > 0 else float("nan")
    urmse = float(np.sqrt(r["rss"] / (n * k)))
    lab = " " * 43

    def row(label, val):
        out.write(f"{lab}{label:<23}={val:>11}\n", "text")
    out.write(f"{'Procrustes analysis (orthogonal)':<43}{'Number of observations':<23}={n:>11}\n", "text")
    row("Model df (df_m)", df_m)
    row("Residual df (df_r)", df_r)
    row("SS(target)", f"{r['ss']:.4f}")
    row("RSS(target)", f"{r['rss']:.4f}")
    row("RMSE = root(RSS/df_r)", f"{rmse:.7g}")
    row("Procrustes = RSS/SS", f"{P:.4f}")
    _proc_matrix(out, "Translation c", "", ["_cons"], tgt, r["trans"].reshape(1, -1))
    _proc_matrix(out, "Rotation & reflection matrix A (orthogonal)", "", src, tgt, r["T"])
    out.write(f"Dilation factor\n\n{'rho =':>18}{r['c']:>10.4f}\n\n", "text")
    ddf = df_r / k if k else float("nan")
    ys = np.vstack([r["ss_i"], r["rss_i"], np.sqrt(r["rss_i"] / ddf), r["rss_i"] / r["ss_i"], r["corr_i"]])
    _proc_matrix(out, "Fit statistics by target variable", "Statistics", ["SS", "RSS", "RMSE", "Procrustes", "Corr_y_yhat"], tgt, ys)
    scalars = [("rho", r["c"]), ("N", n), ("ss", r["ss"]), ("rss", r["rss"]), ("P", P), ("rmse", rmse), ("urmse", urmse),
               ("df_m", df_m), ("df_r", df_r), ("ny", k)]
    macros = [("cmdline", "procrustes " + args.strip()), ("cmd", "procrustes"), ("marginsnotok", "_ALL"),
              ("estat_cmd", "procrustes_estat"), ("predict", "procrustes_p"), ("uniqueA", "1"), ("transform", "orthogonal"),
              ("xlist", " ".join(src)), ("ylist", " ".join(tgt)), ("properties", "nob noV")]
    mats = [("ystats", matrix(ys, ["SS", "RSS", "RMSE", "Procrustes", "Corr_y_yhat"], tgt)),
            ("c", matrix(r["trans"].reshape(1, -1), ["c"], tgt)), ("A", matrix(r["T"], src, tgt))]

    def pred(ss, ee, resid=False, note=False):
        ds = ss.data
        X2 = np.column_stack([ds.get(n_).data.astype(float) for n_ in src])
        ok = np.all(X2 < SYS, axis=1)
        F = r["c"] * X2 @ r["T"] + r["trans"]
        if resid:
            F = np.column_stack([ds.get(n_).data.astype(float) for n_ in tgt]) - F

        def fn(m):
            if note:
                ss.output.write("(fitted assumed)\n", "text")
            return [np.where(ok, F[:, j], SYS) for j in range(k)]
        return k, fn
    post_e(s, "procrustes", scalars, macros, mats, mask, n,
           predict_multi={"": lambda ss, ee: pred(ss, ee, note=True), "fitted": lambda ss, ee: pred(ss, ee),
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
    e1, e2 = (b["sv"][pair[0]] ** 2) / tot, (b["sv"][pair[1]] ** 2) / tot
    s.output.write(f"\nBiplot of {X.shape[0]} observations and {len(names)} variables\n\n"
                   f"    Explained variance by component {pair[0] + 1}  {e1:.4f}\n"
                   f"    Explained variance by component {pair[1] + 1}  {e2:.4f}\n"
                   f"{'Total explained variance':>37}  {e1 + e2:.4f}\n", "text")             # VERIFY
    opts = {"_xtitle": f"Dimension {pair[0] + 1} ({100 * b['sv'][pair[0]] ** 2 / tot:.1f}%)",
            "_ytitle": f"Dimension {pair[1] + 1} ({100 * b['sv'][pair[1]] ** 2 / tot:.1f}%)",
            "title": {"lines": ["Biplot"]}, "legend": {"off": True}}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), o)
