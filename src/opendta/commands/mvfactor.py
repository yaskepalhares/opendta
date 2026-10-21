"""Principal components, factor analysis, rotation and scale reliability ([MV] pca, [MV] factor,
[MV] rotate, [MV] alpha, [MV] screeplot, [MV] scoreplot) with their estat subcommands and predict.

The statistics are in core/multivariate.py. Layouts follow the manuals' examples; what they do not
show carries a VERIFY mark and is settled by the reference log of compat case 1906."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import multivariate as mv
from ..core.dataset import Variable
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from .mvcommon import chi2_p, fill_e, g as gnum, lower_matrix, matrix, numeric_names, post_e, sample
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
RULE = "    "


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

def _header(out, title: str, lines: list[tuple[str, str]]) -> None:
    """Title and indented lines on the left, `label = value` on the right (column 50)."""
    out.write("\n", "text")
    for i, (left, right) in enumerate(lines):
        if i == 0:
            left = title
        out.write(f"{left:<49}{right}\n", "text")


def _right(label: str, value: str) -> str:
    return f"{label:<17}= {value:>10}"


def _g6(x: float) -> str:
    """Six significant digits without a leading zero (.716578), as in the eigenvalue tables."""
    if x == 0:
        return "0"
    t = f"{x:.6g}"
    if "e" in t:
        m, e = t.split("e")
        t = f"{m}e{e[0]}{e[1:].lstrip('0').rjust(2, '0')}"
    return t.replace("0.", ".", 1) if t.startswith(("0.", "-0.")) else t


def _g4(x: float) -> str:
    """Four significant digits without a leading zero (.3973), as in %9.4g."""
    if x == 0:
        return "0"
    t = f"{x:.4g}"
    return t.replace("0.", ".", 1) if t.startswith(("0.", "-0.")) else t


def _eig_table(out, name: str, ev, trace: float, head: str = "Eigenvalue", fixed: bool = False,
               oblique: str = "") -> None:
    """Eigenvalue (or variance) table. `fixed` prints the values with 5 decimals (factor); `oblique`
    drops the Difference column and writes the note in its place."""
    rule1 = RULE + "-" * 74
    fac = name == "Factor"
    stub = (lambda x: f"{x:>11}  |") if fac else (lambda x: f"{x:>12} |")
    out.write("\n" + rule1 + "\n", "text")
    if oblique:
        out.write(f"{RULE}{stub(name)}{head:>13}{'Proportion':>13}{oblique:>{34 if fac else 33}}\n", "text")
    else:
        out.write(f"{RULE}{stub(name)}{head:>13}{'Difference':>13}{'Proportion':>{18 if fac else 19}}{'Cumulative':>13}\n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 60}\n", "text")
    cum = 0.0
    p = len(ev)
    for j in range(p):
        diff = ev[j] - ev[j + 1] if j + 1 < p else None
        prop = ev[j] / trace if trace else float("nan")
        cum += prop
        f_ = (lambda v: f"{v:.5f}") if fixed else _g6
        if oblique:
            out.write(f"{RULE}{stub(_stub(name, j + 1))}{f_(ev[j]):>13}{prop:>13.4f}\n", "text")
            continue
        d = f"{(f_(diff) if diff is not None else '.'):>13}"
        out.write(f"{RULE}{stub(_stub(name, j + 1))}{f_(ev[j]):>13}{d}{prop:>{18 if fac else 19}.4f}{cum:>13.4f}\n", "text")
    out.write(rule1 + "\n", "text")


def _stub(name: str, j: int) -> str:
    return {"Component": f"Comp{j}", "Factor": f"Factor{j}"}.get(name, f"{name}{j}")


def _table(out, title: str, rows: list[str], cols: list[str], M_, fmt, extra=None, sub: str = "", stub: str = "Variable",
           blanks: float = 0.0, cw: int = 9, tri: bool = False, foot=None, note: str = "") -> None:
    """`title`, then a table whose stub column is headed `stub` and whose cells are fmt(value) right
    aligned in 9 characters followed by a space. `extra` = (heading, [row text], width): a column after a
    second bar, whose rows are pre-formatted."""
    k = len(cols)
    w = (cw + 1) * k
    top = RULE + "-" * (14 + w + (extra[2] + 1 if extra else 0))
    sep = f"{RULE}{'-' * 13}+{'-' * w}" + (f"+{'-' * extra[2]}" if extra else "")
    head = f"{RULE}{stub:>12} |" + "".join(f"{c:>{cw}} " for c in cols) + (f"|{extra[0]}" if extra else "")
    lines = []
    for i, nm in enumerate(rows):
        cells = "".join(f"{'':>{cw}} " if (tri and j > i) or abs(M_[i][j]) < blanks else f"{fmt(M_[i][j]):>{cw}} "
                        for j in range(k))
        lines.append(f"{RULE}{_abbr(nm):>12} |{cells}" + (f"|{extra[1][i]}" if extra else ""))
    if foot:
        lines.append(sep)
        for lab, vals in foot:
            lines.append(f"{RULE}{lab:>12} |" + "".join(f"{fmt(v):>{cw}} " for v in vals))
    out.write("\n" + title + "\n" + (sub + "\n" if sub else "") + "\n" + top + "\n" + head + "\n" + sep + "\n"
              + "\n".join(lines) + "\n" + top + "\n" + (note + "\n" if note else ""), "text")


def _f4(x: float) -> str:
    return f"{x:.4f}"


def _unexpl_extra(unexpl) -> tuple[str, list[str], int]:
    return (f"{'Unexplained':>12} ", [f"{(_g4(0.0 if abs(u) < 1e-12 else u)):>12} " for u in unexpl], 13)


def _uniq_extra(psi) -> tuple[str, list[str], int]:
    return (f"{'Uniqueness':>13} ", [f"{u:>12.4f}  " for u in psi], 14)


def _abbr(name: str, w: int = 12) -> str:
    return name if len(name) <= w else name[: w - 1] + "~"


def _print_matrix(out, title: str, names: list[str], M_: np.ndarray, cols: list[str] | None = None,
                  fmt=_f4) -> None:
    cols = cols or [_abbr(n, 8) for n in names]
    _table(out, title, names, cols, M_, fmt, stub="")


# ---------------------------------------------------------------------------
# pca
# ---------------------------------------------------------------------------

_PCA_OPTS = {"components": 4, "mineigen": 4, "correlation": 4, "covariance": 3, "vce": 3, "means": 4,
             "blanks": 5, "tol": 3, "ignore": 6, "norotated": 5, "level": 3}


@command("pca", "pca")
def cmd_pca(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _PCA_OPTS) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    if o.get("vce") and str(o["vce"]).strip() != "none":
        raise StataError(198, "option vce(normal) is not available yet in OpenDTA")
    if o.get("means"):
        raise StataError(198, "option means is not available yet in OpenDTA")
    X, w, wtype, mask = sample(s, p, names)
    means, S, R, sw = mv.weighted_moments(X, w)
    cov = bool(o.get("covariance"))
    if cov and o.get("correlation"):
        raise StataError(198, "options correlation and covariance may not be combined")    # VERIFY
    C = S if cov else R
    fit = mv.pca_fit(C)
    ev, V = fit["ev"], fit["V"]
    pdim = len(names)
    comp = int(float(o["components"])) if o.get("components") else None
    mine = float(o["mineigen"]) if o.get("mineigen") else None
    if comp is not None and (comp < 1 or comp > pdim):
        raise StataError(198, "components() must be between 1 and the number of variables")  # VERIFY
    # VERIFY: with neither option every component is kept (mineigen defaults to 1e-5 of the trace)
    k = comp if comp is not None else (int(np.sum(ev >= mine)) if mine is not None else int(np.sum(ev > 1e-5 * fit["trace"])))
    k = max(k, 1)
    blanks = float(o["blanks"]) if o.get("blanks") else 0.0
    N = int(round(sw))
    unexpl = mv.pca_unexplained(C, ev, V, k)
    rho = float(ev[:k].sum() / fit["trace"])
    out = s.output
    kind = "covariance" if cov else "correlation"
    _header(out, f"Principal components/{kind}", [
        ("", _right("Number of obs", f"{N}")),
        ("", _right("Number of comp.", f"{k}")),
        ("", _right("Trace", gnum(fit["trace"]))),
        ("    Rotation: (unrotated = principal)", _right("Rho", f"{rho:.4f}"))])
    _eig_table(out, "Component", ev, fit["trace"])
    A = V[:, :k]
    _table(out, "Principal components (eigenvectors) ", names, [f"Comp{j + 1}" for j in range(k)], A, _f4,
           extra=_unexpl_extra(unexpl), blanks=blanks)
    out.write("\n", "text")
    _post_pca(s, p, args, names, mask, N, C, means, S, ev, V, k, fit["trace"], rho, kind, wtype, w, cov)


def _post_pca(s, p, args, names, mask, N, C, means, S, ev, V, k, trace, rho, kind, wtype, w, cov, rotated=None):
    pdim = len(names)
    sign, lndet = np.linalg.slogdet(C)
    scalars = [("N", N), ("f", k), ("rho", rho), ("trace", trace), ("lndet", lndet if sign > 0 else SYS),
               ("cond", float(math.sqrt(ev[0] / ev[-1])) if ev[-1] > 0 else SYS)]
    cmdline = "pca " + args.strip()
    macros = [("cmdline", cmdline), ("cmd", "pca"), ("title", "Principal components"), ("marginsnotok", "_ALL"),
              ("estat_cmd", "pca_estat"), ("rotate_cmd", "pca_rotate"), ("predict", "pca_p"), ("Ctype", kind),
              ("properties", "nob noV eigen"), ("wtype", wtype), ("wexp", f"= {p.weight[1]}" if p.weight else "")]
    unexpl = mv.pca_unexplained(C, ev, V, k)
    comps = [f"Comp{j + 1}" for j in range(pdim)]
    mats = [("sds", matrix(np.sqrt(np.diag(S)).reshape(1, -1), ["sds"], names)),
            ("means", matrix(means.reshape(1, -1), ["means"], names)),
            ("C", matrix(C, names, names)),
            ("Psi", matrix(unexpl.reshape(1, -1), ["Psi"], names)),
            ("Ev", matrix(ev.reshape(1, -1), ["Ev"], comps)),
            ("L", matrix(V[:, :k], names, comps[:k]))]
    mvd = {"kind": "pca", "names": names, "C": C, "means": means, "sds": np.sqrt(np.diag(S)), "S": S, "ev": ev,
           "V": V, "k": k, "N": N, "trace": trace, "cov": cov, "rot": rotated}
    post_e(s, "pca", scalars, macros, mats, mask, N,
           estat={"loadings": (3, _estat_loadings), "residuals": (3, _estat_residuals), "smc": (3, _estat_smc),
                  "summarize": (3, _estat_summarize), "kmo": (3, _estat_kmo), "anti": (3, _estat_anti)},
           predict_multi={"score": lambda ss, ee: _pred_scores(ss, ee, show=True), "": lambda ss, ee: _pred_scores(ss, ee, show=True)},
           predict_opts={"fit": 3, "q": 1, "comp": 4}, extra={"mv": mvd})


# ---------------------------------------------------------------------------
# factor
# ---------------------------------------------------------------------------

_FACTOR_OPTS = {"pf": 2, "pcf": 3, "ipf": 3, "ml": 2, "factors": 7, "mineigen": 4, "blanks": 5,
                "heywood": 7, "tolerance": 3, "iterate": 4, "protect": 4, "random": 4, "altdivisor": 3,
                "citerate": 5, "ctolerance": 4, "covariance": 3, "norotated": 5}


@command("factor", "fac")
def cmd_factor(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _FACTOR_OPTS) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 3)
    methods = [m for m in ("pf", "pcf", "ipf", "ml") if o.get(m)]
    if len(methods) > 1:
        raise StataError(198, "only one of pf, pcf, ipf and ml may be specified")      # VERIFY
    method = methods[0] if methods else "pf"
    for opt in ("random", "altdivisor", "protect", "covariance"):
        if o.get(opt):
            raise StataError(198, f"option {opt} is not available yet in OpenDTA")
    X, w, wtype, mask = sample(s, p, names)
    means, S, R, sw = mv.weighted_moments(X, w)
    pdim = len(names)
    N = int(round(sw))
    kreq = int(float(o["factors"])) if o.get("factors") else None
    mine = float(o["mineigen"]) if o.get("mineigen") else None
    heyw = bool(o.get("heywood"))
    tol = float(o["tolerance"]) if o.get("tolerance") else (1e-4 if method == "ipf" else 1e-6)
    # VERIFY: pf and ipf keep the factors with positive eigenvalue, pcf those with eigenvalue >= 1
    if method == "pf":
        fit = mv.factor_pf(R, kreq, 0.0 if mine is None else mine)
    elif method == "ipf":
        fit = mv.factor_ipf(R, kreq, 0.0 if mine is None else mine, tol=tol, heywood=heyw)
    elif method == "pcf":
        fit = mv.factor_pcf(R, kreq, 1.0 if mine is None else mine)
    else:
        if kreq is None:
            # VERIFY: the default number of factors of ml is the number of positive eigenvalues of
            # the principal-factor solution, reduced until the model is identified
            kreq = max(1, mv.factor_pf(R, None)["L"].shape[1])
            while kreq > 1 and ((pdim - kreq) ** 2 - (pdim + kreq)) < 0:
                kreq -= 1
        if ((pdim - kreq) ** 2 - (pdim + kreq)) < 0:
            raise StataError(198, f"{kreq} factors will be too many for {pdim} variables")   # VERIFY
        fit = mv.factor_ml(R, kreq, heywood=heyw, tol=tol)
    L = fit["L"]
    k = L.shape[1]
    if k == 0:
        raise StataError(198, "no factors retained")                                    # VERIFY
    psi = fit["psi"]
    blanks = float(o["blanks"]) if o.get("blanks") else 0.0
    chi_i, df_i = mv.lr_independence(R, N)
    nparams = pdim * k - k * (k - 1) // 2
    out = s.output
    ll = ll0 = SYS
    if method == "ml":
        sign0, ld0 = np.linalg.slogdet(R)
        ll = float(-0.5 * N * fit["F"])
        ll0 = float(0.5 * N * ld0)
    hd = {"N": N, "method": method, "nparams": nparams, "ll": ll}
    out.write(f"(obs={N})\n", "text")
    if method == "ml":
        # the optimizer's own path: only the first and the last values are Stata's (VERIFY)
        for i_, F_ in enumerate(fit["history"]):
            out.write(f"Iteration {i_}:   log likelihood = {gnum(-0.5 * N * F_, '%10.0g')}\n", "text")
    _factor_head(out, hd, "(unrotated)", k)
    ev_show = fit["ev"] if method != "ml" else (L ** 2).sum(0)
    trace = float(fit["ev"].sum()) if method != "ml" else float(ev_show.sum())
    _eig_table(out, "Factor", ev_show, trace if trace else 1.0, fixed=True)
    _lr_line(out, "independent", chi_i, df_i, chi2_p(chi_i, df_i))
    chi_k = df_k = pk = None
    if method == "ml":
        chi_k, df_k = mv.lr_factors(fit["F"], N, pdim, k)
        pk = chi2_p(chi_k, df_k)
        if df_k > 0:
            _lr_line(out, f"{k} factor{'s' if k != 1 else ''}", chi_k, df_k, pk)
        if fit["heywood"]:
            out.write("\n", "text")
            out.write("Beware: solution is a Heywood case (i.e., invalid or boundary values of uniqueness)\n", "text")
    _table(out, "Factor loadings (pattern matrix) and unique variances", names, [f"Factor{j + 1}" for j in range(k)], L, _f4,
           extra=_uniq_extra(psi), blanks=blanks)
    out.write("\n", "text")
    dfr = ((pdim - k) ** 2 - (pdim + k)) / 2.0
    scalars = [("f", k), ("N", N), ("df_m", nparams), ("df_r", dfr)]
    if method == "ml":
        scalars += [("ll", ll), ("aic", -2 * ll + 2 * nparams), ("bic", -2 * ll + nparams * math.log(N))]
    scalars += [("chi2_i", chi_i), ("df_i", df_i), ("p_i", chi2_p(chi_i, df_i))]
    if method == "ml":
        scalars += [("chi2_1", chi_k), ("df_1", df_k), ("ll_0", ll0)]
    scalars.append(("evsum", trace))
    # VERIFY: the macros of the other methods follow the maximum likelihood list of the reference log
    macros = [("cmdline", "factor " + args.strip()), ("cmd", "factor"), ("marginsnotok", "_ALL"),
              ("properties", "nob noV eigen"), ("title", "Factor analysis"), ("predict", "factor_p"),
              ("estat_cmd", "factor_estat"), ("rotate_cmd", "factor_rotate"),
              ("factors", f"factors({o['factors']})" if o.get("factors") else ""),
              ("mtitle", _METHOD[method]), ("method", method), ("wtype", wtype)]
    cn = [f"Factor{j + 1}" for j in range(k)]
    mats = [("sds", matrix(np.sqrt(np.diag(S)).reshape(1, -1), ["sds"], names)),
            ("means", matrix(means.reshape(1, -1), ["means"], names)),
            ("C", matrix(R, names, names)), ("Phi", matrix(np.eye(k), cn, cn)), ("L", matrix(L, names, cn)),
            ("Psi", matrix(psi.reshape(1, -1), ["Psi"], names)),
            ("Ev", matrix(fit["ev"].reshape(1, -1), ["Ev"], [f"Factor{j + 1}" for j in range(len(fit["ev"]))]))]
    mvd = {"kind": "factor", "names": names, "C": R, "means": means, "sds": np.sqrt(np.diag(S)), "S": S,
           "ev": fit["ev"], "L": L, "psi": psi, "k": k, "N": N, "method": method, "Phi": np.eye(k), "rot": None,
           "chi_i": chi_i, "df_i": df_i, "nparams": nparams, "ml": (chi_k, df_k, pk, ll, ll0) if method == "ml" else None, "ll": ll,
           "evsum": trace}
    post_e(s, "factor", scalars, macros, mats, mask, N,
           estat={"kmo": (3, _estat_kmo), "smc": (3, _estat_smc), "anti": (3, _estat_anti), "summarize": (3, _estat_summarize),
                  "residuals": (3, _estat_residuals), "loadings": (3, _estat_loadings),
                  "structure": (3, _estat_structure), "common": (3, _estat_common)},
           predict_multi={"": lambda ss, ee: _pred_scores(ss, ee, show=True),
                          "regression": lambda ss, ee: _pred_scores(ss, ee, show=True),
                          "bartlett": lambda ss, ee: _pred_scores(ss, ee, "bartlett", show=True)},
           predict_opts={"regression": 3, "bartlett": 3, "fitted": 3, "residuals": 3}, extra={"mv": mvd})


def _lr_line(out, what: str, chi: float, df: int, p: float) -> None:
    out.write(f"{RULE}LR test: {what:>11} vs. saturated:  {f'chi2({df})':<8} = {chi:7.2f} Prob>chi2 = {p:6.4f}\n", "text")


def _factor_head(out, d: dict, desc: str, k: int) -> None:
    lines = [("", _right("Number of obs", f"{d['N']}")),
             (f"    Method: {_METHOD[d['method']]}", _right("Retained factors", f"{k}")),
             (f"    Rotation: {desc}", _right("Number of params", f"{d['nparams']}"))]
    if d["method"] == "ml":
        N, npar, ll = d["N"], d["nparams"], d["ll"]
        lines.append(("", _right("Schwarz's BIC", f"{-2 * ll + npar * math.log(N):.4f}")))
        lines.append((f"    Log likelihood = {gnum(ll, '%9.0g')}", _right("(Akaike's) AIC", f"{-2 * ll + 2 * npar:.4f}")))
    _header(out, "Factor analysis/correlation", lines)


_METHOD = {"pf": "principal factors", "pcf": "principal-component factors", "ipf": "iterated principal factors",
           "ml": "maximum likelihood"}


# ---------------------------------------------------------------------------
# rotate
# ---------------------------------------------------------------------------

_ROT = {"varimax": 3, "quartimax": 6, "equamax": 3, "parsimax": 3, "promax": 3, "oblimin": 3, "quartimin": 6,
        "normalize": 4, "factors": 6, "components": 4, "blanks": 5, "norotated": 5, "oblique": 3, "orthogonal": 4}


@command("rotate", "rot")
def cmd_rotate(s: "Session", args: str) -> None:
    est = getattr(s, "_estimates", None)
    d = (getattr(est, "extra", None) or {}).get("mv") if est is not None else None
    if d is None or d["kind"] not in ("pca", "factor"):
        raise StataError(301, "last estimates not found")
    o = match_options(args.strip().lstrip(","), _ROT) if args.strip().lstrip(",").strip() else {}
    methods = [m for m in ("varimax", "quartimax", "equamax", "parsimax", "promax", "oblimin", "quartimin") if o.get(m)]
    if len(methods) > 1:
        raise StataError(198, "only one rotation method may be specified")                   # VERIFY
    method = methods[0] if methods else "varimax"
    norm = bool(o.get("normalize"))
    kind = d["kind"]
    k = d["k"]
    kk = int(float(o.get("factors") or o.get("components") or k))
    if kk < 2:
        raise StataError(198, "rotation requires at least 2 factors/components")             # VERIFY
    kk = min(kk, k)
    # pca rotates the eigenvectors (unit length); the variances are those of the scaled eigenvectors
    A = d["V"][:, :kk] if kind == "pca" else d["L"][:, :kk]
    As = A * np.sqrt(np.maximum(d["ev"][:kk], 0)) if kind == "pca" else A
    pdim = A.shape[0]
    Phi = np.eye(kk)
    if method in ("varimax", "quartimax", "equamax", "parsimax"):
        gamma = {"varimax": 1.0, "quartimax": 0.0, "equamax": kk / 2.0, "parsimax": pdim * (kk - 1) / (pdim + kk - 2.0)}[method]
        Lr, T, _ = mv.orthomax(A, gamma, normalize=norm)
        desc = f"orthogonal {method}"
    elif method == "promax":
        power = float(o["promax"]) if isinstance(o["promax"], str) else 3.0
        Lr, T, Phi = mv.promax(A, power, normalize=norm)
        desc = "oblique promax"
        ctitle, crit = f"promax({power:g})", "promax"
    else:
        gamma = float(o["oblimin"]) if isinstance(o.get("oblimin"), str) else 0.0
        Lr, T, Phi, _ = mv.oblimin(A, gamma, normalize=norm)
        T = np.linalg.inv(T).T                      # the rotation matrix shown is the one with Lr = A T
        desc = "oblique quartimin" if method == "quartimin" else "oblique oblimin"
        ctitle, crit = f"oblimin({gamma:g})", "oblimin"
    if method in ("varimax", "quartimax", "equamax", "parsimax"):
        ctitle = crit = method
    oblique = method in ("promax", "oblimin", "quartimin")
    desc += f" (Kaiser {'on' if norm else 'off'})"
    blanks = float(o["blanks"]) if o.get("blanks") else 0.0
    out = s.output
    names = d["names"]
    N = d["N"]
    var = ((As @ T) ** 2).sum(0)
    if kind == "factor" and oblique:
        var = ((Lr @ Phi) ** 2).sum(0)                  # variance of the structure (VERIFY)
    if kind == "pca":
        unexpl = mv.pca_unexplained(d["C"], d["ev"], d["V"], kk)
        rho = float(d["ev"][:kk].sum() / d["trace"])
        _header(out, f"Principal components/{'covariance' if d['cov'] else 'correlation'}", [
            ("", _right("Number of obs", f"{N}")), ("", _right("Number of comp.", f"{kk}")),
            ("", _right("Trace", gnum(d["trace"]))),
            (f"    Rotation: {desc}", _right("Rho", f"{rho:.4f}"))])
        if oblique:
            _eig_table(out, "Component", var, d["trace"], head="Variance", oblique="Rotated comp. are correlated")
        else:
            _eig_table(out, "Component", var, d["trace"], head="Variance")
        _table(out, "Rotated components", names, [f"Comp{j + 1}" for j in range(kk)], Lr, _f4,
               extra=_unexpl_extra(unexpl), blanks=blanks)
        resid = unexpl
    else:
        _factor_head(out, d, desc, kk)
        if oblique:
            _eig_table(out, "Factor", var, d["evsum"], head="Variance", fixed=True, oblique="Rotated factors are correlated")
        else:
            _eig_table(out, "Factor", var, d["evsum"], head="Variance", fixed=True)
        _lr_line(out, "independent", d["chi_i"], d["df_i"], chi2_p(d["chi_i"], d["df_i"]))
        if d["ml"] is not None:
            chi_k, df_k, pk = d["ml"][:3]
            if df_k > 0:
                _lr_line(out, f"{kk} factor{'s' if kk != 1 else ''}", chi_k, df_k, pk)
        _table(out, "Rotated factor loadings (pattern matrix) and unique variances", names, [f"Factor{j + 1}" for j in range(kk)],
               Lr, _f4, extra=_uniq_extra(d["psi"]), blanks=blanks)
        resid = d["psi"]
    cn = ("Comp" if kind == "pca" else "Factor")
    cols = [f"{cn}{j + 1}" for j in range(kk)]
    # pca shows the matrix with rotated = unrotated T; factor shows its inverse transpose (VERIFY)
    Tshow = np.linalg.inv(T).T if kind == "factor" and oblique else T
    _table(out, ("Component" if kind == "pca" else "Factor") + " rotation matrix", cols, cols, Tshow, _f4, stub="",
           cw=9 if kind == "pca" else 8)
    d2 = dict(d)
    d2["rot"] = {"L": Lr, "T": T, "Phi": Phi, "method": method, "desc": desc, "k": kk, "oblique": oblique,
                 "ctitle": ctitle, "crit": crit, "norm": norm}
    est.extra["mv"] = d2
    base_s, base_m, base_t = est.extra["elist"]
    # VERIFY: the layout of the r_ results follows the oblimin case of the reference log; r_fmin
    # (the value of the rotation criterion) is stored only for oblimin
    rs = [("r_fmin", mv._oblimin_crit(Lr, float(ctitle[8:-1]) if ctitle.startswith("oblimin") else 0.0)[0]),
          ("r_f", kk)] if method in ("oblimin", "quartimin") else [("r_f", kk)]
    rm = [("r_normalization", "kaiser" if norm else "none"), ("r_ctitle", ctitle), ("r_criterion", crit),
          ("r_class", "oblique" if oblique else "orthogonal")]
    if kind == "factor":
        rm = [rm[0], rm[3], rm[2], rm[1]]
    rt = [("r_Ev", matrix(var.reshape(1, -1), ["r_Ev"], cols))]
    if kind == "factor":
        rt.append(("r_Phi", matrix(Phi, cols, cols)))
    rt += [("r_T", matrix(Tshow, cols, cols)), ("r_L", matrix(Lr, names, cols))]
    base_m = [x for x in base_m if not x[0].startswith("r_")]
    fill_e(s, [x for x in base_s if not x[0].startswith("r_")] + rs, rm + base_m, rt + [x for x in base_t if not x[0].startswith("r_")])


# ---------------------------------------------------------------------------
# estat
# ---------------------------------------------------------------------------

def _mvd(est) -> dict:
    return est.extra["mv"]


def _split(rest: str):
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    return body.strip(), opt


def _estat_loadings(s, est, rest):
    d = _mvd(est)
    body, opt = _split(rest)
    o = match_options(opt, {"cnorm": 2, "blanks": 5, "format": 3}) if opt.strip() else {}
    cn = str(o.get("cnorm", "unit")).strip()
    blanks = float(o["blanks"]) if o.get("blanks") else 0.0
    out = s.output
    if d["kind"] == "pca":
        k = d["rot"]["k"] if d["rot"] else d["k"]
        # components scaled to unit length (default), to length sqrt(eigenvalue) or its inverse
        A = d["rot"]["L"] if d["rot"] else d["V"][:, :k] * np.sqrt(np.maximum(d["ev"][:k], 0))
        norms = np.sqrt((A ** 2).sum(0))
        if cn.startswith("unit"):
            A = A / norms
            how = "1"
        elif cn.startswith("inv"):
            A = A / norms ** 2
            how = "1/eigenvalue"                                                         # VERIFY
        else:
            how = "eigenvalue"                                                           # VERIFY
        _table(out, f"Principal component loadings ({'rotated' if d['rot'] else 'unrotated'})", d["names"],
               [f"Comp{j + 1}" for j in range(k)], A, _g4, stub="", blanks=blanks,
               sub=f"    component normalization: sum of squares(column) = {how}")
    else:
        L = d["rot"]["L"] if d["rot"] else d["L"]
        _table(out, "Factor loadings", d["names"], [f"Factor{j + 1}" for j in range(L.shape[1])], L, _g4, stub="",
               blanks=blanks)


def _estat_smc(s, est, rest):
    d = _mvd(est)
    sm = mv.smc(d["C"])
    _table(s.output, "Squared multiple correlations of variables with all other variables", d["names"], ["smc"],
           sm.reshape(-1, 1), _f4, cw=8)
    s.r = {"smc": matrix(sm.reshape(-1, 1), d["names"], ["smc"])}


def _estat_kmo(s, est, rest):
    d = _mvd(est)
    k_, ov = mv.kmo(d["C"])
    _table(s.output, "Kaiser-Meyer-Olkin measure of sampling adequacy", d["names"], ["kmo"], k_.reshape(-1, 1), _f4,
           cw=8, foot=[("Overall", [ov])])
    s.r = {"kmo": ov, "kmo_j": matrix(k_.reshape(-1, 1), d["names"], ["kmo"])}


def _estat_anti(s, est, rest):
    d = _mvd(est)
    body, opt = _split(rest)
    o = match_options(opt, {"nocovariance": 5, "nocorrelation": 5, "format": 3}) if opt.strip() else {}
    corr, cov = mv.anti_image(d["C"])
    names = d["names"]
    if not o.get("nocorrelation"):
        _print_matrix(s.output, "Anti-image correlation matrix", names, corr)
    if not o.get("nocovariance"):
        _print_matrix(s.output, "Anti-image covariance matrix", names, cov)
    s.r = {"Cov": matrix(cov, names, names), "Corr": matrix(corr, names, names)}


def _estat_summarize(s, est, rest):
    d = _mvd(est)
    out = s.output
    ds = s.data
    mask = ds.esample
    out.write(f"\n  Estimation sample {est.cmd:<23}Number of obs = {d['N']:>10}\n\n", "text")
    rule = "-" * 53
    out.write(f"  {'-' * 13}{'-' * 54}\n", "text")
    out.write(f"  {'Variable':>12} |{'Mean':>13}{'Std. Dev.':>15}{'Min':>12}{'Max':>13}\n", "text")
    out.write(f"  {'-' * 13}+{rule}\n", "text")
    for n_, m_, sd in zip(d["names"], d["means"], d["sds"]):
        x = ds.get(n_).data[mask]
        out.write(f"  {_abbr(n_):>12} |{gnum(m_, '%9.0g'):>13}{gnum(sd, '%9.0g'):>14}"
                  f"{gnum(x.min(), '%9.0g'):>13}{gnum(x.max(), '%9.0g'):>13}\n", "text")
    out.write(f"  {'-' * 13}{'-' * 54}\n", "text")
    if est.cmd == "pca":
        out.write("(PCA of covariance matrix)\n" if d["cov"] else
                  "(PCA of correlation matrix, i.e., based on standardized variables)\n", "text")      # VERIFY covariance


def _fitted(d):
    if d["kind"] == "pca":
        k = d["k"]
        L = d["V"][:, :k] * np.sqrt(np.maximum(d["ev"][:k], 0))
        return L @ L.T
    L = d["rot"]["L"] if d["rot"] else d["L"]
    Phi = d["rot"]["Phi"] if d["rot"] else np.eye(L.shape[1])
    return L @ Phi @ L.T + np.diag(d["psi"])


def _estat_residuals(s, est, rest):
    d = _mvd(est)
    body, opt = _split(rest)
    o = match_options(opt, {"fitted": 3, "obs": 3, "format": 3}) if opt.strip() else {}
    F = _fitted(d)
    names = d["names"]
    what = "covariance" if d.get("cov") else "correlation"
    g7 = lambda v: gnum(v, "%9.0g")
    if o.get("fitted"):
        _table(s.output, f"Fitted {what} matrix", names, names, F, g7, cw=10, tri=True)             # VERIFY
        s.r = {"Fit": matrix(F, names, names)}
    else:
        Rm = d["C"] - F
        _table(s.output, f"Residual {what} matrix", names, names, Rm, g7, cw=10, tri=True)
        iu = np.triu_indices(len(names), 1)
        s.r = {"rmse_off": float(math.sqrt((Rm[iu] ** 2).mean())), "Rm": matrix(Rm, names, names)}


def _estat_structure(s, est, rest):
    d = _mvd(est)
    if d["kind"] != "factor":
        raise StataError(321, "subcommand estat structure is unrecognized")
    L = d["rot"]["L"] if d["rot"] else d["L"]
    Phi = d["rot"]["Phi"] if d["rot"] else np.eye(L.shape[1])
    cols = [f"Factor{j + 1}" for j in range(L.shape[1])]
    _table(s.output, "Structure matrix: correlations between variables and common factors", d["names"], cols, L @ Phi, _f4)
    s.r = {"Lstructure": matrix(L @ Phi, d["names"], cols)}


def _estat_common(s, est, rest):
    d = _mvd(est)
    if d["kind"] != "factor":
        raise StataError(321, "subcommand estat common is unrecognized")
    Phi = d["rot"]["Phi"] if d["rot"] else np.eye(d["k"])
    cols = [f"Factor{j + 1}" for j in range(Phi.shape[0])]
    _table(s.output, "Correlation matrix of the common factors", cols, cols, Phi, lambda v: gnum(v, "%9.0g"), stub="Factors",
           tri=True)
    s.r = {"Phi": matrix(Phi, cols, cols)}


# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

def _score_weights(d, how: str) -> np.ndarray:
    """Scoring coefficients of the standardized variables (pca: eigenvectors, or the rotated ones
    through the inverse transpose of the rotation matrix; factor: regression or Bartlett)."""
    rot = d["rot"]
    if d["kind"] == "pca":
        k = rot["k"] if rot else d["k"]
        B = d["V"][:, :k]
        return B @ np.linalg.inv(rot["T"]).T if rot else B
    L = rot["L"] if rot else d["L"]
    Phi = rot["Phi"] if rot else np.eye(L.shape[1])
    if how == "bartlett":
        Pi = np.diag(1.0 / d["psi"])
        return Pi @ L @ np.linalg.inv(L.T @ Pi @ L)
    return np.linalg.solve(d["C"], L @ Phi)


def _show_weights(s, d, W, how: str) -> None:
    rot = d["rot"]
    cols = [f"{'Comp' if d['kind'] == 'pca' else 'Factor'}{j + 1}" for j in range(W.shape[1])]
    if d["kind"] == "pca":
        if rot:
            _table(s.output, f"Scoring coefficients for {'oblique' if rot['oblique'] else 'orthogonal'} {rot['ctitle']} rotation",
                   d["names"], cols, W, _f4, stub="Variable", sub="    sum of squares(column-loading) = 1")  # VERIFY orthogonal
        return
    based = f"; based on {rot['ctitle']} rotated factors" if rot else ""
    _table(s.output, f"Scoring coefficients (method = {how}{based})", d["names"], cols, W, lambda v: f"{v:.5f}")


def _pred_scores(s, est, how: str = "regression", show: bool = False):
    d = _mvd(est)
    o = getattr(est, "_pred_opts", {}) or {}
    names = d["names"]
    k = d["k"]
    if o.get("comp"):
        cols = [int(float(x)) - 1 for x in str(o["comp"]).replace(",", " ").split()]
    else:
        cols = list(range(d["rot"]["k"] if d["rot"] else k))
    nreq = len(cols)
    ds = s.data

    def fn(mask):
        W = _score_weights(d, how)
        if show:
            if d["kind"] == "factor" and not (o.get("regression") or o.get("bartlett")):
                s.output.write("(regression scoring assumed)\n", "text")
            _show_weights(s, d, W, how)
        X = np.column_stack([ds.get(n).data.astype(float) for n in names])
        ok = mask & np.all(X < SYS, axis=1)
        Z = np.where(ok[:, None], X - d["means"], np.nan)
        if not (d["kind"] == "pca" and d["cov"]):
            Z = Z / d["sds"]
        Sc = Z @ W
        return [np.where(ok, Sc[:, j], SYS) for j in cols]

    return nreq, fn, f"too few variables specified\nyou must specify {nreq} new variable(s)"


# ---------------------------------------------------------------------------
# alpha
# ---------------------------------------------------------------------------

_ALPHA = {"asis": 2, "casewise": 4, "detail": 1, "generate": 1, "item": 1, "label": 1, "min": 3, "reverse": 7,
          "std": 3}


@command("alpha", "alpha")
def cmd_alpha(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, _ALPHA) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 2)
    ds = s.data
    from ._util import touse
    mask = touse(s, p)
    cols = np.column_stack([ds.get(n).data for n in names]).astype(float)
    miss = cols >= SYS
    cols = np.where(miss, np.nan, cols)
    casewise = bool(o.get("casewise"))
    X = cols[mask]
    if casewise:
        X = X[~np.isnan(X).any(1)]
    if X.shape[0] == 0:
        raise StataError(2000, "no observations")
    k = len(names)
    std = bool(o.get("std"))
    rev = set(s.expand_varlist(str(o["reverse"]))) if isinstance(o.get("reverse"), str) else set()
    if not set(rev) <= set(names):
        raise StataError(198, "reverse() variables must be in varlist")                   # VERIFY
    sign = np.array([-1.0 if n in rev else 1.0 for n in names])
    if not o.get("asis") and not rev:
        # VERIFY: signs from the first principal component of the correlation matrix (Stata 14 default)
        C0, _ = mv.pairwise_cov(X)
        sd0 = np.sqrt(np.diag(C0))
        R0 = C0 / np.outer(sd0, sd0)
        v = mv.sorted_eigh(np.nan_to_num(R0))[1][:, 0]
        if v.sum() < 0:
            v = -v
        sign = np.where(v < 0, -1.0, 1.0)
    X = X * sign
    C, Nm = mv.pairwise_cov(X)
    res = mv.alpha_from_cov(C, std)
    out = s.output
    out.write(f"\nTest scale = mean({'standardized' if std else 'unstandardized'} items)\n", "text")
    avg = res["rbar"] if std else res["cbar"]
    if o.get("item"):
        _alpha_items(s, names, X, sign, C, Nm, std, res)
    else:
        out.write("\n", "text")
        for lab, val in (("Average interitem " + ("correlation:" if std else "covariance:"), f"{avg:.4f}" if std else gnum(avg, "%9.0g")),
                         ("Number of items in the scale:", f"{k}"),
                         ("Scale reliability coefficient:", f"{res['alpha']:.4f}")):
            out.write(lab + val.rjust(42 - len(lab)) + "\n", "text")
        out.write("\n", "text")
    s.r = {"rho" if std else "cov": avg, "k": float(k), "alpha": res["alpha"]}
    if o.get("detail"):
        sd = np.sqrt(np.diag(C))
        M_ = np.nan_to_num(C / np.outer(sd, sd) if std else C)
        lo, hi = int(Nm.min()), int(Nm.max())
        out.write(f"\nInteritem {'correlations' if std else 'covariances'} "
                  f"(obs={lo if lo == hi else f'{lo} to {hi}'}{' in all pairs' if lo == hi else ''})\n\n", "text")
        lower_matrix(out, names, M_)
    if isinstance(o.get("generate"), str):
        gname = str(o["generate"]).strip()
        if ds.has(gname):
            raise StataError(110, f"{gname} already defined")
        full = cols * sign
        if std:
            full = (full - np.nanmean(X, 0)) / np.nanstd(X, 0, ddof=1)
        mn = int(float(o["min"])) if o.get("min") else k
        cnt = (~np.isnan(full)).sum(1)
        sc = np.where(mask & (cnt >= mn), np.nansum(full, 1) / np.maximum(cnt, 1), SYS)
        var = Variable(gname, "float", np.full(ds.nobs, SYS))
        var.data = sc
        ds.add(var)
        nm = int((sc >= SYS).sum())
        if nm:
            out.write(f"({nm:,} missing value{'s' if nm != 1 else ''} generated)\n", "text")
        s.notify_state()


def _alpha_items(s, names, X, sign, C, Nm, std, res) -> None:
    out = s.output
    k = len(names)
    ok = ~np.isnan(X).any(1)
    Xc = X[ok]
    tot = Xc.sum(1)
    sd = np.sqrt(np.diag(C))
    R = C / np.outer(sd, sd)
    off = lambda M_: float(M_[~np.eye(len(M_), dtype=bool)].mean())
    rows = []
    for j in range(k):
        it = Xc[:, j]
        rest = tot - it
        itc = float(np.corrcoef(it, tot)[0, 1])
        irc = float(np.corrcoef(it, rest)[0, 1])
        idx = [i for i in range(k) if i != j]
        sub = C[np.ix_(idx, idx)]
        a_j = mv.alpha_from_cov(sub, std)["alpha"] if k > 2 else float("nan")
        inter = off(R[np.ix_(idx, idx)] if std else sub) if k > 2 else float("nan")
        rows.append((names[j], int(Nm[j, j]), "-" if sign[j] < 0 else "+", itc, irc, inter, a_j))
    avg = res["rbar"] if std else res["cbar"]
    lab = "correlation" if std else "covariance"
    val = (lambda v: f"{v:.4f}") if std else gnum
    out.write("\n", "text")
    out.write(f"{'average':>67}\n{'item-test':>38}{'item-rest':>14}{'interitem':>16}\n", "text")
    out.write(f"{'Item':<13}|{'Obs':>5}{'Sign':>6}{'correlation':>14}{'correlation':>14}{lab:>15}{'alpha':>11}\n", "text")
    out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")
    for n_, ob, sg, itc, irc, inter, a_j in rows:
        out.write(f"{n_[:12]:<13}|{ob:>5}{sg:>5}{itc:>13.4f}{irc:>14.4f}{val(inter):>16}{a_j:>12.4f}\n", "text")
    out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")
    out.write(f"{'Test scale':<13}|{val(avg):>53}{res['alpha']:>12.4f}\n", "text")
    out.write("-" * 79 + "\n", "text")


# ---------------------------------------------------------------------------
# screeplot, scoreplot, loadingplot
# ---------------------------------------------------------------------------

@command("screeplot", "screeplot")
def cmd_screeplot(s: "Session", args: str) -> None:
    from ..graphics.model import Graph, Layer, Panel
    from .graphcmd import finish
    est = getattr(s, "_estimates", None)
    p = parse_standard(args)
    d = (getattr(est, "extra", None) or {}).get("mv") if est is not None else None
    if p.varlist.strip():
        saved = (s.e, s.data.esample, est)
        with s.output.quietly():
            cmd_pca(s, args if "," not in args else args)
        est2 = s._estimates
        d = est2.extra["mv"]
        s.e, s.data.esample, s._estimates = saved
    if d is None:
        raise StataError(301, "last estimates not found")
    ev = d["ev"]
    o = {}
    mean = False
    if p.options.strip():
        o = match_options(p.options, {"mean": 2, "ci": 2, "yline": 5, "title": 3, "ytitle": 3, "xtitle": 3, "name": 3,
                                     "saving": 3, "nodraw": 3, "components": 4, "correlation": 4, "covariance": 3,
                                     "mineigen": 4, "blanks": 5, "vce": 3, "note": 3, "subtitle": 3, "caption": 3})
        mean = bool(o.get("mean"))
    layers = [Layer("connected", list(range(1, len(ev) + 1)), [float(x) for x in ev], {}, {}, "Eigenvalues")]
    if mean:
        m = float(np.mean(ev))
        layers.append(Layer("line", [1.0, float(len(ev))], [m, m], {}, {"lpattern": "dash"}, "Mean"))
    opts = {"_xtitle": "Number", "_ytitle": "Eigenvalues", "title": {"lines": ["Scree plot of eigenvalues after " + est.cmd if est is not None else "Scree plot of eigenvalues"]}}
    opts["legend"] = {"off": True}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), {k: v for k, v in o.items() if k in ("name", "saving", "nodraw")})


def _plot_xy(s, args: str, scores: bool) -> None:
    from ..graphics.model import Graph, Layer, Panel
    from .graphcmd import finish
    est = getattr(s, "_estimates", None)
    d = (getattr(est, "extra", None) or {}).get("mv") if est is not None else None
    if d is None:
        raise StataError(301, "last estimates not found")
    p = parse_standard(args)
    o = match_options(p.options, {"name": 3, "saving": 3, "nodraw": 3, "title": 3}) if p.options.strip() else {}
    k = d["rot"]["k"] if d["rot"] else d["k"]
    if k < 2:
        raise StataError(198, "at least two components/factors are required")           # VERIFY
    word = "Component" if d["kind"] == "pca" else "Factor"
    if scores:
        nreq, fn, _ = _pred_scores(s, est)
        ds = s.data
        sc = fn(ds.esample)
        ok = (sc[0] < SYS) & (sc[1] < SYS)
        x, y, lab = sc[0][ok].tolist(), sc[1][ok].tolist(), None
        title = f"Score variables ({word.lower()}s)"
    else:
        A = d["rot"]["L"] if d["rot"] else (d["V"][:, :k] * np.sqrt(np.maximum(d["ev"][:k], 0)) if d["kind"] == "pca" else d["L"])
        x, y, lab = A[:, 0].tolist(), A[:, 1].tolist(), list(d["names"])
        title = f"{word} loadings"
    extra = {"mlabel": lab} if lab else {}
    opts = {"_xtitle": f"{word} 1", "_ytitle": f"{word} 2", "title": {"lines": [title]}, "legend": {"off": True}}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", [Layer("scatter", x, y, extra, {}, "")])]), o)


@command("loadingplot", "loadingplot")
def cmd_loadingplot(s: "Session", args: str) -> None:
    _plot_xy(s, args, False)


@command("scoreplot", "scoreplot")
def cmd_scoreplot(s: "Session", args: str) -> None:
    _plot_xy(s, args, True)
