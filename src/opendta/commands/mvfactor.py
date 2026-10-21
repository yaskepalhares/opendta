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
from .mvcommon import chi2_p, f as fnum, g as gnum, matrix, numeric_names, post_e, sample
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
RULE = "    "


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------

def _header(out, title: str, lines: list[tuple[str, str]]) -> None:
    """Title and indented lines on the left, `label = value` on the right (column 49)."""
    out.write("\n", "text")
    for i, (left, right) in enumerate(lines):
        if i == 0:
            left = title
        out.write(f"{left:<48}{right}\n", "text")


def _right(label: str, value: str) -> str:
    return f"{label:<17}= {value:>10}"


def _eig_table(out, name: str, ev, trace: float, head: str = "Eigenvalue") -> None:
    rule1 = RULE + "-" * 74
    out.write("\n" + rule1 + "\n", "text")
    out.write(f"{RULE}{name:>12} |{head:>13}{'Difference':>13}{'Proportion':>19}{'Cumulative':>13}\n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 60}\n", "text")
    cum = 0.0
    p = len(ev)
    for j in range(p):
        diff = ev[j] - ev[j + 1] if j + 1 < p else None
        prop = ev[j] / trace if trace else float("nan")
        cum += prop
        d = f"{gnum(diff, '%9.0g'):>13}" if diff is not None else " " * 13
        out.write(f"{RULE}{_stub(name, j + 1):>12} |{gnum(ev[j], '%9.0g'):>13}{d}"
                  f"{prop:>19.4f}{cum:>13.4f}\n", "text")
    out.write(rule1 + "\n", "text")


def _stub(name: str, j: int) -> str:
    return {"Component": f"Comp{j}", "Factor": f"Factor{j}"}.get(name, f"{name}{j}")


def _load_table(out, title: str, names: list[str], A: np.ndarray, colname: str, extra: np.ndarray | None,
                extra_name: str, blanks: float = 0.0, dec: int = 4) -> None:
    k = A.shape[1]
    head = f"{RULE}{'Variable':>12} |" + "".join(f"{colname + str(j + 1):>10}" for j in range(k))
    sep = f"{RULE}{'-' * 13}+{'-' * (10 * k + 1)}"
    if extra is not None:
        head += f" |{extra_name:>14} "
        sep += f"+{'-' * 14}"
    else:
        head += " "
    rows = []
    for i, nm in enumerate(names):
        cells = "".join(f"{'':>10}" if abs(A[i, j]) < blanks else f"{A[i, j]:>10.{dec}f}" for j in range(k))
        rows.append(f"{RULE}{_abbr(nm):>12} |{cells}" + (f" |{extra[i]:>14.4f}" if extra is not None else ""))
    rule = RULE + "-" * (max(len(head), len(sep)) - len(RULE))
    out.write(f"\n{title}\n\n{rule}\n{head}\n{sep}\n" + "\n".join(rows) + f"\n{rule}\n", "text")


def _abbr(name: str, w: int = 12) -> str:
    return name if len(name) <= w else name[: w - 1] + "~"


def _print_matrix(out, title: str, names: list[str], M_: np.ndarray, cols: list[str] | None = None,
                  dec: int = 4) -> None:
    cols = cols or [_abbr(n, 8) for n in names]
    out.write(f"\n{title}\n\n{RULE}{'':>12} |" + "".join(f"{c:>10}" for c in cols) + "\n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * (10 * len(cols) + 1)}\n", "text")
    for i, nm in enumerate(names):
        out.write(f"{RULE}{_abbr(nm):>12} |" + "".join(f"{M_[i, j]:>10.{dec}f}" for j in range(M_.shape[1])) + "\n", "text")


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
    _load_table(out, "Principal components (eigenvectors) ", names, A, "Comp",
                unexpl if k < pdim else None, "Unexplained", blanks)
    out.write("\n", "text")
    _post_pca(s, p, args, names, mask, N, C, means, S, ev, V, k, fit["trace"], rho, kind, wtype, w, cov)


def _post_pca(s, p, args, names, mask, N, C, means, S, ev, V, k, trace, rho, kind, wtype, w, cov, rotated=None):
    pdim = len(names)
    sign, lndet = np.linalg.slogdet(C)
    scalars = [("N", N), ("f", k), ("rho", rho), ("trace", trace), ("lndet", lndet if sign > 0 else SYS)]
    cmdline = "pca " + args.strip()
    macros = [("cmdline", cmdline), ("cmd", "pca"), ("marginsnotok", "_ALL"), ("predict", "pca_p"),
              ("estat_cmd", "pca_estat"), ("rotate_cmd", "pca_rotate"), ("properties", "nob noV eigen"),
              ("title", "Principal components"), ("Ctype", kind), ("wtype", wtype), ("wexp", f"= {p.weight[1]}" if p.weight else "")]
    unexpl = mv.pca_unexplained(C, ev, V, k)
    mats = [("Ev", matrix(ev.reshape(1, -1), ["Ev"], [f"Comp{j + 1}" for j in range(pdim)])),
            ("L", matrix(V, names, [f"Comp{j + 1}" for j in range(pdim)])),
            ("Psi", matrix(unexpl.reshape(1, -1), ["Psi"], names)),
            ("C", matrix(C, names, names))]
    mats.insert(0, ("means", matrix(means.reshape(1, -1), ["means"], names)))
    mats.insert(1, ("sds", matrix(np.sqrt(np.diag(S)).reshape(1, -1), ["sds"], names)))
    mvd = {"kind": "pca", "names": names, "C": C, "means": means, "sds": np.sqrt(np.diag(S)), "S": S, "ev": ev,
           "V": V, "k": k, "N": N, "trace": trace, "cov": cov, "rot": rotated}
    post_e(s, "pca", scalars, macros, mats, mask, N,
           estat={"loadings": (3, _estat_loadings), "residuals": (3, _estat_residuals), "smc": (3, _estat_smc),
                  "summarize": (3, _estat_summarize), "correlation": (3, _estat_correlation),
                  "covariance": (3, _estat_correlation), "kmo": (3, _estat_kmo), "anti": (3, _estat_anti)},
           predict_multi={"score": lambda ss, ee: _pred_scores(ss, ee), "": lambda ss, ee: _pred_scores(ss, ee)},
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
    tol = float(o["tolerance"]) if o.get("tolerance") else 1e-6
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
    nparams = pdim * k + pdim - k * (k - 1) // 2                                       # VERIFY
    out = s.output
    lines = [("", _right("Number of obs", f"{N}")),
             (f"    Method: {_METHOD[method]}", _right("Retained factors", f"{k}")),
             ("    Rotation: (unrotated)", _right("Number of params", f"{nparams}"))]
    extra_lines = []
    ll = ll0 = SYS
    if method == "ml":
        sign0, ld0 = np.linalg.slogdet(R)
        ll = float(-0.5 * N * (fit["F"] + ld0 + pdim * math.log(2 * math.pi) + pdim))
        ll0 = float(-0.5 * N * (pdim * math.log(2 * math.pi) + pdim))
        bic = -2 * ll + nparams * math.log(N)
        aic = -2 * ll + 2 * nparams
        lines.append((f"Log likelihood = {gnum(ll, '%10.0g')}", _right("Schwarz's BIC", gnum(bic, "%10.0g"))))
        lines.append(("", _right("(Akaike's) AIC", gnum(aic, "%10.0g"))))
    _header(out, "Factor analysis/correlation", lines)
    ev_show = fit["ev"] if method != "ml" else np.maximum(fit["ev"][:k] - 1.0, 0.0)
    trace = float(fit["ev"].sum()) if method != "ml" else float(ev_show.sum())
    _eig_table(out, "Factor", ev_show, trace if trace else 1.0)
    out.write(f"{RULE}LR test: independent vs. saturated:  chi2({df_i}) = {chi_i:8.2f} Prob>chi2 = {chi2_p(chi_i, df_i):6.4f}\n", "text")
    chi_k = df_k = pk = None
    if method == "ml":
        chi_k, df_k = mv.lr_factors(fit["F"], N, pdim, k)
        pk = chi2_p(chi_k, df_k)
        if df_k > 0:
            out.write(f"{RULE}LR test: {k} factor{'s' if k != 1 else ''} vs. saturated:  chi2({df_k}) = {chi_k:8.2f} Prob>chi2 = {pk:6.4f}\n", "text")
        if fit["heywood"]:
            out.write("\n", "text")
            out.write("Beware: solution is a Heywood case (i.e., invalid or boundary values of uniqueness)\n", "text")
    _load_table(out, "Factor loadings (pattern matrix) and unique variances", names, L, "Factor", psi, "Uniqueness", blanks)
    out.write("\n", "text")
    scalars = [("N", N), ("f", k), ("evsum", float(fit["ev"].sum())), ("df_m", nparams), ("df_r", pdim * (pdim + 1) / 2.0 - nparams),
               ("chi2_i", chi_i), ("df_i", df_i), ("p_i", chi2_p(chi_i, df_i))]
    if method == "ml":
        scalars += [("ll", ll), ("ll_0", ll0), ("chi2_1", chi_k), ("p", pk if pk is not None else SYS)]
    scalars.append(("lndet", float(np.linalg.slogdet(R)[1])))
    scalars.append(("rho", float(fit["ev"][:k].sum() / fit["ev"].sum()) if method != "ml" else SYS))
    macros = [("cmdline", "factor " + args.strip()), ("cmd", "factor"), ("method", method), ("title", "Factor analysis"),
              ("rotate_cmd", "factor_rotate"), ("estat_cmd", "factor_estat"), ("predict", "factor_p"),
              ("marginsnotok", "_ALL"), ("properties", "nob noV eigen"), ("wtype", wtype)]
    cn = [f"Factor{j + 1}" for j in range(k)]
    mats = [("Phi", matrix(np.eye(k), cn, cn)), ("C", matrix(R, names, names)),
            ("Ev", matrix(fit["ev"].reshape(1, -1), ["Ev"], [f"Factor{j + 1}" for j in range(len(fit["ev"]))])),
            ("Psi", matrix(psi.reshape(1, -1), ["Psi"], names)), ("L", matrix(L, names, cn))]
    mvd = {"kind": "factor", "names": names, "C": R, "means": means, "sds": np.sqrt(np.diag(S)), "S": S,
           "ev": fit["ev"], "L": L, "psi": psi, "k": k, "N": N, "method": method, "Phi": np.eye(k), "rot": None,
           "chi_i": chi_i, "df_i": df_i, "nparams": nparams, "ml": (chi_k, df_k, pk, ll, ll0) if method == "ml" else None}
    post_e(s, "factor", scalars, macros, mats, mask, N,
           estat={"kmo": (3, _estat_kmo), "smc": (3, _estat_smc), "anti": (3, _estat_anti), "summarize": (3, _estat_summarize),
                  "correlation": (3, _estat_correlation), "covariance": (3, _estat_correlation),
                  "residuals": (3, _estat_residuals), "loadings": (3, _estat_loadings),
                  "structure": (3, _estat_structure), "common": (3, _estat_common)},
           predict_multi={"": lambda ss, ee: _pred_scores(ss, ee), "regression": lambda ss, ee: _pred_scores(ss, ee),
                          "bartlett": lambda ss, ee: _pred_scores(ss, ee, "bartlett")},
           predict_opts={"regression": 3, "bartlett": 3, "fitted": 3, "residuals": 3}, extra={"mv": mvd})


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
    # VERIFY: after pca the eigenvectors scaled by sqrt(eigenvalue) are rotated
    A = d["V"][:, :kk] * np.sqrt(np.maximum(d["ev"][:kk], 0)) if kind == "pca" else d["L"][:, :kk]
    pdim = A.shape[0]
    Phi = np.eye(kk)
    if method in ("varimax", "quartimax", "equamax", "parsimax"):
        gamma = {"varimax": 1.0, "quartimax": 0.0, "equamax": kk / 2.0, "parsimax": pdim * (kk - 1) / (pdim + kk - 2.0)}[method]
        Lr, T, _ = mv.orthomax(A, gamma, normalize=norm)
        desc = f"orthogonal {method}"
    elif method == "promax":
        power = float(o["promax"]) if isinstance(o["promax"], str) else 3.0
        Lr, T, Phi = mv.promax(A, power, normalize=norm)
        desc = f"oblique promax({power:g})"
    else:
        gamma = float(o["oblimin"]) if isinstance(o.get("oblimin"), str) else 0.0
        Lr, T, Phi, _ = mv.oblimin(A, gamma, normalize=norm)
        desc = "oblique quartimin" if method == "quartimin" or gamma == 0 else f"oblique oblimin({gamma:g})"
    oblique = method in ("promax", "oblimin", "quartimin")
    desc += f" (Kaiser {'on' if norm else 'off'})"
    blanks = float(o["blanks"]) if o.get("blanks") else 0.0
    out = s.output
    names = d["names"]
    N = d["N"]
    if kind == "pca":
        unexpl = np.diag(d["C"]) - (Lr ** 2).sum(1) if not oblique else np.diag(d["C"]) - np.einsum("ij,jk,ik->i", Lr, Phi, Lr)
        _header(out, f"Principal components/{'covariance' if d['cov'] else 'correlation'}", [
            ("", _right("Number of obs", f"{N}")), ("", _right("Number of comp.", f"{kk}")),
            ("", _right("Trace", gnum(d["trace"]))),
            (f"    Rotation: {desc}", _right("Rho", f"{float((Lr ** 2).sum() / d['trace']):.4f}"))])
        var = (Lr ** 2).sum(0)
        _eig_table(out, "Component", var, d["trace"], head="Variance")
        disp = Lr / np.sqrt((Lr ** 2).sum(0))
        _load_table(out, "Rotated components", names, disp, "Comp", unexpl if kk < len(names) else None, "Unexplained", blanks)
        resid = unexpl
    else:
        nparams = pdim * kk + pdim - kk * (kk - 1) // 2                                # VERIFY
        _header(out, "Factor analysis/correlation", [
            ("", _right("Number of obs", f"{N}")),
            (f"    Method: {_METHOD[d['method']]}", _right("Retained factors", f"{kk}")),
            (f"    Rotation: {desc}", _right("Number of params", f"{nparams}"))])
        var = (Lr ** 2).sum(0)
        _eig_table(out, "Factor", var, float(var.sum()), head="Variance")
        out.write(f"{RULE}LR test: independent vs. saturated:  chi2({d['df_i']}) = {d['chi_i']:8.2f} Prob>chi2 = {chi2_p(d['chi_i'], d['df_i']):6.4f}\n", "text")
        _load_table(out, "Rotated factor loadings (pattern matrix) and unique variances", names, Lr, "Factor", d["psi"], "Uniqueness", blanks)
        resid = d["psi"]
    cn = ("Comp" if kind == "pca" else "Factor")
    cols = [f"{cn}{j + 1}" for j in range(kk)]
    _print_matrix(out, ("Component" if kind == "pca" else "Factor") + " rotation matrix", cols, T, cols)
    if oblique:
        _print_matrix(out, "Factor correlations" if kind == "factor" else "Component correlations", cols, Phi, cols)
    out.write("\n", "text")
    d2 = dict(d)
    d2["rot"] = {"L": Lr, "T": T, "Phi": Phi, "method": method, "desc": desc, "k": kk, "oblique": oblique}
    est.extra["mv"] = d2
    mats = [("r_L", matrix(Lr, names, cols)), ("r_T", matrix(T, cols, cols)), ("r_Phi", matrix(Phi, cols, cols))]
    s.e.update({k_: m for k_, m in mats})
    s.e["r_criterion"] = float((Lr ** 2).sum())
    s.e["r_class"] = "oblique" if oblique else "orthogonal"
    s.e["r_ctitle"] = method
    s.e["r_normalize"] = "on" if norm else "off"


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
        elif cn.startswith("inv"):
            A = A / norms ** 2
        out.write(f"\nPrincipal components/{'covariance' if d['cov'] else 'correlation'}\n", "text")
        _load_table(out, "Rotated components" if d["rot"] else "Principal components (eigenvectors) ",
                    d["names"], A, "Comp", None, "", blanks)
    else:
        L = d["rot"]["L"] if d["rot"] else d["L"]
        _load_table(out, "Factor loadings", d["names"], L, "Factor", None, "", blanks)
    out.write("\n", "text")


def _estat_correlation(s, est, rest):
    d = _mvd(est)
    body, opt = _split(rest)
    o = match_options(opt, {"covariance": 3, "format": 3}) if opt.strip() else {}
    names = d["names"]
    M_ = d["S"] if o.get("covariance") else d["C"]
    title = "Covariance matrix" if o.get("covariance") else "Correlation matrix"
    _print_matrix(s.output, title, names, M_)
    s.output.write("\n", "text")
    s.r = {"C": matrix(M_, names, names)}
    s.r["C"].colnames = list(names)


def _estat_smc(s, est, rest):
    d = _mvd(est)
    sm = mv.smc(d["C"])
    out = s.output
    out.write("\nSquared multiple correlations of variables with all other variables\n\n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 10}\n", "text")
    out.write(f"{RULE}{'Variable':>12} |{'smc':>9} \n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 10}\n", "text")
    for n_, v in zip(d["names"], sm):
        out.write(f"{RULE}{_abbr(n_):>12} |{v:>9.4f} \n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 10}\n", "text")
    s.r = {"smc": matrix(sm.reshape(-1, 1), d["names"], ["smc"])}


def _estat_kmo(s, est, rest):
    d = _mvd(est)
    k_, ov = mv.kmo(d["C"])
    out = s.output
    out.write("\nKaiser-Meyer-Olkin measure of sampling adequacy\n\n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 10}\n", "text")
    out.write(f"{RULE}{'Variable':>12} |{'kmo':>9} \n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 10}\n", "text")
    for n_, v in zip(d["names"], k_):
        out.write(f"{RULE}{_abbr(n_):>12} |{v:>9.4f} \n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 10}\n", "text")
    out.write(f"{RULE}{'Overall':>12} |{ov:>9.4f} \n", "text")
    out.write(f"{RULE}{'-' * 13}+{'-' * 10}\n", "text")
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
    s.output.write("\n", "text")
    s.r = {"Cov": matrix(cov, names, names), "Corr": matrix(corr, names, names)}


def _estat_summarize(s, est, rest):
    d = _mvd(est)
    out = s.output
    ds = s.data
    mask = ds.esample
    out.write(f"\nEstimation sample {est.cmd}{'':>26}Number of obs = {d['N']:>10}\n\n", "text")
    out.write(f"  {'Variable':>13} |{'Mean':>11}{'Std. Dev.':>12}{'Min':>11}{'Max':>11}\n", "text")
    out.write(f"  {'-' * 13}-+{'-' * 46}\n", "text")
    for n_, m_, sd in zip(d["names"], d["means"], d["sds"]):
        x = ds.get(n_).data[mask]
        out.write(f"  {_abbr(n_):>13} |{gnum(m_, '%9.0g'):>11}{gnum(sd, '%9.0g'):>12}"
                  f"{gnum(x.min(), '%9.0g'):>11}{gnum(x.max(), '%9.0g'):>11}\n", "text")
    out.write("\n", "text")


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
    if o.get("fitted"):
        _print_matrix(s.output, "Fitted correlation matrix", names, F)
        s.r = {"Fit": matrix(F, names, names)}
    else:
        Rm = d["C"] - F
        _print_matrix(s.output, "Raw residuals of correlations (observed-fitted)", names, Rm)
        iu = np.triu_indices(len(names), 1)
        s.r = {"rmse_off": float(math.sqrt((Rm[iu] ** 2).mean())), "Rm": matrix(Rm, names, names)}
    s.output.write("\n", "text")


def _estat_structure(s, est, rest):
    d = _mvd(est)
    if d["kind"] != "factor":
        raise StataError(321, "subcommand estat structure is unrecognized")
    L = d["rot"]["L"] if d["rot"] else d["L"]
    Phi = d["rot"]["Phi"] if d["rot"] else np.eye(L.shape[1])
    _load_table(s.output, "Structure matrix: correlations between variables and common factors", d["names"], L @ Phi,
                "Factor", None, "")
    s.output.write("\n", "text")
    s.r = {"Lstructure": matrix(L @ Phi, d["names"], [f"Factor{j + 1}" for j in range(L.shape[1])])}


def _estat_common(s, est, rest):
    d = _mvd(est)
    if d["kind"] != "factor":
        raise StataError(321, "subcommand estat common is unrecognized")
    Phi = d["rot"]["Phi"] if d["rot"] else np.eye(d["k"])
    cols = [f"Factor{j + 1}" for j in range(Phi.shape[0])]
    _print_matrix(s.output, "Correlation matrix of the common factors", cols, Phi, cols)
    s.output.write("\n", "text")
    s.r = {"Phi": matrix(Phi, cols, cols)}


# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

def _pred_scores(s, est, how: str = "regression"):
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
        X = np.column_stack([ds.get(n).data.astype(float) for n in names])
        ok = mask & np.all(X < SYS, axis=1)
        Z = np.where(ok[:, None], X - d["means"], np.nan)
        if not (d["kind"] == "pca" and d["cov"]):
            Z = Z / d["sds"]
        if d["kind"] == "pca":
            if d["rot"]:
                Lr = d["rot"]["L"]
                Sc = Z @ (Lr @ np.linalg.inv(Lr.T @ Lr))      # VERIFY: scoring of rotated components
            else:
                Sc = Z @ d["V"][:, :k]
        else:
            L = d["rot"]["L"] if d["rot"] else d["L"]
            Phi = d["rot"]["Phi"] if d["rot"] else np.eye(L.shape[1])
            if how == "bartlett":
                Pi = np.diag(1.0 / d["psi"])
                W = Pi @ L @ np.linalg.inv(L.T @ Pi @ L)
            else:
                W = np.linalg.solve(d["C"], L @ Phi)
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
        for lab, val in (("Average interitem " + ("correlation:" if std else "covariance:"), gnum(avg, "%9.0g")),
                         ("Number of items in the scale:", f"{k}"),
                         ("Scale reliability coefficient:", f"{res['alpha']:.4f}")):
            out.write(lab + val.rjust(42 - len(lab)) + "\n", "text")
        out.write("\n", "text")
    s.r = {"alpha": res["alpha"], "k": float(k), "rbar" if std else "cbar": avg}
    if o.get("detail"):
        sd = np.sqrt(np.diag(C))
        _print_matrix(out, "Interitem correlations", names, np.nan_to_num(C / np.outer(sd, sd)))
        out.write("\n", "text")
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
    rows = []
    for j in range(k):
        it = Xc[:, j]
        rest = tot - it
        itc = float(np.corrcoef(it, tot)[0, 1])
        irc = float(np.corrcoef(it, rest)[0, 1])
        idx = [i for i in range(k) if i != j]
        sub = C[np.ix_(idx, idx)]
        a_j = mv.alpha_from_cov(sub, std)["alpha"] if k > 2 else float("nan")
        inter = float(np.mean([R[j, i] if std else C[j, i] for i in idx]))
        rows.append((names[j], int(Nm[j, j]), "-" if sign[j] < 0 else "+", itc, irc, inter, a_j))
    avg = res["rbar"] if std else res["cbar"]
    lab = "correlation" if std else "covariance"
    out.write("\n", "text")
    out.write(f"{'average':>64}\n{'item-test':>38}{'item-rest':>13}{'interitem':>13}\n", "text")
    out.write(f"{'Item':<13}| {'Obs':>4}{'Sign':>5}{'correlation':>14}{'correlation':>13}{lab:>13}{'alpha':>8}\n", "text")
    out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")
    for n_, ob, sg, itc, irc, inter, a_j in rows:
        out.write(f"{n_[:12]:<13}| {ob:>4}{sg:>5}{itc:>14.4f}{irc:>13.4f}{inter:>13.4f}{a_j:>8.4f}\n", "text")
    out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")
    out.write(f"{'Test scale':<13}|{'':>37}{avg:>13.4f}{res['alpha']:>8.4f}\n", "text")
    out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")


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
