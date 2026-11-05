"""Canonical correlations ([MV] canon) and its estat subcommands ([MV] canon postestimation).

The statistics are in core/mvtests.py."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import multivariate as mv
from ..core import rotation as R
from ..core import mvtests as T
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from .mvcommon import matrix, numeric_names, post_e, sample
from .mvfactor import _g6
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
from .mvtests import _find_comma

# ---------------------------------------------------------------------------
# canon
# ---------------------------------------------------------------------------

def _split_sets(text: str) -> tuple[str, str, str]:
    t = text.strip()
    if not t.startswith("("):
        raise StataError(198, "invalid syntax")
    depth, parts, cur, i = 0, [], "", 0
    while i < len(t):
        ch = t[i]
        if ch == "(":
            depth += 1
            if depth == 1:
                cur = ""
                i += 1
                continue
        elif ch == ")":
            depth -= 1
            if depth == 0:
                parts.append(cur)
                i += 1
                continue
        if depth >= 1:
            cur += ch
        elif depth == 0 and ch not in " \t" and ch != ",":
            break
        i += 1
    rest = t[i:]
    if len(parts) != 2:
        raise StataError(198, "invalid syntax")
    return parts[0], parts[1], rest


@command("canon", "canon")
def cmd_canon(s: "Session", args: str) -> None:
    c = _find_comma(args)
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    set1, set2, rest = _split_sets(body)
    p = parse_standard("x " + rest + ("," + opt if opt else ""))
    o = match_options(p.options, {"lc": 2, "noconstant": 3, "stdcoef": 3, "test": 3, "level": 1}) if p.options.strip() else {}
    n1 = numeric_names(s, set1, 1)
    n2 = numeric_names(s, set2, 1)
    X, w, wtype, mask = sample(s, p, n1 + n2)
    p1, p2 = len(n1), len(n2)
    noconstant = bool(o.get("noconstant"))
    if noconstant:
        n_ = X.shape[0]
        sw = float(n_ if w is None else w.sum())
        W = np.ones(n_) if w is None else w
        C = (X * W[:, None]).T @ X / (sw - 1.0)
        means = np.zeros(p1 + p2)
    else:
        means, C, _, sw = mv.weighted_moments(X, w)
    N = int(round(sw))
    S11, S22, S12 = C[:p1, :p1], C[p1:, p1:], C[:p1, p1:]
    cc = T.canonical(S11, S22, S12)
    r, A1, A2 = cc["r"], cc["A1"], cc["A2"]
    # sign convention: the largest standardized coefficient (in absolute value) of each second-set
    # column is positive; it reproduces the auto example of [MV] canon and compat 1906 (VERIFY)
    k = len(r)
    sd1, sd2 = np.sqrt(np.diag(S11)), np.sqrt(np.diag(S22))
    for j in range(k):
        z2 = A2[:, j] * sd2
        if z2[np.argmax(np.abs(z2))] < 0:
            A1[:, j] *= -1
            A2[:, j] *= -1
    R11, R22 = S11 / np.outer(sd1, sd1), S22 / np.outer(sd2, sd2)
    R12 = S12 / np.outer(sd1, sd2)
    Z1, Z2 = A1 * sd1[:, None], A2 * sd2[:, None]
    shown1, shown2 = (Z1, Z2) if o.get("stdcoef") else (A1, A2)
    out = s.output
    out.write(f"\nCanonical correlation analysis{'':>22}Number of obs ={N:>11}\n", "text")
    cn = [str(j + 1) for j in range(k)]
    for title, nm, A_ in (("first", n1, shown1), ("second", n2, shown2)):
        out.write(f"\n{'Standardized' if o.get('stdcoef') else 'Raw'} coefficients for the {title} variable set\n\n", "text")
        out.write(f"{'':>17}|" + "".join(f"{c_:>9} " for c_ in cn) + "\n", "text")
        out.write(f"{'':>4}{'-' * 13}+{'-' * 10 * k}\n", "text")
        for i, n_ in enumerate(nm):
            out.write(f"{n_[:15]:>16} |" + "".join(f"{A_[i, j]:>9.4f} " for j in range(k)) + "\n", "text")
        out.write(f"{'':>4}{'-' * (14 + 10 * k)}\n", "text")
    rule = "-" * 76
    out.write(f"\n{rule}\nCanonical correlations:\n  " + "  ".join(f"{x:.4f}" for x in r) + "\n", "text")
    out.write(f"\n{rule}\nTests of significance of all canonical correlations\n\n", "text")
    out.write(f"{'Statistic':>34}{'df1':>9}{'df2':>9}{'F':>13}{'Prob>F':>11}\n", "text")
    lam = r ** 2 / (1.0 - r ** 2)
    v = N - 1.0 - p2 if not noconstant else N - p2
    stats = T.mv_statistics(lam, p1, p2, v)
    names_full = {"W": "Wilks' lambda", "P": "Pillai's trace", "L": "Lawley-Hotelling trace", "R": "Roy's largest root"}
    for r_ in stats:
        out.write(f"{names_full[r_['name']]:>22}{_g6(r_['stat']):>12}{_g6(r_['df1']):>9}{_g6(r_['df2']):>9}"
                  f"{_g6(r_['F']):>13}{r_['p']:>11.4f} {r_['flag']}\n", "text")
    out.write(f"{rule}\n{'':>28}e = exact, a = approximate, u = upper bound on F\n", "text")
    # e(b), e(V), e(stat_#) and e(stat_m) have the dimensions of Stata's; their contents are not documented (VERIFY)
    type_code = {"e": 1.0, "a": 2.0, "u": 3.0}
    cols = ["statistic", "df1", "df2", "F", "p", "type"]
    def _row(r_):
        return [r_["stat"], r_["df1"], r_["df2"], r_["F"], r_["p"], type_code[r_["flag"]]]
    stat_m = np.array([_row(r_) for r_ in stats])
    stat_j = []
    for j in range(k):
        w_ = T.mv_statistics(lam[j:], p1 - j, p2 - j, v)[0]
        stat_j.append(np.array([[r[j]] + _row(w_)[:5]]))
    b = np.concatenate([A1.T.ravel(), A2.T.ravel()]).reshape(1, -1)
    sc = [("df_r", N - 1.0), ("N", N), ("df", N - 1.0), ("n_cc", k), ("df1", stats[0]["df1"]),
          ("df2", stats[0]["df2"]), ("rank", b.shape[1])]
    mats = [("b", matrix(b, ["y1"], [f"c{i + 1}" for i in range(b.shape[1])])),
            ("V", matrix(np.zeros((b.shape[1], b.shape[1])), [f"c{i + 1}" for i in range(b.shape[1])],
                         [f"c{i + 1}" for i in range(b.shape[1])]))]
    mats += [(f"stat_{j + 1}", matrix(stat_j[j], ["stat"], cols)) for j in reversed(range(k))]
    mats += [("stat_m", matrix(stat_m, ["W", "P", "L", "R"], cols)),
             ("canload11", matrix(R11 @ Z1, n1, cn)), ("canload22", matrix(R22 @ Z2, n2, cn)),
             ("canload12", matrix(R12 @ Z2, n1, cn)), ("canload21", matrix(R12.T @ Z1, n2, cn)),
             ("rawcoef_var1", matrix(A1, n1, cn)), ("rawcoef_var2", matrix(A2, n2, cn)),
             ("stdcoef_var1", matrix(Z1, n1, cn)), ("stdcoef_var2", matrix(Z2, n2, cn)),
             ("ccorr", matrix(r.reshape(1, -1), ["ccorr"], cn)),
             ("corr_var1", matrix(R11, n1, n1)), ("corr_var2", matrix(R22, n2, n2)),
             ("corr_mixed", matrix(R12, n1, n2))]
    macros = [("cmdline", "canon " + args.strip()), ("marginsnotok", "_ALL"), ("cmd", "canon"),
              ("predict", "canon_p"), ("estat_cmd", "canon_estat"), ("properties", "b V")]
    post_e(s, "canon", sc, macros, mats, mask, N,
           estat={"correlations": (3, _estat_corr), "loadings": (3, _estat_loadings),
                  "rotatecompare": (7, _estat_rotatecompare), "rotate": (3, _estat_rotate)},
           extra={"canon": {"n1": n1, "n2": n2, "R11": R11, "R22": R22, "R12": R12, "A1": A1, "A2": A2,
                            "Z1": Z1, "Z2": Z2, "cn": cn, "rot": None}})


# ---------------------------------------------------------------------------
# estat after canon ([MV] canon postestimation)
# ---------------------------------------------------------------------------

def _fmt_opt(rest: str, extra=None) -> tuple[str, dict]:
    c = _find_comma(rest)
    opt = "" if c == -1 else rest[c + 1:]
    if c != -1 and rest[:c].strip() or c == -1 and rest.strip():
        raise StataError(198, "invalid syntax")
    o = match_options(opt, dict({"format": 1}, **(extra or {}))) if opt.strip() else {}
    return (str(o["format"]).strip() if isinstance(o.get("format"), str) else "%8.4f"), o


def _show(s, title: str, M_, rows, cols, fmt: str) -> None:
    """A titled matrix listed by matlist with the given format (VERIFY: border style)."""
    from .matrix import _store
    store = _store(s)
    tmp = "__canon_show"
    store[tmp] = matrix(M_, rows, cols)
    try:
        s.output.write(f"\n{title}\n", "text")
        s.interp._execute_expanded(f"matlist {tmp}, format({fmt}) border(bottom)")
    finally:
        store.pop(tmp, None)


def _canon(est) -> dict:
    return est.extra["canon"]


def _estat_corr(s, est, rest):
    fmt, _ = _fmt_opt(rest)
    c = _canon(est)
    _show(s, "Correlations for variable list 1", c["R11"], c["n1"], c["n1"], fmt)
    _show(s, "Correlations for variable list 2", c["R22"], c["n2"], c["n2"], fmt)
    _show(s, "Correlations between variable lists 1 and 2", c["R12"], c["n1"], c["n2"], fmt)
    s.r = {"corr_mixed": matrix(c["R12"], c["n1"], c["n2"]), "corr_var2": matrix(c["R22"], c["n2"], c["n2"]),
           "corr_var1": matrix(c["R11"], c["n1"], c["n1"])}


def _loadings(c) -> dict:
    return {"canload11": c["R11"] @ c["Z1"], "canload22": c["R22"] @ c["Z2"],
            "canload12": c["R12"] @ c["Z2"], "canload21": c["R12"].T @ c["Z1"]}


def _estat_loadings(s, est, rest):
    fmt, _ = _fmt_opt(rest)
    c = _canon(est)
    L = _loadings(c)
    cn = c["cn"]
    _show(s, "Canonical loadings for variable list 1", L["canload11"], c["n1"], cn, fmt)
    _show(s, "Canonical loadings for variable list 2", L["canload22"], c["n2"], cn, fmt)
    _show(s, "Correlation between variable list 1 and canonical variates from list 2", L["canload12"], c["n1"], cn, fmt)
    _show(s, "Correlation between variable list 2 and canonical variates from list 1", L["canload21"], c["n2"], cn, fmt)
    s.r = {k_: matrix(L[k_], c["n1"] if k_ in ("canload11", "canload12") else c["n2"], cn)
           for k_ in ("canload12", "canload21", "canload22", "canload11")}


_WHAT = {"rawcoefs": ("raw coefficients", "rawcoefs"), "stdcoefs": ("standardized coefficients", "stdcoefs"),
         "loadings": ("canonical loadings", "loadings")}


def _estat_rotate(s, est, rest):
    """Varimax rotation computed on the stacked canonical loadings (Cliff and Krus 1976), applied
    to the requested matrix ([MV] canon postestimation, Methods and formulas)."""
    fmt, o = _fmt_opt(rest, {"rawcoefs": 1, "stdcoefs": 1, "loadings": 1})
    which = [w_ for w_ in ("rawcoefs", "stdcoefs", "loadings") if o.get(w_)]
    if len(which) > 1:
        raise StataError(198, "only one of rawcoefs, stdcoefs and loadings may be specified")   # VERIFY
    what = which[0] if which else "rawcoefs"
    c = _canon(est)
    L = _loadings(c)
    k = len(c["cn"])
    if k < 2:
        raise StataError(198, "rotation requires at least 2 canonical correlations")         # VERIFY
    stacked = np.vstack([L["canload11"], L["canload22"]])
    # rotatemat's varimax, with its column orientation (core/rotation.orient)
    Lv, T, _ = mv.orthomax(stacked, 1.0)
    _, T = R.orient(Lv, T)
    base = {"rawcoefs": np.vstack([c["A1"], c["A2"]]), "stdcoefs": np.vstack([c["Z1"], c["Z2"]]),
            "loadings": stacked}[what]
    AT = base @ T
    rows = c["n1"] + c["n2"]
    out = s.output
    out.write("\n", "text")
    for lab, val in (("Criterion", "varimax"), ("Rotation class", "orthogonal"), ("Normalization", "none")):
        out.write(f"    {lab:<19}{val}\n", "text")                                        # VERIFY: spacing
    title = "Rotated " + _WHAT[what][0]
    _show(s, title, AT, rows, c["cn"], fmt)
    _show(s, "Rotation matrix", T, c["cn"], c["cn"], fmt)
    c["rot"] = {"AT": AT, "T": T, "what": what, "base": base, "rows": rows}
    s.r = {"T": matrix(T, c["cn"], c["cn"]), "AT": matrix(AT, rows, c["cn"]), "criterion": "varimax",
           "class": "orthogonal", "coefficients": _WHAT[what][1]}


def _estat_rotatecompare(s, est, rest):
    fmt, _ = _fmt_opt(rest)
    c = _canon(est)
    rot = c["rot"]
    if rot is None:
        raise StataError(321, "estat rotate must be run first")                               # VERIFY
    w_ = _WHAT[rot["what"]][0]
    _show(s, f"Rotated {w_} orthogonal varimax", rot["AT"], rot["rows"], c["cn"], fmt)
    _show(s, f"Unrotated {w_}", rot["base"], rot["rows"], c["cn"], fmt)
