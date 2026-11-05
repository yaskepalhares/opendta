"""rotate after pca and factor, and rotatemat ([MV] rotate, [MV] rotatemat).

The rotation criteria and the gradient projection algorithm are in core/rotation.py; varimax and
promax keep the dedicated routines of core/multivariate.py, which the manual says ignore the
optimize options."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core import multivariate as mv
from ..core import rng
from ..core import rotation as R
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options
from .mvcommon import chi2_p, fill_e, g as gnum, matrix
from .mvfactor import _eig_table, _f4, _factor_head, _header, _lr_line, _right, _table, _uniq_extra, _unexpl_extra
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

# ---------------------------------------------------------------------------
# options shared by rotate and rotatemat
# ---------------------------------------------------------------------------

_METHODS = ("varimax", "vgpf", "quartimax", "equamax", "parsimax", "entropy", "tandem1", "tandem2", "promax",
            "oblimin", "cf", "bentler", "oblimax", "quartimin", "target", "partial")
_COMMON = {"varimax": 3, "vgpf": 4, "quartimax": 6, "equamax": 3, "parsimax": 3, "entropy": 7, "tandem1": 7,
           "tandem2": 7, "promax": 3, "oblimin": 3, "cf": 2, "bentler": 7, "oblimax": 7, "quartimin": 6,
           "target": 6, "partial": 7, "orthogonal": 4, "oblique": 3, "normalize": 4, "blanks": 5, "format": 3,
           "noloading": 5, "norotation": 6, "iterate": 4, "log": 3, "trace": 3, "tolerance": 3,
           "gtolerance": 4, "ltolerance": 4, "protect": 4, "maxstep": 4, "init": 4, "random": 4}
_ROT = dict(_COMMON, factors=6, components=4, norotated=5, detail=3)
_ROTMAT = dict(_COMMON, nodisplay=5, matname=7, colnames=8)


def _spec(s, o: dict, p: int, k: int) -> dict:
    """The rotation requested by the options: method, its argument, class and optimize options."""
    from .matrix import get_matrix
    methods = [m for m in _METHODS if o.get(m)]
    if len(methods) > 1:
        raise StataError(198, "only one rotation method may be specified")                   # VERIFY
    method = methods[0] if methods else "varimax"
    if o.get("orthogonal") and o.get("oblique"):
        raise StataError(198, "options orthogonal and oblique may not be combined")          # VERIFY
    oblique = bool(o.get("oblique")) or method == "promax"
    if method == "promax" and o.get("orthogonal"):
        raise StataError(198, "promax rotation is oblique; option orthogonal not allowed")   # VERIFY
    if oblique and method in ("varimax",) + R.ORTHOGONAL_ONLY:
        raise StataError(198, f"{method} rotation is orthogonal only; option oblique not allowed")   # VERIFY
    arg, carg = None, ""
    val = o.get(method)
    if method == "promax":
        arg = float(val) if isinstance(val, str) else 3.0
        if arg <= 1:
            raise StataError(198, "promax() power must be greater than 1")                  # VERIFY
        carg = f"{arg:g}"
    elif method == "oblimin":
        arg = float(val) if isinstance(val, str) else 0.0
        carg = f"{arg:g}"
    elif method == "cf":
        if not isinstance(val, str):
            raise StataError(198, "option cf() requires an argument")                       # VERIFY
        arg = float(val)
        if not 0 <= arg <= 1:
            raise StataError(198, "cf() must be in [0,1]")                                  # VERIFY
        carg = f"{arg:g}"
    elif method in ("target", "partial"):
        names = str(val).split() if isinstance(val, str) else []
        if len(names) != (1 if method == "target" else 2):
            raise StataError(198, f"option {method}() incorrectly specified")                # VERIFY
        mats = [get_matrix(s, n_).data for n_ in names]
        if any(m_.shape != (p, k) for m_ in mats):
            raise StataError(503, "conformability error")
        arg = mats[0] if method == "target" else (mats[0], mats[1])
        carg = " ".join(names)
    if o.get("init") and o.get("random"):
        raise StataError(198, "options init() and random may not be combined")              # VERIFY
    T0 = get_matrix(s, str(o["init"])).data if isinstance(o.get("init"), str) else None
    if T0 is not None and T0.shape != (k, k):
        raise StataError(503, "conformability error")

    def num(name, default):
        return float(o[name]) if isinstance(o.get(name), str) else default
    return {"method": method, "arg": arg, "carg": carg, "oblique": oblique, "norm": bool(o.get("normalize")),
            "iterate": int(num("iterate", 1000)), "tol": num("tolerance", 1e-6), "gtol": num("gtolerance", 1e-6),
            "ltol": num("ltolerance", 1e-6), "maxstep": int(num("maxstep", 20)),
            "protect": int(num("protect", 0)), "random": bool(o.get("random")), "T0": T0}


def _ctitle(spec: dict) -> tuple[str, str]:
    """(criterion, ctitle) as rotatemat stores them in r(criterion) and r(ctitle)."""
    m = spec["method"]
    if m in ("promax", "oblimin", "cf"):
        return m, f"{m}({spec['carg']})"
    if m == "vgpf":
        return "varimax", "varimax"                                                          # VERIFY
    return m, m


def _random_T(k: int, oblique: bool) -> np.ndarray:
    """A random start for random/protect(): normal deviates from Stata's uniform stream (VERIFY:
    Stata's draw of the random orthogonal or normal matrix is not documented)."""
    u = np.clip(rng.RNG.uniform(k * k), 1e-16, 1 - 1e-16)
    return R.random_start(k, oblique, st.norm.ppf(u).reshape(k, k))


def rotate_matrix(s, A: np.ndarray, spec: dict) -> dict:
    """Rotate A as the spec says. Returns {"L", "T" (manual's T: L = A T, or A inv(T') when
    oblique), "f", "iter", "trials"}. Kaiser normalization is applied around the rotation; the
    columns get positive sums and decreasing sums of squares (core/rotation.orient)."""
    m, obl = spec["method"], spec["oblique"]
    B, h = R.kaiser(A) if spec["norm"] else (A, np.ones(A.shape[0]))
    trials = []
    if m == "varimax":
        L, T, _ = mv.orthomax(B, 1.0)
        f, it = SYS, SYS
    elif m == "promax":
        L, Tp, _ = mv.promax(B, spec["arg"])
        T = np.linalg.inv(Tp).T
        f, it = SYS, SYS
    else:
        opts = dict(oblique=obl, maxiter=spec["iterate"], tol=spec["tol"], gtol=spec["gtol"], ltol=spec["ltol"],
                    maxstep=spec["maxstep"])
        name = "oblimin" if m == "oblimin" else m
        k = A.shape[1]
        if spec["protect"] > 0:
            best = None
            for _ in range(spec["protect"]):
                r_ = R.gpf(B, name, spec["arg"], T0=_random_T(k, obl), **opts)
                trials.append(r_)
                if best is None or r_["f"] < best["f"]:
                    best = r_
            r_ = best
        else:
            T0 = spec["T0"] if spec["T0"] is not None else (_random_T(k, obl) if spec["random"] else None)
            r_ = R.gpf(B, name, spec["arg"], T0=T0, **opts)
        L, T, f, it = r_["L"], r_["T"], r_["f"], r_["iter"]
    L, T = R.orient(L, T)
    return {"L": L * h[:, None], "T": T, "f": f, "iter": it, "trials": trials}


def _trials(out, trials) -> None:
    """protect(): one line per trial (VERIFY: spacing; the manual's example shows the values in %9.0g)."""
    for i, r_ in enumerate(trials, start=1):
        out.write(f"Trial {i:>2} : min criterion {gnum(r_['f']).strip()}\n", "text")


# ---------------------------------------------------------------------------
# rotate
# ---------------------------------------------------------------------------

@command("rotate", "rot")
def cmd_rotate(s: "Session", args: str) -> None:
    est = getattr(s, "_estimates", None)
    d = (getattr(est, "extra", None) or {}).get("mv") if est is not None else None
    if d is None or d["kind"] not in ("pca", "factor"):
        raise StataError(301, "last estimates not found")
    o = match_options(args.strip().lstrip(","), _ROT) if args.strip().lstrip(",").strip() else {}
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
    spec = _spec(s, o, pdim, kk)
    method, norm = spec["method"], spec["norm"]
    Phi = np.eye(kk)
    fmin = SYS
    if method in ("varimax", "quartimax", "equamax", "parsimax") and not spec["protect"]:
        gamma = {"varimax": 1.0, "quartimax": 0.0, "equamax": kk / 2.0, "parsimax": pdim * (kk - 1) / (pdim + kk - 2.0)}[method]
        Lr, T, _ = mv.orthomax(A, gamma, normalize=norm)
        desc = f"orthogonal {method}"
    elif method == "promax":
        Lr, T, Phi = mv.promax(A, spec["arg"], normalize=norm)
        desc = "oblique promax"
    elif method in ("oblimin", "quartimin") and spec["oblique"] and not spec["protect"] and spec["T0"] is None \
            and not spec["random"]:
        gamma = spec["arg"] if method == "oblimin" else 0.0
        Lr, T, Phi, _ = mv.oblimin(A, gamma, normalize=norm)
        T = np.linalg.inv(T).T                      # the rotation matrix shown is the one with Lr = A T
        desc = f"oblique {method}"
        fmin = mv._oblimin_crit(Lr, gamma)[0]
    else:
        r_ = rotate_matrix(s, A, spec)
        _trials(s.output, r_["trials"])
        Lr, Tm = r_["L"], r_["T"]
        if spec["oblique"]:
            T, Phi = np.linalg.inv(Tm).T, Tm.T @ Tm
        else:
            T = Tm
        # rotate sorts the columns by decreasing explained variance ([MV] rotate, technical note)
        ev_ = ((Lr @ Phi) ** 2).sum(0) if kind == "factor" and spec["oblique"] else ((As @ T) ** 2).sum(0)
        order = np.argsort(-ev_, kind="stable")
        Lr, T, Phi = Lr[:, order], T[:, order], Phi[np.ix_(order, order)]
        fmin = r_["f"]
        desc = f"{'oblique' if spec['oblique'] else 'orthogonal'} {_ctitle(spec)[0]}"
    crit, ctitle = _ctitle(spec)
    oblique = spec["oblique"]
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
    # (the value of the rotation criterion) is stored only for the criteria that have one
    rs = ([("r_fmin", fmin)] if fmin != SYS else []) + [("r_f", kk)]
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
# rotatemat
# ---------------------------------------------------------------------------

def _fmt_cell(fmt: str):
    return lambda x: format_value(x, fmt, pad=False)


@command("rotatemat", "rotatemat")
def cmd_rotatemat(s: "Session", args: str) -> None:
    """rotatemat matname [, options] ([MV] rotatemat): rotation of a Stata matrix. The layout
    follows the manual's examples; the spacing of the header lines is VERIFY."""
    from .matrix import get_matrix
    body, _, opt = args.partition(",")
    name = body.strip()
    if not name or len(name.split()) != 1:
        raise StataError(198, "invalid syntax")
    M_ = get_matrix(s, name)
    A = M_.data
    if (np.abs(A) >= SYS).any():
        raise StataError(504, f"matrix {name} has missing values")
    p, k = A.shape
    o = match_options(opt, _ROTMAT) if opt.strip() else {}
    spec = _spec(s, o, p, k)
    r_ = rotate_matrix(s, A, spec)
    crit, ctitle = _ctitle(spec)
    cls = "oblique" if spec["oblique"] else "orthogonal"
    fmt = str(o["format"]).strip() if isinstance(o.get("format"), str) else "%9.5f"
    blanks = float(o["blanks"]) if isinstance(o.get("blanks"), str) else 0.0
    rows = M_.rownames
    cols = M_.colnames
    out = s.output
    if not o.get("nodisplay"):
        _trials(out, r_["trials"])
        label = str(o["matname"]) if isinstance(o.get("matname"), str) else f"{name}[{p},{k}]"
        out.write(f"\nRotation of {label}\n\n", "text")
        lines = [("Criterion", ctitle), ("Rotation class", cls), ("Kaiser normalization", "on" if spec["norm"] else "off")]
        if r_["f"] != SYS:
            lines += [("Criterion value", gnum(r_["f"]).strip()), ("Number of iterations", str(r_["iter"]))]
        for lab, val in lines:
            out.write(f"    {lab:<22}{val}\n", "text")
        cw = max(len(format_value(0.0, fmt)), 8)
        cname = str(o["colnames"]) if isinstance(o.get("colnames"), str) else "factors"
        if not o.get("noloading"):
            title = f"Rotated {cname}" + (f" (blanks represent abs()<{gnum(blanks).strip()})" if blanks else "")
            _table(out, title, rows, cols, r_["L"], _fmt_cell(fmt), stub="", blanks=blanks, cw=cw)
        if not o.get("norotation"):
            _table(out, f"{'Oblique' if spec['oblique'] else 'Orthogonal'} rotation", cols, cols, r_["T"],
                   _fmt_cell(fmt), stub="", cw=cw)
    res = {"AT": matrix(r_["L"], rows, cols), "T": matrix(r_["T"], cols, cols)}
    if r_["trials"]:
        res["fmin"] = matrix(np.array([[t_["f"] for t_ in r_["trials"]]]))
    res.update({"carg": spec["carg"], "normalization": "kaiser" if spec["norm"] else "none", "class": cls,
                "criterion": crit, "ctitle12": ctitle[:12], "ctitle": ctitle, "cmd": "rotatemat"})
    if r_["trials"]:
        res["nnconv"] = float(sum(not t_["converged"] for t_ in r_["trials"]))
    res.update({"rc": 0.0, "iter": float(r_["iter"]) if r_["iter"] != SYS else SYS, "f": float(r_["f"])})
    s.r = {k_: v_ for k_, v_ in res.items() if not (isinstance(v_, float) and v_ == SYS)}
