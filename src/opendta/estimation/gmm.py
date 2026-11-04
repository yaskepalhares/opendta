"""gmm ([R] gmm): generalized method of moments with interactive or program moment equations,
one-step, two-step and iterated estimators and the weighting matrices of the manual."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.varlist import resolve_name
from ..lang.syntax import match_options, parse_standard
from .mleq import aux_row
from .models2 import g
from .nonlin import ProgFunc, SubExpr, _int_option, _num_jac, _options_tail, _param_names, _program, _split_equations
from .regress import cmdline, comma
from .results import Estimates, coef_table, post

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def cmd_gmm(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "gmm", args)
        return
    ds = s.data
    eqnames: list[str] = []
    evalprog = ""
    pm = re.match(r"^\s*([A-Za-z_]\w*)(.*)$", args, re.S)
    if pm and not args.lstrip().startswith("("):
        # moment-evaluator program: called as `prog residual_vars if touse, at(b)` ([R] gmm)
        prog = _program(s, pm.group(1))
        tail = pm.group(2).strip()
        ot = _options_tail(tail)
        mq = re.search(r"\beq\w*\s*\(([^)]*)\)", ot)
        if mq and not re.match(r"neq", ot[mq.start():mq.start() + 3]):
            eqnames = mq.group(1).split()
        neq = _int_option(ot, "neq") or len(eqnames)
        if not neq:
            raise StataError(198, "option nequations() or equations() required")   # VERIFY
        eqnames = eqnames or [str(q) for q in range(1, neq + 1)]
        pf = ProgFunc(s, prog, neq, _param_names(ot), [], noisy_first=True)
        ses = [pf.eq(q) for q in range(neq)]
        eqs = ["" for _ in range(neq)]
        evalprog = prog
        tail = re.sub(r"\beq\w*\s*\([^)]*\)", " ", tail) if mq else tail
    else:
        eqs, tail = _split_equations(args)
        if not eqs:
            raise StataError(198, "equation required")   # VERIFY
        ses = [SubExpr(s, body) for body in eqs]
    names: list[str] = []
    init: dict = {}
    for se in ses:
        for nm in se.names:
            if nm not in names:
                names.append(nm)
            if se.init.get(nm) is not None:
                init[nm] = se.init[nm]
    p = parse_standard(tail)
    # instruments() may repeat, so it is cut from the option text before matching the others
    opt_text = p.options
    inst_specs = re.findall(r"instruments\s*\(((?:[^()]|\([^()]*\))*)\)", opt_text)
    opt_text = re.sub(r"instruments\s*\(((?:[^()]|\([^()]*\))*)\)", " ", opt_text)
    # derivative() may repeat too: one analytic derivative of an equation's residual per parameter.
    # The plural derivatives() belongs to moment-evaluator programs: with expressions Stata refuses
    # it (compat 1920)
    plural = re.findall(r"derivatives\s*\((?:[^()]|\([^()]*\))*\)", opt_text)
    if plural and not evalprog:
        raise StataError(198, " ".join(plural) + " not allowed with interactive version")
    opt_text = re.sub(r"derivatives\s*\((?:[^()]|\([^()]*\))*\)", " ", opt_text)
    deriv_specs = re.findall(r"\bderivative\s*\(((?:[^()]|\([^()]*\))*)\)", opt_text)
    opt_text = re.sub(r"\bderivative\s*\(((?:[^()]|\([^()]*\))*)\)", " ", opt_text)
    o = match_options(opt_text, {"onestep": 3, "twostep": 3, "igmm": 4, "winitial": 4, "wmatrix": 2, "vce": 3,
                                 "level": 1, "nolog": 4, "log": 3, "eps": 3, "iterate": 4, "from": 3, "robust": 1,
                                 "cluster": 2, "noheader": 6, "notable": 5, "nocommonesample": 5, "quickderivatives": 5,
                                 "derivatives": 3, "variables": 3, "hasconstant": 3, "parameters": 5, "nequations": 3,
                                 "haslfderivatives": 5, "conv_maxiter": 5, "technique": 4, "nparameters": 5,
                                 "igmmiterate": 8, "igmmeps": 7, "igmmweps": 8}) \
        if opt_text.strip() else {}
    if (o.get("parameters") or o.get("nparameters")) and not evalprog:
        # parameters() and nparameters() belong to moment-evaluator programs (compat 1923)
        mt = re.search(r"\bn?param\w*\s*\([^)]*\)", opt_text)
        raise StataError(198, f"{mt.group(0) if mt else 'parameters()'} not allowed with interactive version")
    mask = touse(s, p)
    for se in ses:
        for v in se.variables():
            mask &= ds.get(v).data < SYS
    # instruments per equation (eqno: list), or one list for all
    inst_by_eq: dict[int, tuple[list[str], bool]] = {}
    for spec in inst_specs:
        mm = re.match(r"^\s*(?:(\w+)\s*:)?\s*(.*)$", spec, re.S)
        if mm.group(1) and not mm.group(1).isdigit():
            if mm.group(1) not in eqnames:
                raise StataError(198, f"equation {mm.group(1)} not found")   # VERIFY
            eqno = eqnames.index(mm.group(1))
        else:
            eqno = int(mm.group(1)) - 1 if mm.group(1) else None
        body = mm.group(2)
        noc = False
        if "," in body:
            body, opt = body.split(",", 1)
            noc = "noconstant".startswith(opt.strip()) and bool(opt.strip())
        vs = [resolve_name(ds, v) for v in body.split()]
        for q in range(len(ses)) if eqno is None else [eqno]:
            inst_by_eq[q] = (vs, noc)
    for q in range(len(ses)):
        if q not in inst_by_eq:
            raise StataError(198, "option instruments() required")   # VERIFY
    for vs, _ in inst_by_eq.values():
        for v in vs:
            mask &= ds.get(v).data < SYS
    theta0 = np.array([init.get(nm, 0.0) for nm in names])
    if o.get("from"):
        for mm in re.finditer(r"(\w+)\s+([-+0-9.eE]+)", str(o["from"])):
            if mm.group(1) in names:
                theta0[names.index(mm.group(1))] = float(mm.group(2))
    for se in ses:
        with np.errstate(all="ignore"):
            f0 = se.evaluate({nm: theta0[names.index(nm)] for nm in se.names})
        mask &= (f0 < SYS) & np.isfinite(f0)
    n = int(mask.sum())
    if n == 0:
        raise StataError(2000, "no observations")
    m = len(ses)
    Zs, znames = [], []
    for q in range(m):
        vs, noc = inst_by_eq[q]
        cols = [ds.get(v).data[mask].astype(float) for v in vs] + ([] if noc else [np.ones(n)])
        Zs.append(np.column_stack(cols))
        znames.append(" ".join(vs + ([] if noc else ["_cons"])))
    offs = np.concatenate([[0], np.cumsum([Z.shape[1] for Z in Zs])]).astype(int)
    n_mom = int(offs[-1])
    k = len(names)
    if n_mom < k:
        raise StataError(481, "equation not identified; number of moments is less than the number of parameters")   # VERIFY
    out = s.output
    wm = str(o.get("wmatrix", "robust")).split()
    wkind = wm[0].lower()[:3]
    cl = None
    clustvar = ""
    if wkind == "clu":
        clustvar = resolve_name(ds, wm[1])
        cv = ds.get(clustvar).data[mask]
        _, cl = np.unique(cv, return_inverse=True)
    if o.get("cluster"):
        clustvar = resolve_name(ds, str(o["cluster"]))
        _, cl = np.unique(ds.get(clustvar).data[mask], return_inverse=True)
        wkind = "clu"

    def resid(t):
        cols = []
        with np.errstate(all="ignore"):
            for se in ses:
                cols.append(se.evaluate({nm: t[names.index(nm)] for nm in se.names})[mask])
        return np.column_stack(cols)

    def moments(t):
        e = resid(t)
        return np.column_stack([Zs[q] * e[:, [q]] for q in range(m)])

    def mbar(t):
        return moments(t).mean(axis=0)

    # derivatives(): (equation, parameter) -> expression for d residual / d parameter
    derivs: dict[tuple[int, str], SubExpr] = {}
    for spec in deriv_specs:
        mm = re.match(r"^\s*(?:(\d+)\s*:|(\w*)\s*/)?\s*(\w+)\s*=\s*(.+)$", spec, re.S)
        if not mm or mm.group(3) not in names:
            raise StataError(198, "invalid derivatives() specification")   # VERIFY
        eqno = int(mm.group(1)) - 1 if mm.group(1) else 0
        if not 0 <= eqno < m:
            raise StataError(198, "invalid derivatives() specification")   # VERIFY
        derivs[(eqno, mm.group(3))] = SubExpr(s, mm.group(4))
    if derivs:
        # every parameter of an equation needs its derivative (VERIFY: Stata's rule when only some are given)
        for q, se in enumerate(ses):
            for nm in se.names:
                if (q, nm) not in derivs:
                    raise StataError(198, f"derivatives() required for parameter {nm} of equation {q + 1}")   # VERIFY

    def jac_mbar(t):
        if not derivs:
            return _num_jac(mbar, t)
        J = np.zeros((n_mom, k))
        for (q, nm), de in derivs.items():
            with np.errstate(all="ignore"):
                d = de.evaluate({a: t[names.index(a)] for a in de.names})[mask]
            J[offs[q]:offs[q + 1], names.index(nm)] = (Zs[q] * d[:, None]).mean(axis=0)
        return J

    def meat(t, kind):
        G_ = moments(t)
        if kind == "clu" and cl is not None:
            u = np.zeros((cl.max() + 1, G_.shape[1]))
            np.add.at(u, cl, G_)
            return u.T @ u / n
        return G_.T @ G_ / n

    def step(t, W, label):
        out.write(f"\n{label}\n" if label else "", "text")
        q_old = None
        for it in range(100):
            gb = mbar(t)
            G = jac_mbar(t)
            q = float(gb @ W @ gb)
            if not o.get("nolog"):
                out.write(f"Iteration {it}:   GMM criterion Q(b) = {g(q, '%10.0g'):>10}\n", "text")
            if (q_old is not None and abs(q_old - q) < 1e-10 * (1 + abs(q))) or q < 1e-30:
                break
            d = np.linalg.solve(G.T @ W @ G + 1e-12 * np.eye(k), G.T @ W @ gb)
            lam, tn = 1.0, t
            while lam > 1e-8:
                tn = t - lam * d
                gn = mbar(tn)
                qn = float(gn @ W @ gn)
                if np.isfinite(qn) and qn <= q + 1e-12:
                    break
                lam /= 2
            q_old, t = q, tn
        return t

    # initial weight matrix: unadjusted, Z'Z/n blocks scaled by the residual covariance
    e0 = resid(theta0)
    S0 = np.zeros((n_mom, n_mom))
    sig0 = e0.T @ e0 / n
    for a in range(m):
        for b in range(m):
            S0[offs[a]:offs[a + 1], offs[b]:offs[b + 1]] = sig0[a, b] * (Zs[a].T @ Zs[b]) / n
    winit = str(o.get("winitial", "unadjusted")).strip().lower()
    if winit.startswith("iden"):
        W1 = np.eye(n_mom)
        winit_label = "Identity"
    else:
        W1 = np.linalg.inv(S0 + 1e-12 * np.eye(n_mom))
        winit_label = "Unadjusted"
    theta = step(theta0, W1, "Step 1")
    n_steps = 0
    if o.get("onestep"):
        W = W1
        wm_label = winit_label
    else:
        S = meat(theta, wkind) if wkind in ("rob", "clu") else None
        if S is None:
            e1 = resid(theta)
            sg = e1.T @ e1 / n
            S = np.zeros((n_mom, n_mom))
            for a in range(m):
                for b in range(m):
                    S[offs[a]:offs[a + 1], offs[b]:offs[b + 1]] = sg[a, b] * (Zs[a].T @ Zs[b]) / n
        W = np.linalg.inv(S)
        theta = step(theta, W, "Step 2")
        wm_label = {"rob": "Robust", "una": "Unadjusted", "clu": f"Cluster ({clustvar})"}[wkind]
        if o.get("igmm"):
            # iterated GMM: re-estimate the weight matrix until the parameters (igmmeps()) or the
            # weight matrix (igmmweps()) stop changing (VERIFY the defaults and the log)
            eps = float(o["igmmeps"]) if o.get("igmmeps") else 1e-6
            weps = float(o["igmmweps"]) if o.get("igmmweps") else 1e-6
            maxit = int(float(o["igmmiterate"])) if o.get("igmmiterate") else 300
            for it_ in range(3, maxit + 3):
                if wkind in ("rob", "clu"):
                    S = meat(theta, wkind)
                else:
                    e1 = resid(theta)
                    sg = e1.T @ e1 / n
                    S = np.zeros((n_mom, n_mom))
                    for a in range(m):
                        for b in range(m):
                            S[offs[a]:offs[a + 1], offs[b]:offs[b + 1]] = sg[a, b] * (Zs[a].T @ Zs[b]) / n
                W_new = np.linalg.inv(S)
                t_new = step(theta, W_new, f"Step {it_}")
                dt = float(np.max(np.abs(t_new - theta)) / (1 + np.max(np.abs(theta))))
                dw = float(np.max(np.abs(W_new - W)) / (1 + np.max(np.abs(W))))
                theta, W = t_new, W_new
                n_steps = it_
                if dt < eps or dw < weps:
                    if not o.get("nolog"):
                        out.write("note: iterative GMM parameter vector converged\n" if dt < eps else
                                  "note: iterative GMM weight matrix converged\n", "text")   # compat 1920 (VERIFY the second)
                    break
    gb = mbar(theta)
    G = jac_mbar(theta)
    Q = float(gb @ W @ gb)
    Ainv = np.linalg.inv(G.T @ W @ G)
    vkind = "clu" if (cl is not None) else ("rob" if wkind == "rob" or str(o.get("vce", "")).startswith("rob")
                                              or o.get("robust") else "una")
    if vkind == "una":
        V = Ainv / n
    else:
        V = Ainv @ (G.T @ W @ meat(theta, vkind) @ W @ G) @ Ainv / n
    V = (V + V.T) / 2
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    est = Estimates("gmm", "", list(names), theta, V, n, "z", None)
    est.level = level
    est.vcetype = "Robust" if vkind in ("rob", "clu") else ""
    est.clustvar, est.n_clust = clustvar, (int(cl.max() + 1) if cl is not None else None)
    est.eqnames = ["/"] * k
    est.rows = [aux_row("/" + nm, float(theta[i]), float(np.sqrt(max(V[i, i], 0))), level) for i, nm in enumerate(names)]
    est.extra = dict(N=n, n_par=k, n_mom=n_mom, wm_label=wm_label, winit=winit_label, Q=Q, znames=znames,
                     onestep=bool(o.get("onestep")), eqn=eqnames if eqnames and not eqnames[0].isdigit() else None)
    est.sexps, est.pnames = ses, names

    def display(s_, e, *, header=True, table=True, level=None, **_kw):
        o_ = s_.output
        x = e.extra
        if header:
            o_.write("\nGMM estimation \n\n", "text")
            o_.write(f"Number of parameters = {x['n_par']:>3}\n", "text")
            o_.write(f"Number of moments    = {x['n_mom']:>3}\n", "text")
            o_.write(f"{'Initial weight matrix: ' + x['winit']:<50}{'Number of obs':<16}=", "text")
            o_.write(f"{comma(x['N']):>11}\n", "result")
            if not x["onestep"]:
                o_.write(f"GMM weight matrix:     {x['wm_label']}\n", "text")
            o_.write("\n", "text")
        if not table:
            return
        coef_table(s_, "", e.rows, stat="z", df=None, level=level or e.level, vcetype=e.vcetype)
        for i, zn in enumerate(x["znames"], 1):
            label = x["eqn"][i - 1] if x.get("eqn") else str(i)    # equations() names (compat 1923)
            o_.write(f"Instruments for equation {label}: {zn}\n", "text")
    est.display = display
    # e() as Stata lists it (compat 1920); each parameter is an equation of its own (b1:_cons)
    eqn = eqnames or [str(q) for q in range(1, m + 1)]
    pq = [[nm for nm in names if nm in se.names] for se in ses]
    scalars = [("rank", float(np.linalg.matrix_rank(V))), ("N", n), ("Q", Q), ("J", n * Q), ("J_df", n_mom - k)]
    if not evalprog:
        scalars += [(f"k_{q + 1}", float(len(pq[q]))) for q in range(m)]
    scalars += [("converged", 1.0), ("has_xtinst", 0.0), ("type", 2.0 if eqnames else 1.0),   # VERIFY type 2
                ("n_eq", float(m)), ("k", float(k)), ("n_moments", float(n_mom)), ("k_eq", float(k)),
                ("k_aux", float(k)), ("k_eq_model", 0.0)]
    if n_steps:
        scalars.append(("ic", float(n_steps)))
    estimator = "onestep" if o.get("onestep") else ("igmm" if o.get("igmm") else "twostep")
    vce_name = {"rob": "robust", "clu": "cluster", "una": "unadjusted"}[vkind]
    macros = [("cmdline", cmdline("gmm", args)), ("cmd", "gmm"), ("estat_cmd", "gmm_estat"), ("predict", "gmm_p"),
              ("marginsnotok", "_ALL"), ("eqnames", " ".join(eqn)), ("technique", "gn"), ("winit", winit_label),
              ("estimator", estimator), ("wmatrix", wm_label.lower().split()[0] if not o.get("onestep") else
                                         winit_label.lower()),
              ("vce", vce_name)]
    if vkind != "una":
        macros.append(("vcetype", "Robust"))
    if clustvar:
        macros.append(("clustvar", clustvar))
    if evalprog:
        # compat 1923: the program, the parameter names as given, no per-equation parameters
        macros += [("evalprog", evalprog), ("params", " ".join(names))]
    else:
        macros.append(("params", " ".join(f"{nm}:_cons" for nm in names)))
    for q in range(m):
        macros.append((f"inst_{q + 1}", znames[q]))
        if not evalprog:
            macros.append((f"params_{q + 1}", " ".join(f"{nm}:_cons" for nm in pq[q])))
        if getattr(ses[q], "text", ""):
            sx = re.sub(r"\s+", "", ses[q].text)
            sx = re.sub(r"\{(\w+)(=[^}]*)?\}", lambda mm_: "{" + mm_.group(1) + ":_cons}", sx)   # VERIFY {lc:} terms
            macros.append((f"sexp_{q + 1}", sx))
    macros.append(("properties", "b V"))
    from ..commands.matrix import Matrix
    cn = [f"{nm}:_cons" for nm in names]

    def _mat(a, rows, cols):
        mt = Matrix(np.atleast_2d(a), rows, [c.split(":")[1] for c in cols])
        mt.coleq = [c.split(":")[0] for c in cols]
        return mt
    S_fin = meat(theta, vkind if vkind != "una" else "rob")
    mats = [("init", _mat(theta0.reshape(1, -1), ["y1"], cn)), ("V_modelbased", _mat(Ainv / n, cn, cn)),
            ("S", Matrix(S_fin)), ("W", Matrix(W))]
    post(s, est, scalars, macros, mask, extra_mats=mats)
    display(s, est, header=not o.get("noheader"), table=not o.get("notable"))
