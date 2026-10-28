"""ucm: unobserved-components models ([TS] ucm): trend, seasonal and stochastic-cycle components
estimated by maximum likelihood with the (diffuse) Kalman filter of ssm.py."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from . import fvars
from .ml import maximize, options as ml_options, print_log
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .ssm import kalman_loglik
from .tsops import require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

# (level variance, trend variance, irregular, fixed drift, fixed level) of each trend model
MODELS = {
    "none": (False, False, False, False, False), "ntrend": (False, False, True, False, False),
    "dconstant": (False, False, True, False, True), "llevel": (True, False, True, False, False),
    "rwalk": (True, False, False, False, False), "dtrend": (False, False, True, True, True),
    "lldtrend": (True, False, True, True, False), "rwdrift": (True, False, False, True, False),
    "lltrend": (True, True, True, False, False), "strend": (False, True, True, False, False),
    "rtrend": (False, True, False, False, False),
}
TITLE = {"none": "", "ntrend": "", "dconstant": "deterministic constant", "llevel": "local level",
         "rwalk": "random walk", "dtrend": "deterministic trend", "lldtrend": "local level with deterministic trend",
         "rwdrift": "random walk with drift", "lltrend": "local linear trend", "strend": "smooth trend",
         "rtrend": "random trend"}


@dataclass
class UCSpec:
    model: str = "rwalk"
    seasonal: int = 0
    depname: str = ""
    cycles: list[int] = field(default_factory=list)       # orders
    nreg: int = 0


def param_layout(spec: UCSpec, depname: str) -> list[tuple[str, str]]:
    """[(name, transform)] in the order of e(b); transform: id, var, freq, damp."""
    lv, tv, irr, drift, _lvl = MODELS[spec.model]
    out = [(f"reg{i}", "id") for i in range(spec.nreg)]
    nc = len(spec.cycles)
    for i in range(nc):
        sfx = "" if nc == 1 else str(i + 1)
        out += [(f"frequency{sfx}", "freq"), (f"damping{sfx}", "damp")]
    if drift:
        out.append(("beta", "id"))
    if lv:
        out.append(("var(level)", "var"))
    if tv:
        out.append(("var(slope)", "var"))
    if spec.seasonal:
        out.append(("var(seasonal)", "var"))
    for i in range(nc):
        out.append((f"var(cycle{i + 1})" if nc > 1 else "var(cycle1)", "var"))
    if irr:
        out.append((f"var({depname})", "var"))
    return out


def to_nat(u: float, kind: str) -> float:
    if kind == "var":
        return math.exp(min(u, 700))
    if kind == "freq":
        return math.pi / (1 + math.exp(-u))
    if kind == "damp":
        return 1 / (1 + math.exp(-u))
    return u


def to_unc(x: float, kind: str) -> float:
    if kind == "var":
        return math.log(max(x, 1e-12))
    if kind == "freq":
        return math.log(x / (math.pi - x))
    if kind == "damp":
        return math.log(x / (1 - x))
    return x


def build_state_space(spec: UCSpec, nat: dict):
    """(D, A, S, H, c, diffuse) for the natural-scale parameters in `nat`."""
    lv, tv, irr, drift, fixlvl = MODELS[spec.model]
    blocks_A, blocks_D, vars_, diff, cvec = [], [], [], [], []
    shocks = []
    # trend
    if spec.model not in ("none", "ntrend"):
        slope = spec.model in ("lltrend", "strend", "rtrend")
        n = 2 if slope else 1
        A = np.eye(n)
        if slope:
            A[0, 1] = 1.0
        S = np.zeros((n, n))
        if lv:
            S[0, 0] = nat["var(level)"]
        if tv:
            S[n - 1, n - 1] = nat["var(slope)"]
        c = np.zeros(n)
        if drift:
            c[0] = nat["beta"]
        D = np.zeros(n); D[0] = 1.0
        blocks_A.append(A); blocks_D.append(D); shocks.append(S); diff.append(np.ones(n, dtype=bool)); cvec.append(c)
    # seasonal
    if spec.seasonal:
        ns = spec.seasonal - 1
        A = np.zeros((ns, ns))
        A[0, :] = -1.0
        if ns > 1:
            A[1:, :-1] = np.eye(ns - 1)
        S = np.zeros((ns, ns)); S[0, 0] = nat["var(seasonal)"]
        D = np.zeros(ns); D[0] = 1.0
        blocks_A.append(A); blocks_D.append(D); shocks.append(S); diff.append(np.ones(ns, dtype=bool))
        cvec.append(np.zeros(ns))
    # cycles
    nc = len(spec.cycles)
    for i, order in enumerate(spec.cycles):
        sfx = "" if nc == 1 else str(i + 1)
        lam, rho = nat[f"frequency{sfx}"], nat[f"damping{sfx}"]
        var = nat[f"var(cycle{i + 1})"]
        Rm = rho * np.array([[math.cos(lam), math.sin(lam)], [-math.sin(lam), math.cos(lam)]])
        n = 2 * order
        A = np.zeros((n, n))
        for j in range(order):
            A[2 * j:2 * j + 2, 2 * j:2 * j + 2] = Rm
            if j + 1 < order:
                A[2 * j:2 * j + 2, 2 * j + 2:2 * j + 4] = np.eye(2)
        S = np.zeros((n, n)); S[n - 2:, n - 2:] = var * np.eye(2)
        D = np.zeros(n); D[0] = 1.0
        blocks_A.append(A); blocks_D.append(D); shocks.append(S); diff.append(np.zeros(n, dtype=bool))
        cvec.append(np.zeros(n))
    from scipy.linalg import block_diag
    if not blocks_A:
        A = np.zeros((0, 0)); S = np.zeros((0, 0)); D = np.zeros((1, 0)); c = np.zeros(0); dm = np.zeros(0, dtype=bool)
    else:
        A, S = block_diag(*blocks_A), block_diag(*shocks)
        D = np.concatenate(blocks_D).reshape(1, -1)
        c, dm = np.concatenate(cvec), np.concatenate(diff)
    H = np.array([[nat.get(f"var({spec.depname})", 0.0) if irr else 0.0]])
    return D, A, S, H, c, dm


_OPTS = {"model": 3, "seasonal": 4, "cycle": 3, "constraints": 3, "collinear": 4, "vce": 3, "robust": 1,
         "level": 1, "nocnsreport": 5, "nolog": 3, "iterate": 4, "tolerance": 3, "ltolerance": 4,
         "nrtolerance": 4, "technique": 4, "difficult": 4, "from": 3, "trace": 5}


def cmd_ucm(s: "Session", args: str) -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    from .tsmodels import _check_gaps, _numderiv, _sample_text
    if not args.strip() or args.strip().startswith(","):
        replay(s, "ucm", args)
        return
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    # cycle() may be repeated: collect all of them before the option matcher keeps only the last
    import re
    cyc = [int(float(m_)) for m_ in re.findall(r"\bcycle\(\s*(\d+)", opts)]
    o = match_options(re.sub(r"\bcycle\((?:[^()]|\([^()]*\))*\)", "", opts), _OPTS) if opts.strip() else {}
    for bad in ("constraints", "from", "vce"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() not implemented in OpenDTA")   # OpenDTA
    # cycle(#, frequency(#f)): the starting value of that cycle's frequency, in radians
    freq0 = [float(m_.group(1)) if m_ else None for m_ in
             (re.search(r"frequency\(\s*([-+.\deE]+)\s*\)", c_) for c_ in
              re.findall(r"\bcycle\(((?:[^()]|\([^()]*\))*)\)", opts))]
    for f_ in freq0:
        if f_ is not None and not 0 < f_ < math.pi:
            raise StataError(198, "cycle frequency must be in (2*_pi/_N, _pi) and must not be too close to _pi")
    spec = UCSpec()
    spec.model = str(o["model"]).strip().lower() if isinstance(o.get("model"), str) else "rwalk"
    if spec.model not in MODELS:
        raise StataError(198, "option model() incorrectly specified")
    spec.seasonal = int(float(o["seasonal"])) if isinstance(o.get("seasonal"), str) else 0
    if spec.seasonal == 1 or spec.seasonal < 0:
        raise StataError(198, "seasonal() must be at least 2")   # VERIFY
    spec.cycles = cyc
    for q in cyc:
        if q not in (1, 2, 3):
            raise StataError(198, "cycle() order must be 1, 2 or 3")
    if len(spec.cycles) > 3:
        raise StataError(198, "at most three cycles")
    p = parse_standard(body)
    text = p.varlist + (f" if {p.if_}" if p.if_ else "") + (f" in {p.in_}" if p.in_ else "")
    smp = Sample(s, text + ", noconstant" + (f" level({o['level']})" if isinstance(o.get("level"), str) else ""),
                 opts_spec={}, weights=())
    t0, t1 = _check_gaps(s, smp.mask)
    fvars.drop_collinear(smp.cols, False)
    est_cols = [col for col in smp.cols if col.values is not None]
    n = smp.n
    info = require_ts(s.data)
    order = np.argsort(s.data.get(info.tvar).data[smp.mask], kind="stable")
    X = np.column_stack([col.values for col in est_cols])[order] if est_cols else np.zeros((n, 0))
    y = smp.y.astype(float)[order]
    spec.nreg = len(est_cols)
    spec.depname = smp.depname
    layout = param_layout(spec, smp.depname)
    names = [nm for nm, _k in layout]
    kinds = [k_ for _nm, k_ in layout]
    nreg = spec.nreg

    def natural(u):
        return {nm: to_nat(v, kd) for nm, v, kd in zip(names, u, kinds)}

    def contrib(u):
        nat = natural(u)
        D, A, S, H, cc, dm = build_state_space(spec, nat)
        beta = np.array([nat[f"reg{i}"] for i in range(nreg)])
        dd = (X @ beta).reshape(-1, 1) if nreg else None
        return kalman_loglik(y.reshape(-1, 1), D, A, S, H, cc, dd, dm)

    def ll(u):
        if any(abs(v) > 25 for v, k in zip(u, kinds) if k == "var"):
            return -np.inf            # variances beyond exp(+-25) are numerically meaningless
        c_ = contrib(u)
        # VERIFY (case 1904): the reference log is ln(s) below this value for seasonal(s) (s = 4 seen)
        return -np.inf if c_ is None else float(c_.sum()) - (math.log(spec.seasonal) if spec.seasonal else 0.0)

    def fun(u, todo):
        if todo == 0:
            return ll(u), None, None
        return _numderiv(ll, u)

    vy = float(np.var(np.diff(y))) or 1.0
    start = {"var(level)": 0.3 * vy, "var(slope)": 0.05 * vy, "var(seasonal)": 0.05 * vy, "beta": float(np.mean(np.diff(y))),
             f"var({smp.depname})": 0.5 * vy}
    u0 = []
    beta0 = np.linalg.lstsq(X, y, rcond=None)[0] if nreg else np.zeros(0)
    ncyc = 0
    for nm, kd in layout:
        if nm.startswith("reg"):
            u0.append(float(beta0[int(nm[3:])]))
        elif nm.startswith("frequency"):
            f_ = freq0[ncyc] if ncyc < len(freq0) else None
            ncyc += 1
            u0.append(to_unc(f_ if f_ is not None else 2 * math.pi / max(8.0, n / 6), kd))
        elif nm.startswith("damping"):
            u0.append(to_unc(0.7, kd))
        elif nm.startswith("var(cycle"):
            u0.append(to_unc(0.3 * vy, kd))
        else:
            u0.append(to_unc(start[nm], kd) if kd == "var" else start[nm])
    u0 = np.array(u0)
    if not np.isfinite(ll(u0)):
        raise StataError(430, "initial values not feasible")
    r = maximize(fun, u0, **ml_options(o))
    out = s.output
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, r, "log likelihood")
    u = r.b
    _l, _g, Hm = _numderiv(ll, u)
    try:
        Vu = np.linalg.inv(-Hm)
    except np.linalg.LinAlgError:
        Vu = np.linalg.pinv(-Hm)
    # natural scale by the delta method
    nat_vals = np.array([to_nat(v, kd) for v, kd in zip(u, kinds)])
    J = np.zeros(len(u))
    for i, (v, kd) in enumerate(zip(u, kinds)):
        h = 1e-6 * (abs(v) + 1)
        J[i] = (to_nat(v + h, kd) - to_nat(v - h, kd)) / (2 * h)
    V = (J[:, None] * Vu * J[None, :])
    V = (V + V.T) / 2
    se = np.sqrt(np.maximum(np.diag(V), 0))
    D, A, S, H, cc, dm = build_state_space(spec, natural(u))
    stationary = not dm.any()
    # the Wald test covers the regressors and the frequency and damping of the cycles (compat 1904
    # and 1915)
    idx = list(range(nreg)) + [i for i, nm in enumerate(names) if nm.startswith(("frequency", "damping"))]
    chi2 = SYS
    if idx:
        try:
            chi2 = float(nat_vals[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], nat_vals[idx]))
        except np.linalg.LinAlgError:
            chi2 = SYS
    pv = float(chi2_tail(chi2, len(idx))) if idx and chi2 < SYS else SYS
    blocks: list[list[CoefRow]] = []
    if nreg:
        bmap = {col.name: nat_vals[i] for i, col in enumerate(est_cols)}
        semap = {col.name: se[i] for i, col in enumerate(est_cols)}
        blocks.append(rows_from_columns(est_cols, bmap, semap, s.data, constant=False))
    enames, eqs = [], []
    others: list[CoefRow] = []
    vrows: list[CoefRow] = []
    for i, (nm, kd) in enumerate(layout):
        if nm.startswith("reg"):
            enames.append(est_cols[int(nm[3:])].name); eqs.append(smp.depname)
            continue
        (vrows if kd == "var" else others).append(CoefRow("coef", nm, nat_vals[i], se[i], name=nm, tight=kd == "var", zero_lo=kd == "var"))
        enames.append(nm); eqs.append("")
    blocks += [b_ for b_ in (others, vrows) if b_]
    rows: list[CoefRow] = []
    for b_ in blocks:
        if rows:
            rows.append(CoefRow("sep"))
        rows += b_
    est = Estimates("ucm", smp.depname, enames, nat_vals, V, n, "z", None)
    est.eqnames = [e_ or "_" for e_ in eqs]
    est.rows = rows
    est.level = float(o["level"]) if isinstance(o.get("level"), str) else smp.level
    est.vcetype = "OIM"
    comps = [TITLE[spec.model]] if TITLE[spec.model] else []
    if spec.seasonal:
        comps.append(f"seasonal({spec.seasonal})")
    if spec.cycles:
        comps.append("order " + " ".join(map(str, spec.cycles)) + (" cycle" if len(spec.cycles) == 1 else " cycles")
                     if len(spec.cycles) == 1 else f"{len(spec.cycles)} cycles of order " + " ".join(map(str, spec.cycles)))
    est.extra = {"N": float(n), "ll": r.ll, "chi2": chi2, "p": pv, "df_m": len(idx), "sample": _sample_text(s, t0, t1),
                 "components": ", ".join(comps), "stationary": stationary, "nreg": nreg}
    est.ucm = {"spec": spec, "layout": layout, "u": u, "y": y, "X": X, "order": order, "mask": smp.mask.copy()}
    est.display = display_ucm
    est.predict = _predict_ucm
    scalars = [("N", float(n)), ("k", float(len(u))), ("k_aux", float(len(u) - nreg)), ("k_eq", float(1 + (nreg > 0))),
               ("k_dv", 1.0), ("k_cycles", float(len(spec.cycles))), ("df_m", float(len(idx))), ("ll", r.ll),
               ("chi2", chi2), ("p", pv), ("tmin", t0), ("tmax", t1), ("stationary", 1.0 if stationary else 0.0),
               ("rank", float(len(u))), ("ic", float(r.iterations)), ("rc", 0.0), ("converged", 1.0 if r.converged else 0.0)]
    macros = [("cmd", "ucm"), ("cmdline", cmdline("ucm", args)), ("depvar", smp.depname),
              ("model", spec.model), ("title", "Unobserved-components model"), ("tvar", info.tvar),
              ("vce", "oim"), ("vcetype", "OIM"), ("chi2type", "Wald"), ("properties", "b V")]
    post(s, est, scalars, macros, smp.mask)
    display_ucm(s, est)


def display_ucm(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, **_kw) -> None:
    from .tsmodels import _g, _right
    out = s.output
    x = est.extra
    if header:
        out.write("\nUnobserved-components model\n", "text")
        out.write(f"Components: {x['components']}".rstrip() + "\n\n", "text")
        _right(out, x["sample"].replace("Sample:  ", "Sample: "), "Number of obs", f"{int(x['N']):,}")
        if x["df_m"]:
            _right(out, "", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f"))
        ll = f"Log likelihood = {_g(x['ll'], '%10.0g'):>10}"
        if x["df_m"]:
            _right(out, ll, "Prob > chi2", _g(x["p"], "%10.4f"))
        else:
            out.write(ll + "\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype)
        if not x["stationary"]:
            out.write("Note: Model is not stationary.\n", "text")
        out.write("Note: Tests of variances against zero are one sided, and the two-sided\n"
                  "      confidence intervals are truncated at zero.\n", "text")


def _predict_ucm(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    """One-step predictions (xb, the default), residuals and standardized residuals over the
    estimation sample (VERIFY: trend, seasonal, cycle and the other smethod() options)."""
    from .ssm import kalman_onestep
    U = est.ucm
    spec, layout, u = U["spec"], U["layout"], U["u"]
    nat = {nm: to_nat(v, kd) for (nm, kd), v in zip(layout, u)}
    D, A, S, H, cc, dm = build_state_space(spec, nat)
    beta = np.array([nat[f"reg{i}"] for i in range(spec.nreg)])
    dd = (U["X"] @ beta).reshape(-1, 1) if spec.nreg else None
    yhat, fv = kalman_onestep(U["y"].reshape(-1, 1), D, A, S, H, cc, dd, dm)
    rows = np.nonzero(U["mask"])[0][U["order"]]
    # the first observations, one per diffuse state, have no prediction (compat 1915: 299 of 300
    # with a random walk; VERIFY with other models)
    nd = int(np.asarray(dm).sum())
    yhat = yhat.copy()
    yhat[:nd] = np.nan
    out = np.full(s.data.nobs, SYS)
    if stat in ("", "xb"):
        out[rows] = np.where(np.isnan(yhat[:, 0]), SYS, yhat[:, 0])
        return out, "" if stat else "(option xb assumed; fitted values)"
    e = U["y"] - yhat[:, 0]
    if stat == "residuals":
        out[rows] = np.where(np.isnan(e), SYS, e)
    elif stat == "rstandard":
        z = e / np.sqrt(fv[:, 0])
        out[rows] = np.where(np.isnan(z), SYS, z)
    else:
        raise StataError(198, f"option {stat} not allowed")                      # VERIFY
    return out, ""


def cycle_density(est: Estimates, which: int, omega: np.ndarray, normalized: bool = True) -> np.ndarray:
    """Spectral density of the cycle number `which` (1, 2, 3) of a ucm fit, from the state-space
    form of that cycle alone: (1/2pi) D (I - A z)^-1 S (I - A z)^-H D', z = exp(-i w); divided by
    the cycle's variance unless pspectrum (VERIFY: psdensity's normalization after ucm)."""
    from .ssm import stationary_cov
    U = est.ucm
    spec = U["spec"]
    if not spec.cycles:
        raise StataError(321, "no cycle in the model")                               # VERIFY
    if not 1 <= which <= len(spec.cycles):
        raise StataError(198, f"cycle({which}) not found")                          # VERIFY
    nat = {nm: to_nat(v, kd) for (nm, kd), v in zip(U["layout"], U["u"])}
    nc = len(spec.cycles)
    sfx = "" if nc == 1 else str(which)
    one = UCSpec(model="none", depname=spec.depname, cycles=[spec.cycles[which - 1]])
    sub = {"frequency": nat[f"frequency{sfx}"], "damping": nat[f"damping{sfx}"],
           "var(cycle1)": nat["var(cycle1)" if nc == 1 else f"var(cycle{which})"]}
    D, A, S, _H, _c, _d = build_state_space(one, sub)
    r = A.shape[0]
    f = np.empty(len(omega))
    for i, w in enumerate(omega):
        G = np.linalg.inv(np.eye(r) - A * np.exp(-1j * w))
        f[i] = float(np.real(D @ G @ S @ G.conj().T @ D.T)[0, 0]) / (2 * math.pi)
    if normalized:
        P = stationary_cov(A, S)
        f = f / float((D @ P @ D.T)[0, 0])
    return f
