"""xtfrontier: stochastic frontier models for panel data ([XT] xtfrontier).

The likelihood is in core/panelfrontier.py. Models: ti (time-invariant inefficiency, the default)
and tvd (time-varying decay). The layout of the table follows the manual's description; details
not settled by it carry a VERIFY mark and the case is in compat 1905."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import panelfrontier as pf
from ..core.errors import StataError
from ..core.panelquad import Groups
from .ml import maximize, print_log
from .panel2 import _finish, _Panel
from .results import CoefRow, crit
from .tsops import ts_info
from .xtre import _numeric_fun, _vce

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def cmd_xtfrontier(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtfrontier", args)
        return
    P = _Panel(s, "xtfrontier", args, extra_opts={"ti": 2, "tvd": 3, "cost": 2, "t": 1}, constant=None)
    o, pr, out = P.pr.o, P.pr, s.output
    if o.get("ti") and o.get("tvd"):
        raise StataError(198, "only one of ti or tvd may be specified")           # VERIFY
    if o.get("fe") or o.get("re") or o.get("pa"):
        raise StataError(198, "option fe not allowed")                             # VERIFY
    tvd = bool(o.get("tvd"))
    cost = bool(o.get("cost"))
    ds = s.data
    info = ts_info(ds)
    tvar = str(o["t"]).strip() if isinstance(o.get("t"), str) else (info.tvar if info is not None else "")
    if not tvar:
        raise StataError(459, "time variable not set; use xtset panelvar timevar ...")
    delta = 1.0 if isinstance(o.get("t"), str) or info is None else info.delta
    t = ds.get(tvar).data[pr.smp.mask].astype(float)
    tidx = np.rint((t - t.min()) / (delta or 1.0)).astype(int)
    pr.write_notes(out)
    y = P.y
    X = np.column_stack([pr.X, np.ones(len(y))]) if pr.constant else pr.X
    k = X.shape[1]
    grp = Groups(P.inv)
    d1 = np.ones(len(y))

    def unpack(par):
        b = par[:k]
        mu = par[k]
        j = k + 1
        eta = 0.0
        if tvd:
            eta = par[j]
            j += 1
        l2, lg = par[j], par[j + 1]
        s2 = math.exp(min(max(l2, -700), 700))
        g = 1 / (1 + math.exp(-lg))
        return b, mu, eta, s2, g

    def llg(par):
        b, mu, eta, s2, g = unpack(par)
        su2, sv2 = g * s2, (1 - g) * s2
        if su2 <= 0 or sv2 <= 0:
            return np.full(grp.G, -np.inf)
        d = pf.decay(tidx, grp, eta) if tvd else d1
        return pf.group_loglik(y - X @ b, d, grp, mu, su2, sv2, cost)

    def ll(par):
        v = float(llg(par).sum())
        return v if np.isfinite(v) else -np.inf
    b0 = np.linalg.lstsq(X, y, rcond=None)[0]
    s2_0 = float(np.mean((y - X @ b0) ** 2)) or 1.0
    if pr.constant:
        b0[-1] += (-1 if cost else 1) * math.sqrt(2 / math.pi) * math.sqrt(0.5 * s2_0)
    start = np.concatenate([b0, [0.0], [0.0] if tvd else [], [math.log(s2_0), 0.0]])
    res = pr.maximize(_numeric_fun(ll), start)
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, res, "log likelihood")
    V, robust = _vce(P, llg, res, grp)
    b, mu, eta, s2, g = unpack(res.b)
    n = len(res.b)
    j_mu = k
    j_eta = k + 1 if tvd else None
    j_l2 = n - 2
    j_lg = n - 1
    se = np.sqrt(np.maximum(np.diag(V), 0))
    z = crit("z", pr.smp.level, None)

    def delta_se(grad):
        gr = np.zeros(n)
        for idx, v in grad.items():
            gr[idx] = v
        return float(math.sqrt(max(gr @ V @ gr, 0)))
    su2, sv2 = g * s2, (1 - g) * s2
    dg = g * (1 - g)
    aux = [CoefRow("sep"), CoefRow("aux", "/mu", mu, se[j_mu], name="mu")]
    names_extra = ["/mu"]
    if tvd:
        aux.append(CoefRow("aux", "/eta", eta, se[j_eta], name="eta"))
        names_extra.append("/eta")
    aux += [CoefRow("aux", "/lnsigma2", res.b[j_l2], se[j_l2], name="lnsigma2"),
            CoefRow("aux", "/ilgtgamma", res.b[j_lg], se[j_lg], name="ilgtgamma"), CoefRow("sep")]
    names_extra += ["/lnsigma2", "/ilgtgamma"]
    se_s2 = delta_se({j_l2: s2})
    se_g = delta_se({j_lg: dg})
    se_su = delta_se({j_l2: su2, j_lg: s2 * dg})
    se_sv = delta_se({j_l2: sv2, j_lg: -s2 * dg})
    for nm, val, sev in (("sigma2", s2, se_s2), ("gamma", g, se_g), ("sigma_u2", su2, se_su), ("sigma_v2", sv2, se_sv)):
        aux.append(CoefRow("aux", nm, val, sev, name=nm))
    title = "Time-varying decay inefficiency model" if tvd else "Time-invariant inefficiency model"
    sub = f"Time variable: {tvar}"
    extra_sc = [("sigma2", s2), ("gamma", g), ("sigma_u2", su2), ("sigma_v2", sv2), ("mu", mu)]
    if tvd:
        extra_sc.append(("eta", eta))
    _finish(s, P, "xtfrontier", args, res.b, V, res, names_extra, title, sub, np.ones(len(y), bool),
            aux_rows=aux, robust=robust, extra_scalars=extra_sc,
            extra_macros=[("function", "cost" if cost else "production"), ("tvar", tvar),
                          ("model", "tvd" if tvd else "ti")])
