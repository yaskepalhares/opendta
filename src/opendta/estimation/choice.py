"""Choice models: clogit, rologit, asclogit, slogit, mprobit, asmprobit and asroprobit
([R] clogit, [R] rologit, [R] asclogit, [R] slogit, [R] mprobit, [R] asmprobit, [R] asroprobit).
The likelihoods are in core/choice.py; the layout of the equations in mleq.py.

Layouts come from the manuals; the lines marked VERIFY wait for the Stata 14 log of case 1908."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import choice as C
from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..core.varlist import resolve_name
from . import fvars
from .mleq import (Block, Engine, aux_row, display_eq, finish_model, lab, names_in, peek, std_right)
from .models2 import _Prep, g
from .regress import comma

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _replay(s, cmd, args) -> bool:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return True
    return False


def _group_values(pr, name: str) -> np.ndarray:
    ds = pr.s.data
    var = ds.get(name)
    vals = var.data[pr.smp.mask]
    if var.is_string:
        codes = {v: i for i, v in enumerate(sorted(set(vals)))}
        return np.array([codes[v] for v in vals], dtype=float)
    return vals.astype(float)


def _restrict(pr, keep: np.ndarray) -> None:
    """Drops observations from the sample of a prepared command (groups that carry no information)."""
    if keep.all():
        return
    smp = pr.smp
    full = np.flatnonzero(smp.mask)
    smp.mask = smp.mask.copy()
    smp.mask[full[~keep]] = False
    smp.n = int(smp.mask.sum())
    pr.w = pr.w[keep]
    pr.y = pr.y[keep]
    pr.X = pr.X[keep]
    pr.N = float(pr.w.sum()) if smp.wtype in ("fweight", "iweight") else float(smp.n)
    if getattr(smp, "cluster", None) is not None and len(np.atleast_1d(smp.cluster)) == len(keep):
        smp.cluster = smp.cluster[keep]
    for c in smp.cols:
        if c.values is not None:
            c.values = c.values[keep]


def _notes_dropped(s, name: str, ngroups: int, nobs: int, why: str) -> None:
    if ngroups:
        s.output.write(f"note: {ngroups} group{'s' if ngroups != 1 else ''} ({nobs} obs) dropped because of "
                       f"{why}\n", "text")


def _group_setup(s, args: str, spec: dict, cmd: str):
    o = peek(args, spec)
    if o.get("group") in (None, True):
        raise StataError(100, "option group() required")
    gname = resolve_name(s.data, str(o["group"]).strip())
    pr = _Prep(s, args, spec, constant=False, extra_vars=names_in(gname, o.get("offset")))
    pr.write_notes(s.output)
    gid = _group_values(pr, gname)
    return pr, o, gname, gid


def _cluster_by_group(pr, gname: str, gid: np.ndarray) -> None:
    """vce(robust) in these commands clusters on the group variable."""
    smp = pr.smp
    if pr.robust and smp.vce == "robust":
        smp.vce, smp.clustvar = "cluster", gname
        smp.cluster = gid.copy()


def _block(pr, name: str, X: np.ndarray) -> Block:
    return Block(name, pr.smp.cols, pr.smp.terms, False, X)


# ---------------------------------------------------------------------------
# clogit
# ---------------------------------------------------------------------------

def cmd_clogit(s: "Session", args: str) -> None:
    if _replay(s, "clogit", args):
        return
    spec = {"group": 1, "or": 2, "offset": 3, "strata": 4}
    pr, o, gname, gid = _group_setup(s, args, spec, "clogit")
    y = (pr.y != 0).astype(float)
    groups = C.group_index(gid)
    # groups with only positive or only negative outcomes carry no information
    drop_g = [idx for idx in groups if y[idx].min() == y[idx].max()]
    if drop_g:
        keep = np.ones(len(y), dtype=bool)
        for idx in drop_g:
            keep[idx] = False
        _notes_dropped(s, "clogit", len(drop_g), int(sum(len(i) for i in drop_g)),
                       "all positive or all negative outcomes.")
        if not keep.any():
            raise StataError(2000, "no observations")
        _restrict(pr, keep)
        y, gid = y[keep], gid[keep]
        groups = C.group_index(gid)
    if any(y[idx].sum() > 1 for idx in groups):
        s.output.write("note: multiple positive outcomes within groups encountered.\n", "text")
    _cluster_by_group(pr, gname, gid)
    main = _block(pr, pr.smp.depname, pr.X)
    w = pr.w

    def lik(etas, th):
        ll, d = C.clogit(etas[0], y, groups)
        return ll, [d], None

    eng = Engine(lik, [main], w, 0)
    res = pr.maximize(eng.fun, np.zeros(main.k))
    from .mleq import fit_log
    V, n_clust = pr.vce(res.H, _group_scores(eng, res.b, groups) if pr.robust else None)
    ll0 = C.clogit_null(y, groups)
    ngroups = len(groups)
    nalt = np.array([len(i) for i in groups])

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [est.title, "", "", f"{lab(est)} = {ll4:>9}"]
        right = std_right(x)
        if est.robust:
            right = [r for r in right if r[0] != "Pseudo R2"] + [("Pseudo R2", g(1 - x["ll"] / ll0, "%10.4f"))]
        return left, right

    est = finish_model(s, pr, "clogit", args, title="Conditional (fixed-effects) logistic regression",
                       blocks=[main], eng=eng, res=res, V=V, n_clust=n_clust, ll0=ll0, aux_names=[],
                       aux_rows_fn=lambda lev: [], hdr_fn=hdr, k_eq=1, r2p=1 - res.ll / ll0 if ll0 else None,
                       extra_scalars=[("N_g", float(ngroups)), ("g_min", float(nalt.min())),
                                      ("g_avg", float(nalt.mean())), ("g_max", float(nalt.max())),
                                      ("k_out", float(sum(1 for idx in groups if y[idx].sum() > 0)))],
                       extra_macros=[("group", gname)])
    est.coef_title, est.eform_on = ("Odds Ratio", True) if o.get("or") else ("Coef.", False)
    display_eq(s, est)
    est.coef_title, est.eform_on = "Coef.", False


def _group_scores(eng, b, groups) -> np.ndarray:
    """Scores of the groups (one row each), the ones the sandwich needs when the groups are the
    clusters; the observation scores summed within the group."""
    S = eng.scores(b)
    return S


# ---------------------------------------------------------------------------
# rologit
# ---------------------------------------------------------------------------

def cmd_rologit(s: "Session", args: str) -> None:
    if _replay(s, "rologit", args):
        return
    spec = {"group": 1, "ties": 2, "reverse": 3, "or": 2, "offset": 3, "incomplete": 3, "notestrg": 6}
    pr, o, gname, gid = _group_setup(s, args, spec, "rologit")
    ties = str(o.get("ties", "exactp")).strip().lower()
    for full in ("none", "exactp", "breslow", "efron"):
        if full.startswith(ties) and ties:
            ties = full
            break
    else:
        raise StataError(198, "option ties() invalid")   # VERIFY
    rank = pr.y.astype(float)
    if o.get("reverse"):
        rank = -rank
    groups = C.group_index(gid)
    drop_g = [idx for idx in groups if len(idx) < 2 or rank[idx].min() == rank[idx].max()]
    if drop_g:
        keep = np.ones(len(rank), dtype=bool)
        for idx in drop_g:
            keep[idx] = False
        _notes_dropped(s, "rologit", len(drop_g), int(sum(len(i) for i in drop_g)),
                       "all alternatives tied or a single alternative")   # VERIFY
        if not keep.any():
            raise StataError(2000, "no observations")
        _restrict(pr, keep)
        rank, gid = rank[keep], gid[keep]
        groups = C.group_index(gid)
    if ties == "none" and any(len(np.unique(rank[i])) < len(i) for i in groups):
        raise StataError(459, "ties found in the rankings; specify ties()")   # VERIFY
    _cluster_by_group(pr, gname, gid)
    main = _block(pr, pr.smp.depname, pr.X)
    w = pr.w

    def lik(etas, th):
        ll, d = C.rologit(etas[0], rank, groups, ties)
        return ll, [d], None

    eng = Engine(lik, [main], w, 0)
    res = pr.maximize(eng.fun, np.zeros(main.k))
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ll0 = float(sum(C.rologit_group(np.zeros(len(i)), rank[i], ties)[0] for i in groups))
    nalt = np.array([len(i) for i in groups])
    ngroups = len(groups)

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [est.title, f"Group variable: {gname}", "", f"Ties handled via the {ties} method", "", f"{lab(est)} = {ll4:>9}"]
        right = [("Number of obs", comma(x["N"])), ("Number of groups", comma(ngroups)),
                 ("Alts per group: min", comma(nalt.min())), ("avg", g(nalt.mean(), "%9.1f")),
                 ("max", comma(nalt.max()))] + std_right(x)[1:]   # VERIFY layout
        return left, right

    est = finish_model(s, pr, "rologit", args, title="Rank-ordered logistic regression", blocks=[main],
                       eng=eng, res=res, V=V, n_clust=n_clust, ll0=ll0, aux_names=[], aux_rows_fn=lambda lev: [],
                       hdr_fn=hdr, k_eq=1, r2p=None, extra_scalars=[("N_g", float(ngroups)),
                       ("g_min", float(nalt.min())), ("g_avg", float(nalt.mean())), ("g_max", float(nalt.max()))],
                       extra_macros=[("group", gname), ("ties", ties)])
    est.coef_title, est.eform_on = ("Odds Ratio", True) if o.get("or") else ("Coef.", False)
    display_eq(s, est)
    est.coef_title, est.eform_on = "Coef.", False


# ---------------------------------------------------------------------------
# asclogit
# ---------------------------------------------------------------------------

def _alt_names(ds, altvar: str, levels: np.ndarray) -> list[str]:
    var = ds.get(altvar)
    labels = ds.value_labels.get(var.value_label, {}) if var.value_label else {}
    out = []
    for v in levels:
        if var.is_string:
            out.append(str(v))
        else:
            iv = int(v)
            out.append(labels.get(iv, str(iv)))
    return out


def _alt_setup(s, args: str, spec: dict, cmd: str):
    o = peek(args, spec)
    for key in ("case", "alternatives"):
        if o.get(key) in (None, True):
            raise StataError(100, f"option {key}() required")
    ds = s.data
    case = resolve_name(ds, str(o["case"]).strip())
    alt = resolve_name(ds, str(o["alternatives"]).strip())
    casevars = str(o.get("casevars", "")).strip()
    pr = _Prep(s, args, spec, constant=False, extra_vars=names_in(case, alt, casevars))
    pr.write_notes(s.output)
    cid = _group_values(pr, case)
    mask = pr.smp.mask
    avar = ds.get(alt)
    araw = avar.data[mask]
    if avar.is_string:
        levels = np.array(sorted(set(araw)))
        acode = np.array([int(np.searchsorted(levels, v)) for v in araw])
        lev_names = [str(v) for v in levels]
    else:
        levels = np.unique(araw.astype(float))
        acode = np.searchsorted(levels, araw.astype(float))
        lev_names = _alt_names(ds, alt, levels)
    return pr, o, case, alt, casevars, cid, acode, levels, lev_names


class _Alt:
    """Design of the commands on case/alternative data (long form)."""


def _asc_design(s, args: str, spec: dict, cmd: str, *, ranking: bool = False) -> _Alt:
    from .mleq import block_from_text
    pr, o, case, alt, casevars, cid, acode, levels, lev_names = _alt_setup(s, args, spec, cmd)
    d = _Alt()
    d.pr, d.o, d.case, d.alt, d.cid, d.acode, d.levels, d.lev_names = pr, o, case, alt, cid, acode, levels, lev_names
    y = (pr.y != 0).astype(float) if not ranking else pr.y.astype(float)
    cases = C.group_index(cid)
    if ranking:
        bad = [idx for idx in cases if len(np.unique(y[idx])) < len(idx)]
        why = "ties in the rankings"
    else:
        bad = [idx for idx in cases if y[idx].sum() != 1]
        why = ""
    if bad:
        keep = np.ones(len(y), dtype=bool)
        for idx in bad:
            keep[idx] = False
        if ranking:
            s.output.write(f"note: {len(bad)} case{'s' if len(bad) != 1 else ''} ({int(sum(len(i) for i in bad))} obs) "
                           "dropped because of ties or missing rankings\n", "text")   # VERIFY text
        else:
            n_no = sum(1 for idx in bad if y[idx].sum() == 0)
            if n_no:
                s.output.write(f"note: {n_no} case{'s' if n_no != 1 else ''} "
                               f"({int(sum(len(i) for i in bad if y[i].sum() == 0))} obs) "
                               "dropped because of no positive outcome\n", "text")   # VERIFY text
            n_mu = len(bad) - n_no
            if n_mu:
                s.output.write(f"note: {n_mu} case{'s' if n_mu != 1 else ''} "
                               f"({int(sum(len(i) for i in bad if y[i].sum() > 1))} obs) "
                               "dropped because of multiple positive outcomes\n", "text")   # VERIFY text
        if not keep.any():
            raise StataError(2000, "no observations")
        _restrict(pr, keep)
        y, cid, acode = y[keep], cid[keep], acode[keep]
        d.cid, d.acode = cid, acode
        cases = C.group_index(cid)
    d.y, d.cases = y, cases
    counts = np.bincount(acode, minlength=len(levels))
    base_opt = o.get("basealternative")
    if base_opt in (None, True):
        base = int(np.argmax(counts))
    else:
        b_txt = str(base_opt).strip()
        names_l = [n.lower() for n in lev_names]
        if b_txt.lower() in names_l:
            base = names_l.index(b_txt.lower())
        else:
            try:
                base = int(np.flatnonzero(levels.astype(float) == float(b_txt))[0])
            except (ValueError, IndexError):
                raise StataError(198, f"alternative {b_txt} not found")   # VERIFY
    d.base = base
    _cluster_by_group(pr, case, cid)
    smp = pr.smp
    d.n = smp.n
    d.main = Block(alt, smp.cols, smp.terms, False, pr.X)
    d.blocks = [d.main]
    casevars = str(o.get("casevars", "")).strip()
    for j in range(len(levels)):
        if j == base:
            continue
        ind = (acode == j).astype(float)
        blk = block_from_text(pr, lev_names[j], casevars, constant=not o.get("noconstant"))
        blk.X = blk.X * ind[:, None]
        d.blocks.append(blk)
    d.ncase = len(cases)
    d.nalt = np.array([len(i) for i in cases])
    return d


def _asc_header(d: _Alt, title: str):
    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        llab = getattr(est, "ll_label", None) or lab(est)
        left = [est.title, f"Case variable: {d.case}", "", f"Alternative variable: {d.alt}", "", "", "", f"{llab} = {ll4:>9}"]
        right = [("Number of obs", comma(x["N"])), ("Number of cases", comma(d.ncase)), ("", ""),
                 ("Alts per case: min", comma(d.nalt.min())), ("avg", g(d.nalt.mean(), "%9.1f")),
                 ("max", comma(d.nalt.max())), ("", ""),
                 (f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")), ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left, right
    return hdr


def _asc_finish(s, args, cmd, title, d: _Alt, eng, res, V, n_clust, ll0, *, aux_names=(), aux_rows_fn=None,
                extra_scalars=(), extra_macros=()):
    idx = []
    for j, blk in enumerate(d.blocks):
        idx += list(range(eng.offs[j], eng.offs[j] + blk.nslopes))

    def footer(s_, est):
        s_.output.write(f"(base alternative is {d.lev_names[d.base]})\n", "text")

    est = finish_model(s, d.pr, cmd, args, title=title, blocks=d.blocks, eng=eng, res=res, V=V, n_clust=n_clust,
                       ll0=ll0, aux_names=list(aux_names), aux_rows_fn=aux_rows_fn or (lambda lev: []),
                       hdr_fn=_asc_header(d, title), wald=True, test_idx=idx, df_m=len(idx),
                       k_eq=len(d.blocks), footer=footer, heads=True, noeform_after_first=False,
                       extra_scalars=[("N_case", float(d.ncase)), ("k_alt", float(len(d.levels))),
                                      ("alt_min", float(d.nalt.min())), ("alt_avg", float(d.nalt.mean())),
                                      ("alt_max", float(d.nalt.max())), ("i_base", float(d.base + 1))]
                       + list(extra_scalars),
                       extra_macros=[("case", d.case), ("altvar", d.alt), ("base", d.lev_names[d.base])]
                       + list(extra_macros))
    return est


def cmd_asclogit(s: "Session", args: str) -> None:
    if _replay(s, "asclogit", args):
        return
    spec = {"case": 4, "alternatives": 3, "casevars": 5, "basealternative": 3, "altwise": 3, "or": 2,
            "offset": 3, "noconstant": 4}
    d = _asc_design(s, args, spec, "asclogit")
    pr, o = d.pr, d.o

    def lik(etas, th):
        eta = np.zeros(d.n)
        for e in etas:
            eta = eta + e
        ll, dd = C.clogit(eta, d.y, d.cases)
        return ll, [dd] * len(etas), None

    eng = Engine(lik, d.blocks, pr.w, 0)
    res = pr.maximize(eng.fun, np.zeros(eng.K))
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    est = _asc_finish(s, args, "asclogit", "Alternative-specific conditional logit", d, eng, res, V, n_clust,
                      C.clogit_null(d.y, d.cases))
    est.coef_title, est.eform_on = ("Odds Ratio", True) if o.get("or") else ("Coef.", False)
    display_eq(s, est)
    est.coef_title, est.eform_on = "Coef.", False


# ---------------------------------------------------------------------------
# slogit
# ---------------------------------------------------------------------------

def cmd_slogit(s: "Session", args: str) -> None:
    """Stereotype logistic regression with one dimension: ln{P(y=k)/P(y=b)} = theta_k + phi_k x'beta,
    phi_b = 0, theta_b = 0 and phi = 1 for the lowest non-base outcome."""
    if _replay(s, "slogit", args):
        return
    spec = {"dimension": 3, "baseoutcome": 3, "or": 2, "nocoefficient": 6}
    from scipy import optimize as opt
    pr = _Prep(s, args, spec, constant=False)
    pr.write_notes(s.output)
    o = pr.o
    if o.get("dimension") and str(o["dimension"]).strip() != "1":
        raise StataError(198, "only dimension(1) is available")   # VERIFY
    y = pr.y
    levels = np.unique(y)
    K = len(levels)
    if K < 3:
        raise StataError(148, "outcome does not take at least three values")   # VERIFY
    code = np.searchsorted(levels, y)
    if o.get("baseoutcome"):
        bv = float(o["baseoutcome"])
        if bv not in levels:
            raise StataError(459, f"baseoutcome({o['baseoutcome']}) is not an outcome")   # VERIFY
        base = int(np.flatnonzero(levels == bv)[0])
    else:
        base = int(np.argmax(np.bincount(code)))
    first = 0 if base != 0 else 1                         # phi = 1
    free_phi = [k for k in range(K) if k != base and k != first]
    free_th = [k for k in range(K) if k != base]
    nth = len(free_phi) + len(free_th)
    main = _block(pr, pr.smp.depname, pr.X)
    w = pr.w
    n = pr.smp.n
    Y = np.zeros((n, K))
    Y[np.arange(n), code] = 1.0

    def parts(th):
        phi = np.zeros(K)
        phi[first] = 1.0
        for i, k in enumerate(free_phi):
            phi[k] = th[i]
        theta = np.zeros(K)
        for i, k in enumerate(free_th):
            theta[k] = th[len(free_phi) + i]
        return phi, theta

    def lik(etas, th):
        eta = etas[0]
        phi, theta = parts(np.asarray(th))
        S = theta[None, :] + eta[:, None] * phi[None, :]
        lse = np.logaddexp.reduce(S, axis=1)
        P = np.exp(S - lse[:, None])
        ll = (Y * S).sum(axis=1) - lse
        d_eta = (Y - P) @ phi
        dth = np.zeros((n, nth))
        for i, k in enumerate(free_phi):
            dth[:, i] = eta * (Y[:, k] - P[:, k])
        for i, k in enumerate(free_th):
            dth[:, len(free_phi) + i] = Y[:, k] - P[:, k]
        return ll, [d_eta], dth

    eng = Engine(lik, [main], w, nth)

    def neg(b):
        v = eng.fun(b, 1)
        return -v[0], -v[1]

    best = None
    rng = np.random.default_rng(1)
    Xs = pr.X / (pr.X.std(axis=0) + 1e-12) if pr.X.size else pr.X
    score = np.linalg.lstsq(np.column_stack([pr.X, np.ones(n)]), code.astype(float), rcond=None)[0][:-1]
    for trial in range(4):
        th0 = np.zeros(nth)
        for i, k in enumerate(free_phi):
            th0[i] = 1.0 - (k - first) / max(K - 1 - first, 1) if k > first else 1.0
        b0 = np.concatenate([score * (1 if trial == 0 else rng.normal(0.5, 0.5)), th0]) if trial < 2 else \
            np.concatenate([rng.normal(0, 0.5, main.k), th0 + rng.normal(0, 0.3, nth)])
        try:
            r = opt.minimize(neg, b0, jac=True, method="BFGS", options={"gtol": 1e-8, "maxiter": 2000})
        except Exception:
            continue
        if best is None or r.fun < best.fun:
            best = r
    if best is None:
        raise StataError(430, "convergence not achieved")
    res = pr.maximize(eng.fun, best.x)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    nk = np.bincount(code, minlength=K).astype(float)
    ll0 = float(np.sum(nk * np.log(nk / nk.sum())))
    off = main.k
    from .results import CoefRow
    phi_all, th_all = parts(res.b[off:])
    se = np.sqrt(np.maximum(np.diag(V), 0))
    lev_txt = [g(float(v), "%9.0g") for v in levels]

    def aux_rows(level):
        rows = [CoefRow("sep")]
        for k in range(K):
            nm = f"/phi1_{k + 1}"
            if k in free_phi:
                i = free_phi.index(k)
                rows.append(aux_row(nm, float(res.b[off + i]), float(se[off + i]), level))
            else:
                r = CoefRow("value", nm, float(phi_all[k]))
                r.eq = f"{phi_all[k]:>10g}  (constrained)"
                rows.append(r)
        rows.append(CoefRow("sep"))
        for k in range(K):
            nm = f"/theta{k + 1}"
            if k in free_th:
                i = len(free_phi) + free_th.index(k)
                rows.append(aux_row(nm, float(res.b[off + i]), float(se[off + i]), level))
            else:
                r = CoefRow("value", nm, 0.0)
                r.eq = f"{0:>10g}  (base outcome)"
                rows.append(r)
        return rows

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [est.title, "", "", f"{lab(est)} = {ll4:>9}"]
        right = std_right(x)
        return left, right

    names = [("phi1_%d" % (k + 1), "_cons") for k in free_phi] + [("theta%d" % (k + 1), "_cons") for k in free_th]
    est = finish_model(s, pr, "slogit", args, title="Stereotype logistic regression", blocks=[main], eng=eng,
                       res=res, V=V, n_clust=n_clust, ll0=ll0, aux_names=names, aux_rows_fn=aux_rows, hdr_fn=hdr,
                       k_eq=1 + len(names), test_idx=list(range(main.k)), df_m=main.k,
                       extra_scalars=[("k_out", float(K)), ("i_base", float(base + 1)), ("dim", 1.0)],
                       extra_macros=[("baseoutcome", lev_txt[base])])
    est.coef_title, est.eform_on = "Coef.", False
    display_eq(s, est)
