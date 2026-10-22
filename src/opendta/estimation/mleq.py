"""Engine for maximum-likelihood models with several linear indexes (equations) and ancillary
parameters: zip, hetprobit, heckman, biprobit... ([R] ml).

A model gives `lik(etas, theta) -> (ll_obs, [d ll / d eta_j], d ll / d theta (n x a) or None)`
(see likelihoods.py). `Engine` builds the weighted log likelihood, its score and Hessian over the
stacked parameters [beta_1, ..., beta_m, theta]; `Block` is one equation (its columns, which
include factor variables, and its constant); `assemble` lays out e(b), e(V) and the table rows
the way Stata shows several equations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from ..core import missing as M
from . import fvars
from .models2 import _expand, _header, g, num_hess
from .results import CoefRow, crit, coef_table, rows_from_columns

SYS = M.SYSMISS


@dataclass
class Block:
    """One equation: columns of its varlist (free ones form X), constant and offset."""
    name: str
    cols: list
    terms: list
    constant: bool
    X: np.ndarray
    offset: np.ndarray | None = None
    notes: list = field(default_factory=list)
    off_name: str = ""

    @property
    def k(self) -> int:
        return self.X.shape[1]

    @property
    def nslopes(self) -> int:
        return self.k - (1 if self.constant else 0)


def block_from_pr(pr, name: str) -> Block:
    """The main equation of a `_Prep`: its columns plus the constant."""
    n = pr.smp.n
    X = np.column_stack([pr.X, np.ones(n)]) if pr.constant else pr.X
    return Block(name, pr.smp.cols, pr.smp.terms, pr.constant, X)


def block_from_text(pr, name: str, text: str, constant: bool = True) -> Block:
    """An extra equation from a varlist (factor variables allowed), on the sample of `pr`."""
    ds = pr.s.data
    mask = pr.smp.mask
    with fvars.caching():
        terms = fvars.expand_fv(ds, text) if text.strip() else []
        cols = fvars.build_design(ds, terms, mask, constant=constant)
    notes = fvars.drop_collinear(cols, constant, pr.w)
    free = [c.values for c in cols if c.values is not None]
    n = pr.smp.n
    parts = free + ([np.ones(n)] if constant else [])
    X = np.column_stack(parts) if parts else np.zeros((n, 0))
    return Block(name, cols, terms, constant, X, notes=notes)


def constant_only(blk: Block) -> Block:
    n = blk.X.shape[0]
    return Block(blk.name, [], [], True, np.ones((n, 1)), blk.offset)


class Engine:
    """Weighted log likelihood over [beta_1 ... beta_m, theta]."""

    def __init__(self, lik: Callable, blocks: list[Block], w: np.ndarray, n_theta: int = 0):
        self.lik = lik
        self.blocks = blocks
        self.w = w
        self.n_theta = n_theta
        self.ks = [b.k for b in blocks]
        self.offs = np.concatenate([[0], np.cumsum(self.ks)]).astype(int)
        self.K = int(self.offs[-1]) + n_theta

    def split(self, b: np.ndarray):
        etas = []
        for j, blk in enumerate(self.blocks):
            e = blk.X @ b[self.offs[j]:self.offs[j + 1]]
            if blk.offset is not None:
                e = e + blk.offset
            etas.append(e)
        return etas, b[self.offs[-1]:]

    def ll_obs(self, b: np.ndarray) -> np.ndarray:
        etas, th = self.split(b)
        return self.lik(etas, th)[0]

    def scores(self, b: np.ndarray) -> np.ndarray:
        etas, th = self.split(b)
        _, d_eta, d_th = self.lik(etas, th)
        cols = [blk.X * d[:, None] for blk, d in zip(self.blocks, d_eta)]
        if self.n_theta:
            cols.append(np.asarray(d_th).reshape(len(d_eta[0]), self.n_theta))
        return np.column_stack(cols)

    def grad(self, b: np.ndarray) -> np.ndarray:
        return self.scores(b).T @ self.w

    def fun(self, b: np.ndarray, todo: int):
        with np.errstate(all="ignore"):
            ll = float(self.w @ self.ll_obs(b))
            if not np.isfinite(ll):
                return -np.inf, None, None
            if todo == 0:
                return ll, None, None
            gr = self.grad(b)
            if not np.all(np.isfinite(gr)):
                return -np.inf, None, None
            return ll, gr, num_hess(self.grad, b)


def add_base_row(est, label: str, text: str, before: str | None = None) -> None:
    """The row `label |  text` of the base alternative/outcome, between two equation blocks: ahead
    of the block named `before`, or after the last coefficient block when there is none."""
    old = est.rows_fn

    def with_base(level: float):
        rows = list(old(level))
        i = len(rows)
        for k, r in enumerate(rows):
            if before is not None and r.kind == "eq" and r.label == before:
                i = k - 1
                break
        rows[i:i] = [CoefRow("sep"), CoefRow("eq", label, eq=text)]
        return rows

    est.rows_fn = with_base
    est.rows = with_base(est.level)


def aux_row(label: str, value: float, se: float, level: float, *, ci=None, name: str = "_cons"):
    """Row of an ancillary parameter (/athrho, rho...). An estimated parameter (`ci` None) has z and
    P>|z|; one derived by a transformation (`ci` = (lo, hi), from the estimate's own interval)
    has neither."""
    if ci is None:
        return CoefRow("coef", label, value, se, name=name, noeform=True)
    r = CoefRow("aux", label, value, se, name=name)
    r.eq = f"{float(ci[0])!r} {float(ci[1])!r}"
    return r


def transformed_row(label: str, f: Callable, df_: Callable, b: float, se: float, level: float):
    """Row for f(b) with delta-method standard error and the interval of b transformed."""
    c = crit("z", level, None)
    return aux_row(label, float(f(b)), float(abs(df_(b)) * se), level,
                   ci=(float(f(b - c * se)), float(f(b + c * se))))


def assemble(pr, blocks: list[Block], b: np.ndarray, V: np.ndarray, aux_names: list[tuple[str, str]],
             *, heads: bool | None = None):
    """Names, equation names, e(b), e(V) and the table rows of the coefficient blocks. `aux_names`
    = [(equation, name)] for the ancillary parameters at the end of b (in e(b) they come after
    the blocks). `heads`: write an equation header before each block (default: when there is more
    than one). The rows of the ancillary parameters are left to the caller."""
    ds = pr.s.data
    names, eqs, layout, rows = [], [], [], []
    start = 0
    pos_all: list[int] = []
    heads = len(blocks) > 1 if heads is None else heads
    for j, blk in enumerate(blocks):
        bmap, semap = {}, {}
        slot = start
        for c in blk.cols:
            nm = c.name if not c.omitted else fvars.omitted_name(c.name)
            names.append(nm)
            eqs.append(blk.name)
            if c.values is not None:
                layout.append((nm, True))
                pos_all.append(slot)
                bmap[c.name] = float(b[slot])
                semap[c.name] = float(np.sqrt(max(V[slot, slot], 0)))
                slot += 1
            else:
                layout.append((nm, False))
                bmap[c.name], semap[c.name] = 0.0, 0.0
        if blk.constant:
            names.append("_cons")
            eqs.append(blk.name)
            layout.append(("_cons", True))
            pos_all.append(slot)
            bmap["_cons"] = float(b[slot])
            semap["_cons"] = float(np.sqrt(max(V[slot, slot], 0)))
            slot += 1
        if heads:
            if j:
                rows.append(CoefRow("sep"))
            rows.append(CoefRow("eq", blk.name))
        rows.extend(rows_from_columns(blk.cols, bmap, semap, ds, constant=blk.constant))
        start += blk.k
    for eq, nm in aux_names:
        names.append(nm)
        eqs.append(eq)
        layout.append((nm, True))
        pos_all.append(start)
        start += 1
    kk = len(names)
    b_full = np.zeros(kk)
    V_full = np.zeros((kk, kk))
    sel = [i for i, (_, free) in enumerate(layout) if free]
    b_full[sel] = b[pos_all]
    V_full[np.ix_(sel, sel)] = V[np.ix_(pos_all, pos_all)]
    return names, eqs, b_full, V_full, rows


def display_eq(s, est, *, header: bool = True, table: bool = True, level: float | None = None, **_kw) -> None:
    """Header (left lines, right pairs), table and footer of a model made with `Engine`."""
    lev = level if level is not None else est.level
    if header:
        hd = est.hdr_fn(est)
        if callable(hd):
            hd(s)
        else:
            _header(s, *hd)
    if table:
        if est.vcetype and est.clustvar and getattr(est, "group_cluster", False):
            s.output.write(f"{'(Std. Err. adjusted for clustering on ' + est.clustvar + ')':>78}\n", "text")
        elif est.vcetype and est.clustvar:
            s.output.write(f"{'(Std. Err. adjusted for ' + format(int(est.n_clust), ',') + ' clusters in ' + est.clustvar + ')':>78}\n",
                           "text")
        rows = est.rows_fn(lev) if getattr(est, "rows_fn", None) else est.rows
        coef_table(s, est.depvar_cell, rows, stat="z", df=None, level=lev, vcetype=est.vcetype,
                   coef_title=getattr(est, "coef_title", "Coef."), eform=getattr(est, "eform_on", False))
        foot = getattr(est, "footer", None)
        if foot:
            foot(s, est)


def fit_log(s, pr, title: str, res, label: str | None = None) -> None:
    """The log of a preliminary fit under its heading ("Fitting full model:"), as the ado-file
    commands print them."""
    if pr.o.get("nolog"):
        return
    from .models2 import _log
    s.output.write(f"\n{title}\n", "text")
    _log(s, pr, res, label)


def lab(est) -> str:
    return "Log pseudolikelihood" if est.robust else "Log likelihood"


def std_right(x: dict, extra: list[tuple[str, str]] = ()) -> list[tuple[str, str]]:
    """Number of obs, [extra], chi2, Prob > chi2, Pseudo R2 (when there is one)."""
    from .regress import comma
    right = [("Number of obs", comma(x["N"]))] + list(extra)
    right += [(f"{x['chi2type']} chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
              ("Prob > chi2", g(x["p"], "%10.4f"))]
    if x.get("r2_p") is not None:
        right.append(("Pseudo R2", g(x["r2_p"], "%10.4f")))
    return right


def finish_model(s, pr, cmd: str, args: str, *, title: str, blocks: list[Block], eng: Engine, res, V,
                 n_clust, ll0: float, aux_names: list[tuple[str, str]], aux_rows_fn: Callable,
                 hdr_fn: Callable, depvar_cell: str | None = None, footer: Callable | None = None,
                 wald: bool = False, extra_scalars=(), extra_macros=(), mats=(), order=None,
                 test_idx: list[int] | None = None, df_m: int | None = None, k_eq: int | None = None,
                 heads: bool | None = None, noeform_after_first: bool = True, r2p: float | None = None,
                 predict_opts: dict | None = None):
    """Lays out e(b)/e(V), tests the slopes of the first equation, posts the results and returns
    the Estimates (display is `display_eq`). `aux_rows_fn(level)` gives the rows after the
    coefficient blocks (ancillary parameters and their transformations)."""
    from .models2 import _finish, _model_test
    names, eqs, b_full, V_full, rows = assemble(pr, blocks, res.b, V, aux_names, heads=heads)
    idx = list(range(blocks[0].nslopes)) if test_idx is None else test_idx
    if wald:
        saved = pr.robust
        pr.robust = True
        kind, chi2, pv = _model_test(pr, res.b, V, idx, res.ll, ll0)
        pr.robust = saved
    else:
        kind, chi2, pv = _model_test(pr, res.b, V, idx, res.ll, ll0)
    if r2p is None and ll0 < SYS and ll0 and not wald and not pr.robust:
        r2p = 1 - res.ll / ll0
    if noeform_after_first and len(blocks) > 1:
        seen_first = False
        for r in rows:
            if r.kind == "sep" and not seen_first:
                seen_first = True
            if seen_first:
                r.noeform = True

    def rows_fn(level: float):
        return rows + list(aux_rows_fn(level))

    est = _finish(s, pr, cmd, args, title=title, names=names, eqnames=eqs, b=b_full, V=V_full,
                  rows=rows_fn(pr.smp.level), res=res, ll0=ll0, df_m=len(idx) if df_m is None else df_m,
                  chi2kind=kind, chi2=chi2, p=pv, r2p=r2p, display=display_eq, n_clust=n_clust,
                  k_eq=len(blocks) + len({e for e, _ in aux_names}) if k_eq is None else k_eq,
                  k_aux=len(aux_names), extra_scalars=extra_scalars, extra_macros=extra_macros,
                  constant=blocks[0].constant, mats=mats, order=order, k_free=len(res.b))
    est.rows_fn = rows_fn
    est.hdr_fn = hdr_fn
    est.depvar_cell = pr.smp.depname if depvar_cell is None else depvar_cell
    est.footer = footer
    est.blocks = blocks
    est.predict_opts = predict_opts or {}
    est.n_coef = sum(len(b.cols) + (1 if b.constant else 0) for b in blocks)
    return est


# ---------------------------------------------------------------------------
# options and postestimation helpers
# ---------------------------------------------------------------------------

_COMMON_OPTS = {"noconstant": 4, "level": 1, "vce": 3, "robust": 1, "cluster": 2, "nolog": 4, "log": 3,
                "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "difficult": 4,
                "technique": 4, "noheader": 6}


def peek(args: str, spec: dict[str, int]) -> dict:
    """Options of a command line before the sample is built (to find the variables that
    option arguments name)."""
    from ..lang.syntax import match_options, parse_standard
    p = parse_standard(args)
    full = dict(spec)
    full.update(_COMMON_OPTS)
    return match_options(p.options, full) if p.options.strip() else {}


def sub_options(text: str, spec: dict[str, int]) -> tuple[str, dict]:
    """'z1 z2, offset(v)' -> ('z1 z2', {'offset': 'v'})."""
    from ..lang.syntax import match_options
    text = text.strip()
    if "," not in text:
        return text, {}
    head, tail = text.split(",", 1)
    return head.strip(), match_options(tail, spec)


def names_in(*texts: str) -> str:
    """Variable names to add to the sample from option arguments (blank-separated)."""
    return " ".join(t for t in texts if t and t.strip())


def offset_values(pr, off_opt, exp_opt):
    """(offset array on the sample, label for e()/the table) from offset() or exposure()."""
    ds = pr.s.data
    if exp_opt:
        nm = str(exp_opt).strip()
        v = ds.get(nm).data[pr.smp.mask].astype(float)
        if np.any(v <= 0):
            from ..core.errors import StataError
            raise StataError(459, f"{nm} <= 0 in sample")   # VERIFY
        return np.log(v), f"ln({nm})"
    if off_opt:
        nm = str(off_opt).strip()
        return ds.get(nm).data[pr.smp.mask].astype(float), nm
    return None, ""


def _block_columns(s, blk):
    """Yields (position in the block's free coefficients, column for all observations)."""
    ds = s.data
    n = ds.nobs
    for c in blk.cols:
        if c.values is None:
            continue
        t = blk.terms[c.term]
        col = np.ones(n)
        kk = 0
        for comp in t.comps:
            x = fvars.comp_values(ds, comp)
            if comp.factor:
                col = col * (x == c.cell[kk])
                kk += 1
            else:
                col = col * np.where(x < SYS, x, 0.0)
        yield c, col


def block_index(s, est, j: int):
    """Linear index of equation j for all observations (offset included) and the mask of
    observations where its variables are not missing."""
    ds = s.data
    blk = est.blocks[j]
    pos = sum(len(b.cols) + (1 if b.constant else 0) for b in est.blocks[:j])
    eta = np.zeros(ds.nobs)
    with fvars.caching():
        ok = ~fvars.term_missing(ds, blk.terms) if blk.terms else np.ones(ds.nobs, bool)
        idx = {id(c): i for i, c in enumerate(blk.cols)}
        for c, col in _block_columns(s, blk):
            bj = est.b[pos + idx[id(c)]]
            if bj != 0:
                eta += bj * col
        if blk.constant:
            eta += est.b[pos + len(blk.cols)]
    if getattr(blk, "off_name", ""):
        nm = blk.off_name
        if nm.startswith("ln("):
            v = ds.get(nm[3:-1]).data
            eta += np.where(v < SYS, np.log(np.where(v > 0, v, 1)), 0.0)
        else:
            v = ds.get(nm).data
            eta += np.where(v < SYS, v, 0.0)
    return eta, ok


def make_predict(fn: Callable, default: str = ""):
    """predict for a model made with Engine. fn(stat, etas, theta, est) returns the values for
    `stat` or None when the statistic is not available; stat "" is replaced by `default` with
    Stata's note."""
    from ..core.errors import StataError

    def predict(s, est, stat, mask):
        idxs = [block_index(s, est, j) for j in range(len(est.blocks))]
        ok = np.all([o for _, o in idxs], axis=0)
        etas = [e for e, _ in idxs]
        th = est.b[est.n_coef:] if hasattr(est, "n_coef") else np.zeros(0)
        what = stat or default
        if what == "xb":
            vals, note = etas[0], ""
        else:
            vals = fn(what, etas, th, est)
            note = f"(option {what} assumed; predicted {est.pred_label})" if not stat and getattr(est, "pred_label", "") else ""
            if vals is None:
                raise StataError(198, f"option {stat} not allowed")
        return np.where(ok, vals, SYS), note
    return predict
