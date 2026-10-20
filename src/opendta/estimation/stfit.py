"""What stcox and streg share: the estimation sample of st data, the header of their output
and the scalars that describe the survival data ([ST] stcox, [ST] streg)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import find_top, parse_standard
from . import fvars
from .results import Estimates
from .st import n_subjects, require_st, st_data

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def fmt(x: float, fmt_: str = "%10.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt_, pad=False).strip()


class STSample:
    """Regressors and sample of stcox/streg: the observations with _st == 1, the weights given to
    stset, _t, _t0 and _d of the sample, and the design matrix after collinearity checks."""

    def __init__(self, s: "Session", args: str, spec: dict, extra_vars: str = ""):
        from .regress import Sample
        info = require_st(s)
        ds = s.data
        self.info = info
        if parse_standard(args).weight is not None:
            raise StataError(101, "weights not allowed; specify weights with stset")   # VERIFY
        # the regressors come after a dummy depvar (_t) so that Sample builds the design matrix
        c = find_top(args, ",")
        body, opts = (args, None) if c == -1 else (args[:c], args[c + 1:])
        wt = f" [{info.wtype}={info.wvar}]" if info.wtype else ""
        text = f"_t {body.strip()}{wt}" + (f", {opts}" if opts is not None else "")
        extra = info.idvar if info.idvar and not ds.get(info.idvar).is_string else ""
        extra = f"{extra} {extra_vars}".strip()
        smp = Sample(s, text, opts_spec=spec, weights=("fweight", "pweight", "iweight"),
                     extra_vars=extra, restrict=ds.get("_st").data == 1)
        self.smp = smp
        self.o = smp.o
        self.mask = smp.mask
        self.sd = sd = st_data(s, smp.mask)
        self.t, self.t0, self.d = sd.t, sd.t0, sd.d
        self.w = smp.w[smp.mask].astype(np.float64) if smp.w is not None else np.ones(smp.n)
        self.wtype = smp.wtype or ""
        self.fweight = self.wtype == "fweight"
        self.n = smp.n
        self.notes = fvars.drop_collinear(smp.cols, False, self.w)
        self.est_cols = [c for c in smp.cols if c.values is not None]
        self.X = np.column_stack([c.values for c in self.est_cols]) if self.est_cols else np.zeros((smp.n, 0))
        self.robust = smp.vce in ("robust", "cluster") or self.wtype == "pweight"

    def scalars(self) -> dict:
        """N, N_sub, N_fail and risk (frequency weights count each observation w times)."""
        fw = self.w if self.fweight else np.ones(self.n)
        return {"N_sub": n_subjects(self.sd), "N_fail": float(np.sum(self.d * fw)),
                "risk": float(np.sum((self.t - self.t0) * fw)),
                "N": float(fw.sum()) if self.fweight else float(self.n)}

    def clusters(self) -> tuple[np.ndarray | None, str]:
        """Clusters of the robust variance: vce(cluster) or, with robust and stset id(), the
        subjects."""
        smp = self.smp
        if smp.vce == "cluster":
            return smp.cluster, smp.clustvar
        if smp.vce == "robust" and self.sd.ids is not None:
            return self.sd.ids, self.info.idvar
        return None, ""

    def notes_and_header(self, s: "Session") -> None:
        from .st import header
        if not self.o.get("noshow"):
            header(s, self.info)
        for nm in self.notes:
            s.output.write(f"note: {nm} omitted because of collinearity\n", "text")


def header_block(s: "Session", title: str, est: Estimates) -> None:
    """No. of subjects, failures and time at risk at the left; N, the model test and its p-value
    at the right."""
    x = est.extra
    out = s.output
    out.write(f"\n{title}\n\n", "text")
    llab = "Log pseudolikelihood" if est.robust else "Log likelihood"
    lw = len(llab) if est.robust else 15
    chi = f"{x['chi2type']} chi2({int(x['df_m'])})"
    left = [("No. of subjects", fmt(x["N_sub"], "%12.0gc")), ("No. of failures", fmt(x["N_fail"], "%12.0gc")),
            ("Time at risk", fmt(x["risk"], "%12.0g")), None, (llab, fmt(x["ll"], "%10.0g"))]
    right = [("Number of obs", fmt(x["N"], "%12.0gc")), None, None,
             (chi, f"{x['chi2']:.2f}" if x["chi2"] < SYS else "."),
             ("Prob > chi2", f"{x['p']:.4f}" if x["p"] < SYS else ".")]
    for lft, rgt in zip(left, right):
        line = f"{lft[0]:<{lw}} ={lft[1]:>13}" if lft is not None else ""
        if rgt is not None:
            out.write(f"{line:<48}{rgt[0]:<17}=", "text")
            out.write(f"{rgt[1]:>12}\n", "result")
        else:
            out.write(line + "\n", "text")
    out.write("\n", "text")


def cluster_note(s: "Session", est: Estimates) -> None:
    if est.clustvar:
        msg = f"(Std. Err. adjusted for {int(est.n_clust):,} clusters in {est.clustvar})"
        s.output.write(f"{msg:>78}\n", "text")


def place(smp, b: np.ndarray, V: np.ndarray, *, constant: bool, aux: int = 0):
    """Coefficients and variance of the estimated columns placed in the full e(b) layout (omitted
    columns at zero), then _cons and `aux` ancillary parameters. Returns (names, b, V, positions,
    bmap, semap) with bmap/semap keyed by column name, as rows_from_columns expects."""
    names = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols]
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None]
    if constant:
        pos.append(len(names))
        names.append("_cons")
    pos += [len(names) + j for j in range(aux)]
    names += ["_cons"] * aux
    k = len(names)
    bf, Vf = np.zeros(k), np.zeros((k, k))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = V
    se = np.sqrt(np.maximum(np.diag(Vf), 0))
    bmap = {c.name: bf[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(se[i]) for i, c in enumerate(smp.cols)}
    if constant:
        j = len(smp.cols)
        bmap["_cons"], semap["_cons"] = bf[j], float(se[j])
    return names, bf, Vf, pos, bmap, semap
