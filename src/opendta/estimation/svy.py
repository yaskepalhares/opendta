"""Complex surveys ([SVY] svyset, svydescribe, svy estimation) and the
statistics mean, proportion, total and ratio ([R] mean, [R] proportion...).

Variance by Taylor linearization (design with strata, PSUs and fpc
in one stage):

    V = Σ_h (1 - f_h) n_h/(n_h - 1) Σ_j (z_hj - z̄_h)(z_hj - z̄_h)'

z_hj is the sum, in PSU j of stratum h, of the linearized scores of each
observation. Design degrees of freedom = number of PSUs - number of strata.
The declaration is kept in the _dta characteristics (_svy_*), as in Stata.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from .results import g9

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# svyset
# ---------------------------------------------------------------------------

@dataclass
class Design:
    psu: str = ""            # "" or "_n" = each observation is a PSU
    weight: str = ""
    wtype: str = ""          # pweight / iweight
    strata: str = ""
    fpc: str = ""
    singleunit: str = "missing"
    vce: str = "linearized"
    poststrata: str = ""
    postweight: str = ""
    # stages 2, 3...: [{"su": ..., "strata": ..., "fpc": ...}]
    stages: list = field(default_factory=list)
    # replicate weights and options of the replication methods
    brrweight: str = ""
    fay: float = 0.0
    bsrweight: str = ""
    bsn: float = 1.0
    jkrweight: str = ""
    sdrweight: str = ""
    mse: bool = False
    dof: float | None = None

    @property
    def is_set(self) -> bool:
        return True

    @property
    def repweights(self) -> list[str]:
        txt = {"brr": self.brrweight, "bootstrap": self.bsrweight, "jackknife": self.jkrweight,
               "sdr": self.sdrweight}.get(self.vce, "")
        return txt.split() if txt else []


def get_design(ds) -> Design | None:
    ch = ds.chars.get("_dta", {})
    if ch.get("_svy_version", "") == "":
        return None
    stages = []
    k = 2
    while f"_svy_su{k}" in ch:
        stages.append({"su": ch.get(f"_svy_su{k}", ""), "strata": ch.get(f"_svy_strata{k}", ""),
                       "fpc": ch.get(f"_svy_fpc{k}", "")})
        k += 1
    dof = ch.get("_svy_dof", "")
    return Design(psu=ch.get("_svy_su1", ""), weight=ch.get("_svy_wvar", ""), wtype=ch.get("_svy_wtype", ""),
                  strata=ch.get("_svy_strata1", ""), fpc=ch.get("_svy_fpc1", ""),
                  singleunit=ch.get("_svy_singleunit", "missing") or "missing",
                  vce=ch.get("_svy_vce", "linearized") or "linearized",
                  poststrata=ch.get("_svy_poststrata", ""), postweight=ch.get("_svy_postweight", ""),
                  stages=stages, brrweight=ch.get("_svy_brrweight", ""), fay=float(ch.get("_svy_fay", "0") or 0),
                  bsrweight=ch.get("_svy_bsrweight", ""), bsn=float(ch.get("_svy_bsn", "1") or 1),
                  jkrweight=ch.get("_svy_jkrweight", ""), sdrweight=ch.get("_svy_sdrweight", ""),
                  mse=ch.get("_svy_mse", "") == "mse", dof=float(dof) if dof else None)


def _show_design(s: "Session", d: Design) -> None:
    out = s.output
    rows = []
    if d.weight:
        rows.append((d.wtype, d.weight))
    else:
        rows.append(("pweight", "<none>"))
    rows += [("VCE", d.vce)]
    if d.vce != "linearized":
        # replication methods show MSE: off/on (compat 0602)
        rows.append(("MSE", "on" if d.mse else "off"))
    for lab, val in (("brrweight", d.brrweight), ("bsrweight", d.bsrweight), ("jkrweight", d.jkrweight),
                     ("sdrweight", d.sdrweight)):
        if val:
            rows.append((lab, val))
    if d.vce == "brr" and d.fay:
        rows.append(("Fay", f"{d.fay:g}"))
    if d.vce == "bootstrap" and d.bsn != 1:
        rows.append(("bsn", f"{d.bsn:g}"))
    rows += [("Single unit", d.singleunit), ("Strata 1", d.strata or "<one>"),
             ("SU 1", d.psu or "<observations>"), ("FPC 1", d.fpc or "<zero>")]
    for k, st_ in enumerate(d.stages, start=2):
        rows += [(f"Strata {k}", st_["strata"] or "<one>"), (f"SU {k}", st_["su"] or "<observations>"),
                 (f"FPC {k}", st_["fpc"] or "<zero>")]
    if d.poststrata:
        rows[1:1] = [("Poststrata", d.poststrata), ("Postweight", d.postweight)]
    out.write("\n", "text")
    for lab, val in rows:
        out.write(f"{lab:>13}: ", "text")
        out.write(val + "\n", "result")


_SVYSET_OPTS = {"strata": 3, "fpc": 3, "singleunit": 6, "vce": 3, "poststrata": 5, "postweight": 5,
                "clear": 5, "brrweight": 3, "fay": 3, "bsrweight": 3, "bsn": 3, "jkrweight": 3,
                "sdrweight": 3, "mse": 3, "dof": 3, "noclear": 7}


def _split_stages(text: str) -> list[str]:
    parts, depth, cur, i = [], 0, "", 0
    while i < len(text):
        c = text[i]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
        if depth == 0 and text.startswith("||", i):
            parts.append(cur)
            cur = ""
            i += 2
            continue
        cur += c
        i += 1
    parts.append(cur)
    return parts


def cmd_svyset(s: "Session", args: str) -> None:
    ds = s.data
    t = args.strip()
    if not t:
        d = get_design(ds)
        if d is None:
            s.output.write("no survey characteristics are set\n", "text")   # VERIFY
            return
        _show_design(s, d)
        return
    if t.startswith(",") and re.search(r"\bclear\b", t):
        for k in list(ds.chars.get("_dta", {})):
            if k.startswith("_svy_"):
                del ds.chars["_dta"][k]
        return
    from ..core.varlist import resolve_name
    pieces = _split_stages(t)
    # the general options (replicate weights, vce...) may come after the last stage
    last = pieces[-1]
    p = parse_standard(pieces[0])
    o = match_options(p.options, _SVYSET_OPTS) if p.options.strip() else {}
    stages = []
    for piece in pieces[1:]:
        ps = parse_standard(piece)
        so = match_options(ps.options, _SVYSET_OPTS) if ps.options.strip() else {}
        su = ps.varlist.strip()
        stage = {"su": "" if su in ("_n", "") else resolve_name(ds, su),
                 "strata": resolve_name(ds, str(so["strata"]).strip()) if so.get("strata") else "",
                 "fpc": resolve_name(ds, str(so["fpc"]).strip()) if so.get("fpc") else ""}
        stages.append(stage)
        for k_ in ("vce", "brrweight", "fay", "bsrweight", "bsn", "jkrweight", "sdrweight", "mse", "dof",
                   "singleunit"):
            if k_ in so:
                o[k_] = so[k_]
    del last
    words = p.varlist.split()
    psu = ""
    if words:
        psu = "_n" if words[0] == "_n" else resolve_name(ds, words[0])
    d = Design(psu=psu if psu != "_n" else "", stages=stages)
    if p.weight:
        wtype, wexp = p.weight
        if wtype not in ("pweight", "iweight", "weight"):
            raise StataError(101, f"{wtype}s not allowed")   # VERIFY
        d.wtype = "pweight" if wtype == "weight" else wtype
        d.weight = resolve_name(ds, wexp.strip())
    if o.get("strata"):
        d.strata = resolve_name(ds, str(o["strata"]).strip())
    if o.get("fpc"):
        d.fpc = resolve_name(ds, str(o["fpc"]).strip())
    if o.get("singleunit"):
        su = str(o["singleunit"]).strip()
        if su not in ("missing", "certainty", "scaled", "centered"):
            raise StataError(198, "invalid singleunit() option")   # VERIFY
        d.singleunit = su
    for key in ("brrweight", "bsrweight", "jkrweight", "sdrweight"):
        if o.get(key):
            main = str(o[key]).split(",")[0]
            setattr(d, key, " ".join(s.expand_varlist(main)))
    if o.get("fay"):
        d.fay = float(o["fay"])
        if not 0 <= d.fay < 2 or d.fay == 1:
            raise StataError(198, "fay() must be in [0, 2) and not 1")   # VERIFY
    if o.get("bsn"):
        d.bsn = float(o["bsn"])
    d.mse = bool(o.get("mse"))
    if o.get("dof"):
        d.dof = float(o["dof"])
    if o.get("vce"):
        v = str(o["vce"]).strip().lower()
        v = {"linear": "linearized", "jack": "jackknife", "jackknife": "jackknife", "bs": "bootstrap",
             "boot": "bootstrap"}.get(v, v)
        d.vce = v
    elif d.brrweight:
        d.vce = "brr"
    elif d.bsrweight:
        d.vce = "bootstrap"
    elif d.jkrweight:
        d.vce = "jackknife"
    elif d.sdrweight:
        d.vce = "sdr"
    if d.vce not in ("linearized", "jackknife", "brr", "bootstrap", "sdr"):
        raise StataError(198, f"vce({d.vce}) not allowed")   # VERIFY
    if d.vce == "bootstrap" and not d.bsrweight:
        raise StataError(198, "bsrweight() required for vce(bootstrap)")   # VERIFY
    if d.vce == "sdr" and not d.sdrweight:
        raise StataError(198, "sdrweight() required for vce(sdr)")   # VERIFY
    if o.get("poststrata"):
        d.poststrata = resolve_name(ds, str(o["poststrata"]).strip())
        d.postweight = resolve_name(ds, str(o.get("postweight", "")).strip())
    ch = ds.chars.setdefault("_dta", {})
    for k in list(ch):
        if k.startswith("_svy_"):
            del ch[k]
    ch.update({"_svy_version": "2", "_svy_su1": d.psu, "_svy_wvar": d.weight, "_svy_wtype": d.wtype,
               "_svy_strata1": d.strata, "_svy_fpc1": d.fpc, "_svy_singleunit": d.singleunit,
               "_svy_vce": d.vce, "_svy_poststrata": d.poststrata, "_svy_postweight": d.postweight,
               "_svy_brrweight": d.brrweight, "_svy_fay": f"{d.fay:g}", "_svy_bsrweight": d.bsrweight,
               "_svy_bsn": f"{d.bsn:g}", "_svy_jkrweight": d.jkrweight, "_svy_sdrweight": d.sdrweight,
               "_svy_mse": "mse" if d.mse else "", "_svy_dof": f"{d.dof:g}" if d.dof else ""})
    for k, st_ in enumerate(stages, start=2):
        ch[f"_svy_su{k}"] = st_["su"]
        ch[f"_svy_strata{k}"] = st_["strata"]
        ch[f"_svy_fpc{k}"] = st_["fpc"]
    _show_design(s, d)
    r = [("settings", t), ("singleunit", d.singleunit), ("vce", d.vce)]
    if d.strata:
        r.append(("strata1", d.strata))
    r.append(("su1", d.psu or "_n"))
    for k, st_ in enumerate(stages, start=2):
        if st_["strata"]:
            r.append((f"strata{k}", st_["strata"]))
        r.append((f"su{k}", st_["su"] or "_n"))
    if d.weight:
        r += [("wvar", d.weight), ("wtype", d.wtype)]
    s.r = dict(r)
    s.r = dict(reversed(list(s.r.items())))


def require_design(s: "Session") -> Design:
    d = get_design(s.data)
    if d is None:
        raise StataError(119, "data not set up for svy, use svyset")
    return d


# ---------------------------------------------------------------------------
# structure of the design in the sample
# ---------------------------------------------------------------------------

def _hadamard(r: int) -> np.ndarray:
    """Hadamard matrix of order r (power of 2, Sylvester construction).
    VERIFY: Stata also uses orders 12, 20, 24...; here, the smallest
    power of 2 that covers the strata."""
    H = np.array([[1.0]])
    while H.shape[0] < r:
        H = np.block([[H, H], [H, -H]])
    return H


class DesignData:
    """Weights, strata and PSUs of the observations in `mask` (design sample)."""

    def __init__(self, s: "Session", d: Design, mask: np.ndarray, vce: str | None = None, mse: bool | None = None,
                 hadamard=False):
        ds = s.data
        self.d = d
        self.s = s
        n_all = ds.nobs
        w = ds.get(d.weight).data.astype(float) if d.weight else np.ones(n_all)
        mask = mask & (w < SYS)
        strata = ds.get(d.strata).data if d.strata else np.zeros(n_all)
        psu = ds.get(d.psu).data if d.psu else np.arange(n_all, dtype=float)
        mask = mask & (strata < SYS) & (psu < SYS)
        self.vce = (vce or d.vce or "linearized")
        self.mse = d.mse if mse is None else mse
        self.rep_vars = d.repweights if self.vce == d.vce else []
        self.hadamard = hadamard or None
        if hadamard is not False and self.vce == "brr" and not self.rep_vars and not hadamard:
            # BRR computed by OpenDTA needs the Hadamard matrix, as in Stata (compat 0602): svy
            # brr, hadamard(matrix) ...
            raise StataError(198, "option hadamard() required when brrweights are not svyset")
        R = None
        if self.rep_vars:
            R = np.column_stack([ds.get(v).data.astype(float) for v in self.rep_vars])
            mask = mask & np.all(R < SYS, axis=1)
        self.mask = mask
        self.w = w[mask]
        self.rep = R[mask] if R is not None else None
        self.strata = strata[mask]
        self.psu = psu[mask]
        self.fpc = ds.get(d.fpc).data[mask].astype(float) if d.fpc else None
        self.n = int(mask.sum())
        # following stages (PSU within stratum within the previous unit)
        self.stage_cols = []
        for st_ in d.stages:
            su = ds.get(st_["su"]).data[mask] if st_["su"] else np.arange(n_all, dtype=float)[mask]
            sh = ds.get(st_["strata"]).data[mask] if st_["strata"] else np.zeros(self.n)
            fp = ds.get(st_["fpc"]).data[mask].astype(float) if st_["fpc"] else None
            self.stage_cols.append((sh, su, fp))
        # PSU identifiers within stratum
        self.h_ids, self.h_inv = np.unique(self.strata, return_inverse=True)
        pairs = np.column_stack([self.h_inv, self.psu])
        _, self.u_inv = np.unique(pairs, axis=0, return_inverse=True)
        self.u_inv = self.u_inv.ravel()
        self.n_strata = len(self.h_ids)
        self.n_psu = int(self.u_inv.max()) + 1 if self.n else 0
        # design degrees of freedom (compat 0603): replicate weights without a design, number of
        # replicates - 1; BRR, also number of replicates - 1 (4 Hadamard replicates: df 3);
        # bootstrap and SDR replicate weights: no design df (normal-based intervals)
        self.normal = bool(self.rep is not None and self.vce in ("bootstrap", "sdr"))
        if self.rep is not None and not d.psu and not d.strata:
            self.df = self.rep.shape[1] - 1
        elif self.vce == "brr":
            nrep = None
            if self.rep is not None:
                nrep = self.rep.shape[1]
            elif self.hadamard:
                from ..commands.matrix import get_matrix
                try:
                    nrep = np.asarray(get_matrix(s, str(self.hadamard).strip()).data).shape[0]
                except StataError:
                    nrep = None
            self.df = (nrep - 1) if nrep else self.n_psu - self.n_strata     # VERIFY with brrweight()
        else:
            self.df = self.n_psu - self.n_strata
        if d.dof:
            self.df = int(d.dof)

    # -- linearization ------------------------------------------------------
    def variance(self, z: np.ndarray) -> np.ndarray:
        """V from the scores z (n x k) of the sample observations."""
        z = np.atleast_2d(z.T).T if z.ndim == 1 else z
        k = z.shape[1]
        U = np.zeros((self.n_psu, k))
        np.add.at(U, self.u_inv, z)
        psu_h = np.zeros(self.n_psu, dtype=int)
        psu_h[self.u_inv] = self.h_inv
        V = np.zeros((k, k))
        single = False
        grand = U.mean(axis=0) if self.n_psu else np.zeros(k)
        for h in range(self.n_strata):
            idx = np.nonzero(psu_h == h)[0]
            nh = len(idx)
            f = 0.0
            if self.fpc is not None:
                fv = self.fpc[self.h_inv == h][0]
                f = fv if fv <= 1 else nh / fv
            if nh < 2:
                mode = self.d.singleunit
                if mode == "missing":
                    single = True
                elif mode == "centered":
                    dev = U[idx] - grand
                    V += (1 - f) * dev.T @ dev
                elif mode == "scaled":
                    single = True   # VERIFY: scaling by the mean of the other strata
                continue
            dev = U[idx] - U[idx].mean(axis=0)
            V += (1 - f) * nh / (nh - 1) * dev.T @ dev
            # following stages: they only count if the 1st stage is without replacement (f_h > 0),
            # multiplied by f_h ([SVY] variance estimation)
            if f > 0 and self.stage_cols:
                for u in idx:
                    rows = np.nonzero(self.u_inv == u)[0]
                    V += f * self._stage_var(z, rows, 0)
        if single:
            return np.full((k, k), np.nan)
        return V

    def _stage_var(self, z: np.ndarray, rows: np.ndarray, level: int) -> np.ndarray:
        k = z.shape[1]
        sh, su, fp = self.stage_cols[level]
        V = np.zeros((k, k))
        for h in np.unique(sh[rows]):
            rh = rows[sh[rows] == h]
            units = np.unique(su[rh])
            nh = len(units)
            if nh < 2:
                continue
            U = np.array([z[rh[su[rh] == u]].sum(axis=0) for u in units])
            f = 0.0
            if fp is not None:
                fv = fp[rh][0]
                f = fv if fv <= 1 else nh / fv
            dev = U - U.mean(axis=0)
            V += (1 - f) * nh / (nh - 1) * dev.T @ dev
            if f > 0 and level + 1 < len(self.stage_cols):
                for u in units:
                    V += f * self._stage_var(z, rh[su[rh] == u], level + 1)
        return V

    # -- replication ----------------------------------------------------------
    @property
    def replicated(self) -> bool:
        return self.vce in ("jackknife", "brr", "bootstrap", "sdr")

    @property
    def vcetype(self) -> str:
        if self.mse and self.vce != "linearized":
            # with mse Stata abbreviates and marks with *: "Jknife *" (compat 0602); the other
            # names with * are VERIFY
            return {"jackknife": "Jknife *", "brr": "BRR *", "bootstrap": "Bootstrap *",
                    "sdr": "SDR *"}[self.vce]
        return {"linearized": "Linearized", "jackknife": "Jackknife", "brr": "BRR", "bootstrap": "Bootstrap",
                "sdr": "SDR"}[self.vce]

    def replicates(self):
        """Yields (replicate weights, stratum, multiplier) for the method."""
        w = self.w
        if self.rep is not None:
            for r in range(self.rep.shape[1]):
                yield self.rep[:, r], 0, None
            return
        psu_h = np.zeros(self.n_psu, dtype=int)
        psu_h[self.u_inv] = self.h_inv
        if self.vce == "jackknife":
            for h in range(self.n_strata):
                units = np.nonzero(psu_h == h)[0]
                nh = len(units)
                if nh < 2:
                    continue
                in_h = self.h_inv == h
                for u in units:
                    ww = w.copy()
                    ww[in_h] *= nh / (nh - 1)
                    ww[self.u_inv == u] = 0.0
                    yield ww, h, (nh - 1) / nh
            return
        if self.vce == "brr":
            counts = np.bincount(psu_h, minlength=self.n_strata)
            if np.any(counts != 2):
                raise StataError(459, "BRR requires 2 PSUs per stratum")   # VERIFY
            from ..commands.matrix import get_matrix
            H = np.asarray(get_matrix(self.s, str(self.hadamard).strip()).data, dtype=float)
            if H.shape[1] < self.n_strata:
                raise StataError(198, "Hadamard matrix has too few columns")   # VERIFY
            R = H.shape[0]
            fay = self.d.fay
            for r in range(R):
                ww = w.copy()
                for h in range(self.n_strata):
                    units = np.nonzero(psu_h == h)[0]
                    first, second = units[0], units[1]
                    a, b = (2 - fay, fay) if H[r, h] > 0 else (fay, 2 - fay)
                    ww[self.u_inv == first] *= a
                    ww[self.u_inv == second] *= b
                yield ww, 0, None
            return
        raise StataError(198, f"vce({self.vce}) requires replicate weights")   # VERIFY

    def rep_variance(self, estimator, theta0: np.ndarray) -> tuple[np.ndarray, int]:
        """V by replication; `estimator(w)` returns the vector of estimates."""
        thetas, strata, mults = [], [], []
        reps = list(self.replicates())
        out = self.s.output
        name = {"jackknife": "Jackknife", "brr": "BRR", "bootstrap": "Bootstrap",
                "sdr": "SDR"}.get(self.vce, self.vce)
        # replication dots, as with the jackknife prefix; VERIFY for brr, bootstrap and sdr (in
        # 0602 only the jackknife got to run in Stata)
        out.write(f"\n{name} replications ({len(reps)})\n"
                  "----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 \n", "text")
        for k, (ww, h, m) in enumerate(reps):
            try:
                th = np.asarray(estimator(ww), dtype=float)
            except StataError:
                th = np.full(len(theta0), np.nan)
            thetas.append(th)
            strata.append(h)
            mults.append(m)
            out.write("x" if np.isnan(th).any() else ".", "text")
            if (k + 1) % 50 == 0:
                out.write(f"{k + 1:>6}\n", "text")
        if len(reps) % 50:
            out.write("\n", "text")
        # the estimates of each replicate are stored as float, as in Stata (svy
        # jackknife/brr/bootstrap/sdr only use double with the double option; [SVY] svy jackknife
        # and [R] jackknife) — compat 0602
        T = np.array(thetas).astype(np.float32).astype(float)
        R = len(T)
        k = len(theta0)
        ok = ~np.isnan(T).any(axis=1)
        if self.vce == "jackknife":
            V = np.zeros((k, k))
            hs = np.array(strata)
            for h in np.unique(hs):
                sel = (hs == h) & ok
                if not sel.any():
                    continue
                Th = T[sel]
                c = theta0 if self.mse else Th.mean(axis=0)
                m = mults[int(np.nonzero(hs == h)[0][0])]
                if m is None:
                    m = (R - 1) / R
                f = 0.0
                if self.fpc is not None and self.rep is None:
                    fv = self.fpc[self.h_inv == h][0]
                    nh = int(sel.sum())
                    f = fv if fv <= 1 else nh / fv
                D = Th - c
                V += (1 - f) * m * D.T @ D
            return V, R
        Tg = T[ok]
        c = theta0 if self.mse else Tg.mean(axis=0)
        D = Tg - c
        S = D.T @ D
        if self.vce == "brr":
            return S / (len(Tg) * (1 - self.d.fay) ** 2), R
        if self.vce == "bootstrap":
            return self.d.bsn / len(Tg) * S, R
        if self.vce == "sdr":
            return 4 / len(Tg) * S, R
        raise StataError(198, "invalid vce")


# ---------------------------------------------------------------------------
# estimadores: mean, proportion, total, ratio
# ---------------------------------------------------------------------------

def _over_groups(s, over: str, mask: np.ndarray):
    if not over:
        return [("", np.ones(int(mask.sum()), bool), "")]
    from ..core.varlist import expand
    ds = s.data
    names = expand(ds, over)
    cols = [ds.get(n).data[mask] for n in names]
    keys = sorted(set(zip(*[c.tolist() for c in cols])))
    groups = []
    for key in keys:
        sel = np.ones(int(mask.sum()), bool)
        labs = []
        for n, c, k_ in zip(names, cols, key):
            sel &= c == k_
            v = ds.get(n)
            lab = ds.value_labels.get(v.value_label, {}).get(int(k_)) if v.value_label else None
            labs.append(lab if lab else (f"{int(k_)}" if k_ == int(k_) else f"{k_:g}"))
        groups.append((" ".join(labs), sel, key))
    return groups


def _estimate(stat: str, Y: list[np.ndarray], w: np.ndarray, groups, names: list[str], X=None):
    """Values and scores (n x k) of each statistic, by over() group."""
    vals, scores, labels = [], [], []
    n = len(w)
    for j, y in enumerate(Y):
        for gname, sel, _ in groups:
            ws = w * sel
            if stat == "mean":
                W = ws.sum()
                m = float(ws @ y) / W if W > 0 else SYS
                z = ws * (y - m) / W if W > 0 else np.zeros(n)
            elif stat == "total":
                m = float(ws @ y)
                z = ws * y
            elif stat == "ratio":
                x = X[j]
                Wx = float(ws @ x)
                m = float(ws @ y) / Wx if Wx else SYS
                z = ws * (y - m * x) / Wx if Wx else np.zeros(n)
            else:   # proportion: y is already an indicator
                W = ws.sum()
                m = float(ws @ y) / W if W > 0 else SYS
                z = ws * (y - m) / W if W > 0 else np.zeros(n)
            vals.append(m)
            scores.append(z)
            labels.append((names[j], gname))
    return np.array(vals), np.column_stack(scores) if scores else np.zeros((n, 0)), labels


_TITLES = {"mean": "Mean estimation", "proportion": "Proportion estimation", "total": "Total estimation",
           "ratio": "Ratio estimation"}


def _summary_table(s, stat_title: str, labels, vals, se, df, level, vcetype: str, over_names: str):
    out = s.output
    lv = f"{level:g}"
    width = 62
    crit = st.t.ppf(1 - (1 - level / 100) / 2, df) if df and df > 0 else st.norm.ppf(1 - (1 - level / 100) / 2)
    out.write("-" * width + "\n", "text")
    if vcetype == "Bootstrap" and not (df and df > 0):
        # bootstrap replicate weights: observed value and normal-based interval (compat 0603)
        out.write(f"{'':>12} |   Observed   Bootstrap         Normal-based\n", "text")
    elif vcetype:
        start = 28 + (9 - len(vcetype)) // 2
        out.write(f"{'':>12} |" + " " * (start - 14) + vcetype + "\n", "text")
    head = "Over" if over_names else ""
    out.write(f"{head:>12} | {stat_title:>10}   Std. Err.     [{lv}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 48 + "\n", "text")
    last = None
    for (var, grp), b, e in zip(labels, vals, se):
        if grp:
            if var != last:
                # variable header on the left; between over() blocks, a separator line (compat
                # 0601)
                if last is not None and over_names:
                    out.write("-" * 13 + "+" + "-" * 48 + "\n", "text")
                out.write(f"{var[:12]:<12} |" + "\n", "text")
                last = var
            lab = grp
        else:
            lab = var
        out.write(f"{lab[:12]:>12} |", "text")
        if e is None or not np.isfinite(e):
            out.write(f"  {g9(b):>9}  {'.':>9}     {'.':>9}   {'.':>9}\n", "result")
            continue
        if stat_title == "Proportion" and 0 < b < 1 and vcetype in ("Linearized", ""):
            # CI on the logit scale (proportion's default in Stata 14, VERIFY)
            from scipy.special import expit, logit
            half = crit * e / (b * (1 - b))
            lo, hi = expit(logit(b) - half), expit(logit(b) + half)
        else:
            lo, hi = b - crit * e, b + crit * e
        out.write(f"  {g9(b):>9}  {g9(e):>9}     {g9(lo):>9}   {g9(hi):>9}\n", "result")
    out.write("-" * width + "\n", "text")


def _proportion_vars(s, names, mask):
    """proportion: one indicator per category of each variable."""
    ds = s.data
    Y, labs = [], []
    for n in names:
        v = ds.get(n)
        x = v.data[mask]
        for lv in sorted(set(x.tolist())):
            Y.append((x == lv).astype(float))
            lab = ds.value_labels.get(v.value_label, {}).get(int(lv)) if v.value_label else None
            labs.append((n, lab or (f"{int(lv)}" if lv == int(lv) else f"{lv:g}")))
    return Y, labs


def summary_stat(s: "Session", stat: str, args: str, *, svy: bool = False, svyopts: str = "") -> None:
    from ..commands._util import touse
    from ..core.varlist import expand
    p = parse_standard(args)
    o = match_options(p.options, {"over": 2, "level": 1, "vce": 3, "nolegend": 4, "noheader": 6,
                                  "missing": 4, "percent": 3, "citype": 3}) if p.options.strip() else {}
    ds = s.data
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    mask = touse(s, p)
    if stat == "ratio":
        pairs = re.findall(r"(?:\(?\s*(\w+)\s*:)?\s*(\w+)\s*/\s*(\w+)\s*\)?", p.varlist)
        if not pairs:
            raise StataError(198, "invalid syntax")
        rnames = [nm or f"_ratio_{k + 1}" for k, (nm, _a, _b) in enumerate(pairs)]
        names = list(rnames)
        num = [expand(ds, a)[0] for _, a, _b in pairs]
        den = [expand(ds, b)[0] for _, _a, b in pairs]
        allv = num + den
    else:
        names = expand(ds, p.varlist, allow_empty=False)
        allv = names
    for n in allv:
        x = ds.get(n).data
        mask &= x < SYS
    over = str(o.get("over") or "").split(",")[0].strip()
    if over:
        for n in expand(ds, over):
            mask &= ds.get(n).data < SYS
    out = s.output
    d = require_design(s) if svy else None
    subpop = None
    if svy:
        so = parse_svy_opts(svyopts)
        dd = DesignData(s, d, np.ones(ds.nobs, bool) if not p.if_ and not p.in_ else touse(s, p),
                        vce=so.get("vce"), mse=so.get("mse"), hadamard=so.get("hadamard"))
        if so.get("subpop"):
            subpop = _subpop_mask(s, str(so["subpop"]))
        est_mask = mask & dd.mask
        out.write(f"\n(running {stat} on estimation sample)\n", "text")
        # scores outside the sample (or subpop) are zero, but the design counts all PSUs
        full = dd.mask
        inside = est_mask[full] & (subpop[full] if subpop is not None else True)
        w = dd.w * inside
        rows_idx = np.nonzero(full)[0]
        n_obs = int(inside.sum())
    else:
        full = mask
        w_all = np.ones(ds.nobs)
        if p.weight:
            from ..commands.summarize import weights as get_weights
            w_all, wtype, mask = get_weights(s, p, mask, ("pweight", "aweight", "fweight", "iweight"))
            full = mask
        w = (w_all[full] if p.weight else np.ones(int(full.sum()))).astype(float)
        rows_idx = np.nonzero(full)[0]
        n_obs = int(full.sum()) if not (p.weight and p.weight[0].startswith("f")) else int(w.sum())
    if stat == "proportion":
        Y, labs_p = _proportion_vars(s, names, full)
        names_eff = [f"{a}" for a, _ in labs_p]
    else:
        Y = [ds.get(n).data[full].astype(float) for n in (num if stat == "ratio" else names)]
        Y = [np.where(y < SYS, y, 0.0) for y in Y]
    X = [np.where(ds.get(n).data[full] < SYS, ds.get(n).data[full], 0.0) for n in den] if stat == "ratio" else None
    groups = _over_groups(s, over, full) if over else [("", np.ones(len(w), bool), "")]
    legend: list[str] = []
    if over:
        ovars = over.split()
        seen = []
        for gname, _sel, key in groups:
            if gname in seen:
                continue
            seen.append(gname)
            desc = ", ".join(f"{v} = {lab}" for v, lab in zip(ovars, gname.split(" "))) \
                if len(ovars) > 1 else f"{ovars[0]} = {gname}"
            legend.append(f"{gname[:13]:>13}: {desc}")
    def compute(wv):
        if stat == "proportion":
            vals_, Z_all, labels_ = [], [], []
            for i, y in enumerate(Y):
                for gname, sel, _ in groups:
                    v_, Z_, _l = _estimate("proportion", [y], wv, [(gname, sel, None)], [names_eff[i]])
                    vals_.extend(v_.tolist())
                    Z_all.append(Z_)
                    labels_.append((labs_p[i][0], f"{labs_p[i][1]}" + (f" {gname}" if over else "")))
            return np.array(vals_), (np.column_stack(Z_all) if Z_all else np.zeros((len(wv), 0))), labels_
        return _estimate(stat, Y, wv, groups, names, X)

    vals, Z, labels = compute(w)
    # over() group of each column (for the "analytic" variance without svy)
    colsel = [sel for _y in Y for _g, sel, _k in groups]
    if svy:
        n_rep = None
        if dd.replicated:
            V, n_rep = dd.rep_variance(lambda ww: compute(ww * inside)[0], vals)
        else:
            V = dd.variance(Z)
        df = None if dd.normal else dd.df
        vcetype = dd.vcetype
        # variances under SRS for estat effects ([SVY] estat, Methods and formulas)
        _srs = _srs_variances(Z, w, dd, int(est_mask[full].sum()))
    else:
        # V = N/(N-1) Σ z z' (with fweight, each score counts w times)
        robust = str(o.get("vce", "")).startswith("r")
        pw = bool(p.weight and p.weight[0] == "pweight")
        # with pweight the standard error is linearized, but Stata shows no label
        vcetype = "Robust" if robust else ""
        fw = bool(p.weight and p.weight[0] == "fweight")
        cnt = w if fw else np.ones(len(w))
        nn = float(cnt.sum())
        Zr = Z / np.where(w > 0, w, 1)[:, None] if fw else Z
        if stat == "total" or stat == "ratio" or robust or pw:
            # linearization: centered scores, V = n/(n-1) Σ (z - z̄)(z - z̄)'
            zbar = (cnt @ Zr) / nn if stat == "total" else np.zeros(Z.shape[1])
            Zc = Zr - zbar
            V = nn / (nn - 1) * (Zc * cnt[:, None]).T @ Zc
        else:
            # "analytic" mean/proportion: s²/n within each over() group, with no covariance between
            # groups (compat 0601: .1209974 = sd/sqrt(n_g))
            V = (Zr * cnt[:, None]).T @ Zr
            Wg = np.array([float(cnt @ sel) for sel in colsel])
            c = np.sqrt(np.where(Wg > 1, Wg / (Wg - 1), np.nan))
            V = V * np.outer(c, c)
        df = int(nn) - 1
        if stat == "ratio":
            vcetype = "Linearized"
    se = np.sqrt(np.maximum(np.diag(V), 0)) if V.size else np.zeros(0)
    if np.isnan(V).any():
        se = np.full(len(vals), np.nan)
    # header
    if not o.get("noheader"):
        if svy:
            pop = float(dd.w[est_mask[full]].sum())
            out.write("\n", "text")
            right = [("Number of obs", f"{int(est_mask[full].sum()):,}"), ("Population size", _num(pop))]
            if subpop is not None:
                right += [("Subpop. no. obs", f"{n_obs:,}"), ("Subpop. size", _num(float(w.sum())))]
            if n_rep:
                right.append(("Replications", f"{n_rep:,}"))
            if df is not None:
                right.append(("Design df", f"{df:,}"))
            # replicate weights (compat 0603): jackknife shows the strata but not the PSUs;
            # bootstrap and SDR show neither, and the title shares the line with Number of obs
            if dd.rep is not None and dd.normal:
                left = []
            elif dd.rep is not None:
                left = [f"Number of strata = {dd.n_strata:>7,}", ""]
            else:
                left = [f"Number of strata = {dd.n_strata:>7,}", f"Number of PSUs   = {dd.n_psu:>7,}"]
            if left:
                out.write(f"Survey: {_TITLES[stat]}\n\n", "text")
            else:
                left = [f"Survey: {_TITLES[stat]}"]
            for i, (rl, rv) in enumerate(right):
                lft = left[i] if i < len(left) else ""
                out.write(f"{lft:<34}{rl:<16}= ", "text")
                out.write(f"{rv:>10}\n", "result")
        else:
            out.write("\n", "text")
            out.write(f"{_TITLES[stat]:<34}{'Number of obs':<16}= ", "text")
            out.write(f"{n_obs:>10,}\n", "result")
        if over:
            # group legend: "         masc: sexo = masc" (compat 0601)
            out.write("\n", "text")
            for ln in legend:
                out.write(ln + "\n", "text")
        if stat == "ratio":
            out.write("\n", "text")
            for rn, (num_, den_) in zip(rnames, zip(num, den)):
                out.write(f"{rn[:13]:>13}: {num_}/{den_}\n", "text")
        out.write("\n", "text")
    stat_title = {"mean": "Mean", "proportion": "Proportion", "total": "Total", "ratio": "Ratio"}[stat]
    _summary_table(s, stat_title, labels, vals, se, df, level, vcetype, over)
    # e()
    from ..commands.matrix import Matrix
    from .results import Estimates, post
    cn = [f"{a}:{b}" if b else a for a, b in labels]
    est = Estimates(stat, "", [c.split(":")[-1] if False else c for c in cn], vals, V, n_obs,
                    "t" if svy else "t", df)
    est.rows = []
    est.level = level
    est.display = lambda ss, ee, **kw: _summary_table(ss, stat_title, labels, vals, se, df, level, vcetype, over)
    if svy:
        est.svy_effects = dict(V=V, srs=_srs, labels=labels, vals=vals, title=stat_title, vcetype=vcetype,
                               over=over, level=level, legend=legend)
    if svy:
        # e() of svy in Stata 14, in ereturn list order (compat 0601)
        k = len(vals)
        n_over = len(groups) if over else 1
        scalars = [("df_r", float(df) if df is not None else SYS)]
        if subpop is not None:
            scalars += [("N_subpop", float(w.sum())), ("N_sub", float(n_obs))]
        scalars += [("N_strata_omit", 0.0), ("singleton", 0.0), ("census", 0.0),
                    ("N_pop", float(dd.w[est_mask[full]].sum())), ("N_psu", float(dd.n_psu)),
                    ("N_strata", float(dd.n_strata)), ("N_over", float(n_over)),
                    ("N", float(int(est_mask[full].sum()))), ("stages", 1.0), ("k_eq", 1.0),
                    ("rank", float(np.linalg.matrix_rank(V)) if V.size else 0.0)]
        wvar = d.weight
        cmdtxt = f"{stat} {args.strip()}".strip()
        command = cmdtxt if "," in cmdtxt else cmdtxt + ","
        macros = [("cmdline", f"svy{', ' + svyopts.strip() if svyopts.strip() else ' '}: {cmdtxt}"),
                  ("cmd", stat), ("prefix", "svy"), ("cmdname", stat), ("command", command),
                  ("title", f"Survey: {_TITLES[stat]}"), ("vcetype", vcetype), ("vce", dd.vce),
                  ("estat_cmd", "svy_estat"), ("varlist", " ".join(allv))]
        if over:
            macros += [("over", over), ("over_labels", " ".join(g[0] for g in groups)),
                       ("over_namelist", " ".join(g[0] for g in groups))]   # VERIFY
        if subpop is not None:
            macros.append(("subpop", str(so["subpop"]).strip()))
        macros += [("marginsnotok", "_ALL")]
        if wvar:
            macros += [("wtype", d.wtype or "pweight"), ("wvar", wvar), ("wexp", f"= {wvar}")]
        macros.append(("singleunit", d.singleunit))
        if d.psu and d.psu != "_n":
            macros.append(("su1", d.psu))
        if d.strata:
            macros.append(("strata1", d.strata))
        if d.fpc:
            macros.append(("fpc1", d.fpc))
        macros += [("properties", "b V"), ("depvar", stat_title)]
        Nk = np.array([[float(np.count_nonzero(Z[:, j]) or n_obs) for j in range(k)]])
        Wk = np.array([[float(np.sum(w * (np.abs(Z[:, j]) > 0))) or float(w.sum()) for j in range(k)]])
        nn_ = max(n_obs, 2)
        Vsrs = nn_ / (nn_ - 1) * Z.T @ Z
        mats = [("_N_subp", Matrix(Wk, ["r1"], list(cn))), ("V_srs", Matrix(Vsrs, list(cn), list(cn))),
                ("_N", Matrix(Nk, ["r1"], list(cn))),
                ("_N_strata_certain", Matrix(np.zeros((1, 1)))), ("_N_strata_single", Matrix(np.zeros((1, 1)))),
                ("_N_strata", Matrix(np.array([[float(dd.n_strata)]]))),
                ("error", Matrix(np.zeros((1, k)), ["r1"], list(cn)))]
        post(s, est, scalars, macros, np.isin(np.arange(ds.nobs), rows_idx), extra_mats=mats)
        return
    scalars = [("N", float(n_obs)), ("df_r", float(df) if df is not None else SYS)]
    macros = [("cmd", stat), ("cmdline", f"{stat} {args.strip()}".strip()), ("properties", "b V"),
              ("varlist", " ".join(allv))]
    post(s, est, scalars, macros, np.isin(np.arange(ds.nobs), rows_idx))


def _num(x: float) -> str:
    return f"{x:,.0f}" if abs(x - round(x)) < 1e-6 else format_value(x, "%11.0g", pad=False).strip()


def parse_svy_opts(text: str) -> dict:
    """Options of the svy prefix: subpop(), vce(), mse, dof()."""
    so = match_options(text, {"subpop": 3, "vce": 3, "mse": 3, "dof": 3, "nodots": 5, "dots": 4,
                              "level": 1, "hadamard": 3, "fay": 3}) if text.strip() else {}
    out = {}
    if so.get("subpop"):
        out["subpop"] = str(so["subpop"])
    if so.get("vce"):
        v = str(so["vce"]).strip().lower()
        out["vce"] = {"linear": "linearized", "jack": "jackknife", "bs": "bootstrap",
                      "boot": "bootstrap"}.get(v, v)
    if so.get("mse"):
        out["mse"] = True
    out["hadamard"] = str(so["hadamard"]) if so.get("hadamard") else None
    return out


def _srs_variances(Z: np.ndarray, w: np.ndarray, dd: "DesignData", m: int) -> tuple:
    """V_srswor and V_srswr from [SVY] estat (total formula with the scores
    u = z/w): V = (1 - f) M̂/(m - 1) Σ w (u - ū)(u - ū)'."""
    with np.errstate(divide="ignore", invalid="ignore"):
        U = np.where(w[:, None] > 0, Z / np.where(w > 0, w, 1)[:, None], 0.0)
    Mhat = float(dd.w.sum())
    ubar = (w @ U) / Mhat if Mhat else np.zeros(Z.shape[1])
    D = U - ubar
    S = (D * w[:, None]).T @ D
    Vwr = Mhat / max(m - 1, 1) * S
    f = 0.0
    if dd.fpc is not None:
        fv = float(dd.fpc[0])
        f = fv if fv <= 1 else m / fv   # VERIFY: f = m/M̂ with fpc in the 1st stage
    return (1 - f) * Vwr, Vwr


def _subpop_mask(s: "Session", text: str) -> np.ndarray:
    from ..commands._util import eval_vector
    t = text.strip()
    ds = s.data
    if t.startswith("if "):
        v = np.asarray(eval_vector(s, t[3:]), dtype=float)
        return np.broadcast_to(v, (ds.nobs,)) != 0
    m = re.match(r"^(\w+)(?:\s+if\s+(.*))?$", t)
    if not m:
        raise StataError(198, "invalid subpop() option")
    from ..core.varlist import resolve_name
    x = ds.get(resolve_name(ds, m.group(1))).data
    sel = (x != 0) & (x < SYS)
    if m.group(2):
        v = np.asarray(eval_vector(s, m.group(2)), dtype=float)
        sel &= np.broadcast_to(v, (ds.nobs,)) != 0
    return sel


# ---------------------------------------------------------------------------
# svydescribe
# ---------------------------------------------------------------------------

def cmd_svydescribe(s: "Session", args: str) -> None:
    d = require_design(s)
    ds = s.data
    from ..commands._util import touse
    p = parse_standard(args)
    mask = touse(s, p)
    dd = DesignData(s, d, mask)
    out = s.output
    out.write("\nSurvey: Describing stage 1 sampling units\n", "text")
    _show_design(s, d)
    out.write("\n" + " " * 38 + "#Obs per Unit\n", "text")
    out.write(" " * 30 + "-" * 28 + "\n", "text")
    out.write("Stratum    #Units     #Obs      min       mean      max   \n", "text")
    out.write("  ".join(["-" * 8] * 6) + "\n", "text")
    psu_h = np.zeros(dd.n_psu, dtype=int)
    psu_h[dd.u_inv] = dd.h_inv
    counts = np.bincount(dd.u_inv, minlength=dd.n_psu)
    for h, hv in enumerate(dd.h_ids):
        c = counts[psu_h == h]
        lab = f"{int(hv)}" if hv == int(hv) else f"{hv:g}"
        out.write(f"{lab:>8}  {len(c):>8,}  {int(c.sum()):>8,}  {int(c.min()):>8,}  {c.mean():>8.1f}  "
                  f"{int(c.max()):>8,}\n", "result")
    out.write("  ".join(["-" * 8] * 6) + "\n", "text")
    out.write(f"{dd.n_strata:>8,}  {dd.n_psu:>8,}  {dd.n:>8,}  {int(counts.min()):>8,}  {counts.mean():>8.1f}  "
              f"{int(counts.max()):>8,}\n", "result")


# ---------------------------------------------------------------------------
# prefixo svy:
# ---------------------------------------------------------------------------

def cmd_svy(s: "Session", args: str) -> None:
    """svy [, subpop() vce()]: comando"""
    head, colon, cmdtext = args.partition(":")
    if not colon:
        raise StataError(198, "invalid syntax")
    opts = head.strip().lstrip(",").strip()
    cmdtext = cmdtext.strip()
    word, _, rest = cmdtext.partition(" ")
    from ..commands.registry import lookup
    spec = lookup(word)
    name = spec.name if spec else word
    if name in ("mean", "proportion", "total", "ratio"):
        summary_stat(s, name, rest, svy=True, svyopts=opts)
        return
    if name in ("regress", "logit", "logistic", "probit", "poisson"):
        from .svyreg import svy_model
        svy_model(s, name, rest, opts)
        return
    if name in ("ologit", "oprobit", "mlogit", "nbreg", "glm", "tobit"):
        from .svyreg import svy_generic
        svy_generic(s, name, rest, opts)
        return
    if name in ("tabulate",):
        from .svyreg import svy_tabulate
        svy_tabulate(s, rest, opts)
        return
    raise StataError(199, f"svy: {word} not yet supported by OpenDTA")


def cmd_mean(s, args):
    summary_stat(s, "mean", args)


def cmd_proportion(s, args):
    summary_stat(s, "proportion", args)


def cmd_total(s, args):
    summary_stat(s, "total", args)


def cmd_ratio(s, args):
    summary_stat(s, "ratio", args)


# ---------------------------------------------------------------------------
# estat effects
# ---------------------------------------------------------------------------

def _g6(x: float) -> str:
    """DEFF and DEFT with six significant digits, without the leading zero
    (.64512, 1.05667, .803194; compat 0602)."""
    if x is None or x >= SYS or not np.isfinite(x):
        return "."
    t = f"{x:.6g}"
    if t.startswith("0."):
        t = t[1:]
    elif t.startswith("-0."):
        t = "-" + t[2:]
    return t


def estat_effects(s: "Session", est, rest: str) -> None:
    """DEFF = V/V_srswor and DEFT = sqrt(V/V_srswr) ([SVY] estat). Layout
    checked in compat 0602."""
    info = getattr(est, "svy_effects", None)
    if info is None:
        raise StataError(321, "estat effects not allowed after this command")   # VERIFY
    o = match_options(rest.lstrip(", "), {"deff": 4, "deft": 4, "srssubpop": 3, "meff": 4, "meft": 4}) \
        if rest.strip(" ,") else {}
    V = np.asarray(info["V"], dtype=float)
    Vwor, Vwr = info["srs"]
    out = s.output
    width = 58
    if info.get("legend"):
        # same group legend as the command (compat 0602)
        out.write("\n" + "\n".join(info["legend"]) + "\n", "text")
    out.write("\n" + "-" * width + "\n", "text")
    vt = info["vcetype"]
    start = 28 + (9 - len(vt)) // 2
    out.write(f"{'':>12} |" + " " * (start - 14) + vt + "\n", "text")
    head = "Over" if info.get("over") else info.get("depvar", "")
    out.write(f"{head[:12]:>12} | {info['title']:>10}   Std. Err.       DEFF      DEFT\n", "text")
    out.write("-" * 13 + "+" + "-" * 44 + "\n", "text")
    last = None
    deffs, defts = [], []
    for j, ((var, grp), b) in enumerate(zip(info["labels"], info["vals"])):
        v = V[j, j]
        deff = v / Vwor[j, j] if Vwor[j, j] > 0 else SYS
        deft = float(np.sqrt(v / Vwr[j, j])) if Vwr[j, j] > 0 else SYS
        deffs.append(deff)
        defts.append(deft)
        if grp:
            if var != last:
                if last is not None:
                    out.write("-" * 13 + "+" + "-" * 44 + "\n", "text")
                out.write(f"{var[:12]:<12} |\n", "text")
                last = var
            lab = grp
        else:
            lab = var
        out.write(f"{lab[:12]:>12} |", "text")
        out.write(f"  {g9(b):>9}  {g9(float(np.sqrt(v))):>9}{_g6(deff):>12}{_g6(deft):>10}\n", "result")
    out.write("-" * width + "\n", "text")
    from ..commands.matrix import Matrix
    cn = [f"{a}:{b}" if b else a for a, b in info["labels"]]
    s.r = {"deft": Matrix(np.array([defts]), ["r1"], cn), "deff": Matrix(np.array([deffs]), ["r1"], cn)}
