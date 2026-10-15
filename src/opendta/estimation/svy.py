"""Complex surveys ([SVY] svyset, svydescribe, the svy prefix and estat effects).

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
# mean, proportion, ratio and total under svy (the estimators are in summaries.py)
# ---------------------------------------------------------------------------

class SurveyEngine:
    """svy: mean/proportion/ratio/total: design weights, design variance (linearized or
    replicated), the survey header and the e() of svy."""

    def __init__(self, s: "Session", req, svyopts: str = ""):
        from ..commands._util import touse
        self.design = d = require_design(s)
        self.svyopts = svyopts
        self.so = so = parse_svy_opts(svyopts)
        p = req.p
        all_obs = np.ones(s.data.nobs, bool) if not p.if_ and not p.in_ else touse(s, p)
        self.dd = dd = DesignData(s, d, all_obs, vce=so.get("vce"), mse=so.get("mse"),
                                  hadamard=so.get("hadamard"))
        self.subpop = subpop_mask(s, str(so["subpop"])) if so.get("subpop") else None
        s.output.write(f"\n(running {req.stat} on estimation sample)\n", "text")
        # scores outside the sample (or subpop) are zero, but the design counts all PSUs
        self.rows = dd.mask
        self.in_sample = (req.mask & dd.mask)[dd.mask]
        self.inside = self.in_sample & (self.subpop[dd.mask] if self.subpop is not None else True)
        self.w = dd.w * self.inside
        self.n_obs = int(self.inside.sum())
        self.n_sample = int(self.in_sample.sum())
        self.n_rep = None

    def variance(self, req, fit, refit):
        dd = self.dd
        if dd.replicated:
            V, self.n_rep = dd.rep_variance(lambda ww: refit(ww * self.inside), fit.vals)
        else:
            V = dd.variance(fit.Z)
        # variances under SRS for estat effects ([SVY] estat, Methods and formulas)
        self.srs = srs_variances(fit.Z, self.w, dd, self.n_sample)
        return V, None if dd.normal else dd.df, dd.vcetype

    def header(self, s, req, df):
        dd, out = self.dd, s.output
        title = f"Survey: {req.info.title}"
        out.write("\n", "text")
        right = [("Number of obs", f"{self.n_sample:,}"), ("Population size", format_size(self.population))]
        if self.subpop is not None:
            right += [("Subpop. no. obs", f"{self.n_obs:,}"), ("Subpop. size", format_size(float(self.w.sum())))]
        if self.n_rep:
            right.append(("Replications", f"{self.n_rep:,}"))
        if df is not None:
            right.append(("Design df", f"{df:,}"))
        # replicate weights (compat 0603): jackknife shows the strata but not the PSUs; bootstrap
        # and SDR show neither, and the title shares the line with Number of obs
        if dd.rep is not None and dd.normal:
            left = [title]
        else:
            out.write(f"{title}\n\n", "text")
            left = [f"Number of strata = {dd.n_strata:>7,}",
                    "" if dd.rep is not None else f"Number of PSUs   = {dd.n_psu:>7,}"]
        for i, (rl, rv) in enumerate(right):
            lft = left[i] if i < len(left) else ""
            out.write(f"{lft:<34}{rl:<16}= ", "text")
            out.write(f"{rv:>10}\n", "result")

    @property
    def population(self) -> float:
        return float(self.dd.w[self.in_sample].sum())

    def post(self, s, est, req, fit, df):
        """e() of svy in Stata 14, in ereturn list order (compat 0601)."""
        from ..commands.matrix import Matrix
        from .results import post
        d, dd, V, Z, w = self.design, self.dd, est.V, fit.Z, self.w
        stat, over, cn = req.stat, req.over, fit.colnames
        est.svy_effects = dict(V=V, srs=self.srs, labels=fit.labels, vals=fit.vals, title=req.info.column,
                               vcetype=dd.vcetype, over=over, level=req.level, legend=fit.legend)
        k = len(fit.vals)
        scalars = [("df_r", float(df) if df is not None else SYS)]
        if self.subpop is not None:
            scalars += [("N_subpop", float(w.sum())), ("N_sub", float(self.n_obs))]
        scalars += [("N_strata_omit", 0.0), ("singleton", 0.0), ("census", 0.0),
                    ("N_pop", self.population), ("N_psu", float(dd.n_psu)),
                    ("N_strata", float(dd.n_strata)), ("N_over", float(len(fit.groups) if over else 1)),
                    ("N", float(self.n_sample)), ("stages", 1.0), ("k_eq", 1.0),
                    ("rank", float(np.linalg.matrix_rank(V)) if V.size else 0.0)]
        cmdtxt = f"{stat} {req.args.strip()}".strip()
        opts = self.svyopts.strip()
        macros = [("cmdline", f"svy{', ' + opts if opts else ' '}: {cmdtxt}"),
                  ("cmd", stat), ("prefix", "svy"), ("cmdname", stat),
                  ("command", cmdtxt if "," in cmdtxt else cmdtxt + ","),
                  ("title", f"Survey: {req.info.title}"), ("vcetype", dd.vcetype), ("vce", dd.vce),
                  ("estat_cmd", "svy_estat"), ("varlist", " ".join(req.allv))]
        if over:
            glabs = " ".join(g[0] for g in fit.groups)
            macros += [("over", over), ("over_labels", glabs), ("over_namelist", glabs)]   # VERIFY
        if self.subpop is not None:
            macros.append(("subpop", str(self.so["subpop"]).strip()))
        macros += [("marginsnotok", "_ALL")]
        if d.weight:
            macros += [("wtype", d.wtype or "pweight"), ("wvar", d.weight), ("wexp", f"= {d.weight}")]
        macros.append(("singleunit", d.singleunit))
        if d.psu and d.psu != "_n":
            macros.append(("su1", d.psu))
        if d.strata:
            macros.append(("strata1", d.strata))
        if d.fpc:
            macros.append(("fpc1", d.fpc))
        macros += [("properties", "b V"), ("depvar", req.info.column)]
        Nk = np.array([[float(np.count_nonzero(Z[:, j]) or self.n_obs) for j in range(k)]])
        Wk = np.array([[float(np.sum(w * (np.abs(Z[:, j]) > 0))) or float(w.sum()) for j in range(k)]])
        nn_ = max(self.n_obs, 2)
        Vsrs = nn_ / (nn_ - 1) * Z.T @ Z
        mats = [("_N_subp", Matrix(Wk, ["r1"], list(cn))), ("V_srs", Matrix(Vsrs, list(cn), list(cn))),
                ("_N", Matrix(Nk, ["r1"], list(cn))),
                ("_N_strata_certain", Matrix(np.zeros((1, 1)))), ("_N_strata_single", Matrix(np.zeros((1, 1)))),
                ("_N_strata", Matrix(np.array([[float(dd.n_strata)]]))),
                ("error", Matrix(np.zeros((1, k)), ["r1"], list(cn)))]
        post(s, est, scalars, macros, dd.mask.copy(), extra_mats=mats)


def format_size(x: float) -> str:
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


def srs_variances(Z: np.ndarray, w: np.ndarray, dd: "DesignData", m: int) -> tuple:
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


def subpop_mask(s: "Session", text: str) -> np.ndarray:
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
    from ..commands._util import touse
    d = require_design(s)
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
        from .summaries import summary_stat
        summary_stat(s, name, rest, engine=lambda s_, req: SurveyEngine(s_, req, opts))
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
    if rest.strip(" ,"):
        # validated only: the table always shows DEFF and DEFT
        match_options(rest.lstrip(", "), {"deff": 4, "deft": 4, "srssubpop": 3, "meff": 4, "meft": 4})
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
