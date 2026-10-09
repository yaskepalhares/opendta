"""Amostras complexas ([SVY] svyset, svydescribe, svy estimation) e as
estatísticas mean, proportion, total e ratio ([R] mean, [R] proportion...).

Variância por linearização de Taylor (desenho com estratos, UPAs e fpc
num estágio):

    V = Σ_h (1 - f_h) n_h/(n_h - 1) Σ_j (z_hj - z̄_h)(z_hj - z̄_h)'

z_hj é a soma, na UPA j do estrato h, dos escores linearizados de cada
observação. Graus de liberdade do desenho = nº de UPAs - nº de estratos.
A declaração fica nas características de _dta (_svy_*), como no Stata.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from .results import g9

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# svyset
# ---------------------------------------------------------------------------

@dataclass
class Design:
    psu: str = ""            # "" ou "_n" = cada observação é uma UPA
    weight: str = ""
    wtype: str = ""          # pweight / iweight
    strata: str = ""
    fpc: str = ""
    singleunit: str = "missing"
    vce: str = "linearized"
    poststrata: str = ""
    postweight: str = ""

    @property
    def is_set(self) -> bool:
        return True


def get_design(ds) -> Design | None:
    ch = ds.chars.get("_dta", {})
    if ch.get("_svy_version", "") == "":
        return None
    return Design(psu=ch.get("_svy_su1", ""), weight=ch.get("_svy_wvar", ""), wtype=ch.get("_svy_wtype", ""),
                  strata=ch.get("_svy_strata1", ""), fpc=ch.get("_svy_fpc1", ""),
                  singleunit=ch.get("_svy_singleunit", "missing") or "missing",
                  vce=ch.get("_svy_vce", "linearized") or "linearized",
                  poststrata=ch.get("_svy_poststrata", ""), postweight=ch.get("_svy_postweight", ""))


def _show_design(s: "Session", d: Design) -> None:
    out = s.output
    rows = []
    if d.weight:
        rows.append((d.wtype, d.weight))
    else:
        rows.append(("pweight", "<none>"))
    rows += [("VCE", d.vce), ("Single unit", d.singleunit), ("Strata 1", d.strata or "<one>"),
             ("SU 1", d.psu or "<observations>"), ("FPC 1", d.fpc or "<zero>")]
    if d.poststrata:
        rows[1:1] = [("Poststrata", d.poststrata), ("Postweight", d.postweight)]
    out.write("\n", "text")
    for lab, val in rows:
        out.write(f"{lab:>13}: ", "text")
        out.write(val + "\n", "result")


def cmd_svyset(s: "Session", args: str) -> None:
    ds = s.data
    t = args.strip()
    if not t:
        d = get_design(ds)
        if d is None:
            s.output.write("no survey characteristics are set\n", "text")   # VERIFICAR
            return
        _show_design(s, d)
        return
    if t.startswith(",") and re.search(r"\bclear\b", t):
        for k in list(ds.chars.get("_dta", {})):
            if k.startswith("_svy_"):
                del ds.chars["_dta"][k]
        return
    p = parse_standard(t)
    o = match_options(p.options, {"strata": 3, "fpc": 3, "singleunit": 6, "vce": 3, "poststrata": 5,
                                  "postweight": 5, "clear": 5}) if p.options.strip() else {}
    from ..core.varlist import resolve_name
    words = p.varlist.split()
    psu = ""
    if words:
        psu = "_n" if words[0] == "_n" else resolve_name(ds, words[0])
    d = Design(psu=psu if psu != "_n" else "")
    if p.weight:
        wtype, wexp = p.weight
        if wtype not in ("pweight", "iweight", "weight"):
            raise StataError(101, f"{wtype}s not allowed")   # VERIFICAR
        d.wtype = "pweight" if wtype == "weight" else wtype
        d.weight = resolve_name(ds, wexp.strip())
    if o.get("strata"):
        d.strata = resolve_name(ds, str(o["strata"]).strip())
    if o.get("fpc"):
        d.fpc = resolve_name(ds, str(o["fpc"]).strip())
    if o.get("singleunit"):
        su = str(o["singleunit"]).strip()
        if su not in ("missing", "certainty", "scaled", "centered"):
            raise StataError(198, "invalid singleunit() option")   # VERIFICAR
        d.singleunit = su
    if o.get("vce"):
        d.vce = str(o["vce"]).strip()
    if o.get("poststrata"):
        d.poststrata = resolve_name(ds, str(o["poststrata"]).strip())
        d.postweight = resolve_name(ds, str(o.get("postweight", "")).strip())
    ch = ds.chars.setdefault("_dta", {})
    for k in list(ch):
        if k.startswith("_svy_"):
            del ch[k]
    ch.update({"_svy_version": "2", "_svy_su1": d.psu, "_svy_wvar": d.weight, "_svy_wtype": d.wtype,
               "_svy_strata1": d.strata, "_svy_fpc1": d.fpc, "_svy_singleunit": d.singleunit,
               "_svy_vce": d.vce, "_svy_poststrata": d.poststrata, "_svy_postweight": d.postweight})
    _show_design(s, d)
    s.r = {"settings": t, "singleunit": d.singleunit, "vce": d.vce}
    if d.strata:
        s.r["strata1"] = d.strata
    s.r["su1"] = d.psu or "_n"
    if d.weight:
        s.r["wvar"], s.r["wtype"] = d.weight, d.wtype
    s.r = dict(reversed(list(s.r.items())))


def require_design(s: "Session") -> Design:
    d = get_design(s.data)
    if d is None:
        raise StataError(119, "data not set up for svy, use svyset")
    return d


# ---------------------------------------------------------------------------
# estrutura do desenho na amostra
# ---------------------------------------------------------------------------

class DesignData:
    """Pesos, estratos e UPAs das observações em `mask` (amostra do desenho)."""

    def __init__(self, s: "Session", d: Design, mask: np.ndarray):
        ds = s.data
        self.d = d
        n_all = ds.nobs
        w = ds.get(d.weight).data.astype(float) if d.weight else np.ones(n_all)
        mask = mask & (w < SYS)
        strata = ds.get(d.strata).data if d.strata else np.zeros(n_all)
        psu = ds.get(d.psu).data if d.psu else np.arange(n_all, dtype=float)
        mask = mask & (strata < SYS) & (psu < SYS)
        self.mask = mask
        self.w = w[mask]
        self.strata = strata[mask]
        self.psu = psu[mask]
        self.fpc = ds.get(d.fpc).data[mask].astype(float) if d.fpc else None
        self.n = int(mask.sum())
        # identificadores de UPA dentro de estrato
        self.h_ids, self.h_inv = np.unique(self.strata, return_inverse=True)
        pairs = np.column_stack([self.h_inv, self.psu])
        _, self.u_inv = np.unique(pairs, axis=0, return_inverse=True)
        self.u_inv = self.u_inv.ravel()
        self.n_strata = len(self.h_ids)
        self.n_psu = int(self.u_inv.max()) + 1 if self.n else 0
        self.df = self.n_psu - self.n_strata

    def variance(self, z: np.ndarray) -> np.ndarray:
        """V a partir dos escores z (n x k) das observações da amostra."""
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
                    single = True   # VERIFICAR: escalonamento pela média dos demais estratos
                continue
            dev = U[idx] - U[idx].mean(axis=0)
            V += (1 - f) * nh / (nh - 1) * dev.T @ dev
        if single:
            return np.full((k, k), np.nan)
        return V


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
    """Valores e escores (n x k) de cada estatística, por grupo de over()."""
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
            else:   # proportion: y já é indicador
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
    if vcetype:
        start = 28 + (9 - len(vcetype)) // 2
        out.write(f"{'':>12} |" + " " * (start - 14) + vcetype + "\n", "text")
    head = "Over" if over_names else ""
    out.write(f"{head:>12} | {stat_title:>10}   Std. Err.     [{lv}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 48 + "\n", "text")
    last = None
    for (var, grp), b, e in zip(labels, vals, se):
        if grp:
            if var != last:
                # VERIFICAR: alinhamento do nome da variável (à esquerda com over, à direita sem)
                out.write((f"{var[:12]:<12} |" if over_names else f"{var[:12]:>12} |") + "\n", "text")
                last = var
            lab = grp
        else:
            lab = var
        out.write(f"{lab[:12]:>12} |", "text")
        if e is None or not np.isfinite(e):
            out.write(f"  {g9(b):>9}  {'.':>9}     {'.':>9}   {'.':>9}\n", "result")
            continue
        if stat_title == "Proportion" and 0 < b < 1:
            # IC na escala logit (padrão do proportion no Stata 14, VERIFICAR)
            from scipy.special import expit, logit
            half = crit * e / (b * (1 - b))
            lo, hi = expit(logit(b) - half), expit(logit(b) + half)
        else:
            lo, hi = b - crit * e, b + crit * e
        out.write(f"  {g9(b):>9}  {g9(e):>9}     {g9(lo):>9}   {g9(hi):>9}\n", "result")
    out.write("-" * width + "\n", "text")


def _proportion_vars(s, names, mask):
    """proportion: um indicador por categoria de cada variável."""
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
        names = [f"{a}/{b}" for _, a, b in pairs]
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
        so = match_options(svyopts, {"subpop": 3, "vce": 3}) if svyopts.strip() else {}
        dd = DesignData(s, d, np.ones(ds.nobs, bool) if not p.if_ and not p.in_ else touse(s, p))
        if so.get("subpop"):
            subpop = _subpop_mask(s, str(so["subpop"]))
        est_mask = mask & dd.mask
        out.write(f"\n(running {stat} on estimation sample)\n", "text")
        # escores fora da amostra (ou da subpop) são zero, mas o desenho conta todas as UPAs
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
    if stat == "proportion":
        vals, Z, labels = [], [], []
        for i, y in enumerate(Y):
            for gname, sel, _ in groups:
                v_, Z_, _l = _estimate("proportion", [y], w, [(gname, sel, None)], [names_eff[i]])
                vals.extend(v_.tolist())
                Z.append(Z_)
                labels.append((labs_p[i][0], f"{labs_p[i][1]}" + (f" {gname}" if over else "")))
        vals = np.array(vals)
        Z = np.column_stack(Z) if Z else np.zeros((len(w), 0))
    else:
        vals, Z, labels = _estimate(stat, Y, w, groups, names, X)
    if svy:
        V = dd.variance(Z)
        df = dd.df
        vcetype = "Linearized"
    else:
        # V = N/(N-1) Σ z z' (com fweight, cada escore conta w vezes)
        robust = str(o.get("vce", "")).startswith("r")
        vcetype = "Robust" if robust else ""
        if p.weight and p.weight[0] == "fweight":
            nn = float(w.sum())
            V = nn / (nn - 1) * (Z / np.where(w > 0, w, 1)[:, None]).T @ Z
        else:
            nn = float(len(w))
            V = nn / (nn - 1) * Z.T @ Z
        df = int(nn) - 1
    se = np.sqrt(np.maximum(np.diag(V), 0)) if V.size else np.zeros(0)
    if np.isnan(V).any():
        se = np.full(len(vals), np.nan)
    # cabeçalho
    if not o.get("noheader"):
        if svy:
            pop = float(dd.w[est_mask[full]].sum())
            out.write("\n", "text")
            out.write(f"Survey: {_TITLES[stat]}\n\n", "text")
            right = [("Number of obs", f"{int(est_mask[full].sum()):,}"), ("Population size", _num(pop))]
            if subpop is not None:
                right += [("Subpop. no. obs", f"{n_obs:,}"), ("Subpop. size", _num(float(w.sum())))]
            right.append(("Design df", f"{df:,}"))
            left = [f"Number of strata = {dd.n_strata:>7,}", f"Number of PSUs   = {dd.n_psu:>7,}"]
            for i, (rl, rv) in enumerate(right):
                lft = left[i] if i < len(left) else ""
                out.write(f"{lft:<33}{rl:<18}= ", "text")
                out.write(f"{rv:>11}\n", "result")
        else:
            out.write("\n", "text")
            out.write(f"{_TITLES[stat]:<34}{'Number of obs':<16}= ", "text")
            out.write(f"{n_obs:>10,}\n", "result")
        if over:
            out.write("\n" + f"Over: {over}\n", "text")   # VERIFICAR legenda dos grupos
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
    scalars = [("N", float(n_obs))]
    if svy:
        scalars += [("N_strata", float(dd.n_strata)), ("N_psu", float(dd.n_psu)), ("N_pop", float(dd.w[est_mask[full]].sum())),
                    ("df_r", float(df))]
        if subpop is not None:
            scalars += [("N_sub", float(n_obs)), ("N_subpop", float(w.sum()))]
    else:
        scalars += [("df_r", float(df))]
    macros = [("cmd", stat), ("cmdline", f"{stat} {args.strip()}".strip()), ("properties", "b V"),
              ("varlist", " ".join(allv))]
    if svy:
        macros = [("prefix", "svy"), ("cmdname", stat)] + macros
        macros.append(("vcetype", "Linearized"))
    post(s, est, scalars, macros, np.isin(np.arange(ds.nobs), rows_idx))


def _num(x: float) -> str:
    return f"{x:,.0f}" if abs(x - round(x)) < 1e-6 else format_value(x, "%11.0g", pad=False).strip()


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
    out.write("\n" + " " * 36 + "#Obs per Unit\n", "text")
    out.write(" " * 30 + "-" * 28 + "\n", "text")
    out.write("Stratum   #Units     #Obs       min      mean       max\n", "text")
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
