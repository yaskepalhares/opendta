"""jackknife, permute and simulate prefixes ([R] jackknife, [R] permute, [R] simulate).

They share the expression list, the evaluation of statistics and the replication loop with
bootstrap.py. The layouts follow the manuals; what the manuals do not show has the VERIFY mark
(compat/do/1903)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail as _chi2_tail
from ..lang.syntax import find_top, match_options
from ..lazy import lazy
from .bootstrap import (_evaluate, _hdr, _quiet, _subset, _warn_rclass, parse_exp_list)
from .results import CoefRow, Estimates, abbrev, coef_table, g9, post

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_COMMON = {"saving": 3, "nodots": 5, "dots": 4, "noisily": 3, "trace": 5, "nolegend": 6,
           "verbose": 1, "double": 6, "level": 1, "notable": 5, "noheader": 6, "title": 5,
           "seed": 4, "reps": 1, "force": 5}
_JK = dict(_COMMON, eclass=1, rclass=1, n=1, keep=1, mse=1, cluster=2, idcluster=3, strata=3, dof=3)
_PM = dict(_COMMON, left=2, right=2, strata=3)
_SIM = dict(_COMMON, nolegend=6)


def _split(args: str) -> tuple[str, str, str]:
    """`exp_list [, options] : command` -> (exp_list, options, command)."""
    colon = find_top(args, ":")
    if colon == -1:
        raise StataError(198, "invalid syntax")
    head, command = args[:colon].strip(), args[colon + 1:].strip()
    if not command:
        raise StataError(198, "invalid syntax")
    c = find_top(head, ",")
    return (head, "", command) if c == -1 else (head[:c], head[c + 1:], command)


def _named(items: list[tuple[str, str]], stem: str) -> list[tuple[str, str]]:
    """Default names _jk_1, _jk_2... for the unnamed expressions (not for `_b` and `_se`)."""
    out, k = [], 0
    for n, e in items:
        if not n and e not in ("_b", "_se"):
            k += 1
            n = f"{stem}{k}"
        out.append((n, e))
    return out


def _dots(s: "Session", title: str, n: int, o: dict) -> bool:
    on = not o.get("nodots")
    if on:
        s.output.write(f"\n{title} ({n})\n", "text")
        s.output.write("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 \n", "text")
    return on


def _dot(s: "Session", ok: bool, r: int, n: int) -> None:
    s.output.write("." if ok else "x", "text")
    if (r + 1) % 50 == 0:
        s.output.write(f"{r + 1:>6}\n", "text")
    elif r + 1 == n:
        s.output.write("\n", "text")


def _save(s: "Session", spec: str, names: list[str], R: np.ndarray, double: bool) -> None:
    """saving(filename[, replace]): the replicates as a dataset."""
    file, _, opt = spec.partition(",")
    base = s.data
    try:
        s.data = _dataset(names, R, double)
        s.interp._execute_expanded(f'save "{file.strip()}"' + (", replace" if "replace" in opt else ""))
    finally:
        s.data = base


def _dataset(names: list[str], R: np.ndarray, double: bool):
    from ..core.dataset import Dataset, Variable
    ds = Dataset()
    ds.set_obs(R.shape[0])
    for j, nm in enumerate(names):
        col = np.where(np.isnan(R[:, j]), SYS, R[:, j])
        ds.add(Variable(nm, "double" if double else "float", col))
    return ds


def _run0(s: "Session", command: str, noisy: bool):
    """Runs the command on the data in memory; returns (cmdname, is e-class, e0, est0)."""
    from ..commands.registry import lookup
    spec = lookup(command.split(None, 1)[0])
    cmdname = spec.name if spec else command.split(None, 1)[0]
    _quiet(s, command, noisy)
    est0 = getattr(s, "_estimates", None)
    eclass = "b" in s.e and est0 is not None and getattr(est0, "cmd", "") == cmdname
    return cmdname, eclass, dict(s.e), est0


def _sample(s: "Session", eclass: bool):
    base = s.data
    if eclass and base.esample is not None and len(base.esample) == base.nobs:
        return base.esample.copy()
    return np.ones(base.nobs, dtype=bool)


def _legend(s: "Session", command: str, names: list[str], exps: list[str], extra=()) -> None:
    out = s.output
    out.write(f"{'command':>13}:  {command}\n", "text")
    for n, ex in zip(names, exps):
        out.write(f"{abbrev(n, 12):>13}:  {ex}\n", "text")
    for lab, val in extra:
        out.write(f"{lab:>13}:  {val}\n", "text")


# ---------------------------------------------------------------------------
# jackknife
# ---------------------------------------------------------------------------

def cmd_jackknife(s: "Session", args: str) -> None:
    exp_text, opt_text, command = _split(args)
    o = match_options(opt_text, _JK) if opt_text.strip() else {}
    items = _named(parse_exp_list(exp_text), "_jk_")
    out = s.output
    level = float(o["level"]) if isinstance(o.get("level"), str) else float(s.settings.get("level", 95))
    noisy = bool(o.get("noisily"))
    cmdname, eclass, e0, est0 = _run0(s, command, noisy)
    if o.get("rclass"):
        eclass = False
    if not items:
        if not eclass:
            raise StataError(198, "expression list required")   # VERIFY
        items = [("", "_b")]
    eclass_b = eclass and items == [("", "_b")]
    names, observed = _evaluate(s, items)
    base = s.data
    sample = _sample(s, eclass)
    if eclass:
        out.write(f"(running {cmdname} on estimation sample)\n", "text")
    else:
        _warn_rclass(s, cmdname, "jackknife")
    r0 = dict(s.r)
    work = _subset(base, np.flatnonzero(sample))
    nw = work.nobs
    if isinstance(o.get("cluster"), str):
        cv = work.get(s.expand_varlist(o["cluster"])[0])
        cl = np.asarray(cv.raw if cv.is_string else cv.data)
        units = [np.flatnonzero(cl == v) for v in np.unique(cl)]
    else:
        units = [np.array([i]) for i in range(nw)]
    nu = len(units)
    keep_all = np.ones(nw, dtype=bool)
    R = np.full((nu, len(names)), np.nan)
    dots = _dots(s, "Jackknife replications", nu, o)
    mis = 0
    try:
        for r, u in enumerate(units):
            keep_all[u] = False
            s.data = _subset(work, np.flatnonzero(keep_all))
            keep_all[u] = True
            ok = True
            try:
                _quiet(s, command, noisy)
                _, vals = _evaluate(s, items)
                if len(vals) != len(names):
                    ok = False
                else:
                    R[r] = vals if o.get("double") else np.asarray(vals, dtype=np.float32)
            except StataError:
                ok = False
            mis += 0 if ok else 1
            if dots:
                _dot(s, ok, r, nu)
    finally:
        s.data = base
    s.r = r0
    good = R[~np.isnan(R).any(axis=1)]
    n = len(good)
    if n < 2:
        raise StataError(2000, "insufficient observations to compute jackknife standard errors")   # VERIFY
    centre = observed if o.get("mse") else good.mean(axis=0)
    D = good - centre
    V = (n - 1) / n * (D.T @ D)
    if isinstance(o.get("saving"), str):
        _save(s, o["saving"], names, R, bool(o.get("double")))
    df = float(o["dof"]) if isinstance(o.get("dof"), str) else None
    _post_jk(s, names, observed, V, command, items, eclass_b, est0, e0, cmdname, sample, n, mis, level, df, o)
    _display_jk(s, o)


def _post_jk(s, names, observed, V, command, items, eclass_b, est0, e0, cmdname, sample, n, mis, level,
             df, o) -> None:
    from ..commands.matrix import Matrix
    cols = [x.split(":")[-1] for x in names]
    eqs = [x.split(":")[0] if ":" in x else "_" for x in names]
    se = np.sqrt(np.diag(V))
    if eclass_b and est0 is not None:
        est = est0
        stat = "t" if est.stat == "t" else "z"
        if df is None and stat == "t":
            df = float(n - 1)          # VERIFY: degrees of freedom of the t statistics (clusters - 1)
        if df is not None:
            stat = "t"
        est.V = V
        est.vcetype = "Jackknife"
        est.title = getattr(est0, "title", "") or e0.get("title", "")
        est.stat = stat
        est.df_r = df
        scalars = [(k, v) for k, v in e0.items() if isinstance(v, float) and k not in ("N_reps", "N_misreps")]
        scalars += [("N_reps", float(n)), ("N_misreps", float(mis)), ("level", level)]
        if df is not None:
            scalars = [(k, v) for k, v in scalars if k != "df_r"] + [("df_r", df)]
        macros = [(k, v) for k, v in e0.items() if isinstance(v, str)
                  and k not in ("vce", "vcetype", "prefix", "command")]
        macros += [("prefix", "jackknife"), ("vce", "jackknife"), ("vcetype", "Jackknife"),
                   ("command", command)]
        keep = [j for j, x in enumerate(names) if not x.endswith("_cons") and V[j, j] > 0]
        est.jk = {"names": names, "observed": observed, "chi2": None, "df_m": len(keep)}
        if keep:
            try:
                est.jk["chi2"] = float(observed[keep] @ np.linalg.solve(V[np.ix_(keep, keep)], observed[keep]))
            except np.linalg.LinAlgError:
                pass
        est.jk.update(n=n, level=level, eclass_b=True, command=command, names_exp=[e for _, e in items])
        post(s, est, list(reversed(scalars)), list(reversed(macros)), sample)
        return
    est = Estimates(cmd="jackknife", depvar="", names=cols, b=np.asarray(observed, dtype=float), V=V,
                    N=int(sample.sum()), stat="t" if df else "z")
    est.eqnames = eqs
    est.df_r = df
    est.rows = [CoefRow("coef", label=c, b=float(b), se=float(v), name=nm)
                for nm, c, b, v in zip(names, cols, observed, se)]
    est.vcetype = "Jackknife"
    est.title = "Jackknife results"
    est.level = level
    est.jk = {"names": names, "observed": observed, "chi2": None, "df_m": 0, "n": n, "level": level,
              "eclass_b": False, "command": command, "names_exp": [e for _, e in items]}
    scalars = [("N", float(int(sample.sum()))), ("N_reps", float(n)), ("N_misreps", float(mis)),
               ("k_eq", 1.0), ("k_exp", float(len(names))), ("level", level)]
    if df is not None:
        scalars.insert(1, ("df_r", df))
    macros = [("cmd", "jackknife"), ("command", command), ("prefix", "jackknife"), ("vce", "jackknife"),
              ("vcetype", "Jackknife"), ("title", "Jackknife results"), ("properties", "b V")]
    for j, (_, ex) in enumerate(items, start=1):
        macros.append((f"exp{j}", ex))
    post(s, est, scalars, macros, sample)


def _display_jk(s: "Session", o: dict) -> None:
    est = s._estimates
    x = est.jk
    out = s.output
    if not o.get("noheader"):
        out.write("\n", "text")
        items = [("Number of obs", f"{est.N:,}"), ("Replications", f"{x['n']:,}")]
        if x["eclass_b"]:
            k = x["df_m"]
            if x["chi2"] is not None:
                if est.df_r is not None:
                    items.append((f"F({k}, {int(est.df_r)})", f"{x['chi2'] / k:.2f}"))
                    items.append(("Prob > F", f"{float(st.f.sf(x['chi2'] / k, k, est.df_r)):.4f}"))
                else:
                    items.append((f"Wald chi2({k})", f"{x['chi2']:.2f}"))
                    items.append(("Prob > chi2", f"{_chi2_tail(x['chi2'], k):.4f}"))
            for key, lab, fmt in (("r2", "R-squared", "%9.4f"), ("r2_a", "Adj R-squared", "%9.4f"),
                                  ("r2_p", "Pseudo R2", "%9.4f"), ("rmse", "Root MSE", "%9.4f")):
                v = s.e.get(key)
                if isinstance(v, float) and v < SYS:
                    items.append((lab, g9(v, fmt)))
        _hdr(s, est.title if x["eclass_b"] else "Jackknife results", items)
        if not x["eclass_b"]:
            out.write("\n", "text")
            _legend(s, x["command"], [n for n in x["names"]], x["names_exp"])
        out.write("\n", "text")
    if o.get("notable"):
        return
    se = np.sqrt(np.diag(est.V))
    if x["eclass_b"]:
        by = dict(zip(x["names"], se))
        short = {n.split(":")[-1]: v for n, v in by.items()}
        rows = []
        for r in est.rows:
            if r.kind in ("coef", "aux"):
                key = r.name or r.label
                rows.append(CoefRow(r.kind, r.label, r.b, float(by.get(key, short.get(key.split(":")[-1], r.se))),
                                    r.eq, r.name))
            else:
                rows.append(r)
        dep = est.depvar
    else:
        rows = [CoefRow("coef", label=n.split(":")[-1], b=float(b), se=float(v), name=n)
                for n, b, v in zip(x["names"], x["observed"], se)]
        dep = ""
    coef_table(s, dep, rows, stat=est.stat, df=est.df_r, level=x["level"],
               top_line=f"{'':>12} |" + " " * 14 + "Jackknife")


# ---------------------------------------------------------------------------
# permute
# ---------------------------------------------------------------------------

def cmd_permute(s: "Session", args: str) -> None:
    head, opt_text, command = _split(args)
    o = match_options(opt_text, _PM) if opt_text.strip() else {}
    parts = head.split(None, 1)
    if not parts:
        raise StataError(198, "permutation variable required")
    pvar = s.expand_varlist(parts[0])[0]
    items = _named(parse_exp_list(parts[1] if len(parts) > 1 else ""), "_pm_")
    if not items:
        raise StataError(198, "expression list required")   # VERIFY
    out = s.output
    reps = int(o["reps"]) if isinstance(o.get("reps"), str) else 100
    if reps < 1:
        raise StataError(198, "reps() must be a positive integer")
    level = float(o["level"]) if isinstance(o.get("level"), str) else float(s.settings.get("level", 95))
    if isinstance(o.get("seed"), str):
        s.interp._execute_expanded(f"set seed {o['seed'].strip()}")
    noisy = bool(o.get("noisily"))
    _run0(s, command, noisy)
    names, observed = _evaluate(s, items)
    base = s.data
    r0 = dict(s.r)
    n = base.nobs
    var = base.get(pvar)
    original = var.raw.copy()
    strata = None
    if isinstance(o.get("strata"), str):
        sv = base.get(s.expand_varlist(o["strata"])[0])
        strata = np.asarray(sv.raw if sv.is_string else sv.data)
    groups = [np.arange(n)] if strata is None else [np.flatnonzero(strata == v) for v in np.unique(strata)]
    from ..core import rng
    R = np.full((reps, len(names)), np.nan)
    dots = _dots(s, "Permutation replications", reps, o)
    mis = 0
    try:
        for r in range(reps):
            perm = np.arange(n)
            for g in groups:
                key = rng.RNG.uniform(len(g))
                perm[g] = g[np.argsort(key, kind="stable")]
            var.raw = original[perm]
            ok = True
            try:
                _quiet(s, command, noisy)
                _, vals = _evaluate(s, [(a, b) for a, b in zip(names, [e for _, e in items])])
                R[r] = vals if o.get("double") else np.asarray(vals, dtype=np.float32)
            except StataError:
                ok = False
            mis += 0 if ok else 1
            if dots:
                _dot(s, ok, r, reps)
    finally:
        var.raw = original
    s.r = r0
    if isinstance(o.get("saving"), str):
        _save(s, o["saving"], names, R, bool(o.get("double")))
    side = "left" if o.get("left") else "right" if o.get("right") else "two"
    out.write("\n", "text")
    _hdr(s, "Monte Carlo permutation results", [("Number of obs", f"{n:,}")])
    out.write("\n", "text")
    _legend(s, command, names, [e for _, e in items], [("permute var", pvar)])
    out.write("\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(f"{'T':<12} |     T(obs)       c       n   p=c/n   SE(p) [{level:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    a = 1 - level / 100
    cs, ns = [], []
    for j, nm in enumerate(names):
        col = R[:, j][~np.isnan(R[:, j])]
        obs = observed[j]
        # the observed statistic is stored as float, as the replicates are
        t = float(np.float32(obs)) if not o.get("double") else obs
        c = int((col <= t).sum()) if side == "left" else int((col >= t).sum()) if side == "right" \
            else int((np.abs(col) >= abs(t)).sum())
        m = len(col)
        p = c / m if m else float("nan")
        se = math.sqrt(p * (1 - p) / m) if m else float("nan")
        lo = 0.0 if c == 0 else float(st.beta.ppf(a / 2, c, m - c + 1))
        hi = 1.0 if c == m else float(st.beta.ppf(1 - a / 2, c + 1, m - c))
        out.write(f"{abbrev(nm, 12):>12} |", "text")
        out.write(f"{g9(obs, '%10.0g'):>11}{c:>8}{m:>8}{p:>8.4f}{se:>8.4f}{lo:>10.4f}{hi:>10.4f}\n", "result")
        cs.append(c)
        ns.append(m)
    out.write("-" * 78 + "\n", "text")
    out.write(" Note:  confidence interval is with respect to p=c/n.\n", "text")
    note = {"left": "c = #{T <= T(obs)}", "right": "c = #{T >= T(obs)}", "two": "c = #{|T| >= |T(obs)|}"}[side]
    out.write(f" Note:  {note}\n", "text")      # VERIFY: layout of the notes
    from ..commands.matrix import Matrix
    s.r = {"N_reps": float(reps)}
    s.r["N"] = float(n)
    s.r["c"] = Matrix(np.array(cs, dtype=float).reshape(1, -1), ["c"], list(names))
    s.r["p"] = Matrix(np.array([c / m for c, m in zip(cs, ns)]).reshape(1, -1), ["p"], list(names))


# ---------------------------------------------------------------------------
# simulate
# ---------------------------------------------------------------------------

def cmd_simulate(s: "Session", args: str) -> None:
    exp_text, opt_text, command = _split(args)
    o = match_options(opt_text, _SIM) if opt_text.strip() else {}
    if not isinstance(o.get("reps"), str):
        raise StataError(198, "option reps() required")
    reps = int(o["reps"])
    if reps < 1:
        raise StataError(198, "reps() must be a positive integer")
    items = _named(parse_exp_list(exp_text), "_sim_")
    if not items:
        raise StataError(198, "expression list required")   # VERIFY
    names = [n for n, _ in items]
    exps = [e for _, e in items]
    out = s.output
    if isinstance(o.get("seed"), str):
        s.interp._execute_expanded(f"set seed {o['seed'].strip()}")
    noisy = bool(o.get("noisily"))
    pairs = list(zip(names, exps))
    out.write("\n", "text")
    _legend(s, command, names, exps)
    R = np.full((reps, len(names)), np.nan)
    dots = _dots(s, "Simulations", reps, o)
    r0 = dict(s.r)
    for r in range(reps):
        ok = True
        try:
            _quiet(s, command, noisy)
            _, vals = _evaluate(s, pairs)
            R[r] = vals
        except StataError:
            ok = False
        if dots:
            _dot(s, ok, r, reps)
    s.r = r0
    double = bool(o.get("double"))
    ds = _dataset(names, R, double)
    if isinstance(o.get("saving"), str):
        _save(s, o["saving"], names, R, double)
    s.data = ds
