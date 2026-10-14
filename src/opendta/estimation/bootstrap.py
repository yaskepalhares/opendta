"""bootstrap prefix, vce(bootstrap) and estat bootstrap ([R] bootstrap,
[R] bootstrap postestimation).

The draws use OpenDTA's generator (mt64) resampling observations,
or clusters within strata. The sequence of draws of Stata's bsample is
undocumented: with the same seed, the replicates (and therefore the
standard errors) differ from Stata's, although they are statistically
equivalent (VERIFY with compat/do/0803)."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core.stats import chi2_tail as _chi2_tail

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from .results import CoefRow, Estimates, LINE, abbrev, coef_table, g9

if TYPE_CHECKING:
    from ..core.dataset import Dataset
    from ..session import Session

SYS = M.SYSMISS

_OPTS = {"reps": 1, "seed": 4, "strata": 3, "cluster": 2, "idcluster": 3, "size": 4, "nodots": 5,
         "dots": 4, "saving": 3, "level": 1, "notable": 5, "noheader": 6, "nolegend": 6, "force": 5,
         "reject": 3, "bca": 3, "mse": 3, "noisily": 1, "trace": 5, "title": 5, "nowarn": 6,
         "double": 6}


@dataclass
class BootInfo:
    reps: np.ndarray                     # replicates × statistics (NaN = lost replicate)
    names: list[str]
    observed: np.ndarray
    command: str
    exps: list[str]
    n_reps: int
    n_mis: int
    level: float
    eclass_b: bool
    title: str = "Bootstrap results"
    N: int = 0
    strata: str = ""
    cluster: str = ""
    extra: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# expression list
# ---------------------------------------------------------------------------

def _tokens(text: str) -> list[str]:
    out, i, n = [], 0, len(text)
    while i < n:
        if text[i].isspace():
            i += 1
            continue
        if text[i] == "(":
            depth, j = 0, i
            while j < n:
                if text[j] == "(":
                    depth += 1
                elif text[j] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
            continue
        j = i
        depth = 0
        while j < n and (depth or not text[j].isspace()):
            if text[j] in "([":
                depth += 1
            elif text[j] in ")]":
                depth -= 1
            j += 1
        out.append(text[i:j])
        i = j
    return out


def parse_exp_list(text: str) -> list[tuple[str, str]]:
    """[(name or "", expression)]; `_b` and `_se` stay as they are."""
    items = []
    for tok in _tokens(text):
        t = tok
        if t.startswith("(") and t.endswith(")"):
            t = t[1:-1].strip()
            m = re.match(r"^([A-Za-z_]\w*)\s*:\s*(.+)$", t)
            if m:
                items.append((m.group(1), m.group(2).strip()))
                continue
        m = re.match(r"^([A-Za-z_]\w*)\s*=(?!=)\s*(.+)$", t)
        if m:
            items.append((m.group(1), m.group(2).strip()))
        else:
            items.append(("", t))
    return items


def _evaluate(s: "Session", items: list[tuple[str, str]]) -> tuple[list[str], np.ndarray]:
    names: list[str] = []
    vals: list[float] = []
    k = 0
    for name, ex in items:
        if ex in ("_b", "_se"):
            bm = s.e.get("b")
            if bm is None:
                raise StataError(301, "last estimates not found")
            b = np.asarray(bm.data, dtype=float).ravel()
            cols = list(bm.colnames)
            eqs = list(getattr(bm, "coleq", []) or [])
            if ex == "_se":
                b = np.sqrt(np.diag(np.asarray(s.e["V"].data, dtype=float)))
            for j, c in enumerate(cols):
                eq = eqs[j] if j < len(eqs) and eqs[j] not in ("", "_") else ""
                names.append((eq + ":" if eq else "") + c)
                vals.append(float(b[j]))
            continue
        k += 1
        v = s.eval(ex)
        if isinstance(v, str):
            raise StataError(109, "type mismatch")
        names.append(name or f"_bs_{k}")
        vals.append(float(v))
    return names, np.asarray(vals)


# ---------------------------------------------------------------------------
# reamostragem
# ---------------------------------------------------------------------------

def _subset(ds: "Dataset", idx: np.ndarray) -> "Dataset":
    new = ds.copy()
    for v in new.vars:
        v.raw = v.raw[idx]
    new.nobs = len(idx)
    new.sortlist = []
    new.esample = None
    return new


def _draw(n: int) -> np.ndarray:
    from ..core import rng
    u = rng.RNG.uniform(n)
    return np.minimum((u * n).astype(np.int64), n - 1)


def _sample_index(ds: "Dataset", strata: np.ndarray | None, clusters: np.ndarray | None,
                  size: int | None) -> np.ndarray:
    n = ds.nobs
    groups = [np.arange(n)] if strata is None else [np.flatnonzero(strata == v) for v in np.unique(strata)]
    out = []
    for g in groups:
        if clusters is None:
            m = size if (size and strata is None) else len(g)
            out.append(g[_draw(len(g))[:m] if m <= len(g) else _draw_many(len(g), m)])
        else:
            ids = np.unique(clusters[g])
            pick = ids[_draw(len(ids))]
            for c in pick:
                out.append(g[clusters[g] == c])
    return np.concatenate(out) if out else np.arange(0)


def _draw_many(n: int, m: int) -> np.ndarray:
    from ..core import rng
    u = rng.RNG.uniform(m)
    return np.minimum((u * n).astype(np.int64), n - 1)


# ---------------------------------------------------------------------------
# prefixo
# ---------------------------------------------------------------------------

def cmd_bootstrap(s: "Session", args: str) -> None:
    colon = find_top(args, ":")
    if colon == -1:
        raise StataError(198, "invalid syntax")
    head, command = args[:colon].strip(), args[colon + 1:].strip()
    if not command:
        raise StataError(198, "invalid syntax")
    c = find_top(head, ",")
    exp_text, opt_text = (head, "") if c == -1 else (head[:c], head[c + 1:])
    o = match_options(opt_text, _OPTS) if opt_text.strip() else {}
    run(s, command, parse_exp_list(exp_text), o)


def run(s: "Session", command: str, items: list[tuple[str, str]], o: dict) -> None:
    from ..commands.registry import lookup
    out = s.output
    reps = int(o["reps"]) if isinstance(o.get("reps"), str) else 50
    if reps < 2:
        raise StataError(198, "reps() must be an integer greater than 1")   # VERIFY
    level = float(o["level"]) if isinstance(o.get("level"), str) else float(s.settings.get("level", 95))
    if isinstance(o.get("seed"), str):
        s.interp._execute_expanded(f"set seed {o['seed'].strip()}")
    word = command.split(None, 1)[0]
    spec = lookup(word)
    cmdname = spec.name if spec else word
    noisy = bool(o.get("noisily"))

    # 1) run with the original data
    _quiet(s, command, noisy)
    eclass = "b" in s.e and getattr(s, "_estimates", None) is not None and \
        getattr(s._estimates, "cmd", "") == cmdname
    if not items:
        if not eclass:
            raise StataError(198, "expression list required")   # VERIFY
        items = [("", "_b")]
    eclass_b = eclass and items == [("", "_b")]
    names, observed = _evaluate(s, items)
    base = s.data
    sample = np.ones(base.nobs, dtype=bool)
    if eclass and base.esample is not None and len(base.esample) == base.nobs:
        sample = base.esample.copy()
        if not o.get("_from_vce"):
            # with vce(bootstrap) Stata does not show this line (compat 0803)
            out.write(f"(running {cmdname} on estimation sample)\n", "text")
    if not eclass:
        _warn_rclass(s, cmdname)
    est0 = getattr(s, "_estimates", None) if eclass else None
    e0 = dict(s.e)
    r0 = dict(s.r)
    work = _subset(base, np.flatnonzero(sample))
    strata = clusters = None
    if isinstance(o.get("strata"), str):
        sv = work.get(s.expand_varlist(o["strata"])[0])
        strata = np.asarray(sv.raw if sv.is_string else sv.data)
    if isinstance(o.get("cluster"), str):
        cv = work.get(s.expand_varlist(o["cluster"])[0])
        clusters = np.asarray(cv.raw if cv.is_string else cv.data)
    size = int(o["size"]) if isinstance(o.get("size"), str) else None

    # 2) replicates
    dots = not o.get("nodots")
    if dots:
        out.write(f"\nBootstrap replications ({reps})\n", "text")
        out.write("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 \n", "text")
    R = np.full((reps, len(names)), np.nan)
    mis = 0
    try:
        for r in range(reps):
            idx = _sample_index(work, strata, clusters, size)
            s.data = _subset(work, idx)
            ok = True
            try:
                _quiet(s, command, noisy)
                nm, vals = _evaluate(s, items)
                if len(vals) != len(names):
                    ok = False
                elif o.get("double"):
                    R[r] = vals
                else:
                    # replicates stored as float, as Stata does without the double option
                    R[r] = np.asarray(vals, dtype=np.float32)
            except StataError:
                ok = False
            if not ok:
                mis += 1
            if dots:
                out.write("." if ok else "x", "text")
                if (r + 1) % 50 == 0:
                    out.write(f"{r + 1:>6}\n", "text")
        if dots and reps % 50:
            out.write("\n", "text")
    finally:
        s.data = base
    s.r = r0

    good = R[~np.isnan(R).any(axis=1)]
    nrep = len(good)
    if nrep < 2:
        raise StataError(2000, "insufficient observations to compute bootstrap standard errors")   # VERIFY
    V = np.cov(good, rowvar=False, ddof=1).reshape(len(names), len(names))
    info = BootInfo(reps=R, names=names, observed=observed, command=command, exps=[e for _, e in items],
                    n_reps=nrep, n_mis=mis, level=level, eclass_b=eclass_b, N=int(sample.sum()),
                    strata=o.get("strata", "") if isinstance(o.get("strata"), str) else "",
                    cluster=o.get("cluster", "") if isinstance(o.get("cluster"), str) else "")
    _post(s, info, V, est0, e0, cmdname, sample)
    display(s, s._estimates, header=not o.get("noheader"), table=not o.get("notable"))


def _wrap(text: str, width: int, first: str, rest: str) -> list[str]:
    """Wraps keeping the two spaces after a period, with indentation."""
    parts = re.split(r"( +)", text)
    lines, cur = [], first
    started = False
    for k in range(0, len(parts), 2):
        word = parts[k]
        sep = parts[k - 1] if k > 0 else ""
        if started and len(cur) + len(sep) + len(word) > width:
            lines.append(cur)
            cur = rest + word
        else:
            cur = cur + (sep if started else "") + word
        started = True
    lines.append(cur)
    return lines


def _warn_rclass(s: "Session", cmdname: str) -> None:
    """bootstrap warning for a command that does not store e(sample) (compat 0803)."""
    width = int(s.settings.get("linesize", 80)) - 1
    out = s.output
    p1 = (f"Because {cmdname} is not an estimation command or does not set e(sample), bootstrap has no "
          "way to determine which observations are used in calculating the statistics and so assumes "
          "that all observations are used.  This means that no observations will be excluded from the "
          "resampling because of missing values or other reasons.")
    p2 = ("If the assumption is not true, press Break, save the data, and drop the observations that "
          "are to be excluded.  Be sure that the dataset in memory contains only the relevant data.")
    out.write("\n", "text")
    for ln in _wrap(p1, width, "Warning:  ", " " * 10):
        out.write(ln + "\n", "text")
    out.write("\n", "text")
    for ln in _wrap(p2, width, " " * 10, " " * 10):
        out.write(ln + "\n", "text")


def _quiet(s: "Session", command: str, noisy: bool) -> None:
    out = s.output
    if not noisy:
        out.quiet_depth += 1
    try:
        s.interp._execute_expanded(command)
    finally:
        if not noisy:
            out.quiet_depth -= 1


# ---------------------------------------------------------------------------
# e() and display
# ---------------------------------------------------------------------------

def _post(s: "Session", info: BootInfo, V: np.ndarray, est0, e0: dict, cmdname: str,
          sample: np.ndarray) -> None:
    from ..commands.matrix import Matrix
    from .results import post
    good = info.reps[~np.isnan(info.reps).any(axis=1)]
    se = np.sqrt(np.diag(V))
    b_bs = good.mean(axis=0)
    z = st.norm.ppf(0.5 + info.level / 200)
    a = (100 - info.level) / 2
    pct_lo = np.percentile(good, a, axis=0, method="linear")   # VERIFY definition of the percentile
    pct_hi = np.percentile(good, 100 - a, axis=0, method="linear")
    # bias-corrected: z0 = Φ⁻¹(proportion of replicates below the observed value)
    with np.errstate(divide="ignore", invalid="ignore"):
        prop = np.clip((good < info.observed).mean(axis=0), 1e-12, 1 - 1e-12)
    z0 = st.norm.ppf(prop)
    p1 = st.norm.cdf(2 * z0 - z) * 100
    p2 = st.norm.cdf(2 * z0 + z) * 100
    bc_lo = np.array([np.percentile(good[:, j], p1[j]) for j in range(len(info.names))])
    bc_hi = np.array([np.percentile(good[:, j], p2[j]) for j in range(len(info.names))])
    info.extra = {"se": se, "b_bs": b_bs, "bias": b_bs - info.observed,
                  "normal": (info.observed - z * se, info.observed + z * se),
                  "percentile": (pct_lo, pct_hi), "bc": (bc_lo, bc_hi)}
    names = info.names
    cols = [n.split(":")[-1] for n in names]
    eqs = [n.split(":")[0] if ":" in n else "_" for n in names]

    def mat(row: np.ndarray, rown: str = "y1") -> Matrix:
        m = Matrix(np.asarray(row, dtype=float).reshape(1, -1), [rown], cols)
        m.coleq = eqs
        return m

    def ci(pair) -> Matrix:
        m = Matrix(np.vstack(pair), ["ll", "ul"], cols)
        m.coleq = eqs
        return m

    extra_mats = [("b_bs", mat(b_bs)), ("reps", mat(np.full(len(names), info.n_reps, dtype=float))),
                  ("bias", mat(b_bs - info.observed)), ("se", mat(se)),
                  ("z0", mat(z0)), ("ci_normal", ci(info.extra["normal"])),
                  ("ci_percentile", ci(info.extra["percentile"])), ("ci_bc", ci(info.extra["bc"]))]
    if info.eclass_b and est0 is not None:
        est = est0
        est.V = V
        est.vcetype = "Bootstrap"
        est.stat = "z"
        est.boot = info
        est.title = getattr(est0, "title", "") or e0.get("title", "")
        scalars = [(k, v) for k, v in e0.items() if isinstance(v, float) and k not in ("N_reps", "N_misreps")]
        scalars += [("N_reps", float(info.n_reps)), ("N_misreps", float(info.n_mis)), ("level", info.level)]
        k = int(np.count_nonzero(np.abs(np.diag(V)) > 0))
        _wald(info, V)
        scalars.append(("chi2", info.extra.get("chi2", SYS)))
        macros = [(k_, v) for k_, v in e0.items() if isinstance(v, str)
                  and k_ not in ("vce", "vcetype", "prefix")]
        macros += [("prefix", "bootstrap"), ("vce", "bootstrap"), ("vcetype", "Bootstrap"),
                   ("command", info.command)]
        # Stata lists in the reverse order of storage; here the lists already come in that order
        post(s, est, list(reversed(scalars)), list(reversed(macros)), sample, extra_mats=extra_mats)
        return
    est = Estimates(cmd="bootstrap", depvar="", names=cols, b=np.asarray(info.observed, dtype=float),
                    V=V, N=info.N, stat="z")
    est.eqnames = eqs
    est.rows = [CoefRow("coef", label=n.split(":")[-1], b=float(b), se=float(v), name=n)
                for n, b, v in zip(names, info.observed, se)]
    est.vcetype = "Bootstrap"
    est.boot = info
    est.title = "Bootstrap results"
    est.level = info.level
    scalars = [("N", float(info.N)), ("N_reps", float(info.n_reps)), ("N_misreps", float(info.n_mis)),
               ("k_eq", 1.0), ("k_exp", float(len(names))), ("level", info.level)]
    macros = [("cmd", "bootstrap"), ("command", info.command), ("prefix", "bootstrap"),
              ("vce", "bootstrap"), ("vcetype", "Bootstrap"), ("title", "Bootstrap results"),
              ("properties", "b V")]
    for j, ex in enumerate(info.exps, start=1):
        macros.append((f"exp{j}", ex))
    post(s, est, scalars, macros, sample, extra_mats=extra_mats)


def _wald(info: BootInfo, V: np.ndarray) -> None:
    keep = [j for j, n in enumerate(info.names) if not n.endswith("_cons") and V[j, j] > 0]
    if not keep:
        return
    b = info.observed[keep]
    Vk = V[np.ix_(keep, keep)]
    try:
        info.extra["chi2"] = float(b @ np.linalg.solve(Vk, b))
        info.extra["df"] = len(keep)
    except np.linalg.LinAlgError:
        pass


def _hdr(s: "Session", left: str, items: list[tuple[str, str]]) -> None:
    out = s.output
    for k, (lab, val) in enumerate(items):
        lft = left if k == 0 else ""
        out.write(f"{lft:<48}{lab:<18}= ", "text")
        out.write(f"{val:>10}\n", "result")


def display(s: "Session", est, *, header: bool = True, table: bool = True, level: float | None = None,
            **_) -> None:
    info: BootInfo = est.boot
    out = s.output
    lev = level if level is not None else info.level
    if header:
        out.write("\n", "text")
        items = [("Number of obs", f"{info.N:,}"), ("Replications", f"{info.n_reps:,}")]
        if info.eclass_b:
            if "chi2" in info.extra:
                df = info.extra["df"]
                items.append((f"Wald chi2({df})", f"{info.extra['chi2']:.2f}"))
                items.append(("Prob > chi2", f"{_chi2_tail(info.extra['chi2'], df):.4f}"))
            for key, lab, fmt in (("r2", "R-squared", "%9.4f"), ("r2_a", "Adj R-squared", "%9.4f"),
                                  ("r2_p", "Pseudo R2", "%9.4f"), ("rmse", "Root MSE", "%9.4f")):
                v = s.e.get(key)
                if isinstance(v, float) and v < SYS:
                    items.append((lab, g9(v, fmt)))
        _hdr(s, est.title if info.eclass_b else "Bootstrap results", items)
        if info.cluster:
            out.write(f"{'':>48}(Replications based on clusters in {info.cluster})\n", "text")   # VERIFY
        if not info.eclass_b:
            out.write("\n", "text")
            out.write(f"{'command':>13}:  {info.command}\n", "text")
            k = 0
            for n, ex in zip(info.names, info.exps):
                out.write(f"{abbrev(n, 12):>13}:  {ex}\n", "text")
                k += 1
        out.write("\n", "text")
    if not table:
        return
    se = info.extra["se"]
    rows = []
    if info.eclass_b:
        se_by = dict(zip(info.names, se))
        se_short = {n.split(":")[-1]: v for n, v in se_by.items()}
        for r in est.rows:
            if r.kind in ("coef", "aux"):
                key = r.name or r.label
                v = se_by.get(key, se_short.get(key.split(":")[-1], r.se))
                rows.append(CoefRow(r.kind, r.label, r.b, float(v), r.eq, r.name))
            else:
                rows.append(r)
    else:
        rows = [CoefRow("coef", label=n.split(":")[-1], b=float(b), se=float(v), name=n)
                for n, b, v in zip(info.names, info.observed, se)]
    dep = est.depvar if info.eclass_b else ""
    coef_table(s, dep, rows, stat="z", df=None, level=lev,
               vcetype="", top_line=f"{'':>12} |   Observed   Bootstrap                         Normal-based")


# ---------------------------------------------------------------------------
# estat bootstrap
# ---------------------------------------------------------------------------

def estat_bootstrap(s: "Session", est, rest: str) -> None:
    info: BootInfo | None = getattr(est, "boot", None)
    if info is None:
        raise StataError(321, "estat bootstrap requires bootstrap results")   # VERIFY
    c = find_top(rest, ",")
    opt = "" if c == -1 else rest[c + 1:]
    o = match_options(opt, {"bc": 2, "bca": 3, "normal": 1, "percentile": 1, "all": 3, "noheader": 5,
                            "nolegend": 5, "verbose": 1}) if opt.strip() else {}
    kinds = [k for k in ("normal", "percentile", "bc") if o.get(k) or o.get("all")]
    if not kinds:
        kinds = ["percentile", "bc"] if not o.get("bc") else ["bc"]
        if o.get("percentile"):
            kinds = ["percentile"]
    out = s.output
    lv = f"{info.level:g}"
    out.write("\n", "text")
    _hdr(s, est.title if info.eclass_b else "Bootstrap results",
         [("Number of obs", f"{info.N:,}"), ("Replications", f"{info.n_reps:,}")])
    out.write("\n" + LINE + "\n", "text")
    out.write(f"{'':>12} |    Observed               Bootstrap\n", "text")
    dep = est.depvar if info.eclass_b else ""
    out.write(f"{abbrev(dep):>12} |       Coef.       Bias    Std. Err.  [{lv}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    tags = {"normal": "(N)", "percentile": "(P)", "bc": "(BC)"}
    ex = info.extra
    for j, n in enumerate(info.names):
        first = True
        for k in kinds:
            lo, hi = ex[k][0][j], ex[k][1][j]
            # columns (compat 0803): 12, 11, 12, 12, 11 and the label in 6
            if first:
                out.write(f"{abbrev(n.split(':')[-1]):>12} |", "text")
                out.write(f"{g9(info.observed[j], '%10.0g'):>12}{g9(ex['bias'][j]):>11}"
                          f"{g9(ex['se'][j], '%10.0g'):>12}{g9(lo):>12}{g9(hi):>11}{tags[k]:>6}\n", "result")
                first = False
            else:
                out.write(f"{'':>12} |" + " " * 35, "text")
                out.write(f"{g9(lo):>12}{g9(hi):>11}{tags[k]:>6}\n", "result")
    out.write(LINE + "\n", "text")
    legend = {"normal": "(N)    normal confidence interval",
              "percentile": "(P)    percentile confidence interval",
              "bc": "(BC)   bias-corrected confidence interval"}
    for k in kinds:
        out.write(legend[k] + "\n", "text")


# ---------------------------------------------------------------------------
# vce(bootstrap) in the estimation commands
# ---------------------------------------------------------------------------

_VCE = re.compile(r"\bvce\s*\(\s*(boot(?:s|st|str|stra|strap)?)\b", re.I)


def split_vce_bootstrap(args: str) -> tuple[str, str] | None:
    """If the options have vce(bootstrap [, options]), returns (args without the option,
    bootstrap options)."""
    c = find_top(args, ",")
    if c == -1:
        return None
    opts = args[c + 1:]
    m = _VCE.search(opts)
    if not m:
        return None
    start = m.start()
    depth, j = 0, opts.index("(", start)
    k = j
    while k < len(opts):
        if opts[k] == "(":
            depth += 1
        elif opts[k] == ")":
            depth -= 1
            if depth == 0:
                break
        k += 1
    inner = opts[j + 1:k]
    bopts = inner.partition(",")[2].strip()
    rest = (opts[:start] + opts[k + 1:]).strip()
    new_args = args[:c] + ((", " + rest) if rest else "")
    return new_args, bopts


def wrap_vce_bootstrap(name: str) -> None:
    """Replaces the registered command with one that understands vce(bootstrap)."""
    from ..commands.registry import REGISTRY, CommandSpec
    spec = REGISTRY[name]
    inner = spec.fn

    def fn(s: "Session", args: str, _inner=inner, _name=name) -> None:
        sp = split_vce_bootstrap(args)
        if sp is None:
            _inner(s, args)
            return
        new_args, bopts = sp
        o = match_options(bopts, _OPTS) if bopts else {}
        o["_from_vce"] = True
        run(s, f"{_name} {new_args}", [], o)

    REGISTRY[name] = CommandSpec(spec.name, spec.min_abbrev, fn, spec.prefix, spec.byable)
