"""Least-squares systems and variants: constraint, cnsreg, eivreg, vwls, sureg, reg3 and rreg
([R] constraint, [R] cnsreg, [R] eivreg, [R] vwls, [R] sureg, [R] reg3, [R] rreg). The numerics are
in core/syseq.py; here the syntax, the data, Stata's output and e().

Layouts come from the manuals; the lines marked VERIFY wait for the Stata 14 log of case 1908."""

from __future__ import annotations

import re
from types import SimpleNamespace
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import syseq
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from . import fvars
from .mleq import Block, assemble, block_from_text
from .models2 import _header, g
from .regress import cmdline, comma, sig5
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# constraint
# ---------------------------------------------------------------------------

def _constraints(s: "Session") -> dict:
    if not hasattr(s, "constraints"):
        s.constraints = {}
    return s.constraints


def cmd_constraint(s: "Session", args: str) -> None:
    """constraint define|dir|list|drop ([R] constraint). The expression is kept as text and is
    read against the coefficients of the command that uses it (cnsreg, ...)."""
    cons = _constraints(s)
    m = re.match(r"^(\w+)\s*(.*)$", args.strip(), re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    sub, rest = m.group(1).lower(), m.group(2).strip()
    if "define".startswith(sub) and len(sub) >= 1 and sub not in ("dir", "di", "d") or sub == "def":
        mm = re.match(r"^(\d+)\s+(.+)$", rest, re.S)
        if not mm:
            raise StataError(198, "invalid syntax")
        n = int(mm.group(1))
        if not 1 <= n <= 1999:
            raise StataError(125, "constraint number must be between 1 and 1999")   # VERIFY
        cons[n] = mm.group(2).strip()
        return
    if sub in ("dir", "di", "list", "l", "li", "lis"):
        wanted = _numlist(rest) if rest and rest != "_all" else sorted(cons)
        for n in sorted(wanted):
            if n in cons:
                s.output.write(f"{n:>6}:  {cons[n]}\n", "text")
        return
    if sub == "drop":
        if rest == "_all":
            cons.clear()
            return
        for n in _numlist(rest):
            if n not in cons:
                raise StataError(111, f"constraint {n} not found")   # VERIFY
            del cons[n]
        return
    raise StataError(199, f"unrecognized command: constraint {sub}")


def _numlist(text: str) -> list[int]:
    out: list[int] = []
    for tok in text.replace(",", " ").split():
        if "-" in tok:
            a, b = tok.split("-", 1)
            out += list(range(int(a), int(b) + 1))
        else:
            out.append(int(tok))
    return out


# ---------------------------------------------------------------------------
# sample of the commands of this module
# ---------------------------------------------------------------------------

class _Sample:
    """If/in, weights and options of a command whose varlist this module reads itself."""

    def __init__(self, s: "Session", tail: str, spec: dict[str, int], weights=("aweight", "fweight", "iweight"),
                 var_terms=()):
        from ..commands._util import touse
        from ..commands.summarize import weights as get_weights
        self.s = s
        p = parse_standard(tail)
        spec = dict(spec)
        spec.update({"level": 1, "noconstant": 4})
        self.o = match_options(p.options, spec) if p.options.strip() else {}
        ds = s.data
        mask = touse(s, p)
        bad = fvars.term_missing(ds, list(var_terms)) if var_terms else np.zeros(ds.nobs, dtype=bool)
        mask &= ~bad
        self.w, self.wtype, mask = get_weights(s, p, mask, weights)
        self.mask = mask
        self.n = int(mask.sum())
        if self.n == 0:
            raise StataError(2000, "no observations")
        self.level = float(self.o["level"]) if self.o.get("level") else float(s.settings.get("level", "95"))
        if self.level < 10 or self.level > 99.99:
            raise StataError(198, "level() must be between 10 and 99.99 inclusive")
        self.constant = not self.o.get("noconstant")
        self.vce = "ols"
        self.clustvar = ""
        self.depname = ""


def _equations(text: str) -> tuple[list[tuple[str, str, str]], str]:
    """'(a b c) (n: d = e f) if x, opts' -> [(name, depvar, rhs)], 'if x, opts'."""
    t = text.strip()
    eqs, depth, start, end = [], 0, None, 0
    for i, ch in enumerate(t):
        if ch == "(":
            if depth == 0:
                start = i
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0 and start is not None:
                eqs.append(t[start + 1:i])
                end = i + 1
        elif depth == 0 and ch not in " \t":
            break
    out = []
    for body in eqs:
        name = ""
        mm = re.match(r"^\s*(\w+)\s*:\s*(.*)$", body, re.S)
        if mm:
            name, body = mm.group(1), mm.group(2)
        if "=" in body:
            dep, rhs = body.split("=", 1)
        else:
            toks = fvars._split_tokens(body)
            if not toks:
                raise StataError(100, "varlist required")
            dep, rhs = toks[0], " ".join(toks[1:])
        out.append((name, dep.strip(), rhs.strip()))
    if not out:
        raise StataError(198, "equation required")   # VERIFY
    return out, t[end:].strip()


def _dep(s: "Session", name: str):
    ds = s.data
    terms = fvars.expand_fv(ds, name)
    if len(terms) != 1 or terms[0].is_interaction or terms[0].comps[0].factor:
        raise StataError(198, "depvar may not be a factor variable")   # VERIFY
    return terms[0], terms[0].comps[0]


def _fake_pr(s, mask, n):
    return SimpleNamespace(s=s, smp=SimpleNamespace(mask=mask, n=n), w=np.ones(n))


def _est_names(blk: Block) -> list[str]:
    names = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in blk.cols if c.values is not None]
    if blk.constant:
        names.append("_cons")
    return names


# ---------------------------------------------------------------------------
# cnsreg
# ---------------------------------------------------------------------------

def cmd_cnsreg(s: "Session", args: str) -> None:
    from .postest import _constraint, replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "cnsreg", args)
        return
    from .regress import Sample, wexp_macros
    smp = Sample(s, args, opts_spec={"constraints": 3, "collinear": 4, "noheader": 6, "notable": 5},
                 weights=("aweight", "fweight", "pweight", "iweight"))
    o = smp.o
    if not o.get("constraints"):
        raise StataError(100, "option constraints() required")
    ds = s.data
    w, N = smp.weights_for_fit()
    notes = fvars.drop_collinear(smp.cols, smp.constant, w)
    for nm in notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    est_cols = [c for c in smp.cols if c.values is not None]
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(smp.n)] if smp.constant else []))
    names = [c.name for c in est_cols] + (["_cons"] if smp.constant else [])
    fake = SimpleNamespace(names=names, eqnames=[])
    cons = _constraints(s)
    R, r, texts = [], [], []
    for n in _numlist(str(o["constraints"])):
        if n not in cons:
            raise StataError(111, f"constraint {n} not found")
        coef, rhs = _constraint(s, fake, cons[n])
        row = np.zeros(len(names))
        for k, v in coef.items():
            row[k] = v
        R.append(row)
        r.append(rhs)
        from .postest import _ctext
        texts.append(_ctext(fake, coef, rhs))
    R, r = np.array(R).reshape(-1, len(names)), np.array(r)
    beta, A, rk = syseq.restricted_ls(X, smp.y, w, R, r)
    resid = smp.y - X @ beta
    rss = float(np.sum(w * resid ** 2))
    k = X.shape[1]
    df_r = N - k + rk
    if df_r <= 0:
        raise StataError(2001, "insufficient observations")
    s2 = rss / df_r
    Vb = s2 * A
    if smp.constant:
        ybar = float(np.sum(w * smp.y) / np.sum(w))
        tss = float(np.sum(w * (smp.y - ybar) ** 2))
    else:
        tss = float(np.sum(w * smp.y ** 2))
    slopes = list(range(len(est_cols)))
    df_m = len(slopes) - rk
    F = SYS
    if slopes and df_m > 0:
        bb, VV = beta[slopes], Vb[np.ix_(slopes, slopes)]
        F = float(bb @ np.linalg.pinv(VV) @ bb) / df_m
    rmse = float(np.sqrt(s2))
    kk = len(smp.cols) + (1 if smp.constant else 0)
    names_full = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols] \
        + (["_cons"] if smp.constant else [])
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([kk - 1] if smp.constant else [])
    b_full, V_full = np.zeros(kk), np.zeros((kk, kk))
    b_full[pos] = beta
    V_full[np.ix_(pos, pos)] = Vb
    bmap = {c.name: b_full[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(V_full[i, i], 0))) for i, c in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b_full[-1], float(np.sqrt(max(V_full[-1, -1], 0)))
    est = Estimates("cnsreg", smp.depname, names_full, b_full, V_full, smp.n, "t", df_r)
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level, est.vcetype, est.clustvar, est.n_clust = smp.level, "", "", None
    est.extra = {"N": N, "df_m": df_m, "df_r": df_r, "F": F, "rmse": rmse, "rss": rss, "tss": tss,
                 "r2": 1 - rss / tss if tss > 0 else SYS, "ll": SYS, "rank": k - rk}
    est.wtype = smp.wtype
    est.cns_text = texts
    est.display = _display_cnsreg
    scalars = [("N", N), ("df_m", df_m), ("df_r", df_r), ("F", F), ("rmse", rmse), ("rss", rss),
               ("ll", -0.5 * N * (np.log(2 * np.pi) + np.log(rss / N) + 1)), ("r2", est.extra["r2"]),
               ("rank", k - rk)]
    macros = [("cmdline", cmdline("cnsreg", args)), ("title", "Constrained linear regression"),
              ("marginsok", "XB default"), ("depvar", smp.depname), ("cmd", "cnsreg"),
              ("properties", "b V"), ("predict", "regres_p")] + wexp_macros(smp)
    post(s, est, scalars, macros, smp.mask)
    _display_cnsreg(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _flabel(df_m, df_r, w_m=3, w_r=10) -> str:
    return f"F({int(df_m):>{w_m}},{int(df_r):>{w_r}})"


def _display_cnsreg(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        out = s.output
        pF = g(st.f.sf(x["F"], x["df_m"], x["df_r"]) if x["F"] < SYS else SYS, "%10.4f")
        out.write("\n", "text")
        for lft, lab, val in (("Constrained linear regression", "Number of obs", comma(x["N"])),
                              ("", _flabel(x["df_m"], x["df_r"], 4, 7), g(x["F"], "%10.2f")),
                              ("", "Prob > F", pF), ("", "Root MSE", sig5(x["rmse"]))):
            out.write(f"{lft:<48}{lab:<18}= ", "text")
            out.write(f"{val:>10}\n", "result")
        out.write("\n", "text")
        for i, t in enumerate(est.cns_text, 1):
            s.output.write(f" ({i:>2})  {t}\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=lev)


# ---------------------------------------------------------------------------
# eivreg
# ---------------------------------------------------------------------------

def cmd_eivreg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "eivreg", args)
        return
    from .regress import Sample, wexp_macros
    smp = Sample(s, args, opts_spec={"reliab": 1, "noheader": 6, "notable": 5}, weights=("aweight", "fweight"))
    o = smp.o
    ds = s.data
    est_cols = [c for c in smp.cols if c.values is not None]
    rel = np.ones(len(est_cols))
    if o.get("reliab"):
        for part in re.findall(r"\(([^)]*)\)", str(o["reliab"])) or [str(o["reliab"])]:
            mm = re.match(r"^\s*(\S+)\s+([0-9.eE+-]+)\s*$", part)
            if not mm:
                raise StataError(198, "option reliab() invalid")   # VERIFY
            nm = mm.group(1)
            val = float(mm.group(2))
            if not 0 < val <= 1:
                raise StataError(198, "reliability must be in (0, 1]")   # VERIFY
            idx = [i for i, c in enumerate(est_cols) if c.name == nm]
            if not idx:
                raise StataError(111, f"{nm} not found")
            rel[idx[0]] = val
    if not np.all(np.isfinite(rel)) or smp.n < len(est_cols) + 2:
        raise StataError(2001, "insufficient observations")
    Xs = np.column_stack([c.values for c in est_cols]) if est_cols else np.zeros((smp.n, 0))
    beta, a, V, see = syseq.eivreg(Xs, smp.y, rel, smp.constant)
    N = float(smp.n)
    p = len(est_cols)
    df_r = N - p - (1 if smp.constant else 0)
    r2 = max(1 - see * df_r / ((N - 1) * float(np.var(smp.y, ddof=1))), 0.0) if smp.constant else SYS   # VERIFY: matches the 1908 log
    kk = len(smp.cols) + (1 if smp.constant else 0)
    names_full = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols] \
        + (["_cons"] if smp.constant else [])
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([kk - 1] if smp.constant else [])
    b_vec = np.concatenate([beta, [a]]) if smp.constant else beta
    b_full, V_full = np.zeros(kk), np.zeros((kk, kk))
    b_full[pos] = b_vec
    V_full[np.ix_(pos, pos)] = V
    bmap = {c.name: b_full[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(V_full[i, i], 0))) for i, c in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b_full[-1], float(np.sqrt(max(V_full[-1, -1], 0)))
    est = Estimates("eivreg", smp.depname, names_full, b_full, V_full, smp.n, "t", df_r)
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level, est.vcetype, est.clustvar, est.n_clust = smp.level, "", "", None
    sl = list(range(p))
    F = float(beta @ np.linalg.solve(V[np.ix_(sl, sl)], beta)) / p if p else SYS
    rmse = float(np.sqrt(see)) if df_r > 0 else SYS   # VERIFY: Root MSE is sqrt(sigma_ee)
    relrows = [(c.name, float(rel[i])) for i, c in enumerate(est_cols) if o.get("reliab") and rel[i] != 1.0 or False]
    if len(relrows) < len(est_cols):
        relrows.append(("*", 1.0))
    est.extra = {"N": N, "df_m": p, "df_r": df_r, "F": F, "r2": r2, "rmse": rmse, "rel": relrows}
    est.display = _display_eivreg
    scalars = [("N", N), ("df_m", p), ("df_r", df_r), ("F", F), ("r2", r2), ("rmse", rmse)]
    macros = [("cmdline", cmdline("eivreg", args)), ("title", "Errors-in-variables regression"),
              ("depvar", smp.depname), ("cmd", "eivreg"), ("properties", "b V"), ("predict", "regres_p")] \
        + wexp_macros(smp)
    post(s, est, scalars, macros, smp.mask)
    _display_eivreg(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _display_eivreg(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        out = s.output
        left = [" " * 19 + "assumed", "    variable     reliability", "-" * 28]
        left += [f"{nm:>12}{v:>13.4f}" for nm, v in x["rel"]]
        right = [("", ""), ("", ""), ("Number of obs", comma(x["N"])),
                 (_flabel(x["df_m"], x["df_r"], 3, 10), g(x["F"], "%10.2f")),
                 ("Prob > F", g(st.f.sf(x["F"], x["df_m"], x["df_r"]) if x["F"] < SYS else SYS, "%10.4f")),
                 ("R-squared", g(x["r2"], "%10.4f")), ("Root MSE", g(x["rmse"], "%8.0g"))]
        out.write("\n", "text")
        out.write(f"{left[0]:<48}Errors-in-variables regression\n", "text")
        left, right = left[1:], right[1:]
        m = max(len(left), len(right))
        left += [""] * (m - len(left))
        right += [("", "")] * (m - len(right))
        for lft, (rl, rv) in zip(left, right):
            if rl:
                out.write(f"{lft:<48}{rl:<18}= ", "text")
                out.write(f"{rv:>10}\n", "result")
            else:
                out.write(lft + "\n", "text")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=lev)


# ---------------------------------------------------------------------------
# vwls
# ---------------------------------------------------------------------------

def cmd_vwls(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "vwls", args)
        return
    from .regress import Sample, wexp_macros
    smp = Sample(s, args, opts_spec={"sd": 2, "noheader": 6, "notable": 5}, weights=(), extra_vars="")
    o = smp.o
    ds = s.data
    if not o.get("sd"):
        raise StataError(100, "option sd() required")
    from ..core.varlist import resolve_name
    sdname = resolve_name(ds, str(o["sd"]).strip())
    sd = ds.get(sdname).data[smp.mask].astype(float)
    if np.any(sd <= 0) or np.any(sd >= SYS):
        raise StataError(198, "sd() must be positive and not missing")   # VERIFY
    est_cols = [c for c in smp.cols if c.values is not None]
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(smp.n)] if smp.constant else []))
    beta, V, chi2_fit = syseq.vwls(X, smp.y, sd)
    N = float(smp.n)
    k = X.shape[1]
    df_gof = N - k
    slopes = list(range(len(est_cols)))
    df_m = len(slopes)
    chi2_m = float(beta[slopes] @ np.linalg.solve(V[np.ix_(slopes, slopes)], beta[slopes])) if slopes else SYS
    kk = len(smp.cols) + (1 if smp.constant else 0)
    names_full = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols] \
        + (["_cons"] if smp.constant else [])
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([kk - 1] if smp.constant else [])
    b_full, V_full = np.zeros(kk), np.zeros((kk, kk))
    b_full[pos] = beta
    V_full[np.ix_(pos, pos)] = V
    bmap = {c.name: b_full[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(V_full[i, i], 0))) for i, c in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b_full[-1], float(np.sqrt(max(V_full[-1, -1], 0)))
    est = Estimates("vwls", smp.depname, names_full, b_full, V_full, smp.n, "z", None)
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level, est.vcetype, est.clustvar, est.n_clust = smp.level, "", "", None
    est.extra = {"N": N, "df_gof": df_gof, "chi2_gof": chi2_fit, "df_m": df_m, "chi2": chi2_m}
    est.display = _display_vwls
    scalars = [("N", N), ("df_m", df_m), ("chi2", chi2_m), ("chi2_gf", chi2_fit), ("df_gf", df_gof)]
    macros = [("cmdline", cmdline("vwls", args)), ("title", "Variance-weighted least-squares regression"),
              ("depvar", smp.depname), ("cmd", "vwls"), ("properties", "b V"), ("sdvar", sdname),
              ("predict", "regres_p")]
    post(s, est, scalars, macros, smp.mask)
    _display_vwls(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _display_vwls(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        out = s.output
        out.write("\n", "text")
        out.write(f"{'Variance-weighted least-squares regression':<48}{'Number of obs':<18}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        p_gof = chi2_tail(x["chi2_gof"], x["df_gof"])
        p_m = chi2_tail(x["chi2"], x["df_m"]) if x["df_m"] else SYS
        lab = f"Goodness-of-fit chi2({int(x['df_gof'])})"
        out.write(f"{lab:<27}={g(x['chi2_gof'], '%8.2f'):>8}".ljust(48)
                  + f"{'Model chi2(' + str(int(x['df_m'])) + ')':<18}= {g(x['chi2'], '%10.2f'):>10}\n", "text")
        out.write(f"{'Prob > chi2':<27}={g(p_gof, '%8.4f'):>8}".ljust(48)
                  + f"{'Prob > chi2':<18}= {g(p_m, '%10.4f'):>10}\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=lev)


# ---------------------------------------------------------------------------
# sureg and reg3
# ---------------------------------------------------------------------------

def cmd_sureg(s: "Session", args: str) -> None:
    _system(s, args, "sureg")


def cmd_reg3(s: "Session", args: str) -> None:
    _system(s, args, "reg3")


def _system(s: "Session", args: str, cmd: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return
    eqs, tail = _equations(args)
    tail = re.sub(r"(?<![\w])2sls(?![\w])", "twosls", re.sub(r"(?<![\w])3sls(?![\w])", "threesls", tail))
    ds = s.data
    deps = [_dep(s, d) for _, d, _ in eqs]
    terms = [fvars.expand_fv(ds, rhs) if rhs.strip() else [] for _, _, rhs in eqs]
    spec = {"corr": 4, "isure": 3, "iterate": 4, "small": 2, "dfk": 3, "dfk2": 4, "notable": 5, "noheader": 6,
            "constraints": 3, "nolog": 4, "tolerance": 3, "level": 1, "noconstant": 4, "log": 3,
            "nodfk": 4, "first": 3}
    if cmd == "reg3":
        spec.update({"ols": 3, "twosls": 4, "threesls": 4, "sure": 2, "ireg3": 4, "endog": 3, "exog": 3, "inst": 3,
                     "allexog": 3})
        spec.pop("isure")
    smp = _Sample(s, tail, spec, var_terms=[t for ds_ in deps for t in [ds_[0]]] + [t for tt in terms for t in tt])
    o = smp.o
    mask, n = smp.mask, smp.n
    pr = _fake_pr(s, mask, n)
    blocks, ys, depnames = [], [], []
    for (name, dep, rhs), (dterm, dcomp) in zip(eqs, deps):
        ys.append(fvars.comp_values(ds, dcomp)[mask].astype(float))
        depnames.append(dcomp.tsname)
        blocks.append(block_from_text(pr, name or dcomp.tsname, rhs, smp.constant))
    for blk in blocks:
        for nm in blk.notes:
            s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    Xs = [blk.X for blk in blocks]
    ks = [X.shape[1] for X in Xs]
    cons, cns_text = _system_constraints(s, o, blocks) if o.get("constraints") else (None, [])
    dfk = "small" if o.get("small") and cmd == "sureg" else ""
    if o.get("dfk"):
        dfk = "dfk"
    iterate = bool(o.get("isure") or o.get("ireg3"))
    method = "sure"
    if cmd == "reg3":
        method = "3sls"
        for key, opt in (("ols", "ols"), ("twosls", "2sls"), ("threesls", "3sls"), ("sure", "sure")):
            if o.get(key):
                method = opt
    if cmd == "sureg":
        fit = syseq.sure(Xs, ys, iterate=iterate, dfk=dfk, cons=cons)
        Zs = Xs
    else:
        rhs_names = []
        for blk in blocks:
            for c in blk.cols:
                if c.values is not None and c.name not in rhs_names:
                    rhs_names.append(c.name)
        all_vars = []
        for dep_, blk in zip(depnames, blocks):
            for nm in [dep_] + [c.name for c in blk.cols if c.values is not None]:
                if nm not in all_vars:
                    all_vars.append(nm)
        if o.get("inst"):
            # inst() is the complete list of exogenous variables; everything else is endogenous
            exog_list = list(dict.fromkeys(str(o["inst"]).split()))
            endog_list = [v for v in all_vars if v not in exog_list]
        else:
            endog_list = list(dict.fromkeys(list(depnames) + str(o.get("endog") or "").split()))
            exog_list = [v for v in rhs_names if v not in endog_list]
            exog_list += [v for v in dict.fromkeys(str(o.get("exog") or "").split()) if v not in exog_list]
        endog_names = set(endog_list)
        inst_cols = []
        for nm in exog_list:
            col = next((c.values for b_ in blocks for c in b_.cols if c.name == nm and c.values is not None), None)
            inst_cols.append(col if col is not None else ds.get(nm).data[mask].astype(float))
        if smp.constant:
            inst_cols.append(np.ones(n))
        Zfull = np.column_stack(inst_cols)
        P = Zfull @ np.linalg.pinv(Zfull.T @ Zfull) @ Zfull.T
        Zs = []
        for blk in blocks:
            cols = []
            for c in [c for c in blk.cols if c.values is not None]:
                cols.append(P @ c.values if c.name in endog_names else c.values)
            if blk.constant:
                cols.append(np.ones(n))
            Zs.append(np.column_stack(cols))
        # Stata scales the 2sls/ols covariance by the equation degrees of freedom (observed in the 1908 log)
        # with constraints the 2sls/ols covariance is not rescaled and the tests are z (compat 1918)
        fit = syseq.three_sls(Xs, Zs, ys, method=method, iterate=iterate,
                              dfk="small" if method in ("2sls", "ols") and not dfk and cons is None else dfk, cons=cons)
    foot = (endog_list, exog_list) if cmd == "reg3" else None
    if cmd == "sureg":
        exog = []
        for dep_, blk in zip(depnames, blocks):
            for nm in [dep_] + [c.name for c in blk.cols if c.values is not None]:
                if nm not in exog:
                    exog.append(nm)
    else:
        exog = exog_list
    _post_system(s, args, cmd, eqs, blocks, ys, Xs, fit, smp, o, method if cmd == "reg3" else "sure", dfk, foot,
                 cns_text=cns_text, cons=cons, exog=exog)


def _system_constraints(s, o, blocks):
    """(R, r) over the free coefficients of the equations, and the constraints as Stata lists
    them ([eq]name)."""
    from .postest import _constraint, _ctext
    names, eqn = [], []
    for blk in blocks:
        for c in blk.cols:
            if c.values is not None:
                names.append(c.name)
                eqn.append(blk.name)
        if blk.constant:
            names.append("_cons")
            eqn.append(blk.name)
    fake = SimpleNamespace(names=names, eqnames=eqn)
    cons = _constraints(s)
    R, r, texts = [], [], []
    for k in _numlist(str(o["constraints"])):
        if k not in cons:
            raise StataError(111, f"constraint {k} not found")
        coef, rhs = _constraint(s, fake, cons[k])
        row = np.zeros(len(names))
        for j, v in coef.items():
            row[j] = v
        R.append(row)
        r.append(rhs)
        texts.append(_ctext(fake, coef, rhs))
    return (np.array(R).reshape(-1, len(names)), np.array(r, dtype=float)), texts


def _post_system(s, args, cmd, eqs, blocks, ys, Xs, fit, smp, o, method, dfk, foot=None, cns_text=(), cons=None,
                 exog=()) -> None:
    ds = s.data
    n = smp.n
    m = len(blocks)
    b_all = np.concatenate(fit.b)
    names, eqnames, b_full, V_full, rows = assemble(SimpleNamespace(s=s), blocks, b_all, fit.V, [])
    use_t = bool(o.get("small")) or (cmd == "reg3" and method in ("2sls", "ols") and not cns_text)
    if cmd == "sureg" and o.get("isure") and not o.get("nolog"):
        s.output.write("\n", "text")
        for it, d in enumerate(fit.trace, 1):
            s.output.write(f"Iteration {it}:   tolerance = {g(d, '%10.0g'):>10}\n", "text")
    eq_info = []
    for i, blk in enumerate(blocks):
        e = fit.resid[:, i]
        rss = float(e @ e)
        y = ys[i]
        tss = float(np.sum((y - y.mean()) ** 2)) if blk.constant else float(y @ y)
        k = blk.k
        denom = n - k if use_t or dfk == "small" else n
        rmse = float(np.sqrt(rss / denom))
        r2 = 1 - rss / tss if tss > 0 else SYS
        off = int(sum(b.k for b in blocks[:i]))
        sl = list(range(off, off + blk.nslopes))
        if sl:
            bb, VV = b_all[sl], fit.V[np.ix_(sl, sl)]
            chi2 = float(bb @ np.linalg.solve(VV, bb))
        else:
            chi2 = SYS
        eq_info.append({"obs": n, "parms": blk.nslopes, "rmse": rmse, "r2": r2, "chi2": chi2, "df": len(sl), "rss": rss,
                        "mss": tss - rss})
    level = smp.level
    est = Estimates(cmd, blocks[0].name, names, b_full, V_full, n, "t" if use_t else "z",
                    float(m * n - sum(b.k for b in blocks)) if use_t else None)
    est.foot = foot
    est.cns_text = list(cns_text)
    est.iterated = bool(o.get("isure") or o.get("ireg3"))
    est.eqnames = eqnames
    est.rows = rows
    est.level, est.vcetype, est.clustvar, est.n_clust = level, "", "", None
    est.constant = smp.constant
    est.blocks = blocks
    est.depvar_cell = ""
    est.fit = fit
    est.eq_info = eq_info
    est.method = method
    est.names_eq = [b.name for b in blocks]
    est.corr_opt = bool(o.get("corr"))
    est.cons = smp.constant
    est.extra = {"N": n}
    est.display = _display_system
    est.depvar = " ".join(b.name for b in blocks)
    title = {"sureg": "Seemingly unrelated regression",
             "reg3": {"3sls": "Three-stage least-squares regression", "2sls": "Two-stage least-squares regression",
                      "ols": "Ordinary least-squares regression",
                      "sure": "Seemingly unrelated regression"}.get(method, "Three-stage least-squares regression")
             }[cmd]
    est.title = title + (", iterated " if cmd == "sureg" and o.get("isure") else "")
    # e() as Stata lists it after sureg (compat 1918; reg3 follows it, VERIFY)
    scalars = []
    for i, (blk, e) in enumerate(zip(blocks, eq_info), 1):
        p_i = chi2_tail(e["chi2"], e["df"]) if e["chi2"] < SYS else SYS
        scalars += [(f"cons_{i}", 1.0 if blk.constant else 0.0), (f"chi2_{i}", e["chi2"]), (f"p_{i}", p_i),
                    (f"df_m{i}", e["df"]), (f"rss_{i}", e["rss"]), (f"rmse_{i}", e["rmse"]), (f"mss_{i}", e["mss"]),
                    (f"r2_{i}", e["r2"])]
    S = fit.resid.T @ fit.resid / n
    ll = -0.5 * n * (m * (1 + np.log(2 * np.pi)) + float(np.linalg.slogdet(S)[1]))
    scalars += [("k", float(len(b_all))), ("N", n), ("rank", np.linalg.matrix_rank(fit.V)), ("ll", ll),
                ("dfk2_adj", float(n)), ("k_eq", m), ("ic", fit.iterations if est.iterated else 0)]
    macros = [("cmd", cmd), ("cmdline", cmdline(cmd, args)), ("method", method),
              ("marginsdefault", " ".join(f"predict(xb equation({b.name}))" for b in blocks)),
              ("marginsok", "XB default"), ("marginsnotok", "stdp Residuals stddp"), ("predict", "reg3_p"),
              ("corr", "unstructured"), ("exog", " ".join(exog)), ("depvar", est.depvar),
              ("eqnames", " ".join(b.name for b in blocks)), ("properties", "b V")]
    from ..commands.matrix import Matrix
    mats = []
    if cons is not None:
        mats.append(("Cns", Matrix(np.column_stack([cons[0], cons[1]]), [str(i + 1) for i in range(len(cons[1]))],
                                   list(names) + ["_r"])))
    sig_names = [b.name for b in blocks]
    mats.append(("Sigma", Matrix(fit.sigma, sig_names, sig_names)))
    post(s, est, scalars, macros, smp.mask, extra_mats=mats)
    _display_system(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _display_system(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    out = s.output
    lev = level if level is not None else est.level
    if header:
        out.write(f"\n{est.title}\n", "text")
        bar = "-" * 74
        out.write(bar + "\n", "text")
        statcol = "F-Stat" if est.stat == "t" else "chi2"
        out.write(f'Equation             Obs   Parms        RMSE    "R-sq"{statcol:>11}{"P":>9}\n', "text")
        out.write(bar + "\n", "text")
        for nm, e in zip(est.names_eq, est.eq_info):
            if est.stat == "t":
                v = e["chi2"] / e["df"] if e["chi2"] < SYS else SYS
                p = st.f.sf(v, e["df"], est.df_r) if v < SYS else SYS
            else:
                v = e["chi2"]
                p = chi2_tail(v, e["df"]) if v < SYS else SYS
            out.write(f"{nm:<20}{int(e['obs']):>4}{int(e['parms']):>8}{g(e['rmse'], '%9.0g'):>12}"
                      f"{g(e['r2'], '%10.4f'):>10}{g(v, '%11.2f'):>11}{g(p, '%9.4f'):>9}\n", "text")
        out.write(bar + "\n", "text")
        if est.corr_opt:
            _corr_block(s, est)
        out.write("\n", "text")
        for i, t in enumerate(getattr(est, "cns_text", []), 1):      # VERIFY the place of the list
            out.write(f" ({i:>2})  {t}\n", "text")
        if getattr(est, "cns_text", []):
            out.write("\n", "text")
    if table:
        coef_table(s, "", est.rows, stat=est.stat, df=est.df_r, level=lev)
        if est.foot:
            out.write(f"Endogenous variables:  {' '.join(est.foot[0])} \n", "text")
            out.write(f"Exogenous variables:   {' '.join(est.foot[1])} \n", "text")
            out.write("-" * 78 + "\n", "text")


def _corr_block(s, est) -> None:
    out = s.output
    c = est.fit.corr
    nm = est.names_eq
    out.write("\nCorrelation matrix of residuals:\n\n", "text")
    out.write(" " * 12 + "".join(f"{x[:8]:>9}" for x in nm) + "\n", "text")
    for i, x in enumerate(nm):
        out.write(f"{x[:12]:>12}" + "".join(f"{c[i, j]:>9.4f}" for j in range(i + 1)) + "\n", "text")
    out.write(f"\nBreusch-Pagan test of independence: chi2({est.fit.df_bp}) = {est.fit.chi2_bp:>9.3f}, "
              f"Pr = {chi2_tail(est.fit.chi2_bp, est.fit.df_bp):.4f}\n", "text")   # VERIFY


# ---------------------------------------------------------------------------
# rreg
# ---------------------------------------------------------------------------

def cmd_rreg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "rreg", args)
        return
    from .regress import Sample
    smp = Sample(s, args, opts_spec={"tune": 2, "tolerance": 3, "iterate": 2, "genwt": 3, "nolog": 4,
                                     "graph": 5, "noheader": 6, "notable": 5}, weights=())
    o = smp.o
    ds = s.data
    est_cols = [c for c in smp.cols if c.values is not None]
    notes = fvars.drop_collinear(smp.cols, smp.constant, np.ones(smp.n))
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(smp.n)] if smp.constant else []))
    tune = float(o["tune"]) if o.get("tune") else 7.0
    tol = float(o["tolerance"]) if o.get("tolerance") else 0.01
    b, V, wts, keep, log, scale = syseq.rreg(X, smp.y, tune=tune, tol=tol)
    out = s.output
    ndrop = int((~keep).sum())
    if ndrop:
        out.write(f"\n   Cook's D > 1 dropped {ndrop} observations\n", "text")   # VERIFY
    if not o.get("nolog"):
        out.write("\n", "text")
        for kind, it, diff in log:
            out.write(f"{kind:>8} iteration {it}:  maximum difference in weights = {g(diff, '%9.0g').strip()}\n",
                      "text")
    if o.get("genwt"):
        from ..core.dataset import Variable
        name = str(o["genwt"]).strip()
        if ds.has(name):
            raise StataError(110, f"{name} already defined")
        var = Variable(name, "double", np.full(ds.nobs, SYS))
        full = np.full(ds.nobs, SYS)
        full[np.flatnonzero(smp.mask)] = wts
        var.data = full
        ds.add(var)
    N = int(keep.sum())
    k = X.shape[1]
    slopes = list(range(len(est_cols)))
    df_m = len(slopes)
    df_r = N - k
    F = float(b[slopes] @ np.linalg.solve(V[np.ix_(slopes, slopes)], b[slopes])) / df_m if slopes else SYS
    kk = len(smp.cols) + (1 if smp.constant else 0)
    names_full = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols] \
        + (["_cons"] if smp.constant else [])
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([kk - 1] if smp.constant else [])
    b_full, V_full = np.zeros(kk), np.zeros((kk, kk))
    b_full[pos] = b
    V_full[np.ix_(pos, pos)] = V
    bmap = {c.name: b_full[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(V_full[i, i], 0))) for i, c in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b_full[-1], float(np.sqrt(max(V_full[-1, -1], 0)))
    est = Estimates("rreg", smp.depname, names_full, b_full, V_full, N, "t", df_r)
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level, est.vcetype, est.clustvar, est.n_clust = smp.level, "", "", None
    est.extra = {"N": N, "df_m": df_m, "df_r": df_r, "F": F}
    est.display = _display_rreg
    scalars = [("N", N), ("df_m", df_m), ("df_r", df_r), ("F", F), ("rmse", scale)]
    macros = [("cmdline", cmdline("rreg", args)), ("title", "Robust regression"), ("depvar", smp.depname),
              ("cmd", "rreg"), ("properties", "b V"), ("predict", "regres_p")]
    post(s, est, scalars, macros, smp.mask)
    _display_rreg(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _display_rreg(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        right = [("Number of obs", comma(x["N"])),
                 (_flabel(x["df_m"], x["df_r"], 3, 10), g(x["F"], "%10.2f")),
                 ("Prob > F", g(st.f.sf(x["F"], x["df_m"], x["df_r"]) if x["F"] < SYS else SYS, "%10.4f"))]
        _header(s, ["Robust regression"], right)
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=lev)
