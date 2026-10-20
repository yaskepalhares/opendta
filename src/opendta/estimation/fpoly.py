"""Fractional polynomials: the fp and mfp prefixes ([R] fp, [R] mfp).

fp fits the command with every first- and second-degree fractional polynomial of the variable
(powers -2 -1 -.5 0 .5 1 2 3 by default), compares them by deviance and leaves the best model of
the requested dimension. mfp runs the multivariable selection (Royston and Altman's closed test
procedure) over the variables of the model.

Stata's scaling rule, the layout of the comparison tables and the order of the variables in mfp
are only partly documented: the marked lines have the VERIFY tag (compat/do/1903)."""

from __future__ import annotations

import itertools
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options
from ..lazy import lazy
from .results import g9

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_POWERS = (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0)


def fp_terms(x: np.ndarray, powers: tuple[float, ...]) -> np.ndarray:
    """Columns of the fractional polynomial of x (> 0) with the given powers (sorted):
    repeated powers multiply the previous term by ln x."""
    cols = []
    lx = np.log(x)
    for j, p in enumerate(sorted(powers)):
        if j and p == sorted(powers)[j - 1]:
            cols.append(cols[-1] * lx)
        else:
            cols.append(lx if p == 0 else x ** p)
    return np.column_stack(cols)


def fp_sets(degree: int, powers: tuple[float, ...]) -> list[tuple[float, ...]]:
    return list(itertools.combinations_with_replacement(sorted(powers), degree))


def scaling(x: np.ndarray) -> tuple[float, float]:
    """Shift a and scale b so that (x + a) / b is positive and of moderate size (VERIFY: Stata's
    rule is only described in words)."""
    u = np.unique(x)
    a = 0.0
    if u[0] <= 0:
        gaps = np.diff(u)
        delta = float(gaps[gaps > 0].min()) if np.any(gaps > 0) else 1.0
        a = delta - float(u[0])
    rng = float(u[-1] + a)
    b = 10.0 ** int(np.trunc(np.log10(rng))) if rng > 0 else 1.0
    return a, b


def _deviance(s: "Session"):
    est = getattr(s, "_estimates", None)
    ll = s.e.get("ll")
    if est is None or not isinstance(ll, float) or ll >= SYS:
        raise StataError(322, "the command does not store the log likelihood e(ll)")     # VERIFY
    return -2.0 * ll, s.e.get("rmse")


class _Fitter:
    def __init__(self, s: "Session", command: str, var: str, stub: str, a: float, b: float, center: float):
        self.s, self.command, self.var, self.stub = s, command, var, stub
        self.a, self.b, self.center = a, b, center
        self.x = s.data.get(var).data
        self.ok = self.x < SYS
        self.cache: dict = {}

    def names(self, k: int) -> list[str]:
        return [f"{self.stub}_{j + 1}" for j in range(k)]

    def write(self, powers: tuple[float, ...]) -> list[str]:
        from ..core.dataset import Variable
        ds = self.s.data
        xs = (self.x[self.ok] + self.a) / self.b
        cols = fp_terms(xs, powers)
        if self.center:
            cols = cols - self.center_values(powers)
        names = self.names(len(powers))
        for j, nm in enumerate(names):
            full = np.full(ds.nobs, SYS)
            full[self.ok] = cols[:, j]
            if ds.has(nm):
                ds.get(nm).data = full
            else:
                ds.add(Variable(nm, "double", full))
        return names

    def center_values(self, powers):
        xs = (self.x[self.ok] + self.a) / self.b
        return fp_terms(xs, powers).mean(axis=0)

    def run(self, replacement: str, noisy: bool = False):
        text = self.command.replace(f"<{self.var}>", replacement)
        out = self.s.output
        if not noisy:
            out.quiet_depth += 1
        try:
            self.s.interp._execute_expanded(text)
        finally:
            if not noisy:
                out.quiet_depth -= 1
        return _deviance(self.s)

    def fit(self, powers, noisy: bool = False):
        if powers == "omit":
            return self.run("", noisy)
        if powers == "linear":
            return self.run(self.var, noisy)
        names = self.write(powers)
        return self.run(" ".join(names), noisy)

    def drop(self) -> None:
        ds = self.s.data
        ds.drop_vars([n for n in self.names(4) if ds.has(n)])


def best_fp(fitter: _Fitter, degree: int, powers: tuple[float, ...], progress=None):
    best = None
    for ps in fp_sets(degree, powers):
        dev, rmse = fitter.fit(ps)
        if progress:
            progress()
        if np.isfinite(dev) and (best is None or dev < best[0] - 1e-10):
            best = (dev, rmse, ps)
    return best


def _ftest(s: "Session", dev0: float, dev1: float, df: int, n: int, rmse_big, dfr: float | None):
    """p-value for the deviance difference dev0 - dev1 (df extra parameters); an F test with the
    residual degrees of freedom when the model is a linear regression, a chi2 otherwise."""
    diff = max(dev0 - dev1, 0.0)
    if dfr:
        # regress: deviance = n (1 + ln 2pi + ln(rss/n)); F from the rss ratio
        rss1 = np.exp(dev1 / n - 1 - np.log(2 * np.pi)) * n
        rss0 = np.exp(dev0 / n - 1 - np.log(2 * np.pi)) * n
        f = ((rss0 - rss1) / df) / (rss1 / dfr)
        return diff, float(st.f.sf(f, df, dfr))
    return diff, float(chi2_tail(diff, df))


def _split(args: str):
    colon = find_top(args, ":")
    if colon == -1:
        raise StataError(198, "invalid syntax")
    head, command = args[:colon].strip(), args[colon + 1:].strip()
    c = find_top(head, ",")
    return (head, "") if c == -1 else (head[:c].strip(), head[c + 1:])


_OPT = {"df": 2, "powers": 2, "dimension": 3, "scale": 2, "noscale": 3, "center": 3, "nocenter": 3,
        "replace": 3, "all": 2, "zero": 3, "catzero": 4, "alpha": 3, "select": 3, "cycles": 3,
        "level": 3, "classic": 3}


def _powers(o: dict) -> tuple[float, ...]:
    if isinstance(o.get("powers"), str):
        return tuple(sorted(float(v) for v in o["powers"].replace(",", " ").split()))
    return _POWERS


def cmd_fp(s: "Session", args: str) -> None:
    command = args[find_top(args, ":") + 1:].strip() if find_top(args, ":") != -1 else ""
    head, opt = _split(args)
    m = re.fullmatch(r"<(\w+)>", head.strip())
    if not m:
        raise StataError(198, "variable must be specified in angle brackets")
    var = s.expand_varlist(m.group(1))[0]
    o = match_options(opt, _OPT) if opt.strip() else {}
    if f"<{var}>" not in command:
        raise StataError(198, f"<{var}> must appear in the command")
    powers = _powers(o)
    degree = int(o["dimension"]) if isinstance(o.get("dimension"), str) else 2
    if isinstance(o.get("df"), str):
        degree = max(int(float(o["df"])) // 2, 1)
    ds = s.data
    x = ds.get(var).data
    ok = x < SYS
    a = b = 0.0
    a, b = (0.0, 1.0)
    if o.get("scale"):
        a, b = scaling(x[ok])
        if isinstance(o["scale"], str):
            a, b = (float(v) for v in o["scale"].split())
    elif x[ok].min() <= 0 and not o.get("zero"):
        raise StataError(411, f"{var} has nonpositive values; use the scale option")      # VERIFY
    stub = var
    for j in (1, 2):
        if ds.has(f"{stub}_{j}") and not o.get("replace"):
            raise StataError(110, f"variable {stub}_{j} already defined")    # VERIFY
    fitter = _Fitter(s, command, var, stub, a, b, 0.0)
    out = s.output
    total = sum(len(fp_sets(d, powers)) for d in range(1, degree + 1))
    out.write(f"(fitting {total} models)\n", "text")
    done = [0]
    marks = iter(range(10, 101, 10))
    nxt = [next(marks)]
    out.write("(", "text")

    def progress():
        done[0] += 1
        while nxt[0] is not None and done[0] * 100 >= nxt[0] * total:
            out.write("....", "text")
            out.write(f"{nxt[0]}%", "text")
            nxt[0] = next(marks, None)
    models = {}
    models["omit"] = fitter.fit("omit")
    models["linear"] = fitter.fit("linear")
    for d in range(1, degree + 1):
        models[d] = best_fp(fitter, d, powers, progress)
    out.write(")\n", "text")
    n = int(s.e.get("N", 0))
    # N and the residual df of the best largest model
    big = models[degree]
    fitter.fit(big[2])
    dfr = s.e.get("df_r")
    dfr = float(dfr) if isinstance(dfr, float) and dfr < SYS else None
    has_rmse = big[1] is not None and isinstance(big[1], float)
    _compare_table(s, var, models, degree, n, dfr, has_rmse)
    # final model: the best of the requested dimension
    s.macros.globals.pop("_fp_tmp", None)
    names = fitter.write(big[2])
    _label(s, names, var, big[2], a, b)
    fitter.run(" ".join(names), noisy=True)
    s.r = {}


def _label(s: "Session", names: list[str], var: str, powers, a: float, b: float) -> None:
    base = var if (a == 0 and b == 1) else f"({var}{'+' + format(a, 'g') if a else ''})/{b:g}"
    for j, (nm, p) in enumerate(zip(names, sorted(powers))):
        if j and p == sorted(powers)[j - 1]:
            lab = s.data.get(names[j - 1]).label + "*ln(X)"
        else:
            lab = "ln(X)" if p == 0 else f"X^{p:g}"
        s.data.get(nm).label = lab.replace("X", base)


def _compare_table(s: "Session", var: str, models: dict, degree: int, n: int, dfr, has_rmse: bool) -> None:
    out = s.output
    big = models[degree]
    rows = [("Omitted", 0, models["omit"], "."), ("Linear", 1, models["linear"], "1")]
    for d in range(1, degree + 1):
        rows.append((f"m = {d}", 2 * d, models[d], " ".join(f"{p:g}" for p in models[d][2])))
    width = 71 if has_rmse else 60
    out.write("\nFractional polynomial comparisons:\n", "text")
    out.write("-" * width + "\n", "text")
    head = f"{var[:12]:>12} |   df    Deviance" + ("   Res. s.d." if has_rmse else "") + \
        "   Dev. dif.   P(*)    Powers"
    out.write(head + "\n", "text")
    out.write("-" * 13 + "+" + "-" * (width - 13) + "\n", "text")
    for lab, df, (dev, rmse, *rest), pw in rows:
        out.write(f"{lab:<12} |{df:>5}{dev:>12.3f}", "result")
        if has_rmse:
            out.write(f"{rmse:>12.5g}" if isinstance(rmse, float) else f"{'':>12}", "result")
        dd = dev - big[0]
        if lab.startswith("m =") and df == 2 * degree:
            p = "     "
            out.write(f"{dd:>12.3f}{'':>7}  {pw}\n", "result")
        else:
            ddf = 2 * degree - df
            diff, p = _ftest(s, dev, big[0], ddf, n, None, dfr)
            out.write(f"{dd:>12.3f}{p:>7.3f}  {pw}\n", "result")
    out.write("-" * width + "\n", "text")
    out.write(f"(*) P-value from deviance difference comparing reported model with m = {degree} model\n", "text")


# ---------------------------------------------------------------------------
# mfp
# ---------------------------------------------------------------------------

def cmd_mfp(s: "Session", args: str) -> None:
    colon = find_top(args, ":")
    if colon == -1:
        raise StataError(198, "invalid syntax")
    head, command = args[:colon].strip(), args[colon + 1:].strip()
    o = match_options(head.lstrip(", "), _OPT) if head.strip() else {}
    alpha = float(o["alpha"]) if isinstance(o.get("alpha"), str) else 0.05
    select = float(o["select"]) if isinstance(o.get("select"), str) else 1.0
    powers = _powers(o)
    dfdef = int(float(o["df"])) if isinstance(o.get("df"), str) else 4
    cycles = int(o["cycles"]) if isinstance(o.get("cycles"), str) else 5
    m = re.match(r"(\S+)\s+(\S+)\s+(.*)$", command, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    cmd, dep, rest = m.groups()
    cm = re.match(r"([^,\[]*?)(\s+if\s.*|\s+in\s.*|\s*\[.*|\s*,.*)?$", rest, re.S)
    varlist, tail = (cm.group(1) if cm else rest), (cm.group(2) or "") if cm else ""
    vars_ = s.expand_varlist(varlist)
    ds = s.data
    out = s.output

    def info(v):
        x = ds.get(v).data
        ok = x < SYS
        nu = len(np.unique(x[ok]))
        return nu

    # model state per variable: ("fp", powers) | ("lin",) | ("out",)
    state = {}
    dfs = {}
    for v in vars_:
        nu = info(v)
        dfs[v] = dfdef if nu >= 4 else 1 if nu > 1 else 1
        if nu >= 5 and dfs[v] >= 4:
            dfs[v] = dfdef
        elif nu < 5:
            dfs[v] = 1
        state[v] = ("lin",)
    scal = {}
    for v in vars_:
        x = ds.get(v).data
        ok = x < SYS
        scal[v] = scaling(x[ok]) if x[ok].min() <= 0 or x[ok].max() > 100 else (0.0, 1.0)

    def text(st_: dict, tmp: dict) -> str:
        parts = []
        for v in vars_:
            if st_[v][0] == "lin":
                parts.append(v)
            elif st_[v][0] == "fp":
                parts += tmp[v]
        return f"{cmd} {dep} {' '.join(parts)}{tail}"

    created: dict[str, list[str]] = {}

    def write_fp(v, powers_):
        from ..core.dataset import Variable
        x = ds.get(v).data
        ok = x < SYS
        a, b = scal[v]
        cols = fp_terms((x[ok] + a) / b, powers_)
        names = [f"I{v}__{j + 1}" for j in range(len(powers_))]
        for j, nm in enumerate(names):
            full = np.full(ds.nobs, SYS)
            full[ok] = cols[:, j]
            if ds.has(nm):
                ds.get(nm).data = full
            else:
                ds.add(Variable(nm, "double", full))
        return names

    def fit(st_: dict):
        tmp = {v: write_fp(v, st_[v][1]) for v in vars_ if st_[v][0] == "fp"}
        for v in tmp:
            created[v] = tmp[v]
        out.quiet_depth += 1
        try:
            s.interp._execute_expanded(text(st_, tmp))
        finally:
            out.quiet_depth -= 1
        dev, _ = _deviance(s)
        dfr = s.e.get("df_r")
        return dev, int(s.e.get("N", 0)), (float(dfr) if isinstance(dfr, float) and dfr < SYS else None)

    try:
        dev_lin, n, dfr = fit(state)
        out.write(f"Deviance for model with all terms untransformed = {dev_lin:.3f}, {n} observations\n", "text")
        # order: increasing p-value of the Wald test of each variable in the linear model
        est = s._estimates
        pv = {}
        for v in vars_:
            idx = [i for i, nm in enumerate(est.names) if nm.split(":")[-1] == v]
            if idx:
                b_, V_ = np.asarray(est.b)[idx], np.asarray(est.V)[np.ix_(idx, idx)]
                chi = float(b_ @ np.linalg.solve(V_, b_))
                pv[v] = float(chi2_tail(chi, len(idx)))
            else:
                pv[v] = 1.0
        order = sorted(vars_, key=lambda v: pv[v])
        history = []
        for cyc in range(1, cycles + 1):
            before = dict(state)
            for v in order:
                state = _select_one(state, v, dfs[v], powers, alpha, select, fit)
            if state == before:
                break
        out.write(f"\nFractional polynomial fitting algorithm converged after {cyc} cycles.\n", "text")
        out.write("\nTransformations of covariates:\n\n", "text")
        for v in vars_:
            if state[v][0] == "fp":
                a, b = scal[v]
                out.write(f"-> {v}: powers {' '.join(f'{p:g}' for p in state[v][1])}"
                          f" (shift {a:g}, scale {b:g})\n", "text")   # VERIFY: layout
        out.write(f"\nFinal multivariable fractional polynomial model for {dep}\n", "text")
        out.write("-" * 76 + "\n", "text")
        out.write("    Variable    |-----Initial-----   |----------Final----------\n", "text")
        out.write("                |   df     Select   Alpha    Status    df    Powers\n", "text")
        out.write("-" * 76 + "\n", "text")
        for v in vars_:
            stt = state[v]
            status = "out" if stt[0] == "out" else "in"
            fdf = 0 if stt[0] == "out" else 1 if stt[0] == "lin" else 2 * len(stt[1])
            pw = "" if stt[0] == "out" else "1" if stt[0] == "lin" else " ".join(f"{p:g}" for p in stt[1])
            out.write(f"{v:>15} |{dfs[v]:>5}{select:>11.4f}{alpha:>8.4f}{status:>10}{fdf:>6}  {pw}\n", "result")
        out.write("-" * 76 + "\n", "text")
        tmp = {v: write_fp(v, state[v][1]) for v in vars_ if state[v][0] == "fp"}
        for v, nm in tmp.items():
            created[v] = nm
        s.interp._execute_expanded(text(state, tmp))
    finally:
        pass


def _select_one(state: dict, v: str, df: int, powers, alpha: float, select: float, fit) -> dict:
    """Closed test procedure for one variable with the others fixed."""
    n_deg = max(df // 2, 1) if df > 1 else 0
    cands = {}
    cur = dict(state)
    cur[v] = ("out",)
    dev_out, n, dfr = fit(cur)
    cur[v] = ("lin",)
    dev_lin, _, _ = fit(cur)
    if n_deg == 0:
        diff, p = _ftest(None, dev_out, dev_lin, 1, n, None, dfr)
        keep = select >= 1 or p < select
        new = dict(state)
        new[v] = ("lin",) if keep else ("out",)
        return new
    best_by_deg = {}
    for d in range(1, n_deg + 1):
        bestd = None
        for ps in fp_sets(d, powers):
            cur[v] = ("fp", ps)
            dev, _, _ = fit(cur)
            if bestd is None or dev < bestd[0]:
                bestd = (dev, ps)
        best_by_deg[d] = bestd
    dev_big, ps_big = best_by_deg[n_deg]
    new = dict(state)
    # 1) the variable in the model at all (4 df test for FP2): only when select < 1
    _, p_in = _ftest(None, dev_out, dev_big, 2 * n_deg, n, None, dfr)
    if select < 1 and p_in >= select:
        new[v] = ("out",)
        return new
    # 2) nonlinearity
    _, p_lin = _ftest(None, dev_lin, dev_big, 2 * n_deg - 1, n, None, dfr)
    if p_lin >= alpha:
        new[v] = ("lin",)
        return new
    # 3) degree: FP(m) against FP(m-1)
    chosen = n_deg
    for d in range(n_deg - 1, 0, -1):
        _, p_d = _ftest(None, best_by_deg[d][0], dev_big, 2 * (n_deg - d), n, None, dfr)
        if p_d >= alpha:
            chosen = d
        else:
            break
    new[v] = ("fp", best_by_deg[chosen][1])
    return new
