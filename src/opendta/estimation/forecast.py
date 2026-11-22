"""forecast ([TS] forecast): models made of stored estimation results, identities and exogenous
variables, solved period by period (Gauss-Seidel within each period) into forecast variables.

The layouts and the default begin() follow compat 2008. Estimation equations use the linear
prediction of single-equation models; stochastic simulation (simulate()) is not available yet."""

from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Variable
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lang.words import strip_outer_quotes

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_SUBS = [("create", 2), ("estimates", 3), ("identity", 2), ("coefvector", 2), ("exogenous", 2),
         ("solve", 1), ("adjust", 2), ("describe", 1), ("list", 1), ("clear", 2), ("drop", 2),
         ("query", 1)]
_SUB_TEXT = ("forecast subcommands are create, estimates, identity, coefvector, exogenous, solve, "
             "adjust, describe, list, clear, drop, and query.")


@dataclass
class Model:
    name: str
    equations: list = field(default_factory=list)      # ("estimates", depvar, Estimates) | ("identity", var, exp)
    exog: list[str] = field(default_factory=list)
    adjust: dict[str, list[str]] = field(default_factory=dict)
    commands: list[str] = field(default_factory=list)
    solved: bool = False
    created: list[str] = field(default_factory=list)   # variables made by forecast solve

    def endogenous(self) -> list[str]:
        return [eq[1] for eq in self.equations]


def _model(s: "Session") -> Model | None:
    return getattr(s, "_forecast_model", None)


def _need(s: "Session", sub: str) -> Model:
    m = _model(s)
    if m is None:
        # compat 2008 (forecast describe); VERIFY: the other subcommands
        raise StataError(111, "forecast problem not found\n    Before you can issue forecast "
                              f"{sub} you must begin the forecast problem by issuing forecast "
                              "create.")
    return m


def cmd_forecast(s: "Session", args: str) -> None:
    t = args.strip()
    if not t:
        tail = "" if _model(s) else ("  You have not started a forecast model yet, so you must begin "
                                     "with forecast create.")
        raise StataError(198, "missing forecast subcommand\n    " + _SUB_TEXT + tail)   # compat 2005
    head, _, rest = t.partition(" ")
    head_word = head.split(",")[0]
    sub = next((n for n, k in _SUBS if len(head_word) >= k and n.startswith(head_word)), None)
    if sub is None:
        raise StataError(198, f"invalid forecast subcommand {head_word}")   # VERIFY
    rest = (head[len(head_word):] + " " + rest).strip()
    {"create": _create, "estimates": _estimates, "identity": _identity, "coefvector": _coefvector,
     "exogenous": _exogenous, "solve": _solve, "adjust": _adjust, "describe": _describe,
     "list": _list, "clear": _clear, "drop": _drop, "query": _query}[sub](s, rest, t)


def _split(text: str) -> tuple[str, str]:
    c = find_top(text, ",")
    return (text, "") if c == -1 else (text[:c], text[c + 1:])


def _create(s, rest, full) -> None:
    body, opt = _split(rest)
    o = match_options(opt, {"replace": 7})
    if _model(s) is not None and not o.get("replace"):
        raise StataError(110, f"forecast model {_model(s).name} already exists; specify replace")  # VERIFY
    name = body.strip() or "forecast"
    s._forecast_model = Model(name)
    s._forecast_model.commands.append("forecast " + full)
    s.output.write(f"  Forecast model {name} started.\n", "text")   # compat 2008


def _estimates(s, rest, full) -> None:
    from .postest import _stored
    m = _need(s, "estimates")
    body, opt = _split(rest)
    o = match_options(opt, {"predict": 4, "names": 2, "advise": 3})
    if o.get("predict"):
        raise StataError(198, "option predict() is not available yet in OpenDTA")
    name = body.strip()
    snap = _stored(s).get(name)
    if snap is None:
        raise StataError(111, f"estimation result {name} not found")
    est = snap["est"]
    if est is None or getattr(est, "eqnames", None):
        raise StataError(198, "forecast estimates of multiple-equation models is not available yet in "
                              "OpenDTA")
    depvar = (o.get("names") and str(o["names"]).split(",")[0].strip()) or est.depvar
    if depvar in m.endogenous():
        raise StataError(110, f"{depvar} is already an endogenous variable in the model")   # VERIFY
    est = copy.deepcopy(est)
    est._stored_name = name   # forecast describe lists it
    m.equations.append(("estimates", depvar, est))
    m.commands.append("forecast " + full)
    s.output.write(f"  Added estimation results from {est.cmd}.\n"
                   f"  Forecast model {m.name} now contains {len(m.equations)} endogenous "
                   f"variable{'s' if len(m.equations) > 1 else ''}.\n", "text")


def _identity(s, rest, full) -> None:
    m = _need(s, "identity")
    body, opt = _split(rest)
    o = match_options(opt, {"generate": 3})
    var, eq, exp = body.partition("=")
    var = var.strip()
    words = var.split()
    if len(words) == 2 and words[0] in ("byte", "int", "long", "float", "double"):
        var = words[1]
    if not eq or not re.fullmatch(r"[A-Za-z_]\w*", var):
        raise StataError(198, "invalid syntax")
    if not s.data.has(var):
        if not o.get("generate"):
            raise StataError(111, f"variable {var} not found")
        # generate fills the variable with the expression (compat 2008: f_z holds y + x before
        # the forecast)
        from ..commands._util import eval_vector
        vec = np.broadcast_to(np.asarray(eval_vector(s, exp.strip()), dtype=float), (s.data.nobs,))
        ok = vec < SYS
        col = np.full(s.data.nobs, SYS)
        col[ok] = vec[ok].astype(np.float32)
        s.data.add(Variable(var, "float", col))
    m.equations.append(("identity", var, exp.strip()))
    m.commands.append("forecast " + full)
    s.output.write(f"  Forecast model {m.name} now contains {len(m.equations)} endogenous "
                   f"variable{'s' if len(m.equations) > 1 else ''}.\n", "text")


def _coefvector(s, rest, full) -> None:
    raise StataError(198, "forecast coefvector is not available yet in OpenDTA")


def _exogenous(s, rest, full) -> None:
    m = _need(s, "exogenous")
    names = s.expand_varlist(rest)
    if not names:
        raise StataError(100, "varlist required")
    for n in names:
        if n not in m.exog:
            m.exog.append(n)
    m.commands.append("forecast " + full)
    s.output.write(f"  Forecast model {m.name} now contains {len(m.exog)} declared exogenous "
                   f"variable{'s' if len(m.exog) > 1 else ''}.\n", "text")


def _adjust(s, rest, full) -> None:
    m = _need(s, "adjust")
    var, eq, exp = rest.partition("=")
    var = var.strip()
    if not eq or var not in m.endogenous():
        raise StataError(198, f"{var} is not an endogenous variable of the model")   # VERIFY
    m.adjust.setdefault(var, []).append(exp.strip())
    m.commands.append("forecast " + full)
    k = len(m.adjust[var])
    s.output.write(f"  Endogenous variable {var} now has {k} adjustment{'s' if k > 1 else ''}.\n",
                   "text")   # compat 2008


def _clear(s, rest, full) -> None:
    s._forecast_model = None


def _query(s, rest, full) -> None:
    m = _model(s)
    if m is None:
        s.output.write("  No forecast model exists.\n", "text")   # compat 2008
    else:
        s.output.write(f"  Forecast model {m.name} exists.\n", "text")


def _list(s, rest, full) -> None:
    m = _need(s, "list")
    for c in m.commands:
        s.output.write(c + "\n", "text")


def _describe(s, rest, full) -> None:
    m = _need(s, "describe")
    # layout of compat 2008. VERIFY: the wording of the empty and solved cases
    out = s.output
    rule = "-" * 65 + "\n"

    def center(x: str) -> str:
        pad = max(21 - len(x), 0)
        return " " * (pad // 2) + x + " " * (pad - pad // 2)

    def head(n: int, one: str, many: str) -> None:
        what = f"no {many}" if n == 0 else f"{n} {one if n == 1 else many}"
        out.write(f"\nForecast model {m.name} contains {what}{'.' if n == 0 else ':'}\n", "text")

    ests = [(v, o) for k, v, o in m.equations if k == "estimates"]
    ids = [v for k, v, _ in m.equations if k == "identity"]
    head(len(ests), "estimation result", "estimation results")
    if ests:
        out.write(rule + "       Estimation    |                     |\n"
                  "         result      |       Command       |   # LHS variables\n" + rule, "text")
        for i, (v, o) in enumerate(ests, 1):
            name = getattr(o, "_stored_name", "") or v
            out.write(f"  {f'{i}. {name}':<20}{center(o.cmd)}{1:>13}\n", "result")
        out.write(rule, "text")
    head(len(ids), "identity", "identities")
    if ids:
        out.write(rule + "       Endogenous variable\n" + rule, "text")
        for i, v in enumerate(ids, 1):
            out.write(f"  {i}. {v}\n", "result")
        out.write(rule, "text")
    head(0, "coefficient vector", "coefficient vectors")
    head(len(m.exog), "declared exogenous variable", "declared exogenous variables")
    if m.exog:
        out.write(rule, "text")
        for i, v in enumerate(m.exog, 1):
            out.write(f"  {i}. {v}\n", "result")
        out.write(rule, "text")
    head(len(m.equations), "endogenous variable", "endogenous variables")
    if m.equations:
        out.write(rule + "       Variable      |       Source        |    # adjustments\n" + rule,
                  "text")
        for i, (k, v, _) in enumerate(m.equations, 1):
            out.write(f"  {f'{i}. {v}':<20}{center(k)}{len(m.adjust.get(v, [])):>13}\n", "result")
        out.write(rule, "text")
    adj = [v for v in m.endogenous() if m.adjust.get(v)]
    head(len(adj), "variable with adjustments", "variables with adjustments")
    if adj:
        out.write(rule, "text")
        for i, v in enumerate(adj, 1):
            out.write(f"  {i}. {v}\n", "result")
        out.write(rule, "text")
    out.write(f"\nForecast model {m.name} has {'' if m.solved else 'not '}been solved.\n", "text")


def _drop(s, rest, full) -> None:
    m = _need(s, "drop")
    _, opt = _split(rest)
    o = match_options(opt, {"prefix": 3, "suffix": 3})
    names = [n for n in m.created if s.data.has(n)]
    if o.get("prefix") or o.get("suffix"):
        pre, suf = str(o.get("prefix") or ""), str(o.get("suffix") or "")
        names = [pre + v + suf for v in m.endogenous() if s.data.has(pre + v + suf)]
    if names:
        s.data.drop_vars(names)
    m.created = [n for n in m.created if n not in names]
    s.output.write(f"  ({len(names)} forecast variable{'s' if len(names) != 1 else ''} dropped)\n", "text")


# ---------------------------------------------------------------------------
# solve
# ---------------------------------------------------------------------------

def solve_period(t: int, model_eval, names: list[str], cols: dict[str, np.ndarray],
                 maxiter: int = 500, tol: float = 1e-9) -> int:
    """Gauss-Seidel at observation t: each equation in turn, until no value changes by more
    than tol (relative). Returns the number of iterations."""
    for it in range(1, maxiter + 1):
        delta = 0.0
        for k, v in enumerate(names):
            new = model_eval(k, t)
            old = cols[v][t]
            cols[v][t] = new
            if new < SYS and old < SYS:
                delta = max(delta, abs(new - old) / max(1.0, abs(old)))
            elif (new < SYS) != (old < SYS):
                delta = np.inf
        if delta <= tol:
            return it
    raise StataError(430, f"convergence not achieved in period {t + 1}")   # VERIFY


def _period_label(s, tvar: str, i: int) -> str:
    from ..core.formats import format_value
    v = s.data.get(tvar)
    fmt = v.fmt if v.fmt.startswith("%t") else "%9.0g"
    return format_value(float(v.data[i]), fmt, pad=False).strip()


def _time_index(s, tvar: str, text: str) -> int:
    val = float(s.eval(text))
    hits = np.flatnonzero(s.data.get(tvar).data == val)
    if not len(hits):
        raise StataError(198, f"{text} is not a period in the data")   # VERIFY
    return int(hits[0])


def _solve(s, rest, full) -> None:
    from .postest import xb_all
    from .tsops import ts_info
    m = _need(s, "solve")
    _, opt = _split(rest)
    o = match_options(opt, {"prefix": 3, "suffix": 3, "begin": 3, "end": 3, "periods": 3,
                            "static": 4, "actuals": 3, "log": 3, "simulate": 3, "double": 6,
                            "vtolerance": 4, "ztolerance": 4, "maxiter": 4})
    if o.get("simulate"):
        raise StataError(198, "option simulate() is not available yet in OpenDTA")
    if not m.equations:
        raise StataError(198, f"forecast model {m.name} has no endogenous variables")   # VERIFY
    info = ts_info(s.data)
    if info is None:
        raise StataError(111, "time variable not set, use tsset varname ...")
    if info.panel:
        raise StataError(198, "forecast with panel data is not available yet in OpenDTA")
    pre = str(o["prefix"]) if o.get("prefix") else ""
    suf = str(o["suffix"]) if o.get("suffix") else ""
    if not pre and not suf:
        pre = "f_"
    ds = s.data
    endo = m.endogenous()
    actual = {v: ds.get(v).data.astype(float).copy() for v in endo}
    n = ds.nobs
    if o.get("begin"):
        begin = _time_index(s, info.tvar, str(o["begin"]))
    else:
        # the period after the last one in which every endogenous variable that has data is
        # nonmissing (compat 2008: 26). VERIFY: models whose variables end at different times
        full = np.ones(n, dtype=bool)
        for v in endo:
            if (actual[v] < SYS).any():
                full &= actual[v] < SYS
        last = np.flatnonzero(full)
        begin = min(int(last[-1]) + 1, n - 1) if len(last) else 0
    end = _time_index(s, info.tvar, str(o["end"])) if o.get("end") else n - 1
    if o.get("periods"):
        end = min(begin + int(o["periods"]) - 1, n - 1)
    if end < begin:
        raise StataError(198, "end() must not be before begin()")   # VERIFY
    names = [pre + v + suf for v in endo]
    for nm in names:
        if ds.has(nm):
            raise StataError(110, f"variable {nm} already defined")
    cols = {v: np.where(np.arange(n) < begin, actual[v], SYS) for v in endo}
    # the forecast variables hold the actual values before the first period (compat 2008)
    static, use_actuals = bool(o.get("static")), bool(o.get("actuals"))
    work = ds.copy()

    def view(t: int) -> None:
        for v in endo:
            col = cols[v].copy()
            if static:
                col[:t] = actual[v][:t]
            work.get(v).data = col

    def model_eval(k: int, t: int) -> float:
        kind, var, obj = m.equations[k]
        if use_actuals and actual[var][t] < SYS:
            return float(actual[var][t])
        view(t)
        saved = s.data
        s.data = work
        try:
            if kind == "estimates":
                xb, ok = xb_all(s, obj)
                val = float(xb[t]) if ok[t] else SYS
            else:
                from ..commands._util import eval_vector
                vec = np.broadcast_to(np.asarray(eval_vector(s, obj), dtype=float), (n,))
                val = float(vec[t])
            for exp in m.adjust.get(var, []):
                work.get(var).data = np.where(np.arange(n) == t, val, work.get(var).data)
                from ..commands._util import eval_vector
                vec = np.broadcast_to(np.asarray(eval_vector(s, exp), dtype=float), (n,))
                val = float(vec[t])
        finally:
            s.data = saved
        return val if val < SYS else SYS

    out = s.output
    kind = "static" if static else "dynamic"
    title = f"Computing {kind} forecasts for model {m.name}."
    out.write(f"\n{title}\n" + "-" * len(title) + "\n", "text")
    out.write(f"Starting period:  {_period_label(s, info.tvar, begin)}\n", "text")
    out.write(f"Ending period:    {_period_label(s, info.tvar, end)}\n", "text")
    out.write(f"Forecast {'prefix' if pre else 'suffix'}:  {pre or suf}\n\n", "text")
    show = str(o.get("log") or "brief") not in ("off", "quietly")
    for t in range(begin, end + 1):
        if show:
            out.write(f"{_period_label(s, info.tvar, t)}: ", "text")
        its = solve_period(t, model_eval, endo, cols)
        bad = [v for v in endo if not cols[v][t] < SYS]
        if bad:
            # compat 2008 (a static forecast reaching missing lagged values)
            if show:
                out.write("\n", "text")
            listed = bad[0] if len(bad) == 1 else (" and ".join(bad) if len(bad) == 2 else
                                                     ", ".join(bad[:-1]) + ", and " + bad[-1])
            raise StataError(416, "missing values encountered\n    Missing values were encountered "
                                  "while attempting to solve the model at time "
                                  f"{_period_label(s, info.tvar, t)}.  "
                                  f"Variable{'s' if len(bad) > 1 else ''} {listed} evaluate"
                                  f"{'' if len(bad) > 1 else 's'} to missing.")
        if show:
            # VERIFY: Stata's solver takes more iterations (one dot each) than this one
            out.write(" " + "." * its + "\n", "text")
    vtype = "double" if o.get("double") else "float"
    m.solved = True
    for v, nm in zip(endo, names):
        col = np.where(np.arange(n) <= end, cols[v], SYS)
        ds.add(Variable(nm, vtype, col))
        m.created.append(nm)
    foot = (f"Forecast {len(endo)} variable{'s' if len(endo) != 1 else ''} spanning "
            f"{end - begin + 1} period{'s' if end > begin else ''}.")
    out.write(f"\n{foot}\n" + "-" * len(foot) + "\n", "text")
    s.notify_state()
