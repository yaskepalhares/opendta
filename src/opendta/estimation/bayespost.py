"""Postestimation after bayesmh and the bayes prefix ([BAYES] bayesstats, bayestest, bayesgraph):
summaries, information criteria and interval tests computed from the saved draws."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import mcmc
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lazy import lazy
from .bayes import _BRACE, _commas, _current, _srow
from .results import g9

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")


def _eval_draws(est, text: str) -> np.ndarray:
    """A function of the parameters, evaluated draw by draw: {name} refer to parameters."""
    names = list(est.extra["pnames"])
    draws = est.extra["draws"]
    py = text
    env = {"exp": np.exp, "ln": np.log, "log": np.log, "log10": np.log10, "sqrt": np.sqrt,
           "abs": np.abs, "invlogit": lambda v: 1 / (1 + np.exp(-v)), "logit": lambda v: np.log(v / (1 - v)),
           "normal": lambda v: special.ndtr(v), "invnormal": lambda v: special.ndtri(v)}

    def sub(m):
        body = m.group(1).strip()
        nm = body if body in names else None
        if nm is None:
            raise StataError(198, f"parameter {{{body}}} not found")   # VERIFY
        env[f"_p{names.index(nm)}"] = draws[:, names.index(nm)]
        return f"_p{names.index(nm)}"
    py = _BRACE.sub(sub, py).replace("^", "**")
    try:
        return np.asarray(eval(py, {"__builtins__": {}}, env), float) * np.ones(draws.shape[0])
    except Exception:
        raise StataError(198, f"invalid expression {text}")   # VERIFY


def _param_items(est, text: str):
    """[(label, draws, index or None)] for a list of {parameters} and (expressions)."""
    names = list(est.extra["pnames"])
    draws = est.extra["draws"]
    text = text.strip()
    if not text:
        return [(n, draws[:, i], i) for i, n in enumerate(names)]
    items = []
    i = 0
    while i < len(text):
        ch = text[i]
        if ch.isspace():
            i += 1
        elif ch == "(":
            depth, j = 1, i + 1
            while j < len(text) and depth:
                depth += {"(": 1, ")": -1}.get(text[j], 0)
                j += 1
            e = text[i + 1:j - 1]
            lm = re.match(r"^\s*(\w+)\s*:(?!\w*\})", e)
            lab, body = (lm.group(1), e[lm.end():]) if lm and "{" not in e[:lm.end()] else (e.strip(), e)
            items.append((lab, _eval_draws(est, body), None))
            i = j
        elif ch == "{":
            j = text.index("}", i)
            body = text[i + 1:j]
            i = j + 1
            if ":" in body:
                label, rest = body.split(":", 1)
                label = label.strip()
                sel = [k for k, n in enumerate(names) if n.startswith(label + ":")]
                if rest.strip():
                    want = set(rest.replace(",", " ").split())
                    sel = [k for k in sel if names[k].split(":", 1)[1] in want]
                for k in sel:
                    items.append((names[k], draws[:, k], k))
            else:
                for nm in body.split():
                    if nm not in names:
                        raise StataError(198, f"parameter {{{nm}}} not found")   # VERIFY
                    k = names.index(nm)
                    items.append((nm, draws[:, k], k))
        else:
            raise StataError(198, f"invalid parameter list: {text}")   # VERIFY
    if not items:
        raise StataError(198, "no parameters selected")   # VERIFY
    return items


def cmd_bayesstats(s: "Session", args: str) -> None:
    est = _current(s)
    m = re.match(r"\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "bayesstats requires summary, ess or ic")   # VERIFY
    sub = m.group(1).lower()
    rest = m.group(2)
    subs = [w for w in ("summary", "ess", "ic") if w.startswith(sub)]
    if len(subs) != 1:
        raise StataError(198, f"bayesstats {sub}: invalid subcommand")   # VERIFY
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    o = match_options(opt, {"hpd": 3, "clevel": 4, "skip": 4, "nolog": 3, "all": 3}) if opt.strip() else {}
    level = float(o["clevel"]) if o.get("clevel") else est.extra["level"]
    if subs[0] == "ic":
        _stats_ic(s, est, body, o)
        return
    items = _param_items(est, body)
    draws = np.column_stack([d for _l, d, _i in items])
    if subs[0] == "summary":
        summ = mcmc.summarize(draws, level, bool(o.get("hpd")))
        s.output.write("\nPosterior summary statistics" + " " * 22 + f"MCMC sample size = {_commas(draws.shape[0]):>9}\n \n", "text")
        _flat_table(s, [l for l, _d, _i in items], summ, level, bool(o.get("hpd")))
    else:
        s.output.write("\nEfficiency summaries    " + f"MCMC sample size = {_commas(draws.shape[0]):>9}\n \n", "text")
        out = s.output
        W = 12
        out.write("-" * (W + 40) + "\n", "text")
        out.write(f"{'':>{W}} |        ESS   Corr. time    Efficiency\n", "text")
        out.write("-" * (W + 1) + "+" + "-" * 38 + "\n", "text")

        def one(lab, d):
            e = mcmc.ess(d)
            out.write(f"{lab:>{W}} |", "text")
            out.write(f"{e:>11.2f}{len(d) / e:>13.2f}{e / len(d):>14.4f}\n", "result")
        _grouped(out, [l for l, _d, _i in items], [d for _l, d, _i in items], one, W, 38)
        out.write("-" * (W + 40) + "\n", "text")


def _grouped(out, labels, payload, one, W, width):
    """Rows of eq:name parameters under an equation header, the others after a separator, as
    Stata's bayesstats tables (case 1909)."""
    seen, done = [], False
    for lab in labels:
        if ":" in lab and lab.split(":", 1)[0] not in seen:
            seen.append(lab.split(":", 1)[0])
    for k, e in enumerate(seen):
        if k:
            out.write("-" * (W + 1) + "+" + "-" * width + "\n", "text")
        out.write(f"{e:<{W}} |\n", "text")
        for lab, pl in zip(labels, payload):
            if ":" in lab and lab.split(":", 1)[0] == e:
                one(lab.split(":", 1)[1], pl)
        done = True
    rest = [(lab, pl) for lab, pl in zip(labels, payload) if ":" not in lab]
    if rest and done:
        out.write("-" * (W + 1) + "+" + "-" * width + "\n", "text")
    for lab, pl in rest:
        one(lab, pl)


def _flat_table(s, labels, summ, level, use_hpd) -> None:
    out = s.output
    W = 12
    out.write("-" * (W + 66) + "\n", "text")
    out.write(f"{'':>{W}} |" + " " * 48 + ("         HPD" if use_hpd else "Equal-tailed") + "\n", "text")
    out.write(f"{'':>{W}} |      Mean   Std. Dev.     MCSE     Median  [{level:g}% Cred. Interval]\n", "text")
    out.write("-" * (W + 1) + "+" + "-" * 64 + "\n", "text")

    def one(lab, r):
        out.write(f"{lab[-W:]:>{W}} |" + _srow(r) + "\n", "result")
    _grouped(out, labels, summ, one, W, 64)
    out.write("-" * (W + 66) + "\n", "text")


def _stored_models(s, names: str):
    if not names.strip():
        est = _current(s)
        return [("active", est)]
    out = []
    store = getattr(s, "_stored_estimates", {})
    for n in names.split():
        if n == ".":
            out.append(("active", _current(s)))
            continue
        snap = store.get(n)
        if snap is None:
            raise StataError(111, f"estimation result {n} not found")
        est = snap["est"]
        if est is None or est.cmd != "bayesmh":
            raise StataError(198, f"estimation result {n} is not from bayesmh")   # VERIFY
        out.append((n, est))
    return out


def _stats_ic(s, est, body: str, o) -> None:
    models = _stored_models(s, body)
    out = s.output
    out.write("\nBayesian information criteria\n\n", "text")
    W = 12
    out.write("-" * (W + 36) + "\n", "text")
    out.write(f"{'':>{W}} |   Avg. log DIC  log(ML)   log(BF)\n".replace("Avg. log DIC", "       DIC"), "text")
    out.write("-" * (W + 1) + "+" + "-" * 34 + "\n", "text")
    base = models[0][1].extra["lml"]
    for name, e in models:
        x = e.extra
        bf = x["lml"] - base
        out.write(f"{name:>{W}} |  {g9(x['dic']):>9}  {g9(x['lml']):>8}  {g9(bf):>8}\n", "result")
    out.write("-" * (W + 36) + "\n", "text")
    out.write("Note: Marginal likelihood (ML) is computed using Laplace-Metropolis approximation.\n", "text")


def cmd_bayestest(s: "Session", args: str) -> None:
    est = _current(s)
    m = re.match(r"\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "bayestest requires interval or model")   # VERIFY
    sub = m.group(1).lower()
    rest = m.group(2)
    subs = [w for w in ("interval", "model") if w.startswith(sub)]
    if len(subs) != 1:
        raise StataError(198, f"bayestest {sub}: invalid subcommand")   # VERIFY
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    out = s.output
    if subs[0] == "interval":
        o = match_options(opt, {"lower": 3, "upper": 3, "prior": 5}) if opt.strip() else {}
        items = _param_items(est, body)
        n = len(items[0][1])
        out.write("\nInterval tests     " + f"MCMC sample size = {_commas(n):>9}\n\n", "text")
        lo = float(o["lower"]) if o.get("lower") else -np.inf
        hi = float(o["upper"]) if o.get("upper") else np.inf
        res = []
        for k, (lab, d, _i) in enumerate(items, 1):
            cond = (f"{g9(lo)} < {{{lab}}} < {g9(hi)}" if np.isfinite(lo) and np.isfinite(hi) else
                    f"{{{lab}}} > {g9(lo)}" if np.isfinite(lo) else f"{{{lab}}} < {g9(hi)}")
            out.write(f"{'prob' + str(k):>12} : {cond}\n", "text")
            ind = ((d > lo) & (d < hi)).astype(float)
            p = float(ind.mean())
            sd = float(np.std(ind, ddof=1))
            e = mcmc.ess(ind) if sd > 0 else float(len(ind))
            res.append((f"prob{k}", p, sd, sd / math.sqrt(e) if e > 0 else 0.0))
        W = 12
        out.write(" \n" + "-" * (W + 35) + "\n", "text")
        out.write(f"{'':>{W}} |      Mean    Std. Dev.      MCSE\n", "text")
        out.write("-" * (W + 1) + "+" + "-" * 33 + "\n", "text")
        for lab, p, sd, mc in res:
            out.write(f"{lab:>{W}} |", "text")
            out.write(f"{g9(p):>10}{sd:>12.5f}{g9(mc):>11}\n", "result")
        out.write("-" * (W + 35) + "\n", "text")
        return
    models = _stored_models(s, body)
    if len(models) < 2:
        raise StataError(198, "bayestest model requires at least two models")   # VERIFY
    lml = np.array([e.extra["lml"] for _n, e in models])
    if not np.all(np.isfinite(lml)):
        raise StataError(498, "log marginal likelihood is not available for these models")   # VERIFY
    pri = np.full(len(models), 1.0 / len(models))
    post_ = pri * np.exp(lml - lml.max())
    post_ /= post_.sum()
    out.write("\nBayesian model tests\n\n", "text")
    W = 12
    out.write("-" * (W + 36) + "\n", "text")
    out.write(f"{'':>{W}} |    log(ML)    P(M)  P(M|y)\n", "text")
    out.write("-" * (W + 1) + "+" + "-" * 34 + "\n", "text")
    for (n, e), pr, po in zip(models, pri, post_):
        out.write(f"{n:>{W}} |  {g9(e.extra['lml']):>9}  {pr:>6.4f}  {po:>6.4f}\n", "result")
    out.write("-" * (W + 36) + "\n", "text")


def cmd_bayesgraph(s: "Session", args: str) -> None:
    est = _current(s)
    m = re.match(r"\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "bayesgraph requires a graph type")   # VERIFY
    sub = m.group(1).lower()
    kinds = ("diagnostics", "trace", "histogram", "kdensity", "ac", "cusum")
    hit = [k for k in kinds if k.startswith(sub)]
    if len(hit) != 1:
        raise StataError(198, f"bayesgraph {sub}: invalid graph type")   # VERIFY
    kind = hit[0]
    rest = m.group(2)
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    items = _param_items(est, body)
    if len(items) > 1 and kind != "diagnostics":
        pass
    for lab, d, _i in items:
        if kind == "diagnostics":
            names = []
            for k in ("trace", "histogram", "ac", "kdensity"):
                nm = f"__bg_{k}"
                _one_graph(s, k, lab, d, opt, name=nm)
                names.append(nm)
            s.run(f"graph combine {' '.join(names)}")
        else:
            _one_graph(s, kind, lab, d, opt)


def _one_graph(s, kind: str, lab: str, d: np.ndarray, opt: str, name: str | None = None) -> None:
    from ..core.dataset import Variable
    n = len(d)
    ttl = lab
    s.run("preserve")
    try:
        s.run("quietly clear")
        if kind == "trace":
            s.run(f"quietly set obs {n}")
            s.data.add(Variable("iter", "long", np.arange(1, n + 1, dtype=float)))
            s.data.add(Variable("theta", "double", np.asarray(d, float)))
            cmd = f'twoway line theta iter, title("{ttl}") xtitle("Iteration number") ytitle("") legend(off)'
        elif kind == "histogram":
            s.run(f"quietly set obs {n}")
            s.data.add(Variable("theta", "double", np.asarray(d, float)))
            cmd = f'histogram theta, title("{ttl}") xtitle("") legend(off)'
        elif kind == "kdensity":
            s.run(f"quietly set obs {n}")
            s.data.add(Variable("theta", "double", np.asarray(d, float)))
            cmd = f'kdensity theta, title("{ttl}") xtitle("") legend(off) note("")'
        elif kind == "ac":
            ml = min(40, n - 1)
            ac = mcmc.autocorr(np.asarray(d, float), ml)
            s.run(f"quietly set obs {ml}")
            s.data.add(Variable("lag", "long", np.arange(1, ml + 1, dtype=float)))
            s.data.add(Variable("ac", "double", ac[1:]))
            cmd = f'twoway (spike ac lag), title("{ttl}") xtitle("Lag") ytitle("Autocorrelation") legend(off)'
        else:
            s.run(f"quietly set obs {n}")
            cs = np.cumsum(d - np.mean(d))
            s.data.add(Variable("iter", "long", np.arange(1, n + 1, dtype=float)))
            s.data.add(Variable("cusum", "double", cs))
            cmd = f'twoway line cusum iter, title("{ttl}") xtitle("Iteration number") ytitle("Cusum") legend(off)'
        if name:
            cmd += f" name({name}, replace)"
        s.run(cmd)
    finally:
        s.run("restore")
