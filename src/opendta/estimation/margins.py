"""margins ([R] margins) e nlcom ([R] nlcom): previsões médias, efeitos
marginais e combinações não lineares, com erros-padrão pelo método delta.

As derivadas em relação a b e a x são numéricas (diferenças centrais)."""

from __future__ import annotations

import itertools
import re
from typing import TYPE_CHECKING

import numpy as np
from scipy import special, stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from . import fvars
from .results import LINE, CoefRow, coef_table, current, crit, g9, pvalue

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# previsão em função de b e de valores fixados das variáveis
# ---------------------------------------------------------------------------

def _xb(s: "Session", est, b: np.ndarray, rows: np.ndarray, overrides: dict) -> np.ndarray:
    """xb nas linhas `rows`, com variáveis trocadas por `overrides` (var -> vetor)."""
    ds = s.data
    n = len(rows)
    out = np.zeros(n)
    for j, c in enumerate(est.cols):
        if c.values is None and not c.base:
            continue
        if b[j] == 0 and c.base:
            continue
        t = est.terms[c.term]
        y = np.ones(n)
        k = 0
        for comp in t.comps:
            if comp.var in overrides and not comp.ts:
                x = overrides[comp.var]
            else:
                x = fvars.comp_values(ds, comp)[rows]
            if comp.factor:
                y = y * (x == c.cell[k])
                k += 1
            else:
                y = y * x
        out += b[j] * y
    if est.constant:
        out += b[len(est.cols)]
    return out


def _transform(est, predict: str):
    """Função de xb que dá a previsão padrão do modelo, e o texto do Expression."""
    cmd = getattr(est, "ml_cmd", est.cmd)
    dep = est.depvar
    p = predict.strip()
    if p in ("xb", "xb()"):
        return (lambda xb: xb), "Linear prediction, predict(xb)"
    if est.cmd == "regress" or cmd in ("regress", "areg", "xtreg", "ivregress"):
        return (lambda xb: xb), "Linear prediction, predict()"
    if cmd in ("logit", "logistic"):
        return special.expit, f"Pr({dep}), predict()"
    if cmd == "probit":
        return st.norm.cdf, f"Pr({dep}), predict()"
    if cmd in ("poisson", "nbreg"):
        return (lambda xb: np.exp(xb)), "Predicted number of events, predict()"
    if cmd == "glm":
        return est.glm_ginv, f"Predicted mean {dep}, predict()"
    raise StataError(322, f"margins not yet supported by OpenDTA after {cmd}")


def _jacobian(f, b: np.ndarray) -> np.ndarray:
    k = b.size
    f0 = np.atleast_1d(f(b))
    J = np.zeros((f0.size, k))
    for i in range(k):
        h = 1e-6 * (abs(b[i]) + 1e-3)
        a, c = b.copy(), b.copy()
        a[i] += h
        c[i] -= h
        J[:, i] = (np.atleast_1d(f(a)) - np.atleast_1d(f(c))) / (2 * h)
    return J


def _parse_at(s: "Session", text: str) -> list[dict]:
    """at(x=(1 2) w=3) -> lista de dicionários {var: valor | "mean"}."""
    specs: list[tuple[str, list]] = []
    pairs = []
    i = 0
    while i < len(text):
        m = re.compile(r"\s*(\(\w+\))?\s*([\w.]+)\s*=\s*").match(text, i)
        if not m:
            break
        j = m.end()
        if j < len(text) and text[j] == "(":
            depth, k = 0, j
            while k < len(text):
                if text[k] == "(":
                    depth += 1
                elif text[k] == ")":
                    depth -= 1
                    if depth == 0:
                        break
                k += 1
            pairs.append((m.group(1), m.group(2), text[j:k + 1]))
            i = k + 1
        else:
            m2 = re.compile(r"\S+").match(text, j)
            pairs.append((m.group(1), m.group(2), m2.group(0) if m2 else ""))
            i = m2.end() if m2 else len(text)
    for stat, var, val in pairs:
        from ..core.varlist import resolve_name
        var = resolve_name(s.data, var)
        if val.startswith("("):
            vals = [float(v) for v in _numlist(val[1:-1])]
        else:
            vals = [float(val)]
        specs.append((var, vals))
    for m in re.finditer(r"\((mean|median)\)\s+([\w ]+?)(?=\s*\(|$)", text):
        for var in m.group(2).split():
            specs.append((var, ["mean"]))
    if not specs:
        return [{}]
    names = [v for v, _ in specs]
    return [dict(zip(names, combo)) for combo in itertools.product(*[vals for _, vals in specs])]


def _numlist(t: str) -> list[float]:
    out = []
    for part in t.replace(",", " ").split():
        m = re.fullmatch(r"(-?[\d.]+)\((-?[\d.]+)\)(-?[\d.]+)", part)
        if m:
            a, st_, b = float(m.group(1)), float(m.group(2)), float(m.group(3))
            x = a
            while (st_ > 0 and x <= b + 1e-9) or (st_ < 0 and x >= b - 1e-9):
                out.append(round(x, 12))
                x += st_
            continue
        if "/" in part:
            a, b = part.split("/")
            out.extend(float(v) for v in range(int(float(a)), int(float(b)) + 1))
            continue
        out.append(float(part))
    return out


def cmd_margins(s: "Session", args: str) -> None:
    est = current(s)
    p = parse_standard(args)
    o = match_options(p.options, {"dydx": 4, "eyex": 4, "dyex": 4, "eydx": 4, "at": 2, "atmeans": 3,
                                  "predict": 4, "level": 1, "post": 4, "noatlegend": 6, "nose": 4,
                                  "vce": 3, "over": 4, "asobserved": 3}) if p.options.strip() else {}
    ds = s.data
    from .postest import esample
    from ..commands._util import touse
    sample = esample(s) & touse(s, p)
    rows = np.nonzero(sample)[0]
    if not len(rows):
        raise StataError(2000, "no observations")
    f, expr = _transform(est, str(o.get("predict") or ""))
    b = est.b.copy()
    V = est.V
    lev = float(o["level"]) if o.get("level") else getattr(est, "level", 95.0)
    terms = fvars.expand_fv(ds, p.varlist) if p.varlist.strip() else []
    ats = _parse_at(s, str(o.get("at") or "")) if o.get("at") else [{}]
    atmeans = bool(o.get("atmeans"))
    means = {}
    if atmeans:
        for c in est.terms:
            for comp in c.comps:
                if not comp.factor:
                    means[comp.var] = float(np.mean(fvars.comp_values(ds, comp)[rows]))
    dydx_vars = []
    dydx_kind = ""
    for key in ("dydx", "eyex", "dyex", "eydx"):
        if o.get(key):
            dydx_kind = key
            txt = str(o[key]).strip()
            if txt in ("*", "_all"):
                dydx_vars = sorted({c.var for t in est.terms for c in t.comps}, key=lambda v: ds.index(v))
            else:
                dydx_vars = txt.split()
    results = []   # (rótulo do grupo, rótulo da linha, valor, J)

    def base_overrides(at: dict) -> dict:
        ov = {}
        for var, val in means.items():
            ov[var] = np.full(len(rows), val)
        for var, val in at.items():
            if val == "mean":
                ov[var] = np.full(len(rows), float(np.mean(ds.get(var).data[rows])))
            else:
                ov[var] = np.full(len(rows), float(val))
        return ov

    def margin_fn(ov):
        return lambda bb: float(np.mean(f(_xb(s, est, bb, rows, ov))))

    factor_vars = {c.var for t in est.terms for c in t.comps if c.factor}
    for ai, at in enumerate(ats):
        ov0 = base_overrides(at)
        at_lab = f"{ai + 1}" if len(ats) > 1 else ""
        if dydx_vars:
            for var in dydx_vars:
                from ..core.varlist import resolve_name
                var = resolve_name(ds, var.split(".")[-1])
                if var in factor_vars:
                    lv = sorted(set(ds.get(var).data[rows].tolist()))
                    basev = lv[0]
                    for c in est.cols:
                        if est.terms[c.term].comps[0].var == var and len(est.terms[c.term].comps) == 1 and c.base:
                            basev = c.cell[0]
                    for level in lv:
                        if level == basev:
                            continue
                        ov1 = dict(ov0, **{var: np.full(len(rows), level)})
                        ovb = dict(ov0, **{var: np.full(len(rows), basev)})

                        def fn(bb, ov1=ov1, ovb=ovb):
                            return margin_fn(ov1)(bb) - margin_fn(ovb)(bb)
                        results.append((var, f"{int(level)}", fn, at_lab))
                    continue
                x = ov0.get(var, ds.get(var).data[rows].astype(float))
                h = 1e-5 * (np.abs(x).mean() + 1e-3)

                def fn(bb, var=var, x=x, h=h):
                    up = dict(ov0, **{var: x + h})
                    dn = dict(ov0, **{var: x - h})
                    d = (f(_xb(s, est, bb, rows, up)) - f(_xb(s, est, bb, rows, dn))) / (2 * h)
                    if dydx_kind == "eyex":
                        d = d * x / f(_xb(s, est, bb, rows, ov0))
                    elif dydx_kind == "dyex":
                        d = d * x
                    elif dydx_kind == "eydx":
                        d = d / f(_xb(s, est, bb, rows, ov0))
                    return float(np.mean(d))
                results.append(("", var, fn, at_lab))
        elif terms:
            for t in terms:
                comps = t.comps
                for c in comps:
                    if c.var in factor_vars:
                        c.factor = True   # margins g: g é fator no modelo
                if not all(c.factor for c in comps):
                    raise StataError(198, "only factor variables and their interactions are allowed")   # VERIFICAR
                lvls = [sorted(set(ds.get(c.var).data[rows].tolist())) for c in comps]
                for cell in itertools.product(*lvls):
                    ov = dict(ov0, **{c.var: np.full(len(rows), lv) for c, lv in zip(comps, cell)})
                    lab = "#".join(f"{int(v)}" for v in cell)
                    head = "#".join(c.var for c in comps)
                    results.append((head, lab, margin_fn(ov), at_lab))
        else:
            results.append(("", "_cons", margin_fn(ov0), at_lab))
    vals = np.array([fn(b) for _, _, fn, _ in results])
    J = np.vstack([_jacobian(fn, b) for _, _, fn, _ in results])
    Vm = J @ V @ J.T
    se = np.sqrt(np.maximum(np.diag(Vm), 0))
    # exibição
    out = s.output
    if dydx_vars:
        title = "Conditional marginal effects" if atmeans else "Average marginal effects"
    elif atmeans or (o.get("at") and _all_fixed(est, ats)):
        title = "Adjusted predictions"
    else:
        title = "Predictive margins"
    vce_label = {"regress": "OLS"}.get(est.cmd, "Robust" if est.vcetype else "OIM")
    if est.vcetype:
        vce_label = "Robust"
    out.write("\n", "text")
    out.write(f"{title:<48}{'Number of obs':<18}= ", "text")
    out.write(f"{len(rows):>10,}\n", "result")
    out.write(f"Model VCE    : {vce_label}\n", "text")
    out.write("\n", "text")
    out.write(f"Expression   : {expr}\n", "text")
    if dydx_vars:
        lab = {"dydx": "dy/dx", "eyex": "ey/ex", "dyex": "dy/ex", "eydx": "ey/dx"}[dydx_kind]
        out.write(f"{lab} w.r.t. : {' '.join(dydx_vars)}\n", "text")
    if atmeans or o.get("at"):
        _at_legend(s, est, ats, means, rows, atmeans)
    out.write("\n", "text")
    stat = "t" if est.df_r else "z"
    rws = []
    last_head = None
    multi_at = len(ats) > 1
    for i, (head, lab, _, at_lab) in enumerate(results):
        if multi_at and not head and not dydx_vars:
            if last_head != "_at":
                rws.append(CoefRow("header", "_at"))
                last_head = "_at"
            rws.append(CoefRow("coef", f"{at_lab} ", vals[i], se[i]))
            continue
        if head:
            if head != last_head:
                if rws and rws[-1].kind != "blank":
                    rws.append(CoefRow("blank"))
                rws.append(CoefRow("header", head))
                last_head = head
            rws.append(CoefRow("coef", f"{lab} ", vals[i], se[i]))
        else:
            rws.append(CoefRow("coef", lab, vals[i], se[i]))
    col = "Margin" if not dydx_vars else {"dydx": "dy/dx", "eyex": "ey/ex", "dyex": "dy/ex",
                                          "eydx": "ey/dx"}[dydx_kind]
    _margins_table(s, rws, stat, est.df_r, lev, col)
    if dydx_vars and any(h for h, *_ in results):
        out.write("Note: dy/dx for factor levels is the discrete change from the base level.\n", "text")
    s.r = {}
    from ..commands.matrix import Matrix
    names = [f"{h}:{lab}" if h else lab for h, lab, *_ in results]
    s.r["b"] = Matrix(vals.reshape(1, -1), ["r1"], [n.replace(":", ".") for n in names])
    s.r["V"] = Matrix(Vm, list(s.r["b"].colnames), list(s.r["b"].colnames))
    s.r["N"] = float(len(rows))
    if o.get("post"):
        from .results import Estimates, post
        e2 = Estimates("margins", est.depvar, list(s.r["b"].colnames), vals, Vm, len(rows), stat, est.df_r)
        e2.rows = rws
        e2.level = lev
        e2.display = lambda ss, ee, **kw: _margins_table(ss, ee.rows, stat, est.df_r, lev, col)
        post(s, e2, [("N", float(len(rows)))], [("cmd", "margins"), ("properties", "b V")], sample)


def _all_fixed(est, ats) -> bool:
    fixed = set().union(*[set(a) for a in ats]) if ats else set()
    allv = {c.var for t in est.terms for c in t.comps}
    return allv <= fixed


def _at_legend(s, est, ats, means, rows, atmeans):
    out = s.output
    ds = s.data
    for i, at in enumerate(ats):
        lead = f"{i + 1}._at" if len(ats) > 1 else "at"
        items = []
        if atmeans:
            for var, val in means.items():
                items.append((var, "(mean)", val))
        for var, val in at.items():
            v = float(np.mean(ds.get(var).data[rows])) if val == "mean" else float(val)
            items.append((var, "(mean)" if val == "mean" else "", v))
        for j, (var, tag, val) in enumerate(items):
            first = f"{lead:<13}: " if j == 0 else " " * 15
            out.write(f"{first}{var:<16}{tag:<8}= {format_value(val, '%10.0g', pad=False).strip():>10}\n", "text")
        if len(ats) > 1 and i < len(ats) - 1:
            out.write("\n", "text")


def _margins_table(s, rws, stat, df, lev, col):
    out = s.output
    lvt = f"{lev:g}"
    out.write(LINE + "\n", "text")
    out.write(f"{'':>12} |            Delta-method\n", "text")
    stat_h = "t    P>|t|" if stat == "t" else "z    P>|z|"
    out.write(f"{'':>12} | {col:>10}   Std. Err.      {stat_h}     {'[' + lvt + '% Conf. Interval]':>20}\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    coef_rows_only(s, rws, stat, df, lev)
    out.write(LINE + "\n", "text")


def coef_rows_only(s, rws, stat, df, lev):
    """Corpo da tabela (sem cabeçalho e rodapé)."""
    out = s.output
    c = crit(stat, lev, df)
    for r in rws:
        if r.kind == "header":
            out.write(f"{r.label:>12} |\n", "text")
            continue
        if r.kind == "blank":
            out.write(f"{'':>12} |\n", "text")
            continue
        b, se = r.b, r.se
        out.write(f"{r.label[:12]:>12} |", "text")
        if se > 0:
            tv = b / se
            p = pvalue(stat, tv, df)
            out.write(f"  {g9(b):>9}  {g9(se):>9}  {tv:>7.2f}  {p:>6.3f}    {g9(b - c * se):>9}   "
                      f"{g9(b + c * se):>9}\n", "result")
        else:
            out.write(f"  {g9(b):>9}  {'(omitted)':>9}\n", "result")


# ---------------------------------------------------------------------------
# nlcom
# ---------------------------------------------------------------------------

def cmd_nlcom(s: "Session", args: str) -> None:
    est = current(s)
    head, comma, opts = args.partition(",")
    o = match_options(opts, {"level": 1, "post": 4, "iterate": 4}) if comma and opts.strip() else {}
    text = head.strip()
    if text.startswith("("):
        from .postest import _split_parens
        exprs = _split_parens(text)
    else:
        exprs = [text]
    items = []
    for i, e in enumerate(exprs):
        m = re.match(r"^\s*(\w+)\s*:\s*(.+)$", e)
        if m and not e.strip().startswith("_b["):
            items.append((m.group(1), m.group(2).strip()))
        else:
            items.append((f"_nl_{i + 1}", e.strip()))
    b0 = est.b.copy()

    def evaluator(expr):
        from .postest import _coef_index

        def fn(bb):
            def rep(m):
                kind, name = m.group(1), m.group(2)
                idx = _coef_index(est, s, name)
                if kind == "_se":
                    return repr(float(np.sqrt(max(est.V[idx, idx], 0))))
                return f"({float(bb[idx])!r})"
            t = re.sub(r"(_b|_coef|_se)\[([^\]]+)\]", rep, expr)
            v = s.eval(t)
            return float(v)
        return fn

    fns = [evaluator(e) for _, e in items]
    vals = np.array([f(b0) for f in fns])
    J = np.vstack([_jacobian(f, b0) for f in fns])
    Vm = J @ est.V @ J.T
    se = np.sqrt(np.maximum(np.diag(Vm), 0))
    out = s.output
    out.write("\n", "text")
    for name, e in items:
        out.write(f"{name:>12}:  {e}\n", "text")
    out.write("\n", "text")
    lev = float(o["level"]) if o.get("level") else getattr(est, "level", 95.0)
    rws = [CoefRow("coef", name, vals[i], se[i]) for i, (name, _) in enumerate(items)]
    stat = "t" if est.df_r else "z"
    coef_table(s, est.depvar, rws, stat=stat, df=est.df_r, level=lev)
    from ..commands.matrix import Matrix
    names = [n for n, _ in items]
    s.r = {"b": Matrix(vals.reshape(1, -1), ["r1"], names), "V": Matrix(Vm, names, names)}
    s.r = dict(reversed(list(s.r.items())))
    if o.get("post"):
        from .results import Estimates, post
        e2 = Estimates("nlcom", est.depvar, names, vals, Vm, est.N, stat, est.df_r)
        e2.rows = rws
        e2.level = lev
        e2.display = lambda ss, ee, **kw: coef_table(ss, ee.depvar, ee.rows, stat=stat, df=est.df_r, level=lev)
        post(s, e2, [("N", float(est.N))], [("cmd", "nlcom"), ("properties", "b V")], s.data.esample)
