"""mlexp ([R] mlexp): maximum likelihood with the log likelihood of each observation given as an
expression, fitted by the ml engine of mlcmd.py."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from .mlcmd import _OPT, _balanced, _init, _maximize
from .mlmodel import Eq, Model, _full

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_BRACE = re.compile(r"\{([^{}]*)\}")


def cmd_mlexp(s: "Session", args: str) -> None:
    t = args.strip()
    if not t.startswith("("):
        raise StataError(198, "invalid syntax")
    j = _balanced(t, 0)
    expr, tail = t[1:j].strip(), t[j + 1:]
    p = parse_standard(tail)
    o = match_options(p.options, dict(_OPT, initial=4, variables=3)) if p.options.strip() else {}
    if o.get("constraints") or o.get("svy"):
        raise StataError(198, "option constraints() and svy are not available yet in OpenDTA")
    ds = s.data
    combos: dict[str, Eq] = {}
    scalars: dict[str, float] = {}
    order: list[str] = []

    def sub(mt):
        body = mt.group(1).strip()
        if ":" in body:
            name, _, rhs = body.partition(":")
            name = name.strip()
            lst, _, opt = rhs.partition(",")
            if name not in combos:
                words = lst.split()
                has_cons = "_cons" in words      # VERIFY: the manual lists _cons explicitly
                vars_ = s.expand_varlist(" ".join(w for w in words if w != "_cons")) if words else []
                combos[name] = Eq(name, [], vars_, has_cons or "noconstant" not in opt)
                order.append(name)
            return f"\x00c{name}\x00"
        name, _, init = body.partition("=")
        name = name.strip()
        if name not in scalars:
            scalars[name] = float(s.eval(init)) if init.strip() else 0.0
            order.append("/" + name)
        return f"\x00s{name}\x00"
    template = _BRACE.sub(sub, expr)
    if not order:
        raise StataError(198, "no parameters in the expression")
    from ..commands._util import touse
    from ..commands.summarize import weights
    mask = touse(s, p)
    w, wtype, mask = weights(s, p, mask, ("fweight", "aweight", "pweight", "iweight"))
    used = {v for e in combos.values() for v in e.indep}
    used |= {x for x in re.findall(r"[A-Za-z_]\w*", _BRACE.sub("", expr)) if ds.has(x)}
    clustvar = ""
    vce = o.get("vce")
    if isinstance(o.get("cluster"), str):
        clustvar = o["cluster"]
    elif isinstance(vce, str) and vce.split()[0].lower().startswith("cl"):
        clustvar = vce.split()[1]
    if clustvar:
        used.add(clustvar)
    for v in used:
        var = ds.get(v)
        mask = mask & (np.array([x != "" for x in var.raw]) if var.is_string else var.data < SYS)
    eqs: list[Eq] = []
    pos = 0
    for key in order:
        if key.startswith("/"):
            e = Eq(key[1:], [], [], True, scalar=True)
        else:
            e = combos[key]
        cols = [ds.get(v).data[mask] for v in e.indep]
        if e.cons:
            cols.append(np.ones(int(mask.sum())))
        e.X = np.column_stack(cols) if cols else np.zeros((int(mask.sum()), 0))
        e.off = np.zeros(int(mask.sum()))
        e.start = pos
        pos += e.k
        eqs.append(e)
    m = Model("mlexp", "", eqs, mask, None if w is None else w[mask], wtype, o, f"mlexp {args}")
    m.N = int(mask.sum())
    if m.N == 0:
        raise StataError(2000, "no observations")
    if clustvar:
        cv = ds.get(clustvar)
        m.cluster = np.asarray(cv.raw if cv.is_string else cv.data)[mask]
        m.clustvar = clustvar
    for e in eqs:
        if e.scalar:
            m.names.append(("/", e.name))
        else:
            m.names += [(e.name, v) for v in e.indep] + ([(e.name, "_cons")] if e.cons else [])
    m.template, m.lexp = template, expr
    m.b0 = np.zeros(m.nb)
    for e in eqs:
        if e.scalar:
            m.b0[e.start] = scalars[e.name]
    if isinstance(o.get("from"), str):
        s._ml = m
        _init(s, o["from"])
    elif isinstance(o.get("initial"), str):
        s._ml = m
        _init(s, " ".join(f"{k}={v}" if "=" not in k else k for k, v in
                          zip(o["initial"].split()[::2], o["initial"].split()[1::2])))
    elif not m.b0.any():
        m.b0 = None          # no starting values given: the search steps apply
    s._ml = m
    _maximize(s, "")


def _mlexp_obs(s: "Session", m: Model, th: np.ndarray) -> np.ndarray:
    """mlexp: the expression evaluated for the given linear predictors (one column per
    parameter group, like lf); returns the log likelihood per observation (NaN when missing)."""
    from ..commands._util import eval_vector
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    ds = s.data
    names = getattr(m, "_names", None)
    if names is None or any(not ds.has(n) for n in names.values()):
        names = {e.name: _temp_name(s) for e in m.eqs if not e.scalar}
        for nm in names.values():
            ds.add(Variable(nm, "double", np.full(ds.nobs, SYS)))
        m._names = names
        m._tmp = list(names.values())
    text = m.template
    for j, e in enumerate(m.eqs):
        if e.scalar:
            text = text.replace(f"\x00s{e.name}\x00", f"({float(th[0, j])!r})")
        else:
            ds.get(names[e.name]).data = _full(m, th[:, j], ds.nobs)
            text = text.replace(f"\x00c{e.name}\x00", names[e.name])
    v = np.broadcast_to(np.asarray(eval_vector(s, text), dtype=float), (ds.nobs,))[m.touse]
    return np.where(v >= SYS, np.nan, v)
