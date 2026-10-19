"""unab, tsunab, fvunab, unabcmd, fvexpand, fvrevar, tsrevar and _rmcoll
([P] unab, [P] unabcmd, [R] fvrevar, [P] fvexpand, [TS] tsrevar, [P] _rmcoll).

They answer programs that receive a varlist: unabbreviate it, expand its factor-variable
and time-series operators, replace them by (temporary) variables and flag the collinear
ones. The terms and their columns come from estimation/fvars.py, the same ones the models use.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import NAME_RE, Variable
from ..core.errors import StataError
from ..estimation import fvars, tsops
from ..lang.adopath import find_ado
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .preserve import current_scope
from .registry import command, lookup
from .summarize import weights
from .syntaxcmd import _temp_name

if TYPE_CHECKING:
    from ..session import Session

# commands the interpreter runs itself (they are not in the registry)
_INTERPRETER_WORDS = ("if", "else", "foreach", "forvalues", "while", "by", "bysort", "input",
                      "mata")


# ---------------------------------------------------------------------------
# unab, tsunab, fvunab
# ---------------------------------------------------------------------------

def _count_error(n: int, lo: int, hi: int, name: str) -> StataError:
    if lo == hi:
        need = f"{lo} variable{'' if lo == 1 else 's'}"
    elif hi - lo == 1:
        need = f"{lo} or {hi} variables"      # VERIFY: plural of the 2-value range
    else:
        need = f"{lo} - {hi} variables"
    if n > hi:
        return StataError(103, f"too many variables specified\n{name}: {need} required")   # VERIFY
    return StataError(102, f"too few variables specified\n{name}: {need} required")        # VERIFY


def _unab(s: "Session", args: str, kind: str) -> None:
    head, colon, rest = args.partition(":")
    lname = head.strip()
    if not colon or not NAME_RE.match(lname):
        raise StataError(198, "invalid syntax")
    text, comma, opts = rest.partition(",")
    o = match_options(opts, {"min": 3, "max": 3, "name": 4}) if comma else {}
    lo = int(float(s.eval(str(o["min"])))) if o.get("min") not in (None, True) else 1
    hi = int(float(s.eval(str(o["max"])))) if o.get("max") not in (None, True) else 32000
    label = str(o["name"]).strip().strip('"') if o.get("name") not in (None, True) else text.strip()
    ds = s.data
    abbrev = s.settings.get("varabbrev", "on") == "on"
    text = text.strip()
    if kind == "unab":
        if fvars.has_fv_or_ts(text):
            raise StataError(101, "factor variables and time-series operators not allowed")   # VERIFY
        names = s.expand_varlist(text) if text else []
    else:
        names = fvars.unabbreviate(ds, text, abbrev, fv=kind == "fvunab")
    if not lo <= len(names) <= hi:
        raise _count_error(len(names), lo, hi, label)
    s.macros.set_local(lname, " ".join(names))


@command("unab")
def cmd_unab(s: "Session", args: str) -> None:
    _unab(s, args, "unab")


@command("tsunab")
def cmd_tsunab(s: "Session", args: str) -> None:
    _unab(s, args, "tsunab")


@command("fvunab")
def cmd_fvunab(s: "Session", args: str) -> None:
    _unab(s, args, "fvunab")


@command("unabcmd")
def cmd_unabcmd(s: "Session", args: str) -> None:
    word = args.strip()
    if not word or len(word.split()) != 1:
        raise StataError(198, "invalid syntax")
    spec = lookup(word)
    if spec is not None:
        full = spec.name
    elif word in _INTERPRETER_WORDS or word in s.programs or find_ado(s, word) is not None:
        full = word
    else:
        raise StataError(111, f"command {word} not found as either built-in or ado-file")
    s.r = {"cmd": full}


# ---------------------------------------------------------------------------
# terms for fvexpand, fvrevar and _rmcoll
# ---------------------------------------------------------------------------

def _terms(s: "Session", text: str) -> list[fvars.Term]:
    """The terms of the varlist, with the factors written before the continuous variables inside
    an interaction (c.x#i.g is named g#c.x, as in the examples of [R] fvrevar)."""
    if not text.strip():
        raise StataError(100, "varlist required")
    abbrev = s.settings.get("varabbrev", "on") == "on"
    terms = fvars.expand_fv(s.data, text, abbrev)
    for t in terms:
        if t.is_interaction:
            t.comps.sort(key=lambda c: not c.factor)
    if any(c.ts for t in terms for c in t.comps):
        tsops.require_ts(s.data)
    return terms


@command("fvexpand")
def cmd_fvexpand(s: "Session", args: str) -> None:
    p = parse_standard(args)
    terms = _terms(s, p.varlist)
    cols = fvars.build_design(s.data, terms, touse(s, p))
    s.r = {"varlist": " ".join(c.name for c in cols)}


def _indicator(ds, term: fvars.Term, cell: tuple) -> np.ndarray:
    """Values of the column of `term` at the factor levels `cell`, in all the observations."""
    y = np.ones(ds.nobs)
    levels = iter(cell)
    for c in term.comps:
        x = fvars.comp_values(ds, c)
        if c.factor:
            y = y * (x == next(levels))
        else:
            y = np.where((y == 0) | (x >= M.SYSMISS), np.where(y == 0, 0.0, M.SYSMISS), y * x)
    return y


def _column_type(ds, term: fvars.Term) -> str:
    if all(c.factor and not c.ts for c in term.comps):
        return "byte"
    return "double" if any(ds.get(c.var).vtype == "double" for c in term.comps) else "float"


def _new_variable(s: "Session", term: fvars.Term, values: np.ndarray, stub: str,
                  made: list[str]) -> str:
    """A temporary variable (or stub#) with the values; the one made earlier by another call with
    the same values is used again."""
    ds = s.data
    cache = getattr(s, "_revar_cache", None)
    if cache is None:
        cache = s._revar_cache = {}
    if not stub:
        for name in cache.get(id(ds), []):
            if ds.has(name) and np.array_equal(ds.get(name).data, values):
                return name
    if stub:
        name = f"{stub}{len(made) + 1}"
        if ds.has(name):
            raise StataError(110, f"variable {name} already defined")
    else:
        name = _temp_name(s)
        current_scope(s).tempvars.append(name)
        cache.setdefault(id(ds), []).append(name)
    ds.add(Variable(name, _column_type(ds, term), values))
    made.append(name)
    return name


def _revar(s: "Session", args: str, command_name: str) -> None:
    p = parse_standard(args)
    spec = {"substitute": 3, "list": 1}
    if command_name == "fvrevar":
        spec.update({"tsonly": 2, "stub": 4})
    o = match_options(p.options, spec) if p.options.strip() else {}
    if o.get("list") and (o.get("substitute") or o.get("tsonly") or o.get("stub")):
        raise StataError(198, "options list and substitute, tsonly or stub() may not be combined")
    ds = s.data
    terms = _terms(s, p.varlist)
    if command_name == "tsrevar" and any(t.is_interaction or t.has_factor or t.comps[0].explicit_c
                                         for t in terms):
        raise StataError(101, "factor variables not allowed")      # VERIFY
    if o.get("list"):
        names = sorted(fvars.base_vars(terms), key=ds.index)
        s.r = {"varlist": " ".join(names)}
        return
    stub = str(o.get("stub", "")) if o.get("stub") not in (None, True) else ""
    if stub and not NAME_RE.match(stub):
        raise StataError(198, "option stub() incorrectly specified")
    mask = touse(s, p)
    made: list[str] = []
    out: list[str] = []
    for ti, t in enumerate(terms):
        plain = not t.is_interaction and not t.has_factor and not t.comps[0].explicit_c
        if plain and not t.comps[0].ts:
            out.append(t.comps[0].var)
            continue
        if o.get("tsonly") and not plain:
            cols = fvars.build_design(ds, [t], mask)
            out.extend(c.name for c in cols)
            continue
        for c in fvars.build_design(ds, [t], mask):
            out.append(_new_variable(s, t, _indicator(ds, t, c.cell), stub, made))
    s.r = {"varlist": " ".join(out)}


@command("fvrevar")
def cmd_fvrevar(s: "Session", args: str) -> None:
    _revar(s, args, "fvrevar")


@command("tsrevar")
def cmd_tsrevar(s: "Session", args: str) -> None:
    _revar(s, args, "tsrevar")


# ---------------------------------------------------------------------------
# _rmcoll
# ---------------------------------------------------------------------------

@command("_rmcoll")
def cmd_rmcoll(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"noconstant": 3, "collinear": 4, "expand": 3,
                                  "forcedrop": 5}) if p.options.strip() else {}
    ds = s.data
    terms = _terms(s, p.varlist)
    if o.get("forcedrop") and any(t.is_interaction or t.has_factor for t in terms):
        raise StataError(198, "option forcedrop not allowed with factor variables")   # VERIFY
    mask = touse(s, p) & ~fvars.term_missing(ds, terms)
    w, _, mask = weights(s, p, mask, ("fweight", "aweight", "iweight", "pweight"))
    constant = not o.get("noconstant")
    cols = fvars.build_design(ds, terms, mask, constant=constant)
    notes: dict[int, str] = {}
    for k, c in enumerate(cols):
        if c.values is not None and not np.any(c.values):
            c.omitted, c.values = True, None
            notes[k] = f"note: {c.name} identifies no observations in the sample"
    if not o.get("collinear"):
        flagged = fvars.drop_collinear(cols, constant, None if w is None else w[mask])
        for k, c in enumerate(cols):
            if c.omitted and c.name in flagged:
                notes[k] = f"note: {c.name} omitted because of collinearity"
    for k in sorted(notes):
        s.output.write(notes[k] + "\n", "text")
    if o.get("forcedrop"):
        names = [c.name for c in cols if not c.omitted]
    elif o.get("expand"):
        names = [fvars.omitted_name(c.name) if c.omitted and not c.base else c.name for c in cols]
    else:
        names = []
        for ti, t in enumerate(terms):
            mine = [c for c in cols if c.term == ti]
            if t.is_interaction or t.has_factor:
                names.append(fvars.specific_name(ds, t, mask, mine))
            else:
                names.append(fvars.omitted_name(mine[0].name) if mine[0].omitted else mine[0].name)
    s.r = {"k_omitted": float(sum(1 for c in cols if c.omitted and not c.base)),
           "varlist": " ".join(names)}
