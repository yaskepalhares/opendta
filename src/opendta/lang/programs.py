"""User-defined programs ([P] program) and execution scopes.

    program [define] name [, rclass eclass sclass nclass byable(recall|onecall)
                            properties(...) sortpreserve plugin]
        ...
    end

When `name arguments` is called, the body runs in a new scope of local macros
with `0' = everything after the name and `1', `2', ... = the words.

Each program (and each do-file) opens a Scope, which keeps what must be
undone on exit: temporary variables, names and files, and the data
kept by preserve. rclass/sclass programs accumulate the results of
`return`/`sreturn` and only publish them in r()/s() when they end without an error.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import ExitRequest, StataError
from .lexer import LogicalLine
from .words import split_words

if TYPE_CHECKING:
    from ..session import Session

_CLASSES = ("rclass", "eclass", "sclass", "nclass")


@dataclass
class Program:
    name: str
    lines: list[LogicalLine]
    kind: str = "nclass"                 # rclass | eclass | sclass | nclass
    byable: str = ""                     # recall | onecall | ""
    properties: str = ""
    source: str = ""                     # source .ado file, if any
    options: str = ""
    sortpreserve: bool = False           # restores the sort order of the data at the end
    plugin: str = ""                     # file named by `program name, plugin using()`


@dataclass
class Scope:
    kind: str                            # program | do
    name: str = ""
    program: Program | None = None
    pending: dict = field(default_factory=dict)      # return/sreturn not yet published
    tempvars: list[str] = field(default_factory=list)
    tempnames: list[str] = field(default_factory=list)
    tempfiles: list[Path] = field(default_factory=list)
    preserved: object = None             # ver commands/preserve.py


def parse_definition(rest: str) -> tuple[str, dict[str, str]]:
    """Text after `program [define]` -> (name, options)."""
    head, comma, opts = rest.partition(",")
    words = head.split()
    if words and len(words[0]) >= 2 and "define".startswith(words[0]):
        words = words[1:]
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    name = words[0]
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,31}", name):
        raise StataError(198, f"{name} invalid name")
    options: dict[str, str] = {}
    for m in re.finditer(r"([A-Za-z_]+)(?:\(([^)]*)\))?", opts if comma else ""):
        options[m.group(1).lower()] = m.group(2) or ""
    return name, options


def make_program(name: str, options: dict[str, str], lines: list[LogicalLine],
                 source: str = "") -> Program:
    kind = next((c for c in _CLASSES if c in options), "nclass")
    byable = options.get("byable", "")
    if "byable" in options and byable not in ("recall", "onecall"):
        raise StataError(198, "option byable() must be recall or onecall")
    return Program(name, lines, kind, byable, options.get("properties", ""), source,
                   sortpreserve="sortpreserve" in options)


def define_plugin(s: "Session", name: str, options: dict[str, str]) -> None:
    """program name, plugin [using(file)]: the file is looked for as Stata does, but OpenDTA cannot
    load a compiled Stata plugin."""
    from .adopath import find_file
    fname = options.get("using", "").strip().strip('"') or f"{name}.plugin"
    if not Path(fname).expanduser().is_file() and not find_file(s, fname):
        raise StataError(601, f"file {fname} not found")
    raise StataError(198, f"plugin {fname} cannot be loaded: OpenDTA does not run compiled "
                          "Stata plugins")


def call_program(s: "Session", prog: Program, args: str) -> None:
    """Runs the program with the given arguments."""
    argv = split_words(args)
    frame = {str(k): v for k, v in enumerate(argv, start=1)}
    if args.strip():
        frame["0"] = args.strip()
    if len(s.scopes) > 60:
        raise StataError(1000, "system limit exceeded: nested program calls")   # VERIFY
    s.macros.push_frame(frame)
    scope = Scope("program", prog.name, prog)
    s.scopes.append(scope)
    index = None
    if prog.sortpreserve:
        index = _mark_order(s)
    ok = False
    try:
        try:
            s.interp.run_program(prog)
        except ExitRequest as e:
            if e.clear_all:
                raise
            if e.rc:
                raise StataError(e.rc, "")
        ok = True
    finally:
        s.scopes.pop()
        s.macros.pop_frame()
        s.close_scope(scope)       # restores the preserved data first
        if index is not None:
            _restore_order(s, index, ok)
    if ok:
        if prog.kind == "rclass":
            s.r = scope.pending
        elif prog.kind == "sclass":
            s.sret = scope.pending


def _mark_order(s: "Session") -> tuple[str, list[str]]:
    """sortpreserve: a variable with the order of the observations, whose name the program finds in
    the local macro _sortindex, and the sort variables."""
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    import numpy as np
    ds = s.data
    name = _temp_name(s)
    ds.add(Variable(name, "long", np.arange(1, ds.nobs + 1, dtype=np.float64)))
    s.macros.set_local("_sortindex", name)
    return name, list(ds.sortlist)


def _restore_order(s: "Session", index: tuple[str, list[str]], ok: bool) -> None:
    name, sortlist = index
    ds = s.data
    if not ds.has(name):
        if ok:
            # VERIFY: return code of this error
            raise StataError(111, "could not restore sort order because variables were dropped")
        return
    import numpy as np
    v = ds.get(name)
    order = np.argsort(v.data, kind="stable")
    if not np.array_equal(order, np.arange(ds.nobs)):
        ds.reorder_obs(order)
    changed = ds.changed
    ds.drop_vars([name])
    ds.changed = changed
    ds.sortlist = [x for x in sortlist if ds.has(x)]
    s.notify_state()


def current_program_scope(s: "Session") -> Scope | None:
    for sc in reversed(s.scopes):
        if sc.kind == "program":
            return sc
    return None
