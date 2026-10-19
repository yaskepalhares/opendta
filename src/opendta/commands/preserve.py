"""preserve/restore and tempvar/tempname/tempfile ([P] preserve, [P] macro).

Temporaries and preserve belong to the current scope (program or
do-file) and are undone when it ends, with or without an error. Outside
programs and do-files, they belong to the interactive session.

preserve keeps a copy of the data in memory (the native vectors of
each variable), which is fast; Stata writes a temporary file.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.dataset import NAME_RE
from ..core.errors import StataError
from ..lang.programs import Scope
from ..lang.syntax import match_options
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def current_scope(s: "Session") -> Scope:
    if s.scopes:
        return s.scopes[-1]
    top = getattr(s, "top_scope", None)
    if top is None:
        top = s.top_scope = Scope("interactive")
    return top


def restore_preserved(s: "Session", scope: Scope) -> None:
    """End of the program or do-file: back to the preserved data."""
    if scope.preserved is not None:
        s.data = scope.preserved
        scope.preserved = None
        s.notify_state()


@command("preserve")
def cmd_preserve(s: "Session", args: str) -> None:
    head, comma, opts = args.partition(",")
    if head.strip():
        raise StataError(198, "invalid syntax")
    if comma:
        match_options(opts, {"changed": 6})
    scope = current_scope(s)
    if scope.preserved is not None:
        raise StataError(621, "already preserved")
    scope.preserved = s.data.copy()


@command("restore")
def cmd_restore(s: "Session", args: str) -> None:
    head, comma, opts = args.partition(",")
    if head.strip():
        raise StataError(198, "invalid syntax")
    o = match_options(opts, {"not": 3, "preserve": 8}) if comma else {}
    scope = next((sc for sc in reversed(s.scopes) if sc.preserved is not None), None)
    if scope is None:
        top = getattr(s, "top_scope", None)
        scope = top if top is not None and top.preserved is not None else None
    if scope is None:
        raise StataError(622, "nothing to restore")
    if o.get("not"):
        scope.preserved = None
        return
    snapshot = scope.preserved
    s.data = snapshot.copy() if o.get("preserve") else snapshot
    if not o.get("preserve"):
        scope.preserved = None
    s.notify_state()


def _counter(s: "Session") -> str:
    k = getattr(s, "_temp_counter", 0)
    while True:
        name = f"__{k:06X}"
        k += 1
        if not s.data.has(name) and name not in s.scalars and name not in getattr(s, "matrices", {}):
            s._temp_counter = k
            return name


def _names(args: str) -> list[str]:
    names = args.split()
    if not names:
        raise StataError(100, "something required")
    for n in names:
        if not NAME_RE.match(n):
            raise StataError(198, f"{n} invalid name")
    return names


@command("tempvar")
def cmd_tempvar(s: "Session", args: str) -> None:
    scope = current_scope(s)
    for n in _names(args):
        t = _counter(s)
        scope.tempvars.append(t)
        s.macros.set_local(n, t)


@command("tempname")
def cmd_tempname(s: "Session", args: str) -> None:
    scope = current_scope(s)
    for n in _names(args):
        t = _counter(s)
        scope.tempnames.append(t)
        s.macros.set_local(n, t)


@command("tempfile")
def cmd_tempfile(s: "Session", args: str) -> None:
    scope = current_scope(s)
    folder = Path(tempfile.gettempdir()) / f"opendta-{os.getpid()}"
    folder.mkdir(exist_ok=True)
    for n in _names(args):
        t = folder / f"S{_counter(s)}"
        scope.tempfiles.append(t)
        s.macros.set_local(n, str(t))
