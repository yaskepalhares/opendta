"""scalar — manual [P] scalar."""

from __future__ import annotations

import fnmatch
import re
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.formats import format_value
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("scalar", "sca")
def cmd_scalar(s: "Session", args: str) -> None:
    t = args.strip()
    sub, _, rest = t.partition(" ")
    if sub in ("list", "l", "li", "lis", "dir"):
        _list(s, rest.split())
        return
    if sub == "drop":
        names = rest.split()
        if names == ["_all"]:
            s.scalars.clear()
            return
        for pat in names:
            matches = [k for k in s.scalars if fnmatch.fnmatchcase(k, pat)]
            if not matches:
                raise StataError(111, f"scalar {pat} not found")
            for k in matches:
                del s.scalars[k]
        return
    if sub in ("define", "de", "def", "defi", "defin"):
        t = rest
    m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$", t, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    name = m.group(1)
    if len(name) > 32:
        raise StataError(198, f"{name} invalid name")
    s.scalars[name] = s.eval(m.group(2))


def _list(s: "Session", names: list[str]) -> None:
    out = s.output
    items = sorted(s.scalars.items())
    if names and names != ["_all"]:
        items = [(k, v) for k, v in items if any(fnmatch.fnmatchcase(k, p) for p in names)]
    for name, value in items:
        if isinstance(value, str):
            shown = value
        else:
            shown = format_value(value, "%10.0g", pad=True).strip().rjust(10)
        out.write(f"{name:>10} = ", "text")
        out.write(f"{shown}\n", "result")
