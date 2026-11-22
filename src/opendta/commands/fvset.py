"""fvset ([R] fvset): default base levels and designs of factor variables, kept as the
characteristics _fv_base and _fv_design of each variable."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.varlist import expand
from ..lang.syntax import find_top, match_options
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _base_spec(text: str) -> str:
    t = text.strip()
    for full, k in (("default", 3), ("first", 1), ("last", 1), ("frequent", 4), ("none", 1)):
        if len(t) >= k and full.startswith(t):
            return full
    if re.fullmatch(r"\d+", t):
        return t
    raise StataError(198, "invalid base specification")   # compat 2008


def _design_spec(text: str) -> str:
    t = text.strip()
    for full, k in (("default", 3), ("asbalanced", 5), ("asobserved", 5)):
        if len(t) >= k and full.startswith(t):
            return full
    raise StataError(198, f"invalid design specification '{t}'")   # VERIFY


def _set(s: "Session", names: list[str], key: str, value: str) -> None:
    for n in names:
        chars = s.data.chars.setdefault(n, {})
        if value == "default":
            chars.pop(key, None)
        else:
            chars[key] = value
    s.data.changed = True


@command("fvset")
def cmd_fvset(s: "Session", args: str) -> None:
    """fvset base|design|clear|report. VERIFY: the layout of fvset report and what fvset
    without a subcommand does (compat 2005: rc 0 and no output without data)."""
    t = args.strip()
    if not t:
        return
    sub, _, rest = t.partition(" ")
    ds = s.data
    if "report".startswith(sub) and len(sub) >= 3:
        c = find_top(rest, ",")
        body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
        o = match_options(opt, {"base": 1, "design": 1})
        names = expand(ds, body) if body.strip() else [v.name for v in ds.vars]
        rows = []
        for n in names:
            ch = ds.chars.get(n, {})
            b, d = ch.get("_fv_base", ""), ch.get("_fv_design", "")
            if o.get("base") and b != _base_spec(str(o["base"])):
                continue
            if o.get("design") and d != _design_spec(str(o["design"])):
                continue
            if b or d:
                rows.append((n, b, d))
        if rows:
            # compat 2008: each column as wide as its longest entry plus 4 (the last without
            # the 4), between rules of dashes. VERIFY: long variable names
            out = s.output
            wn = max([8] + [len(r[0]) for r in rows]) + 4
            wb = max([4] + [len(r[1]) for r in rows]) + 4
            wd = max([6] + [len(r[2]) for r in rows])
            rule = "  " + "-" * (wn + wb + wd) + "\n"
            out.write(f"\n  {'Variable':<{wn}}{'Base':<{wb}}Design\n" + rule, "text")
            for n, b, d in rows:
                out.write(f"  {n:<{wn}}{b:<{wb}}{d}".rstrip() + "\n", "result")
            out.write(rule, "text")
        return
    if sub == "clear":
        names = expand(ds, rest)
        for n in names:
            ch = ds.chars.get(n, {})
            ch.pop("_fv_base", None)
            ch.pop("_fv_design", None)
        ds.changed = True
        return
    if sub in ("base", "design"):
        spec, _, vl = rest.strip().partition(" ")
        if not vl.strip():
            raise StataError(100, "varlist required")
        names = expand(ds, vl)
        for n in names:
            if ds.get(n).is_string:
                raise StataError(109, f"{n}:  string variables may not be used as factor variables")
        if sub == "base":
            _set(s, names, "_fv_base", _base_spec(spec))
        else:
            _set(s, names, "_fv_design", _design_spec(spec))
        return
    raise StataError(198, f"unrecognized subcommand {sub}")   # VERIFY
