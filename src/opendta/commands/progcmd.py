"""program drop/dir/list and return, ereturn, sreturn ([P] program, [P] return).

r() of an rclass program only changes when it ends without an error: `return`
writes to an accumulator of the scope, published at the end (lang/programs.py).
e() changes immediately, as in Stata, and s() follows the same rule as r() for sclass.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.programs import current_program_scope
from .macro import _assign
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("program", "pr")
def cmd_program(s: "Session", args: str) -> None:
    words = args.split()
    if not words:
        raise StataError(198, "invalid syntax")
    sub, rest = words[0], words[1:]
    if sub == "drop":
        if not rest:
            raise StataError(198, "invalid syntax")
        if rest in (["_all"], ["_allado"]):
            if rest == ["_all"]:
                s.programs.clear()
            else:
                for k in [k for k, p in s.programs.items() if p.source.endswith(".ado")]:
                    del s.programs[k]
            return
        for name in rest:
            if name not in s.programs:
                raise StataError(111, f"program {name} not found")   # VERIFY
            del s.programs[name]
        return
    if sub == "dir":
        out = s.output
        # VERIFY: layout of program dir
        for name, p in s.programs.items():
            size = sum(len(ln.text) + 1 for ln in p.lines)
            out.write(f"  {'ado' if p.source.endswith('.ado') else '   '} {size:>8,}  {name}\n", "text")
        return
    if sub in ("list", "li", "l"):
        names = rest or list(s.programs)
        if names == ["_all"]:
            names = list(s.programs)
        out = s.output
        for name in names:
            p = s.programs.get(name)
            if p is None:
                raise StataError(111, f"program {name} not found")
            # VERIFY: layout of program list
            out.write(f"\n{name}{', ' + p.kind if p.kind != 'nclass' else ''}:\n", "text")
            for k, ln in enumerate(p.lines, start=1):
                out.write(f"  {k}.  {ln.text}\n", "text")
        return
    # a prefix (capture, quietly...) sends `program name, plugin` here
    from ..lang.programs import define_plugin, parse_definition
    name, options = parse_definition(args)
    if "plugin" in options:
        define_plugin(s, name, options)
    raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# return / sreturn
# ---------------------------------------------------------------------------

def _scalar_value(s: "Session", text: str) -> float:
    m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$", text, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    v = s.eval(m.group(2))
    if isinstance(v, str):
        raise StataError(109, "type mismatch")
    return m.group(1), float(v)


def _pending(s: "Session", kind: str) -> dict:
    sc = current_program_scope(s)
    if sc is None or sc.program is None or sc.program.kind != kind:
        letter = {"rclass": "r", "sclass": "s"}[kind]
        raise StataError(151, f"non {letter}-class program may not set {letter}()")   # VERIFY
    return sc.pending


def _store_cmd(s: "Session", args: str, kind: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    if sub in ("list", "li", "l"):
        _results_list(s, s.r if kind == "rclass" else s.sret, "r" if kind == "rclass" else "s")
        return
    if sub == "clear":
        if kind == "rclass":
            sc = current_program_scope(s)
            if sc is not None and sc.program is not None and sc.program.kind == "rclass":
                sc.pending.clear()
            s.r = {}
        else:
            s.sret = {}
            sc = current_program_scope(s)
            if sc is not None and sc.program is not None and sc.program.kind == "sclass":
                sc.pending.clear()
        return
    store = _pending(s, kind)
    if sub in ("local", "loc"):
        _assign(s, rest, lambda name, value: store.__setitem__(name, value))
        return
    if sub in ("scalar", "sca") and kind == "rclass":
        name, value = _scalar_value(s, rest)
        store[name] = value
        return
    if sub == "add" and kind == "rclass":
        for k, v in s.r.items():
            store.setdefault(k, v)
        return
    if sub in ("matrix", "mat") and kind == "rclass":
        from .matrix import take_matrix
        store_matrix(store, take_matrix(s, rest))
        return
    raise StataError(198, "invalid syntax")


def store_matrix(store: dict, item) -> None:
    name, mat = item
    store[name] = mat


@command("return", "ret")
def cmd_return(s: "Session", args: str) -> None:
    _store_cmd(s, args, "rclass")


@command("sreturn", "sret")
def cmd_sreturn(s: "Session", args: str) -> None:
    _store_cmd(s, args, "sclass")


@command("ereturn", "eret")
def cmd_ereturn(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    if sub in ("list", "li", "l"):
        _results_list(s, s.e, "e")
        return
    if sub == "clear":
        s.e = {}
        return
    if sub in ("local", "loc"):
        _assign(s, rest, lambda name, value: s.e.__setitem__(name, value) if value
                else s.e.pop(name, None))
        return
    if sub in ("scalar", "sca"):
        name, value = _scalar_value(s, rest)
        s.e[name] = value
        return
    if sub in ("matrix", "mat", "post", "repost", "display", "di"):
        from .matrix import ereturn_matrix
        ereturn_matrix(s, sub, rest)
        return
    raise StataError(198, "invalid syntax")


def _results_list(s: "Session", store: dict, prefix: str) -> None:
    """return list / ereturn list / sreturn list."""
    from .matrix import Matrix
    out = s.output
    scalars = {k: v for k, v in store.items() if isinstance(v, float) or isinstance(v, int)}
    macros = {k: v for k, v in store.items() if isinstance(v, str)}
    mats = {k: v for k, v in store.items() if isinstance(v, Matrix)}
    # Stata lists from the last stored result to the first (observed in
    # compat/expected/0201_programas.log) and aligns the names in 22 columns
    scalars = dict(reversed(list(scalars.items())))
    macros = dict(reversed(list(macros.items())))
    mats = dict(reversed(list(mats.items())))
    if scalars:
        out.write("\nscalars:\n", "text")
        for k, v in scalars.items():
            out.write(f"{prefix + '(' + k + ')':>22} =  ", "text")
            # values with up to 16 digits (%18.0g), as in Stata 14 (compat 0501, 0504)
            out.write(format_value(float(v), "%18.0g", pad=False).strip() + "\n", "result")
    if macros:
        out.write("\nmacros:\n", "text")
        for k, v in macros.items():
            out.write(f"{prefix + '(' + k + ')':>22} : ", "text")
            out.write(f'"{v}"\n', "result")
    if mats:
        out.write("\nmatrices:\n", "text")
        for k, v in mats.items():
            out.write(f"{prefix + '(' + k + ')':>22} :  ", "text")
            out.write(f"{v.rows} x {v.cols}\n", "result")
    if prefix == "e" and s.data.esample is not None and ("b" in store or getattr(getattr(s, "_estimates", None), "no_b", False)):
        # VERIFY spaces after e(sample)
        out.write("\nfunctions:\n", "text")
        out.write(f"{'e(sample)':>22}\n", "text")
