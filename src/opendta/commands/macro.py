"""local, global, macro e funções estendidas de macro (manual [P] macro)."""

from __future__ import annotations

import fnmatch
import re
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.formats import number_to_macro
from ..lang.words import split_words, strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_NAME = re.compile(r"^\s*([A-Za-z_0-9][A-Za-z0-9_]*)(.*)$", re.S)


def _assign(s: "Session", args: str, setter) -> None:
    inc = re.match(r"^\s*(\+\+|--)([A-Za-z_][A-Za-z0-9_]*)\s*$", args)
    if inc and setter == s.macros.set_local:
        # local ++i / local --i
        name = inc.group(2)
        current = s.macros.get_local(name) or "0"
        try:
            x = float(current)
        except ValueError:
            raise StataError(198, "invalid syntax")
        x += 1 if inc.group(1) == "++" else -1
        setter(name, number_to_macro(x))
        return
    m = _NAME.match(args)
    if not m:
        raise StataError(198, "invalid syntax")
    name, rest = m.group(1), m.group(2)
    if len(name) > 31:
        raise StataError(198, f"{name} invalid name")
    rest_s = rest.strip()
    if rest_s.startswith("="):
        v = s.eval(rest_s[1:])
        value = v if isinstance(v, str) else number_to_macro(v)
    elif rest_s.startswith(":"):
        value = extended_function(s, rest_s[1:].strip())
    else:
        value = strip_outer_quotes(rest_s) if rest_s else ""
    setter(name, value)


@command("local", "loc")
def cmd_local(s: "Session", args: str) -> None:
    _assign(s, args, s.macros.set_local)


@command("global", "gl")
def cmd_global(s: "Session", args: str) -> None:
    _assign(s, args, s.macros.set_global)


@command("macro", "ma")
def cmd_macro(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    rest = rest.strip()
    out = s.output
    if sub in ("list", "l", "li", "lis", "dir", "di"):
        names = rest.split()
        for name, value in sorted(s.macros.globals.items()):
            if names and name not in names:
                continue
            out.write(f"{name + ':':<16}{value}\n", "text")
        for name, value in sorted(s.macros.locals.items()):
            if names and "_" + name not in names:
                continue
            out.write(f"{'_' + name + ':':<16}{value}\n", "text")
        return
    if sub == "drop":
        for pat in rest.split():
            if pat == "_all":
                s.macros.clear_globals()
                continue
            target = s.macros.locals if pat.startswith("_") else s.macros.globals
            key = pat[1:] if pat.startswith("_") else pat
            for name in [k for k in target if fnmatch.fnmatchcase(k, key)]:
                del target[name]
        return
    if sub == "shift":
        loc = s.macros.locals
        k = 1
        while str(k + 1) in loc:
            loc[str(k)] = loc[str(k + 1)]
            k += 1
        loc.pop(str(k), None)
        return
    raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# Funções estendidas: local x : <função>
# ---------------------------------------------------------------------------

def extended_function(s: "Session", text: str) -> str:
    t = text.strip()

    m = re.match(r"^word\s+count\s+(.*)$", t, re.S)
    if m:
        return str(len(split_words(m.group(1))))
    m = re.match(r"^word\s+(\S+)\s+of\s+(.*)$", t, re.S)
    if m:
        k = int(float(s.eval(m.group(1))))
        words = split_words(m.group(2))
        return words[k - 1] if 1 <= k <= len(words) else ""
    m = re.match(r"^(?:length|strlen)\s+(local|global)\s+(\S+)$", t)
    if m:
        store = s.macros.get_local if m.group(1) == "local" else s.macros.get_global
        return str(len(store(m.group(2))))
    m = re.match(r"^di(?:s|sp|spl|spla|splay)?\s+(.*)$", t, re.S)
    if m:
        return _capture_display(s, m.group(1))
    m = re.match(r"^list\s+(.*)$", t, re.S)
    if m:
        return _list_function(s, m.group(1).strip())
    m = re.match(r"^env(?:ironment)?\s+(\S+)$", t)
    if m:
        import os
        return os.environ.get(m.group(1), "")
    m = re.match(r"^piece\s+(\S+)\s+(\S+)\s+of\s+(.*)$", t, re.S)
    if m:
        return _piece(int(float(m.group(1))), int(float(m.group(2))), strip_outer_quotes(m.group(3)))
    raise StataError(198, f"função estendida ainda não implementada: {t.split()[0] if t else ''}")


def _capture_display(s: "Session", args: str) -> str:
    from ..core.output import Capture
    from .display import cmd_display
    cap = Capture()
    out = s.output
    saved = out._listeners, out.column, out.quiet_depth, out.noisy_depth
    out._listeners = [cap]
    out.quiet_depth = out.noisy_depth = 0
    try:
        cmd_display(s, args + " _continue")
    finally:
        out._listeners, out.column, out.quiet_depth, out.noisy_depth = saved
    return cap.text


def _piece(n: int, length: int, text: str) -> str:
    words = text.split()
    pieces: list[str] = []
    cur = ""
    for w in words:
        cand = (cur + " " + w).strip()
        if len(cand) <= length or not cur:
            cur = cand
        else:
            pieces.append(cur)
            cur = w
    if cur:
        pieces.append(cur)
    return pieces[n - 1] if 1 <= n <= len(pieces) else ""


def _list_function(s: "Session", spec: str) -> str:
    get = s.macros.get_local
    words = spec.split()
    if not words:
        raise StataError(198, "invalid syntax")
    op = words[0]
    if op in ("uniq", "dups", "sort", "clean", "retokenize", "sizeof"):
        items = split_words(get(words[1]), keep_quotes=True)
        if op == "uniq":
            seen: list[str] = []
            for it in items:
                if it not in seen:
                    seen.append(it)
            return " ".join(seen)
        if op == "dups":
            seen, dups = [], []
            for it in items:
                if it in seen:
                    dups.append(it)
                seen.append(it)
            return " ".join(dups)
        if op == "sort":
            return " ".join(sorted(items))
        if op in ("clean", "retokenize"):
            return " ".join(items)
        return str(len(items))
    if op == "posof":
        m = re.match(r'^posof\s+(.+?)\s+in\s+(\S+)$', spec)
        if not m:
            raise StataError(198, "invalid syntax")
        target = strip_outer_quotes(m.group(1))
        items = split_words(get(m.group(2)))
        return str(items.index(target) + 1) if target in items else "0"
    if len(words) == 3 and words[1] in ("|", "&", "-", "==", "===", "in"):
        a = split_words(get(words[0]), keep_quotes=True)
        b = split_words(get(words[2]), keep_quotes=True)
        o = words[1]
        if o == "|":
            return " ".join(a + [x for x in b if x not in a])
        if o == "&":
            return " ".join(x for x in a if x in b)
        if o == "-":
            return " ".join(x for x in a if x not in b)
        if o == "==":
            return "1" if sorted(a) == sorted(b) else "0"
        if o == "===":
            return "1" if a == b else "0"
        return "1" if all(x in b for x in a) else "0"
    raise StataError(198, "invalid syntax")
