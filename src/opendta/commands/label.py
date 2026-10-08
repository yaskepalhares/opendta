"""label — manual [D] label."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.varlist import expand, resolve_name, unique
from ..lang.syntax import find_top, match_options
from ..lang.words import split_words, strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("label", "la")
def cmd_label(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    rest = rest.strip()
    ds = s.data
    out = s.output

    if "variable".startswith(sub) and len(sub) >= 3:
        name, _, text = rest.partition(" ")
        var = ds.get(resolve_name(ds, name))
        text = strip_outer_quotes(text) if text.strip() else ""
        if len(text) > 80:
            text = text[:80]
        var.label = text
        ds.changed = True
    elif "data".startswith(sub) and len(sub) >= 1:
        ds.label = strip_outer_quotes(rest) if rest else ""
        ds.changed = True
    elif "define".startswith(sub) and len(sub) >= 3:
        _define(s, rest)
    elif "values".startswith(sub) and len(sub) >= 3:
        words = rest.split()
        if not words:
            raise StataError(100, "varlist required")
        if len(words) == 1:
            vnames, lbl = words, ""          # sem nome: remove o rótulo
        else:
            *vnames, lbl = words
        if lbl == ".":
            lbl = ""
        for n in unique(expand(ds, " ".join(vnames), allow_empty=False)):
            var = ds.get(n)
            if var.is_string:
                raise StataError(181, "may not label strings")
            var.value_label = lbl
        ds.changed = True
    elif "list".startswith(sub) and len(sub) >= 1:
        names = rest.split() or sorted(ds.value_labels)
        for name in names:
            if name not in ds.value_labels:
                raise StataError(111, f"value label {name} not found")
            out.write(f"{name}:\n", "text")
            for k, v in sorted(ds.value_labels[name].items()):
                out.write(f"{k:>12} {v}\n", "text")
    elif sub == "dir":
        for name in sorted(ds.value_labels):
            out.write(name + "\n", "text")
    elif sub == "drop":
        names = rest.split()
        if names == ["_all"]:
            ds.value_labels.clear()
        else:
            for name in names:
                if name not in ds.value_labels:
                    raise StataError(111, f"value label {name} not found")
                del ds.value_labels[name]
        ds.changed = True
    else:
        raise StataError(198, "invalid syntax")
    s.notify_state()


def _define(s: "Session", text: str) -> None:
    ds = s.data
    comma = find_top(text, ",")
    opts = {}
    if comma != -1:
        opts = match_options(text[comma + 1:], {"add": 3, "modify": 3, "replace": 7, "nofix": 5})
        text = text[:comma]
    name, _, body = text.strip().partition(" ")
    if not re.match(r"^[A-Za-z_][A-Za-z0-9_]{0,31}$", name):
        raise StataError(198, f"{name} invalid name")
    tokens = split_words(body)
    if len(tokens) % 2:
        raise StataError(198, "invalid syntax")
    pairs: dict[int, str] = {}
    for k in range(0, len(tokens), 2):
        try:
            val = float(tokens[k])
        except ValueError:
            raise StataError(198, f"{tokens[k]} not an integer")
        if val != int(val):
            raise StataError(198, "may not label non-integers")
        pairs[int(val)] = tokens[k + 1]
    exists = name in ds.value_labels
    if exists and not (opts.get("add") or opts.get("modify") or opts.get("replace")):
        raise StataError(110, f"label {name} already defined")
    if opts.get("replace") or not exists:
        ds.value_labels[name] = {}
    current = ds.value_labels[name]
    for k, v in pairs.items():
        if k in current and opts.get("add") and not opts.get("modify"):
            raise StataError(180, f"invalid attempt to modify label")
        current[k] = v
    ds.changed = True
