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
            vnames, lbl = words, ""          # no name: removes the label
        else:
            *vnames, lbl = words
        if lbl == ".":
            lbl = ""
        for n in unique(expand(ds, " ".join(vnames), allow_empty=False)):
            var = ds.get(n)
            if var.is_string:
                raise StataError(181, "may not label strings")
            var.value_label = lbl
            _widen(var, ds.value_labels.get(lbl))
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


@command("numlabel")
def cmd_numlabel(s: "Session", args: str) -> None:
    """numlabel [lblnamelist], add|remove [mask(str) force] ([D] labelbook): puts the value in
    front of each value label, or takes it away. A label whose texts all have the prefix
    already is left as it is, even with force (compat 2008). VERIFY: what force changes."""
    comma = find_top(args, ",")
    body, opt = (args, "") if comma == -1 else (args[:comma], args[comma + 1:])
    o = match_options(opt, {"add": 3, "remove": 6, "mask": 4, "force": 5, "detail": 3})
    if bool(o.get("add")) == bool(o.get("remove")):
        if o.get("add"):
            raise StataError(198, "options add and remove must not be both specified")
        raise StataError(198, "option add or remove should be specified")   # compat 2005
    if o.get("detail"):
        raise StataError(198, "option detail is not available yet in OpenDTA")
    mask = strip_outer_quotes(str(o["mask"])) if o.get("mask") else "#. "
    if "#" not in mask:
        raise StataError(198, "mask() must contain #")   # VERIFY
    ds = s.data
    names = split_words(body)
    if not names or names == ["_all"]:
        names = list(ds.value_labels)
    for n in names:
        if n not in ds.value_labels:
            raise StataError(111, f"value label {n} not found")
    changed = False
    for n in names:
        lab = ds.value_labels[n]
        has = {v: text.startswith(mask.replace("#", str(v))) for v, text in lab.items()}
        if o.get("add") and all(has.values()) or o.get("remove") and not any(has.values()):
            continue
        for v, text in lab.items():
            prefix = mask.replace("#", str(v))
            if o.get("add") and not has[v]:
                lab[v] = (prefix + text)[:32000]
            elif o.get("remove") and has[v]:
                lab[v] = text[len(prefix):]
        changed = True
    if not changed:
        s.output.write("(no value label to be modified)\n", "text")
        return
    ds.changed = True
    s.notify_state()


def _widen(var, labels) -> None:
    """label values widens a %#.#g display format to the length of the longest label
    (compat 2007 and 2008: %9.0g becomes %46.0g). VERIFY: other formats and labels defined
    after label values."""
    if not labels:
        return
    m = re.fullmatch(r"%(-?)(\d+)\.(\d+)g", var.fmt)
    longest = max(len(v) for v in labels.values())
    if m and int(m.group(2)) < longest:
        var.fmt = f"%{m.group(1)}{longest}.{m.group(3)}g"
