"""`class exit` ([P] class exit) and `classutil` ([P] classutil).

    class exit [rvalue]
    classutil drop instance [instance ...]
    classutil describe object [, recurse newok]
    classutil dir [pattern] [, all detail]
    classutil cdir [pattern]
    classutil which classname [, all]

`cutil` is a synonym of classutil. The object model is in lang/classes.py and lang/classeval.py.
The layout of classutil describe, dir and cdir follows Stata 14's output (compat/do/1928_classes.do).
"""

from __future__ import annotations

import fnmatch
import re
from typing import TYPE_CHECKING

from ..core.errors import ExitRequest, StataError
from ..lang.classes import ClassDef, Instance, system, type_of
from ..lang.classeval import evaluator, parse_path
from ..lang.words import split_words
from ..lang.syntax import match_options
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("class")
def cmd_class(s: "Session", args: str) -> None:
    word, _, rest = args.strip().partition(" ")
    if word != "exit":
        raise StataError(198, "invalid syntax")
    cs = system(s)
    if not cs.ctx:
        raise StataError(198, "command -class exit- can only be used in a class program")
    frame = cs.ctx[-1]
    if rest.strip():
        frame.retval = evaluator(s).rvalue(rest)
    raise ExitRequest(0)


def _dotted(name: str) -> str:
    name = name.strip()
    return name if name.startswith(".") else "." + name


def _drop(s: "Session", rest: str) -> None:
    cs = system(s)
    names = rest.split()
    if not names:
        raise StataError(198, "invalid syntax")
    for n in names:
        n = n.lstrip(".")
        if "." in n or "[" in n:
            raise StataError(459, f"{n}:  can only drop top-level objects")
        if n not in cs.globals:
            raise StataError(111, f".{n}:  object not found")
        cs.drop(n)


def _describe(s: "Session", rest: str) -> None:
    head, comma, opts = rest.partition(",")
    o = match_options(opts, {"recurse": 7, "newok": 5}) if comma else {}
    target = head.strip()
    if not target or len(target.split()) != 1:
        raise StataError(198, "invalid syntax")
    ev = evaluator(s)
    segs, end = parse_path(_dotted(target))
    if not segs or _dotted(target)[end:].strip():
        raise StataError(198, "invalid syntax")
    r = ev.resolve(segs, global_context=True)
    if r[0] == "class":
        value, name = r[1], target.lstrip(".")
    elif r[0] == "var":
        value, name = r[1].value, _dotted(target)
    else:
        raise StataError(198, f"{_dotted(target)} is not an object")       # VERIFY
    if isinstance(value, (Instance, ClassDef)):
        s.output.write("\n")
    _describe_value(s, name, value, bool(o.get("recurse")), set())
    if isinstance(value, (Instance, ClassDef)):
        typ = value.cls.name if isinstance(value, Instance) else value.name
        bitype = "class" if isinstance(value, ClassDef) else "instance"
    else:
        typ = bitype = type_of(value)
    s.r = {"type": typ, "bitype": bitype}


def _line(typ: str, name: str, v: object) -> str:
    """One member: type, name and, for doubles and strings, the value."""
    out = f"{typ:<12} {name}"
    if isinstance(v, (float, str)):
        out = f"{typ:<12} {name:<14}{_short(v)}"
    return out


def _describe_value(s: "Session", name: str, value: object, recurse: bool, seen: set) -> None:
    out = s.output
    if isinstance(value, ClassDef):
        out.write(f"class {value.name} {{\n")
        for scope in ("classwide", "instancespecific"):
            decls = [d for d in value.effective().values() if d.scope == scope]
            if decls:
                out.write(f"   {scope}:\n")
                for d in decls:
                    out.write(f"       {d.type:<13} {d.name}\n")
        out.write("}\n")
        return
    if not isinstance(value, Instance):
        out.write(_line(type_of(value), name, value) + "\n")
        return
    # VERIFY: classwide and dynamic members of an instance; arrays
    out.write(f"{value.cls.name} {name}:\n")
    nested = []
    names = [n for n, d in value.cls.effective().items() if d.scope != "classwide"]
    for n in names + value.dynamic:
        hit = value.member(n)
        if hit is None:
            continue
        v = hit[1].value
        decl = value.cls.effective().get(n)
        typ = decl.type if decl is not None else type_of(v)
        out.write("        " + _line(typ, "." + n, v) + "\n")
        if isinstance(v, Instance):
            nested.append((f"{name}.{n}", v))
    if recurse:
        for n, v in nested:
            if id(v) not in seen:
                seen.add(id(v))
                out.write("\n")
                _describe_value(s, n, v, recurse, seen)


def _short(v: object) -> str:
    if isinstance(v, float):
        from ..core.formats import number_to_macro
        return f" = {number_to_macro(v)}"
    if isinstance(v, str):
        return f' = "{v}"'
    return ""


def _columns(s: "Session", names: list[str]) -> None:
    """Names side by side, each column two wider than the longest name."""
    if not names:
        return
    w = max(len(n) for n in names) + 2
    per = max(1, (int(float(s.settings.get("linesize", 80))) + 2) // w)   # VERIFY: where Stata wraps
    for i in range(0, len(names), per):
        s.output.write("".join(f"{n:<{w}}" for n in names[i:i + per]).rstrip() + "\n")


def _dir(s: "Session", rest: str) -> None:
    head, comma, opts = rest.partition(",")
    o = match_options(opts, {"all": 3, "detail": 3}) if comma else {}
    pattern = head.strip() or "*"
    cs = system(s)
    names = [("." + n, c.value) for n, c in cs.globals.items()
             if fnmatch.fnmatchcase(n, pattern.lstrip("."))]
    if o.get("all"):
        names += [(n, d) for n, d in cs.defs.items() if fnmatch.fnmatchcase(n, pattern)]
    if o.get("detail"):
        for n, v in names:
            kind = "class" if isinstance(v, ClassDef) else type_of(v)
            s.output.write(_line(kind, n, v) + "\n")
    else:
        _columns(s, [n for n, _ in names])
    s.r = {"list": " ".join(n for n, _ in names)}


def _class_files(s: "Session") -> list[str]:
    from ..lang.adopath import adopath, resolve_entry
    found: list[str] = []
    for entry in adopath(s):
        folder = resolve_entry(s, entry)
        if not folder.is_dir():
            continue
        for f in sorted(folder.glob("*.class")) + sorted(folder.glob("?/*.class")):
            if f.stem not in found:
                found.append(f.stem)
    return found


def _cdir(s: "Session", rest: str) -> None:
    pattern = rest.strip() or "*"
    names = sorted(n for n in _class_files(s) if fnmatch.fnmatchcase(n, pattern))
    _columns(s, names)
    s.r = {"list": " ".join(names)}


def _which(s: "Session", rest: str) -> None:
    from ..lang.adopath import find_file
    head, comma, opts = rest.partition(",")
    o = match_options(opts, {"all": 3}) if comma else {}
    words = split_words(head)
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    fname = f"{words[0].lstrip('.')}.class"
    hits = find_file(s, fname, all_=bool(o.get("all")))
    if not hits:
        raise StataError(601, f'file "{fname}" not found')
    for h in hits:
        s.output.write(f"{h}\n")
        for line in h.read_text(encoding="utf-8", errors="replace").splitlines():
            if line.startswith("*!"):
                s.output.write(line + "\n")
    s.r = {"fn": " ".join(f'"{h}"' for h in hits) if o.get("all") else str(hits[0])}


_SUBCOMMANDS = {"drop": _drop, "describe": _describe, "dir": _dir, "cdir": _cdir,
                "which": _which}


@command("classutil")
def cmd_classutil(s: "Session", args: str) -> None:
    m = re.match(r"\s*([A-Za-z]+)(.*)$", args, re.S)
    word, rest = (m.group(1), m.group(2)) if m else ("", "")
    # VERIFY: abbreviations of the subcommands (d for describe?)
    fn = _SUBCOMMANDS.get(word)
    if fn is None:
        raise StataError(198, "invalid syntax")
    fn(s, rest)


@command("cutil")
def cmd_cutil(s: "Session", args: str) -> None:
    cmd_classutil(s, args)
