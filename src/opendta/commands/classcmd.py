"""`class exit` ([P] class exit) and `classutil` ([P] classutil).

    class exit [rvalue]
    classutil drop instance [instance ...]
    classutil describe object [, recurse newok]
    classutil dir [pattern] [, all detail]
    classutil cdir [pattern]
    classutil which classname [, all]

`cutil` is a synonym of classutil. The object model is in lang/classes.py and lang/classeval.py.
VERIFY: the manual does not show the output of classutil describe and classutil dir; the layout
here is OpenDTA's own until it is checked against Stata (compat/do/1928_classes.do).
"""

from __future__ import annotations

import fnmatch
import re
from typing import TYPE_CHECKING

from ..core.errors import ExitRequest, StataError
from ..lang.classes import Array, ClassDef, Instance, system, type_of
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
        # VERIFY: Stata's message outside a member program
        raise StataError(198, "class exit may be used only from class-member programs")
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
            raise StataError(198, f".{n} is not a top-level instance")     # VERIFY
        if n not in cs.globals:
            raise StataError(111, f".{n} not found")                         # VERIFY
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
    out = s.output
    out.write("\n")
    _describe_value(s, name, value, bool(o.get("recurse")), set())
    typ = "classname" if isinstance(value, (Instance, ClassDef)) else type_of(value)
    bitype = ("class" if isinstance(value, ClassDef) else "instance") if typ == "classname" else typ
    s.r = {"type": typ, "bitype": bitype}


def _describe_value(s: "Session", name: str, value: object, recurse: bool, seen: set) -> None:
    out = s.output
    if isinstance(value, ClassDef):
        cls, inst = value, None
        out.write(f"{name} is a class\n")
    elif isinstance(value, Instance):
        cls, inst = value.cls, value
        out.write(f"{name} is an instance of {cls.name}\n")
    else:
        out.write(f"{name} is {_article(type_of(value))}{_short(value)}\n")
        return
    if cls.parents:
        out.write(f"  inherits from {' '.join(p.name for p in cls.parents)}\n")
    nested = []
    for scope in ("classwide", "instancespecific"):
        decls = [d for d in cls.effective().values() if d.scope == scope]
        if not decls:
            continue
        out.write(f"  {scope}:\n")
        for d in decls:
            v = d.value
            if inst is not None and inst.member(d.name):
                v = inst.member(d.name)[1].value
            out.write(f"    {d.type:<12} .{d.name}{_short(v)}\n")
            if isinstance(v, Instance):
                nested.append((f"{name}.{d.name}", v))
    if inst is not None and inst.dynamic:
        out.write("  dynamic:\n")
        for n in inst.dynamic:
            v = inst.members[n].value
            out.write(f"    {type_of(v):<12} .{n}{_short(v)}\n")
            if isinstance(v, Instance):
                nested.append((f"{name}.{n}", v))
    programs = sorted({p for c in cls.mro() for p in c.programs})
    if programs:
        out.write("  programs: " + " ".join("." + p for p in programs) + "\n")
    if recurse:
        for n, v in nested:
            if id(v) not in seen:
                seen.add(id(v))
                out.write("\n")
                _describe_value(s, n, v, recurse, seen)


def _article(t: str) -> str:
    return ("an " if t[:1] in "aeiou" else "a ") + t


def _short(v: object) -> str:
    if isinstance(v, float):
        from ..core.formats import number_to_macro
        return f" = {number_to_macro(v)}"
    if isinstance(v, str):
        return f' = "{v}"'
    if isinstance(v, Array):
        return f" ({len(v.items)} element{'s' if len(v.items) != 1 else ''})"
    return ""


def _dir(s: "Session", rest: str) -> None:
    head, comma, opts = rest.partition(",")
    o = match_options(opts, {"all": 3, "detail": 3}) if comma else {}
    pattern = head.strip() or "*"
    cs = system(s)
    names = [("." + n, c.value) for n, c in cs.globals.items()
             if fnmatch.fnmatchcase(n, pattern.lstrip("."))]
    if o.get("all"):
        names += [(n, d) for n, d in cs.defs.items() if fnmatch.fnmatchcase(n, pattern)]
    for n, v in names:
        if o.get("detail"):
            kind = "class" if isinstance(v, ClassDef) else type_of(v)
            s.output.write(f"  {n:<20} {kind}{_short(v)}\n")
        else:
            s.output.write(f"  {n}\n")
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
    for n in names:
        s.output.write(f"  {n}\n")
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
