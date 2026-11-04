"""Stata's class system ([P] class): the object model and the .class files.

    begin coordinate.class
    version 14
    class coordinate {
        [classwide: | instancespecific:]
        double x [= rvalue]
        string s [= rvalue]
        array  a [= {el, el, ...}]
        othercls o [= rvalue]
        name = rvalue                   (type taken from the value)
    } [, inherit(parent ...)]
    program .set
        ...
    end
    end coordinate.class

Every member value lives in a Cell (the manual's "index card"). Identifiers that share an object
through `.a.ref = .b.ref` point to the same Cell; `nref` counts them, and when it reaches 0 the
value is released (instances run their .destructor). Classwide variables live in the ClassDef that
declares them while instances of it exist. Initial values are computed once, when the definition is
read, in the global context ("prerecorded"), and each new instance gets a copy of them.

The statements, the rvalues and the references are in classeval.py; the commands `class exit` and
`classutil` in commands/classcmd.py.
"""

from __future__ import annotations

import itertools
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.missing import SYSMISS
from .lexer import split_commands
from .programs import Program, make_program, parse_definition

if TYPE_CHECKING:
    from ..session import Session

_ids = itertools.count(1)
_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]{0,31}")
PRIMITIVE = ("double", "string", "array")


class Cell:
    """One member value; the identifiers that share it (.ref) point to the same Cell."""
    __slots__ = ("value", "nref", "uid", "dead")

    def __init__(self, value: object, nref: int = 1):
        self.value = value
        self.nref = nref
        self.uid = next(_ids)
        self.dead = False


class Array:
    """Arrays are sparse and untyped: each element has its own type."""

    def __init__(self, items: dict[int, Cell] | None = None):
        self.items: dict[int, Cell] = items if items is not None else {}

    def nels(self) -> int:
        return max(self.items, default=0)


@dataclass(eq=False)
class Decl:
    scope: str                  # classwide | instancespecific
    type: str                   # double | string | array | a class name
    name: str
    value: object               # prerecorded initial value
    owner: "ClassDef | None" = None


@dataclass(eq=False)
class ClassDef:
    name: str
    path: str = ""
    parents: list["ClassDef"] = field(default_factory=list)
    decls: list[Decl] = field(default_factory=list)
    programs: dict[str, Program] = field(default_factory=dict)
    classwide: dict[str, Cell] = field(default_factory=dict)
    ninst: int = 0              # instances of this class or of a class that inherits from it

    def mro(self) -> list["ClassDef"]:
        """The class, then its parents; with multiple inheritance the rightmost parent wins."""
        out = [self]
        for p in reversed(self.parents):
            for c in p.mro():
                if c not in out:
                    out.append(c)
        return out

    def effective(self) -> dict[str, Decl]:
        """Member variables, the inherited ones first; the class's own definition wins."""
        out: dict[str, Decl] = {}
        for p in self.parents:
            out.update(p.effective())
        for d in self.decls:
            out[d.name] = d
        return out

    def find_program(self, name: str) -> tuple["ClassDef", Program] | None:
        for c in self.mro():
            if name in c.programs:
                return c, c.programs[name]
        return None

    def isof(self, name: str) -> bool:
        return any(c.name == name for c in self.mro())


class Instance:
    __slots__ = ("cls", "members", "dynamic")

    def __init__(self, cls: ClassDef):
        self.cls = cls
        self.members: dict[str, Cell] = {}       # instance-specific and dynamic
        self.dynamic: list[str] = []

    def member(self, name: str) -> tuple[dict, Cell] | None:
        if name in self.members:
            return self.members, self.members[name]
        d = self.cls.effective().get(name)
        if d is not None and d.scope == "classwide" and name in d.owner.classwide:
            return d.owner.classwide, d.owner.classwide[name]
        return None


def type_of(v: object) -> str:
    """.objtype: double, string, array or the class name."""
    if isinstance(v, float):
        return "double"
    if isinstance(v, str):
        return "string"
    if isinstance(v, Array):
        return "array"
    if isinstance(v, Instance):
        return v.cls.name
    return v.name           # a class definition     VERIFY


def same_type(old: object, new: object) -> bool:
    if isinstance(old, Instance):
        return isinstance(new, Instance) and new.cls.isof(old.cls.name)
    return type(old) is type(new) and not isinstance(old, ClassDef)


def mismatch() -> StataError:
    return StataError(109, "type mismatch")


@dataclass(eq=False)
class Frame:
    """A running member program: its implied context and the class that defines the program."""
    cell: Cell
    owner: ClassDef
    retval: object = None       # what `class exit` returned (see classeval.RV)
    oncopy_src: str = ""


class ClassSystem:
    def __init__(self, s: "Session"):
        self.s = s
        self.defs: dict[str, ClassDef] = {}
        self.globals: dict[str, Cell] = {}
        self.keys: dict[str, Cell] = {}
        self.ctx: list[Frame] = []
        self._loading: set[str] = set()

    # -- class definitions ----------------------------------------------------------------------
    def load(self, name: str) -> ClassDef:
        if name in self.defs:
            return self.defs[name]
        from .adopath import find_file
        hits = find_file(self.s, f"{name}.class")
        if not hits:
            raise StataError(111, f"class {name} not found")          # VERIFY
        if name in self._loading:
            raise StataError(198, f"class {name} refers to itself")   # VERIFY
        self._loading.add(name)
        try:
            cls = self._parse(name, str(hits[0]), hits[0].read_text(encoding="utf-8",
                                                                     errors="replace"))
        finally:
            self._loading.discard(name)
        self.defs[name] = cls
        return cls

    def _parse(self, name: str, path: str, source: str) -> ClassDef:
        lines = [ln for ln in split_commands(source) if ln.kind == "cmd" and ln.text.strip()]
        cls = ClassDef(name, path)
        body: list[str] | None = None
        i, seen_class = 0, False
        while i < len(lines):
            t = lines[i].text.strip()
            i += 1
            word, _, rest = t.partition(" ")
            if body is not None:
                if t.startswith("}"):
                    self._header_end(cls, t[1:])
                    self._declarations(cls, body)
                    body = None
                else:
                    body.append(t)
                continue
            if word == "version":
                continue
            if word == "class" and not seen_class:
                seen_class = True
                head, brace, after = rest.partition("{")
                if not brace:
                    raise StataError(198, "class: { expected")        # VERIFY
                cname = head.strip()
                if cname and cname != name:
                    raise StataError(198, f"class {cname} defined in {name}.class")   # VERIFY
                after = after.strip()
                if after.startswith("}"):
                    self._header_end(cls, after[1:])
                else:
                    body = [after] if after else []
                continue
            if len(word) >= 3 and "program".startswith(word):
                j = i
                while j < len(lines) and lines[j].text.strip() != "end":
                    j += 1
                if j >= len(lines):
                    raise StataError(198, "program define: end not found")    # VERIFY
                pname, opts = self._program_header(rest)
                cls.programs[pname] = make_program("." + pname, opts, list(lines[i:j]), path)
                i = j + 1
                continue
            raise StataError(198, f"{t}: invalid line in {name}.class")         # VERIFY
        if body is not None:
            raise StataError(198, "class: } expected")                         # VERIFY
        return cls

    @staticmethod
    def _program_header(rest: str) -> tuple[str, dict[str, str]]:
        words = rest.split(",")[0].split()
        if words and len(words[0]) >= 2 and "define".startswith(words[0]):
            rest = rest.split(None, 1)[1] if len(rest.split(None, 1)) > 1 else ""
        rest = rest.strip()
        if not rest.startswith("."):
            raise StataError(198, "member programs are defined as program .name")   # VERIFY
        return parse_definition(rest[1:])

    def _header_end(self, cls: ClassDef, tail: str) -> None:
        tail = tail.strip()
        if not tail:
            return
        m = re.fullmatch(r",\s*inherit\s*\(([^)]*)\)\s*", tail)
        if not m:
            raise StataError(198, "invalid syntax")
        cls.parents = [self.load(p) for p in m.group(1).split()]

    def _declarations(self, cls: ClassDef, body: list[str]) -> None:
        scope = "instancespecific"
        for t in body:
            if t in ("classwide:", "instancespecific:"):
                scope = t[:-1]
                continue
            d = self.declaration(t, scope, cls.name)
            d.owner = cls
            cls.decls.append(d)

    def declaration(self, text: str, scope: str, being_defined: str = "") -> Decl:
        """[type] name [= rvalue], in a class definition or in .Declare."""
        from .classeval import evaluator
        lhs, eq, rhs = text.partition("=")
        words = lhs.split()
        if len(words) == 1 and eq:
            typ, name = "", words[0]
        elif len(words) == 2:
            typ, name = words
        else:
            raise StataError(198, "invalid syntax")
        if not _NAME.fullmatch(name) or (typ and not _NAME.fullmatch(typ)):
            raise StataError(198, f"{name} invalid name")
        if typ and typ == being_defined:
            raise StataError(198, f"class {typ} cannot contain itself")   # VERIFY
        if eq:
            value = evaluator(self.s).value_of(rhs, global_context=True)
            if typ:
                default = self._default(typ)
                ok = same_type(default, value)
                self.release_value(default)
                if not ok:
                    self.release_value(value)
                    raise mismatch()
            else:
                typ = type_of(value)
        else:
            value = self._default(typ)
        return Decl(scope, typ, name, value)

    def _default(self, typ: str) -> object:
        if typ == "double":
            return SYSMISS
        if typ == "string":
            return ""
        if typ == "array":
            return Array()
        return self.new_instance(self.load(typ), "")

    # -- instances --------------------------------------------------------------------------------
    def new_instance(self, cls: ClassDef, args: str) -> Instance:
        inst = Instance(cls)
        self._count(cls, 1)
        for d in cls.effective().values():
            if d.scope == "classwide":
                if d.name not in d.owner.classwide:
                    d.owner.classwide[d.name] = Cell(self.copy_value(d.value, oncopy=False))
            else:
                inst.members[d.name] = Cell(self.copy_value(d.value, oncopy=False))
        found = cls.find_program("new")
        if found is not None:
            try:
                self.drop_rv(self.call(Cell(inst, 0), found[0], found[1], args))
            except BaseException:
                self.destroy(inst)
                raise
        return inst

    def _count(self, cls: ClassDef, delta: int) -> None:
        for c in cls.mro():
            c.ninst += delta
            if c.ninst == 0:
                # the last instance is gone: classwide variables start again from their
                # initial values with the next instance
                for cell in c.classwide.values():
                    self.release_cell(cell)
                c.classwide.clear()

    def copy_value(self, v: object, *, oncopy: bool = True) -> object:
        if isinstance(v, Array):
            return Array({k: Cell(self.copy_value(c.value, oncopy=oncopy))
                          for k, c in v.items.items()})
        if not isinstance(v, Instance):
            return v
        new = Instance(v.cls)
        self._count(v.cls, 1)
        new.members = {k: Cell(self.copy_value(c.value, oncopy=oncopy))
                       for k, c in v.members.items()}
        new.dynamic = list(v.dynamic)
        found = v.cls.find_program("oncopy") if oncopy else None
        if found is not None:
            src = self.objkey(Cell(v, 0))
            self.drop_rv(self.call(Cell(new, 0), found[0], found[1], "", oncopy_src=src))
        return new

    def release_cell(self, cell: Cell) -> None:
        cell.nref -= 1
        if cell.nref <= 0 and not cell.dead:
            cell.dead = True
            self.release_value(cell.value)

    def release_value(self, v: object) -> None:
        if isinstance(v, Array):
            for c in list(v.items.values()):
                self.release_cell(c)
            v.items.clear()
        elif isinstance(v, Instance):
            self.destroy(v)

    def destroy(self, inst: Instance) -> None:
        found = inst.cls.find_program("destructor")
        if found is not None:
            # VERIFY: whether an error in a .destructor is shown; here it is ignored
            try:
                self.drop_rv(self.call(Cell(inst, 0), found[0], found[1], ""))
            except StataError:
                pass
        for c in list(inst.members.values()):
            self.release_cell(c)
        inst.members.clear()
        self._count(inst.cls, -1)

    # -- member programs --------------------------------------------------------------------------
    def call(self, cell: Cell, owner: ClassDef, prog: Program, args: str, *,
             oncopy_src: str = "") -> object:
        from .programs import call_program
        frame = Frame(cell, owner, oncopy_src=oncopy_src)
        self.ctx.append(frame)
        try:
            call_program(self.s, prog, args)
        except BaseException:
            self.drop_rv(frame.retval)
            raise
        finally:
            self.ctx.pop()
        return frame.retval

    def drop_rv(self, rv: object) -> None:
        """Discards a value returned by a member program that nobody used."""
        if rv is None:
            return
        if rv.cell is not None:
            self.release_cell(rv.cell)
        else:
            self.release_value(rv.value)

    # -- keys and unames --------------------------------------------------------------------------
    def objkey(self, cell: Cell) -> str:
        key = f"__objkey_{cell.uid}"          # VERIFY: the form of Stata's keys
        self.keys[key] = cell
        return key

    @staticmethod
    def uname(cell: Cell) -> str:
        return f"__uname{cell.uid:06X}"       # VERIFY: the form of Stata's unames

    # -- dropping ---------------------------------------------------------------------------------
    def drop(self, name: str) -> None:
        cell = self.globals.pop(name)
        self.release_cell(cell)

    def drop_names(self, names: list[str]) -> None:
        """Top-level objects named by tempname are dropped when the program ends."""
        for n in names:
            if n in self.globals:
                self.drop(n)

    def clear(self) -> None:
        for name in list(self.globals):
            if name in self.globals:
                self.drop(name)
        self.defs.clear()
        self.keys.clear()


def system(s: "Session") -> ClassSystem:
    cs = getattr(s, "_classes", None)
    if cs is None:
        cs = s._classes = ClassSystem(s)
    return cs


def clear_classes(s: "Session") -> None:
    """discard and clear all drop every instance and class definition."""
    cs = getattr(s, "_classes", None)
    if cs is not None:
        cs.clear()
