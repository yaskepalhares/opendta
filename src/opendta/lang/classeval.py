"""Statements, references and rvalues of the class system ([P] class, sections 7-10).

    lvalue = rvalue                 .a.b[exp].c = ...   (a copy of the right side)
    lvalue.ref = lvalue.ref         .a.ref = .b.ref     (the same object)
    lvalue.ref = rvalue             makes .a unique again
    .a.b.prog [args]                runs a member program in the context of .a.b
    .a.Declare [type] name [= rvalue], .a.list.Arrpush "x", ...    (built-in modifiers)

rvalue is "string", `"string"', #, exp, (exp), .id[.id...] [program arguments], {} or
{el, el, ...}. A reference `.x` inside a member program means .context.x when the context has
such a member (or program), and the global .x otherwise; .Global and .Local force one of them.
`.a.b' in a macro expands to the printable form of the object; an error gives nothing.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

from ..core.errors import ExitRequest, StataError
from ..core.formats import number_to_macro
from .classes import (Array, Cell, ClassDef, ClassSystem, Frame, Instance, mismatch, same_type,
                      system, type_of)

if TYPE_CHECKING:
    from ..session import Session

_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

# built-ins whose result can be followed by more of the reference (they are member variables)
_MEMBER_BUILTINS = ("classmv", "instancemv", "dynamicmv", "superclass")
_BUILTINS = ("new", "copy", "ref", "ref_n", "objtype", "isa", "classname", "isofclass",
             "objkey", "uname", "arrnels", "arrindexof", "oncopy_src",
             "Declare", "Arrdropel", "Arrdropall", "Arrpop", "Arrpush") + _MEMBER_BUILTINS


@dataclass
class Seg:
    name: str
    index: str | None = None    # name[exp]
    arg: str | None = None      # Super(classname)


@dataclass
class RV:
    """A returned value: a fresh value owned by whoever receives it, or a shared Cell (.ref)."""
    value: object = None
    cell: Cell | None = None


@dataclass
class Node:
    value: object
    cell: Cell | None = None
    slot: tuple | None = None   # (dict, key) where the cell is kept, for assignment


def parse_path(text: str, pos: int = 0) -> tuple[list[Seg], int]:
    segs: list[Seg] = []
    while pos < len(text) and text[pos] == ".":
        m = _NAME.match(text, pos + 1)
        if not m:
            break
        seg = Seg(m.group(0))
        pos = m.end()
        if pos < len(text) and text[pos] == "[":
            j = _closing(text, pos)
            seg.index, pos = text[pos + 1:j], j + 1
        elif seg.name == "Super" and pos < len(text) and text[pos] == "(":
            j = _closing(text, pos)
            seg.arg, pos = text[pos + 1:j].strip(), j + 1
        segs.append(seg)
    return segs, pos


def _closing(text: str, i: int) -> int:
    """Position of the bracket that closes text[i], skipping quotes and nested brackets."""
    pairs = {"(": ")", "[": "]", "{": "}"}
    stack = [pairs[text[i]]]
    j = i + 1
    while j < len(text):
        c = text[j]
        if c == '"':
            k = text.find('"', j + 1)
            j = len(text) if k == -1 else k
        elif text.startswith('`"', j):
            j = _compound_end(text, j) - 1
        elif c in pairs:
            stack.append(pairs[c])
        elif c == stack[-1]:
            stack.pop()
            if not stack:
                return j
        j += 1
    raise StataError(198, "invalid syntax")


def _compound_end(text: str, i: int) -> int:
    """text[i:] starts with `" ; returns the position after the matching "'."""
    depth, j = 0, i
    while j < len(text) - 1:
        if text.startswith('`"', j):
            depth += 1
            j += 2
        elif text.startswith("\"'", j):
            depth -= 1
            j += 2
            if depth == 0:
                return j
        else:
            j += 1
    raise StataError(198, "invalid syntax")


def _quoted(t: str) -> tuple[str, int]:
    """A simple or compound quoted string at the start of t -> (contents, end)."""
    if t.startswith('`"'):
        end = _compound_end(t, 0)
        return t[2:end - 2], end
    k = t.find('"', 1)
    if k == -1:
        raise StataError(198, "invalid syntax")
    return t[1:k], k + 1


def _split_top(text: str) -> list[str]:
    """Splits on the commas that are outside quotes and brackets."""
    out, depth, start, j = [], 0, 0, 0
    while j < len(text):
        c = text[j]
        if c == '"':
            k = text.find('"', j + 1)
            j = len(text) if k == -1 else k
        elif text.startswith('`"', j):
            j = _compound_end(text, j) - 1
        elif c in "([{":
            depth += 1
        elif c in ")]}":
            depth -= 1
        elif c == "," and depth == 0:
            out.append(text[start:j])
            start = j + 1
        j += 1
    out.append(text[start:])
    return out


def _is_ref(t: str) -> bool:
    return len(t) > 1 and t[0] == "." and (t[1].isalpha() or t[1] == "_")


def _dotted(segs: list[Seg]) -> str:
    return "".join(f".{g.name}" + (f"[{g.index}]" if g.index is not None else "") for g in segs)


def _kind(v: object) -> str:
    """Type names of Stata's "may not assign" message."""
    if isinstance(v, (Instance, ClassDef)):
        return "class"
    return type_of(v)


def not_found(segs: list[Seg]) -> StataError:
    return StataError(111, f"{_dotted(segs)} not found")           # VERIFY: message and code


class Evaluator:
    def __init__(self, cs: ClassSystem):
        self.cs = cs
        self.s = cs.s

    # -- references -------------------------------------------------------------------------------
    def resolve(self, segs: list[Seg], *, lvalue: bool = False, global_context: bool = False):
        """-> ("var", Node) | ("new", (dict, key)) | ("prog", Node, owner, Program)
              | ("builtin", Node, name) | ("class", ClassDef)."""
        cs = self.cs
        segs = list(segs)
        force = ""
        if segs[0].name in ("Global", "Local") and segs[0].index is None:
            force = segs.pop(0).name
            if not segs:
                raise StataError(198, "invalid syntax")
        frame = None if global_context or force == "Global" or not cs.ctx else cs.ctx[-1]
        first = segs[0]
        if frame is not None and (first.name in _BUILTINS or first.name == "Super"
                                  or self._in_context(frame, first.name)):
            return self._walk(Node(frame.cell.value, frame.cell), segs, frame, lvalue, segs)
        if force == "Local":
            raise not_found(segs)
        name = first.name
        if name in cs.keys:
            cell = cs.keys[name]
            if cell.dead:
                raise not_found(segs)
            node = Node(cell.value, cell)
        elif name in cs.globals:
            node = Node(cs.globals[name].value, cs.globals[name], (cs.globals, name))
        elif lvalue:
            if len(segs) == 1 and first.index is None:
                return ("new", (cs.globals, name))
            raise not_found(segs[:1])
        else:
            try:
                cls = cs.load(name)
            except StataError as e:
                if len(segs) > 1 and e.message == f"class {name} not found":
                    raise StataError(4023, f"{_dotted(segs)[1:]}: class member function not found")
                raise
            if first.index is not None:
                raise StataError(198, "invalid syntax")
            if len(segs) == 1:
                return ("class", cls)
            return self._walk(Node(cls), segs[1:], None, lvalue, segs)
        node = self._index(node, first, lvalue and len(segs) == 1, segs[:1])
        if isinstance(node, tuple):
            return node
        if len(segs) == 1:
            return ("var", node)
        return self._walk(node, segs[1:], None, lvalue, segs)

    @staticmethod
    def _in_context(frame: Frame, name: str) -> bool:
        inst = frame.cell.value
        return isinstance(inst, Instance) and (inst.member(name) is not None
                                               or inst.cls.find_program(name) is not None)

    def _walk(self, node: Node, segs: list[Seg], frame: Frame | None, lvalue: bool,
              full: list[Seg]):
        i = 0
        while i < len(segs):
            seg = segs[i]
            last = i == len(segs) - 1
            v = node.value
            if seg.name == "Super":
                if last or not isinstance(v, Instance):
                    raise StataError(198, "invalid syntax")
                owner = frame.owner if frame is not None and node.cell is frame.cell else v.cls
                parents = owner.parents
                if seg.arg:
                    parents = [p for c in owner.mro() for p in c.parents if p.name == seg.arg][:1]
                    if not parents:
                        raise StataError(111, f"{seg.arg} is not a parent of {owner.name}")
                pname = segs[i + 1].name
                if i + 1 != len(segs) - 1:
                    raise StataError(198, "invalid syntax")
                for p in parents:
                    found = p.find_program(pname)
                    if found is not None:
                        return ("prog", node, found[0], found[1])
                if pname == "new":
                    return ("noop",)
                raise not_found(full)
            if seg.name in _MEMBER_BUILTINS:
                node = self._index(Node(self._member_array(node, seg.name)), seg, False, full)
                i += 1
                continue
            if seg.name in _BUILTINS:
                if not last:
                    raise StataError(198, "invalid syntax")      # VERIFY
                return ("builtin", node, seg.name)
            if isinstance(v, Instance):
                hit = v.member(seg.name)
                if hit is not None:
                    node = Node(hit[1].value, hit[1], (hit[0], seg.name))
                    node = self._index(node, seg, lvalue and last, full[:len(full) - len(segs) + i + 1])
                    if isinstance(node, tuple):
                        return node
                    i += 1
                    continue
                found = v.cls.find_program(seg.name)
                if found is not None:
                    if not last or seg.index is not None:
                        raise StataError(198, "invalid syntax")
                    return ("prog", node, found[0], found[1])
            if lvalue:
                raise StataError(198, f"{_dotted(full)[1:]} invalid name")
            raise not_found(full[:len(full) - len(segs) + i + 1])
        return ("var", node)

    def _index(self, node: Node, seg: Seg, creatable: bool, segs: list[Seg]):
        if seg.index is None:
            return node
        if not isinstance(node.value, Array):
            raise StataError(198, f"{_dotted(segs)}: not an array")          # VERIFY
        x = self.s.eval(seg.index)
        if isinstance(x, str) or x != x or x >= 8.988465674311579e307:
            raise StataError(198, "invalid syntax")
        k = math.trunc(x)
        items = node.value.items
        if k in items:
            return Node(items[k].value, items[k], (items, k))
        if creatable and k >= 1:
            return ("new", (items, k))
        raise not_found(segs)

    def _member_array(self, node: Node, which: str) -> Array:
        """.classmv, .instancemv, .dynamicmv and .superclass: arrays whose elements are the
        members themselves; with no elements they are an error, as the manual says."""
        v = node.value
        cells: list[Cell] = []
        if isinstance(v, Instance):
            eff = v.cls.effective()
            if which == "classmv":
                cells = [v.member(n)[1] for n, d in eff.items()
                         if d.scope == "classwide" and v.member(n)]
            elif which == "instancemv":
                cells = [v.members[n] for n, d in eff.items()
                         if d.scope != "classwide" and n in v.members]
            elif which == "dynamicmv":
                cells = [v.members[n] for n in v.dynamic]
            else:
                cells = [Cell(p, 0) for p in v.cls.parents]
        elif isinstance(v, ClassDef) and which == "superclass":
            cells = [Cell(p, 0) for p in v.parents]
        if not cells:
            raise StataError(111, f".{which} not found")                     # VERIFY
        return Array({k: c for k, c in enumerate(cells, start=1)})

    # -- rvalues ----------------------------------------------------------------------------------
    def rvalue(self, text: str, *, global_context: bool = False) -> RV:
        t = text.strip()
        if not t:
            raise StataError(198, "invalid syntax")
        if t[0] == '"' or t.startswith('`"'):
            content, end = _quoted(t)
            if t[end:].strip():
                raise StataError(198, "invalid syntax")
            return RV(content)
        if t[0] == "{":
            j = _closing(t, 0)
            if t[j + 1:].strip():
                raise StataError(198, "invalid syntax")
            return RV(self._array(t[1:j], global_context))
        if _is_ref(t):
            segs, end = parse_path(t)
            return self.object_rvalue(segs, t[end:], global_context=global_context)
        return RV(self._expression(t))

    def _expression(self, t: str) -> object:
        v = self.s.eval(t)
        if isinstance(v, str):
            return v
        try:
            return float(v)
        except (TypeError, ValueError):
            raise mismatch()                    # a matrix

    def _array(self, inner: str, global_context: bool) -> Array:
        arr = Array()
        if not inner.strip():
            return arr
        for k, el in enumerate(_split_top(inner), start=1):
            el = el.strip()
            if not el:
                continue
            if el[0] == "{":
                raise StataError(198, "recursive array definitions not allowed")   # VERIFY
            if el[0] == "[" and el[-1] == "]":
                rv = self.rvalue(el[1:-1], global_context=global_context)
            else:
                rv = self.rvalue(el, global_context=global_context)
            arr.items[k] = Cell(self.take(rv))
        return arr

    def object_rvalue(self, segs: list[Seg], args: str, *, global_context: bool = False) -> RV:
        r = self.resolve(segs, global_context=global_context)
        kind = r[0]
        if kind == "prog":
            rv = self.cs.call(r[1].cell or Cell(r[1].value, 0), r[2], r[3], args.strip())
            return rv if rv is not None else RV("")
        if kind == "builtin":
            rv = self.builtin(r[1], r[2], args)
            return rv if rv is not None else RV("")
        if kind == "noop":
            return RV("")
        if args.strip():
            raise StataError(198, "invalid syntax")
        if kind == "class":
            raise StataError(198, f".{r[1].name} is a class, not an object")      # VERIFY
        return RV(self.cs.copy_value(r[1].value))

    def take(self, rv: RV) -> object:
        """The value of rv as a fresh value."""
        if rv.cell is None:
            return rv.value
        v = self.cs.copy_value(rv.cell.value)
        self.cs.release_cell(rv.cell)
        return v

    def value_of(self, text: str, *, global_context: bool = False) -> object:
        return self.take(self.rvalue(text, global_context=global_context))

    # -- statements -------------------------------------------------------------------------------
    def statement(self, text: str) -> None:
        segs, end = parse_path(text)
        if not segs:
            raise StataError(198, "invalid syntax")
        rest = text[end:]
        rs = rest.lstrip()
        if rs.startswith("=") and not rs.startswith("=="):
            ref = segs[-1].name == "ref" and segs[-1].index is None and len(segs) > 1
            self.assign(segs[:-1] if ref else segs, ref, rs[1:])
            return
        r = self.resolve(segs)
        if r[0] == "prog":
            self.cs.drop_rv(self.cs.call(r[1].cell or Cell(r[1].value, 0), r[2], r[3],
                                         rest.strip()))
        elif r[0] == "builtin":
            self.cs.drop_rv(self.builtin(r[1], r[2], rest))
        elif r[0] in ("var", "class") and rest.strip():
            raise StataError(198, "invalid syntax")
        # VERIFY: a bare reference to a variable or a class does nothing here

    def assign(self, lhs: list[Seg], ref: bool, rhs: str) -> None:
        cs = self.cs
        rv = self.rvalue(rhs)
        try:
            r = self.resolve(lhs, lvalue=True)
        except BaseException:
            cs.drop_rv(rv)
            raise
        if r[0] == "new":
            holder, key = r[1]
            if ref and rv.cell is not None:
                holder[key] = rv.cell                  # the reference held by rv moves here
            else:
                holder[key] = Cell(self.take(rv))
            return
        if r[0] != "var" or r[1].slot is None:
            cs.drop_rv(rv)
            raise StataError(198, "invalid syntax")    # VERIFY: assignment to a program/built-in
        node = r[1]
        if not ref and rv.cell is not None:
            cs.drop_rv(rv)
            raise StataError(9999, "= assignment of a reference is not allowed")
        new_value = rv.cell.value if rv.cell is not None else rv.value
        if not same_type(node.value, new_value):
            cs.drop_rv(rv)
            if isinstance(node.value, Instance) and not isinstance(new_value, Instance):
                raise StataError(4026, "rvalue invalid")
            # VERIFY: an instance of another class, and arrays
            raise StataError(4027, "incompatible types\n     may not assign "
                             f"{_kind(new_value)} to {_kind(node.value)}")
        holder, key = node.slot
        if ref:
            if rv.cell is node.cell:
                cs.release_cell(rv.cell)
                return
            holder[key] = rv.cell if rv.cell is not None else Cell(rv.value)
            cs.release_cell(node.cell)
            return
        new_value = self.take(rv)
        old = node.cell.value
        node.cell.value = new_value         # in place: identifiers sharing the cell see it too
        cs.release_value(old)

    # -- built-ins --------------------------------------------------------------------------------
    def builtin(self, node: Node, name: str, args: str) -> RV | None:
        cs = self.cs
        v = node.value
        a = args.strip()
        if name == "new":
            cls = v if isinstance(v, ClassDef) else v.cls if isinstance(v, Instance) else None
            if cls is None:
                raise StataError(198, ".new may be used only with classes and instances")  # VERIFY
            return RV(cs.new_instance(cls, a))
        if name in ("Declare", "Arrdropall", "Arrpop", "Arrpush", "Arrdropel"):
            self._modifier(node, name, a)
            return None
        if a and name not in ("isofclass", "arrindexof"):
            raise StataError(198, "invalid syntax")
        if name == "copy":
            return RV(cs.copy_value(v))
        if name == "ref":
            if node.cell is None:
                raise StataError(198, "invalid syntax")
            node.cell.nref += 1
            return RV(cell=node.cell)
        if name == "ref_n":
            return RV(float(node.cell.nref if node.cell is not None else 1))
        if name == "objtype":
            return RV(type_of(v))
        if name == "isa":
            if isinstance(v, Instance):
                return RV("class")
            return RV("classtype" if isinstance(v, ClassDef) else type_of(v))
        if name == "classname":
            if isinstance(v, (Instance, ClassDef)):
                return RV(v.cls.name if isinstance(v, Instance) else v.name)
            return RV("")
        if name == "isofclass":
            cls = v.cls if isinstance(v, Instance) else v if isinstance(v, ClassDef) else None
            return RV(1.0 if cls is not None and cls.isof(a.lstrip(".")) else 0.0)
        if name == "objkey":
            return RV(cs.objkey(node.cell or Cell(v, 0)))
        if name == "uname":
            return RV(cs.uname(node.cell or Cell(v, 0)))
        if name == "oncopy_src":
            src = next((f.oncopy_src for f in reversed(cs.ctx) if f.oncopy_src), "")
            return RV(src)
        if not isinstance(v, Array):
            raise StataError(4020, "not an array")      # VERIFY: seen for .Arrpush only
        if name == "arrnels":
            return RV(float(v.nels()))
        if name == "arrindexof":
            from .words import strip_outer_quotes
            target = strip_outer_quotes(a)
            for k in sorted(v.items):
                # VERIFY: elements that are not strings never match
                if isinstance(v.items[k].value, str) and v.items[k].value == target:
                    return RV(float(k))
            return RV(0.0)
        raise StataError(198, "invalid syntax")

    def _modifier(self, node: Node, name: str, a: str) -> None:
        cs = self.cs
        v = node.value
        if name == "Declare":
            if not isinstance(v, Instance):
                raise StataError(198, ".Declare may be used only with class instances")  # VERIFY
            d = cs.declaration(a, "instancespecific")
            if v.member(d.name) is not None:
                cs.release_value(d.value)
                raise StataError(4017, f"{d.name}: name already declared in current class")
            v.members[d.name] = Cell(d.value)
            v.dynamic.append(d.name)
            return
        if name == "Arrdropel" and not a:
            # .array[i].Arrdropel: the element itself
            holder = node.slot[0] if node.slot else None
            if not isinstance(holder, dict) or not isinstance(node.slot[1], int):
                raise StataError(198, ".Arrdropel may be used only with array elements")  # VERIFY
            cs.release_cell(holder.pop(node.slot[1]))
            return
        if not isinstance(v, Array):
            raise StataError(4020, "not an array")
        if name == "Arrdropel":
            k = math.trunc(float(self.s.eval(a)))
            if k in v.items:
                cs.release_cell(v.items.pop(k))
        elif name == "Arrdropall":
            cs.release_value(v)
        elif name == "Arrpop":
            n = v.nels()
            if n:
                if n in v.items:
                    cs.release_cell(v.items.pop(n))
                v.top = n - 1
        elif name == "Arrpush":
            v.items[v.nels() + 1] = Cell(self.value_of(a))

    # -- macro expansion --------------------------------------------------------------------------
    def text(self, inner: str) -> str:
        """`.id...' : the printable form; an error gives nothing (the message is not shown)."""
        segs, end = parse_path(inner)
        if not segs:
            return ""
        args = inner[end:]
        try:
            r = self.resolve(segs)
            if r[0] == "var":
                return "" if args.strip() else self.printable(r[1].value)
            if r[0] == "prog":
                rv = self.cs.call(r[1].cell or Cell(r[1].value, 0), r[2], r[3], args.strip())
            elif r[0] == "builtin":
                rv = self.builtin(r[1], r[2], args)
            else:
                return ""
        except (StataError, ExitRequest):
            return ""
        if rv is None:
            return ""
        out = self.printable(rv.cell.value if rv.cell is not None else rv.value)
        self.cs.drop_rv(rv)
        return out

    def printable(self, v: object) -> str:
        if isinstance(v, str):
            return v
        if isinstance(v, float):
            return number_to_macro(v)
        if isinstance(v, Instance):
            found = v.cls.find_program("macroexpand")
            if found is None:
                return ""
            # run on a copy, so that .macroexpand cannot change the object
            copy = self.cs.copy_value(v, oncopy=False)
            try:
                rv = self.cs.call(Cell(copy, 0), found[0], found[1], "")
            except (StataError, ExitRequest):
                return ""
            finally:
                self.cs.destroy(copy)
            if rv is None:
                return ""
            res = rv.cell.value if rv.cell is not None else rv.value
            out = self.printable(res) if isinstance(res, (str, float)) else ""
            self.cs.drop_rv(rv)
            return out
        return ""


def evaluator(s: "Session") -> Evaluator:
    return Evaluator(system(s))


def statement(s: "Session", text: str) -> None:
    evaluator(s).statement(text)


def is_statement(text: str) -> bool:
    return _is_ref(text)


def macro_text(s: "Session", inner: str) -> str:
    return evaluator(s).text(inner)
