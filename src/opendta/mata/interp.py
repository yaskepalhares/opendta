"""Execução do Mata: instruções, expressões, funções do usuário e a
ligação com a sessão (blocos `mata ... end`, `mata: instrução`).

Passagem de argumentos: como no Mata, por referência. Quando o argumento é
uma variável, o que a função fizer com o parâmetro volta para a variável
de quem chamou.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from . import ops
from .display import render
from .lexer import MataSyntaxError
from .parser import (AddrOf, Assign, Binary, Block, Break, Call, Continue, Decl, Deref, DoWhile,
                     Empty, ExprStmt, For, FuncDef, If, IncDec, Member, Name, Num, Parser, RangeSub,
                     Return, Str, StructDef, Subscript, Ternary, Transpose, Unary, While)
from .values import MV, MataError, Pointer, StructVal, SYS, empty, real, string, truth

if TYPE_CHECKING:
    from ..session import Session


class _Break(Exception):
    pass


class _Continue(Exception):
    pass


class _Return(Exception):
    def __init__(self, value):
        self.value = value


class MataLeave(Exception):
    """exit() sem código: sai do bloco mata em silêncio."""


class Frame:
    def __init__(self, name: str = "<istmt>", parent: "Frame | None" = None):
        self.name = name
        self.vars: dict[str, MV] = {}
        self.types: dict[str, tuple[str, str]] = {}
        self.parent = parent           # só para variáveis external (globais)

    def get(self, name: str) -> MV:
        if name in self.vars:
            return self.vars[name]
        raise MataError(3499, f"{name} not found")

    def has(self, name: str) -> bool:
        return name in self.vars


class Ref:
    """Argumento passado por referência a uma função da biblioteca."""

    def __init__(self, frame: Frame, name: str):
        self.frame = frame
        self.name = name

    def set(self, value: MV) -> None:
        self.frame.vars[self.name] = value

    def get(self) -> MV | None:
        return self.frame.vars.get(self.name)


_DEFAULT_ORG_TYPES = {"real": "real", "complex": "complex", "string": "string", "pointer": "pointer",
                      "numeric": "real", "transmorphic": "real"}


def _names_used(node, out: set | None = None) -> set:
    """Nomes citados em qualquer ponto de um trecho da AST (declarações não contam)."""
    import dataclasses
    if out is None:
        out = set()
    if isinstance(node, Name):
        out.add(node.name)
    elif isinstance(node, Call):
        out.add(node.func)
        _names_used(node.args, out)
    elif isinstance(node, Decl):
        pass
    elif isinstance(node, (list, tuple)):
        for x in node:
            _names_used(x, out)
    elif dataclasses.is_dataclass(node):
        for fld in dataclasses.fields(node):
            _names_used(getattr(node, fld.name), out)
    return out


class MataEngine:
    def __init__(self, session: "Session"):
        self.s = session
        self.globals = Frame()
        self.funcs: dict[str, FuncDef] = {}
        self.structs: dict[str, StructDef] = {}
        self.depth = 0
        self.argc: list[int] = []      # nº de argumentos de cada chamada (args())
        from . import library  # noqa: F401  (registra as funções)

    # ------------------------------------------------------------------
    # Entrada a partir do Stata
    # ------------------------------------------------------------------
    def run_statement(self, text: str) -> None:
        """`mata: instrução` ou `mata instrução` numa linha só."""
        try:
            stmts = Parser(text).parse_program()
        except MataSyntaxError as e:
            raise StataError(3000, str(e))
        for st in stmts:
            self._top(st)

    def run_block(self, chunks: list[tuple[str, tuple[str, ...]]], *, echo: bool, stop_on_error: bool) -> None:
        """Bloco mata ... end. `chunks`: (texto sem comentários, linhas para o
        eco). Cada instrução é ecoada com ': ' (continuações com '> ')."""
        out = self.s.output
        k = 0
        n = len(chunks)
        while k < n:
            text, raw = chunks[k]
            if not text.strip():
                k += 1
                continue
            if text.strip().split()[0] == "mata" and len(text.split()) > 1:
                if echo:
                    self._echo(raw, first=True)
                self._guard(lambda t=text: self.subcommand(t.split(None, 1)[1]), stop_on_error, echo)
                k += 1
                continue
            dup = self._redefinition(text)
            if dup:
                # o Stata acusa logo depois do cabeçalho; o corpo vira
                # instruções soltas (observado no Stata 14)
                if echo:
                    self._echo(raw, first=True)
                self._guard(lambda d=dup: (_ for _ in ()).throw(StataError(3000, f"{d}() already exists")),
                            stop_on_error, echo, self._remaining(chunks, k))
                k += 1
                continue
            buf, echo_lines = text, list(raw)
            j = k
            while True:
                try:
                    stmts = Parser(buf).parse_program()
                    break
                except MataSyntaxError as e:
                    if e.incomplete and j + 1 < n:
                        j += 1
                        buf += "\n" + chunks[j][0]
                        echo_lines.extend(chunks[j][1])
                        continue
                    stmts = e
                    break
            if echo:
                for idx, ln in enumerate(echo_lines):
                    out.write((": " if idx == 0 else "> ") + ln + "\n", "command")
            if isinstance(stmts, MataSyntaxError):
                self._guard(lambda e=stmts: (_ for _ in ()).throw(StataError(3000, str(e))),
                            stop_on_error, echo, self._remaining(chunks, j))
            else:
                def run(sts=stmts):
                    for st in sts:
                        self._top(st)
                self._guard(run, stop_on_error, echo, self._remaining(chunks, j))
            k = j + 1

    @staticmethod
    def _remaining(chunks, j: int) -> int:
        return sum(1 for text, _ in chunks[j + 1:] if text.strip())

    _HEADER = None

    def _redefinition(self, text: str) -> str | None:
        """Nome da função (ou struct) que a linha redefine, se já existe."""
        import re
        if MataEngine._HEADER is None:
            types = r"(?:real|string|complex|pointer|transmorphic|numeric|void|struct\s+\w+|class\s+\w+)"
            orgs = r"(?:scalar|vector|rowvector|colvector|matrix)"
            MataEngine._HEADER = re.compile(
                rf"^\s*(?:function\s+(\w+)\s*\(|{types}\s+(?:{orgs}\s+)?(\w+)\s*\(|"
                rf"{orgs}\s+(\w+)\s*\(|struct\s+(\w+)\s*\{{?\s*$)")
        m = MataEngine._HEADER.match(text)
        if not m:
            return None
        name = next(g for g in m.groups() if g)
        return name if (name in self.funcs or name in self.structs) else None

    def _echo(self, raw, first: bool) -> None:
        for idx, ln in enumerate(raw):
            self.s.output.write((": " if idx == 0 and first else "> ") + ln + "\n", "command")

    def _guard(self, action, stop_on_error: bool, echo: bool, remaining: int = 0) -> None:
        out = self.s.output
        try:
            action()
        except StataError as e:
            if stop_on_error:
                # mata: para no erro e conta as linhas que não rodaram (Stata 14)
                if remaining:
                    e.message = (e.message + "\n" if e.message else "") + \
                        f"({remaining} line{'s' if remaining != 1 else ''} skipped)"
                raise
            self.s.report_error(e)
            self.s.set_rc(e.rc)
        if echo:
            out.ensure_line_start()
            out.write("\n", "text")

    # ------------------------------------------------------------------
    # Instruções de nível superior
    # ------------------------------------------------------------------
    def _top(self, st) -> None:
        try:
            if isinstance(st, FuncDef):
                if st.name in self.funcs:
                    raise MataError(3000, f"{st.name}() already exists")   # VERIFICAR
                self.funcs[st.name] = st
                # o compilador avisa os argumentos que o corpo não usa
                # (observado no Stata 14, compat 0406)
                used = _names_used(st.body)
                for prm in st.params:
                    if prm.name not in used:
                        self.s.output.write(f"note: argument {prm.name} unused\n", "text")
                return
            if isinstance(st, StructDef):
                self.structs[st.name] = st
                return
            if isinstance(st, ExprStmt) and not isinstance(st.expr, (Assign, IncDec)):
                v = self.eval(st.expr, self.globals)
                if isinstance(v, MV):
                    for line in render(v):
                        self.s.output.write(line + "\n", "result")
                return
            self.exec(st, self.globals)
        except (_Break, _Continue):
            raise StataError(3000, "break or continue outside of a loop")
        except _Return:
            return
        except MataError as e:
            raise StataError(e.code, self._format_error(e))
        except Exception as e:     # noqa: BLE001
            from .library import MataExit
            if isinstance(e, MataExit):
                if e.rc:
                    raise StataError(e.rc, "")
                raise MataLeave()
            raise

    @staticmethod
    def _format_error(e: MataError) -> str:
        if e.code == 3000 and not e.where:
            return e.message
        lines = []
        stack = e.where or []
        if stack:
            lines.append(f"{stack[0]:>24}:  {e.code:>4}  {e.message}")
            for name in stack[1:]:
                lines.append(f"{name:>24}:     -  function returned error")
            lines.append(f"{'<istmt>':>24}:     -  function returned error")
        else:
            lines.append(f"{'<istmt>':>24}:  {e.code:>4}  {e.message}")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Subcomandos (mata clear, mata describe, mata drop)
    # ------------------------------------------------------------------
    def subcommand(self, text: str) -> None:
        words = text.split()
        if not words:
            return
        sub = words[0]
        out = self.s.output
        if sub == "clear":
            self.globals = Frame()
            self.funcs.clear()
            self.structs.clear()
            return
        if sub in ("describe", "d", "de", "des", "desc", "descr", "descri", "describ"):
            # VERIFICAR: layout do mata describe
            out.write("\n      # bytes   type                        name and extent\n", "text")
            out.write("-" * 79 + "\n", "text")
            for name in sorted(self.funcs):
                f = self.funcs[name]
                out.write(f"{'':>13}   {(f.rettype + ' ' + f.retorg).strip():<27} {name}()\n", "text")
            for name in sorted(self.globals.vars):
                v = self.globals.vars[name]
                out.write(f"{v.a.nbytes:>13,}   {(v.t + ' ' + v.orgtype()):<27} {name}"
                          f"[{v.rows},{v.cols}]\n" if not v.is_scalar else
                          f"{v.a.nbytes:>13,}   {(v.t + ' ' + v.orgtype()):<27} {name}\n", "text")
            out.write("-" * 79 + "\n", "text")
            return
        if sub == "drop":
            for w in words[1:]:
                if w.endswith("()"):
                    if self.funcs.pop(w[:-2], None) is None:
                        raise StataError(3499, f"{w} not found")
                else:
                    if self.globals.vars.pop(w, None) is None:
                        raise StataError(3499, f"{w} not found")
            return
        if sub in ("set", "mlib", "mosave", "which", "memory", "query", "rename"):
            if sub == "rename" and len(words) == 3:
                old, new = words[1], words[2]
                if old.endswith("()"):
                    self.funcs[new[:-2]] = self.funcs.pop(old[:-2])
                else:
                    self.globals.vars[new] = self.globals.vars.pop(old)
            return       # aceitos sem efeito (VERIFICAR)
        self.run_statement(text)

    # ------------------------------------------------------------------
    # Instruções
    # ------------------------------------------------------------------
    def exec(self, st, fr: Frame) -> None:
        if isinstance(st, ExprStmt):
            v = self.eval(st.expr, fr)
            # no nível de cima (fora de funções), expressões mostram o valor,
            # mesmo dentro de laços: for (i=1; i<=3; i++) i
            if fr is self.globals and isinstance(v, MV) and not isinstance(st.expr, (Assign, IncDec)):
                for line in render(v):
                    self.s.output.write(line + "\n", "result")
        elif isinstance(st, Block):
            for s2 in st.body:
                self.exec(s2, fr)
        elif isinstance(st, If):
            if truth(self.eval(st.cond, fr)):
                self.exec(st.then, fr)
            elif st.other is not None:
                self.exec(st.other, fr)
        elif isinstance(st, For):
            if st.init is not None:
                self.eval(st.init, fr)
            while st.cond is None or truth(self.eval(st.cond, fr)):
                try:
                    self.exec(st.body, fr)
                except _Break:
                    break
                except _Continue:
                    pass
                if st.step is not None:
                    self.eval(st.step, fr)
        elif isinstance(st, While):
            while truth(self.eval(st.cond, fr)):
                try:
                    self.exec(st.body, fr)
                except _Break:
                    break
                except _Continue:
                    continue
        elif isinstance(st, DoWhile):
            while True:
                try:
                    self.exec(st.body, fr)
                except _Break:
                    break
                except _Continue:
                    pass
                if not truth(self.eval(st.cond, fr)):
                    break
        elif isinstance(st, Break):
            raise _Break()
        elif isinstance(st, Continue):
            raise _Continue()
        elif isinstance(st, Return):
            raise _Return(None if st.value is None else self.eval(st.value, fr))
        elif isinstance(st, Decl):
            for nm in st.names:
                if st.external:
                    if nm not in self.globals.vars:
                        self.globals.vars[nm] = self._initial(st.eltype, st.org)
                    fr.vars[nm] = self.globals.vars[nm]
                    fr.types[nm] = ("external", st.org)
                    continue
                fr.types[nm] = (st.eltype, st.org)
                fr.vars[nm] = self._initial(st.eltype, st.org)
        elif isinstance(st, (Empty,)):
            pass
        elif isinstance(st, (FuncDef, StructDef)):
            raise MataError(3000, "function definitions must be at the top level")
        else:
            raise MataError(3000, "invalid statement")

    def _initial(self, eltype: str, org: str) -> MV:
        """Valor inicial de uma variável declarada (escalar missing, vetor vazio...)."""
        if eltype.startswith("struct "):
            return MV(np.array([[self.new_struct(eltype[7:])]], dtype=object), "struct")
        t = _DEFAULT_ORG_TYPES.get(eltype, "real")
        if org == "scalar":
            if t == "string":
                return string("")
            if t == "pointer":
                return MV(np.array([[None]], dtype=object), "pointer")
            if t == "complex":
                return MV(np.array([[complex(SYS, 0)]]), "complex")
            return real(SYS)
        r, c = (1, 0) if org == "rowvector" else (0, 1) if org == "colvector" else (0, 0)
        return empty(t, r, c)

    def new_struct(self, sname: str) -> StructVal:
        sd = self.structs.get(sname)
        if sd is None:
            raise MataError(3000, f"struct {sname} not found")
        return StructVal(sname, {nm: self._initial(el, org) for el, org, nm in sd.fields})

    # ------------------------------------------------------------------
    # Expressões
    # ------------------------------------------------------------------
    def eval(self, node, fr: Frame):
        if isinstance(node, Num):
            return real(node.value)
        if isinstance(node, Str):
            return string(node.value)
        if isinstance(node, Name):
            if node.name in fr.vars:
                return fr.vars[node.name]
            if fr is not self.globals and node.name in self.funcs:
                pass
            raise MataError(3499, f"{node.name} not found")
        if isinstance(node, Binary):
            if node.op in ("&&", "||"):
                left = truth(self.eval(node.left, fr))
                if node.op == "&&" and not left:
                    return real(0.0)
                if node.op == "||" and left:
                    return real(1.0)
                return real(1.0 if truth(self.eval(node.right, fr)) else 0.0)
            a = self._value(node.left, fr)
            b = self._value(node.right, fr)
            try:
                return ops.binary(node.op, a, b)
            except MataError as e:
                # o produto de matrizes aparece como "função" (*:  3200 ...); os
                # demais operadores dão o erro em <istmt> (observado no Stata 14)
                if not e.where and node.op == "*":
                    e.where.append(node.op)
                raise
        if isinstance(node, Unary):
            return ops.unary(node.op, self._value(node.operand, fr))
        if isinstance(node, Transpose):
            return ops.transpose(self._value(node.operand, fr))
        if isinstance(node, Ternary):
            return self.eval(node.a if truth(self._value(node.cond, fr)) else node.b, fr)
        if isinstance(node, Assign):
            v = self._value(node.value, fr)
            self.assign(node.target, v, fr)
            return v
        if isinstance(node, IncDec):
            cur = self._value(node.target, fr)
            new = ops.binary("+" if node.op == "++" else "-", cur, real(1.0))
            self.assign(node.target, new, fr)
            return new if node.prefix else cur
        if isinstance(node, Call):
            return self.call(node, fr)
        if isinstance(node, Subscript):
            return self.subscript(self._value(node.base, fr), node, fr)
        if isinstance(node, RangeSub):
            base = self._value(node.base, fr)
            r, c = self._range_index(base, self._value(node.spec, fr))
            return MV(base.a[np.ix_(r, c)].copy(), base.t)
        if isinstance(node, AddrOf):
            return MV(np.array([[self._pointer_to(node.operand, fr)]], dtype=object), "pointer")
        if isinstance(node, Deref):
            p = self._value(node.operand, fr)
            return self._deref(p)
        if isinstance(node, Member):
            return self._member(node, fr)
        raise MataError(3000, "invalid expression")

    def _value(self, node, fr: Frame) -> MV:
        v = self.eval(node, fr)
        if v is None:
            raise MataError(3000, "void function used in expression")   # VERIFICAR
        return v

    # ponteiros e structs --------------------------------------------------
    def _pointer_to(self, target, fr: Frame) -> Pointer:
        if isinstance(target, Name):
            name = target.name
            if name not in fr.vars:
                raise MataError(3499, f"{name} not found")
            return Pointer(lambda: fr.vars[name], lambda v: fr.vars.__setitem__(name, v), name)
        if isinstance(target, Call) and not target.args and target.func in self.funcs:
            return Pointer(lambda: target.func, None, target.func + "()")
        val = self._value(target, fr)
        box = {"v": val.copy()}
        return Pointer(lambda: box["v"], lambda v: box.__setitem__("v", v), "temp")

    def _deref(self, p: MV) -> MV:
        if p.t != "pointer":
            raise MataError(3000, "pointer required")   # VERIFICAR
        ptr = p.scalar()
        if ptr is None:
            raise MataError(3010, "attempt to dereference NULL pointer")
        v = ptr.get()
        if isinstance(v, str):          # ponteiro para função
            raise MataError(3000, "function pointer dereferenced without call")
        return v

    def _member(self, node: Member, fr: Frame) -> MV:
        if fr is self.globals and isinstance(node.base, Name) and not node.arrow:
            v = fr.vars.get(node.base.name)
            if v is None or v.t != "struct":
                # erro de compilação no Stata 14
                raise MataError(3000, "type mismatch:  exp.exp:  transmorphic found where struct expected")
        base = self._value(node.base, fr)
        if node.arrow:
            base = self._deref(base)
        if base.t != "struct":
            raise MataError(3000, "struct required")
        sv = base.scalar()
        if node.name not in sv.fields:
            raise MataError(3000, f"{node.name} not found in struct {sv.sname}")
        return sv.fields[node.name]

    # subscritos -----------------------------------------------------------
    def _index(self, spec, n: int, fr: Frame) -> np.ndarray:
        if spec is None:
            return np.arange(n)
        v = self._value(spec, fr)
        if v.t != "real":
            raise MataError(3253, "nonreal found where real required")
        idx = v.a.ravel()
        if idx.size and (np.any(idx >= SYS) or np.any(idx < 1) or np.any(idx > n)):
            raise MataError(3301, "subscript invalid")
        return idx.astype(np.int64) - 1

    def subscript(self, base: MV, node: Subscript, fr: Frame) -> MV:
        if node.single:
            if base.rows != 1 and base.cols != 1:
                raise MataError(3301, "subscript invalid")
            flat = base.a.ravel()
            idx = self._index(node.rows, flat.size, fr)
            out = flat[idx]
            shape = (1, len(idx)) if base.rows == 1 else (len(idx), 1)
            return MV(out.reshape(shape).copy(), base.t)
        r = self._index(node.rows, base.rows, fr)
        c = self._index(node.cols, base.cols, fr)
        return MV(base.a[np.ix_(r, c)].copy(), base.t)

    def _range_index(self, base: MV, spec: MV) -> tuple[np.ndarray, np.ndarray]:
        s = spec.a
        if spec.t != "real":
            raise MataError(3253, "nonreal found where real required")
        if s.shape == (2, 2):
            r1, c1, r2, c2 = s[0, 0], s[0, 1], s[1, 0], s[1, 1]
        elif s.shape in ((2, 1), (1, 2)) and (base.rows == 1 or base.cols == 1):
            a, b = s.ravel()
            if base.rows == 1:
                r1, r2, c1, c2 = 1, 1, a, b
            else:
                r1, r2, c1, c2 = a, b, 1, 1
        else:
            raise MataError(3301, "subscript invalid")
        r1 = 1 if r1 >= SYS else r1
        c1 = 1 if c1 >= SYS else c1
        r2 = base.rows if r2 >= SYS else r2
        c2 = base.cols if c2 >= SYS else c2
        if not (1 <= r1 <= r2 + 1 and r2 <= base.rows and 1 <= c1 <= c2 + 1 and c2 <= base.cols):
            raise MataError(3301, "subscript invalid")
        return np.arange(int(r1) - 1, int(r2)), np.arange(int(c1) - 1, int(c2))

    # atribuição -------------------------------------------------------------
    def assign(self, target, value: MV, fr: Frame) -> None:
        if isinstance(target, Name):
            declared = fr.types.get(target.name)
            if declared and declared[0] != "external":
                self._check_type(target.name, declared, value)
            new = value.copy()
            if declared and declared[0] == "external":
                self.globals.vars[target.name] = new
            fr.vars[target.name] = new
            return
        if isinstance(target, Deref):
            p = self._value(target.operand, fr)
            if p.t != "pointer" or p.scalar() is None:
                raise MataError(3010, "attempt to dereference NULL pointer")
            ptr = p.scalar()
            if ptr.set is None:
                raise MataError(3000, "invalid lval")
            ptr.set(value.copy())
            return
        if isinstance(target, Member):
            if fr is self.globals and isinstance(target.base, Name) and not target.arrow:
                v = fr.vars.get(target.base.name)
                if v is None or v.t != "struct":
                    raise MataError(3000, "type mismatch:  exp.exp:  transmorphic found where struct expected")
            base = self._value(target.base, fr)
            if target.arrow:
                base = self._deref(base)
            if base.t != "struct":
                raise MataError(3000, "struct required")
            base.scalar().fields[target.name] = value.copy()
            return
        if isinstance(target, (Subscript, RangeSub)):
            base = self._value(target.base, fr)
            if value.t != base.t and not (base.a.size == 0):
                if {value.t, base.t} == {"real", "complex"}:
                    base.a = base.a.astype(np.complex128)
                    base.t = "complex"
                else:
                    raise MataError(3250, "type mismatch")
            if isinstance(target, RangeSub):
                r, c = self._range_index(base, self._value(target.spec, fr))
            elif target.single:
                if base.rows != 1 and base.cols != 1:
                    raise MataError(3301, "subscript invalid")
                n = base.a.size
                idx = self._index(target.rows, n, fr)
                if base.rows == 1:
                    r, c = np.array([0]), idx
                else:
                    r, c = idx, np.array([0])
            else:
                r = self._index(target.rows, base.rows, fr)
                c = self._index(target.cols, base.cols, fr)
            block = (len(r), len(c))
            val = value.a
            if value.is_scalar:
                val = np.broadcast_to(val, block)
            elif target.__class__ is Subscript and target.single and val.size == block[0] * block[1]:
                val = val.reshape(block)
            elif val.shape != block:
                raise MataError(3200, "conformability error")
            base.a[np.ix_(r, c)] = val
            if base.view is not None:
                from .stata_api import write_view
                write_view(self, base, c)
            return
        raise MataError(3000, "invalid lval")

    @staticmethod
    def _check_type(name: str, declared: tuple[str, str], v: MV) -> None:
        el, org = declared
        want = {"real": "real", "string": "string", "complex": "complex", "pointer": "pointer"}.get(el)
        if want and v.t != want and not (want == "complex" and v.t == "real"):
            raise MataError(3250, "type mismatch")
        if org == "scalar" and not v.is_scalar:
            raise MataError(3204, f"{v.orgtype()} found where scalar required")
        if org == "rowvector" and v.rows != 1:
            raise MataError(3203, f"{v.orgtype()} found where rowvector required")
        if org == "colvector" and v.cols != 1:
            raise MataError(3203, f"{v.orgtype()} found where colvector required")
        if org == "vector" and v.rows != 1 and v.cols != 1:
            raise MataError(3203, f"{v.orgtype()} found where vector required")

    # chamadas de função -----------------------------------------------------
    def call(self, node: Call, fr: Frame):
        name = node.func
        # variável que guarda ponteiro para função: (*f)(x) não é sintaxe aceita aqui
        if name in self.funcs:
            return self.call_user(self.funcs[name], node.args, fr)
        from .library import LIBRARY
        entry = LIBRARY.get(name)
        if entry is None:
            if name in fr.vars and fr.vars[name].t == "pointer":
                target = fr.vars[name].scalar().get() if fr.vars[name].scalar() else None
                if isinstance(target, str) and target in self.funcs:
                    return self.call_user(self.funcs[target], node.args, fr)
            raise MataError(3499, f"{name}() not found")
        fn, lo, hi, outs = entry
        n = len(node.args)
        while n and node.args[n - 1] is None:
            n -= 1
        if not lo <= n <= hi:
            raise MataError(3001, f"expected {lo if lo == hi else f'{lo} to {hi}'} arguments "
                                  f"but received {n}")   # VERIFICAR texto
        vals, refs = [], []
        for k, a in enumerate(node.args[:n]):
            if a is None:
                vals.append(None)
                refs.append(None)
                continue
            if k in outs:
                if isinstance(a, Name):
                    refs.append(Ref(fr, a.name))
                    vals.append(fr.vars.get(a.name))
                    continue
                refs.append(None)
                vals.append(self._value(a, fr))
                continue
            refs.append(Ref(fr, a.name) if isinstance(a, Name) else None)
            vals.append(self._value(a, fr))
        try:
            return fn(self, vals, refs)
        except MataError as e:
            if not e.where:
                e.where.append(name + "()")
            raise

    def call_with_values(self, f: FuncDef, values: list[MV]) -> list:
        """Chama f com valores (não expressões) e devolve o valor final de cada
        parâmetro (os "argumentos de saída" de avaliadores, como em optimize)."""
        local = Frame(f.name)
        for p, v in zip(f.params, values):
            local.vars[p.name] = v.copy()
            local.types[p.name] = (p.eltype, p.org)
        self.depth += 1
        self.argc.append(len(values))
        try:
            self.exec(f.body, local)
        except _Return:
            pass
        except MataError as e:
            e.where.append(f.name + "()")
            raise
        finally:
            self.depth -= 1
            self.argc.pop()
        return [local.vars.get(p.name) for p in f.params[:len(values)]]

    def call_user(self, f: FuncDef, args: list, fr: Frame):
        n = len(args)
        required = sum(1 for p in f.params if not p.optional)
        if not required <= n <= len(f.params):
            exp = required if required == len(f.params) else len(f.params)
            raise MataError(3001, f"expected {exp} argument{'s' if exp != 1 else ''} but received {n}")
        if self.depth > 500:
            raise MataError(3000, "stack overflow")
        local = Frame(f.name)
        back: list[tuple[str, str]] = []
        for p, a in zip(f.params, args):
            if a is None:
                continue
            v = self._value(a, fr)
            if p.eltype != "transmorphic" or p.org != "matrix":
                try:
                    self._check_type(p.name, (p.eltype, p.org), v)
                except MataError as e:
                    e.where.insert(0, f.name + "()")
                    raise
            local.vars[p.name] = v
            local.types[p.name] = (p.eltype, p.org)
            if isinstance(a, Name):
                back.append((p.name, a.name))
        for p in f.params[n:]:
            local.types[p.name] = (p.eltype, p.org)
        self.depth += 1
        self.argc.append(n)
        result = None
        try:
            self.exec(f.body, local)
        except _Return as r:
            result = r.value
        except (_Break, _Continue):
            raise MataError(3000, "break or continue outside of a loop")
        except MataError as e:
            e.where.append(f.name + "()")
            raise
        finally:
            self.depth -= 1
            self.argc.pop()
        # passagem por referência: o valor final do parâmetro volta à variável
        for pname, vname in back:
            if pname in local.vars:
                fr.vars[vname] = local.vars[pname]
        if f.rettype == "void":
            return None
        if result is None:
            if f.rettype == "transmorphic" and f.retorg == "matrix":
                return None
            raise MataError(3000, f"{f.name}() returned nothing")   # VERIFICAR
        return result
