"""Analisador do Mata: instruções e expressões ([M-2]).

Precedência (da mais forte para a mais fraca), conforme [M-2] op_intro:

    postfixos   '  [ ]  [| |]  ( )  ++  --  ->  .
    prefixos    -  !  &  *  ++  --
    ^  :^
    *  /  :*  :/  #
    +  -  :+  :-
    ..  ::                      (sequências)
    ,                           (junta colunas)
    \\                          (junta linhas)
    ==  !=  >  >=  <  <=  e as versões com dois-pontos
    &  :&
    |  :|
    &&
    ||
    ?:
    =                           (atribuição, associa à direita)
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .lexer import MataSyntaxError, Token, tokenize

# ---------------------------------------------------------------------------
# AST
# ---------------------------------------------------------------------------


@dataclass
class Num:
    value: float


@dataclass
class Str:
    value: str


@dataclass
class Name:
    name: str


@dataclass
class Unary:
    op: str
    operand: object


@dataclass
class Binary:
    op: str
    left: object
    right: object


@dataclass
class Ternary:
    cond: object
    a: object
    b: object


@dataclass
class Assign:
    target: object
    value: object


@dataclass
class IncDec:
    op: str            # ++ ou --
    target: object
    prefix: bool


@dataclass
class Call:
    func: str
    args: list


@dataclass
class Subscript:
    base: object
    rows: object       # None = todas (.), ou expressão
    cols: object       # None se subscrito de um índice só
    single: bool = False


@dataclass
class RangeSub:
    base: object
    spec: object       # expressão 2x2 (ou 1x2) com os cantos


@dataclass
class Transpose:
    operand: object


@dataclass
class Member:
    base: object
    name: str
    arrow: bool


@dataclass
class AddrOf:
    operand: object


@dataclass
class Deref:
    operand: object


@dataclass
class Missing:
    """O "." usado como subscrito (todas as linhas ou colunas)."""


# instruções

@dataclass
class ExprStmt:
    expr: object


@dataclass
class Block:
    body: list


@dataclass
class If:
    cond: object
    then: object
    other: object = None


@dataclass
class For:
    init: object
    cond: object
    step: object
    body: object


@dataclass
class While:
    cond: object
    body: object


@dataclass
class DoWhile:
    body: object
    cond: object


@dataclass
class Break:
    pass


@dataclass
class Continue:
    pass


@dataclass
class Return:
    value: object = None


@dataclass
class Decl:
    eltype: str
    org: str
    names: list
    external: bool = False


@dataclass
class Param:
    name: str
    eltype: str = "transmorphic"
    org: str = "matrix"
    optional: bool = False


@dataclass
class FuncDef:
    name: str
    rettype: str
    retorg: str
    params: list
    body: Block
    source: str = ""


@dataclass
class StructDef:
    name: str
    fields: list = field(default_factory=list)    # (eltype, org, name)


@dataclass
class Empty:
    pass


ELTYPES = {"real", "complex", "string", "pointer", "transmorphic", "numeric", "struct", "class", "void",
           "function"}
ORGS = {"scalar", "vector", "rowvector", "colvector", "matrix"}


# ---------------------------------------------------------------------------
# Analisador
# ---------------------------------------------------------------------------

class Parser:
    def __init__(self, text: str):
        self.text = text
        self.toks: list[Token] = tokenize(text)
        self.i = 0
        self.no_comma = 0
        self.in_function = False

    # utilidades -----------------------------------------------------------
    @property
    def tok(self) -> Token:
        return self.toks[self.i]

    def peek(self, k: int = 1) -> Token:
        return self.toks[min(self.i + k, len(self.toks) - 1)]

    def at(self, *ops: str) -> bool:
        t = self.tok
        return t.kind == "op" and t.value in ops

    def at_id(self, *names: str) -> bool:
        return self.tok.kind == "id" and self.tok.value in names

    def take(self) -> Token:
        t = self.tok
        self.i += 1
        return t

    def expect(self, op: str) -> None:
        if self.tok.kind == "eof":
            raise MataSyntaxError(f"'{op}' expected", incomplete=True)
        if not self.at(op):
            raise MataSyntaxError("invalid expression")
        self.i += 1

    def skip_semis(self) -> None:
        while self.at(";"):
            self.i += 1

    # programa --------------------------------------------------------------
    def parse_program(self) -> list:
        out = []
        self.skip_semis()
        while self.tok.kind != "eof":
            out.append(self.statement(top=True))
            self.skip_semis()
        return out

    # instruções --------------------------------------------------------------
    def statement(self, top: bool = False):
        t = self.tok
        if t.kind == "eof":
            raise MataSyntaxError("statement expected", incomplete=True)
        if self.at("{"):
            return self.block()
        if self.at("}"):
            raise MataSyntaxError("expression invalid")
        if self.at(";"):
            self.take()
            return Empty()
        if t.kind == "id":
            kw = t.value
            if kw == "if":
                self.take()
                self.expect("(")
                cond = self.expression()
                self.expect(")")
                then = self.statement()
                other = None
                save = self.i
                self.skip_semis()
                if self.at_id("else"):
                    self.take()
                    other = self.statement()
                else:
                    self.i = save
                return If(cond, then, other)
            if kw == "for":
                self.take()
                self.expect("(")
                init = None if self.at(";") else self.expression()
                self.expect(";")
                cond = None if self.at(";") else self.expression()
                self.expect(";")
                step = None if self.at(")") else self.expression()
                self.expect(")")
                return For(init, cond, step, self.statement())
            if kw == "while":
                self.take()
                self.expect("(")
                cond = self.expression()
                self.expect(")")
                return While(cond, self.statement())
            if kw == "do":
                self.take()
                body = self.statement()
                self.skip_semis()
                if self.tok.kind == "eof":
                    raise MataSyntaxError("while expected", incomplete=True)
                if not self.at_id("while"):
                    raise MataSyntaxError("while expected")
                self.take()
                self.expect("(")
                cond = self.expression()
                self.expect(")")
                return DoWhile(body, cond)
            if kw == "break":
                self.take()
                return Break()
            if kw == "continue":
                self.take()
                return Continue()
            if kw == "return":
                if not self.in_function:
                    # fora de função (observado no Stata 14)
                    raise MataSyntaxError("'return' found where almost anything else expected")
                self.take()
                if self.at("("):
                    self.take()
                    if self.at(")"):
                        self.take()
                        return Return(None)
                    v = self.expression()
                    self.expect(")")
                    return Return(v)
                if self.tok.kind == "eof" or self.at(";", "}"):
                    return Return(None)
                return Return(self.expression())
            if kw == "else":
                raise MataSyntaxError("else without if")
            if kw == "struct" and self.peek().kind == "id" and self.peek(2).kind == "op" \
                    and self.peek(2).value == "{" and top:
                return self.struct_def()
            if self._looks_like_function(top):
                return self.func_def()
            save = self.i
            decl = self._try_decl()
            if decl is not None:
                if not self.in_function:
                    # fora de função não há declarações (observado no Stata 14):
                    # no nível de cima o Mata espera o "(" de uma função
                    if not top:
                        raise MataSyntaxError(f"'{kw}' found where almost anything else expected")
                    self.i = save
                    self._type_spec()
                    self.take()
                    if self.tok.kind == "eof":
                        raise MataSyntaxError("'(' expected", incomplete=True)
                    raise MataSyntaxError(f"'{self.tok.value}' found where '(' expected")
                return decl
        return ExprStmt(self.expression())

    def block(self) -> Block:
        self.expect("{")
        body = []
        while True:
            self.skip_semis()
            if self.tok.kind == "eof":
                raise MataSyntaxError("} expected", incomplete=True)
            if self.at("}"):
                self.take()
                return Block(body)
            body.append(self.statement())

    # declarações -----------------------------------------------------------
    def _type_spec(self) -> tuple[str, str] | None:
        """Lê [eltype] [org] (ex.: real scalar, string matrix, pointer(real) scalar,
        struct nome scalar). Devolve None se não há tipo aqui."""
        save = self.i
        el, org = "", ""
        if self.at_id("struct", "class"):
            self.take()
            if self.tok.kind != "id":
                self.i = save
                return None
            el = "struct " + self.take().value
        elif self.at_id(*(ELTYPES - {"struct", "class"})):
            el = self.take().value
            if el == "pointer" and self.at("("):
                depth = 0
                while True:            # pointer(real matrix) — o alvo não importa
                    if self.at("("):
                        depth += 1
                    elif self.at(")"):
                        depth -= 1
                    self.take()
                    if depth == 0:
                        break
        if self.at_id(*ORGS):
            org = self.take().value
        if not el and not org:
            self.i = save
            return None
        return el or "transmorphic", org or "matrix"

    def _looks_like_function(self, top: bool) -> bool:
        if not top:
            return False
        save = self.i
        try:
            if self.at_id("function"):
                self.take()
            else:
                spec = self._type_spec()
                if spec is None:
                    return False
            if self.tok.kind != "id":
                return False
            self.take()
            return self.at("(")
        finally:
            self.i = save

    def _try_decl(self):
        save = self.i
        external = False
        if self.at_id("external"):
            self.take()
            external = True
        spec = self._type_spec()
        if spec is None or self.tok.kind != "id":
            self.i = save
            return None
        # "real x" é declaração; "real(x)" seria chamada de função
        names = []
        while True:
            if self.tok.kind != "id":
                raise MataSyntaxError("invalid declaration")
            names.append(self.take().value)
            if self.at("["):            # real matrix A[2,2]? (não usado no Mata moderno)
                raise MataSyntaxError("invalid declaration")
            if not self.at(","):
                break
            self.take()
        return Decl(spec[0], spec[1], names, external)

    def func_def(self) -> FuncDef:
        start = self.tok.pos
        if self.at_id("function"):
            self.take()
            rt, ro = "transmorphic", "matrix"
        else:
            rt, ro = self._type_spec()
        name = self.take().value
        self.expect("(")
        params = []
        optional = False
        while not self.at(")"):
            if self.tok.kind == "eof":
                raise MataSyntaxError(") expected", incomplete=True)
            if self.at("|"):
                self.take()
                optional = True
                continue
            spec = self._type_spec() or ("transmorphic", "matrix")
            if self.tok.kind != "id":
                raise MataSyntaxError("invalid function declaration")
            params.append(Param(self.take().value, spec[0], spec[1], optional))
            if self.at(","):
                self.take()
        self.take()
        self.skip_semis()
        if self.tok.kind == "eof":
            raise MataSyntaxError("{ expected", incomplete=True)
        self.in_function = True
        try:
            body = self.block()
        finally:
            self.in_function = False
        src = self.text[start:self.toks[self.i - 1].pos + 1]
        return FuncDef(name, rt, ro, params, body, src)

    def struct_def(self) -> StructDef:
        self.take()
        name = self.take().value
        self.expect("{")
        fields = []
        while True:
            self.skip_semis()
            if self.tok.kind == "eof":
                raise MataSyntaxError("} expected", incomplete=True)
            if self.at("}"):
                self.take()
                break
            d = self._try_decl()
            if d is None:
                raise MataSyntaxError("invalid struct definition")
            for nm in d.names:
                fields.append((d.eltype, d.org, nm))
        return StructDef(name, fields)

    # expressões ----------------------------------------------------------------
    def expression(self):
        return self.assignment()

    def assignment(self):
        left = self.ternary()
        if self.at("="):
            self.take()
            if not isinstance(left, (Name, Subscript, RangeSub, Member, Deref)):
                raise MataSyntaxError("invalid lval")
            return Assign(left, self.assignment())
        return left

    def ternary(self):
        cond = self.oror()
        if self.at("?"):
            self.take()
            a = self.ternary()
            self.expect(":")
            b = self.ternary()
            return Ternary(cond, a, b)
        return cond

    def _binary(self, sub, ops):
        left = sub()
        while self.tok.kind == "op" and self.tok.value in ops:
            op = self.take().value
            if self.tok.kind == "eof":
                raise MataSyntaxError("invalid expression", incomplete=True)
            right = sub()
            _literal_types(op, left, right)
            left = Binary(op, left, right)
        return left

    def oror(self):
        return self._binary(self.andand, ("||",))

    def andand(self):
        return self._binary(self.orop, ("&&",))

    def orop(self):
        return self._binary(self.andop, ("|", ":|"))

    def andop(self):
        return self._binary(self.relational, ("&", ":&"))

    def relational(self):
        return self._binary(self.rowjoin, ("==", "!=", ">", ">=", "<", "<=",
                                           ":==", ":!=", ":>", ":>=", ":<", ":<="))

    def rowjoin(self):
        return self._binary(self.coljoin, ("\\",))

    def coljoin(self):
        if self.no_comma:
            return self.range_()
        return self._binary(self.range_, (",",))

    def range_(self):
        return self._binary(self.additive, ("..", "::"))

    def additive(self):
        return self._binary(self.multiplicative, ("+", "-", ":+", ":-"))

    def multiplicative(self):
        return self._binary(self.unary, ("*", "/", ":*", ":/", "#"))

    def unary(self):
        if self.at("-", "!", "+"):
            op = self.take().value
            operand = self.unary()
            return operand if op == "+" else Unary(op, operand)
        if self.at("&"):
            self.take()
            return AddrOf(self.unary())
        if self.at("*"):
            self.take()
            return Deref(self.unary())
        if self.at("++", "--"):
            op = self.take().value
            return IncDec(op, self.unary(), prefix=True)
        return self.power()

    def power(self):
        left = self.postfix()
        while self.at("^", ":^"):
            op = self.take().value
            if self.at("-", "+"):          # 2^-1
                sign = self.take().value
                right = self.postfix()
                right = Unary("-", right) if sign == "-" else right
            else:
                right = self.postfix()
            left = Binary(op, left, right)
        return left

    def postfix(self):
        node = self.primary()
        while True:
            if self.at("'"):
                self.take()
                node = Transpose(node)
            elif self.at("["):
                self.take()
                self.no_comma += 1
                try:
                    r = self._sub_index()
                    if self.at(","):
                        self.take()
                        c = self._sub_index()
                        self.expect("]")
                        node = Subscript(node, r, c)
                    else:
                        self.expect("]")
                        node = Subscript(node, r, None, single=True)
                finally:
                    self.no_comma -= 1
            elif self.at("[|"):
                self.take()
                save = self.no_comma
                self.no_comma = 0
                spec = self.expression()
                self.no_comma = save
                self.expect("|]")
                node = RangeSub(node, spec)
            elif self.at("++", "--") and isinstance(node, (Name, Subscript, Member)):
                node = IncDec(self.take().value, node, prefix=False)
            elif self.at("->"):
                self.take()
                node = Member(node, self.take().value, arrow=True)
            elif self.at(".") and self.peek().kind == "id":
                self.take()
                node = Member(node, self.take().value, arrow=False)
            else:
                return node

    def _sub_index(self):
        if self.tok.kind == "num" and self.tok.value >= 8.98846567431158e307 \
                and self.peek().kind == "op" and self.peek().value in (",", "]"):
            self.take()
            return None
        return self.expression()

    def primary(self):
        t = self.tok
        if t.kind == "num":
            self.take()
            return Num(float(t.value))
        if t.kind == "str":
            self.take()
            return Str(str(t.value))
        if t.kind == "id":
            self.take()
            if self.at("("):
                self.take()
                args = []
                self.no_comma += 1
                try:
                    while not self.at(")"):
                        if self.tok.kind == "eof":
                            raise MataSyntaxError(") expected", incomplete=True)
                        if self.at(","):            # argumento omitido: f(a, , b)
                            self.take()
                            args.append(None)
                            continue
                        args.append(self.expression())
                        if self.at(","):
                            self.take()
                            if self.at(")"):
                                args.append(None)
                        elif not self.at(")"):
                            raise MataSyntaxError("invalid expression")
                finally:
                    self.no_comma -= 1
                self.take()
                return Call(t.value, args)
            return Name(t.value)
        if self.at("("):
            self.take()
            save = self.no_comma
            self.no_comma = 0
            try:
                if self.at(")"):
                    raise MataSyntaxError("invalid expression")
                e = self.expression()
            finally:
                self.no_comma = save
            self.expect(")")
            return e
        if t.kind == "eof":
            raise MataSyntaxError("invalid expression", incomplete=True)
        raise MataSyntaxError("invalid expression")


_ARITH = {"+", "-", "*", "/", "^", ":+", ":-", ":*", ":/", ":^", "#"}


def _literal_types(op: str, a, b) -> None:
    """Constantes de tipos diferentes já dão erro na compilação:
    "a" + 1 → type mismatch:  string + real not allowed (Stata 14)."""
    if op not in _ARITH:
        return
    kinds = []
    for x in (a, b):
        kinds.append("string" if isinstance(x, Str) else "real" if isinstance(x, Num) else None)
    if None in kinds or kinds[0] == kinds[1]:
        return
    if op in ("*", ":*") and "real" in kinds:      # "ab" * 3 é repetição
        return
    raise MataSyntaxError(f"type mismatch:  {kinds[0]} {op} {kinds[1]} not allowed")


def parse(text: str) -> list:
    return Parser(text).parse_program()
