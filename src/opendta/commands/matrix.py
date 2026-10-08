"""matrix ([P] matrix): definição, operadores, funções, listagem e nomes.

Matrizes ficam em s.matrices (e em r()/e() como resultados). Os valores são
float64 com os códigos de missing do OpenDTA; nas contas, missing vira NaN
e volta a '.' no fim (missing estendidos viram '.').

Expressões de matriz:

    A + B   A - B   A * B (produto ou escalar)   A # B (Kronecker)
    A / k   -A   A'   (A, B) junta colunas   (A \\ B) junta linhas
    J(r,c,z) I(n) inv() invsym() cholesky() diag() vecdiag() hadamard()
    corr() vec() nullmat() matuniform(r,c)   A[i..j, k...]   r(nome) e(nome)

Nas expressões comuns: A[i,j], el(A,i,j), rowsof(), colsof(), trace(),
det(), issymmetric(), rownumb(), colnumb(), matmissing(), mreldif().
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value, parse_format
from ..lang.syntax import match_options
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@dataclass
class Matrix:
    data: np.ndarray
    rownames: list[str] = field(default_factory=list)
    colnames: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.data = np.atleast_2d(np.asarray(self.data, dtype=np.float64))
        r, c = self.data.shape
        if len(self.rownames) != r:
            self.rownames = [f"r{i}" for i in range(1, r + 1)]
        if len(self.colnames) != c:
            self.colnames = [f"c{j}" for j in range(1, c + 1)]

    @property
    def rows(self) -> int:
        return int(self.data.shape[0])

    @property
    def cols(self) -> int:
        return int(self.data.shape[1])

    def copy(self) -> "Matrix":
        return Matrix(self.data.copy(), list(self.rownames), list(self.colnames))


def _store(s: "Session") -> dict[str, Matrix]:
    if not hasattr(s, "matrices"):
        s.matrices = {}
    return s.matrices


def _nan(a: np.ndarray) -> np.ndarray:
    out = np.array(a, dtype=np.float64)
    out[out >= M.SYSMISS] = np.nan
    return out


def _miss(a: np.ndarray) -> np.ndarray:
    out = np.array(a, dtype=np.float64)
    bad = ~np.isfinite(out) | (np.abs(out) >= M.SYSMISS)
    out[bad] = M.SYSMISS
    return out


def get_matrix(s: "Session", name: str) -> Matrix:
    name = name.strip()
    m = re.fullmatch(r"([re])\((\w+)\)", name)
    if m:
        store = s.r if m.group(1) == "r" else s.e
        v = store.get(m.group(2))
        if isinstance(v, Matrix):
            return v
        raise StataError(111, f"matrix {name} not found")
    mat = _store(s).get(name)
    if mat is None:
        raise StataError(111, f"matrix {name} not found")
    return mat


def matrix_names(s: "Session", name: str, which: str) -> str:
    mat = get_matrix(s, name)
    if which in ("rownames", "rowfullnames"):
        return " ".join(mat.rownames)
    if which in ("colnames", "colfullnames"):
        return " ".join(mat.colnames)
    return " ".join("_" for _ in (mat.rownames if which == "roweq" else mat.colnames))


# ---------------------------------------------------------------------------
# expressões de matriz
# ---------------------------------------------------------------------------

_TOKEN = re.compile(r"""
    (?P<num>\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)
  | (?P<res>[re]\(\w+\))
  | (?P<name>[A-Za-z_][A-Za-z0-9_]*)
  | (?P<range>\.\.\.|\.\.)
  | (?P<miss>\.[a-z]?)
  | (?P<str>"[^"]*")
  | (?P<op>[-+*#/'\\,()\[\]^])
  | (?P<ws>\s+)
""", re.X)

_MATFUNCS = {"J", "I", "inv", "invsym", "cholesky", "diag", "vecdiag", "hadamard", "corr", "vec",
             "nullmat", "matuniform", "get"}


class _Parser:
    def __init__(self, s: "Session", text: str):
        self.s = s
        self.text = text
        self.toks: list[tuple[str, str, int]] = []
        pos = 0
        while pos < len(text):
            m = _TOKEN.match(text, pos)
            if not m:
                raise StataError(198, "invalid syntax")
            kind = m.lastgroup
            if kind != "ws":
                self.toks.append((kind, m.group(0), pos))
            pos = m.end()
        self.toks.append(("end", "", len(text)))
        self.i = 0

    @property
    def tok(self):
        return self.toks[self.i]

    def take(self):
        t = self.toks[self.i]
        self.i += 1
        return t

    def is_op(self, v: str) -> bool:
        return self.tok[0] == "op" and self.tok[1] == v

    def expect(self, v: str) -> None:
        if not self.is_op(v):
            raise StataError(198, "invalid syntax")
        self.i += 1

    def parse(self) -> Matrix:
        m = self.joined()             # A \ B e A, B também valem fora de parênteses
        if self.tok[0] != "end":
            raise StataError(198, "invalid syntax")
        return m

    def expr(self) -> Matrix:
        left = self.term()
        while self.is_op("+") or self.is_op("-"):
            op = self.take()[1]
            right = self.term()
            left = _add(left, right, op)
        return left

    def term(self) -> Matrix:
        left = self.unary()
        while self.is_op("*") or self.is_op("#") or self.is_op("/"):
            op = self.take()[1]
            right = self.unary()
            if op == "*":
                left = _mul(left, right)
            elif op == "#":
                d = _miss(np.kron(_nan(left.data), _nan(right.data)))
                left = Matrix(d)
            else:
                if right.data.shape != (1, 1):
                    raise StataError(503, "conformability error")
                left = Matrix(_miss(_nan(left.data) / _nan(right.data)[0, 0]),
                              left.rownames, left.colnames)
        return left

    def unary(self) -> Matrix:
        if self.is_op("-"):
            self.take()
            m = self.unary()
            return Matrix(_miss(-_nan(m.data)), m.rownames, m.colnames)
        if self.is_op("+"):
            self.take()
            return self.unary()
        return self.power()

    def power(self) -> Matrix:
        """x^y só entre escalares (1 x 1), como dentro de (`i', `i'^2)."""
        base = self.postfix()
        if self.is_op("^"):
            self.take()
            exp = self.unary()
            if base.data.shape != (1, 1) or exp.data.shape != (1, 1):
                raise StataError(503, "conformability error")
            x = self.s.eval(f"({float(base.data[0, 0])!r})^({float(exp.data[0, 0])!r})") \
                if base.data[0, 0] < M.SYSMISS and exp.data[0, 0] < M.SYSMISS else M.SYSMISS
            return Matrix(np.array([[float(x)]]))
        return base

    def postfix(self) -> Matrix:
        m = self.primary()
        while True:
            if self.is_op("'"):
                self.take()
                m = Matrix(m.data.T.copy(), list(m.colnames), list(m.rownames))
            elif self.is_op("["):
                m = self.submatrix(m)
            else:
                return m

    def primary(self) -> Matrix:
        kind, val, pos = self.tok
        if self.is_op("("):
            self.take()
            m = self.joined()
            self.expect(")")
            return m
        if kind == "num":
            self.take()
            return Matrix(np.array([[float(val)]]))
        if kind == "miss":
            self.take()
            return Matrix(np.array([[M.missing_code(val)]]))
        if kind == "res":
            self.take()
            return get_matrix(self.s, val).copy()
        if kind == "name":
            nxt = self.toks[self.i + 1]
            if nxt[0] == "op" and nxt[1] == "(":
                return self.function()
            if val in _store(self.s):
                self.take()
                return _store(self.s)[val].copy()
            # escalar com nome (scalar, variável na 1ª obs.)
            self.take()
            x = self.s.eval(val)
            if isinstance(x, str):
                raise StataError(109, "type mismatch")
            return Matrix(np.array([[float(x)]]))
        raise StataError(198, "invalid syntax")

    def joined(self) -> Matrix:
        """Dentro de parênteses: ',' junta colunas e '\\' junta linhas.
        nullmat() de matriz inexistente some da junção."""
        rows = [self.row()]
        while self.is_op("\\"):
            self.take()
            rows.append(self.row())
        rows = [r for r in rows if not isinstance(r, _Empty)] or [_Empty()]
        if len(rows) == 1:
            return rows[0]
        cols = {r.cols for r in rows}
        if len(cols) != 1:
            raise StataError(503, "conformability error")
        names = [n for r in rows for n in r.rownames]
        return Matrix(np.vstack([r.data for r in rows]),
                      names if len(set(names)) == len(names) and not _default(names, "r") else [],
                      rows[0].colnames)

    def row(self) -> Matrix:
        parts = [self.expr()]
        while self.is_op(","):
            self.take()
            parts.append(self.expr())
        parts = [p for p in parts if not isinstance(p, _Empty)] or [_Empty()]
        if len(parts) == 1:
            return parts[0]
        nrows = {p.rows for p in parts}
        if len(nrows) != 1:
            raise StataError(503, "conformability error")
        names = [n for p in parts for n in p.colnames]
        return Matrix(np.hstack([p.data for p in parts]), parts[0].rownames,
                      names if len(set(names)) == len(names) and not _default(names, "c") else [])

    def function(self) -> Matrix:
        name = self.take()[1]
        self.expect("(")
        start = self.tok[2]
        depth = 1
        args_raw: list[str] = []
        arg_start = start
        while True:
            kind, val, pos = self.tok
            if kind == "end":
                raise StataError(132, "too many '(' or '['")
            if kind == "op" and val == "(":
                depth += 1
            elif kind == "op" and val == ")":
                depth -= 1
                if depth == 0:
                    args_raw.append(self.text[arg_start:pos])
                    self.take()
                    break
            elif kind == "op" and val == "," and depth == 1:
                args_raw.append(self.text[arg_start:pos])
                arg_start = pos + 1
            self.take()
        args = [a.strip() for a in args_raw if a.strip() or len(args_raw) > 1]
        return _call(self.s, name, args)

    def submatrix(self, m: Matrix) -> Matrix:
        self.expect("[")
        start = self.tok[2]
        depth = 1
        while True:
            kind, val, pos = self.tok
            if kind == "end":
                raise StataError(198, "invalid syntax")
            if kind == "op" and val == "[":
                depth += 1
            if kind == "op" and val == "]":
                depth -= 1
                if depth == 0:
                    inner = self.text[start:pos]
                    self.take()
                    break
            self.take()
        parts = _split_top(inner, ",")
        if len(parts) != 2:
            raise StataError(198, "invalid syntax")
        ri = _index_range(self.s, parts[0], m.rownames)
        ci = _index_range(self.s, parts[1], m.colnames)
        return Matrix(m.data[np.ix_(ri, ci)], [m.rownames[i] for i in ri],
                      [m.colnames[j] for j in ci])


def _default(names: list[str], prefix: str) -> bool:
    return all(re.fullmatch(prefix + r"\d+", n) for n in names)


def _split_top(text: str, sep: str) -> list[str]:
    out, depth, cur, q = [], 0, "", False
    for ch in text:
        if ch == '"':
            q = not q
        if not q:
            if ch in "([":
                depth += 1
            elif ch in ")]":
                depth -= 1
            elif ch == sep and depth == 0:
                out.append(cur)
                cur = ""
                continue
        cur += ch
    out.append(cur)
    return out


def _index_range(s: "Session", text: str, names: list[str]) -> list[int]:
    t = text.strip()
    n = len(names)

    def one(x: str) -> int:
        x = x.strip()
        if x.startswith('"'):
            label = x.strip('"')
            if label not in names:
                raise StataError(111, f"{label} not found")
            return names.index(label) + 1
        v = s.eval(x)
        return int(v)

    if "..." in t:
        a = one(t.split("...")[0])
        b = n
    elif ".." in t:
        a_txt, b_txt = t.split("..", 1)
        a, b = one(a_txt), one(b_txt)
    else:
        a = b = one(t)
    if not (1 <= a <= n and 1 <= b <= n and a <= b):
        raise StataError(503, "conformability error")
    return list(range(a - 1, b))


def _add(a: Matrix, b: Matrix, op: str) -> Matrix:
    if a.data.shape != b.data.shape:
        raise StataError(503, "conformability error")
    d = _nan(a.data) + _nan(b.data) if op == "+" else _nan(a.data) - _nan(b.data)
    return Matrix(_miss(d), a.rownames, a.colnames)


def _mul(a: Matrix, b: Matrix) -> Matrix:
    if a.data.shape == (1, 1) and b.data.shape != (1, 1):
        return Matrix(_miss(_nan(b.data) * _nan(a.data)[0, 0]), b.rownames, b.colnames)
    if b.data.shape == (1, 1) and a.data.shape != (1, 1):
        return Matrix(_miss(_nan(a.data) * _nan(b.data)[0, 0]), a.rownames, a.colnames)
    if a.cols != b.rows:
        raise StataError(503, "conformability error")
    return Matrix(_miss(_nan(a.data) @ _nan(b.data)), a.rownames, b.colnames)


def _scalar(s: "Session", text: str) -> float:
    v = s.eval(text)
    if isinstance(v, str):
        raise StataError(109, "type mismatch")
    return float(v)


def _mat_arg(s: "Session", text: str) -> Matrix:
    return _Parser(s, text).parse()


def _call(s: "Session", name: str, args: list[str]) -> Matrix:
    if name not in _MATFUNCS:
        # função escalar: vira matriz 1 x 1
        x = _scalar(s, f"{name}({', '.join(args)})")
        return Matrix(np.array([[x]]))
    if name == "J":
        if len(args) != 3:
            raise StataError(198, "invalid syntax")
        r, c = int(_scalar(s, args[0])), int(_scalar(s, args[1]))
        return Matrix(np.full((r, c), _scalar(s, args[2])))
    if name == "I":
        n = int(_scalar(s, args[0]))
        return Matrix(np.eye(n))
    if name == "matuniform":
        r, c = int(_scalar(s, args[0])), int(_scalar(s, args[1]))
        return Matrix(np.random.default_rng().random((r, c)))
    if name == "nullmat":
        nm = args[0].strip()
        if nm in _store(s):
            return _store(s)[nm].copy()
        return _Empty()
    a = _mat_arg(s, args[0])
    if name in ("inv", "invsym"):
        if a.rows != a.cols:
            raise StataError(505, "matrix not symmetric" if name == "invsym" else "matrix not square")
        x = _nan(a.data)
        try:
            d = np.linalg.inv(x)
        except np.linalg.LinAlgError:
            if name == "inv":
                raise StataError(504, "matrix has missing values" if np.isnan(x).any()
                                 else "matrix not positive definite")   # VERIFICAR
            d = np.linalg.pinv(x)   # VERIFICAR: o Stata zera linhas/colunas dependentes
        return Matrix(_miss(d), list(a.colnames), list(a.rownames))
    if name == "cholesky":
        try:
            d = np.linalg.cholesky(_nan(a.data))
        except np.linalg.LinAlgError:
            raise StataError(506, "matrix not positive definite")
        return Matrix(d, a.rownames, a.colnames)
    if name == "diag":
        v = a.data.ravel()
        names = a.colnames if a.rows == 1 else a.rownames
        return Matrix(np.diag(v), names, names)
    if name == "vecdiag":
        if a.rows != a.cols:
            raise StataError(505, "matrix not square")
        return Matrix(np.diag(a.data)[None, :], ["r1"], list(a.colnames))
    if name == "hadamard":
        b = _mat_arg(s, args[1])
        if a.data.shape != b.data.shape:
            raise StataError(503, "conformability error")
        return Matrix(_miss(_nan(a.data) * _nan(b.data)), a.rownames, a.colnames)
    if name == "corr":
        d = np.sqrt(np.diag(_nan(a.data)))
        return Matrix(_miss(_nan(a.data) / np.outer(d, d)), a.rownames, a.colnames)
    if name == "vec":
        col = a.data.T.reshape(-1, 1)
        names = [f"{c}:{r}" for c in a.colnames for r in a.rownames]
        return Matrix(col, names, ["c1"])
    raise StataError(133, f"unknown function {name}()")


class _Empty(Matrix):
    """nullmat() de matriz inexistente: elemento neutro da junção."""

    def __init__(self):
        self.data = np.zeros((0, 0))
        self.rownames, self.colnames = [], []


def evaluate_matrix(s: "Session", text: str) -> Matrix:
    m = _Parser(s, text).parse()
    if isinstance(m, _Empty):
        raise StataError(503, "conformability error")
    return m


# ---------------------------------------------------------------------------
# funções escalares sobre matrizes (usadas pelo avaliador de expressões)
# ---------------------------------------------------------------------------

MATRIX_SCALAR_FUNCS = {"rowsof", "colsof", "el", "trace", "det", "issymmetric", "rownumb",
                       "colnumb", "matmissing", "mreldif"}


def matrix_scalar(s: "Session", name: str, args: list, eval_arg) -> float:
    """args: nós da expressão; o 1º (e o 2º em mreldif) é nome de matriz."""
    from ..lang.expr import Name, Call

    def mat(node) -> Matrix:
        if isinstance(node, Name):
            return get_matrix(s, node.name)
        if isinstance(node, Call) and node.name in ("r", "e"):
            return get_matrix(s, f"{node.name}({node.raw})")
        raise StataError(198, "invalid syntax")

    a = mat(args[0])
    if name == "rowsof":
        return float(a.rows)
    if name == "colsof":
        return float(a.cols)
    if name == "el":
        i, j = int(eval_arg(args[1])), int(eval_arg(args[2]))
        if 1 <= i <= a.rows and 1 <= j <= a.cols:
            return float(a.data[i - 1, j - 1])
        return M.SYSMISS
    if name == "trace":
        if a.rows != a.cols:
            raise StataError(505, "matrix not square")
        return float(_miss(np.array([np.trace(_nan(a.data))]))[0])
    if name == "det":
        if a.rows != a.cols:
            raise StataError(505, "matrix not square")
        return float(_miss(np.array([np.linalg.det(_nan(a.data))]))[0])
    if name == "issymmetric":
        return 1.0 if a.rows == a.cols and np.array_equal(a.data, a.data.T) else 0.0
    if name in ("rownumb", "colnumb"):
        label = eval_arg(args[1])
        names = a.rownames if name == "rownumb" else a.colnames
        return float(names.index(label) + 1) if label in names else M.SYSMISS
    if name == "matmissing":
        return 1.0 if (a.data >= M.SYSMISS).any() else 0.0
    if name == "mreldif":
        b = mat(args[1])
        if a.data.shape != b.data.shape:
            return M.SYSMISS
        x, y = _nan(a.data), _nan(b.data)
        return float(np.nanmax(np.abs(x - y) / (np.abs(y) + 1)))
    raise StataError(133, f"unknown function {name}()")


def matrix_element(s: "Session", name: str, i: float, j: float) -> float:
    a = get_matrix(s, name)
    i, j = int(i), int(j)
    if 1 <= i <= a.rows and 1 <= j <= a.cols:
        return float(a.data[i - 1, j - 1])
    raise StataError(503, "conformability error")   # VERIFICAR


# ---------------------------------------------------------------------------
# matrix list
# ---------------------------------------------------------------------------

def list_matrix(s: "Session", title: str, mat: Matrix, fmt: str = "%10.0g",
                *, names: bool = True, header: bool = True, show_title: str = "") -> None:
    out = s.output
    f = parse_format(fmt)
    sym = mat.rows == mat.cols and mat.rows > 1 and np.array_equal(mat.data, mat.data.T)
    cells = [[format_value(float(x), f, pad=False).strip() if x < M.SYSMISS
              else M.missing_name(float(x)) for x in row] for row in mat.data]
    rw = max([len(n) for n in mat.rownames] + [0]) if names else 0
    widths = []
    for j, cname in enumerate(mat.colnames):
        col = [cells[i][j] for i in range(mat.rows) if not sym or i >= j]
        widths.append(max([len(cname) if names else 0] + [len(c) for c in col]))
    if header:
        lead = "symmetric " if sym else ""
        extra = f":  {show_title}" if show_title else ""
        out.write(f"\n{lead}{title}[{mat.rows},{mat.cols}]{extra}\n", "text")
    if names:
        line = " " * rw + "".join("  " + c.rjust(w) for c, w in zip(mat.colnames, widths))
        out.write(line + "\n", "text")
    for i in range(mat.rows):
        out.write((mat.rownames[i].ljust(rw) if names else ""), "text")
        last = i + 1 if sym else mat.cols
        out.write("".join("  " + cells[i][j].rjust(widths[j]) for j in range(last)) + "\n", "result")


# ---------------------------------------------------------------------------
# comando matrix
# ---------------------------------------------------------------------------

_NAME = re.compile(r"[A-Za-z_][A-Za-z0-9_]{0,31}")


def _define(s: "Session", text: str) -> None:
    m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*(\[(.*?)\])?\s*=\s*(.*)$", text, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    name, sub, rhs = m.group(1), m.group(3), m.group(4).strip()
    if sub is not None:
        # matrix A[i,j] = exp
        mat = get_matrix(s, name)
        parts = _split_top(sub, ",")
        if len(parts) != 2:
            raise StataError(198, "invalid syntax")
        i, j = int(_scalar(s, parts[0])), int(_scalar(s, parts[1]))
        if not (1 <= i <= mat.rows and 1 <= j <= mat.cols):
            raise StataError(503, "conformability error")
        x = _scalar(s, rhs)
        mat.data[i - 1, j - 1] = x
        return
    _store(s)[name] = evaluate_matrix(s, rhs)


def _input(s: "Session", text: str) -> None:
    m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*\((.*)\)\s*$", text, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    rows = []
    for row in m.group(2).split("\\"):
        vals = []
        for cell in row.replace(",", " ").split():
            if cell == "." or re.fullmatch(r"\.[a-z]", cell):
                vals.append(M.missing_code(cell))
            else:
                try:
                    vals.append(float(cell))
                except ValueError:
                    raise StataError(198, "invalid syntax")
        rows.append(vals)
    if len({len(r) for r in rows}) != 1:
        raise StataError(503, "conformability error")
    _store(s)[m.group(1)] = Matrix(np.array(rows))


def _names_cmd(s: "Session", which: str, text: str) -> None:
    m = re.match(r"\s*(\S+)\s*=\s*(.*)$", text, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    mat = get_matrix(s, m.group(1))
    names = m.group(2).split()
    n = mat.rows if which.startswith("row") else mat.cols
    if len(names) != n:
        raise StataError(503, "conformability error")
    if which.endswith("eq"):
        return                                   # equações: VERIFICAR (fase 5)
    if which.startswith("row"):
        mat.rownames = names
    else:
        mat.colnames = names


@command("matrix", "mat")
def cmd_matrix(s: "Session", args: str) -> None:
    t = args.strip()
    sub, _, rest = t.partition(" ")
    if re.match(r"\s*(=|\[)", rest) or "=" in sub or "[" in sub:
        _define(s, t)                 # matrix d = ...: "d" é o nome, não define
        return
    if sub in ("define", "def", "de", "defi", "defin", "d"):
        _define(s, rest)
        return
    if sub in ("input", "in", "inp", "inpu"):
        _input(s, rest)
        return
    if sub in ("list", "l", "li", "lis"):
        head, comma, opts = rest.partition(",")
        o = match_options(opts, {"format": 3, "title": 3, "nonames": 7, "noheader": 8,
                                 "nohalf": 6, "nodotz": 6}) if comma else {}
        name = head.strip()
        mat = get_matrix(s, name)
        fmt = str(o.get("format", "%10.0g")).strip()
        list_matrix(s, name, mat, fmt, names=not o.get("nonames"), header=not o.get("noheader"),
                    show_title=str(o.get("title", "")).strip().strip('"'))
        return
    if sub in ("drop", "drop"):
        names = rest.split()
        store = _store(s)
        if names in (["_all"], []):
            store.clear()
            return
        for n in names:
            if n not in store:
                raise StataError(111, f"matrix {n} not found")
            del store[n]
        return
    if sub in ("dir",):
        for n, mat in sorted(_store(s).items()):
            s.output.write(f"{n:>18}[{mat.rows},{mat.cols}]\n", "text")   # VERIFICAR
        return
    if sub in ("rename", "ren"):
        a, b = rest.split()
        store = _store(s)
        store[b] = get_matrix(s, a)
        del store[a]
        return
    if sub in ("rownames", "colnames", "roweq", "coleq", "rown", "coln"):
        _names_cmd(s, {"rown": "rownames", "coln": "colnames"}.get(sub, sub), rest)
        return
    # matrix A = ...
    _define(s, t)


# ---------------------------------------------------------------------------
# return/ereturn matrix
# ---------------------------------------------------------------------------

def take_matrix(s: "Session", text: str) -> tuple[str, Matrix]:
    """return matrix nome = matriz [, copy]: a matriz sai de s.matrices."""
    head, comma, opts = text.partition(",")
    m = re.match(r"\s*([A-Za-z_]\w*)\s*=\s*(\S+)\s*$", head)
    if not m:
        raise StataError(198, "invalid syntax")
    copy = comma and "copy" in opts
    mat = get_matrix(s, m.group(2))
    if not copy and m.group(2) in _store(s):
        del _store(s)[m.group(2)]
    return m.group(1), (mat.copy() if copy else mat)


def ereturn_matrix(s: "Session", sub: str, rest: str) -> None:
    if sub in ("matrix", "mat"):
        name, mat = take_matrix(s, rest)
        s.e[name] = mat
        return
    if sub in ("post", "repost"):
        head, comma, opts = rest.partition(",")
        names = head.split()
        if sub == "post":
            s.e = {}
        if names:
            b = get_matrix(s, names[0])
            s.e["b"] = b.copy()
            _store(s).pop(names[0], None)
        if len(names) > 1:
            V = get_matrix(s, names[1])
            s.e["V"] = V.copy()
            _store(s).pop(names[1], None)
        for m in re.finditer(r"(depname|obs|dof)\(([^)]*)\)", opts if comma else ""):
            key, val = m.group(1), m.group(2).strip()
            if key == "depname":
                s.e["depvar"] = val
            else:
                s.e["N" if key == "obs" else "df_r"] = float(s.eval(val))
        return
    raise StataError(198, "ereturn display arrives with the estimation commands (phase 5)")
