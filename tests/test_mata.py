"""Mata: analisador, operadores, funções do usuário, biblioteca e st_*."""

import numpy as np
import pytest

from opendta.mata.lexer import MataSyntaxError
from opendta.mata.parser import parse

MISS = 8.98846567431158e307


def mata(run, code: str) -> str:
    """Roda um bloco mata ... end e devolve a saída (sem eco)."""
    return run("mata\n" + code + "\nend")


def val(run, name):
    return run.session.mata.globals.vars[name]


def test_parser_precedence():
    (st,) = parse("x = 1 + 2 * 3 ^ 2")
    assert st.expr.value.op == "+" and st.expr.value.right.op == "*"
    (st,) = parse("(1,2\\3,4)")
    assert st.expr.op == "\\" and st.expr.left.op == ","
    (st,) = parse("-2^2")
    assert st.expr.op == "-" and st.expr.operand.op == "^"
    (st,) = parse("A[1, .]")
    assert st.expr.cols is None
    with pytest.raises(MataSyntaxError) as e:
        parse("x = (1,")
    assert e.value.incomplete
    with pytest.raises(MataSyntaxError) as e:
        parse("real scalar f(real x) {")
    assert e.value.incomplete


def test_scalar_and_matrix_display(run):
    out = mata(run, "x = 2\nx\nA = (1,2\\3,4)\nA\n\"hi\"")
    assert "  2\n" in out
    assert ("       1   2\n    +---------+\n  1 |  1   2  |\n  2 |  3   4  |\n"
            "    +---------+\n") in out
    assert "  hi\n" in out
    out = mata(run, "(2,1\\1,2)")
    assert out.startswith("[symmetric]\n")
    assert "  1 |  2      |" in out


def test_arithmetic_and_missing(run):
    mata(run, "a = (1,2\\3,4)\nb = a * a\nc = a :* a\nd = a + 1\ne = (1, .) :+ 1\n"
              "f = a'\ng = 2^3\nh = (1,2) # (1\\1)\ni = 1/0")
    assert val(run, "b").a.tolist() == [[7, 10], [15, 22]]
    assert val(run, "c").a.tolist() == [[1, 4], [9, 16]]
    assert val(run, "d").a.tolist() == [[2, 3], [4, 5]]
    assert val(run, "e").a[0, 1] >= MISS
    assert val(run, "f").a.tolist() == [[1, 3], [2, 4]]
    assert val(run, "g").a[0, 0] == 8
    assert val(run, "h").a.tolist() == [[1, 2], [1, 2]]
    assert val(run, "i").a[0, 0] >= MISS


def test_ranges_joins_subscripts(run):
    mata(run, "r = 1..4\nc = 4::2\nA = (1,2,3\\4,5,6)\nx = A[2,3]\ny = A[.,2]\n"
              "z = A[(1\\2), (3,1)]\nw = r[(2,4)]\nq = A[|1,2 \\ 2,3|]\nA[1,1] = 9\nA[2,.] = (7,8,9)")
    assert val(run, "r").a.tolist() == [[1, 2, 3, 4]]
    assert val(run, "c").a.tolist() == [[4], [3], [2]]
    assert val(run, "x").a[0, 0] == 6
    assert val(run, "y").a.tolist() == [[2], [5]]
    assert val(run, "z").a.tolist() == [[3, 1], [6, 4]]
    assert val(run, "w").a.tolist() == [[2, 4]]
    assert val(run, "q").a.tolist() == [[2, 3], [5, 6]]
    assert val(run, "A").a.tolist() == [[9, 2, 3], [7, 8, 9]]


def test_relational_and_logic(run):
    mata(run, "a = (1,2) == (1,2)\nb = (1,2) :== (1,3)\nc = 3 > 2 & 1 < 0\nd = (2 > 1) ? \"sim\" : \"não\"")
    assert val(run, "a").a[0, 0] == 1
    assert val(run, "b").a.tolist() == [[1, 0]]
    assert val(run, "c").a[0, 0] == 0
    assert val(run, "d").a[0, 0] == "sim"


def test_strings(run):
    mata(run, 's = "ab" + "cd"\nt = "x" * 3\nu = strupper(("a","b"))\nv = substr("hello", 2, 3)\n'
              'w = tokens("um dois  tres")\nn = strlen("abc")\nz = strtoreal("1.5")\n'
              'q = invtokens(("a","b"), "-")')
    assert val(run, "s").a[0, 0] == "abcd"
    assert val(run, "t").a[0, 0] == "xxx"
    assert val(run, "u").a.tolist() == [["A", "B"]]
    assert val(run, "v").a[0, 0] == "ell"
    assert val(run, "w").a.tolist() == [["um", "dois", "tres"]]
    assert val(run, "n").a[0, 0] == 3 and val(run, "z").a[0, 0] == 1.5
    assert val(run, "q").a[0, 0] == "a-b"


def test_control_flow_and_display_in_loops(run):
    out = mata(run, "for (i=1; i<=3; i++) i\nk = 0\nwhile (k < 5) k = k + 2\n"
                    "do { k-- } while (k > 0)\nif (k == 0) \"zero\"\nelse \"outro\"")
    assert "  1\n  2\n  3\n" in out
    assert val(run, "k").a[0, 0] == 0
    assert "  zero\n" in out and "outro" not in out


def test_user_functions_and_references(run):
    out = mata(run, """real scalar dobro(real scalar x)
{
    return(2*x)
}
void incrementa(real scalar x)
{
    x = x + 1
}
real matrix fat(real scalar n)
{
    if (n <= 1) return(1)
    return(n * fat(n-1))
}
function opcional(a, | b)
{
    if (args() == 1) return(a)
    return(a + b)
}
y = 5
incrementa(y)
dobro(y)
fat(5)""")
    assert "  12\n" in out and "  120\n" in out
    assert val(run, "y").a[0, 0] == 6


def test_errors(run):
    out = mata(run, "nada\n(1,2) * (3,4)\nx = 1")
    assert "                 <istmt>:  3499  nada not found\nr(3499);" in out
    assert "                       *:  3200  conformability error" in out
    assert val(run, "x").a[0, 0] == 1             # o bloco continua
    out = run("mata: nada2")
    assert run.rc == 3499
    out = run("mata:\nnada3\nx = 7\nend")
    assert run.rc == 3499                          # mata: para no erro
    out = mata(run, "real scalar f() {\n    return((1,2)*(3,4))\n}\nf()")
    # o operador aparece como a "função" que falhou, depois f() e <istmt>
    assert ("                       *:  3200  conformability error\n"
            "                     f():     -  function returned error\n"
            "                 <istmt>:     -  function returned error\n") in out


def test_library_numeric(run):
    mata(run, """m = mean((1,2\\3,4\\5,.))
v = variance((1,2\\3,5\\5,6))
s = sum((1,.,3))
mx = max((1,7,.))
rs = rowsum((1,2\\3,4))
i = invsym((4,2\\2,3))
d = det((1,2\\3,4))
c = cholesky((4,2\\2,3))
x = lusolve((2,0\\0,4), (2\\8))
o = sort((3,1\\1,2\\2,9), 1)
u = uniqrows((2\\1\\2))
sel = select((1,2,3)', (1\\0\\1))
L = .
V = .
symeigensystem((2,1\\1,2), V, L)
J1 = J(2, 3, 7)
id = I(2)
miss = missing((1,.,.a))""")
    assert val(run, "m").a.tolist() == [[2, 3]]            # linha com missing sai
    X = np.array([[1, 2], [3, 5], [5, 6]], dtype=float)
    assert np.allclose(val(run, "v").a, np.cov(X.T))
    assert val(run, "s").a[0, 0] == 4 and val(run, "mx").a[0, 0] == 7
    assert val(run, "rs").a.tolist() == [[3], [7]]
    assert np.allclose(val(run, "i").a, np.linalg.inv([[4, 2], [2, 3]]))
    assert val(run, "d").a[0, 0] == pytest.approx(-2)
    assert np.allclose(val(run, "c").a, np.linalg.cholesky([[4, 2], [2, 3]]))
    assert val(run, "x").a.tolist() == [[1], [2]]
    assert val(run, "o").a.tolist() == [[1, 2], [2, 9], [3, 1]]
    assert val(run, "u").a.tolist() == [[1], [2]]
    assert val(run, "sel").a.tolist() == [[1], [3]]
    assert np.allclose(val(run, "L").a, [[3, 1]])
    assert val(run, "J1").a.tolist() == [[7, 7, 7], [7, 7, 7]]
    assert val(run, "miss").a[0, 0] == 2


def test_invsym_singular_zeroes_dependent_columns(run):
    mata(run, "i = invsym((1,1,0\\1,1,0\\0,0,2))")
    i = val(run, "i").a
    assert i[1].tolist() == [0, 0, 0] and i[:, 1].tolist() == [0, 0, 0]
    assert i[0, 0] == pytest.approx(1) and i[2, 2] == pytest.approx(0.5)


def test_printf_and_sprintf(run):
    out = mata(run, 'printf("%5.2f|%s|%g\\n", 3.14159, "ab", 2)\n'
                    'printf("{txt}texto {res}%9.0g\\n", 10)\nz = sprintf("%04.1f", 2.5)')
    assert " 3.14|ab|2\n" in out
    assert "texto        10\n" in out


def test_structs_and_pointers(run):
    out = mata(run, """struct ponto {
    real scalar x, y
}
real scalar norma(struct ponto scalar p)
{
    return(sqrt(p.x^2 + p.y^2))
}
real scalar teste()
{
    struct ponto scalar p
    p.x = 3
    p.y = 4
    return(norma(p))
}
teste()
a = 5
p = &a
*p = 9
a""")
    assert "  5\n" in out
    assert val(run, "a").a[0, 0] == 9


def test_stata_interface(run):
    run("set obs 3\ngen x = _n\ngen str2 s = \"a\" + string(_n)")
    mata(run, """X = st_data(., "x")
S = st_sdata(., "s")
st_store(., "x", X :* 10)
k = st_addvar("double", "y")
st_store(., "y", X :^ 2)
st_local("loc", "valor")
st_global("glob", "g")
st_numscalar("sc", 42)
st_matrix("M", (1,2\\3,4))
n = st_nobs()
nv = st_nvar()
nome = st_varname(2)
idx = st_varindex("y")
st_numscalar("r(teste)", 7)""")
    ds = run.session.data
    assert list(ds.get("x").data) == [10, 20, 30]
    assert list(ds.get("y").data) == [1, 4, 9]
    assert val(run, "S").a.ravel().tolist() == ["a1", "a2", "a3"]
    assert val(run, "k").a[0, 0] == 3 and val(run, "n").a[0, 0] == 3
    assert val(run, "nome").a[0, 0] == "s" and val(run, "idx").a[0, 0] == 3
    assert run.session.macros.get_local("loc") == "valor"
    assert run.session.macros.get_global("glob") == "g"
    assert run.session.scalars["sc"] == 42
    assert run.session.matrices["M"].data.tolist() == [[1, 2], [3, 4]]
    assert run.session.r["teste"] == 7
    out = run("mata: st_matrix(\"M\")")
    assert "  2 |  3   4  |" in out


def test_stata_from_mata_and_mata_subcommands(run):
    out = run('mata: stata("display 1+1")')
    assert "2\n" in out
    run("mata: x = 1")
    out = run("mata: x")
    assert out == "  1\n"
    run("mata clear")
    run("mata: x")
    assert run.rc == 3499


def test_block_echo_in_do_file(run, tmp_path):
    do = tmp_path / "m.do"
    do.write_text("mata\nx = 1\nx\nreal scalar f(real scalar a)\n{\n    return(a)\n}\nend\n")
    out = run(f'do "{do}"')
    assert "------------------------------------------------- mata (type end to exit) ------\n" in out
    assert ": x = 1\n\n: x\n  1\n\n" in out
    assert ": real scalar f(real scalar a)\n> {\n>     return(a)\n> }\n\n: end\n" in out
