"""matrix: definição, operadores, funções, listagem e uso em expressões."""

import numpy as np


def mat(run, name):
    return run.session.matrices[name]


def test_define_and_ops(run):
    run("matrix A = (1,2\\3,4)\nmatrix B = A'\nmatrix C = A*B\nmatrix D = A + 2*I(2)\n"
        "matrix E = A # I(2)\nmatrix F = (A, A \\ B, B)\nmatrix G = -A / 2")
    assert mat(run, "B").data.tolist() == [[1, 3], [2, 4]]
    assert mat(run, "C").data.tolist() == [[5, 11], [11, 25]]
    assert mat(run, "D").data.tolist() == [[3, 2], [3, 6]]
    assert mat(run, "E").data.shape == (4, 4) and mat(run, "F").data.shape == (4, 4)
    assert mat(run, "G").data.tolist() == [[-0.5, -1], [-1.5, -2]]
    run("matrix X = A + J(3,3,0)")
    assert run.rc == 503


def test_functions(run):
    run("matrix S = (4,2\\2,3)\nmatrix Si = inv(S)\nmatrix L = cholesky(S)\n"
        "matrix d = vecdiag(S)\nmatrix Dg = diag(d)\nmatrix H = hadamard(S, S)\n"
        "matrix v = vec(S)\nmatrix R = corr(S)\nmatrix Z = nullmat(nada) \\ (1,2)")
    assert np.allclose(mat(run, "Si").data @ mat(run, "S").data, np.eye(2))
    assert np.allclose(mat(run, "L").data @ mat(run, "L").data.T, mat(run, "S").data)
    assert mat(run, "d").data.tolist() == [[4, 3]] and mat(run, "Dg").data.tolist() == [[4, 0], [0, 3]]
    assert mat(run, "H").data.tolist() == [[16, 4], [4, 9]] and mat(run, "v").data.shape == (4, 1)
    assert abs(mat(run, "R").data[0, 1] - 2 / np.sqrt(12)) < 1e-12
    assert mat(run, "Z").data.tolist() == [[1, 2]]


def test_scalar_functions_and_elements(run):
    run("matrix A = (1,2,3\\4,5,6)\nmatrix rownames A = a b\nmatrix colnames A = x y z")
    out = run('display rowsof(A) " " colsof(A) " " el(A,2,3) " " A[1,2] " " colnumb(A,"y") '
              '" " rownumb(A,"q") " " el(A,9,9)')
    assert out == "2 3 6 2 2 . .\n"
    run("matrix Q = (2,0\\0,3)")
    assert run("display det(Q) trace(Q) issymmetric(Q) matmissing(Q)") == "6510\n"
    run("matrix Sub = A[1..2, \"y\"...]")
    assert mat(run, "Sub").colnames == ["y", "z"] and mat(run, "Sub").rownames == ["a", "b"]
    run("matrix A[2,1] = 99")
    assert mat(run, "A").data[1, 0] == 99
    loc = run.session.macros.get_local
    run("local r : rownames A\nlocal c : colnames A\nlocal n : rowsof A")
    assert (loc("r"), loc("c"), loc("n")) == ("a b", "x y z", "2")


def test_list_and_manage(run):
    run("matrix A = (1,2\\3,4)")
    assert run("matrix list A") == "\nA[2,2]\n    c1  c2\nr1   1   2\nr2   3   4\n"
    run("matrix S = (1,.5\\.5,1)")
    assert run("matrix list S") == "\nsymmetric S[2,2]\n    c1  c2\nr1   1\nr2  .5   1\n"
    out = run("matrix list A, format(%5.2f) title(teste)")
    assert "A[2,2]:  teste" in out and "1.00" in out
    run("matrix rename A B\nmatrix drop S")
    assert set(run.session.matrices) == {"B"}
    run("matrix input M = (1 2 . \\ 4 5 6)")
    assert mat(run, "M").data.shape == (2, 3)
    run("matrix drop _all")
    assert run.session.matrices == {}
    run("matrix list nada")
    assert run.rc == 111


def test_return_ereturn(run):
    run("program pm, rclass\n matrix T = (1,2)\n return matrix T = T\nend\npm")
    assert run.session.r["T"].data.tolist() == [[1, 2]] and "T" not in run.session.matrices
    run("matrix b = (1,2)\nmatrix V = I(2)\nereturn post b V, depname(y) obs(10)")
    assert run.session.e["N"] == 10 and run.session.e["depvar"] == "y"
    run("display e(b)")
    assert run.rc == 109
    run("matrix bb = e(b)\nmatrix VV = e(V)")
    assert mat(run, "bb").data.tolist() == [[1, 2]]
    out = run("ereturn list")
    assert "e(b) :  1 x 2" in out


def test_top_level_joins(run):
    run("matrix A = (1,2)\nmatrix B = A, A\nmatrix C = A \\ A")
    assert mat(run, "B").data.tolist() == [[1, 2, 1, 2]] and mat(run, "C").data.shape == (2, 2)
    run("forvalues i = 1/3 {\n matrix R = nullmat(R) \\ (`i', `i'^2)\n}")
    assert mat(run, "R").data.tolist() == [[1, 1], [2, 4], [3, 9]]
