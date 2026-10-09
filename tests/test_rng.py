"""Números aleatórios: MT19937-64, set seed, c(seed) e funções r*()."""

import numpy as np
import pytest

from opendta.core.rng import MT64

MISS = 8.98846567431158e307


def col(run, name):
    return np.array(run.session.data.get(name).data)


def test_mt64_reference_value():
    # valor exigido pelo padrão C++ para std::mt19937_64 com a semente 5489
    g = MT64(5489)
    assert int(g.uint64(10000)[-1]) == 9981545732273789042
    # gerar em pedaços dá a mesma sequência
    g1, g2 = MT64(42), MT64(42)
    a = g1.uint64(1000)
    b = np.concatenate([g2.uint64(7), g2.uint64(500), g2.uint64(493)])
    assert np.array_equal(a, b)


def test_seed_reproducible(run):
    run("set seed 2024\nset obs 1000\ngen a = runiform()\nset seed 2024\ngen b = runiform()")
    assert np.array_equal(col(run, "a"), col(run, "b"))
    run("gen c = runiform()")
    assert not np.array_equal(col(run, "a"), col(run, "c"))
    a = col(run, "a")
    assert (a > 0).all() and (a < 1).all()
    assert run.session.data.get("a").vtype == "float"


def test_scalar_and_vector_share_stream(run):
    run("set seed 7")
    first = run.eval("runiform()")
    run("set seed 7\nset obs 2\ngen double u = runiform()")
    assert col(run, "u")[0] == first


def test_seed_state_roundtrip(run):
    run("set seed 99\nset obs 50\nlocal st = c(seed)\ngen a = runiform()\nset seed `st'\ngen b = runiform()")
    assert run.rc == 0
    assert np.array_equal(col(run, "a"), col(run, "b"))
    assert run.eval("c(seed)").startswith("X")


def test_seed_errors(run):
    run("set seed -1")                  # o Stata 14 aceita
    assert run.rc == 0
    run("set seed abc")
    assert run.rc == 198
    run("set seed X123")
    assert run.rc == 198
    run("set rng kiss32")
    assert run.rc == 198
    run("set rng mt64")
    assert run.rc == 0


def test_distribution_moments(run):
    run("""set seed 1
set obs 200000
gen double n = rnormal(5, 2)
gen double e = rexponential(3)
gen double p = rpoisson(4)
gen double b = rbinomial(10, .3)
gen double c = rchi2(3)
gen double g = rgamma(2, 3)
gen double be = rbeta(2, 5)
gen double t = rt(10)
gen double k = runiformint(1, 6)
gen double ab = runiform(-1, 3)
gen double lg = rlogistic()
""")
    assert run.rc == 0
    tol = 0.03
    assert col(run, "n").mean() == pytest.approx(5, abs=tol)
    assert col(run, "n").std() == pytest.approx(2, abs=tol)
    assert col(run, "e").mean() == pytest.approx(3, rel=tol)
    assert col(run, "p").mean() == pytest.approx(4, rel=tol)
    assert col(run, "b").mean() == pytest.approx(3, rel=tol)
    assert col(run, "c").mean() == pytest.approx(3, rel=tol)
    assert col(run, "g").mean() == pytest.approx(6, rel=tol)
    assert col(run, "be").mean() == pytest.approx(2 / 7, rel=tol)
    assert col(run, "t").var() == pytest.approx(10 / 8, rel=0.05)
    k = col(run, "k")
    assert set(np.unique(k)) == {1, 2, 3, 4, 5, 6}
    ab = col(run, "ab")
    assert ab.min() >= -1 and ab.max() < 3
    assert col(run, "lg").var() == pytest.approx(np.pi ** 2 / 3, rel=0.05)
    p = col(run, "p")
    assert (p == np.trunc(p)).all()


def test_invalid_parameters_give_missing(run):
    run("set obs 3\ngen a = rnormal(0, -1)\ngen b = runiform(2, 1)\ngen c = rpoisson(-1)\n"
        "gen d = runiformint(1.5, 3)\ngen e = rbinomial(10, 2)\ngen f = rnormal(., 1)")
    for v in "abcdef":
        assert (col(run, v) >= MISS).all(), v
    run('gen z = runiform("a")')
    assert run.rc == 109
    run("gen z = rpoisson()")
    assert run.rc == 198


def test_vector_parameters(run):
    run("set obs 4\ngen m = _n * 100\ngen x = rnormal(m, 0)")
    assert list(col(run, "x")) == [100, 200, 300, 400]


def test_runiform_matches_stata14():
    # valores de compat/expected/0307_aleatorios.log (Stata/SE 14.0)
    cases = {123: [0.31320017867847072, 0.55597911939485856], 0: [0.15979336337046079],
             2147483647: [0.16803268731662413]}
    for seed, want in cases.items():
        got = MT64(seed).uniform(len(want))
        assert [float(f"{x:.17f}") for x in got] == want


def test_draws_only_in_sample(run):
    # como no Stata, o if/in decide quais observações consomem a sequência
    run("set obs 4\nset seed 1\ngen double a = runiform() if mod(_n, 2)\n"
        "set seed 1\ngen double b = runiform() in 1/2")
    a, b = col(run, "a"), col(run, "b")
    assert a[0] == b[0] and a[2] == b[1]
    assert a[1] >= MISS and b[3] >= MISS
