"""Block 20d: the program-based forms of nl, nlsur and gmm, and gmm with parameters()."""

import numpy as np

from opendta.session import Session

DATA = """clear
set obs 400
set seed 78
gen x1 = round(invnormal(runiform()), .01)
gen x2 = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen xp = round(runiform()*4 + 1, .01)
gen yp = round(2*xp^.7 + invnormal(runiform())*.2, .01)
gen y = round(exp(.3 + .5*x1 - .4*x2) + invnormal(runiform())*.5, .01)
gen y1 = round(1 + 2*x1 + invnormal(runiform())*.5, .01)
gen y2 = round(-1 + .5*x1 + 3*x2 + invnormal(runiform())*.5, .01)
program nlpow
    syntax varlist(min=2 max=2) if, at(name)
    local y : word 1 of `varlist'
    local x : word 2 of `varlist'
    replace `y' = `at'[1,1]*`x'^`at'[1,2] `if'
end
program nlpow2, rclass
    syntax varlist(min=2 max=2) [if]
    local y : word 1 of `varlist'
    local x : word 2 of `varlist'
    return local eq "`y' = {a}*`x'^{b=1}"
    return local title "Power curve"
end
program gmm_poi
    syntax varlist if, at(name)
    replace `varlist' = y - exp(`at'[1,1]*x1 + `at'[1,2]*x2 + `at'[1,3]) `if'
end
program nlsurlin
    syntax varlist(min=3 max=3) if, at(name)
    tokenize `varlist'
    replace `1' = `at'[1,1] + `at'[1,2]*`3' `if'
    replace `2' = `at'[1,3] + `at'[1,2]*`3' `if'
end
program nlsurlin2, rclass
    syntax varlist(min=3 max=3) [if]
    tokenize `varlist'
    return scalar n_eq = 2
    return local eq_1 "`1' = {a1} + {b}*`3'"
    return local eq_2 "`2' = {a2} + {b}*`3'"
end
"""


def _run(text, s=None):
    s = s or Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def _bV(s):
    e = s._estimates
    return np.asarray(e.b, dtype=float), np.asarray(e.V, dtype=float)


def test_nl_function_evaluator_matches_expression():
    s, _ = _run(DATA + "nl (yp = {a}*xp^{b=1})\n")
    b0, V0 = _bV(s)
    _, out = _run("nl pow @ yp xp, parameters(a b) initial(b 1)\n", s)
    b1, V1 = _bV(s)
    assert np.allclose(b0, b1, atol=1e-6) and np.allclose(V0, V1, rtol=1e-4)
    assert list(s._estimates.names) == ["a", "b"]
    _, out = _run('display "`e(funcprog)\'"\n', s)
    assert "nlpow" in out


def test_nl_nparameters_and_sexp_program():
    s, _ = _run(DATA + "nl pow @ yp xp, nparameters(2) initial(b2 1)\n")
    assert list(s._estimates.names) == ["b1", "b2"]
    b1, _ = _bV(s)
    _, out = _run("nl pow2 : yp xp\n", s)
    assert np.allclose(_bV(s)[0], b1, atol=1e-6)
    assert "Power curve" in out and '"Power curve"' not in out


def test_gmm_moment_program_and_parameters():
    s, _ = _run(DATA + "gmm (y - exp({b1}*x1 + {b2}*x2 + {b0})), instruments(x1 x2 z)\n")
    b0, V0 = _bV(s)
    _run("gmm gmm_poi, nequations(1) instruments(x1 x2 z) parameters(b1 b2 b0)\n", s)
    assert np.allclose(_bV(s)[0], b0, atol=1e-6)
    assert np.allclose(_bV(s)[1], V0, rtol=1e-4)
    _run("gmm gmm_poi, equations(m1) instruments(m1: x1 x2 z) nparameters(3)\n", s)
    assert list(s._estimates.names) == ["b1", "b2", "b3"]
    _run("gmm (y - exp({b1}*x1 + {b2}*x2 + {b0})), instruments(x1 x2 z) parameters(b0 b1 b2)\n", s)
    assert list(s._estimates.names) == ["b0", "b1", "b2"]
    assert np.allclose(_bV(s)[0], b0[[2, 0, 1]], atol=1e-6)


def test_gmm_parameters_must_match():
    s = Session()
    _run(DATA, s)
    rc = s.run_text("gmm (y - exp({b1}*x1 + {b0})), instruments(x1 z) parameters(b0 c)\n")
    assert rc == 198


def test_nlsur_programs():
    s, _ = _run(DATA + "nlsur (y1 = {a1} + {b}*x1) (y2 = {a2} + {b}*x1), nls\n")
    b0, _ = _bV(s)
    _run("nlsur lin @ y1 y2 x1, nequations(2) parameters(a1 b a2) nls\n", s)
    assert np.allclose(_bV(s)[0], b0, atol=1e-6)
    _run("nlsur lin2 : y1 y2 x1, nls\n", s)
    assert np.allclose(_bV(s)[0], b0, atol=1e-6)
    _run("nlsur lin @ y1 y2 x1, nequations(2) parameters(a1 b a2)\n", s)
    assert s._estimates.extra["deps"] == ["y1", "y2"]


def test_unknown_program():
    s = Session()
    _run(DATA, s)
    assert s.run_text("nl nosuch @ yp xp, nparameters(2)\n") == 199
