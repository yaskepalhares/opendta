"""Mata classes and moptimize()."""

import numpy as np

from opendta.session import Session


def _run(text: str):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(text)
    return s, "".join(lines), rc


CLASSES = """mata:
class point {
    real scalar x, y
    static real scalar count
    void set()
    real scalar norm()
    virtual string scalar describe()
}
void point::new() {
    count = (count == . ? 1 : count + 1)
}
void point::set(real scalar a, real scalar b) {
    x = a
    this.y = b
}
real scalar point::norm() {
    return(sqrt(x^2 + y^2))
}
string scalar point::describe() {
    return("point " + strofreal(norm()))
}
class point3 extends point {
    real scalar z
    virtual string scalar describe()
    real scalar norm()
}
real scalar point3::norm() {
    return(sqrt(super.norm()^2 + z^2))
}
string scalar point3::describe() {
    return("point3 " + strofreal(norm()))
}
string scalar test() {
    class point scalar p
    class point3 scalar q
    p.set(3, 4)
    q.set(1, 2)
    q.z = 2
    return(p.describe() + "|" + q.describe() + "|" + strofreal(p.count))
}
st_local("r", test())
p = point()
p.set(6, 8)
st_numscalar("n", p.norm())
end
"""


def test_classes_methods_inheritance_statics():
    s, out, rc = _run(CLASSES + 'display "`r\'"\ndisplay n')
    assert rc == 0, out
    assert "point 5|point3 3|2" in out
    assert "\n10\n" in out


def test_class_undefined_method():
    s, out, rc = _run(CLASSES + "mata: p.nothing()")
    assert rc != 0


MOPT = """clear
set seed 408
set obs 200
gen x1 = rnormal()
gen x2 = runiform()
gen y = runiform() < invlogit(-.5 + x1 - 2*x2)
quietly logit y x1 x2
matrix bl = e(b)
matrix Vl = e(V)
mata:
void lflogit(transmorphic M, real rowvector b, real colvector lnfj) {
    real colvector xb, y
    xb = moptimize_util_xb(M, b, 1)
    y = moptimize_util_depvar(M, 1)
    lnfj = y:*ln(invlogit(xb)) + (1:-y):*ln(invlogit(-xb))
}
void d2logit(transmorphic M, real scalar todo, real rowvector b, lnf, g, H) {
    real colvector xb, y, p
    xb = moptimize_util_xb(M, b, 1)
    y = moptimize_util_depvar(M, 1)
    p = invlogit(xb)
    lnf = moptimize_util_sum(M, y:*ln(p) + (1:-y):*ln(1:-p))
    if (todo >= 1) {
        g = moptimize_util_vecsum(M, 1, y - p, lnf)
        if (todo == 2) H = -moptimize_util_matsum(M, 1, 1, p:*(1:-p), lnf)
    }
}
M = moptimize_init()
moptimize_init_evaluator(M, &lflogit())
moptimize_init_evaluatortype(M, "lf")
moptimize_init_depvar(M, 1, "y")
moptimize_init_eq_indepvars(M, 1, "x1 x2")
moptimize_init_tracelevel(M, "none")
moptimize(M)
st_matrix("b1", moptimize_result_coefs(M))
M2 = moptimize_init()
moptimize_init_evaluator(M2, &d2logit())
moptimize_init_evaluatortype(M2, "d2")
moptimize_init_depvar(M2, 1, "y")
moptimize_init_eq_indepvars(M2, 1, "x1 x2")
moptimize_init_tracelevel(M2, "none")
moptimize(M2)
st_matrix("b2", moptimize_result_coefs(M2))
st_matrix("V2", moptimize_result_V(M2))
moptimize_result_post(M2)
end
"""


def test_moptimize_matches_logit():
    s, out, rc = _run(MOPT)
    assert rc == 0, out[-2000:]
    m = s.matrices if hasattr(s, "matrices") else None
    from opendta.commands.matrix import get_matrix
    bl = get_matrix(s, "bl").data
    b1 = get_matrix(s, "b1").data
    b2 = get_matrix(s, "b2").data
    assert np.allclose(b1, bl, atol=1e-5)
    assert np.allclose(b2, bl, atol=1e-7)
    assert np.allclose(get_matrix(s, "V2").data, get_matrix(s, "Vl").data, rtol=1e-6)
    assert s.e["N"] == 200 and abs(s.e["ll"] - (-104.92824)) < 1e-4
    del m
