"""The vectorized forms of the functions (lang/vexpr.py) give the same values as the
observation-by-observation ones (lang/functions.py), missing values included."""

import numpy as np
import pytest

from opendta.core import missing as M
from opendta.lang import functions as F
from opendta.lang import vexpr as V

SYS = M.SYSMISS
rng = np.random.default_rng(7)
X = np.concatenate([rng.normal(0, 3, 40), [0.0, 1.0, -1.0, 0.5, 2.0, SYS, 1e300]])
P = np.concatenate([rng.uniform(0, 1, 40), [0.0, 1.0, -0.2, 1.2, 0.5, SYS, 0.999999]])
DF = np.concatenate([rng.integers(1, 30, 40).astype(float), [0.0, -1.0, 1.0, 2.5, 3.0, SYS, 5.0]])
K = np.concatenate([rng.integers(0, 12, 40).astype(float), [0.0, 1.0, 2.0, 3.0, 4.0, SYS, 0.0]])

CASES = {
    "sin": (X,), "cos": (X,), "tan": (X,), "asin": (P,), "acos": (P,), "atan": (X,),
    "atan2": (X, P), "sinh": (P,), "cosh": (P,), "tanh": (X,), "sign": (X,),
    "lngamma": (X,), "lnfactorial": (K,), "logit": (P,), "invlogit": (X,), "digamma": (DF,),
    "normal": (X,), "normalden": (X,), "invnormal": (P,), "lnnormal": (X,),
    "chi2": (DF, np.abs(X)), "chi2tail": (DF, np.abs(X)), "invchi2": (DF, P), "invchi2tail": (DF, P),
    "F": (DF, DF + 1, np.abs(X)), "Ftail": (DF, DF + 1, np.abs(X)), "invF": (DF, DF + 2, P),
    "invFtail": (DF, DF + 2, P), "t": (DF, X), "ttail": (DF, X), "invt": (DF, P), "invttail": (DF, P),
    "binomial": (K + 5, K, P), "binomialtail": (K + 5, K, P), "poisson": (np.abs(X) + 0.1, K),
    "poissontail": (np.abs(X) + 0.1, K), "betaden": (DF, DF + 1, P), "ibeta": (DF, DF + 1, P),
    "gammap": (DF, np.abs(X)),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_vector_equals_scalar(name):
    args = CASES[name]
    vec = V._FAST[name](*args)
    scal = np.array([F.call(name, [float(a[i]) for a in args]) for i in range(len(args[0]))])
    both_miss = (vec >= SYS) & (scal >= SYS)
    close = np.isclose(vec, scal, rtol=1e-13, atol=1e-300)
    assert np.all(both_miss | close), (name, vec[~(both_miss | close)], scal[~(both_miss | close)])


def test_elementwise_fallback_caches_repeated_arguments(monkeypatch):
    calls = []
    real = F.call

    def spy(name, args):
        calls.append(args)
        return real(name, args)

    monkeypatch.setattr(F, "call", spy)
    col = np.array([float(i % 5) for i in range(1000)])
    out = V._elementwise("string", [col], 1000)
    assert out[:6] == ["0", "1", "2", "3", "4", "0"]
    assert len(calls) == 5
