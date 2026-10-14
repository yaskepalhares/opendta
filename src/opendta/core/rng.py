"""Random-number generator: 64-bit Mersenne Twister (MT19937-64).

The Stata 14 manual ([R] set seed) says the default generator is the 64-bit
Mersenne Twister (mt64). Here it follows the algorithm published by
Matsumoto and Nishimura (init_genrand64 + genrand64_int64). Undocumented
details are marked VERIFY (case compat/do/0307_aleatorios.do):

* `set seed #` becomes the initial state through init_genrand64(#) and the
  64-bit integer becomes a double in (0,1) from the top 52 bits + 1/2 (like
  genrand64_real3). With this, runiform() reproduces Stata 14's sequence
  exactly (checked in compat/expected/0307_aleatorios.log with the
  seeds 0, 42, 123 and 2147483647).
* rbinomial() (inversion with one uniform), rexponential(), rchi2(2) and
  rgamma(1, b) (-b·ln u) also reproduce Stata (compat 0309).
* rnormal() consumes three uniforms per draw in Stata (observed in
  compat 0309), with an undocumented algorithm; here it uses inversion, so the
  sequence differs. The same holds for runiformint(), rpoisson(), rt() and
  the general forms of rgamma/rchi2 (VERIFY).
"""

from __future__ import annotations

import numpy as np

from . import missing as M

NN = 312
MM = 156
_MATRIX_A = np.uint64(0xB5026F5AA96619E9)
_UM = np.uint64(0xFFFFFFFF80000000)   # 33 bits altos
_LM = np.uint64(0x7FFFFFFF)           # 31 bits baixos
_MASK64 = (1 << 64) - 1


class MT64:
    def __init__(self, seed: int = 123456789) -> None:
        self.mt = np.zeros(NN, dtype=np.uint64)
        self.mti = NN + 1
        self.seed(seed)

    # -- estado ----------------------------------------------------------------
    def seed(self, seed: int) -> None:
        mt = [0] * NN
        mt[0] = seed & _MASK64
        for i in range(1, NN):
            prev = mt[i - 1]
            mt[i] = (6364136223846793005 * (prev ^ (prev >> 62)) + i) & _MASK64
        self.mt = np.array(mt, dtype=np.uint64)
        self.mti = NN

    def state(self) -> str:
        """State as text (c(seed)); `set seed` accepts it back."""
        return "X" + format(self.mti, "03x") + "".join(format(int(x), "016x") for x in self.mt)

    def set_state(self, text: str) -> None:
        body = text[1:]
        if len(body) != 3 + 16 * NN:
            raise ValueError("invalid state")
        self.mti = int(body[:3], 16)
        self.mt = np.array([int(body[3 + 16 * i: 19 + 16 * i], 16) for i in range(NN)], dtype=np.uint64)

    # -- generation -------------------------------------------------------------
    def _twist(self) -> None:
        mt = self.mt
        one = np.uint64(1)

        def mix(upper: np.ndarray, lower: np.ndarray, far: np.ndarray) -> np.ndarray:
            x = (upper & _UM) | (lower & _LM)
            mag = np.where((x & one) == one, _MATRIX_A, np.uint64(0))
            return far ^ (x >> one) ^ mag

        # i = 0..155: usa valores antigos
        a = mix(mt[0:NN - MM], mt[1:NN - MM + 1], mt[MM:NN])
        mt[0:NN - MM] = a
        # i = 156..310: the "distant" term is already new (mt[i-156])
        mt[NN - MM:NN - 1] = mix(mt[NN - MM:NN - 1], mt[NN - MM + 1:NN], mt[0:MM - 1])
        # i = 311
        mt[NN - 1] = mix(mt[NN - 1:NN], mt[0:1], mt[MM - 1:MM])[0]
        self.mti = 0

    def uint64(self, n: int) -> np.ndarray:
        out = np.empty(n, dtype=np.uint64)
        k = 0
        while k < n:
            if self.mti >= NN:
                self._twist()
            take = min(n - k, NN - self.mti)
            y = self.mt[self.mti:self.mti + take].copy()
            y ^= (y >> np.uint64(29)) & np.uint64(0x5555555555555555)
            y ^= (y << np.uint64(17)) & np.uint64(0x71D67FFFEDA60000)
            y ^= (y << np.uint64(37)) & np.uint64(0xFFF7EEE000000000)
            y ^= y >> np.uint64(43)
            out[k:k + take] = y
            k += take
            self.mti += take
        return out

    def uniform(self, n: int) -> np.ndarray:
        """Doubles in (0,1). VERIFY the conversion used by Stata."""
        x = self.uint64(n)
        return ((x >> np.uint64(12)).astype(np.float64) + 0.5) * (1.0 / 4503599627370496.0)


RNG = MT64()


def reset(seed: int = 123456789) -> None:
    RNG.seed(seed)


# ---------------------------------------------------------------------------
# Distributions: (min, max number of arguments, function(u, *args) -> values)
# ---------------------------------------------------------------------------

def _bad(*conds) -> np.ndarray:
    out = np.zeros(np.broadcast(*conds).shape if conds else (), dtype=bool)
    for c in conds:
        out = out | c
    return out


def _clean(v: np.ndarray, bad: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    return np.where(bad | ~np.isfinite(v), M.SYSMISS, v)


def _uniform(u, a=0.0, b=1.0):
    return _clean(a + (b - a) * u, (a >= M.SYSMISS) | (b >= M.SYSMISS) | (a > b))


def _uniformint(u, a, b):
    a2, b2 = np.trunc(a), np.trunc(b)
    bad = (a >= M.SYSMISS) | (b >= M.SYSMISS) | (a2 > b2) | (a != a2) | (b != b2)
    return _clean(a2 + np.floor(u * (b2 - a2 + 1)), bad)


def _normal(u, m=0.0, s=1.0):
    from scipy.special import ndtri
    bad = (m >= M.SYSMISS) | (s >= M.SYSMISS) | (s < 0)
    return _clean(m + s * ndtri(u), bad)


def _ppf(dist, *params_ok):
    def f(u, *args):
        from scipy import stats as st
        bad = _bad(*[np.asarray(a) >= M.SYSMISS for a in args]) if args else np.zeros(np.shape(u), dtype=bool)
        for check in params_ok:
            bad = bad | ~check(*args)
        safe = [np.where(bad, 1.0, a) for a in args]
        with np.errstate(all="ignore"):
            v = getattr(st, dist).ppf(u, *safe)
        return _clean(v, bad)
    return f


def _gamma(u, a, b=1.0):
    from scipy import stats as st
    bad = (np.asarray(a) >= M.SYSMISS) | (np.asarray(b) >= M.SYSMISS) | (a <= 0) | (b <= 0)
    with np.errstate(all="ignore"):
        v = st.gamma.ppf(u, np.where(bad, 1, a), scale=np.where(bad, 1, b))
        # shape 1 = exponential: -b·ln(u), as in Stata 14 (rchi2(2), compat 0309); VERIFY the other
        # shapes
        v = np.where(np.asarray(a) == 1, -np.asarray(b) * np.log(u), v)
    return _clean(v, bad)


def _chi2(u, df):
    """rchi2(df) = 2·rgamma(df/2): with df = 2 gives -2·ln(u), as in Stata 14."""
    return _gamma(u, np.asarray(df, dtype=np.float64) / 2, 2.0)


def _exponential(u, b):
    # -b·ln(u): reproduces Stata 14's rexponential() (compat 0309)
    bad = (np.asarray(b) >= M.SYSMISS) | (b <= 0)
    return _clean(-b * np.log(u), bad)


def _logistic(u, m=0.0, s=1.0):
    bad = (np.asarray(m) >= M.SYSMISS) | (np.asarray(s) >= M.SYSMISS) | (s < 0)
    return _clean(m + s * np.log(u / (1 - u)), bad)


def _weibull(u, a, b, g=0.0):
    bad = (a <= 0) | (b <= 0) | (np.asarray(a) >= M.SYSMISS) | (np.asarray(b) >= M.SYSMISS)
    with np.errstate(all="ignore"):
        v = g + b * (-np.log1p(-u)) ** (1 / np.where(bad, 1, a))
    return _clean(v, bad)


def _nbinomial(u, n, p):
    from scipy import stats as st
    bad = (n <= 0) | (p <= 0) | (p > 1) | (np.asarray(n) >= M.SYSMISS) | (np.asarray(p) >= M.SYSMISS)
    with np.errstate(all="ignore"):
        v = st.nbinom.ppf(u, np.where(bad, 1, n), np.where(bad, 0.5, p))
    return _clean(v, bad)


def _hyper(u, N, K, n):
    from scipy import stats as st
    bad = (K > N) | (n > N) | (N <= 0) | (K < 0) | (n < 0)
    with np.errstate(all="ignore"):
        v = st.hypergeom.ppf(u, np.where(bad, 1, N), np.where(bad, 0, K), np.where(bad, 0, n))
    return _clean(v, bad)


DISTRIBUTIONS = {
    "runiform": (0, 2, _uniform),
    "runiformint": (2, 2, _uniformint),
    "rnormal": (0, 2, _normal),
    "rbinomial": (2, 2, _ppf("binom", lambda n, p: (n >= 1) & (n == np.trunc(n)) & (p >= 0) & (p <= 1))),
    "rpoisson": (1, 1, _ppf("poisson", lambda m: m > 0)),
    "rchi2": (1, 1, _chi2),
    "rt": (1, 1, _ppf("t", lambda df: df > 0)),
    "rbeta": (2, 2, _ppf("beta", lambda a, b: (a > 0) & (b > 0))),
    "rgamma": (2, 2, _gamma),
    "rexponential": (1, 1, _exponential),
    "rlogistic": (0, 2, _logistic),
    "rweibull": (2, 3, _weibull),
    "rnbinomial": (2, 2, _nbinomial),
    "rhypergeometric": (3, 3, _hyper),
}


def draw(name: str, args: list, n: int) -> np.ndarray:
    """n draws from distribution `name`; scalar args or vectors of length n."""
    lo, hi, fn = DISTRIBUTIONS[name]
    u = RNG.uniform(n)
    arrs = [np.broadcast_to(np.asarray(a, dtype=np.float64), (n,)) for a in args]
    with np.errstate(all="ignore"):
        return fn(u, *arrs)
