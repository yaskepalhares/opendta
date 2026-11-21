"""Functions available in expressions (manual [FN]).

Each function receives the list of already evaluated arguments (float or str) and
returns float or str. The FUNCTIONS table maps name -> (implementation,
minimum number of arguments, maximum number of arguments).

Covered: mathematics, strings, logical, distributions, dates and random
numbers. Full Stata 14 coverage is tracked in the ROADMAP.
"""

from __future__ import annotations

import datetime as _dt
import fnmatch
import math
import os
import re
from typing import Any, Callable

import numpy as np

from ..core import missing as M
from ..core import stdate
from ..core.errors import StataError, type_mismatch
from ..core.formats import STATA_EPOCH, format_value, parse_format
from ..lazy import lazy

# scipy is imported when a distribution function is first called (start-up time)
special = lazy("scipy.special")
stats = lazy("scipy.stats")

Value = Any
Impl = Callable[..., Value]
FUNCTIONS: dict[str, tuple[Impl, int, int]] = {}


def register(name: str, min_args: int, max_args: int | None = None):
    def deco(fn: Impl) -> Impl:
        FUNCTIONS[name] = (fn, min_args, min_args if max_args is None else max_args)
        return fn
    return deco


def alias(new: str, old: str) -> None:
    FUNCTIONS[new] = FUNCTIONS[old]


# -- utilidades -------------------------------------------------------------

def _n(v: Value) -> float:
    if isinstance(v, str):
        raise type_mismatch()
    return v


def _s(v: Value) -> str:
    if not isinstance(v, str):
        raise type_mismatch()
    return v


def _miss(*xs: float) -> bool:
    return any(M.is_missing(x) for x in xs)


def _num(fn: Callable[..., float]) -> Impl:
    """Wraps a numeric function: missing in -> missing out; errors -> missing."""
    def wrapped(*args: Value) -> float:
        xs = [_n(a) for a in args]
        if _miss(*xs):
            return M.SYSMISS
        try:
            r = fn(*xs)
        except (ValueError, OverflowError, ZeroDivisionError):
            return M.SYSMISS
        if r is None:
            return M.SYSMISS
        return M.normalize(float(r))
    return wrapped


def call(name: str, args: list[Value]) -> Value:
    entry = FUNCTIONS.get(name)
    if entry is None:
        raise StataError(133, f"unknown function {name}()")
    fn, lo, hi = entry
    if not (lo <= len(args) <= hi):
        raise StataError(198, "invalid syntax")
    return fn(*args)


# -- mathematics --------------------------------------------------------------

register("abs", 1)(_num(abs))
register("ceil", 1)(_num(math.ceil))
register("floor", 1)(_num(math.floor))
register("int", 1)(_num(lambda x: float(math.trunc(x))))
alias("trunc", "int")
register("exp", 1)(_num(math.exp))
register("ln", 1)(_num(lambda x: math.log(x) if x > 0 else None))
alias("log", "ln")
register("log10", 1)(_num(lambda x: math.log10(x) if x > 0 else None))
register("sqrt", 1)(_num(lambda x: math.sqrt(x) if x >= 0 else None))
register("sin", 1)(_num(math.sin))
register("cos", 1)(_num(math.cos))
register("tan", 1)(_num(math.tan))
register("asin", 1)(_num(math.asin))
register("acos", 1)(_num(math.acos))
register("atan", 1)(_num(math.atan))
register("atan2", 2)(_num(math.atan2))
register("sinh", 1)(_num(math.sinh))
register("cosh", 1)(_num(math.cosh))
register("tanh", 1)(_num(math.tanh))
register("sign", 1)(_num(lambda x: (x > 0) - (x < 0)))
register("lngamma", 1)(_num(lambda x: math.lgamma(x) if x > 0 else None))
register("lnfactorial", 1)(_num(lambda n: math.lgamma(n + 1) if n >= 0 and n == int(n) else None))
register("comb", 2)(_num(lambda n, k: float(math.comb(int(n), int(k)))
                         if n == int(n) and k == int(k) and 0 <= k <= n else None))
register("logit", 1)(_num(lambda p: math.log(p / (1 - p)) if 0 < p < 1 else None))
register("invlogit", 1)(_num(lambda x: 1 / (1 + math.exp(-x)) if x > -700 else 0.0))
register("float", 1)(_num(lambda x: float(np.float32(x))))
register("digamma", 1)(_num(lambda x: float(special.digamma(x))))


@register("round", 1, 2)
def _round(x: Value, y: Value = 1.0) -> float:
    """round(x, y) = y * floor(x/y + .5): halves rounded up
    (round(-2.5) = -2), confirmed in the Stata logs."""
    x, y = _n(x), _n(y)
    if _miss(x, y):
        return M.SYSMISS
    if y == 0:
        return x
    return M.normalize(math.floor(x / y + 0.5) * y)


@register("mod", 2)
def _mod(x: Value, y: Value) -> float:
    x, y = _n(x), _n(y)
    if _miss(x, y) or y == 0:
        return M.SYSMISS
    return M.normalize(x - y * math.floor(x / y))


@register("min", 1, 10_000)
def _min(*xs: Value) -> float:
    vals = [_n(x) for x in xs if not M.is_missing(_n(x))]
    return min(vals) if vals else M.SYSMISS


@register("max", 1, 10_000)
def _max(*xs: Value) -> float:
    vals = [_n(x) for x in xs if not M.is_missing(_n(x))]
    return max(vals) if vals else M.SYSMISS


# -- logical / programming ------------------------------------------------------

@register("cond", 3, 4)
def _cond(x: Value, a: Value, b: Value, c: Value | None = None) -> Value:
    x = _n(x)
    if c is not None and M.is_missing(x):
        return c
    return a if x != 0 else b


@register("inlist", 2, 255)
def _inlist(z: Value, *items: Value) -> float:
    for it in items:
        if isinstance(it, str) != isinstance(z, str):
            raise type_mismatch()
    return 1.0 if z in items else 0.0


@register("inrange", 3)
def _inrange(z: Value, a: Value, b: Value) -> float:
    if isinstance(z, str):
        return 1.0 if _s(a) <= z <= _s(b) else 0.0
    z, a, b = _n(z), _n(a), _n(b)
    if M.is_missing(z):
        return 1.0 if M.is_missing(b) and a <= z else 0.0
    lo = -math.inf if M.is_missing(a) else a
    return 1.0 if lo <= z <= b else 0.0


@register("missing", 1, 10_000)
def _missing(*xs: Value) -> float:
    for x in xs:
        if isinstance(x, str):
            if x == "":
                return 1.0
        elif M.is_missing(x):
            return 1.0
    return 0.0


alias("mi", "missing")


# -- strings -------------------------------------------------------------------

def _str1(fn: Callable[[str], Value]) -> Impl:
    return lambda s: fn(_s(s))


register("strlen", 1)(_str1(lambda s: float(len(s.encode("utf-8")))))
alias("length", "strlen")
register("ustrlen", 1)(_str1(lambda s: float(len(s))))
register("strupper", 1)(_str1(lambda s: s.upper()))
register("strlower", 1)(_str1(lambda s: s.lower()))
def _proper(s: str) -> str:
    """Uppercases ASCII letters that come after something that is not an ASCII letter.
    As in Stata, accented characters count as "not a letter":
    proper("joão") = "JoãO" (use ustrtitle() for Unicode text)."""
    out = []
    prev_letter = False
    for ch in s:
        is_letter = ch.isascii() and ch.isalpha()
        if is_letter:
            out.append(ch.lower() if prev_letter else ch.upper())
        else:
            out.append(ch)
        prev_letter = is_letter
    return "".join(out)


register("strproper", 1)(_str1(_proper))
register("ustrtitle", 1)(_str1(lambda s: s.title()))
register("ustrupper", 1)(_str1(lambda s: s.upper()))
register("ustrlower", 1)(_str1(lambda s: s.lower()))
alias("upper", "strupper")
alias("lower", "strlower")
alias("proper", "strproper")
register("strtrim", 1)(_str1(lambda s: s.strip(" ")))
register("strltrim", 1)(_str1(lambda s: s.lstrip(" ")))
register("strrtrim", 1)(_str1(lambda s: s.rstrip(" ")))
register("stritrim", 1)(_str1(lambda s: re.sub(r" {2,}", " ", s)))
alias("trim", "strtrim")
alias("ltrim", "strltrim")
alias("rtrim", "strrtrim")
alias("itrim", "stritrim")
register("strreverse", 1)(_str1(lambda s: s[::-1]))
alias("reverse", "strreverse")
register("wordcount", 1)(_str1(lambda s: float(len(s.split()))))


@register("word", 2)
def _word(s: Value, n: Value) -> str:
    words = _s(s).split()
    n = _n(n)
    if M.is_missing(n):
        return ""
    k = int(n)
    if k > 0 and k <= len(words):
        return words[k - 1]
    if k < 0 and -k <= len(words):
        return words[k]
    return ""


@register("substr", 3)
def _substr(s: Value, n1: Value, n2: Value) -> str:
    s, n1, n2 = _s(s), _n(n1), _n(n2)
    if M.is_missing(n1):
        return ""
    L = len(s)
    start = int(n1)
    if start < 0:
        start = L + start + 1
    if start < 1 or start > L:
        return ""
    length = L - start + 1 if M.is_missing(n2) else int(n2)
    if length <= 0:
        return ""
    return s[start - 1:start - 1 + length]


@register("strpos", 2)
def _strpos(s1: Value, s2: Value) -> float:
    return float(_s(s1).find(_s(s2)) + 1)


@register("subinstr", 4)
def _subinstr(s: Value, old: Value, new: Value, n: Value) -> str:
    s, old, new, n = _s(s), _s(old), _s(new), _n(n)
    if old == "":
        return s
    return s.replace(old, new) if M.is_missing(n) else s.replace(old, new, int(n))


@register("subinword", 4)
def _subinword(s: Value, old: Value, new: Value, n: Value) -> str:
    s, old, new, n = _s(s), _s(old), _s(new), _n(n)
    count = 0 if M.is_missing(n) else int(n)
    pattern = re.compile(r"(?<!\S)" + re.escape(old) + r"(?!\S)")
    return pattern.sub(new, s, count=count)


@register("strmatch", 2)
def _strmatch(s: Value, pattern: Value) -> float:
    return 1.0 if fnmatch.fnmatchcase(_s(s), _s(pattern)) else 0.0


@register("char", 1)
def _char(n: Value) -> str:
    n = _n(n)
    if M.is_missing(n) or not (1 <= n <= 255):
        return ""
    return chr(int(n))


@register("real", 1)
def _real(s: Value) -> float:
    t = _s(s).strip()
    if t == "" or t == ".":
        return M.SYSMISS
    try:
        return M.missing_code(t)
    except ValueError:
        pass
    try:
        return M.normalize(float(t))
    except ValueError:
        return M.SYSMISS


@register("string", 1, 2)
def _string(x: Value, fmt: Value = "%9.0g") -> str:
    """VERIFY: string(n) without a format uses %9.0g."""
    x = _n(x)
    return format_value(x, parse_format(_s(fmt)), pad=False, dp=False).strip()


alias("strofreal", "string")


@register("abbrev", 2)
def _abbrev(s: Value, n: Value) -> str:
    """abbrev("variavel_muito_longa", 10) = "variavel~a": the first n-2
    characters, '~' and the last character (observed in Stata). Minimum n: 5."""
    s, n = _s(s), int(_n(n))
    n = max(n, 5)
    if len(s) <= n:
        return s
    return s[: n - 2] + "~" + s[-1]


# regular expressions: regexm keeps the captures for regexs()
_last_match: list[re.Match[str] | None] = [None]


def _stata_regex(p: str) -> str:
    return p


@register("regexm", 2)
def _regexm(s: Value, pattern: Value) -> float:
    m = re.search(_stata_regex(_s(pattern)), _s(s))
    _last_match[0] = m
    return 1.0 if m else 0.0


@register("regexr", 3)
def _regexr(s: Value, pattern: Value, repl: Value) -> str:
    return re.sub(_stata_regex(_s(pattern)), _s(repl).replace("\\", "\\\\"), _s(s), count=1)


@register("regexs", 1)
def _regexs(n: Value) -> str:
    m = _last_match[0]
    k = int(_n(n))
    if m is None or k > (m.re.groups):
        raise StataError(198, "regexs(): subexpression out of range")
    return m.group(k) or ""


@register("strrpos", 2)
def _strrpos(s1: Value, s2: Value) -> float:
    return float(_s(s1).rfind(_s(s2)) + 1) if _s(s2) else 0.0


@register("indexnot", 2)
def _indexnot(s1: Value, s2: Value) -> float:
    """Position of the first character of s1 not found in s2; 0 if there is none."""
    s1, s2 = _s(s1), _s(s2)
    return float(next((i + 1 for i, ch in enumerate(s1) if ch not in s2), 0))


@register("plural", 2, 3)
def _plural(n: Value, s: Value, p: Value | None = None) -> str:
    """plural(n, "cat") = "cats"; with p, "+x" adds x, "-x" removes x, anything else replaces."""
    n, s = _n(n), _s(s)
    if n in (1.0, -1.0):
        return s
    if p is None:
        return s + "s"
    p = _s(p)
    if p.startswith("+"):
        return s + p[1:]
    if p.startswith("-"):
        return s[: len(s) - len(p) + 1] if p[1:] and s.endswith(p[1:]) else s
    return p


@register("strtoname", 1, 2)
def _strtoname(s: Value, p: Value = 1.0) -> str:
    """A valid Stata name: other characters become _, a leading digit gets _ in front (unless p
    is 0), and the result has at most 32 characters. VERIFY: strtoname("") is "_"."""
    s, p = _s(s), _n(p)
    out = re.sub(r"[^A-Za-z0-9_]", "_", s)
    if out == "":
        out = "_"
    elif out[0].isdigit() and p != 0:
        out = "_" + out
    return out[:32]


_SOUNDEX = {**dict.fromkeys("bfpv", "1"), **dict.fromkeys("cgjkqsxz", "2"),
            **dict.fromkeys("dt", "3"), "l": "4", **dict.fromkeys("mn", "5"), "r": "6"}


def _soundex_code(s: str, nara: bool) -> str:
    """Soundex code; with nara, letters with the same code separated by h or w count once (the
    rule of the U.S. National Archives). VERIFY: the variant soundex() uses."""
    letters = [c for c in s.lower() if c.isascii() and c.isalpha()]
    if not letters:
        return ""
    out, last = [letters[0].upper()], _SOUNDEX.get(letters[0], "")
    for c in letters[1:]:
        code = _SOUNDEX.get(c, "")
        if code and code != last:
            out.append(code)
        if not (nara and c in "hw"):
            last = code
    return ("".join(out) + "000")[:4]


register("soundex", 1)(_str1(lambda s: _soundex_code(s, False)))
register("soundex_nara", 1)(_str1(lambda s: _soundex_code(s, True)))


# -- Unicode strings (characters, not bytes) ---------------------------------------------------

@register("ustrpos", 2, 3)
def _ustrpos(s1: Value, s2: Value, n: Value = 1.0) -> float:
    s1, s2, n = _s(s1), _s(s2), _n(n)
    start = 0 if M.is_missing(n) else max(int(n), 1) - 1
    return float(s1.find(s2, start) + 1) if s2 else 0.0


@register("ustrrpos", 2, 3)
def _ustrrpos(s1: Value, s2: Value, n: Value | None = None) -> float:
    s1, s2 = _s(s1), _s(s2)
    end = len(s1) if n is None or M.is_missing(_n(n)) else int(_n(n))
    return float(s1.rfind(s2, 0, end) + 1) if s2 else 0.0


FUNCTIONS["usubstr"] = FUNCTIONS["substr"]


@register("ustrleft", 2)
def _ustrleft(s: Value, n: Value) -> str:
    s, n = _s(s), _n(n)
    return s if M.is_missing(n) else s[: max(int(n), 0)]


@register("ustrright", 2)
def _ustrright(s: Value, n: Value) -> str:
    s, n = _s(s), _n(n)
    if M.is_missing(n):
        return s
    return s[len(s) - int(n):] if int(n) > 0 else ""


register("ustrtrim", 1)(_str1(lambda s: s.strip()))
register("ustrltrim", 1)(_str1(lambda s: s.lstrip()))
register("ustrrtrim", 1)(_str1(lambda s: s.rstrip()))
register("ustrreverse", 1)(_str1(lambda s: s[::-1]))
register("ustrwordcount", 1, 2)(lambda s, loc=None: float(len(_s(s).split())))


@register("ustrword", 2, 3)
def _ustrword(s: Value, n: Value, loc: Value | None = None) -> str:
    """nth word; a negative n counts from the end."""
    words, n = _s(s).split(), _n(n)
    if M.is_missing(n) or n == 0:
        return ""
    k = int(n) - 1 if n > 0 else len(words) + int(n)
    return words[k] if 0 <= k < len(words) else ""


@register("uchar", 1)
def _uchar(n: Value) -> str:
    n = _n(n)
    if M.is_missing(n) or not 0 <= n <= 0x10FFFF:
        return ""
    return chr(int(n))


def _uflags(noc: Value | None) -> int:
    """The case-insensitive flag of ustrregex*(); missing counts as not 0, like true."""
    return re.IGNORECASE if noc is not None and _n(noc) != 0 else 0


_ulast_match: list[re.Match[str] | None] = [None]


@register("ustrregexm", 2, 3)
def _ustrregexm(s: Value, pattern: Value, noc: Value | None = None) -> float:
    m = re.search(_s(pattern), _s(s), _uflags(noc))
    _ulast_match[0] = m
    return 1.0 if m else 0.0


@register("ustrregexrf", 3, 4)
def _ustrregexrf(s: Value, pattern: Value, repl: Value, noc: Value | None = None) -> str:
    return re.sub(_s(pattern), _s(repl).replace("\\", "\\\\"), _s(s), count=1, flags=_uflags(noc))


@register("ustrregexra", 3, 4)
def _ustrregexra(s: Value, pattern: Value, repl: Value, noc: Value | None = None) -> str:
    return re.sub(_s(pattern), _s(repl).replace("\\", "\\\\"), _s(s), flags=_uflags(noc))


@register("ustrregexs", 1)
def _ustrregexs(n: Value) -> str:
    m, k = _ulast_match[0], int(_n(n))
    if m is None or not 0 <= k <= m.re.groups:
        return ""
    return m.group(k) or ""


# -- distributions ---------------------------------------------------------------

register("normal", 1)(_num(lambda z: stats.norm.cdf(z)))
register("normalden", 1, 3)(_num(lambda z, m=0.0, s=1.0: stats.norm.pdf(z, m, s) if s > 0 else None))
register("invnormal", 1)(_num(lambda p: stats.norm.ppf(p) if 0 < p < 1 else None))
register("lnnormal", 1)(_num(lambda z: stats.norm.logcdf(z)))
register("chi2", 2)(_num(lambda df, x: stats.chi2.cdf(x, df) if df > 0 else None))
register("chi2tail", 2)(_num(lambda df, x: stats.chi2.sf(x, df) if df > 0 else None))
register("invchi2", 2)(_num(lambda df, p: stats.chi2.ppf(p, df) if 0 <= p < 1 else None))
register("invchi2tail", 2)(_num(lambda df, p: stats.chi2.isf(p, df) if 0 < p <= 1 else None))
register("F", 3)(_num(lambda d1, d2, f: stats.f.cdf(f, d1, d2)))
register("Ftail", 3)(_num(lambda d1, d2, f: stats.f.sf(f, d1, d2)))
register("invF", 3)(_num(lambda d1, d2, p: stats.f.ppf(p, d1, d2) if 0 <= p < 1 else None))
register("invFtail", 3)(_num(lambda d1, d2, p: stats.f.isf(p, d1, d2) if 0 < p <= 1 else None))
register("t", 2)(_num(lambda df, t: stats.t.cdf(t, df)))
register("ttail", 2)(_num(lambda df, t: stats.t.sf(t, df)))
register("invt", 2)(_num(lambda df, p: stats.t.ppf(p, df) if 0 < p < 1 else None))
register("invttail", 2)(_num(lambda df, p: stats.t.isf(p, df) if 0 < p < 1 else None))
register("binomial", 3)(_num(lambda n, k, p: stats.binom.cdf(k, n, p)))
register("binomialtail", 3)(_num(lambda n, k, p: stats.binom.sf(k - 1, n, p)))
register("poisson", 2)(_num(lambda m, k: stats.poisson.cdf(k, m)))
register("poissontail", 2)(_num(lambda m, k: stats.poisson.sf(k - 1, m)))
register("betaden", 3)(_num(lambda a, b, x: stats.beta.pdf(x, a, b)))
register("ibeta", 3)(_num(lambda a, b, x: special.betainc(a, b, x)))
register("gammap", 2)(_num(lambda a, x: special.gammainc(a, x)))

register("binomialp", 3)(_num(lambda n, k, p: stats.binom.pmf(k, n, p)))
register("invbinomial", 3)(_num(lambda n, k, P: 1 - special.betaincinv(n - k, k + 1, P)
                                if 0 <= k < n and 0 < P < 1 else None))
register("invbinomialtail", 3)(_num(lambda n, k, P: special.betaincinv(k, n - k + 1, P)
                                    if 0 < k <= n and 0 < P < 1 else None))
register("poissonp", 2)(_num(lambda m, k: stats.poisson.pmf(k, m) if m > 0 else None))
register("invpoisson", 2)(_num(lambda k, P: special.gammainccinv(k + 1, P) if 0 < P < 1 else None))
register("invpoissontail", 2)(_num(lambda k, q: special.gammaincinv(k, q)
                                   if k >= 1 and 0 < q < 1 else None))
register("nbinomial", 3)(_num(lambda n, k, p: stats.nbinom.cdf(k, n, p)))
register("nbinomialp", 3)(_num(lambda n, k, p: stats.nbinom.pmf(k, n, p)))
register("nbinomialtail", 3)(_num(lambda n, k, p: stats.nbinom.sf(k - 1, n, p)))
register("invnbinomial", 3)(_num(lambda n, k, q: special.betaincinv(n, k + 1, q)
                                 if 0 < q < 1 else None))
register("invnbinomialtail", 3)(_num(lambda n, k, q: special.betaincinv(n, k, 1 - q)
                                     if k >= 1 and 0 < q < 1 else None))
register("hypergeometricp", 4)(_num(lambda N, K, n, k: stats.hypergeom.pmf(k, N, K, n)))
register("hypergeometric", 4)(_num(lambda N, K, n, k: stats.hypergeom.cdf(k, N, K, n)))
register("tden", 2)(_num(lambda df, t: stats.t.pdf(t, df) if df > 0 else None))
register("ntden", 3)(_num(lambda df, np_, t: stats.nct.pdf(t, df, np_)))
register("nt", 3)(_num(lambda df, np_, t: stats.nct.cdf(t, df, np_)))
register("nttail", 3)(_num(lambda df, np_, t: stats.nct.sf(t, df, np_)))
register("invnttail", 3)(_num(lambda df, np_, p: stats.nct.isf(p, df, np_) if 0 < p < 1 else None))
register("Fden", 3)(_num(lambda d1, d2, f: stats.f.pdf(f, d1, d2)))
register("nFden", 4)(_num(lambda d1, d2, np_, f: stats.ncf.pdf(f, d1, d2, np_)))
register("nF", 4)(_num(lambda d1, d2, np_, f: stats.ncf.cdf(f, d1, d2, np_)))
register("nFtail", 4)(_num(lambda d1, d2, np_, f: stats.ncf.sf(f, d1, d2, np_)))
register("invnFtail", 4)(_num(lambda d1, d2, np_, p: stats.ncf.isf(p, d1, d2, np_)
                              if 0 < p < 1 else None))
register("chi2den", 2)(_num(lambda df, x: stats.chi2.pdf(x, df) if df > 0 else None))
register("nchi2", 3)(_num(lambda df, np_, x: stats.ncx2.cdf(x, df, np_)))
register("nchi2tail", 3)(_num(lambda df, np_, x: stats.ncx2.sf(x, df, np_)))
register("invnchi2", 3)(_num(lambda df, np_, p: stats.ncx2.ppf(p, df, np_) if 0 < p < 1 else None))
register("invnchi2tail", 3)(_num(lambda df, np_, p: stats.ncx2.isf(p, df, np_)
                                 if 0 < p < 1 else None))


def _npnchi2(df: float, x: float, p: float) -> float | None:
    """Noncentrality np such that nchi2(df, np, x) = p (decreasing in np)."""
    from scipy.optimize import brentq
    if not (0 < p < 1) or df <= 0 or x <= 0 or stats.chi2.cdf(x, df) < p:
        return None
    hi = 1.0
    while stats.ncx2.cdf(x, df, hi) > p:
        hi *= 2
        if hi > 1e7:
            return None
    return brentq(lambda np_: stats.ncx2.cdf(x, df, np_) - p, 0.0, hi, xtol=1e-12)


register("npnchi2", 3)(_num(_npnchi2))
register("gammaden", 4)(_num(lambda a, b, g, x: stats.gamma.pdf(x, a, loc=g, scale=b)))
register("gammaptail", 2)(_num(lambda a, x: special.gammaincc(a, x)))
register("invgammap", 2)(_num(lambda a, p: special.gammaincinv(a, p) if 0 < p < 1 else None))
register("invgammaptail", 2)(_num(lambda a, p: special.gammainccinv(a, p) if 0 < p < 1 else None))
register("ibetatail", 3)(_num(lambda a, b, x: 1 - special.betainc(a, b, x)))
register("invibeta", 3)(_num(lambda a, b, p: special.betaincinv(a, b, p) if 0 <= p <= 1 else None))
register("invibetatail", 3)(_num(lambda a, b, p: special.betaincinv(a, b, 1 - p)
                                 if 0 <= p <= 1 else None))
register("lnnormalden", 1, 3)(_num(lambda z, m=0.0, s=1.0: stats.norm.logpdf(z, m, s) if s > 0 else None))
register("binormal", 3)(_num(lambda h, k, r: stats.multivariate_normal.cdf(
    [h, k], mean=[0, 0], cov=[[1, r], [r, 1]]) if -1 <= r <= 1 else None))
register("trigamma", 1)(_num(lambda x: float(special.polygamma(1, x))))
register("reldif", 2)(_num(lambda x, y: abs(x - y) / (abs(y) + 1)))
register("logistic", 1)(_num(lambda x: 1 / (1 + math.exp(-x)) if x > -700 else 0.0))
register("cloglog", 1)(_num(lambda x: math.log(-math.log(1 - x)) if 0 < x < 1 else None))
register("invcloglog", 1)(_num(lambda x: 1 - math.exp(-math.exp(x)) if x < 700 else 1.0))

# -- datas --------------------------------------------------------------------

def _to_date(d: float) -> _dt.date | None:
    if M.is_missing(d):
        return None
    try:
        return STATA_EPOCH + _dt.timedelta(days=math.floor(d))
    except OverflowError:
        return None


@register("bofd", 2)
def _bofd(cal: Value, td: Value) -> float:
    """Business date of the %td date td in calendar cal (missing if not a business day)."""
    from ..core import bcal
    if not isinstance(cal, str):
        raise StataError(109, "type mismatch")
    x = _n(td)
    if _miss(x):
        return M.SYSMISS
    b = bcal.get(cal).bofd(float(x))
    return M.SYSMISS if b != b else b


@register("dofb", 2)
def _dofb(b: Value, cal: Value) -> float:
    """%td date of the business date b in calendar cal."""
    from ..core import bcal
    if not isinstance(cal, str):
        raise StataError(109, "type mismatch")
    x = _n(b)
    if _miss(x):
        return M.SYSMISS
    d = bcal.get(cal).dofb(float(x))
    return M.SYSMISS if d != d else d


@register("mdy", 3)
def _mdy(m: Value, d: Value, y: Value) -> float:
    m, d, y = _n(m), _n(d), _n(y)
    if _miss(m, d, y):
        return M.SYSMISS
    try:
        return float((_dt.date(int(y), int(m), int(d)) - STATA_EPOCH).days)
    except ValueError:
        return M.SYSMISS


def _date_part(fn: Callable[[_dt.date], int]) -> Impl:
    def wrapped(d: Value) -> float:
        dt = _to_date(_n(d))
        return M.SYSMISS if dt is None else float(fn(dt))
    return wrapped


register("year", 1)(_date_part(lambda d: d.year))
register("month", 1)(_date_part(lambda d: d.month))
register("day", 1)(_date_part(lambda d: d.day))
register("dow", 1)(_date_part(lambda d: (d.weekday() + 1) % 7))
register("doy", 1)(_date_part(lambda d: d.timetuple().tm_yday))
register("quarter", 1)(_date_part(lambda d: (d.month - 1) // 3 + 1))

def _numf(fn: Callable[..., float | None]) -> Impl:
    """Like _num, but the arguments are truncated to integers first (date components)."""
    return _num(lambda *xs: fn(*[math.floor(x) for x in xs]))


def _opt_top(topyear: Value | None) -> float | None:
    if topyear is None:
        return None
    t = _n(topyear)
    return None if M.is_missing(t) else math.floor(t)


def _strdate(fn: Callable[..., float | None]) -> Impl:
    def wrapped(s: Value, mask: Value, topyear: Value | None = None) -> float:
        r = fn(_s(s), _s(mask), _opt_top(topyear))
        return M.SYSMISS if r is None else float(r)
    return wrapped


register("date", 2, 3)(_strdate(stdate.daily))
alias("daily", "date")
register("clock", 2, 3)(_strdate(stdate.clock))
register("Clock", 2, 3)(_strdate(lambda s, m, t: None if (c := stdate.clock(s, m, t)) is None
                                  else stdate.Cofc(c)))
register("monthly", 2, 3)(_strdate(lambda s, m, t: stdate.period(s, m, "M", 12, t)))
register("quarterly", 2, 3)(_strdate(lambda s, m, t: stdate.period(s, m, "Q", 4, t)))
register("halfyearly", 2, 3)(_strdate(lambda s, m, t: stdate.period(s, m, "H", 2, t)))
register("weekly", 2, 3)(_strdate(lambda s, m, t: stdate.period(s, m, "W", 52, t)))
register("yearly", 2, 3)(_strdate(lambda s, m, t: None if (f := stdate.parse(s, m, t)) is None
                                   else f.get("Y")))

# periods since 1960 from their components, and back
register("ym", 2)(_numf(lambda y, m: 12 * (y - 1960) + m - 1 if 1 <= m <= 12 else None))
register("yq", 2)(_numf(lambda y, q: 4 * (y - 1960) + q - 1 if 1 <= q <= 4 else None))
register("yh", 2)(_numf(lambda y, h: 2 * (y - 1960) + h - 1 if 1 <= h <= 2 else None))
register("yw", 2)(_numf(lambda y, w: 52 * (y - 1960) + w - 1 if 1 <= w <= 52 else None))
register("dofm", 1)(_numf(lambda p: stdate.dofperiod(p, 12)))
register("dofq", 1)(_numf(lambda p: stdate.dofperiod(p, 4)))
register("dofh", 1)(_numf(lambda p: stdate.dofperiod(p, 2)))
register("dofw", 1)(_numf(stdate.dofw))
register("dofy", 1)(_numf(lambda y: stdate.mdy(1, 1, y)))
register("week", 1)(_date_part(stdate.week))
register("halfyear", 1)(_date_part(lambda d: 1 if d.month <= 6 else 2))
register("yofd", 1)(_date_part(lambda d: d.year))
register("mofd", 1)(_date_part(lambda d: 12 * (d.year - 1960) + d.month - 1))
register("qofd", 1)(_date_part(lambda d: 4 * (d.year - 1960) + (d.month - 1) // 3))
register("hofd", 1)(_date_part(lambda d: 2 * (d.year - 1960) + (d.month > 6)))
register("wofd", 1)(_date_part(lambda d: 52 * (d.year - 1960) + stdate.week(d) - 1))

# times in milliseconds: %tc without leap seconds, %tC with them
_MS_DAY = stdate.MS_DAY
register("hms", 3)(_num(stdate.hms))
register("dhms", 4)(_num(lambda d, h, m, s: math.floor(d) * _MS_DAY + stdate.hms(h, m, s)))
register("mdyhms", 6)(_num(lambda mo, d, y, h, mi, s: None if (td := stdate.mdy(mo, d, y)) is None
                           else td * _MS_DAY + stdate.hms(h, mi, s)))
register("Cmdyhms", 6)(_num(lambda mo, d, y, h, mi, s: None if (td := stdate.mdy(mo, d, y)) is None
                            else stdate.Cofc(td * _MS_DAY + stdate.hms(h, mi, s))))
register("cofd", 1)(_num(lambda d: math.floor(d) * _MS_DAY))
register("Cofd", 1)(_num(lambda d: stdate.Cofc(math.floor(d) * _MS_DAY)))
register("dofc", 1)(_num(lambda c: math.floor(c / _MS_DAY)))
register("dofC", 1)(_num(lambda c: math.floor(stdate.cofC(c) / _MS_DAY)))
register("Cofc", 1)(_num(stdate.Cofc))
register("cofC", 1)(_num(stdate.cofC))
register("hh", 1)(_num(lambda c: math.floor((c % _MS_DAY) / 3_600_000)))
register("mm", 1)(_num(lambda c: math.floor((c % 3_600_000) / 60_000)))
register("ss", 1)(_num(lambda c: (c % 60_000) / 1000))
# VERIFY: hhC(), mmC() and ssC() inside a leap second (Stata shows 23:59:60)
register("hhC", 1)(_num(lambda c: math.floor((stdate.cofC(c) % _MS_DAY) / 3_600_000)))
register("mmC", 1)(_num(lambda c: math.floor((stdate.cofC(c) % 3_600_000) / 60_000)))
register("ssC", 1)(_num(lambda c: (stdate.cofC(c) % 60_000) / 1000))
register("hours", 1)(_num(lambda ms: ms / 3_600_000))
register("minutes", 1)(_num(lambda ms: ms / 60_000))
register("seconds", 1)(_num(lambda ms: ms / 1000))
register("msofhours", 1)(_num(lambda h: h * 3_600_000))
register("msofminutes", 1)(_num(lambda m: m * 60_000))
register("msofseconds", 1)(_num(lambda s: s * 1000))


# -- constantes -------------------------------------------------------------------

register("byteorder", 0)(lambda: float(__import__("opendta.core.binfmt", fromlist=["x"]).NATIVE_ORDER))
register("maxdouble", 0)(lambda: M.MAXDOUBLE)
register("mindouble", 0)(lambda: M.MINDOUBLE)
register("epsdouble", 0)(lambda: M.EPSDOUBLE)
register("smallestdouble", 0)(lambda: 2.0 ** -1022)
# limits of the storage types (.dta specification); VERIFY: minfloat() is -maxfloat()
_MAXFLOAT = float(np.frombuffer(bytes.fromhex("7effffff"), dtype=">f4")[0])
register("maxbyte", 0)(lambda: 100.0)
register("minbyte", 0)(lambda: -127.0)
register("maxint", 0)(lambda: 32740.0)
register("minint", 0)(lambda: -32767.0)
register("maxlong", 0)(lambda: 2147483620.0)
register("minlong", 0)(lambda: -2147483647.0)
register("maxfloat", 0)(lambda: _MAXFLOAT)
register("minfloat", 0)(lambda: -_MAXFLOAT)


# -- recoding ------------------------------------------------------------------------

@register("recode", 3, 10_000)
def _recode(x: Value, *cuts: Value) -> float:
    """The first cut that x does not exceed, or the last cut."""
    x, cs = _n(x), [_n(c) for c in cuts]
    if M.is_missing(x):
        return M.SYSMISS
    return next((c for c in cs if x <= c), cs[-1])


@register("irecode", 3, 10_000)
def _irecode(x: Value, *cuts: Value) -> float:
    """Index (0 to n) of the interval of x among the cuts."""
    x, cs = _n(x), [_n(c) for c in cuts]
    if M.is_missing(x):
        return M.SYSMISS
    return float(next((i for i, c in enumerate(cs) if x <= c), len(cs)))


@register("autocode", 4)
def _autocode(x: Value, n: Value, lo: Value, hi: Value) -> float:
    """Upper bound of the interval of x among n equal intervals from lo to hi. VERIFY: x
    below lo gives the first bound and x above hi gives hi."""
    x, n, lo, hi = _n(x), _n(n), _n(lo), _n(hi)
    if _miss(x, n, lo, hi) or n < 1 or hi <= lo:
        return M.SYSMISS
    k = int(n)
    for i in range(1, k + 1):
        u = lo + i * (hi - lo) / k
        if x <= u:
            return u
    return hi


# -- files (paths relative to the working folder) -----------------------------------------

register("fileexists", 1)(_str1(lambda f: float(os.path.isfile(f))))
register("filereadable", 1)(_str1(lambda f: float(os.path.isfile(f) and os.access(f, os.R_OK))))


@register("fileread", 1)
def _fileread(f: Value) -> str:
    """Contents of the file; "fileread() error 601" when it cannot be read (VERIFY: code)."""
    try:
        with open(_s(f), "rb") as fh:
            return fh.read().decode("utf-8", errors="replace")
    except OSError:
        return "fileread() error 601"


@register("filewrite", 2, 3)
def _filewrite(f: Value, text: Value, r: Value = 0.0) -> float:
    """Writes text to f (r = 1 replaces, r = 2 appends, otherwise f must not exist); returns
    the bytes written or a negative error code. VERIFY: the error codes."""
    f, data, r = _s(f), _s(text).encode("utf-8"), _n(r)
    if r not in (1.0, 2.0) and os.path.exists(f):
        return -602.0
    try:
        with open(f, "ab" if r == 2.0 else "wb") as fh:
            fh.write(data)
    except OSError:
        return -603.0
    return float(len(data))


# -- random numbers (core/rng.py) --------------------------------------------

def _random(name: str) -> Impl:
    def fn(*args: Value) -> float:
        from ..core import rng
        return float(rng.draw(name, [_n(a) for a in args], 1)[0])
    return fn


def _register_random() -> None:
    from ..core.rng import DISTRIBUTIONS
    for _name, (_lo, _hi, _) in DISTRIBUTIONS.items():
        register(_name, _lo, _hi)(_random(_name))


_register_random()
