"""Funções disponíveis em expressões (manual [FN]).

Cada função recebe a lista de argumentos já avaliados (float ou str) e
devolve float ou str. A tabela FUNCTIONS registra nome -> (implementação,
mínimo de argumentos, máximo de argumentos).

Esta é a primeira leva: matemática, strings, lógicas, algumas
distribuições e datas. A cobertura completa do Stata 14 está no ROADMAP.
"""

from __future__ import annotations

import datetime as _dt
import fnmatch
import math
import re
from typing import Any, Callable

import numpy as np
from scipy import special, stats

from ..core import missing as M
from ..core.errors import StataError, type_mismatch
from ..core.formats import STATA_EPOCH, format_value, parse_format

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
    """Envolve função numérica: missing entra -> missing sai; erros -> missing."""
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


# -- matemática --------------------------------------------------------------

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
    """round(x, y) = y * floor(x/y + .5): meio arredondado para cima
    (round(-2.5) = -2), confirmado nos logs do Stata."""
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


# -- lógicas / programação ------------------------------------------------------

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
    """Maiúscula em letras ASCII que vêm depois de algo que não é letra ASCII.
    Como no Stata, caracteres acentuados contam como "não letra":
    proper("joão") = "JoãO" (use ustrtitle() para texto Unicode)."""
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
    """VERIFICAR: string(n) sem formato usa %9.0g."""
    x = _n(x)
    return format_value(x, parse_format(_s(fmt)), pad=False, dp=False).strip()


alias("strofreal", "string")


@register("abbrev", 2)
def _abbrev(s: Value, n: Value) -> str:
    """abbrev("variavel_muito_longa", 10) = "variavel~a": os primeiros n-2
    caracteres, '~' e o último caractere (observado no Stata). n mínimo: 5."""
    s, n = _s(s), int(_n(n))
    n = max(n, 5)
    if len(s) <= n:
        return s
    return s[: n - 2] + "~" + s[-1]


# expressões regulares: regexm guarda as capturas para regexs()
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


# -- distribuições ---------------------------------------------------------------

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


# -- datas --------------------------------------------------------------------

def _to_date(d: float) -> _dt.date | None:
    if M.is_missing(d):
        return None
    try:
        return STATA_EPOCH + _dt.timedelta(days=math.floor(d))
    except OverflowError:
        return None


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

_MONTH_NAMES = {m: i + 1 for i, m in enumerate(
    ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"])}


@register("date", 2, 3)
def _date(s: Value, mask: Value, topyear: Value | None = None) -> float:
    """Subconjunto de date(): máscaras com D, M e Y, mês numérico ou por nome."""
    text, mask = _s(s).strip().lower(), _s(mask).upper()
    parts = re.findall(r"[a-z]+|\d+", text)
    order = [c for c in mask if c in "DMY"]
    if len(parts) != len(order):
        return M.SYSMISS
    values: dict[str, int] = {}
    for p, c in zip(parts, order):
        if c == "M" and p.isalpha():
            if p[:3] not in _MONTH_NAMES:
                return M.SYSMISS
            values["M"] = _MONTH_NAMES[p[:3]]
        elif p.isdigit():
            values[c] = int(p)
        else:
            return M.SYSMISS
    try:
        return float((_dt.date(values["Y"], values["M"], values["D"]) - STATA_EPOCH).days)
    except (KeyError, ValueError):
        return M.SYSMISS


# -- constantes -------------------------------------------------------------------

register("maxdouble", 0)(lambda: M.MAXDOUBLE)
register("mindouble", 0)(lambda: M.MINDOUBLE)
register("epsdouble", 0)(lambda: M.EPSDOUBLE)
register("smallestdouble", 0)(lambda: 2.0 ** -1022)


# -- números aleatórios (core/rng.py) --------------------------------------------

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
