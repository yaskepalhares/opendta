"""Formatos de exibição (%9.0g, %9.2f, %10.3e, %-12s, %td ...).

Referência: manual [D] format. Pontos que dependem de detalhes não
documentados (por exemplo, quando %g troca a notação fixa pela exponencial)
estão marcados com "VERIFICAR" e têm casos na suíte compat/.
"""

from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass
from decimal import ROUND_HALF_UP, Decimal

from . import missing as M
from .errors import StataError

_NUMFMT = re.compile(r"^%(-)?(0)?(\d+)\.(\d+)([fge])(c)?$")
_STRFMT = re.compile(r"^%(-|~)?(\d+)s$")
_DATEFMT = re.compile(r"^%(-)?t([dcCwmqhy])(.*)$")

STATA_EPOCH = _dt.date(1960, 1, 1)
_MONTHS = ["jan", "feb", "mar", "apr", "may", "jun",
           "jul", "aug", "sep", "oct", "nov", "dec"]


@dataclass(frozen=True)
class Format:
    kind: str          # f, g, e, s, t
    width: int = 9
    decimals: int = 0
    left: bool = False
    zero_pad: bool = False
    comma: bool = False
    sub: str = ""      # para %t: d, c, ...
    text: str = ""

    def __str__(self) -> str:
        return self.text


def parse_format(text: str) -> Format:
    t = text.strip()
    m = _NUMFMT.match(t)
    if m:
        left, zero, w, d, kind, comma = m.groups()
        return Format(kind, int(w), int(d), bool(left), bool(zero), bool(comma), text=t)
    m = _STRFMT.match(t)
    if m:
        just, w = m.groups()
        return Format("s", int(w), left=(just == "-"), text=t)
    m = _DATEFMT.match(t)
    if m:
        left, sub, _rest = m.groups()
        return Format("t", 9, left=bool(left), sub=sub, text=t)
    raise StataError(120, f"invalid %format")


def is_format(text: str) -> bool:
    try:
        parse_format(text)
        return True
    except StataError:
        return False


DEFAULT_NUMERIC = parse_format("%9.0g")
DISPLAY_NUMERIC = parse_format("%10.0g")   # display sem formato explícito
MACRO_NUMERIC = parse_format("%18.0g")     # local x = exp, `=exp'

# set dp comma|period: muda o separador decimal dos formatos de exibição
# (não afeta string() nem macros, confirmado nos logs do Stata).
_DP = {"comma": False}


def set_decimal_comma(on: bool) -> None:
    _DP["comma"] = bool(on)


def decimal_comma() -> bool:
    return _DP["comma"]


def _swap_dp(s: str) -> str:
    return s.translate(str.maketrans({".": ",", ",": "."}))


def _pad(s: str, fmt: Format) -> str:
    if len(s) >= fmt.width:
        return s
    return s.ljust(fmt.width) if fmt.left else s.rjust(fmt.width)


def format_value(value, fmt: Format | str = DEFAULT_NUMERIC, *, pad: bool = True,
                 dp: bool = True, sign_outside_width: bool = False) -> str:
    """Formata número ou string. `pad=False` remove o preenchimento (display
    sem formato explícito). `dp=False` ignora `set dp comma` (string(), macros).
    `sign_outside_width` é aceito por compatibilidade e não tem mais efeito."""
    if isinstance(fmt, str):
        fmt = parse_format(fmt)
    if isinstance(value, str):
        out = value
        return _pad(out, fmt) if pad and fmt.kind == "s" else out
    if fmt.kind == "s":
        fmt = DEFAULT_NUMERIC
    x = float(value)
    if M.is_missing(x):
        out = M.missing_name(x)
    elif fmt.kind == "f":
        out = _fmt_f(x, fmt.decimals, fmt.comma)
        if fmt.zero_pad and not fmt.left:
            out = out.zfill(fmt.width)
    elif fmt.kind == "e":
        out = _fmt_e(x, fmt.decimals)
    elif fmt.kind == "g":
        out = general(x, fmt.width)
    elif fmt.kind == "t":
        out = _date(x, fmt)
    else:
        out = repr(x)
    if dp and _DP["comma"] and fmt.kind in ("f", "g", "e") and not M.is_missing(x):
        out = _swap_dp(out)
    return _pad(out, fmt) if pad else out


def _fmt_f(x: float, decimals: int, comma: bool = False) -> str:
    """Notação fixa com arredondamento de meio para cima sobre o valor binário
    exato (o Stata dá 1.235e+03 para 1234.5 em %10.3e; o Python daria 1.234)."""
    q = Decimal(x).quantize(Decimal(1).scaleb(-decimals), rounding=ROUND_HALF_UP)
    return f"{q:,.{decimals}f}" if comma else f"{q:.{decimals}f}"


def _fmt_e(x: float, decimals: int) -> str:
    if x == 0:
        return f"{0:.{decimals}e}"
    d = Decimal(x)
    exp = d.adjusted()
    mant = (d.scaleb(-exp)).quantize(Decimal(1).scaleb(-decimals), rounding=ROUND_HALF_UP)
    if abs(mant) >= 10:
        exp += 1
        mant = (d.scaleb(-exp)).quantize(Decimal(1).scaleb(-decimals), rounding=ROUND_HALF_UP)
    return f"{mant:.{decimals}f}" + _exp_text(exp)


def _exp_text(e: int) -> str:
    return f"e{'-' if e < 0 else '+'}{abs(e):02d}"


def general(x: float, width: int) -> str:
    """Formato %w.0g, com a regra observada nos logs do Stata 14:

    * uma posição fica reservada para o sinal: o número usa no máximo
      A = w - 1 caracteres, e o '-' ocupa a posição reservada;
    * no máximo P = w - 2 dígitos significativos;
    * notação exponencial quando o expoente (após arredondar para P dígitos)
      é < -4 ou >= P, como o %g da linguagem C; senão, notação fixa sem zeros
      à direita e sem o zero antes do ponto (.5);
    * na exponencial, a mantissa ocupa todo o espaço: 1.235e+09 em %10.0g.
    """
    if x == 0:
        return "0"
    neg = x < 0
    ax = -x if neg else x
    A = max(width - 1, 1)
    P = max(width - 2, 1)

    X = int(_fmt_e(ax, P - 1).split("e")[1])
    out = None
    if -4 <= X < P:
        decimals = max(0, P - 1 - X)
        while decimals >= 0:
            s = _fmt_f(ax, decimals)
            if "." in s:
                s = s.rstrip("0").rstrip(".")
            if s.startswith("0."):
                s = s[1:]
            if len(s) <= A and s not in ("", "0", "."):
                out = s
                break
            decimals -= 1
    if out is None:
        # exponencial: d.ddd...e+XX ocupando A caracteres
        for d in range(max(A - 2 - 4, 0), -1, -1):
            s = _fmt_e(ax, d)
            if len(s) <= A or d == 0:
                out = s
                break
    return "-" + out if neg else out


def _date(x: float, fmt: Format) -> str:
    if fmt.sub != "d":
        # outros formatos %t chegam na fase 1
        return general(x, DEFAULT_NUMERIC.width)
    try:
        d = STATA_EPOCH + _dt.timedelta(days=int(x // 1))
    except OverflowError:
        return "."
    return f"{d.day:02d}{_MONTHS[d.month - 1]}{d.year}"


def number_to_macro(x: float) -> str:
    """Como um número vira texto em `local x = exp` e `=exp': formato %18.0g
    (16 dígitos significativos), confirmado nos logs do Stata."""
    if M.is_missing(x):
        return M.missing_name(x)
    return general(x, MACRO_NUMERIC.width)
