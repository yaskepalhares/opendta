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
        left, sub, rest = m.groups()
        width = _DATE_WIDTH.get(sub, 9) if not rest else 0
        return Format("t", width, left=bool(left), sub=sub, text=t)
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


# Formatos de data e hora ([D] datetime display formats). Sem códigos, cada
# tipo usa o padrão abaixo; com códigos (%tdnn/dd/CCYY), segue os códigos.
_DATE_DEFAULT = {"d": "DDmonCCYY", "c": "DDmonCCYY_HH:MM:SS", "C": "DDmonCCYY_HH:MM:SS",
                 "w": "CCYY!www", "m": "CCYY!mnn", "q": "CCYY!qq", "h": "CCYY!hh", "y": "CCYY"}
_DATE_WIDTH = {"d": 9, "c": 18, "C": 18, "w": 7, "m": 7, "q": 6, "h": 6, "y": 4}
_DATE_CODES = sorted([
    "DAYNAME", "Dayname", "dayname", "Month", "month", "Mon", "mon", "Day", "day", "Da", "da",
    "JJJ", "jjj", "CC", "cc", "YY", "yy", "NN", "nn", "DD", "dd", "WW", "ww",
    "HH", "Hh", "hH", "hh", "MM", "mm", "SS", "ss", ".sss", ".ss", ".s",
    "a.m.", "A.M.", "am", "AM", "h", "q"], key=len, reverse=True)
_DAYNAMES = ["Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"]


def _date_parts(x: float, sub: str):
    """(datetime, semana, trimestre, semestre) do valor x no tipo `sub`."""
    week = None
    if sub == "d":
        d = _dt.datetime.combine(STATA_EPOCH, _dt.time()) + _dt.timedelta(days=int(x // 1))
    elif sub in ("c", "C"):
        d = _dt.datetime(1960, 1, 1) + _dt.timedelta(milliseconds=int(x // 1))
    elif sub == "w":
        k = int(x // 1)
        year, w = 1960 + k // 52, k % 52
        d = _dt.datetime(year, 1, 1) + _dt.timedelta(days=7 * w)
        week = w + 1
    elif sub == "m":
        k = int(x // 1)
        d = _dt.datetime(1960 + k // 12, k % 12 + 1, 1)
    elif sub == "q":
        k = int(x // 1)
        d = _dt.datetime(1960 + k // 4, 3 * (k % 4) + 1, 1)
    elif sub == "h":
        k = int(x // 1)
        d = _dt.datetime(1960 + k // 2, 6 * (k % 2) + 1, 1)
    else:   # y
        d = _dt.datetime(int(x // 1), 1, 1)
    if week is None:
        week = min((d.timetuple().tm_yday - 1) // 7 + 1, 52)
    return d, week


def _date_code(code: str, d: _dt.datetime, week: int, ms: int) -> str:
    y = d.year
    if code == "CC":
        return f"{y // 100:02d}"
    if code == "cc":
        return str(y // 100)
    if code == "YY":
        return f"{y % 100:02d}"
    if code == "yy":
        return str(y % 100)
    if code in ("JJJ", "jjj"):
        j = d.timetuple().tm_yday
        return f"{j:03d}" if code == "JJJ" else str(j)
    if code in ("Month", "month", "Mon", "mon"):
        full = _dt.date(2000, d.month, 1).strftime("%B")
        name = full if code.lower() == "month" else full[:3]
        return name if code[0].isupper() else name.lower()
    if code == "NN":
        return f"{d.month:02d}"
    if code == "nn":
        return str(d.month)
    if code == "DD":
        return f"{d.day:02d}"
    if code == "dd":
        return str(d.day)
    if code in ("DAYNAME", "Dayname", "dayname", "Day", "day", "Da", "da"):
        full = _DAYNAMES[d.weekday()]
        name = full if code.lower() == "dayname" else full[:3] if code.lower() == "day" else full[:2]
        return name.upper() if code == "DAYNAME" else name.lower() if code[0].islower() else name
    if code == "WW":
        return f"{week:02d}"
    if code == "ww":
        return str(week)
    if code == "h":
        return str(1 if d.month <= 6 else 2)
    if code == "q":
        return str((d.month - 1) // 3 + 1)
    h12 = d.hour % 12 or 12
    if code == "HH":
        return f"{d.hour:02d}"
    if code == "hH":
        return str(d.hour)
    if code == "Hh":
        return f"{h12:02d}"
    if code == "hh":
        return str(h12)
    if code == "MM":
        return f"{d.minute:02d}"
    if code == "mm":
        return str(d.minute)
    if code == "SS":
        return f"{d.second:02d}"
    if code == "ss":
        return str(d.second)
    if code.startswith("."):
        digits = len(code) - 1
        return "." + f"{ms:03d}"[:digits]
    pm = d.hour >= 12
    return {"am": "pm" if pm else "am", "AM": "PM" if pm else "AM",
            "a.m.": "p.m." if pm else "a.m.", "A.M.": "P.M." if pm else "A.M."}[code]


def _date(x: float, fmt: Format) -> str:
    sub = fmt.sub
    m = _DATEFMT.match(fmt.text) if fmt.text else None
    codes = (m.group(3) if m else "") or _DATE_DEFAULT.get(sub, "")
    try:
        d, week = _date_parts(x, sub)
    except (OverflowError, ValueError):
        return "."
    ms = int(x // 1) % 1000 if sub in ("c", "C") else 0
    out: list[str] = []
    i = 0
    while i < len(codes):
        if codes[i] == "!" and i + 1 < len(codes):
            out.append(codes[i + 1])
            i += 2
            continue
        for c in _DATE_CODES:
            if codes.startswith(c, i):
                out.append(_date_code(c, d, week, ms))
                i += len(c)
                break
        else:
            ch = codes[i]
            out.append(" " if ch == "_" else ch)
            i += 1
    return "".join(out)


def number_to_macro(x: float) -> str:
    """Como um número vira texto em `local x = exp` e `=exp': formato %18.0g
    (16 dígitos significativos), confirmado nos logs do Stata."""
    if M.is_missing(x):
        return M.missing_name(x)
    return general(x, MACRO_NUMERIC.width)
