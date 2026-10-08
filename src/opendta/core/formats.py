"""Formatos de exibição (%9.0g, %9.2f, %10.3e, %-12s, %td ...).

Referência: manual [D] format. Pontos que dependem de detalhes não
documentados (por exemplo, quando %g troca a notação fixa pela exponencial)
estão marcados com "VERIFICAR" e têm casos na suíte compat/.
"""

from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass

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


def _pad(s: str, fmt: Format) -> str:
    if len(s) >= fmt.width:
        return s
    return s.ljust(fmt.width) if fmt.left else s.rjust(fmt.width)


def format_value(value, fmt: Format | str = DEFAULT_NUMERIC, *, pad: bool = True,
                 sign_outside_width: bool = False) -> str:
    """Formata número ou string. `pad=False` remove o preenchimento (display
    sem formato explícito). `sign_outside_width=True` faz o sinal de menos
    não consumir largura em %g (comportamento observado no display)."""
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
        out = f"{x:,.{fmt.decimals}f}" if fmt.comma else f"{x:.{fmt.decimals}f}"
        if fmt.zero_pad and not fmt.left:
            out = out.zfill(fmt.width)
    elif fmt.kind == "e":
        out = f"{x:.{fmt.decimals}e}"
    elif fmt.kind == "g":
        out = _general(x, fmt, sign_outside_width)
    elif fmt.kind == "t":
        out = _date(x, fmt)
    else:
        out = repr(x)
    return _pad(out, fmt) if pad else out


def _strip_leading_zero(s: str) -> str:
    if s.startswith("0."):
        return s[1:]
    if s.startswith("-0."):
        return "-" + s[2:]
    return s


def _sig_digits(s: str) -> int:
    digits = re.sub(r"[^0-9]", "", s.split("e")[0])
    return len(digits.lstrip("0"))


def _general(x: float, fmt: Format, sign_outside_width: bool) -> str:
    """%w.0g: mostra o máximo de dígitos significativos que cabem em w.

    VERIFICAR: regra exata de escolha entre notação fixa e exponencial.
    Aqui escolhemos a que exibe mais dígitos significativos (empate: fixa).
    """
    neg = x < 0
    ax = -x if neg else x
    w = fmt.width
    if neg and not sign_outside_width:
        w -= 1
    w = max(w, 1)

    # inteiros que cabem
    if ax == int(ax) and len(str(int(ax))) <= w:
        s = str(int(ax))
        return "-" + s if neg else s

    candidates: list[str] = []

    # notação fixa
    int_digits = len(str(int(ax))) if ax >= 1 else 0
    if int_digits <= w:
        decimals = w - int_digits - 1
        if fmt.decimals > 0:
            decimals = min(decimals, fmt.decimals)
        if decimals >= 0:
            s = f"{ax:.{decimals}f}"
            if ax < 1:
                s = _strip_leading_zero(s)
                if len(s) > w:  # sem zero à esquerda sobra uma casa
                    s = s[:w]
            if len(s) <= w and _sig_digits(s) > 0:
                candidates.append(s)

    # notação exponencial: d.ddde+XX
    exp = f"{ax:e}".split("e")[1]
    exp_len = 1 + len(exp)  # 'e' + sinal + dígitos
    mant = w - exp_len - 2  # '1.' + decimais
    if mant >= 0:
        s = f"{ax:.{mant}e}"
        if len(s) <= w:
            candidates.append(s)

    if not candidates:
        s = f"{ax:.0e}"
    else:
        fixed = [c for c in candidates if "e" not in c]
        expo = [c for c in candidates if "e" in c]
        if fixed and (not expo or _sig_digits(fixed[0]) >= _sig_digits(expo[0])):
            s = fixed[0]
            if "." in s:
                s = s.rstrip("0").rstrip(".") or "0"
        else:
            s = expo[0]
    return "-" + s if neg else s


def _date(x: float, fmt: Format) -> str:
    if fmt.sub != "d":
        # outros formatos %t chegam na fase 1
        return _general(x, DEFAULT_NUMERIC, False)
    try:
        d = STATA_EPOCH + _dt.timedelta(days=int(x // 1))
    except OverflowError:
        return "."
    return f"{d.day:02d}{_MONTHS[d.month - 1]}{d.year}"


def number_to_macro(x: float) -> str:
    """Como um número vira texto em `local x = exp` e `=exp'.

    VERIFICAR: precisão exata usada pelo Stata (aqui, a menor representação
    que preserva o double, sem o zero à esquerda).
    """
    if M.is_missing(x):
        return M.missing_name(x)
    if x == int(x) and abs(x) < 1e16:
        return str(int(x))
    s = repr(x)
    if "e" in s:
        mant, exp = s.split("e")
        sign = exp[0] if exp[0] in "+-" else "+"
        digits = exp.lstrip("+-").rjust(2, "0")
        s = f"{mant}e{sign}{digits}"
    return _strip_leading_zero(s)
