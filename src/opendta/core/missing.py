"""Valores missing numéricos.

O OpenDTA representa números como double (float64) e usa, para os missing,
os mesmos padrões de bits que o formato .dta documenta para variáveis double:

    .   = 0x7fe0000000000000  (2**1023)
    .a  = 0x7fe0010000000000
    ...
    .z  = 0x7fe01a0000000000

Com isso, qualquer valor >= 2**1023 é missing, e a ordenação natural dos
floats já reproduz a regra do Stata: todo número < . < .a < ... < .z.
Também facilita ler e gravar .dta sem conversão.
"""

from __future__ import annotations

import math
import struct

_BASE_BITS = 0x7FE0000000000000
_STEP_BITS = 0x0000010000000000


def _from_bits(bits: int) -> float:
    return struct.unpack("<d", struct.pack("<Q", bits))[0]


def _to_bits(x: float) -> int:
    return struct.unpack("<Q", struct.pack("<d", x))[0]


SYSMISS: float = _from_bits(_BASE_BITS)          # .
MISSING_THRESHOLD: float = SYSMISS                # x >= isto  =>  missing
EXTENDED: dict[str, float] = {
    chr(ord("a") + i): _from_bits(_BASE_BITS + (i + 1) * _STEP_BITS) for i in range(26)
}
_BY_VALUE: dict[float, str] = {SYSMISS: "."}
_BY_VALUE.update({v: "." + k for k, v in EXTENDED.items()})

# Constantes c(maxdouble), c(mindouble), c(epsdouble) seguem da representação
# acima: o maior double não-missing é o predecessor de 2**1023.
MAXDOUBLE: float = math.nextafter(SYSMISS, 0.0)
MINDOUBLE: float = -MAXDOUBLE
EPSDOUBLE: float = 2.0 ** -52


def is_missing(x: float) -> bool:
    """True para ., .a–.z e qualquer valor fora do intervalo representável."""
    return x != x or x >= MISSING_THRESHOLD  # x != x captura NaN por segurança


def missing_code(name: str) -> float:
    """'.' -> SYSMISS, '.a' -> .a, ..."""
    if name == ".":
        return SYSMISS
    if len(name) == 2 and name[0] == "." and name[1] in EXTENDED:
        return EXTENDED[name[1]]
    raise ValueError(f"código de missing inválido: {name!r}")


def missing_name(x: float) -> str:
    """Nome de exibição de um missing ('.', '.a', ...). Pressupõe is_missing(x)."""
    if x != x:
        return "."
    return _BY_VALUE.get(x, ".")


def normalize(x: float) -> float:
    """Converte NaN/inf e qualquer valor >= 2**1023 que não seja um código
    conhecido em SYSMISS. Resultados de operações aritméticas passam por aqui."""
    if x != x or x in (math.inf, -math.inf):
        return SYSMISS
    if x >= MISSING_THRESHOLD and x not in _BY_VALUE:
        return SYSMISS
    if x <= -MISSING_THRESHOLD:
        return SYSMISS
    return x
