"""Numeric missing values.

OpenDTA represents numbers as doubles (float64) and uses, for missing values,
the same bit patterns that the .dta format documents for double variables:

    .   = 0x7fe0000000000000  (2**1023)
    .a  = 0x7fe0010000000000
    ...
    .z  = 0x7fe01a0000000000

Thus any value >= 2**1023 is missing, and the natural ordering of
floats already reproduces Stata's rule: every number < . < .a < ... < .z.
It also makes reading and writing .dta files easy, without conversion.
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
MISSING_THRESHOLD: float = SYSMISS                # x >= this  =>  missing
EXTENDED: dict[str, float] = {
    chr(ord("a") + i): _from_bits(_BASE_BITS + (i + 1) * _STEP_BITS) for i in range(26)
}
_BY_VALUE: dict[float, str] = {SYSMISS: "."}
_BY_VALUE.update({v: "." + k for k, v in EXTENDED.items()})

# The constants c(maxdouble), c(mindouble), c(epsdouble) follow from the representation above: the
# largest non-missing double is the predecessor of 2**1023.
MAXDOUBLE: float = math.nextafter(SYSMISS, 0.0)
MINDOUBLE: float = -MAXDOUBLE
EPSDOUBLE: float = 2.0 ** -52


def is_missing(x: float) -> bool:
    """True for ., .a–.z and any value outside the representable range."""
    return x != x or x >= MISSING_THRESHOLD  # x != x catches NaN for safety


def missing_code(name: str) -> float:
    """'.' -> SYSMISS, '.a' -> .a, ..."""
    if name == ".":
        return SYSMISS
    if len(name) == 2 and name[0] == "." and name[1] in EXTENDED:
        return EXTENDED[name[1]]
    raise ValueError(f"código de missing inválido: {name!r}")


def missing_name(x: float) -> str:
    """Display name of a missing value ('.', '.a', ...). Assumes is_missing(x)."""
    if x != x:
        return "."
    return _BY_VALUE.get(x, ".")


def normalize(x: float) -> float:
    """Converts NaN/inf and any value >= 2**1023 that is not a known code
    into SYSMISS. Results of arithmetic operations go through here."""
    if x != x or x in (math.inf, -math.inf):
        return SYSMISS
    if x >= MISSING_THRESHOLD and x not in _BY_VALUE:
        return SYSMISS
    if x <= -MISSING_THRESHOLD:
        return SYSMISS
    return x
