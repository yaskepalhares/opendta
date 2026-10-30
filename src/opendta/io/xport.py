"""Reader of SAS XPORT version 5 files with every member (the public layout of SAS technical paper
TS-140): 80-byte records, a header for the library and, for each member, its descriptors, one
140-byte NAMESTR per variable and the observations. Numbers are IBM floating point, big-endian.

pyreadstat (io/statpkg.py) reads the first member of a file; this module gives the list of members
(import sasxport, describe), any one member (member()) and the formats.xpf catalog of value labels.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path

from ..core.errors import StataError

_REC = 80
_MEMBER = b"HEADER RECORD*******MEMBER  HEADER RECORD!!!!!!!"
_NAMESTR = b"HEADER RECORD*******NAMESTR HEADER RECORD!!!!!!!"
_OBS = b"HEADER RECORD*******OBS     HEADER RECORD!!!!!!!"
_LIBRARY = b"HEADER RECORD*******LIBRARY HEADER RECORD!!!!!!!"


@dataclass
class XVar:
    name: str
    numeric: bool
    length: int
    label: str
    format: str            # the SAS format name, without width or period
    position: int


@dataclass
class Member:
    name: str
    label: str
    created: str = ""      # the date and time the file records, as written there
    variables: list[XVar] = field(default_factory=list)
    rows: list[list] = field(default_factory=list)      # float (None = missing) or str


def ibm_to_float(b: bytes) -> float | None:
    """An IBM floating-point number of 2 to 8 bytes; None for a SAS missing value."""
    b = b.ljust(8, b"\x00")
    if b[1:] == b"\x00" * 7 and b[0] in (0x2E, 0x5F) or (0x41 <= b[0] <= 0x5A and b[1:] == b"\x00" * 7):
        return None
    hi, lo = struct.unpack(">II", b)
    sign = -1.0 if hi & 0x80000000 else 1.0
    exponent = (hi >> 24) & 0x7F
    mantissa = ((hi & 0x00FFFFFF) << 32) | lo
    if mantissa == 0:
        return 0.0
    return sign * mantissa / float(1 << 56) * 16.0 ** (exponent - 64)


def _text(b: bytes) -> str:
    return b.decode("latin-1").rstrip(" \x00")


def read(path: Path) -> list[Member]:
    """Every member of the file, with its variables and observations."""
    data = Path(path).read_bytes()
    if not data.startswith(_LIBRARY):
        raise StataError(610, f"file {path} not SAS XPORT format")
    members: list[Member] = []
    pos = data.find(_MEMBER)
    while pos >= 0:
        namestr_len = int(data[pos + 74:pos + 78] or b"140")
        d1 = data[pos + 2 * _REC:pos + 3 * _REC]          # after the DSCRPTR header
        d2 = data[pos + 3 * _REC:pos + 4 * _REC]
        m = Member(_text(d1[8:16]), _text(d2[32:72]), _text(d1[64:80]))
        nh = data.find(_NAMESTR, pos)
        nvars = int(data[nh + 54:nh + 58])
        start = nh + _REC
        for k in range(nvars):
            r = data[start + k * namestr_len:start + (k + 1) * namestr_len]
            ntype, _, nlng, _ = struct.unpack(">hhhh", r[:8])
            npos = struct.unpack(">i", r[84:88])[0]
            m.variables.append(XVar(_text(r[8:16]), ntype == 1, nlng, _text(r[16:56]),
                                    _text(r[56:64]), npos))
        oh = data.find(_OBS, start)
        body_start = oh + _REC
        nxt = data.find(_MEMBER, body_start)
        body = data[body_start:nxt if nxt >= 0 else len(data)]
        reclen = sum(v.length for v in m.variables)
        n = len(body) // reclen if reclen else 0
        # the last record is padded with blanks: a row made only of blanks that starts inside the
        # last 80 bytes is padding
        while n and (n - 1) * reclen >= len(body) - _REC and \
                body[(n - 1) * reclen:n * reclen].strip(b" ") == b"":
            n -= 1
        for i in range(n):
            rec = body[i * reclen:(i + 1) * reclen]
            row = []
            for v in m.variables:
                cell = rec[v.position:v.position + v.length]
                row.append(ibm_to_float(cell) if v.numeric else _text(cell))
            m.rows.append(row)
        members.append(m)
        pos = nxt
    return members


def value_labels(path: Path) -> dict[str, dict[int, str]]:
    """The numeric formats of a formats.xpf file (a CNTLOUT data set: FMTNAME, START, LABEL, TYPE)
    as value labels: format name (lowercase) -> {value: label}. Ranges and character formats are
    left out (Stata's value labels map single integers)."""
    out: dict[str, dict[int, str]] = {}
    for m in read(path):
        names = [v.name.upper() for v in m.variables]
        if not {"FMTNAME", "START", "LABEL"} <= set(names):
            continue
        col = {n: k for k, n in enumerate(names)}
        for row in m.rows:
            if "TYPE" in col and str(row[col["TYPE"]]).strip().upper() not in ("N", ""):
                continue
            start, end = row[col["START"]], row[col["END"]] if "END" in col else None
            try:
                x = float(start)
            except (TypeError, ValueError):
                continue
            if end not in (None, "") and str(end).strip() != str(start).strip():
                continue
            if x == int(x):
                out.setdefault(str(row[col["FMTNAME"]]).strip().lower(), {})[int(x)] = str(row[col["LABEL"]])
    return out
