"""Leitura e gravação de arquivos .dta.

Implementação própria dos formatos 117 (Stata 13), 118 (Stata 14/15) e 119
(Stata 15+, mais de 32.767 variáveis), escrita a partir da especificação
pública do formato. Formatos anteriores (até 115) são lidos pela ReadStat
(pyreadstat), que também serve, nos testes, de verificação independente
de tudo o que este módulo lê e grava.

Estrutura dos formatos 117+ (seções delimitadas por marcadores de texto):

    <stata_dta><header>...</header><map>...</map><variable_types>...
    <varnames>...<sortlist>...<formats>...<value_label_names>...
    <variable_labels>...<characteristics>...<data>...<strls>...
    <value_labels>...</stata_dta>

VERIFICAR: os detalhes de layout estão cobertos por testes de ida e volta
contra a ReadStat e por arquivos gravados e abertos no Stata 14.
"""

from __future__ import annotations

import datetime as _dt
import struct
from pathlib import Path

import numpy as np

from ..core import missing as M
from ..core.dataset import Dataset, Variable, default_format, str_len
from ..core.errors import StataError

# códigos de tipo (117+)
T_DOUBLE, T_FLOAT, T_LONG, T_INT, T_BYTE = 65526, 65527, 65528, 65529, 65530
T_STRL = 32768
_TYPE_CODE = {"double": T_DOUBLE, "float": T_FLOAT, "long": T_LONG, "int": T_INT, "byte": T_BYTE}
_CODE_TYPE = {v: k for k, v in _TYPE_CODE.items()}
_NP = {"double": "f8", "float": "f4", "long": "i4", "int": "i2", "byte": "i1"}
_SIZE = {"double": 8, "float": 4, "long": 4, "int": 2, "byte": 1}

# códigos de missing nos tipos inteiros: '.' e depois .a ... .z
_INT_MISS = {"byte": 101, "int": 32741, "long": 2147483621}
_FLOAT_MISS_BITS = 0x7F000000
_FLOAT_MISS_STEP = 0x00000800

_MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


class _Layout:
    """Tamanhos de campos que mudam entre as versões 117, 118 e 119."""

    def __init__(self, release: int):
        self.release = release
        self.k_size = 4 if release == 119 else 2
        self.n_size = 4 if release == 117 else 8
        self.label_len_size = 1 if release == 117 else 2
        self.name = 33 if release == 117 else 129
        self.fmt = 49 if release == 117 else 57
        self.varlabel = 81 if release == 117 else 321
        self.sort_size = 4 if release == 119 else 2


# ---------------------------------------------------------------------------
# Conversão de valores numéricos
# ---------------------------------------------------------------------------

def _to_disk(values: np.ndarray, vtype: str) -> np.ndarray:
    """Valores do OpenDTA (float64 com códigos de missing) para o tipo em disco."""
    v = np.asarray(values, dtype=np.float64)
    miss = v >= M.SYSMISS
    if vtype == "double":
        return v.astype("<f8")
    # índice do missing: 0 para '.', 1..26 para .a..z
    idx = np.zeros(len(v), dtype=np.int64)
    if miss.any():
        codes = [M.SYSMISS] + [M.EXTENDED[c] for c in "abcdefghijklmnopqrstuvwxyz"]
        for k, code in enumerate(codes):
            idx[v == code] = k
    if vtype == "float":
        out = np.where(miss, 0.0, v).astype("<f4")
        if miss.any():
            bits = out.view("<u4").copy()
            bits[miss] = _FLOAT_MISS_BITS + idx[miss] * _FLOAT_MISS_STEP
            out = bits.view("<f4")
        return out
    base = _INT_MISS[vtype]
    return np.where(miss, base + idx, np.where(miss, 0, v)).astype("<" + _NP[vtype])


def _from_disk(raw: np.ndarray, vtype: str) -> np.ndarray:
    """Valores em disco (little-endian) para float64 com os códigos de missing
    do OpenDTA. Em double os códigos são os mesmos do formato."""
    codes = np.array([M.SYSMISS] + [M.EXTENDED[c] for c in "abcdefghijklmnopqrstuvwxyz"])
    if vtype == "double":
        return raw.astype(np.float64)
    if vtype == "float":
        bits = np.ascontiguousarray(raw.astype("<f4")).view("<u4")
        v = raw.astype(np.float64)
        # só os positivos acima do maior float válido; com sinal (bit 31) é número
        miss = (bits >= _FLOAT_MISS_BITS) & (bits < 0x80000000)
        if miss.any():
            k = ((bits[miss] - _FLOAT_MISS_BITS) // _FLOAT_MISS_STEP).astype(np.int64).clip(0, 26)
            v[miss] = codes[k]
        return v
    ints = raw.astype(np.int64)
    v = ints.astype(np.float64)
    base = _INT_MISS[vtype]
    miss = ints >= base
    if miss.any():
        v[miss] = codes[(ints[miss] - base).clip(0, 26)]
    return v


# ---------------------------------------------------------------------------
# Gravação (formato 118)
# ---------------------------------------------------------------------------

_ENC = {"enc": "utf-8"}   # formato 117 (Stata 13) não é Unicode: latin-1


def _enc(text: str) -> bytes:
    return text.encode(_ENC["enc"], errors="replace")


def _fixed(text: str, size: int) -> bytes:
    b = _enc(text)[: size - 1]
    return b + b"\x00" * (size - len(b))


def _timestamp(now: _dt.datetime | None = None) -> str:
    now = now or _dt.datetime.now()
    return f"{now.day:02d} {_MONTHS[now.month - 1]} {now.year} {now.hour:02d}:{now.minute:02d}"


def write_dta(ds: Dataset, path: str | Path, *, release: int = 118,
              timestamp: str | None = None) -> str:
    """Grava ds em path e devolve o carimbo de data gravado no cabeçalho."""
    if release not in (117, 118, 119):
        raise StataError(198, f"dta release {release} not supported")
    L = _Layout(release)
    if release == 118 and ds.nvars > 32767:
        release, L = 119, _Layout(119)
    _ENC["enc"] = "latin-1" if release == 117 else "utf-8"
    K, N = ds.nvars, ds.nobs
    out = bytearray()
    offsets: dict[str, int] = {}

    def tag(name: str) -> None:
        out.extend(f"<{name}>".encode())

    def end(name: str) -> None:
        out.extend(f"</{name}>".encode())

    offsets["stata_dta"] = 0
    tag("stata_dta")
    tag("header")
    out.extend(f"<release>{release}</release>".encode())
    out.extend(b"<byteorder>LSF</byteorder>")
    out.extend(b"<K>" + K.to_bytes(L.k_size, "little") + b"</K>")
    out.extend(b"<N>" + N.to_bytes(L.n_size, "little") + b"</N>")
    label = _enc(ds.label)[:320 if release != 117 else 80]
    out.extend(b"<label>" + len(label).to_bytes(L.label_len_size, "little") + label + b"</label>")
    ts_text = timestamp if timestamp is not None else _timestamp()
    ts = ts_text.encode("ascii")
    out.extend(b"<timestamp>" + bytes([len(ts)]) + ts + b"</timestamp>")
    end("header")

    offsets["map"] = len(out)
    tag("map")
    map_pos = len(out)
    out.extend(b"\x00" * 8 * 14)          # preenchido no fim
    end("map")

    offsets["variable_types"] = len(out)
    tag("variable_types")
    for v in ds.vars:
        if v.vtype == "strL":
            code = T_STRL
        elif v.is_string:
            code = int(v.vtype[3:])
        else:
            code = _TYPE_CODE[v.vtype]
        out.extend(code.to_bytes(2, "little"))
    end("variable_types")

    offsets["varnames"] = len(out)
    tag("varnames")
    for v in ds.vars:
        out.extend(_fixed(v.name, L.name))
    end("varnames")

    offsets["sortlist"] = len(out)
    tag("sortlist")
    names = ds.names
    for k in range(K + 1):
        idx = names.index(ds.sortlist[k]) + 1 if k < len(ds.sortlist) else 0
        out.extend(idx.to_bytes(L.sort_size, "little"))
    end("sortlist")

    offsets["formats"] = len(out)
    tag("formats")
    for v in ds.vars:
        out.extend(_fixed(v.fmt, L.fmt))
    end("formats")

    offsets["value_label_names"] = len(out)
    tag("value_label_names")
    for v in ds.vars:
        out.extend(_fixed(v.value_label, L.name))
    end("value_label_names")

    offsets["variable_labels"] = len(out)
    tag("variable_labels")
    for v in ds.vars:
        out.extend(_fixed(v.label, L.varlabel))
    end("variable_labels")

    offsets["characteristics"] = len(out)
    tag("characteristics")
    for owner, chars in getattr(ds, "chars", {}).items():
        for cname, content in chars.items():
            body = _fixed(owner, L.name) + _fixed(cname, L.name) + _enc(content) + b"\x00"
            out.extend(b"<ch>" + len(body).to_bytes(4, "little") + body + b"</ch>")
    end("characteristics")

    # dados: registros de largura fixa, montados coluna a coluna
    offsets["data"] = len(out)
    tag("data")
    strls: list[tuple[int, int, bytes]] = []
    widths = []
    for v in ds.vars:
        if v.vtype == "strL":
            widths.append(8)
        elif v.is_string:
            widths.append(int(v.vtype[3:]))
        else:
            widths.append(_SIZE[v.vtype])
    rec = sum(widths)
    block = np.zeros((N, rec), dtype=np.uint8)
    pos = 0
    for j, (v, w) in enumerate(zip(ds.vars, widths)):
        if v.vtype == "strL":
            for i, text in enumerate(v.data):
                if text == "":
                    continue                       # (0,0) indica string vazia
                vv, oo = j + 1, i + 1
                if release == 117:
                    ref = vv.to_bytes(4, "little") + oo.to_bytes(4, "little")
                else:
                    ref = vv.to_bytes(2, "little") + oo.to_bytes(6, "little")
                block[i, pos:pos + 8] = np.frombuffer(ref, dtype=np.uint8)
                strls.append((vv, oo, _enc(str(text))))
        elif v.is_string:
            for i, text in enumerate(v.data):
                b = _enc(str(text))[:w]
                if b:
                    block[i, pos:pos + len(b)] = np.frombuffer(b, dtype=np.uint8)
        else:
            disk = _to_disk(v.data, v.vtype)
            block[:, pos:pos + w] = disk.view(np.uint8).reshape(N, w) if N else block[:, pos:pos + w]
        pos += w
    out.extend(block.tobytes())
    end("data")

    offsets["strls"] = len(out)
    tag("strls")
    for vv, oo, payload in strls:
        o_bytes = oo.to_bytes(4 if release == 117 else 8, "little")
        out.extend(b"GSO" + vv.to_bytes(4, "little") + o_bytes + bytes([130])
                   + (len(payload) + 1).to_bytes(4, "little") + payload + b"\x00")
    end("strls")

    offsets["value_labels"] = len(out)
    tag("value_labels")
    for lname, mapping in ds.value_labels.items():
        items = sorted(mapping.items())
        txt = bytearray()
        offs, vals = [], []
        for val, text in items:
            offs.append(len(txt))
            vals.append(int(val))
            txt.extend(_enc(text) + b"\x00")
        n = len(items)
        table = (n.to_bytes(4, "little") + len(txt).to_bytes(4, "little")
                 + b"".join(o.to_bytes(4, "little") for o in offs)
                 + b"".join(struct.pack("<i", val) for val in vals) + bytes(txt))
        out.extend(b"<lbl>" + len(table).to_bytes(4, "little") + _fixed(lname, L.name)
                   + b"\x00" * 3 + table + b"</lbl>")
    end("value_labels")
    offsets["stata_dta_end"] = len(out)
    end("stata_dta")
    offsets["eof"] = len(out)

    order = ["stata_dta", "map", "variable_types", "varnames", "sortlist", "formats",
             "value_label_names", "variable_labels", "characteristics", "data", "strls",
             "value_labels", "stata_dta_end", "eof"]
    out[map_pos:map_pos + 8 * 14] = b"".join(offsets[k].to_bytes(8, "little") for k in order)
    Path(path).write_bytes(bytes(out))
    return ts_text


# ---------------------------------------------------------------------------
# Leitura
# ---------------------------------------------------------------------------

class _Reader:
    def __init__(self, data: bytes):
        self.b = data
        self.p = 0
        self.bo = "<"

    def expect(self, text: str) -> None:
        t = text.encode()
        if self.b[self.p:self.p + len(t)] != t:
            raise StataError(610, "file not Stata format")
        self.p += len(t)

    def take(self, n: int) -> bytes:
        chunk = self.b[self.p:self.p + n]
        self.p += n
        return chunk

    def uint(self, n: int) -> int:
        return int.from_bytes(self.take(n), "little" if self.bo == "<" else "big")

    def until(self, text: str) -> bytes:
        t = text.encode()
        j = self.b.find(t, self.p)
        if j == -1:
            raise StataError(610, "file not Stata format")
        chunk = self.b[self.p:j]
        self.p = j
        return chunk


def _cstr(b: bytes, release: int = 118) -> str:
    """Texto terminado em zero. O Stata 13 (formato 117) e anteriores gravam
    Latin-1, mas outros programas gravam UTF-8 nesses formatos: tenta UTF-8
    primeiro e cai para Latin-1."""
    j = b.find(b"\x00")
    raw = b if j == -1 else b[:j]
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def _fix_text(text: str) -> str:
    """Desfaz UTF-8 lido como Latin-1 ('joÃ£o' -> 'joão'), caso da ReadStat
    com arquivos antigos gravados em UTF-8."""
    try:
        fixed = text.encode("latin-1").decode("utf-8")
    except (UnicodeEncodeError, UnicodeDecodeError):
        return text
    return fixed


def read_dta(path: str | Path) -> Dataset:
    p = Path(path)
    if not p.exists():
        raise StataError(601, f"file {p} not found")
    data = p.read_bytes()
    if data.startswith(b"<stata_dta>"):
        return _read_117plus(data, p)
    return _read_legacy(p)


def _read_117plus(data: bytes, path: Path) -> Dataset:
    r = _Reader(data)
    r.expect("<stata_dta><header><release>")
    release = int(r.until("</release>"))
    if release not in (117, 118, 119):
        return _read_legacy(path)
    L = _Layout(release)
    r.expect("</release><byteorder>")
    r.bo = "<" if r.take(3) == b"LSF" else ">"
    r.expect("</byteorder><K>")
    K = r.uint(L.k_size)
    r.expect("</K><N>")
    N = r.uint(L.n_size)
    r.expect("</N><label>")
    label = _cstr(r.take(r.uint(L.label_len_size)), release)
    r.expect("</label><timestamp>")
    ts = r.take(r.uint(1)).decode("ascii", "replace")
    r.expect("</timestamp></header><map>")
    offsets = [r.uint(8) for _ in range(14)]

    def seek(k: int, name: str) -> None:
        r.p = offsets[k]
        r.expect(f"<{name}>")

    seek(2, "variable_types")
    codes = [r.uint(2) for _ in range(K)]
    seek(3, "varnames")
    names = [_cstr(r.take(L.name), release) for _ in range(K)]
    seek(4, "sortlist")
    sort_idx = [r.uint(L.sort_size) for _ in range(K + 1)]
    seek(5, "formats")
    fmts = [_cstr(r.take(L.fmt), release) for _ in range(K)]
    seek(6, "value_label_names")
    vlabels = [_cstr(r.take(L.name), release) for _ in range(K)]
    seek(7, "variable_labels")
    varlabels = [_cstr(r.take(L.varlabel), release) for _ in range(K)]

    seek(8, "characteristics")
    chars: dict[str, dict[str, str]] = {}
    while r.b.startswith(b"<ch>", r.p):
        r.p += 4
        n = r.uint(4)
        body = r.take(n)
        owner = _cstr(body[:L.name], release)
        cname = _cstr(body[L.name:2 * L.name], release)
        chars.setdefault(owner, {})[cname] = _cstr(body[2 * L.name:], release)
        r.expect("</ch>")

    # tipos e larguras
    vtypes, widths = [], []
    for c in codes:
        if c == T_STRL:
            vtypes.append("strL")
            widths.append(8)
        elif 1 <= c <= 2045:
            vtypes.append(f"str{c}")
            widths.append(c)
        elif c in _CODE_TYPE:
            vtypes.append(_CODE_TYPE[c])
            widths.append(_SIZE[_CODE_TYPE[c]])
        else:
            raise StataError(610, f"unknown variable type {c}")

    seek(9, "data")
    rec = sum(widths)
    block = np.frombuffer(r.take(rec * N), dtype=np.uint8).reshape(N, rec) if N else np.zeros((0, rec), np.uint8)

    # strLs
    seek(10, "strls")
    gso: dict[tuple[int, int], str] = {}
    while r.b.startswith(b"GSO", r.p):
        r.p += 3
        v = r.uint(4)
        o = r.uint(4 if release == 117 else 8)
        t = r.uint(1)
        n = r.uint(4)
        payload = r.take(n)
        if t == 130 and payload.endswith(b"\x00"):
            payload = payload[:-1]
        gso[(v, o)] = _cstr(payload + b"\x00", release)

    # rótulos de valor
    seek(11, "value_labels")
    value_labels: dict[str, dict[int, str]] = {}
    while r.b.startswith(b"<lbl>", r.p):
        r.p += 5
        n_table = r.uint(4)
        lname = _cstr(r.take(L.name), release)
        r.take(3)
        table = r.take(n_table)
        n = int.from_bytes(table[0:4], "little")
        txtlen = int.from_bytes(table[4:8], "little")
        offs = np.frombuffer(table[8:8 + 4 * n], dtype="<i4")
        vals = np.frombuffer(table[8 + 4 * n:8 + 8 * n], dtype="<i4")
        txt = table[8 + 8 * n:8 + 8 * n + txtlen]
        value_labels[lname] = {int(v): _cstr(txt[int(o):], release) for o, v in zip(offs, vals)}
        r.expect("</lbl>")

    ds = Dataset()
    ds.nobs = N
    pos = 0
    for j, (name, vtype, w) in enumerate(zip(names, vtypes, widths)):
        col = block[:, pos:pos + w]
        if vtype == "strL":
            refs = col.copy()
            vals = []
            for i in range(N):
                rb = refs[i].tobytes()
                if release == 117:
                    vv, oo = int.from_bytes(rb[:4], "little"), int.from_bytes(rb[4:], "little")
                else:
                    vv, oo = int.from_bytes(rb[:2], "little"), int.from_bytes(rb[2:], "little")
                vals.append(gso.get((vv, oo), "") if vv else "")
            arr = np.array(vals, dtype=object)
        elif vtype.startswith("str"):
            arr = np.array([_cstr(col[i].tobytes(), release) for i in range(N)], dtype=object)
        else:
            raw = np.ascontiguousarray(col).view(r.bo + _NP[vtype]).reshape(N)
            arr = _from_disk(raw.astype("<" + _NP[vtype]), vtype)
        pos += w
        ds.vars.append(Variable(name, vtype, arr, fmt=fmts[j] or default_format(vtype),
                                label=varlabels[j], value_label=vlabels[j]))
    ds.sortlist = [names[k - 1] for k in sort_idx if k > 0 and k <= K]
    ds.label = label
    ds.timestamp = ts
    ds.value_labels = value_labels
    ds.chars = chars
    ds.filename = str(path)
    ds.fullpath = str(Path(path).resolve())
    ds.changed = False
    return ds


def _read_legacy(path: Path) -> Dataset:
    """Formatos até 115 (Stata 12 e anteriores), lidos pela ReadStat."""
    try:
        import pyreadstat
    except ImportError:
        raise StataError(610, "file format not supported: install pyreadstat to read older .dta files")
    try:
        # saída em dicionário de listas: dispensa o pandas
        cols, meta = pyreadstat.read_dta(str(path), user_missing=True, output_format="dict")
    except Exception as e:  # noqa: BLE001
        raise StataError(610, f"file not Stata format ({e})")
    rs_types = {"int8": "byte", "int16": "int", "int32": "long", "float": "float", "double": "double"}
    ds = Dataset()
    ds.nobs = len(next(iter(cols.values()), []))
    for name, column in cols.items():
        rt = meta.readstat_variable_types.get(name, "double")
        fmt = meta.original_variable_types.get(name, "")
        if rt == "string":
            width = meta.variable_storage_width.get(name, 1) or 1
            vtype = f"str{width}"
            arr = np.array(["" if (x is None or x != x) else _fix_text(str(x)) for x in column], dtype=object)
        else:
            vtype = rs_types.get(rt, "double")
            vals = []
            for x in column:
                if isinstance(x, str) and len(x) == 1 and x.isalpha():
                    vals.append(M.EXTENDED[x.lower()])
                elif x is None or x != x:
                    vals.append(M.SYSMISS)
                else:
                    vals.append(float(x))
            arr = np.array(vals, dtype=np.float64)
        ds.vars.append(Variable(name, vtype, arr, fmt=fmt or default_format(vtype),
                                label=_fix_text(meta.column_names_to_labels.get(name) or ""),
                                value_label=meta.variable_to_label.get(name, "") or ""))
    ds.value_labels = {k: {int(a): _fix_text(b) for a, b in v.items()}
                       for k, v in (meta.value_labels or {}).items()}
    ds.label = _fix_text(meta.file_label or "")
    ds.timestamp = ""
    ds.chars = {}
    ds.filename = str(path)
    ds.fullpath = str(Path(path).resolve())
    ds.changed = False
    return ds
