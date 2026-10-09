"""Leitura e gravação de arquivos .dta.

Implementação própria dos formatos 113 a 115 (Stata 8 a 12, só leitura em
113), 117 (Stata 13), 118 (Stata 14/15) e 119 (Stata 15+, mais de 32.767
variáveis), escrita a partir da especificação pública do formato. Formatos
anteriores a 113 são lidos pela ReadStat (pyreadstat), que também serve,
nos testes, de verificação independente de tudo o que este módulo lê e
grava.

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
import os
import struct
from pathlib import Path

import numpy as np

from ..core import missing as M
from ..core.dataset import Dataset, Variable, default_format
from ..core.errors import StataError

# códigos de tipo (117+)
T_DOUBLE, T_FLOAT, T_LONG, T_INT, T_BYTE = 65526, 65527, 65528, 65529, 65530
T_STRL = 32768
_TYPE_CODE = {"double": T_DOUBLE, "float": T_FLOAT, "long": T_LONG, "int": T_INT, "byte": T_BYTE}
_CODE_TYPE = {v: k for k, v in _TYPE_CODE.items()}
_NP = {"double": "f8", "float": "f4", "long": "i4", "int": "i2", "byte": "i1"}
# os vetores nativos das variáveis (core/storage.py) já usam a codificação de
# missing do formato: ler e gravar é copiar bytes
_SIZE = {"double": 8, "float": 4, "long": 4, "int": 2, "byte": 1}


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
    return f"{now.day:2d} {_MONTHS[now.month - 1]} {now.year} {now.hour:02d}:{now.minute:02d}"


CHUNK_BYTES = 32 * 2**20      # dados gravados e lidos em blocos de ~32 MB
_CACHE_MAX = 50_000           # valores distintos por variável string antes de desistir do cache


def _widths(ds: Dataset) -> list[int]:
    out = []
    for v in ds.vars:
        if v.vtype == "strL":
            out.append(8)
        elif v.is_string:
            out.append(int(v.vtype[3:]))
        else:
            out.append(_SIZE[v.vtype])
    return out


def _fill_block(vars_, widths: list[int], rec: int, start: int, stop: int, release: int,
                strls: list[tuple[int, int, bytes]]) -> np.ndarray:
    """Registros das observações start..stop-1, prontos para o disco."""
    n = stop - start
    block = np.zeros((n, rec), dtype=np.uint8)
    pos = 0
    for j, (v, w) in enumerate(zip(vars_, widths)):
        part = v.raw[start:stop]
        if v.vtype == "strL":
            refs = np.zeros((n, 8), dtype=np.uint8)
            for i, text in enumerate(part):
                if text == "":
                    continue                   # (0,0) indica string vazia
                vv, oo = j + 1, start + i + 1
                if release == 117:
                    ref = vv.to_bytes(4, "little") + oo.to_bytes(4, "little")
                else:
                    ref = vv.to_bytes(2, "little") + oo.to_bytes(6, "little")
                refs[i] = np.frombuffer(ref, dtype=np.uint8)
                strls.append((vv, oo, _enc(str(text))))
            block[:, pos:pos + 8] = refs
        elif v.is_string:
            # dtype S{w}: completa com zeros, como no formato
            enc = np.array([_enc(str(x))[:w] for x in part], dtype=f"S{w}")
            block[:, pos:pos + w] = enc.view(np.uint8).reshape(n, w)
        else:
            disk = np.ascontiguousarray(part, dtype="<" + _NP[v.vtype])
            block[:, pos:pos + w] = disk.view(np.uint8).reshape(n, w)
        pos += w
    return block


# ---------------------------------------------------------------------------
# Formatos 114 e 115 (Stata 10 a 12): cabeçalho binário de tamanho fixo
# ---------------------------------------------------------------------------

_OLD_TYPE = {"byte": 251, "int": 252, "long": 253, "float": 254, "double": 255}
_OLD_TYPE_BACK = {v: k for k, v in _OLD_TYPE.items()}
OLD_STR_MAX = 244


def _old_vars(ds: Dataset) -> tuple[list, list[str]]:
    """Variáveis prontas para 114/115: strL e str acima de 244 viram str244
    (cortando o texto). Devolve (variáveis, nomes das que foram cortadas)."""
    out, cut = [], []
    for v in ds.vars:
        if v.is_string and (v.vtype == "strL" or int(v.vtype[3:]) > OLD_STR_MAX):
            texts = [_enc(str(x))[:OLD_STR_MAX].decode("latin-1") for x in v.raw]
            nv = Variable(v.name, f"str{OLD_STR_MAX}", texts, fmt=v.fmt, label=v.label,
                          value_label=v.value_label)
            if any(len(_enc(str(x))) > OLD_STR_MAX for x in v.raw):
                cut.append(v.name)
            out.append(nv)
        else:
            out.append(v)
    return out, cut


def _write_old_body(ds: Dataset, tmp: Path, release: int, ts_text: str) -> None:
    vars_, _ = _old_vars(ds)
    K, N = len(vars_), ds.nobs
    if K > 32767:
        raise StataError(900, "too many variables for this file format")
    with open(tmp, "wb") as f:
        out = bytearray()
        out += bytes([release, 2, 1, 0])               # versão, LOHI, filetype, livre
        out += K.to_bytes(2, "little") + N.to_bytes(4, "little")
        out += _fixed(ds.label, 81)
        out += _fixed(ts_text, 18)
        for v in vars_:
            out.append(int(v.vtype[3:]) if v.is_string else _OLD_TYPE[v.vtype])
        for v in vars_:
            out += _fixed(v.name, 33)
        pos_of = {n: k for k, n in enumerate([v.name for v in vars_])}
        for k in range(K + 1):
            idx = pos_of[ds.sortlist[k]] + 1 if k < len(ds.sortlist) else 0
            out += idx.to_bytes(2, "little")
        for v in vars_:
            out += _fixed(v.fmt, 49)
        for v in vars_:
            out += _fixed(v.value_label, 33)
        for v in vars_:
            out += _fixed(v.label, 81)
        # campos de expansão: características
        for owner, chars in getattr(ds, "chars", {}).items():
            for cname, content in chars.items():
                body = _fixed(owner, 33) + _fixed(cname, 33) + _enc(content) + b"\x00"
                out += bytes([1]) + len(body).to_bytes(4, "little") + body
        out += b"\x00" * 5
        f.write(out)

        widths = [int(v.vtype[3:]) if v.is_string else _SIZE[v.vtype] for v in vars_]
        rec = sum(widths)
        step = max(1, CHUNK_BYTES // max(rec, 1))
        for start in range(0, N, step):
            block = _fill_block(vars_, widths, rec, start, min(N, start + step), release, [])
            f.write(block.data)
            del block

        out = bytearray()
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
            out += len(table).to_bytes(4, "little") + _fixed(lname, 33) + b"\x00" * 3 + table
        f.write(out)


def _read_old(path: Path) -> Dataset:
    """Formatos 113, 114 e 115 (Stata 8 a 12)."""
    with open(path, "rb") as f:
        head = f.read(109)
        release, bo_code = head[0], head[1]
        bo = ">" if bo_code == 1 else "<"
        big = bo == ">"

        def u(b: bytes) -> int:
            return int.from_bytes(b, "big" if big else "little")

        K = u(head[4:6])
        N = u(head[6:10])
        label = _cstr(head[10:91], 117)
        ts = _cstr(head[91:109], 117)
        fmt_len = 12 if release == 113 else 49
        meta_len = K + 33 * K + 2 * (K + 1) + fmt_len * K + 33 * K + 81 * K
        meta = f.read(meta_len)
        p = 0
        codes = list(meta[p:p + K]); p += K
        names = [_cstr(meta[p + 33 * i:p + 33 * (i + 1)], 117) for i in range(K)]; p += 33 * K
        sort_idx = [u(meta[p + 2 * i:p + 2 * i + 2]) for i in range(K + 1)]; p += 2 * (K + 1)
        fmts = [_cstr(meta[p + fmt_len * i:p + fmt_len * (i + 1)], 117) for i in range(K)]
        p += fmt_len * K
        vlabels = [_cstr(meta[p + 33 * i:p + 33 * (i + 1)], 117) for i in range(K)]; p += 33 * K
        varlabels = [_cstr(meta[p + 81 * i:p + 81 * (i + 1)], 117) for i in range(K)]

        chars: dict[str, dict[str, str]] = {}
        while True:
            hdr = f.read(5)
            if len(hdr) < 5:
                raise StataError(610, "file not Stata format")
            dtype, n = hdr[0], u(hdr[1:5])
            if dtype == 0 and n == 0:
                break
            body = f.read(n)
            if dtype == 1 and n >= 66:
                owner, cname = _cstr(body[:33], 117), _cstr(body[33:66], 117)
                chars.setdefault(owner, {})[cname] = _cstr(body[66:], 117)

        vtypes, widths = [], []
        for c in codes:
            if c in _OLD_TYPE_BACK:
                vtypes.append(_OLD_TYPE_BACK[c])
                widths.append(_SIZE[_OLD_TYPE_BACK[c]])
            elif 1 <= c <= OLD_STR_MAX:
                vtypes.append(f"str{c}")
                widths.append(c)
            else:
                raise StataError(610, f"unknown variable type {c}")
        rec = sum(widths)
        cols: list = [np.empty(N, dtype=object) if t.startswith("str") else np.empty(N, dtype=_NP[t])
                      for t in vtypes]
        step = max(1, CHUNK_BYTES // max(rec, 1))
        for start in range(0, N, step):
            stop = min(N, start + step)
            n = stop - start
            buf = f.read(n * rec)
            if len(buf) != n * rec:
                raise StataError(610, "file not Stata format (data section truncated)")
            block = np.frombuffer(buf, dtype=np.uint8).reshape(n, rec)
            pos = 0
            for j, (vtype, w) in enumerate(zip(vtypes, widths)):
                part = block[:, pos:pos + w]
                if vtype.startswith("str"):
                    fixed = np.ascontiguousarray(part).view(f"S{w}").reshape(n).tolist()
                    cols[j][start:stop] = [_cstr(x, 117) for x in fixed]
                else:
                    cols[j][start:stop] = np.ascontiguousarray(part).view(bo + _NP[vtype]).reshape(n)
                pos += w

        value_labels: dict[str, dict[int, str]] = {}
        while True:
            hdr = f.read(4)
            if len(hdr) < 4:
                break
            n_table = u(hdr)
            lname = _cstr(f.read(33), 117)
            f.read(3)
            table = f.read(n_table)
            n = u(table[0:4])
            txtlen = u(table[4:8])
            offs = np.frombuffer(table[8:8 + 4 * n], dtype=bo + "i4")
            vals = np.frombuffer(table[8 + 4 * n:8 + 8 * n], dtype=bo + "i4")
            txt = table[8 + 8 * n:8 + 8 * n + txtlen]
            value_labels[lname] = {int(v): _cstr(txt[int(o):], 117) for o, v in zip(offs, vals)}

    ds = Dataset()
    ds.nobs = N
    for j, (name, vtype) in enumerate(zip(names, vtypes)):
        if vtype.startswith("str"):
            var = Variable(name, vtype, cols[j])
        else:
            var = Variable(name, vtype, np.empty(0))
            var.raw = cols[j]
        var.fmt = fmts[j] or default_format(vtype)
        var.label = varlabels[j]
        var.value_label = vlabels[j]
        ds.vars.append(var)
    ds.sortlist = [names[k - 1] for k in sort_idx if 0 < k <= K]
    ds.label = label
    ds.timestamp = ts
    ds.value_labels = value_labels
    ds.chars = chars
    ds.filename = str(path)
    ds.fullpath = str(Path(path).resolve())
    ds.changed = False
    return ds


def write_dta(ds: Dataset, path: str | Path, *, release: int = 118,
              timestamp: str | None = None) -> str:
    """Grava ds em path e devolve o carimbo de data gravado no cabeçalho.

    Os dados vão para o disco em blocos de linhas, sem montar o arquivo
    inteiro na memória. A gravação é feita num arquivo temporário na mesma
    pasta, que só substitui o destino no fim: uma falha no meio não destrói
    o arquivo antigo."""
    if release not in (114, 115, 117, 118, 119):
        raise StataError(198, f"dta release {release} not supported")
    if release == 118 and ds.nvars > 32767:
        release = 119
    K, N = ds.nvars, ds.nobs
    path = Path(path)
    tmp = path.with_name(path.name + ".opendta-tmp")
    offsets: dict[str, int] = {}
    ts_text = timestamp if timestamp is not None else _timestamp()

    try:
        if release in (114, 115):
            _ENC["enc"] = "latin-1"
            _write_old_body(ds, tmp, release, ts_text)
        else:
            L = _Layout(release)
            _ENC["enc"] = "latin-1" if release == 117 else "utf-8"
            _write_body(ds, tmp, release, L, K, N, ts_text, offsets)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, path)
    return ts_text


def _write_body(ds: Dataset, tmp: Path, release: int, L: "_Layout", K: int, N: int,
                ts_text: str, offsets: dict[str, int]) -> None:
    with open(tmp, "wb") as f:
        out = bytearray()

        def flush() -> None:
            f.write(out)
            out.clear()

        def here() -> int:
            return f.tell() + len(out)

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
        ts = ts_text.encode("ascii")
        out.extend(b"<timestamp>" + bytes([len(ts)]) + ts + b"</timestamp>")
        end("header")

        offsets["map"] = here()
        tag("map")
        map_pos = here()
        out.extend(b"\x00" * 8 * 14)          # preenchido no fim
        end("map")

        offsets["variable_types"] = here()
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

        offsets["varnames"] = here()
        tag("varnames")
        for v in ds.vars:
            out.extend(_fixed(v.name, L.name))
        end("varnames")

        offsets["sortlist"] = here()
        tag("sortlist")
        pos_of = {n: k for k, n in enumerate(ds.names, start=1)}
        for k in range(K + 1):
            idx = pos_of[ds.sortlist[k]] if k < len(ds.sortlist) else 0
            out.extend(idx.to_bytes(L.sort_size, "little"))
        end("sortlist")

        offsets["formats"] = here()
        tag("formats")
        for v in ds.vars:
            out.extend(_fixed(v.fmt, L.fmt))
        end("formats")

        offsets["value_label_names"] = here()
        tag("value_label_names")
        for v in ds.vars:
            out.extend(_fixed(v.value_label, L.name))
        end("value_label_names")

        offsets["variable_labels"] = here()
        tag("variable_labels")
        for v in ds.vars:
            out.extend(_fixed(v.label, L.varlabel))
        end("variable_labels")

        offsets["characteristics"] = here()
        tag("characteristics")
        for owner, chars in getattr(ds, "chars", {}).items():
            for cname, content in chars.items():
                body = _fixed(owner, L.name) + _fixed(cname, L.name) + _enc(content) + b"\x00"
                out.extend(b"<ch>" + len(body).to_bytes(4, "little") + body + b"</ch>")
        end("characteristics")

        # dados: registros de largura fixa, gravados em blocos de linhas
        offsets["data"] = here()
        tag("data")
        flush()
        widths = _widths(ds)
        rec = sum(widths)
        strls: list[tuple[int, int, bytes]] = []
        step = max(1, CHUNK_BYTES // max(rec, 1))
        for start in range(0, N, step):
            stop = min(N, start + step)
            block = _fill_block(ds.vars, widths, rec, start, stop, release, strls)
            f.write(block.data)
            del block
        end("data")

        offsets["strls"] = here()
        tag("strls")
        for vv, oo, payload in strls:
            o_bytes = oo.to_bytes(4 if release == 117 else 8, "little")
            out.extend(b"GSO" + vv.to_bytes(4, "little") + o_bytes + bytes([130])
                       + (len(payload) + 1).to_bytes(4, "little") + payload + b"\x00")
            if len(out) > CHUNK_BYTES:
                flush()
        end("strls")

        offsets["value_labels"] = here()
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
        offsets["stata_dta_end"] = here()
        end("stata_dta")
        offsets["eof"] = here()
        flush()

        order = ["stata_dta", "map", "variable_types", "varnames", "sortlist", "formats",
                 "value_label_names", "variable_labels", "characteristics", "data", "strls",
                 "value_labels", "stata_dta_end", "eof"]
        f.seek(map_pos)
        f.write(b"".join(offsets[k].to_bytes(8, "little") for k in order))


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
    with open(p, "rb") as f:
        head = f.read(11)
        if head == b"<stata_dta>":
            return _read_117plus(f, p)
    if len(head) >= 4 and head[0] in (113, 114, 115) and head[1] in (1, 2) and head[2] == 1:
        return _read_old(p)
    return _read_legacy(p)


def _read_117plus(f, path: Path) -> Dataset:
    """Lê os metadados de uma vez, os dados em blocos de linhas direto para
    os vetores nativos de cada variável, e o fim do arquivo (strLs e rótulos)
    de uma vez."""
    f.seek(0)
    head = f.read(4096)
    r = _Reader(head)
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
    bo = r.bo

    # metadados: tudo antes de <data>
    f.seek(0)
    r = _Reader(f.read(offsets[9]))
    r.bo = bo

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

    # dados, em blocos
    f.seek(offsets[9])
    if f.read(6) != b"<data>":
        raise StataError(610, "file not Stata format")
    rec = sum(widths)
    cols: list[np.ndarray] = []
    for vtype in vtypes:
        if vtype == "strL":
            cols.append(np.zeros((N, 8), dtype=np.uint8))
        elif vtype.startswith("str"):
            cols.append(np.empty(N, dtype=object))
        else:
            cols.append(np.empty(N, dtype=_NP[vtype]))
    # strings repetidas (categorias, municípios...) viram um único objeto str:
    # 8 bytes por célula em vez de ~60. Desliga para variáveis quase únicas.
    caches: list[dict | None] = [{} for _ in vtypes]
    step = max(1, CHUNK_BYTES // max(rec, 1))
    for start in range(0, N, step):
        stop = min(N, start + step)
        n = stop - start
        buf = f.read(n * rec)
        if len(buf) != n * rec:
            raise StataError(610, "file not Stata format (data section truncated)")
        block = np.frombuffer(buf, dtype=np.uint8).reshape(n, rec)
        pos = 0
        for j, (vtype, w) in enumerate(zip(vtypes, widths)):
            part = block[:, pos:pos + w]
            if vtype == "strL":
                cols[j][start:stop] = part
            elif vtype.startswith("str"):
                fixed = np.ascontiguousarray(part).view(f"S{w}").reshape(n).tolist()
                cache = caches[j]
                if cache is None:
                    cols[j][start:stop] = [_cstr(x, release) for x in fixed]
                else:
                    out = []
                    for x in fixed:
                        s = cache.get(x)
                        if s is None:
                            s = cache[x] = _cstr(x, release)
                        out.append(s)
                    cols[j][start:stop] = out
                    if len(cache) > _CACHE_MAX:
                        caches[j] = None
            else:
                raw = np.ascontiguousarray(part).view(bo + _NP[vtype]).reshape(n)
                cols[j][start:stop] = raw            # converte a ordem dos bytes se preciso
            pos += w
        del block, buf

    # fim do arquivo: strLs e rótulos de valor
    f.seek(offsets[10])
    r = _Reader(f.read())
    r.bo = bo
    base = offsets[10]

    def seek_tail(k: int, name: str) -> None:
        r.p = offsets[k] - base
        r.expect(f"<{name}>")

    seek_tail(10, "strls")
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

    seek_tail(11, "value_labels")
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
    for j, (name, vtype) in enumerate(zip(names, vtypes)):
        col = cols[j]
        if vtype == "strL":
            if release == 117:
                vv = col[:, :4].copy().view("<u4").reshape(N).astype(np.int64)
                oo = col[:, 4:].copy().view("<u4").reshape(N).astype(np.int64)
            else:
                vv = col[:, :2].copy().view("<u2").reshape(N).astype(np.int64)
                o8 = np.zeros((N, 8), dtype=np.uint8)
                o8[:, :6] = col[:, 2:]
                oo = o8.view("<u8").reshape(N).astype(np.int64)
            arr = np.array([gso.get((a, b), "") if a else "" for a, b in zip(vv.tolist(), oo.tolist())],
                           dtype=object)
            var = Variable(name, vtype, arr)
        elif vtype.startswith("str"):
            var = Variable(name, vtype, col)
        else:
            var = Variable(name, vtype, np.empty(0))
            var.raw = col
        var.fmt = fmts[j] or default_format(vtype)
        var.label = varlabels[j]
        var.value_label = vlabels[j]
        ds.vars.append(var)
        cols[j] = None
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
