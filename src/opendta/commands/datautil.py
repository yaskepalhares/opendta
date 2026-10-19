"""File and operating-system utilities ([D] copy, dir, mkdir, rmdir, shell,
filefilter, zipfile, checksum, snapshot, unicode). Messages not shown in the
manual are marked VERIFY and checked in compat/do/1501_datamgmt.do."""

from __future__ import annotations

import datetime as _dt
import os
import shutil
import subprocess
import zipfile as _zip
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lang.words import strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _split(args: str) -> tuple[list[str], str]:
    """Words before the comma (quotes respected) and the options."""
    import shlex
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    try:
        words = shlex.split(body, posix=True)
    except ValueError:
        raise StataError(198, "invalid syntax")
    return words, opt


def _path(t: str) -> Path:
    return Path(os.path.expanduser(strip_outer_quotes(t)))


# ---------------------------------------------------------------------------
# copy, mkdir, rmdir, dir/ls
# ---------------------------------------------------------------------------

@command("copy", "copy")
def cmd_copy(s: "Session", args: str) -> None:
    words, opt = _split(args)
    o = match_options(opt, {"public": 3, "text": 4, "replace": 7})
    if len(words) != 2:
        raise StataError(198, "invalid syntax")
    src, dst = words
    target = _path(dst)
    if target.is_dir():
        target = target / Path(src).name
    if target.exists() and not o.get("replace"):
        raise StataError(602, f"file {dst} already exists")
    if src.startswith(("http://", "https://")):
        import urllib.request
        from ..updater import ssl_context
        try:
            with urllib.request.urlopen(src, timeout=30, context=ssl_context()) as r:
                data = r.read()
        except Exception as e:   # noqa: BLE001
            raise StataError(631, f"host not found or file {src} not found ({e})")   # VERIFY
        target.write_bytes(data)
        return
    p = _path(src)
    if not p.exists():
        raise StataError(601, f"file {src} not found")
    if o.get("text"):
        data = p.read_bytes().replace(b"\r\n", b"\n")
        target.write_bytes(data.replace(b"\n", os.linesep.encode()))
    else:
        shutil.copyfile(p, target)


@command("mkdir", "mkdir")
def cmd_mkdir(s: "Session", args: str) -> None:
    words, opt = _split(args)
    match_options(opt, {"public": 3})
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    try:
        os.mkdir(_path(words[0]))
    except OSError:
        raise StataError(693, f"could not create directory {words[0]}")


@command("rmdir", "rmdir")
def cmd_rmdir(s: "Session", args: str) -> None:
    words, opt = _split(args)
    match_options(opt, {})
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    p = _path(words[0])
    if not p.is_dir():
        raise StataError(601, f"directory {words[0]} not found")   # VERIFY
    try:
        p.rmdir()
    except OSError:
        raise StataError(693, f"could not remove directory {words[0]}")


def _dir(s: "Session", args: str) -> None:
    words, opt = _split(args)
    o = match_options(opt, {"wide": 1})
    spec = words[0] if words else "*"
    p = _path(spec)
    if p.is_dir():
        files = sorted(p.iterdir())
    else:
        files = sorted(p.parent.glob(p.name)) if any(c in p.name for c in "*?[") else ([p] if p.exists() else [])
    out = s.output
    if o.get("wide"):
        names = [f.name + ("/" if f.is_dir() else "") for f in files]
        width = max((len(n) for n in names), default=0) + 2
        per = max(1, 79 // max(width, 1))
        for i in range(0, len(names), per):
            out.write("".join(f"{n:<{width}}" for n in names[i:i + per]).rstrip() + "\n", "text")
        return
    for f in files:
        st = f.stat()
        when = _dt.datetime.fromtimestamp(st.st_mtime)
        # VERIFY: Stata's Unix layout (size in k, date m/dd/yy, time)
        size = "<dir>" if f.is_dir() else _ksize(st.st_size)
        out.write(f"{size:>8}  {when.month:>2}/{when:%d/%y} {when:%H:%M}  {f.name}\n", "text")


def _ksize(n: int) -> str:
    return f"{n / 1024:.1f}k" if n < 1024 * 1024 else f"{n / 1024 / 1024:.1f}M"


@command("dir", "dir")
def cmd_dir(s: "Session", args: str) -> None:
    _dir(s, args)


@command("ls", "ls")
def cmd_ls(s: "Session", args: str) -> None:
    _dir(s, args)


# ---------------------------------------------------------------------------
# shell (and the ! prefix)
# ---------------------------------------------------------------------------

@command("shell", "sh")
def cmd_shell(s: "Session", args: str) -> None:
    cmd = args.strip()
    if not cmd:
        # Stata opens an interactive shell; here there is no terminal to give (VERIFY)
        raise StataError(198, "shell without a command is not available in OpenDTA")
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        raise StataError(1, "shell command timed out")
    if r.stdout:
        s.output.write(r.stdout if r.stdout.endswith("\n") else r.stdout + "\n", "text")
    if r.stderr:
        s.output.write(r.stderr if r.stderr.endswith("\n") else r.stderr + "\n", "text")


# ---------------------------------------------------------------------------
# filefilter
# ---------------------------------------------------------------------------

_ESC = {"\\BS": b"\\", "\\r": b"\r", "\\n": b"\n", "\\t": b"\t", "\\Q": b'"', "\\LQ": b"`",
        "\\RQ": b"'", "\\$": b"$", "\\W": b"\r\n", "\\M": b"\r", "\\U": b"\n"}


def _pattern(text: str) -> bytes:
    """filefilter pattern: escape codes (\\BS, \\n, \\W...) and %ddd decimal bytes."""
    t = strip_outer_quotes(text.strip())
    out = bytearray()
    i = 0
    while i < len(t):
        if t[i] == "\\":
            for code in sorted(_ESC, key=len, reverse=True):
                if t.startswith(code, i):
                    out += _ESC[code]
                    i += len(code)
                    break
            else:
                out += t[i].encode()
                i += 1
        elif t[i] == "%" and t[i + 1:i + 4].isdigit() and len(t[i + 1:i + 4]) == 3:
            out.append(int(t[i + 1:i + 4]))
            i += 4
        else:
            out += t[i].encode("utf-8")
            i += 1
    return bytes(out)


@command("filefilter", "filefilter")
def cmd_filefilter(s: "Session", args: str) -> None:
    words, opt = _split(args)
    o = match_options(opt, {"from": 4, "to": 2, "replace": 7})
    if len(words) != 2 or not isinstance(o.get("from"), str) or not isinstance(o.get("to"), str):
        raise StataError(198, "invalid syntax")
    src, dst = _path(words[0]), _path(words[1])
    if not src.exists():
        raise StataError(601, f"file {words[0]} not found")
    if src.resolve() == dst.resolve():
        raise StataError(602, "the old and new files must be different")   # VERIFY
    if dst.exists() and not o.get("replace"):
        raise StataError(602, f"file {words[1]} already exists")
    old, new = _pattern(str(o["from"])), _pattern(str(o["to"]))
    if not old:
        raise StataError(198, "from() may not be empty")   # VERIFY
    data = src.read_bytes()
    n = data.count(old)
    dst.write_bytes(data.replace(old, new))
    s.r = {"occurrences": float(n), "bytes_from": float(n * len(old)), "bytes_to": float(n * len(new))}


# ---------------------------------------------------------------------------
# zipfile, unzipfile
# ---------------------------------------------------------------------------

@command("zipfile", "zipfile")
def cmd_zipfile(s: "Session", args: str) -> None:
    words, opt = _split(args)
    o = match_options(opt, {"saving": 2})
    if not words or not isinstance(o.get("saving"), str):
        raise StataError(198, "option saving() required")
    sv = str(o["saving"])
    c = find_top(sv, ",")
    name, rep = (sv, "") if c == -1 else (sv[:c], sv[c + 1:])
    target = _path(name.strip())
    if target.suffix == "":
        target = target.with_suffix(".zip")
    if target.exists() and "replace" not in rep:
        raise StataError(602, f"file {target} already exists")
    files: list[Path] = []
    for w in words:
        p = _path(w)
        hits = sorted(p.parent.glob(p.name)) if any(ch in p.name for ch in "*?[") else [p]
        for h in hits:
            if not h.exists():
                raise StataError(601, f"file {w} not found")
            files += [h] if h.is_file() else [f for f in sorted(h.rglob("*")) if f.is_file()]
    out = s.output
    with _zip.ZipFile(target, "w", _zip.ZIP_DEFLATED) as z:
        for f in files:
            z.write(f, os.path.relpath(f))
            out.write(f"    adding: {os.path.relpath(f)}\n", "text")   # VERIFY
    out.write(f"\nfile {target} saved\n", "text")


@command("unzipfile", "unzipfile")
def cmd_unzipfile(s: "Session", args: str) -> None:
    words, opt = _split(args)
    o = match_options(opt, {"replace": 7})
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    p = _path(words[0])
    if p.suffix == "" and not p.exists():
        p = p.with_suffix(".zip")
    if not p.exists():
        raise StataError(601, f"file {words[0]} not found")
    out = s.output
    with _zip.ZipFile(p) as z:
        for info in z.infolist():
            dest = Path(info.filename)
            if dest.is_absolute() or ".." in dest.parts:
                raise StataError(198, f"unsafe path {info.filename} in archive")
            if dest.exists() and not info.is_dir() and not o.get("replace"):
                raise StataError(602, f"file {info.filename} already exists")
            z.extract(info)
            if not info.is_dir():
                out.write(f"    inflating: {info.filename}\n", "text")   # VERIFY
    out.write(f"\nsuccessfully unzipped {p} to current directory\n", "text")


# ---------------------------------------------------------------------------
# checksum
# ---------------------------------------------------------------------------

def cksum(data: bytes) -> int:
    """POSIX cksum CRC (polynomial 0x04C11DB7, length appended). VERIFY: [D] checksum does not
    name its algorithm."""
    crc = 0
    table = _CK_TABLE
    for b in data:
        crc = ((crc << 8) & 0xFFFFFFFF) ^ table[((crc >> 24) ^ b) & 0xFF]
    n = len(data)
    while n:
        crc = ((crc << 8) & 0xFFFFFFFF) ^ table[((crc >> 24) ^ (n & 0xFF)) & 0xFF]
        n >>= 8
    return (~crc) & 0xFFFFFFFF


def _ck_table() -> list[int]:
    t = []
    for i in range(256):
        c = i << 24
        for _ in range(8):
            c = ((c << 1) ^ 0x04C11DB7) if c & 0x80000000 else (c << 1)
        t.append(c & 0xFFFFFFFF)
    return t


_CK_TABLE = _ck_table()


@command("checksum", "checksum")
def cmd_checksum(s: "Session", args: str) -> None:
    words, opt = _split(args)
    o = match_options(opt, {"save": 4, "replace": 7, "saving": 6})
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    p = _path(words[0])
    if not p.exists():
        raise StataError(601, f"file {words[0]} not found")
    data = p.read_bytes()
    ck = cksum(data)
    s.output.write(f"Checksum for {words[0]} = {ck}, size = {len(data)}\n", "text")
    if o.get("save") or isinstance(o.get("saving"), str):
        target = _path(str(o["saving"])) if isinstance(o.get("saving"), str) else p.with_suffix(".sum")
        if target.exists() and not o.get("replace"):
            raise StataError(602, f"file {target} already exists")
        target.write_text(f"1 {len(data)} {ck}\n", encoding="ascii")   # VERIFY: .sum layout
    s.r = {"version": 1.0, "filelen": float(len(data)), "checksum": float(ck)}


# ---------------------------------------------------------------------------
# snapshot
# ---------------------------------------------------------------------------

def _snaps(s: "Session") -> dict:
    if not hasattr(s, "_snapshots"):
        s._snapshots = {}
    return s._snapshots


def _stamp() -> str:
    return _dt.datetime.now().strftime("%d %b %Y %H:%M")


@command("snapshot", "snapshot")
def cmd_snapshot(s: "Session", args: str) -> None:
    words, opt = _split(args)
    if not words:
        raise StataError(198, "invalid syntax")
    sub, rest = words[0], words[1:]
    snaps = _snaps(s)
    out = s.output
    if sub == "save":
        o = match_options(opt, {"label": 1})
        k = max(snaps, default=0) + 1
        label = strip_outer_quotes(str(o["label"])) if isinstance(o.get("label"), str) else ""
        snaps[k] = {"data": s.data.copy(), "label": label, "when": _stamp()}
        out.write(f"snapshot {k}{' (' + label + ')' if label else ''} created at {snaps[k]['when']}\n", "text")
        s.r = {"snapshot": float(k)}
        return
    if sub == "label":
        if len(rest) != 2:
            raise StataError(198, "invalid syntax")
        k = _snap_no(snaps, rest[0])
        snaps[k]["label"] = rest[1]
        return
    if sub == "restore":
        if len(rest) != 1:
            raise StataError(198, "invalid syntax")
        k = _snap_no(snaps, rest[0])
        s.data = snaps[k]["data"].copy()
        s.notify_state()
        return
    if sub in ("list", "erase"):
        if not rest or rest == ["_all"]:
            keys = sorted(snaps)
        else:
            keys = [_snap_no(snaps, r) for r in _numlist(rest)]
        for k in keys:
            if sub == "erase":
                snaps.pop(k, None)
            else:
                lab = snaps[k]["label"]
                out.write(f"snapshot {k}{' (' + lab + ')' if lab else ''} created at {snaps[k]['when']}\n", "text")
        return
    raise StataError(198, f"snapshot {sub}: unknown subcommand")   # VERIFY


def _numlist(words: list[str]) -> list[str]:
    out = []
    for w in words:
        if "/" in w:
            a, b = w.split("/")
            out += [str(i) for i in range(int(a), int(b) + 1)]
        else:
            out.append(w)
    return out


def _snap_no(snaps: dict, text: str) -> int:
    try:
        k = int(text)
    except ValueError:
        raise StataError(198, "invalid snapshot number")
    if k not in snaps:
        raise StataError(198, f"snapshot {k} not found")   # VERIFY
    return k


# ---------------------------------------------------------------------------
# unicode
# ---------------------------------------------------------------------------

def _is_utf8(data: bytes) -> bool:
    try:
        data.decode("utf-8")
        return True
    except UnicodeDecodeError:
        return False


_TEXT_EXT = {".do", ".ado", ".txt", ".csv", ".log", ".smcl", ".sthlp", ".hlp", ".dct", ".raw", ".prn"}


@command("unicode", "unicode")
def cmd_unicode(s: "Session", args: str) -> None:
    words, opt = _split(args)
    if not words:
        raise StataError(198, "invalid syntax")
    sub = words[0]
    out = s.output
    if sub == "encoding":
        if len(words) >= 2 and words[1] == "set":
            if len(words) != 3:
                raise StataError(198, "invalid syntax")
            import codecs
            try:
                codecs.lookup(words[2])
            except LookupError:
                raise StataError(198, f"encoding {words[2]} not found")   # VERIFY
            s.settings["unicode_encoding"] = words[2]
            out.write(f"(default encoding now {words[2]})\n", "text")   # VERIFY
            return
        if len(words) >= 2 and words[1] in ("list", "alias"):
            from encodings.aliases import aliases
            names = sorted(set(aliases.values()))
            pat = words[2].lower() if len(words) > 2 else ""
            for n in names:
                if pat.replace("*", "") in n:
                    out.write(f"  {n}\n", "text")
            return
        raise StataError(198, "invalid syntax")
    if sub == "convertfile":
        o = match_options(opt, {"srcencoding": 3, "dstencoding": 3, "srccallback": 4, "dstcallback": 4,
                                "replace": 7})
        if len(words) != 3:
            raise StataError(198, "invalid syntax")
        src, dst = _path(words[1]), _path(words[2])
        if not src.exists():
            raise StataError(601, f"file {words[1]} not found")
        if dst.exists() and not o.get("replace"):
            raise StataError(602, f"file {words[2]} already exists")
        senc = str(o.get("srcencoding")) if isinstance(o.get("srcencoding"), str) else "utf-8"
        denc = str(o.get("dstencoding")) if isinstance(o.get("dstencoding"), str) else "utf-8"
        errors = "replace" if o.get("srccallback") or o.get("dstcallback") else "strict"
        try:
            dst.write_bytes(src.read_bytes().decode(senc, errors).encode(denc, errors))
        except (UnicodeError, LookupError) as e:
            raise StataError(198, f"conversion failed: {e}")   # VERIFY
        return
    if sub in ("analyze", "translate", "retranslate"):
        o = match_options(opt, {"redo": 4, "nodata": 6, "invalid": 3, "transutf8": 5})
        if len(words) != 2:
            raise StataError(198, "invalid syntax")
        enc = s.settings.get("unicode_encoding")
        if sub != "analyze" and not enc:
            raise StataError(198, "you must first specify the encoding with unicode encoding set")
        p = _path(words[1])
        files = sorted(p.parent.glob(p.name)) if any(c in p.name for c in "*?[") else [p]
        files = [f for f in files if f.is_file() and (f.suffix in _TEXT_EXT or f.suffix == ".dta")]
        need = 0
        for f in files:
            if f.suffix == ".dta":
                bad = _dta_needs(f)
            else:
                bad = not _is_utf8(f.read_bytes())
            need += bad
            if sub == "analyze":
                state = "needs translation" if bad else "does not need translation"   # VERIFY wording
                out.write(f"  File {f.name} ({state})\n", "text")
            elif bad or sub == "retranslate":
                _translate(f, enc)
                out.write(f"  File {f.name} successfully translated\n", "text")   # VERIFY
        out.write(f"\n  Summary\n    {len(files)} file(s) specified\n    {need} file(s) needed translation\n",
                  "text")
        return
    raise StataError(198, f"unicode {sub} is not available in OpenDTA")


def _dta_needs(f: Path) -> bool:
    """Old (pre-Unicode) .dta files with bytes above 127 in their strings need translation."""
    with open(f, "rb") as h:
        head = h.read(32)
    if head.startswith(b"<stata_dta><header><release>118") or head.startswith(b"<stata_dta><header><release>119"):
        return False
    from ..io.dta import read_dta
    ds = read_dta(f)
    texts = [ds.label] + [v.label for v in ds.vars] + [t for v in ds.vars if v.is_string for t in v.raw]
    texts += [t for m in ds.value_labels.values() for t in m.values()]
    return any(any(ord(c) > 127 for c in t) for t in texts)


def _translate(f: Path, enc: str) -> None:
    """Translates to UTF-8, keeping the original in bak.stunicode/ (as Stata does)."""
    bak = f.parent / "bak.stunicode"
    bak.mkdir(exist_ok=True)
    shutil.copy2(f, bak / f.name)
    if f.suffix != ".dta":
        f.write_bytes(f.read_bytes().decode(enc, "replace").encode("utf-8"))
        return
    from ..io.dta import read_dta, write_dta

    def fix(t: str) -> str:
        # old files are read as latin-1: recover the bytes and decode them with the encoding
        return t.encode("latin-1", "replace").decode(enc, "replace")
    ds = read_dta(f)
    ds.label = fix(ds.label)
    for v in ds.vars:
        v.label = fix(v.label)
        if v.is_string:
            v.data = [fix(t) for t in v.raw]
    ds.value_labels = {k: {n: fix(t) for n, t in m.items()} for k, m in ds.value_labels.items()}
    write_dta(ds, f)
