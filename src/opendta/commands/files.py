"""use, save, saveold, notes, char.

Mensagens seguem o Stata 14; as marcadas VERIFICAR têm casos em
compat/do/0104_arquivos.do.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from ..core.varlist import expand, resolve_name, unique
from ..io.dta import read_dta, write_dta
from ..lang.syntax import match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import touse
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def dta_path(text: str) -> Path:
    t = strip_outer_quotes(text.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    p = Path(os.path.expanduser(t))
    if p.suffix == "":
        p = p.with_suffix(".dta")
    return p


def _shown(p: Path) -> str:
    """Nome do arquivo como o Stata mostra nas mensagens (o que foi digitado)."""
    return str(p)


# ---------------------------------------------------------------------------
# use
# ---------------------------------------------------------------------------

@command("use", "u")
def cmd_use(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"clear": 5, "nolabel": 5})
    if p.using is not None:
        path = dta_path(p.using)
        want_vars, if_, in_ = p.varlist, p.if_, p.in_
    else:
        path = dta_path(p.varlist)
        want_vars, if_, in_ = "", p.if_, p.in_
    ds = s.data
    if ds.changed and ds.nvars and not opts.get("clear"):
        raise StataError(4, "no; data in memory would be lost")
    if not path.exists():
        raise StataError(601, f"file {_shown(path)} not found")
    new = read_dta(path)
    if opts.get("nolabel"):
        new.value_labels = {}
        for v in new.vars:
            v.value_label = ""
    keep = unique(expand(new, want_vars)) if want_vars else None
    old = s.data
    s.data = new
    try:
        if if_ or in_:
            # VERIFICAR: aqui o if pode usar variáveis fora do varlist
            # (avaliado antes de descartar as colunas).
            from ..lang.syntax import Parsed
            mask = touse(s, Parsed(if_=if_, in_=in_))
            new.keep_obs(mask)
    except Exception:
        s.data = old
        raise
    if keep is not None:
        new.drop_vars([n for n in new.names if n not in set(keep)])
    new.changed = bool(want_vars or if_ or in_)
    if new.label:
        s.output.write(f"({new.label})\n", "text")
    s.notify_state()


# ---------------------------------------------------------------------------
# save / saveold
# ---------------------------------------------------------------------------

def _save(s: "Session", args: str, release: int, extra_opts: dict | None = None,
          old_note: bool = False) -> None:
    head, comma, opts_text = args.partition(",")
    spec = {"replace": 3, "emptyok": 5}
    spec.update(extra_opts or {})
    opts = match_options(opts_text, spec) if comma else {}
    ds = s.data
    if head.strip():
        path = dta_path(head)
        shown = _shown(path)
    elif ds.fullpath or ds.filename:
        path = Path(ds.fullpath or ds.filename)
        shown = ds.filename or _shown(path)   # VERIFICAR: nome exibido em `save, replace`
    else:
        raise StataError(198, "invalid file specification")
    if ds.nvars == 0 and not opts.get("emptyok"):
        raise StataError(2000, "no variables defined")
    exists = path.exists()
    if exists and not opts.get("replace"):
        raise StataError(602, f"file {shown} already exists")
    if opts.get("version"):
        v = int(str(opts["version"]))
        if v == 13:
            release = 117
        elif v in (14, 15):
            release = 118
        else:
            raise StataError(198, f"option version({v}) not supported by OpenDTA (use 13 or 14)")
    if not exists and opts.get("replace"):
        s.output.write(f"(note: file {shown} not found)\n", "text")
    if old_note and release == 117:
        # VERIFICAR: avisos do saveold do Stata 14
        s.output.write("(saving in Stata 13 format)\n", "text")
        if not opts.get("version"):
            s.output.write("(FYI, saveold has options version(12) and version(11) "
                           "that write files in older Stata formats)\n", "text")
    try:
        ds.timestamp = write_dta(ds, path, release=release)
    except OSError as e:
        raise StataError(603, f"file {shown} could not be opened ({e.strerror})")
    if head.strip():
        ds.filename = str(path)
    ds.fullpath = str(path.resolve())
    ds.changed = False
    s.output.write(f"file {shown} saved\n", "text")
    s.notify_state()


@command("save", "sa")
def cmd_save(s: "Session", args: str) -> None:
    _save(s, args, 118)


@command("saveold")
def cmd_saveold(s: "Session", args: str) -> None:
    # saveold do Stata 14 grava para o Stata 13 (formato 117) por padrão
    _save(s, args, 117, {"version": 1}, old_note=True)


# ---------------------------------------------------------------------------
# notes (guardadas como características: note0 = quantidade, note1.. = textos)
# ---------------------------------------------------------------------------

def _notes_of(s: "Session", owner: str) -> list[str]:
    chars = s.data.chars.get(owner, {})
    try:
        n = int(chars.get("note0", "0"))
    except ValueError:
        n = 0
    return [chars.get(f"note{k}", "") for k in range(1, n + 1)]


def _set_notes(s: "Session", owner: str, notes: list[str]) -> None:
    chars = s.data.chars.setdefault(owner, {})
    for k in [k for k in chars if re.match(r"^note\d+$", k)]:
        del chars[k]
    if notes:
        chars["note0"] = str(len(notes))
        for k, text in enumerate(notes, start=1):
            chars[f"note{k}"] = text
    if not chars:
        s.data.chars.pop(owner, None)
    s.data.changed = True


@command("notes", "note")
def cmd_notes(s: "Session", args: str) -> None:
    t = args.strip()
    ds = s.data
    m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)?\s*:\s*(.*)$", t, re.S)
    if m:
        owner = m.group(1) or "_dta"
        if owner != "_dta":
            owner = resolve_name(ds, owner)
        text = m.group(2).strip()
        _set_notes(s, owner, _notes_of(s, owner) + [text])
        return
    words = t.split()
    if words and words[0] == "drop":
        targets = words[1:] or ["_all"]
        owners = ["_dta"] + ds.names if targets == ["_all"] else [
            "_dta" if w == "_dta" else resolve_name(ds, w) for w in targets]
        for o in owners:
            if _notes_of(s, o):
                _set_notes(s, o, [])
        return
    if words and words[0] in ("list", "l", "li", "lis"):
        words = words[1:]
    owners = ["_dta"] + ds.names if not words or words == ["_all"] else [
        "_dta" if w == "_dta" else resolve_name(ds, w) for w in words]
    out = s.output
    # VERIFICAR: layout da listagem de notes
    for o in owners:
        notes = _notes_of(s, o)
        if not notes:
            continue
        out.write(f"\n{o}:\n", "text")
        for k, text in enumerate(notes, start=1):
            out.write(f"{k:>3}.  {text}\n", "text")


# ---------------------------------------------------------------------------
# char
# ---------------------------------------------------------------------------

@command("char")
def cmd_char(s: "Session", args: str) -> None:
    t = args.strip()
    ds = s.data
    sub, _, rest = t.partition(" ")
    if sub in ("define", "def", "de", "defi", "defin") or "[" in sub:
        spec = rest if sub.startswith("de") and "[" not in sub else t
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)\[([A-Za-z_][A-Za-z0-9_]*)\]\s*(.*)$", spec.strip(), re.S)
        if not m:
            raise StataError(198, "invalid syntax")
        owner = m.group(1) if m.group(1) == "_dta" else resolve_name(ds, m.group(1))
        value = strip_outer_quotes(m.group(3)) if m.group(3).strip() else ""
        chars = ds.chars.setdefault(owner, {})
        if value:
            chars[m.group(2)] = value
        else:
            chars.pop(m.group(2), None)
            if not chars:
                ds.chars.pop(owner, None)
        ds.changed = True
        return
    if sub in ("list", "l", "li", "lis", ""):
        targets = rest.split()
        out = s.output
        for owner, chars in ds.chars.items():
            if targets and owner not in targets:
                continue
            for name, value in chars.items():
                out.write(f"  {owner}[{name}]:".ljust(30) + f"{value}\n", "text")
        return
    raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# browse / edit (abrem o Data Browser na interface gráfica)
# ---------------------------------------------------------------------------

def _browse(s: "Session", args: str) -> None:
    from ..lang.syntax import Parsed
    p = parse_standard(args)
    if p.options.strip():
        match_options(p.options, {"nolabel": 3})
    names = unique(expand(s.data, p.varlist)) if p.varlist.strip() else None
    rows = None
    if p.if_ or p.in_:
        rows = np.flatnonzero(touse(s, Parsed(if_=p.if_, in_=p.in_)))
    hook = s.ui_hooks.get("browse")
    if hook is not None:
        # sem interface (modo batch) o Stata também não mostra nada
        hook(names, rows)


@command("browse", "br")
def cmd_browse(s: "Session", args: str) -> None:
    _browse(s, args)


@command("edit", "ed")
def cmd_edit(s: "Session", args: str) -> None:
    # a edição de células chega na fase de ferramentas; por ora abre em Browse
    _browse(s, args)


# ---------------------------------------------------------------------------
# type (mostra um arquivo de texto)
# ---------------------------------------------------------------------------

@command("type", "ty")
def cmd_type(s: "Session", args: str) -> None:
    head, comma, opts_text = args.partition(",")
    o = match_options(opts_text, {"asis": 4, "smcl": 4, "showtabs": 8, "starbang": 8,
                                  "lines": 4}) if comma else {}
    t = strip_outer_quotes(head.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    path = Path(os.path.expanduser(t))
    if not path.exists():
        raise StataError(601, f"file {t} not found")
    raw = path.read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if o.get("lines"):
        lines = lines[:int(str(o["lines"]))]
    out = []
    for ln in lines:
        if not o.get("starbang") and ln.startswith("*!"):
            continue   # VERIFICAR: linhas *! só aparecem com starbang
        # VERIFICAR: tabulações expandidas a cada 8 colunas; showtabs mostra <T>
        ln = ln.replace("\t", "<T>") if o.get("showtabs") else ln.expandtabs(8)
        out.append(ln)
    s.output.write("".join(x + "\n" for x in out), "text")


# ---------------------------------------------------------------------------
# erase / rm
# ---------------------------------------------------------------------------

def _erase(s: "Session", args: str) -> None:
    t = strip_outer_quotes(args.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    path = Path(os.path.expanduser(t))
    if not path.exists():
        raise StataError(601, f"file {t} not found")
    try:
        path.unlink()
    except OSError as e:
        raise StataError(693, f"could not erase file {t} ({e.strerror})")


@command("erase")
def cmd_erase(s: "Session", args: str) -> None:
    _erase(s, args)


@command("rm")
def cmd_rm(s: "Session", args: str) -> None:
    _erase(s, args)
