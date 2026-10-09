"""infile (formato livre e com dicionário) e infix ([D] infile, [D] infix).

Mensagens e regras não documentadas estão marcadas VERIFICAR, com casos em
compat/do/0107_infile_infix.do.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..io.fixed import (FreeVar, parse_dictionary, parse_infix_spec, read_dictionary_free,
                        read_fixed, read_free, read_text)
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import touse
from .registry import command

if TYPE_CHECKING:
    from ..core.dataset import Dataset
    from ..session import Session

_TYPES = ("byte", "int", "long", "float", "double")


def _path(text: str, ext: str) -> Path:
    t = strip_outer_quotes(text.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    p = Path(t).expanduser()
    if p.suffix == "":
        p = p.with_suffix(ext)
    if not p.exists():
        raise StataError(601, f"file {p} not found")
    return p


def _echo_dictionary(s: "Session", text: str) -> None:
    """O infile com dicionário mostra o dicionário, da linha `dictionary`
    até a chave que fecha (observado no Stata 14, compat 0108)."""
    lines = text.splitlines()
    start = next((k for k, ln in enumerate(lines) if re.search(r"\bdictionary\b", ln)), None)
    if start is None:
        return
    for ln in lines[start:]:
        s.output.write(ln.rstrip() + "\n", "text")
        if ln.strip().startswith("}") or ln.rstrip().endswith("}"):
            break


def _check_empty(s: "Session", clear: bool) -> None:
    ds = s.data
    if ds.nvars and not clear:
        # observado no Stata 14 (infile; VERIFICAR o infix)
        raise StataError(18, "you must start with an empty dataset")


def _finish(s: "Session", new: "Dataset", if_: str | None, in_: str | None,
            warnings: list | None = None, partial: bool = False) -> None:
    names = [v.name for v in new.vars]
    if len(set(names)) != len(names):
        dup = next(n for n in names if names.count(n) > 1)
        raise StataError(110, f"variable {dup} already defined")
    old = s.data
    s.data = new
    try:
        mask = touse(s, Parsed(if_=if_, in_=in_)) if (if_ or in_) else None
    except Exception:
        s.data = old
        raise
    # avisos de valores não numéricos: o Stata numera pela posição que o
    # registro ocupa nos dados no momento da leitura (os registros descartados
    # antes pelo if não contam; observado no Stata 14)
    for k, (text, name, rec) in enumerate(warnings or []):
        if k >= 20:
            break
        obs = int(mask[:rec].sum()) + 1 if mask is not None else rec + 1
        s.output.write(f"'{text}' cannot be read as a number for {name}[{obs}]\n", "text")
    if partial:
        s.output.write("(eof not at end of obs)\n", "text")
    if mask is not None:
        new.keep_obs(mask)
    n = new.nobs
    s.output.write(f"({n:,} observation{'s' if n != 1 else ''} read)\n", "text")
    s.notify_state()


def _expand_names(token: str) -> list[str]:
    """x1-x5 -> x1 ... x5 (novas variáveis)."""
    m = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_]*?)(\d+)-\1(\d+)", token)
    if m:
        a, b = int(m.group(2)), int(m.group(3))
        step = 1 if b >= a else -1
        return [f"{m.group(1)}{k}" for k in range(a, b + step, step)]
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,31}", token):
        raise StataError(198, f"{token} invalid name")
    return [token]


def parse_free_spec(text: str, default_type: str) -> list[FreeVar]:
    """O tipo vale só para o nome seguinte (str10 nome idade: idade é
    numérica); type(a b c) aplica o tipo a todos os nomes entre parênteses."""
    spec: list[FreeVar] = []

    def group(m: re.Match) -> str:
        return " ".join(f"{m.group(1)} {n}" for n in m.group(2).split())
    text = re.sub(r"\b(byte|int|long|float|double|str\d*|strL)\s*\(([^)]*)\)", group, text)
    vtype = default_type
    for tok in text.split():
        m = re.fullmatch(r"_skip(?:\((\d+)\))?", tok)
        if m:
            spec.append(FreeVar("", "", int(m.group(1) or 1)))
            continue
        if tok in _TYPES or re.fullmatch(r"str\d*|strL", tok):
            vtype = tok
            continue
        lbl = ""
        if ":" in tok:
            tok, lbl = tok.split(":", 1)
        for name in _expand_names(tok):
            spec.append(FreeVar(name, vtype, value_label=lbl))
        vtype = default_type
    if not any(v.name for v in spec):
        raise StataError(100, "varlist required")
    return spec


@command("infile", "inf")
def cmd_infile(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.using is None:
        raise StataError(100, "using required")
    o = match_options(p.options, {"clear": 5, "automatic": 4, "byvariable": 2, "using": 5})
    _check_empty(s, bool(o.get("clear")))
    if not p.varlist.strip():
        # dicionário
        dpath = _path(p.using, ".dct")
        dtext = read_text(dpath)
        d = parse_dictionary(dtext, "infile")
        _echo_dictionary(s, dtext)
        if o.get("using"):
            data = read_text(_path(str(o["using"]), ".raw"))
        elif d.datafile:
            data = read_text(_path(d.datafile if Path(d.datafile).is_absolute()
                                   else str(dpath.parent / d.datafile), ".raw"))
        elif d.inline is not None:
            data = d.inline
        else:
            raise StataError(198, "dictionary does not specify a data file")
        new = read_dictionary_free(data, d.layout) if d.free else read_fixed(data, d.layout)
        _finish(s, new, p.if_, p.in_)
        return
    if o.get("byvariable"):
        raise StataError(198, "option byvariable() not yet supported by OpenDTA")
    spec = parse_free_spec(p.varlist, s.default_type())
    text = read_text(_path(p.using, ".raw"))
    new, warnings, partial = read_free(text, spec, automatic=bool(o.get("automatic")))
    _finish(s, new, p.if_, p.in_, warnings, partial)


@command("infix")
def cmd_infix(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.using is None:
        raise StataError(100, "using required")
    o = match_options(p.options, {"clear": 5, "using": 5})
    _check_empty(s, bool(o.get("clear")))
    if not p.varlist.strip() or p.varlist.strip() == "dictionary":
        dpath = _path(p.using, ".dct")
        d = parse_dictionary(read_text(dpath), "infix")
        if o.get("using"):
            data = read_text(_path(str(o["using"]), ".raw"))
        elif d.datafile:
            data = read_text(_path(d.datafile if Path(d.datafile).is_absolute()
                                   else str(dpath.parent / d.datafile), ".raw"))
        elif d.inline is not None:
            data = d.inline
        else:
            raise StataError(198, "dictionary does not specify a data file")
        layout = d.layout
    else:
        layout = parse_infix_spec(p.varlist)
        data = read_text(_path(p.using, ".raw"))
    if not layout.fields:
        raise StataError(198, "invalid syntax")
    _finish(s, read_fixed(data, layout), p.if_, p.in_)
