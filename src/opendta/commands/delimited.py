"""import delimited, export delimited, insheet, outsheet ([D] import delimited).

Mensagens e regras não documentadas estão marcadas VERIFICAR, com casos em
compat/do/0105_delimitado.do.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from ..core.varlist import expand, unique
from ..io.delimited import ReadOptions, WriteOptions, read_delimited, write_delimited
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lang.words import parse_numlist, strip_outer_quotes
from ._util import touse
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def text_path(text: str, default_ext: str) -> Path:
    t = strip_outer_quotes(text.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    p = Path(t).expanduser()
    if p.suffix == "":
        p = p.with_suffix(default_ext)
    return p


def _delimiter(arg: str) -> str:
    a = strip_outer_quotes(arg.split(",")[0].strip())
    if a in ("tab", "\\t"):
        return "\t"
    if a in ("comma",):
        return ","
    if len(a) != 1:
        # VERIFICAR: delimitadores com mais de um caractere
        raise StataError(198, "delimiters() must be a single character")
    return a


def _cols(arg: str) -> set[int]:
    a = arg.strip()
    if a == "_all":
        return {-1}
    return {int(x) for x in parse_numlist(a)}


def _range(arg: str) -> tuple[int | None, int | None]:
    a, _, b = arg.strip().partition(":")
    return (int(a) if a.strip() else None, int(b) if b.strip() else None)


def _load(s: "Session", path: Path, opts: ReadOptions, clear: bool, names: list[str]) -> None:
    ds = s.data
    if ds.changed and ds.nvars and not clear:
        raise StataError(4, "no; data in memory would be lost")
    if not path.exists():
        raise StataError(601, f"file {path} not found")
    new = read_delimited(path, opts)
    for k, n in enumerate(names):
        if k < new.nvars:
            new.vars[k].name = n
    s.data = new
    s.output.write(f"({new.nvars} var{'s' if new.nvars != 1 else ''}, "
                   f"{new.nobs:,} ob{'s' if new.nobs != 1 else ''})\n", "text")
    s.notify_state()


def _resolve_all(opts: ReadOptions, path: Path) -> None:
    """stringcols(_all)/numericcols(_all) dependem do número de colunas."""
    if -1 in opts.stringcols or -1 in opts.numericcols:
        probe = read_delimited(path, ReadOptions(delimiter=opts.delimiter, varnames=0,
                                                 encoding=opts.encoding, rowrange=(1, 1)))
        allc = set(range(1, probe.nvars + 1))
        if -1 in opts.stringcols:
            opts.stringcols = allc
        if -1 in opts.numericcols:
            opts.numericcols = allc


def _split_using(args: str) -> tuple[str, str, str]:
    """[namelist] using f, opts | f, opts -> (namelist, arquivo, opções)."""
    p = parse_standard(args)
    if p.using is not None:
        return p.varlist, p.using, p.options
    return "", p.varlist, p.options


def import_delimited(s: "Session", args: str) -> None:
    namelist, file, opt_text = _split_using(args)
    o = match_options(opt_text, {
        "delimiters": 5, "varnames": 7, "case": 4, "asdouble": 8, "asfloat": 7,
        "clear": 5, "rowrange": 8, "colrange": 8, "encoding": 8, "stringcols": 10,
        "numericcols": 11, "bindquote": 9, "stripquotes": 11})
    opts = ReadOptions()
    if o.get("delimiters"):
        opts.delimiter = _delimiter(str(o["delimiters"]))
    if "varnames" in o:
        v = str(o["varnames"]).strip()
        opts.varnames = 0 if v == "nonames" else int(v)
    if o.get("case"):
        c = str(o["case"]).strip()
        if c not in ("preserve", "lower", "upper"):
            raise StataError(198, "option case() invalid")
        opts.case = c
    opts.asdouble = bool(o.get("asdouble"))
    if o.get("rowrange"):
        opts.rowrange = _range(str(o["rowrange"]))
    if o.get("colrange"):
        opts.colrange = _range(str(o["colrange"]))
    if o.get("encoding"):
        opts.encoding = strip_outer_quotes(str(o["encoding"]).strip())
    if o.get("stringcols"):
        opts.stringcols = _cols(str(o["stringcols"]))
    if o.get("numericcols"):
        opts.numericcols = _cols(str(o["numericcols"]))
    if o.get("bindquote"):
        opts.bindquote = str(o["bindquote"]).strip()
    path = text_path(file, ".csv")
    if path.exists():
        _resolve_all(opts, path)
    _load(s, path, opts, bool(o.get("clear")), namelist.split())


def _export(s: "Session", args: str, opts: WriteOptions, spec: dict[str, int], *,
            old: bool) -> None:
    p = parse_standard(args)
    if p.using is not None:
        names_text, file = p.varlist, p.using
    elif old:
        raise StataError(100, "using required")
    else:
        names_text, file = "", p.varlist
    o = match_options(p.options, spec)
    ds = s.data
    if ds.nvars == 0:
        raise StataError(102, "no variables defined")   # VERIFICAR
    names = unique(expand(ds, names_text)) if names_text.strip() else list(ds.names)
    path = text_path(file, ".out" if old and not o.get("comma") else ".csv")
    if old and o.get("comma"):
        opts.delimiter = ","
    if o.get("delimiter"):
        opts.delimiter = _delimiter(str(o["delimiter"]))
    if o.get("novarnames") or o.get("nonames"):
        opts.varnames = False
    opts.nolabel = bool(o.get("nolabel"))
    opts.datafmt = bool(o.get("datafmt"))
    if o.get("quote"):
        opts.quote = True
    if o.get("noquote"):
        opts.quote = False
    if path.exists() and not o.get("replace"):
        raise StataError(602, f"file {path} already exists")
    rows = np.flatnonzero(touse(s, Parsed(if_=p.if_, in_=p.in_)))
    try:
        write_delimited(ds, path, names, rows, opts)
    except OSError as e:
        raise StataError(603, f"file {path} could not be opened ({e.strerror})")
    # VERIFICAR: o Stata 14 não mostra mensagem ao exportar


def export_delimited(s: "Session", args: str) -> None:
    _export(s, args, WriteOptions(), {
        "delimiter": 5, "novarnames": 6, "nolabel": 5, "datafmt": 7, "quote": 5,
        "replace": 7}, old=False)


@command("import")
def cmd_import(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    if len(sub) >= 5 and "delimited".startswith(sub):
        import_delimited(s, rest)
        return
    if len(sub) >= 3 and "excel".startswith(sub):
        from .excel import import_excel
        import_excel(s, rest)
        return
    raise StataError(199, f"import {sub} not yet available in OpenDTA")


@command("export")
def cmd_export(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    if len(sub) >= 5 and "delimited".startswith(sub):
        export_delimited(s, rest)
        return
    if len(sub) >= 3 and "excel".startswith(sub):
        from .excel import export_excel
        export_excel(s, rest)
        return
    raise StataError(199, f"export {sub} not yet available in OpenDTA")


@command("insheet")
def cmd_insheet(s: "Session", args: str) -> None:
    namelist, file, opt_text = _split_using(args)
    o = match_options(opt_text, {"double": 6, "tab": 3, "comma": 5, "delimiter": 9,
                                 "clear": 5, "case": 4, "names": 5, "nonames": 7})
    opts = ReadOptions()
    if o.get("tab"):
        opts.delimiter = "\t"
    if o.get("comma"):
        opts.delimiter = ","
    if o.get("delimiter"):
        opts.delimiter = _delimiter(str(o["delimiter"]))
    opts.case = "preserve" if o.get("case") else "lower"
    if o.get("names"):
        opts.varnames = 1
    if o.get("nonames"):
        opts.varnames = 0
    opts.asdouble = bool(o.get("double"))
    _load(s, text_path(file, ".raw"), opts, bool(o.get("clear")), namelist.split())


@command("outsheet")
def cmd_outsheet(s: "Session", args: str) -> None:
    # outsheet: tabulação e aspas em todas as strings por padrão
    _export(s, args, WriteOptions(delimiter="\t", quote=True), {
        "comma": 5, "delimiter": 9, "nonames": 7, "nolabel": 5, "noquote": 7,
        "replace": 7}, old=True)
