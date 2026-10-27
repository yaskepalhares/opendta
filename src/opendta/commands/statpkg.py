"""import sas, import spss, import sasxport and export sasxport ([D] import sas, [D] import spss,
[D] import sasxport); the readers are in io/statpkg.py.

Messages and undocumented rules are marked VERIFY, with the sasxport cases in
compat/do/1701_datamgmt_options.do (Stata 14 has import sasxport but not import sas or spss).
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..core.dataset import Dataset
from ..core.errors import StataError
from ..core.varlist import expand, unique
from ..io import statpkg
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import touse

if TYPE_CHECKING:
    from ..session import Session


def _file(text: str, ext: str, *, must_exist: bool = True) -> Path:
    t = strip_outer_quotes(text.strip())
    if not t:
        raise StataError(198, "invalid file specification")
    p = Path(t).expanduser()
    if p.suffix == "":
        p = p.with_suffix(ext)
    if must_exist and not p.exists():
        raise StataError(601, f"file {p} not found")
    return p


def _install(s: "Session", new: Dataset, p: Parsed, namelist: list[str], quiet: bool = False) -> None:
    if namelist:
        statpkg.select(new, namelist)
    old = s.data
    s.data = new
    try:
        if p.if_ or p.in_:
            new.keep_obs(touse(s, Parsed(if_=p.if_, in_=p.in_)))
    except Exception:
        s.data = old
        raise
    s.r = {"N": float(new.nobs), "k": float(new.nvars)}
    if not quiet:
        s.output.write(f"({new.nvars} vars, {new.nobs} obs)\n", "text")    # import sasxport is silent (compat 1701)
    s.notify_state()


def _check_replace(s: "Session", clear: bool) -> None:
    if s.data.changed and s.data.nvars and not clear:
        raise StataError(4, "no; data in memory would be lost")


def _namelist(text: str) -> list[str]:
    return text.split()


def import_sas(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"bcat": 4, "case": 4, "clear": 5, "encoding": 8})
    names = _namelist(p.varlist) if p.using is not None else []
    path = _file(p.using if p.using is not None else p.varlist, ".sas7bdat")
    case = str(o.get("case", "preserve")).strip()
    if case not in ("lower", "upper", "preserve"):
        raise StataError(198, "option case() invalid")
    _check_replace(s, bool(o.get("clear")))
    bcat = _file(str(o["bcat"]), ".sas7bcat") if o.get("bcat") else None
    cols, meta = statpkg.read_sas(path, bcat=bcat, encoding=_enc(o))
    _install(s, statpkg.dataset_from_columns(cols, meta, case=case), p, names)


def import_spss(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"case": 4, "clear": 5, "zsav": 4, "encoding": 8})
    names = _namelist(p.varlist) if p.using is not None else []
    zsav = bool(o.get("zsav"))
    path = _file(p.using if p.using is not None else p.varlist, ".zsav" if zsav else ".sav")
    case = str(o.get("case", "preserve")).strip()
    if case not in ("lower", "upper", "preserve"):
        raise StataError(198, "option case() invalid")
    _check_replace(s, bool(o.get("clear")))
    cols, meta = statpkg.read_spss(path, zsav=zsav, encoding=_enc(o))
    # the label of a variable has up to 256 characters in the SPSS format
    _install(s, statpkg.dataset_from_columns(cols, meta, case=case, label_max=256), p, names)


def _enc(o: dict) -> str | None:
    return strip_outer_quotes(str(o["encoding"]).strip()) if o.get("encoding") else None


def _formats_file(folder: Path) -> Path | None:
    for f in ("formats.xpf", "FORMATS.xpf", "FORMATS.XPF"):
        if (folder / f).exists():
            return folder / f
    return None


def _xport_describe(s: "Session", path: Path, member: str) -> None:
    """import sasxport, describe: each member (or the one in member()) with its variables, in the
    layout of compat 1911. VERIFY: files with several members, string variables."""
    from ..io import xport
    members = xport.read(path)
    total = len(members)
    if member:
        members = [_pick(members, member)]
    out = s.output
    rule = "-" * 76 + "\n"
    for m in members:
        n, k = len(m.rows), len(m.variables)
        # size: the observations as stored in the file plus 4 bytes each (compat 1911)
        size = n * (sum(v.length for v in m.variables) + 4)
        out.write(f"\ndata from {path.name}, member({m.name.lower()})\n", "text")
        out.write(f"  obs:{n:>17,}             {m.created.lower()}\n", "text")
        out.write(f" vars:{k:>17,}             (date shown exactly as recorded in file)\n", "text")
        out.write(f" size:{size:>17,}\n" + rule, "text")
        out.write("variable    variable     value\n  name        type       label      variable label\n"
                  + rule, "text")
        for v in m.variables:
            kind = "numeric" if v.numeric else f"char({v.length})"           # VERIFY: strings
            out.write(f"{v.name.lower():<12}{kind:<12}{v.format.lower():<12}{v.label}".rstrip() + "\n",
                      "text")
        out.write(rule, "text")
    # return list shows the results in the reverse order of this dictionary
    s.r = {"members": " ".join(m.name.lower() for m in members), "n_members": float(total),
           "size": float(size), "k": float(k), "N": float(n)}


def _pick(members: list, name: str):
    for m in members:
        if m.name.lower() == name.lower():
            return m
    raise StataError(198, f"member {name} not found")              # VERIFY: return code


def import_sasxport(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"clear": 5, "novallabels": 5, "member": 6, "describe": 8})
    path = _file(p.varlist, ".xpt")
    member = str(o["member"]).strip() if o.get("member") not in (None, True) else ""
    if o.get("describe"):
        _xport_describe(s, path, member)
        return
    _check_replace(s, bool(o.get("clear")))
    fmt_file = _formats_file(path.parent)
    if not member and fmt_file is None:
        cols, meta = statpkg.read_xport(path)
        _install(s, statpkg.dataset_from_columns(cols, meta, case="lower", xport=True), Parsed(), [], quiet=True)
        return
    from types import SimpleNamespace
    from ..io import xport
    try:
        members = xport.read(path)
    except (IndexError, ValueError, struct.error):
        raise StataError(610, f"file {path} not SAS XPORT format")
    m = _pick(members, member) if member else members[0]
    cols = {v.name: [(float("nan") if x is None else x) if v.numeric else x for x in
                     (row[k] for row in m.rows)] for k, v in enumerate(m.variables)}
    meta = SimpleNamespace(column_names=[v.name for v in m.variables],
                           column_labels=[v.label or None for v in m.variables], file_label=m.label,
                           variable_storage_width={v.name: v.length for v in m.variables if not v.numeric})
    ds = statpkg.dataset_from_columns(cols, meta, case="lower", xport=True)
    if fmt_file is not None:
        # the value labels of formats.xpf go to the variables whose SAS format has their name;
        # with novallabels the variables keep the names but the labels are not defined (compat
        # 1911; VERIFY)
        labels = xport.value_labels(fmt_file)
        for var, xv in zip(ds.vars, m.variables):
            fname = xv.format.lower()
            if fname in labels and not var.is_string:
                if not o.get("novallabels"):
                    ds.value_labels[fname] = labels[fname]
                var.value_label = fname
    _install(s, ds, Parsed(), [], quiet=True)


def _sas_names(names: list[str], rename: bool, what: str, out=None) -> dict[str, str]:
    """Names that SAS XPORT accepts: up to 8 characters, the same in upper and lower case. With
    rename, a longer or clashing name is cut, upper-cased and given a number (VERIFY: Stata's rule
    for clashes). out: where the list of renamed variables goes (compat 1911)."""
    if out is not None:
        new = _sas_names(names, True, what)
        changed = [(n, v.upper()) for n, v in new.items() if v != n]
        if changed:
            items = "".join(f"                {n} -> {v}\n" for n, v in changed)
            if not rename:
                out.write("the following variable(s) have names that must be changed to fit into "
                          ".xpt format:  (suggested renamings shown):\n\n" + items + "\n", "text")
                raise StataError(198, "specify option rename to save .xpt file with suggested names")
            out.write("the following variable(s) were renamed in the output file:\n\n" + items + "\n",
                      "text")
        return {n: v.upper() if v != n else v for n, v in new.items()}
    out: dict[str, str] = {}
    taken: set[str] = set()
    for n in names:
        new = n
        if len(n) > 8 or n.upper() in taken:
            if not rename:
                raise StataError(198, f"{what} {n} has more than 8 characters or clashes with another; "
                                      "use option rename")               # VERIFY: Stata's message
            new, k = n[:8], 1
            while new.upper() in taken:
                tail = str(k)
                new, k = n[:8 - len(tail)] + tail, k + 1
        taken.add(new.upper())
        out[n] = new
    return out


def export_sasxport(s: "Session", args: str) -> None:
    """export sasxport [varlist using] filename [if] [in] [, rename replace vallabfile()]"""
    p = parse_standard(args)
    o = match_options(p.options, {"rename": 6, "replace": 7, "vallabfile": 6})
    ds = s.data
    if ds.nvars == 0:
        raise StataError(102, "no variables defined")   # VERIFY
    if p.using is not None:
        names, file = unique(expand(ds, p.varlist)), p.using
    else:
        names, file = list(ds.names), p.varlist
    vallab = str(o.get("vallabfile", "xpf")).strip() if o.get("vallabfile") not in (None, True) else "xpf"
    if vallab not in ("xpf", "sascode", "both", "none"):
        raise StataError(198, "option vallabfile() incorrectly specified")      # VERIFY
    rename = bool(o.get("rename"))
    newnames = _sas_names(names, rename, "variable", s.output)
    # the value labels in the order of the variables that use them (compat 1911)
    used = list(dict.fromkeys(ds.get(n).value_label for n in names if ds.get(n).value_label
                              and ds.get(n).value_label in ds.value_labels)) if vallab != "none" else []
    labnames = _sas_names(used, rename, "value label")
    path = _file(file, ".xpt", must_exist=False)
    if path.exists() and not o.get("replace"):
        raise StataError(602, f"file {path} already exists")
    rows = np.flatnonzero(touse(s, Parsed(if_=p.if_, in_=p.in_)))
    formats = {newnames[n]: labnames[ds.get(n).value_label].upper() for n in names
               if ds.get(n).value_label in labnames}
    statpkg.write_xport(ds, path, names, rows, rename=newnames, formats=formats)
    s.output.write(f"file {path} saved\n", "text")
    sas = path.with_suffix(".sas")
    if vallab in ("xpf", "none") and sas.exists():
        sas.unlink()        # Stata leaves no .sas file of an earlier export (compat 1911; VERIFY)
    if used:
        sets = {labnames[n]: ds.value_labels[n] for n in used}
        if vallab in ("sascode", "both"):
            _write_sascode(sas, sets)
            s.output.write(f"file {sas} saved\n", "text")
        if vallab in ("xpf", "both"):
            statpkg.write_formats_xpf(path.parent / "formats.xpf", sets)
            s.output.write("file formats.xpf saved\n", "text")


def _write_sascode(path: Path, sets: dict[str, dict]) -> None:
    """The value labels as a SAS program with PROC FORMAT (VERIFY: Stata's file)."""
    lines = ["proc format;"]
    for name, lab in sets.items():
        lines.append(f"value {name}")
        for k in sorted(lab):
            text = str(lab[k]).replace('"', '""')
            lines.append(f'  {k} = "{text}"')
        lines.append(";")
    lines.append("run;")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
