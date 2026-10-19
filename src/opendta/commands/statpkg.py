"""import sas, import spss, import sasxport and export sasxport ([D] import sas, [D] import spss,
[D] import sasxport); the readers are in io/statpkg.py.

Messages and undocumented rules are marked VERIFY, with the sasxport cases in
compat/do/1701_datamgmt_options.do (Stata 14 has import sasxport but not import sas or spss).
"""

from __future__ import annotations

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


def import_sasxport(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"clear": 5, "novallabels": 5, "member": 6, "describe": 8})
    if o.get("describe") or o.get("member"):
        raise StataError(198, "options describe and member() are not available yet in OpenDTA")
    path = _file(p.varlist, ".xpt")
    _check_replace(s, bool(o.get("clear")))
    if not o.get("novallabels"):
        for f in ("formats.xpf", "FORMATS.xpf"):
            if (path.parent / f).exists():
                raise StataError(198, f"value labels in {f} are not available yet in OpenDTA; use novallabels")
    cols, meta = statpkg.read_xport(path)
    _install(s, statpkg.dataset_from_columns(cols, meta, case="lower", xport=True), Parsed(), [], quiet=True)


def export_sasxport(s: "Session", args: str) -> None:
    """export sasxport [varlist using] filename [if] [in] [, rename replace vallabfile(none)]"""
    p = parse_standard(args)
    o = match_options(p.options, {"rename": 6, "replace": 7, "vallabfile": 6})
    ds = s.data
    if ds.nvars == 0:
        raise StataError(102, "no variables defined")   # VERIFY
    if p.using is not None:
        names, file = unique(expand(ds, p.varlist)), p.using
    else:
        names, file = list(ds.names), p.varlist
    if str(o.get("vallabfile", "")).strip() != "none" and any(ds.get(n).value_label for n in names):
        raise StataError(198, "formats.xpf is not available yet in OpenDTA; use vallabfile(none)")
    if o.get("rename"):
        raise StataError(198, "option rename is not available yet in OpenDTA")
    bad = [n for n in names if len(n) > 8]
    if bad:
        # VERIFY: Stata's message; the manual says an error is issued and no file is saved
        raise StataError(198, f"variable {bad[0]} has more than 8 characters; use shorter names")
    path = _file(file, ".xpt", must_exist=False)
    if path.exists() and not o.get("replace"):
        raise StataError(602, f"file {path} already exists")
    rows = np.flatnonzero(touse(s, Parsed(if_=p.if_, in_=p.in_)))
    statpkg.write_xport(ds, path, names, rows)
    s.output.write(f"file {path} saved\n", "text")
