"""SAS, SPSS and SAS XPORT files through ReadStat (pyreadstat, optional like for old .dta).

`dataset_from_columns` is the conversion to a Dataset; it takes the columns and the metadata
that pyreadstat returns, so that it can be tested without the package. Rules that the manual does
not give are marked VERIFY.
"""

from __future__ import annotations

import datetime as _dt
from pathlib import Path

import numpy as np

from ..core import missing as M
from ..core.dataset import NAME_RE, RESERVED, Dataset, Variable, str_type_for
from ..core.errors import StataError

_EPOCH = _dt.datetime(1960, 1, 1)
# largest and smallest values of the integer types (the rest are missing values)
_INT_RANGES = (("byte", -127, 100), ("int", -32767, 32740), ("long", -2147483647, 2147483620))


def pyreadstat():
    try:
        import pyreadstat as mod
    except ImportError:
        raise StataError(610, "this file format requires the pyreadstat package (pip install pyreadstat)")
    return mod


def legal_names(raw: list[str], case: str) -> tuple[list[str], dict[str, str]]:
    """Names that Stata accepts; the others become v# (k = position) and the original name is
    returned to be kept as a characteristic."""
    names: list[str] = []
    original: dict[str, str] = {}
    for k, name in enumerate(raw, start=1):
        n = name.lower() if case == "lower" else name.upper() if case == "upper" else name
        if not NAME_RE.match(n) or n in RESERVED or n in names:
            n = f"v{k}"
            original[n] = name
        names.append(n)
    return names, original


def _numeric(values) -> tuple[str, np.ndarray]:
    arr = np.array([M.SYSMISS if (x is None or x != x) else float(x) for x in values], dtype=np.float64)
    ok = arr[arr < M.SYSMISS]
    if ok.size == 0 or not np.all(ok == np.floor(ok)):
        return ("byte" if ok.size == 0 else "double"), arr
    for vtype, lo, hi in _INT_RANGES:       # VERIFY: the smallest integer type that holds the values
        if ok.min() >= lo and ok.max() <= hi:
            return vtype, arr
    return "double", arr


def _date(values) -> tuple[str, str, np.ndarray] | None:
    """%td for dates and %tc for date and time; None when the column holds neither."""
    kinds = {type(x) for x in values if x is not None and x == x}
    if kinds and kinds <= {_dt.date}:
        days = [M.SYSMISS if (x is None or x != x) else float((x - _EPOCH.date()).days) for x in values]
        return "double", "%td", np.array(days)
    if kinds and kinds <= {_dt.datetime} | {_dt.date}:
        def ms(x):
            if x is None or x != x:
                return M.SYSMISS
            d = x if isinstance(x, _dt.datetime) else _dt.datetime(x.year, x.month, x.day)
            return (d - _EPOCH).total_seconds() * 1000.0
        return "double", "%tc", np.array([ms(x) for x in values])
    return None


def dataset_from_columns(cols: dict, meta, *, case: str = "preserve", label_max: int = 80,
                         xport: bool = False) -> Dataset:
    """cols: name -> list of values; meta: the object of pyreadstat (column_names,
    column_labels, variable_value_labels, file_label). xport: Stata 14 keeps what the XPORT file
    declares (every number is a double, a string has the width of the file; compat 1701, VERIFY)."""
    raw_names = list(meta.column_names)
    names, original = legal_names(raw_names, case)
    labels = list(getattr(meta, "column_labels", None) or [None] * len(raw_names))
    value_labels = getattr(meta, "variable_value_labels", None) or {}
    ds = Dataset()
    nobs = len(next(iter(cols.values()), []))
    ds.nobs = nobs
    nlab = 0
    for raw, name, label in zip(raw_names, names, labels):
        col = list(cols[raw])
        label = (label or "")[:label_max]
        if any(isinstance(x, str) for x in col):
            col = ["" if (x is None or x != x) else str(x) for x in col]
            vtype = str_type_for(max([len(x.encode("utf-8")) for x in col] + [1]))
            width = (getattr(meta, "variable_storage_width", None) or {}).get(raw) if xport else None
            if width and 0 < int(width) <= 2045:
                vtype = f"str{int(width)}"
            ds.vars.append(Variable(name, vtype, col, label=label))
        else:
            dated = _date(col)
            if dated:
                vtype, fmt, arr = dated
                ds.vars.append(Variable(name, vtype, arr, fmt=fmt, label=label))
            else:
                vtype, arr = _numeric(col)
                if xport:
                    vtype = "double"
                var = Variable(name, vtype, arr, label=label)
                vl = value_labels.get(raw)
                if vl and all(float(k) == int(k) for k in vl):
                    nlab += 1                  # VERIFY: labels1, labels2 ... in the order of the variables
                    ds.value_labels[f"labels{nlab}"] = {int(k): str(v) for k, v in vl.items()}
                    var.value_label = f"labels{nlab}"
                ds.vars.append(var)
        if name in original:
            ds.chars.setdefault(name, {})["original_name"] = original[name]    # VERIFY: name of the characteristic
    ds.label = (getattr(meta, "file_label", "") or "")[:label_max]
    ds.changed = True
    return ds


def select(ds: Dataset, namelist: list[str]) -> Dataset:
    """Keeps the variables of the namelist, in the order of the file."""
    wanted = {n.lower() for n in namelist}
    missing = wanted - {v.name.lower() for v in ds.vars}
    if missing:
        raise StataError(111, f"variable {sorted(missing)[0]} not found")
    ds.vars = [v for v in ds.vars if v.name.lower() in wanted]
    return ds


def read_sas(path: Path, *, bcat: Path | None = None, encoding: str | None = None) -> tuple[dict, object]:
    rs = pyreadstat()
    kw = {"encoding": encoding} if encoding else {}
    if bcat is not None:
        kw["catalog_file"] = str(bcat)
    try:
        return rs.read_sas7bdat(str(path), output_format="dict", **kw)
    except Exception as e:  # noqa: BLE001
        raise StataError(610, f"file not SAS format ({e})")


def read_spss(path: Path, *, zsav: bool = False, encoding: str | None = None) -> tuple[dict, object]:
    rs = pyreadstat()
    kw = {"encoding": encoding} if encoding else {}
    try:
        return (rs.read_zsav if zsav and hasattr(rs, "read_zsav") else rs.read_sav)(
            str(path), output_format="dict", **kw)
    except Exception as e:  # noqa: BLE001
        raise StataError(610, f"file not SPSS format ({e})")


def read_xport(path: Path) -> tuple[dict, object]:
    rs = pyreadstat()
    try:
        return rs.read_xport(str(path), output_format="dict")
    except Exception as e:  # noqa: BLE001
        raise StataError(610, f"file not SAS XPORT format ({e})")


def write_xport(ds: Dataset, path: Path, names: list[str], rows: np.ndarray) -> None:
    """XPORT version 5 (names of up to 8 characters, no value labels)."""
    rs = pyreadstat()
    try:
        import pandas as pd
    except ImportError:
        raise StataError(610, "export sasxport requires the pandas package (pip install pandas)")
    frame = {}
    for n in names:
        v = ds.get(n)
        if v.is_string:
            # the file keeps the width of the variable (str6 comes back as str6)
            w = int(v.vtype[3:]) if v.vtype[3:].isdigit() else 0
            frame[n] = [str(v.value(int(i))).ljust(w) for i in rows]
        else:
            x = v.data[rows].astype(float)
            frame[n] = np.where(x >= M.SYSMISS, np.nan, x)
    labels = [ds.get(n).label or None for n in names]
    try:
        rs.write_xport(pd.DataFrame(frame), str(path), file_label=ds.label or "", column_labels=labels,
                       file_format_version=5)
    except Exception as e:  # noqa: BLE001
        raise StataError(603, f"file {path} could not be written ({e})")
