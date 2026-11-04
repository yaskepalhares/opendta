"""Multiple imputation: mi set, mi register, mi impute, mi estimate,
mi passive, mi xeq, mi extract, mi describe, mi misstable ([MI] manual).

Storage follows Stata's mlong and flong styles: _mi_m (0 for the original
data), _mi_id (observation identifier) and _mi_miss (1 for observations
incomplete in m=0), with the settings in the _dta characteristics.

mi impute follows [MI] mi impute regress/pmm/logit/ologit/mlogit/poisson,
Methods and formulas: the parameters are drawn from their posterior (scaled
chi-square for the variance, normal for the coefficients) and then the
missing values. The draws come from OpenDTA's random-number generator in
an order Stata does not document, so imputed values differ from Stata's
(like bootstrap; compat/known_differences.txt). mi estimate combines the
estimates with Rubin's rules, with the small-sample degrees of freedom
of Barnard and Rubin (1999) when the command has residual degrees of
freedom, and the model F test of Li, Raghunathan and Rubin (1991).
"""

from __future__ import annotations

import time
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top
from ..lazy import lazy

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..core.dataset import Dataset
    from ..session import Session

SYS = M.SYSMISS
SYSV = ("_mi_m", "_mi_id", "_mi_miss")


# ---------------------------------------------------------------------------
# settings
# ---------------------------------------------------------------------------

def _ch(ds: "Dataset") -> dict:
    return ds.chars.setdefault("_dta", {})


def mi_style(ds: "Dataset") -> str:
    ch = ds.chars.get("_dta", {})
    return ch.get("_mi_style", "") if ch.get("_mi_marker") else ""


def require_mi(s: "Session") -> str:
    st_ = mi_style(s.data)
    if not st_:
        raise StataError(119, "no; data not mi set")   # VERIFY
    return st_


def _list(ds, key: str) -> list[str]:
    return [v for v in _ch(ds).get(key, "").split() if ds.has(v)]


def _set_list(ds, key: str, names: list[str]) -> None:
    ch = _ch(ds)
    if names:
        ch[key] = " ".join(dict.fromkeys(names))
    else:
        ch.pop(key, None)


def nimp(ds) -> int:
    return int(_ch(ds).get("_mi_M", "0") or 0)


def _touch(ds) -> None:
    _ch(ds)["_mi_update"] = str(time.time())


def _num(ds, name) -> np.ndarray:
    return ds.get(name).data.astype(np.float64)


def _update_miss(ds) -> None:
    """_mi_miss: 1 in m=0 for observations with a missing imputed variable."""
    m = _num(ds, "_mi_m")
    ivars = _list(ds, "_mi_ivars")
    miss = np.zeros(ds.nobs, dtype=bool)
    for v in ivars:
        x = ds.get(v)
        miss |= (x.data == "") if x.is_string else (x.data >= SYS)
    vals = np.where(m == 0, miss.astype(float), SYS)
    ds.get("_mi_miss").data = vals


# ---------------------------------------------------------------------------
# the dataset of imputation m
# ---------------------------------------------------------------------------

def dataset_m(ds: "Dataset", m: int) -> "Dataset":
    """A copy of the data of imputation m (without the system variables)."""
    style = mi_style(ds)
    mm = _num(ds, "_mi_m")
    if style == "flong" or m == 0:
        rows = np.nonzero(mm == m)[0]
    else:
        ids = _num(ds, "_mi_id")
        base = np.nonzero(mm == 0)[0]
        imp = np.nonzero(mm == m)[0]
        repl = dict(zip(ids[imp].tolist(), imp.tolist()))
        rows = np.array([repl.get(i, r) for i, r in zip(ids[base].tolist(), base.tolist())], dtype=int)
    new = ds.copy()
    keep = np.zeros(ds.nobs, dtype=bool)
    keep[rows] = True
    # keep the rows in the order of m=0 (rows may come from the end of the data)
    order = rows
    for v in new.vars:
        v.raw = v.raw[order]
    new.nobs = len(order)
    new.drop_vars([v for v in SYSV if new.has(v)])
    new.esample = None
    ch = new.chars.get("_dta", {})
    for k in [k for k in ch if k.startswith("_mi_")]:
        ch.pop(k)
    return new


class _with_data:
    """Runs commands on another dataset and restores the original one."""

    def __init__(self, s: "Session", ds):
        self.s, self.ds = s, ds

    def __enter__(self):
        self.old = self.s.data
        self.s.data = self.ds
        return self.ds

    def __exit__(self, *exc):
        self.s.data = self.old
        return False


# ---------------------------------------------------------------------------
# mi dispatcher
# ---------------------------------------------------------------------------

def cmd_mi(s: "Session", args: str) -> None:
    from .mi_estimate import mi_estimate
    from .mi_impute import mi_impute
    text = args.strip()
    if not text:
        raise StataError(198, "subcommand required")   # VERIFY
    head = text.split(None, 1)[0].rstrip(":,")
    rest = text[len(text.split(None, 1)[0]):].strip() if " " in text else ""
    if text.split(None, 1)[0].endswith(":"):
        rest = ": " + rest
    table = [("set", 3, mi_set), ("register", 3, mi_register), ("unregister", 5, mi_unregister),
             ("describe", 1, mi_describe), ("impute", 3, mi_impute), ("estimate", 3, mi_estimate),
             ("passive", 3, mi_passive), ("xeq", 3, mi_xeq), ("extract", 3, mi_extract),
             ("misstable", 5, mi_misstable), ("query", 1, mi_query)]
    for name, k, fn in table:
        if len(head) >= k and name.startswith(head):
            if name in ("estimate", "passive", "xeq") and ":" in text:
                i = text.index(":")
                fn(s, text[len(text.split(None, 1)[0]):i].strip().lstrip(",").strip()
                   if len(text.split(None, 1)[0]) <= i else "", text[i + 1:].strip())
            elif name in ("estimate", "passive", "xeq"):
                raise StataError(198, "colon required")   # VERIFY
            else:
                fn(s, rest)
            return
    raise StataError(198, f"mi {head}: invalid subcommand")   # VERIFY


def mi_set(s: "Session", args: str) -> None:
    from ..core.dataset import Variable
    ds = s.data
    style = args.strip().lower()
    if style in ("m", "ml", "mlo", "mlon", "mlong"):
        style = "mlong"
    elif style in ("flong",):
        style = "flong"
    elif style in ("wide", "w", "flongsep"):
        raise StataError(198, f"mi set {style}: OpenDTA supports the mlong and flong styles")   # OpenDTA
    else:
        raise StataError(198, "mi set: style required (mlong or flong)")   # VERIFY
    if mi_style(ds):
        raise StataError(119, "data already mi set")   # VERIFY
    n = ds.nobs
    ds.add(Variable("_mi_m", "int", np.zeros(n)))
    ds.add(Variable("_mi_id", "long" if n > 32000 else "int", np.arange(1, n + 1, dtype=float)))
    ds.add(Variable("_mi_miss", "byte", np.zeros(n)))
    ch = _ch(ds)
    ch.update({"_mi_marker": "_mi_ds_1", "_mi_style": style, "_mi_M": "0"})
    _touch(ds)
    s.notify_state()


def _reg(s, args, key):
    ds = s.data
    require_mi(s)
    words = args.split()
    if not words:
        raise StataError(198, "invalid syntax")
    kind = words[0]
    names = s.expand_varlist(" ".join(words[1:]))
    full = next((k for k in ("imputed", "passive", "regular") if len(kind) >= 3 and k.startswith(kind)), None)
    if full is None:
        raise StataError(198, f"mi register {kind}: invalid")   # VERIFY
    return full, names


def mi_register(s: "Session", args: str) -> None:
    full, names = _reg(s, args, "")
    ds = s.data
    for nm in names:
        if nm in SYSV:
            raise StataError(198, f"{nm} is a system variable")   # VERIFY
    for key in ("_mi_ivars", "_mi_pvars", "_mi_rvars"):
        _set_list(ds, key, [v for v in _list(ds, key) if v not in names])
    key = {"imputed": "_mi_ivars", "passive": "_mi_pvars", "regular": "_mi_rvars"}[full]
    _set_list(ds, key, _list(ds, key) + names)
    if full == "imputed":
        before = int(np.sum(_num(ds, "_mi_miss") == 1))
        _update_miss(ds)
        now = int(np.sum(_num(ds, "_mi_miss") == 1)) - before
        if now > 0:
            s.output.write(f"({now} m=0 obs. now marked as incomplete)\n", "text")     # case 1909
    _touch(ds)


def mi_unregister(s: "Session", args: str) -> None:
    require_mi(s)
    ds = s.data
    names = s.expand_varlist(args)
    for key in ("_mi_ivars", "_mi_pvars", "_mi_rvars"):
        _set_list(ds, key, [v for v in _list(ds, key) if v not in names])
    _update_miss(ds)


def mi_query(s: "Session", args: str) -> None:
    ds = s.data
    style = mi_style(ds)
    if not style:
        s.output.write("(data not mi set)\n", "text")
        return
    s.output.write(f"data mi set {style}, M = {nimp(ds)}\n", "text")
    s.r = {"M": float(nimp(ds)), "style": style}


def mi_describe(s: "Session", args: str) -> None:
    style = require_mi(s)
    ds = s.data
    out = s.output
    m = _num(ds, "_mi_m")
    miss = _num(ds, "_mi_miss")
    m0 = m == 0
    n_inc = int(np.sum(miss[m0] == 1))
    n0 = int(m0.sum())
    M_ = nimp(ds)
    ivars, pvars, rvars = _list(ds, "_mi_ivars"), _list(ds, "_mi_pvars"), _list(ds, "_mi_rvars")
    out.write(f"\n  Style:  {style}\n", "text")
    ago = time.time() - float(_ch(ds).get("_mi_update", time.time()))
    stamp = time.strftime("%d%b%Y %H:%M:%S", time.localtime(time.time() - ago)).lower()
    out.write(f"          last mi update {stamp}, {int(ago)} seconds ago\n\n", "text")
    out.write(f"  Obs.:   complete    {n0 - n_inc:>9,}\n", "text")
    out.write(f"          incomplete  {n_inc:>9,}  (M = {M_} imputation{'s' if M_ != 1 else ''})\n", "text")
    out.write("          ---------------------\n", "text")
    out.write(f"          total       {n0:>9,}\n\n", "text")

    def lst(title, names, count=False):
        parts = []
        for v in names:
            if count:
                x = ds.get(v)
                k = int(np.sum(((x.data == "") if x.is_string else (x.data >= SYS))[m0]))
                parts.append(f"{v}({k})")
            else:
                parts.append(v)
        out.write(f"  {'Vars.:  ' if title == 'imputed' else '        '}{title + ':':<9}{len(names):>2}"
                  f"{'; ' + ' '.join(parts) if parts else ''}\n\n", "text")
    lst("imputed", ivars, True)
    lst("passive", pvars, True)      # with the missing counts, as Stata (case 1909)
    lst("regular", rvars)
    out.write(f"          system:   3; {' '.join(SYSV)}\n\n", "text")
    unreg = [v for v in ds.names if v not in ivars + pvars + rvars and v not in SYSV]
    if unreg:
        out.write(f"         (there {'are' if len(unreg) != 1 else 'is'} {len(unreg)} unregistered "
                  f"variable{'s' if len(unreg) != 1 else ''}; {' '.join(unreg)})\n", "text")
    s.r = {"M": float(M_), "N": float(n0), "N_incomplete": float(n_inc), "N_complete": float(n0 - n_inc),
           "ivars": " ".join(ivars)}


def mi_extract(s: "Session", args: str) -> None:
    require_mi(s)
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    m = int(float(body.strip()))
    if m < 0 or m > nimp(s.data):
        raise StataError(198, f"m={m} does not exist")   # VERIFY
    new = dataset_m(s.data, m)
    s.data = new
    new.changed = True
    s.notify_state()


def mi_misstable(s: "Session", args: str) -> None:
    require_mi(s)
    with _with_data(s, dataset_m(s.data, 0)):
        s.run("misstable " + args, echo=False)


def mi_xeq(s: "Session", mlist: str, cmd: str) -> None:
    require_mi(s)
    ds = s.data
    from ..lang.words import parse_numlist
    ms = [int(v) for v in parse_numlist(mlist)] if mlist.strip() else list(range(nimp(ds) + 1))
    out = s.output
    for m in ms:
        out.write(f"\nm={m} data:\n", "text")
        for c in cmd.split(";"):
            if not c.strip():
                continue
            out.write(f"-> {c.strip()}\n", "text")
            with _with_data(s, dataset_m(ds, m)):
                rc = s.run(c.strip(), echo=False)
            if rc:
                raise StataError(rc, "")


def mi_passive(s: "Session", opts: str, cmd: str) -> None:
    """mi passive: generate|replace|egen ... : runs the command on every imputation."""
    require_mi(s)
    ds = s.data
    words = cmd.split()
    if not words or words[0] not in ("generate", "gen", "g", "replace", "egen"):
        raise StataError(198, "mi passive: generate, replace or egen required")   # VERIFY
    target = words[1].split("=")[0]
    if target in ("byte", "int", "long", "float", "double") or target.startswith("str"):
        target = words[2].split("=")[0]
    out = s.output
    mm = _num(ds, "_mi_m")
    M_ = nimp(ds)
    style = mi_style(ds)
    rc = 0
    vals = None
    for m in range(M_ + 1):
        out.write(f"m={m}:\n", "text")
        sub = dataset_m(ds, m)
        if words[0] != "replace" and sub.has(target):
            sub.drop_vars([target])
        with _with_data(s, sub):
            rc = s.run(cmd, echo=False)
        if rc:
            raise StataError(rc, "")
        v = sub.get(target)
        if vals is None:
            from ..core.dataset import Variable
            if not ds.has(target):
                ds.add(Variable(target, v.vtype, np.full(ds.nobs, "" if v.is_string else SYS,
                                                         dtype=object if v.is_string else float)))
            vals = ds.get(target)
        # write back the rows of m
        if m == 0 or style == "flong":
            rows = np.nonzero(mm == m)[0]
            vals_arr = np.array(vals.data, copy=True)
            vals_arr[rows] = v.data
            vals.data = vals_arr
        else:
            ids = _num(ds, "_mi_id")
            imp = np.nonzero(mm == m)[0]
            base_ids = _num(sub, "_mi_id") if sub.has("_mi_id") else None
            pos0 = {int(i): k for k, i in enumerate(ids[mm == 0].tolist())}
            arr = np.array(vals.data, copy=True)
            for r in imp:
                arr[r] = v.data[pos0[int(ids[r])]]
            vals.data = arr
            del base_ids
    _set_list(ds, "_mi_pvars", _list(ds, "_mi_pvars") + [target])
    _touch(ds)
    s.notify_state()


# ---------------------------------------------------------------------------
# mi impute
# ---------------------------------------------------------------------------


_METHODS = {"regress": 3, "pmm": 3, "logit": 4, "ologit": 5, "mlogit": 5, "poisson": 4, "chained": 5,
            "monotone": 4, "mvn": 3, "truncreg": 5, "intreg": 4, "nbreg": 5}
