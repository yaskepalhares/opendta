"""bsample ([R] bsample): a bootstrap sample of the data in memory."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Variable
from ..core.errors import StataError
from ..core.varlist import expand
from ..lang.syntax import match_options, parse_standard
from ._util import eval_vector, touse
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def draw_units(n: int, m: int, u: np.ndarray) -> np.ndarray:
    """m indices drawn with replacement from 0..n-1, one uniform each."""
    return np.minimum((u[:m] * n).astype(np.int64), n - 1)


def _keys(s: "Session", names: list[str], rows: np.ndarray) -> np.ndarray:
    """Integer code of each row's combination of values."""
    if not names:
        return np.zeros(len(rows), dtype=np.int64)
    cols = [s.data.get(n).data[rows] for n in names]
    combo = [tuple(c[i] for c in cols) for i in range(len(rows))]
    codes = {k: j for j, k in enumerate(sorted(set(combo), key=str))}
    return np.array([codes[k] for k in combo], dtype=np.int64)


@command("bsample")
def cmd_bsample(s: "Session", args: str) -> None:
    """bsample [exp] [if] [in] [, strata() cluster() idcluster() weight()]. VERIFY: Stata's
    order of the draws (the sample differs from Stata's for the same seed), and that the
    observations outside if/in are dropped."""
    from ..core import rng
    p = parse_standard(args)
    o = match_options(p.options, {"strata": 6, "cluster": 7, "idcluster": 3, "weight": 6})
    ds = s.data
    rows = np.flatnonzero(touse(s, p)) if ds.nobs else np.arange(0)
    strata = expand(ds, str(o["strata"])) if o.get("strata") else []
    clusters = expand(ds, str(o["cluster"])) if o.get("cluster") else []
    if o.get("idcluster") and not clusters:
        raise StataError(198, "idcluster() requires cluster()")   # VERIFY
    st = _keys(s, strata, rows)
    cl = _keys(s, clusters, rows) if clusters else np.arange(len(rows))
    size = None
    if p.varlist.strip():
        v = float(eval_vector(s, p.varlist))
        if M.is_missing(v) or v < 1 or v != int(v):
            raise StataError(198, "resample size must be a positive integer")   # VERIFY
        size = int(v)
    picked: list[np.ndarray] = []
    ids: list[np.ndarray] = []
    next_id = 1
    for h in np.unique(st):
        in_h = np.flatnonzero(st == h)
        units = np.unique(cl[in_h])
        m = size if size is not None else len(units)
        if m > len(units):
            raise StataError(498, "resample size must not be greater than number of "
                                  + ("clusters" if clusters else "observations"))
        draw = units[draw_units(len(units), m, rng.RNG.uniform(m))]
        for unit in draw:
            members = in_h[cl[in_h] == unit]
            picked.append(rows[members])
            ids.append(np.full(len(members), next_id))
            next_id += 1
    idx = np.concatenate(picked) if picked else np.arange(0, dtype=np.int64)
    if o.get("weight"):
        wname = expand(ds, str(o["weight"]))[0]
        w = ds.get(wname)
        if w.is_string:
            raise StataError(109, "type mismatch")
        counts = np.zeros(ds.nobs)
        np.add.at(counts, idx, 1.0)
        ds.set_numeric(w, counts)
        s.notify_state()
        return
    for v in ds.vars:
        v.raw = v.raw[idx]
    ds.nobs = len(idx)
    ds.sortlist = []
    if o.get("idcluster"):
        name = str(o["idcluster"]).strip()
        if ds.has(name):
            ds.drop_vars([name])
        ds.add(Variable(name, "long", np.concatenate(ids).astype(float) if ids else np.zeros(0)))
    ds.changed = True
    s.notify_state()
