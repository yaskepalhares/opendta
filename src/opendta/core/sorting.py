"""Sorting of observations (sort, gsort, bysort).

Missing values go after all numbers (. < .a < ... < .z), which is already the
natural order of the codes used. Sorting is stable. Stata does not guarantee
the order of ties without the `stable` option; here it is always stable.
"""

from __future__ import annotations

import numpy as np

from .dataset import Dataset
from .missing import SYSMISS


def _rank(data: np.ndarray) -> np.ndarray:
    if data.dtype == object:
        _, inv = np.unique(data.astype(str), return_inverse=True)
        return inv.astype(np.float64)
    return data.astype(np.float64)


def sort_order(ds: Dataset, keys: list[str], descending: list[bool] | None = None,
               *, missing_first: bool = False) -> np.ndarray:
    descending = descending or [False] * len(keys)
    cols = []
    for name, desc in zip(keys, descending):
        var = ds.get(name)
        r = _rank(var.data)
        if desc:
            if var.is_string:
                r = -r
            else:
                miss = r >= SYSMISS
                # gsort -x: missing values go last, unless mfirst
                r = np.where(miss, (-np.inf if missing_first else np.inf), -r)
        cols.append(r)
    if not cols:
        return np.arange(ds.nobs)
    return np.lexsort(cols[::-1])


def is_sorted(ds: Dataset, keys: list[str]) -> bool:
    order = sort_order(ds, keys)
    return bool(np.all(order == np.arange(ds.nobs)))


def sort(ds: Dataset, keys: list[str], descending: list[bool] | None = None, **kw) -> None:
    order = sort_order(ds, keys, descending, **kw)
    if not np.all(order == np.arange(ds.nobs)):
        ds.reorder_obs(order)
    ds.sortlist = list(keys) if not descending or not any(descending) else []
