"""Groups defined by key variables (by(), collapse, contract, merge...).

`group_ids` numbers the groups 0..k-1 in increasing key order, as
sort would leave them: numbers before missing (. < .a < ... < .z) and strings in
byte order. Observations outside the mask get -1.
"""

from __future__ import annotations

import numpy as np

from .dataset import Dataset


def key_columns(ds: Dataset, names: list[str]) -> list[np.ndarray]:
    """float64 columns that sort like the variables (strings become ranks)."""
    cols = []
    for name in names:
        v = ds.get(name)
        if v.is_string:
            _, inv = np.unique(v.raw.astype(str), return_inverse=True)
            cols.append(inv.astype(np.float64).ravel())
        else:
            cols.append(np.asarray(v.data, dtype=np.float64))
    return cols


def group_ids(ds: Dataset, names: list[str], mask: np.ndarray | None = None
              ) -> tuple[np.ndarray, int]:
    """(group id per observation, number of groups)."""
    n = ds.nobs
    mask = np.ones(n, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    ids = np.full(n, -1, dtype=np.int64)
    if not mask.any():
        return ids, 0
    if not names:
        ids[mask] = 0
        return ids, 1
    inv = _dense_ids([c[mask] for c in key_columns(ds, names)])
    ids[mask] = inv
    return ids, int(inv.max()) + 1


def _dense_ids(cols: list[np.ndarray]) -> np.ndarray:
    """0..k-1 in increasing key order. One key: np.unique; several: a lexicographic sort and
    the points where any key changes (np.unique(axis=0) sorts rows as structured records and
    is many times slower)."""
    if len(cols) == 1:
        _, inv = np.unique(cols[0], return_inverse=True)
        return inv.ravel().astype(np.int64)
    m = len(cols[0])
    order = np.lexsort(cols[::-1])
    change = np.zeros(m, dtype=bool)
    if m:
        change[0] = True
    for c in cols:
        cs = c[order]
        change[1:] |= cs[1:] != cs[:-1]
    inv = np.empty(m, dtype=np.int64)
    inv[order] = np.cumsum(change) - 1
    return inv


def first_index(ids: np.ndarray, k: int) -> np.ndarray:
    """Index of the first observation of each group (in data order)."""
    first = np.full(k, -1, dtype=np.int64)
    sel = np.flatnonzero(ids >= 0)
    # walks backward: the last write is the first obs.
    first[ids[sel[::-1]]] = sel[::-1]
    return first


def joint_keys(a: list[np.ndarray], b: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray, int]:
    """Common ids for the keys of two datasets (merge, joinby).
    Each list has one column per key variable, already comparable with each other
    (numbers as float64; strings as str objects)."""
    na = len(a[0]) if a else 0
    cols = []
    for ca, cb in zip(a, b):
        if ca.dtype == object or cb.dtype == object:
            both = np.concatenate([ca.astype(str), cb.astype(str)])
            _, inv = np.unique(both, return_inverse=True)
            cols.append(inv.astype(np.float64).ravel())
        else:
            cols.append(np.concatenate([ca, cb]).astype(np.float64))
    if not cols:
        return np.zeros(na, dtype=np.int64), np.zeros(len(b[0]) if b else 0, dtype=np.int64), 1
    if len(cols[0]) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), 0
    inv = _dense_ids(cols)
    return inv[:na], inv[na:], int(inv.max()) + 1
