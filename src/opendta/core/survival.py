"""Plain functions for survival-time data manipulation ([ST] stsplit, stjoin, stfill, stgen, stvary,
snapspan): they take arrays and return arrays, without a Session, variables or printing."""

from __future__ import annotations

import numpy as np


def factorize(keys: np.ndarray) -> tuple[np.ndarray, int]:
    """Subject number (0, 1, ...) of each element of `keys` (numbers or strings), by sorted value."""
    _u, inv = np.unique(keys, return_inverse=True)
    return inv.astype(np.int64), len(_u)


def time_order(gid: np.ndarray, t0: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Positions that sort the records by subject, then entry time, then exit time (stable)."""
    return np.lexsort((t, t0, gid))


def group_bounds(gid_sorted: np.ndarray) -> np.ndarray:
    """Start position of each run of equal subject numbers in a sorted vector (plus the end)."""
    if len(gid_sorted) == 0:
        return np.zeros(1, dtype=np.int64)
    cut = np.flatnonzero(np.diff(gid_sorted)) + 1
    return np.concatenate([[0], cut, [len(gid_sorted)]]).astype(np.int64)


# ---------------------------------------------------------------------------
# stsplit
# ---------------------------------------------------------------------------

def split_pieces(t0: np.ndarray, t: np.ndarray, cuts: list[np.ndarray]):
    """Splits each record (t0, t] at the cut points of that record that fall strictly inside it.

    `cuts[i]` is the sorted vector of cut points of record i. Returns (parent, a, b, last): the
    record each piece comes from, its limits (a, b] and whether it is the final piece of the record
    (the one that keeps the failure indicator)."""
    parent: list[int] = []
    lo: list[float] = []
    hi: list[float] = []
    last: list[bool] = []
    for i in range(len(t)):
        c = cuts[i]
        inside = c[(c > t0[i]) & (c < t[i])] if len(c) else c
        edges = np.concatenate([[t0[i]], inside, [t[i]]])
        k = len(edges) - 1
        parent.extend([i] * k)
        lo.extend(edges[:-1].tolist())
        hi.extend(edges[1:].tolist())
        last.extend([False] * (k - 1) + [True])
    return (np.array(parent, dtype=np.int64), np.array(lo), np.array(hi), np.array(last, dtype=bool))


def category_labels(b: np.ndarray, cuts: list[np.ndarray], parent: np.ndarray, at: np.ndarray
                    ) -> np.ndarray:
    """The lower bound of the interval of `at` that contains each piece: pieces end at or before
    the next cut point, so the number of cut points strictly below the end of the piece picks the
    interval. The first interval is labelled 0, or one less than the smallest at() value when some
    value is negative ([ST] stsplit)."""
    first = 0.0 if not len(at) or at.min() >= 0 else float(at.min()) - 1
    labels = np.concatenate([[first], at])
    out = np.empty(len(b))
    for j in range(len(b)):
        c = cuts[parent[j]]
        out[j] = labels[int(np.sum(c < b[j]))] if len(c) else labels[0]
    return out


# ---------------------------------------------------------------------------
# stjoin
# ---------------------------------------------------------------------------

def join_chains(order: np.ndarray, gid: np.ndarray, t0: np.ndarray, t: np.ndarray,
                censored: np.ndarray, same: np.ndarray) -> np.ndarray:
    """For records in `order` (sorted by subject and time), True where the record is joined to the
    one before it: same subject, it starts where the previous one ended, the previous one ended
    without an event and `same` (record i has the same values as the one before it in `order`)."""
    n = len(order)
    join = np.zeros(n, dtype=bool)
    if n < 2:
        return join
    g, a, b, c = gid[order], t0[order], t[order], censored[order]
    join[1:] = (g[1:] == g[:-1]) & (a[1:] == b[:-1]) & c[:-1] & same
    return join


# ---------------------------------------------------------------------------
# stfill, stvary, stgen
# ---------------------------------------------------------------------------

def is_missing_values(x: np.ndarray, is_string: bool, sysmiss: float) -> np.ndarray:
    if is_string:
        return np.array([v == "" for v in x], dtype=bool)
    return x >= sysmiss


def fill_forward(values: np.ndarray, miss: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    """Positions of the source value of each element after carrying the last nonmissing value
    forward within each subject (-1 where nothing has been observed yet)."""
    n = len(values)
    src = np.full(n, -1, dtype=np.int64)
    for a, b in zip(bounds[:-1], bounds[1:]):
        cur = -1
        for k in range(a, b):
            if not miss[k]:
                cur = k
            src[k] = cur
    return src


def subject_summary(values: np.ndarray, miss: np.ndarray, bounds: np.ndarray):
    """Per subject: (constant, varying, never missing, always missing, sometimes missing) as in
    stvary. A subject whose values are all missing is neither constant nor varying."""
    cons = vary = never = always = some = 0
    for a, b in zip(bounds[:-1], bounds[1:]):
        m = miss[a:b]
        ok = values[a:b][~m]
        if m.all():
            always += 1
        elif m.any():
            some += 1
        else:
            never += 1
        if len(ok):
            if len(set(ok.tolist())) == 1:
                cons += 1
            else:
                vary += 1
    return cons, vary, never, always, some


def cumulative_before(flag: np.ndarray, bounds: np.ndarray, include_current: bool) -> np.ndarray:
    """Number of records with flag set up to the previous record (or up to this one) of the same
    subject: count() and count0() of stgen."""
    out = np.zeros(len(flag))
    for a, b in zip(bounds[:-1], bounds[1:]):
        c = np.cumsum(flag[a:b].astype(float))
        out[a:b] = c if include_current else c - flag[a:b]
    return out


def gaps(t0: np.ndarray, t: np.ndarray, bounds: np.ndarray):
    """For records sorted by subject and time: the number of gaps before each record and their
    total length. Delayed entry (the first record) is not a gap."""
    n = len(t)
    ng = np.zeros(n)
    gl = np.zeros(n)
    for a, b in zip(bounds[:-1], bounds[1:]):
        cnt, tot = 0.0, 0.0
        for k in range(a, b):
            if k > a and t0[k] > t[k - 1]:
                cnt += 1
                tot += t0[k] - t[k - 1]
            ng[k], gl[k] = cnt, tot
    return ng, gl
