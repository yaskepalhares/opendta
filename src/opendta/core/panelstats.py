"""Panel descriptive tabulation and data transformations ([XT] xttab, [XT] xtdata).

Plain functions on arrays: no Session, no variables, no printing."""

from __future__ import annotations

import numpy as np


def xttab_counts(codes: np.ndarray, inv: np.ndarray, ncat: int) -> dict:
    """Overall, between and within figures of a categorical variable. `codes` is the category
    (0..ncat-1) of every observation and `inv` its panel (0..G-1).

    freq: observations per category; bfreq: panels with at least one observation in it;
    within: mean over those panels of the share (percent) of their observations in the category."""
    G = int(inv.max()) + 1
    freq = np.bincount(codes, minlength=ncat).astype(float)
    tot_i = np.bincount(inv, minlength=G).astype(float)
    cell = np.zeros((G, ncat))
    np.add.at(cell, (inv, codes), 1.0)
    bfreq = (cell > 0).sum(axis=0).astype(float)
    share = 100.0 * cell / tot_i[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        within = np.where(bfreq > 0, (share * (cell > 0)).sum(axis=0) / bfreq, 0.0)
    return {"freq": freq, "bfreq": bfreq, "within": within, "G": G, "N": float(len(codes))}


def panel_transform(x: np.ndarray, inv: np.ndarray, mode: str) -> np.ndarray:
    """fe: x - mean_i + overall mean; be: mean_i."""
    G = int(inv.max()) + 1
    cnt = np.bincount(inv, minlength=G)
    means = np.bincount(inv, weights=x, minlength=G) / cnt
    if mode == "be":
        return means[inv]
    return x - means[inv] + x.mean()
