"""Stochastic frontier for panel data, Battese and Coelli (1988, 1992) ([XT] xtfrontier, Methods and
formulas).

Plain functions on arrays: no Session, no variables, no printing. The composed error of the
production frontier is e_it = v_it - u_it, with v_it ~ N(0, sv2) and u_it = d_it u_i, where u_i is
N(mu, su2) truncated at zero and d_it = exp(-eta (t - T_i)) (d_it = 1 for the time-invariant
model). For a cost frontier the sign of u is reversed."""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy
from .panelquad import Groups

special = lazy("scipy.special")

LOG_2PI = math.log(2 * math.pi)


def decay(tidx: np.ndarray, grp: Groups, eta: float) -> np.ndarray:
    """d_it = exp(-eta (t - T_i)), with T_i the last period of the panel."""
    last = np.full(grp.G, -np.inf)
    np.maximum.at(last, grp.inv, tidx.astype(float))
    return np.exp(-eta * (tidx - last[grp.inv]))


def group_loglik(e: np.ndarray, d: np.ndarray, grp: Groups, mu: float, su2: float, sv2: float,
                 cost: bool = False) -> np.ndarray:
    """Log likelihood of every panel: e are the residuals y - xb and d the inefficiency
    multipliers. The panel's u_i is integrated analytically."""
    ds = -d if cost else d
    T = grp.counts.astype(float)
    sds2 = grp.sums(ds * ds)
    sde = grp.sums(ds * e)
    see = grp.sums(e * e)
    A = sds2 / sv2 + 1.0 / su2
    B = sde / sv2 - mu / su2
    m = -B / A
    return (-0.5 * T * (LOG_2PI + math.log(sv2)) - 0.5 * (LOG_2PI + math.log(su2))
            - special.log_ndtr(mu / math.sqrt(su2)) + 0.5 * (LOG_2PI - np.log(A))
            + special.log_ndtr(m * np.sqrt(A)) + B ** 2 / (2 * A) - see / (2 * sv2) - mu ** 2 / (2 * su2))
