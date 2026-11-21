"""Pharmacokinetic measures and crossover analyses ([R] pk, pkexamine, pksumm, pkcollapse,
pkequiv, pkcross). Plain functions of arrays: no Session, no printing.

VERIFY: the manual does not give every formula; where it does not, the choices are the usual
noncompartmental ones and are marked."""

from __future__ import annotations

import math

import numpy as np

from ..lazy import lazy

st = lazy("scipy.stats")
interp = lazy("scipy.interpolate")


def auc(t: np.ndarray, c: np.ndarray, trapezoid: bool = False) -> float:
    """Area under the concentration curve from the first to the last time: trapezoidal rule, or
    the integral of a cubic spline through the points (the default). VERIFY: the spline Stata
    uses."""
    if len(t) < 2:
        return 0.0
    if trapezoid or len(t) < 3:
        return float(np.sum((c[1:] + c[:-1]) / 2 * np.diff(t)))
    return float(interp.CubicSpline(t, c, bc_type="natural").integrate(t[0], t[-1]))


def measures(t: np.ndarray, c: np.ndarray, fit: int = 3, trapezoid: bool = False) -> dict:
    """The measures of pkexamine: cmax, tomc (time of cmax), tmax (last time), auc on
    [0, tmax], the elimination rate ke and half-life from the regression of ln(c) on t over the
    last `fit` points, and three AUCs on [0, inf): adding c_last/ke (auc_ln), the area under the
    linear fit of c on t (auc_line) and under the exponential fit (auc_exp)."""
    o = np.argsort(t, kind="stable")
    t, c = np.asarray(t, dtype=float)[o], np.asarray(c, dtype=float)[o]
    if t[0] > 0:
        t, c = np.concatenate([[0.0], t]), np.concatenate([[0.0], c])   # VERIFY
    i = int(np.argmax(c))
    out = {"cmax": float(c[i]), "tomc": float(t[i]), "tmax": float(t[-1]),
           "auc": auc(t, c, trapezoid)}
    k = min(fit, len(t))
    tt, cc = t[-k:], c[-k:]
    ok = cc > 0
    ke = half = auc_ln = auc_exp = auc_line = math.nan
    if ok.sum() >= 2:
        b, a = np.polyfit(tt[ok], np.log(cc[ok]), 1)
        if b < 0:
            ke = -b
            half = math.log(2) / ke
            auc_ln = out["auc"] + c[-1] / ke
            auc_exp = out["auc"] + math.exp(a + b * t[-1]) / ke
    if len(tt) >= 2:
        b1, a1 = np.polyfit(tt, cc, 1)
        if b1 < 0:
            c_end = a1 + b1 * t[-1]
            t0 = -a1 / b1
            auc_line = out["auc"] + max(c_end, 0.0) * max(t0 - t[-1], 0.0) / 2
    out.update({"ke": ke, "half": half, "auc_ln": auc_ln, "auc_exp": auc_exp, "auc_line": auc_line})
    return out


def summary(x: np.ndarray) -> dict:
    """Mean, median, variance, skewness, kurtosis and the p-value of the Shapiro-Wilk test of
    a measure across subjects (pksumm). VERIFY: the normality test pksumm reports."""
    x = x[np.isfinite(x)]
    n = len(x)
    d = x - x.mean()
    m2 = (d ** 2).mean()
    skew = (d ** 3).mean() / m2 ** 1.5 if m2 > 0 else math.nan
    kurt = (d ** 4).mean() / m2 ** 2 if m2 > 0 else math.nan
    p = float(st.shapiro(x).pvalue) if n >= 3 and m2 > 0 else math.nan
    return {"mean": float(x.mean()), "median": float(np.median(x)),
            "var": float(x.var(ddof=1)) if n > 1 else math.nan, "skew": skew, "kurt": kurt, "p": p}


def equivalence(y: np.ndarray, treat: np.ndarray, period: np.ndarray, seq: np.ndarray, ident: np.ndarray,
                test, ref, level: float = 90.0, limit: float = 0.2) -> dict:
    """Classic confidence interval for average bioequivalence in a two-period crossover: the
    difference test - reference of the treatment means, with the residual variance of the
    crossover ANOVA (sequence, subject within sequence, period, treatment); the limits are
    +/- limit times the reference mean. VERIFY: the degrees of freedom and the probability
    reported by pkequiv."""
    yt, yr = y[treat == test], y[treat == ref]
    mt, mr = yt.mean(), yr.mean()
    a = anova_crossover(y, seq, ident, period, treat, None)
    mse, df = a["mse"], a["df_r"]
    n_t, n_r = len(yt), len(yr)
    se = math.sqrt(mse * (1 / n_t + 1 / n_r))
    tcrit = float(st.t.ppf(1 - (1 - level / 100) / 2, df))
    diff = mt - mr
    lo, hi = diff - tcrit * se, diff + tcrit * se
    eq = limit * mr
    p_in = float(st.t.cdf((eq - diff) / se, df) - st.t.cdf((-eq - diff) / se, df))
    return {"diff": diff, "lo": lo, "hi": hi, "eq_lo": -eq, "eq_hi": eq, "ratio": mt / mr,
            "rlo": 1 + lo / mr, "rhi": 1 + hi / mr, "p": p_in, "se": se, "df": df,
            "mean_t": mt, "mean_r": mr}


def _dummies(x: np.ndarray) -> np.ndarray:
    levels = np.unique(x)
    return np.column_stack([(x == lv).astype(float) for lv in levels[1:]]) if len(levels) > 1 \
        else np.zeros((len(x), 0))


def anova_crossover(y, seq, ident, period, treat, carry) -> dict:
    """Sequential sums of squares of the crossover model y = sequence + id|sequence + period +
    treatment [+ carryover] (pkcross). The id terms are the subjects within sequence."""
    n = len(y)
    terms = [("sequence", _dummies(seq)), ("id|sequence", _dummies(ident)), ("period", _dummies(period)),
             ("treatment", _dummies(treat))]
    if carry is not None:
        terms.append(("carryover", _dummies(carry)))
    X = np.ones((n, 1))
    rss_prev = float(((y - y.mean()) ** 2).sum())
    rank_prev = 1
    out_terms = []
    for name, D in terms:
        Xn = np.column_stack([X, D]) if D.shape[1] else X
        beta, *_ = np.linalg.lstsq(Xn, y, rcond=None)
        e = y - Xn @ beta
        rss = float(e @ e)
        rank = int(np.linalg.matrix_rank(Xn))
        out_terms.append((name, rss_prev - rss, rank - rank_prev))
        X, rss_prev, rank_prev = Xn, rss, rank
    df_r = n - rank_prev
    mse = rss_prev / df_r if df_r > 0 else math.nan
    rows = []
    for name, ss, df in out_terms:
        ms = ss / df if df else math.nan
        F = ms / mse if df and mse > 0 else math.nan
        rows.append({"term": name, "ss": ss, "df": df, "ms": ms, "F": F,
                     "p": float(st.f.sf(F, df, df_r)) if df and math.isfinite(F) else math.nan})
    # the sequence effect is tested against the subjects within sequence (between-subject error)
    sub = rows[1]
    if rows[0]["df"] and sub["df"]:
        rows[0]["F"] = rows[0]["ms"] / sub["ms"] if sub["ms"] > 0 else math.nan
        rows[0]["p"] = float(st.f.sf(rows[0]["F"], rows[0]["df"], sub["df"])) \
            if math.isfinite(rows[0]["F"]) else math.nan
    return {"rows": rows, "rss": rss_prev, "df_r": df_r, "mse": mse,
            "tss": float(((y - y.mean()) ** 2).sum()), "n": n}
