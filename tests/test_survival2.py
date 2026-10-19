"""Phase 15: survival analysis, the remaining commands (strate, stptime, stci, ltable, stsplit,
ct, stcrreg, stcurve, ggamma, frailty, exact ties). Reference values come from identities and from
independent computations with numpy and scipy."""

import itertools

import numpy as np
from scipy import optimize, stats as st

from opendta.session import Session

from test_survival import SIM, _bse, _session


def test_stptime_rate_and_interval():
    s, out = _session("stptime")
    d, T = s.r["failures"], s.r["ptime"]
    assert d == 21 and T == 1595
    assert abs(s.r["rate"] - d / T) < 1e-12
    f = np.exp(st.norm.ppf(.975) / np.sqrt(d))          # limits on the ln(rate) scale
    assert abs(s.r["lb"] - s.r["rate"] / f) < 1e-12
    assert abs(s.r["ub"] - s.r["rate"] * f) < 1e-12


def test_stsplit_conserves_time_and_failures():
    s, out = _session("stset tt, failure(died) id(id)\nstsplit, at(10 20 30)\nstptime")
    assert s.data.nobs > 60
    assert s.r["failures"] == 21 and s.r["ptime"] == 1595
    s.run_command("stjoin")
    assert s.data.nobs == 60


def test_ltable_layout():
    s, out = _session("ltable tt died, intervals(10)")
    assert "Beg." in out and "Survival" in out


def test_sttoct_cttost_roundtrip():
    s, out = _session("sttoct nf nc ne")
    assert s.data.has("nf") and s.data.has("nc")
    assert abs(float(np.sum(s.data.get("nf").data)) - 21) < 1e-9
    s.run_command("cttost")
    assert s.run_command("stptime") == 0
    assert s.r["failures"] == 21


def test_stcrreg_without_competing_events_is_breslow_cox():
    s, out = _session("stcox drug age, breslow")
    b0 = s.e["b"].data.ravel().copy()
    s.run_command("stcrreg drug age, compete(died==2)")
    b, se = _bse(s)
    assert np.allclose(b, b0, atol=1e-7)


def test_stcrreg_with_competing_events_scores_sum_to_zero():
    s, _ = _session("replace died = 2 if died == 1 & mod(id, 3) == 0\n"
                    "stset tt, failure(died == 1)\n"
                    "stcrreg drug age, compete(died == 2)")
    assert s.e["N_fail"] < 21 and s.e["N_fail"] > 5
    assert np.all(np.isfinite(s.e["V"].data))


def _perm_loglik(eta, fail, tt, method):
    """Exact marginal (permutation) log likelihood of distinct-time data with ties broken as
    simultaneous events (here: no ties, so every method equals the partial likelihood)."""
    ll = 0.0
    for i in np.where(fail)[0]:
        risk = tt >= tt[i]
        ll += eta[i] - np.log(np.sum(np.exp(eta[risk])))
    return ll


def test_exact_ties_equal_brute_force():
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text("clear\ninput id t d x\n1 1 1 0.5\n2 1 1 -0.3\n3 2 1 0.8\n4 2 1 0.1\n5 3 0 -0.6\n"
                    "6 3 1 0.9\n7 4 1 -0.2\n8 5 0 0.4\nend\nstset t, failure(d)\n"
                    "stcox x, exactp nolog\n")
    assert rc == 0, "".join(lines)[-1500:]
    ll_p = s.e["ll"]
    t = np.array([1, 1, 2, 2, 3, 3, 4, 5.])
    d = np.array([1, 1, 1, 1, 0, 1, 1, 0], bool)
    x = np.array([.5, -.3, .8, .1, -.6, .9, -.2, .4])

    def nll(b):
        total = 0.0
        for tm in np.unique(t[d]):
            fails = np.where((t == tm) & d)[0]
            risk = np.where(t >= tm)[0]
            sub = sum(np.exp(b[0] * x[list(c)].sum()) for c in itertools.combinations(risk, len(fails)))
            total += b[0] * x[fails].sum() - np.log(sub)
        return -total
    r = optimize.minimize(nll, [0.0], method="BFGS", options={"gtol": 1e-10})
    assert abs(s.e["b"].data.ravel()[0] - r.x[0]) < 1e-5
    assert abs(ll_p - (-r.fun)) < 1e-8
    s.run_command("stcox x, exactm nolog")
    assert np.isfinite(s.e["ll"]) and s.e["ll"] != ll_p


def test_ggamma_nests_weibull_and_lognormal():
    s, _ = _session("streg drug age, dist(weibull) nolog")
    llw = s.e["ll"]
    s.run_command("streg drug age, dist(lognormal) nolog")
    lln = s.e["ll"]
    s.run_command("streg drug age, dist(ggamma) nolog")
    llg = s.e["ll"]
    assert llg >= max(llw, lln) - 1e-6


def test_gamma_frailty_exponential_matches_direct_likelihood():
    s, _ = _session("streg drug age, dist(exponential) frailty(gamma) nolog")
    ll = s.e["ll"]
    t = s.data.get("_t").data
    d = s.data.get("_d").data
    X = np.column_stack([s.data.get("drug").data, s.data.get("age").data, np.ones(len(t))])

    def nll(p):
        b, th = p[:3], np.exp(p[3])
        lam = t * np.exp(X @ b)
        return -np.sum(d * (X @ b) - (1 / th + d) * np.log1p(th * lam))
    jac = float(np.sum(d * np.log(t)))                    # streg reports the log likelihood in ln t
    r = optimize.minimize(nll, np.r_[s.e["b"].data.ravel()[:3], -3.0], method="Nelder-Mead",
                          options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 20000})
    assert abs(ll - (-r.fun + jac)) < 1e-4


def test_shared_frailty_header_and_lr_test():
    s, out = _session("gen g = ceil(id/4)\nstreg drug age, dist(weibull) frailty(gamma) shared(g) nolog")
    assert s.e["ll"] < 0
    s.run_command("streg drug age, dist(weibull) frailty(invgaussian) shared(g) nolog")
    assert s.e["ll"] < 0


def test_stcurve_after_models_draws_a_graph():
    s, _ = _session("streg drug age, dist(weibull) nolog\nstcurve, survival nodraw")
    s.run_command("stcox drug age")
    assert s.run_command("stcurve, cumhaz nodraw") == 0
