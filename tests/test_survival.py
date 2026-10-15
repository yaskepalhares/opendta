"""Phase 10: survival-time data (stset, sts, stcox, streg).

Reference values for Cox, Kaplan-Meier and the log-rank and Wilcoxon tests
come from statsmodels 0.15 (PHReg, SurvfuncRight, survdiff) on the same
data; the parametric models are checked against an independent
maximization of the likelihood written with scipy.stats."""

import numpy as np
import pytest
from scipy import optimize, stats as st

from opendta.session import Session

SIM = """clear
set obs 60
set seed 1234
gen id = _n
gen drug = 1 + (runiform() < .5)
gen age = 40 + floor(20*runiform())
gen double tt = ceil(-ln(runiform())/(0.02*exp(0.6*(drug==1) + 0.03*(age-50))))
gen died = runiform() < .75
replace tt = 40 if tt > 40
replace died = 0 if tt == 40
stset tt, failure(died)
"""


def _session(extra: str = ""):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(SIM + extra)
    assert rc == 0, "".join(lines)[-2000:]
    return s, "".join(lines)


def _bse(s):
    b = s.e["b"].data.ravel()
    V = s.e["V"].data
    return b, np.sqrt(np.diag(V))


def test_stset_variables_and_report():
    s, out = _session()
    ds = s.data
    assert ds.has("_st") and ds.has("_d") and ds.has("_t") and ds.has("_t0")
    assert ds.get("_t").vtype == "double"           # the time variable is double
    assert "21  failures in single-record/single-failure data" in out
    assert "1595  total analysis time at risk and under observation" in out
    assert np.all(ds.get("_t0").data == 0)


def test_cox_breslow_efron_robust():
    s, out = _session("stcox drug age, nohr")
    b, se = _bse(s)
    assert np.allclose(b, [-0.09546213, 0.04079989], atol=1e-7)
    assert np.allclose(se, [0.43848997, 0.03959882], rtol=1e-6)
    assert abs(s.e["ll"] - (-78.12700055261848)) < 1e-8
    assert "Cox regression -- Breslow method for ties" in out
    s.run_command("stcox drug age, efron")
    b, se = _bse(s)
    assert np.allclose(b, [-0.09338175, 0.04064022], atol=1e-7)
    assert np.allclose(se, [0.43850156, 0.03962899], rtol=1e-6)
    assert abs(s.e["ll"] - (-78.06471393075837)) < 1e-8
    s.run_command("stcox drug age, vce(robust)")
    b, se = _bse(s)
    # statsmodels' sandwich without Stata's small-sample factor n/(n-1)
    n = s.e["N"]
    assert np.allclose(np.exp(b) * se, np.array([0.38938351, 0.04435619]) * np.sqrt(n / (n - 1)),
                       rtol=1e-6)


def test_logrank_and_wilcoxon():
    s, out = _session("sts test drug")
    assert abs(s.r["chi2"] - 0.06973395780522201) < 1e-10
    s.run_command("sts test drug, wilcoxon")
    assert abs(s.r["chi2"] - 0.016542568246011188) < 1e-10


def test_km_generate_and_ci():
    s, _ = _session("sts generate km = s se = se(s) lb = lb(s) ub = ub(s), by(drug)")
    ds = s.data
    t = ds.get("_t").data
    g = ds.get("drug").data
    i = np.nonzero((t == 4) & (g == 1))[0][0]
    assert abs(ds.get("km").data[i] - 0.966667) < 1e-6
    assert abs(ds.get("se").data[i] - 0.032773) < 1e-6
    # ln(-ln S) interval
    S, sg = 29 / 30, np.sqrt(1 / (30 * 29)) / abs(np.log(29 / 30))
    z = st.norm.ppf(.975)
    assert abs(ds.get("lb").data[i] - S ** np.exp(z * sg)) < 1e-6
    assert abs(ds.get("ub").data[i] - S ** np.exp(-z * sg)) < 1e-6


def _weibull_ll(theta, t, d, X):
    b, lnp = theta[:-1], theta[-1]
    p = np.exp(lnp)
    xb = X @ b
    return -np.sum(d * (xb + lnp + (p - 1) * np.log(t)) - np.exp(xb) * t ** p)


def test_streg_weibull_matches_direct_maximization():
    s, out = _session("streg drug age, distribution(weibull)")
    ds = s.data
    t = ds.get("_t").data
    d = ds.get("_d").data
    X = np.column_stack([ds.get("drug").data, ds.get("age").data, np.ones(len(t))])
    r = optimize.minimize(_weibull_ll, np.array([0, 0, -5, 0.0]), args=(t, d, X), method="BFGS",
                          options={"gtol": 1e-10})
    b, _ = _bse(s)
    assert np.allclose(b, r.x, atol=1e-5)
    # Stata reports the log likelihood of ln t: sum(d ln t) more than that of t
    assert abs(s.e["ll"] - np.sum(d * np.log(t)) + r.fun) < 1e-8
    assert "Weibull regression -- log relative-hazard form" in out
    # AFT: same likelihood, coefficients -b/p
    s.run_command("streg drug age, distribution(weibull) time")
    b2, _ = _bse(s)
    p = np.exp(r.x[-1])
    assert np.allclose(b2[:3], -r.x[:3] / p, atol=1e-5)


@pytest.mark.parametrize("dist", ["exponential", "gompertz", "lognormal", "loglogistic"])
def test_streg_other_distributions_converge(dist):
    s, out = _session(f"streg drug age, distribution({dist})")
    assert s.e["converged"] == 1
    # Stata's convergence criterion: g V g' below nrtolerance()
    g = s.e["gradient"].data.ravel()
    assert float(g @ s.e["V"].data @ g) < 1e-5


def test_streg_lognormal_matches_direct_maximization():
    s, _ = _session("streg drug age, distribution(lognormal)")
    ds = s.data
    t = ds.get("_t").data
    d = ds.get("_d").data
    X = np.column_stack([ds.get("drug").data, ds.get("age").data, np.ones(len(t))])

    def nll(th):
        sig = np.exp(th[-1])
        z = (np.log(t) - X @ th[:-1]) / sig
        return -np.sum(np.where(d > 0, st.norm.logpdf(z) - np.log(sig) - np.log(t), st.norm.logsf(z)))
    r = optimize.minimize(nll, np.array([0, 0, 3, 0.0]), method="BFGS", options={"gtol": 1e-10})
    b, _ = _bse(s)
    assert np.allclose(b, r.x, atol=1e-5)


def test_multi_record_stset():
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    s.run_text("""clear
input id t0 t1 died
1 0 5 0
1 5 9 1
2 0 4 0
2 4 12 0
3 2 6 1
4 0 3 0
4 3 3 1
5 0 . 1
6 1 8 1
7 0 10 1
7 10 14 1
end
stset t1, id(id) failure(died) time0(t0)""")
    out = "".join(lines)
    assert "observation begins on or after (first) failure" in out
    assert "6  subjects" in out and "4  failures in single-failure-per-subject data" in out
    st_ = s.data.get("_st").data
    assert st_.tolist() == [1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0]
    s.run_command("stset t1, id(id) failure(died) time0(t0) exit(time .)")
    assert s.data.get("_st").data[10] == 1


def test_phtest_and_predict():
    s, out = _session("stcox drug age\nestat phtest, detail\npredict cs, csnell\npredict hr")
    assert "global test" in out
    ds = s.data
    # sum of the Cox-Snell residuals equals the number of failures (Breslow estimator)
    assert abs(np.sum(ds.get("cs").data) - 21) < 1e-4
    assert np.allclose(ds.get("hr").data, np.exp(s.e["b"].data.ravel() @
                                                 np.vstack([ds.get("drug").data, ds.get("age").data])),
                       rtol=1e-6)


def test_not_st():
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    s.run_text("set obs 3\ngen t = _n\ncapture stcox t\ndisplay _rc")
    assert "119" in "".join(lines)


def test_stset_exclusion_lines():
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    s.run_text("""clear
input id t1 died
1 5 0
1 9 1
4 3 0
4 3 1
5 . 1
7 10 1
7 14 1
end
stset t1, id(id) failure(died)""")
    out = "".join(lines)
    # both records at the same instant are excluded; the reason continues on a second line
    assert "          2  multiple records at same instant" + " " * 19 + "PROBABLE ERROR\n" in out
    assert "             (t1[_n-1]==t1)\n" in out
    assert "observation begins on or after (first) failure" in out
    assert "          3  observations remaining, representing" in out
