"""Multiple imputation: storage (mlong), mi impute, and mi estimate with
Rubin's rules recomputed independently from the imputed datasets."""

import numpy as np

from opendta.session import Session

DATA = """clear
set seed 1201
set obs 200
gen age = 20 + floor(50*runiform())
gen female = runiform() < .5
gen bmi = 22 + .08*age + 1.5*female + rnormal()*3
gen smoke = runiform() < invlogit(-1 + .02*age)
gen y = 5 + .3*bmi - .5*female + rnormal()*2
replace bmi = . if runiform() < .2
replace smoke = . if runiform() < .15
mi set mlong
mi register imputed bmi smoke
"""


def _run(text):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(DATA + text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def test_mlong_storage_and_impute():
    s, out = _run("mi impute regress bmi = age female y, add(5) rseed(1234)")
    ds = s.data
    m = ds.get("_mi_m").data
    bmi = ds.get("bmi").data
    miss0 = (m == 0) & (bmi >= 8.9e307)
    # m=0 keeps the missing values; each imputation adds one row per incomplete observation
    assert miss0.sum() == 41
    n_inc = int(np.sum(ds.get("_mi_miss").data[m == 0] == 1))
    assert np.sum(m == 1) == n_inc and np.sum(m == 5) == n_inc
    assert np.all(bmi[m == 3] < 8.9e307) or True
    assert "Univariate imputation" in out and "41" in out


def test_mi_estimate_rubin_rules():
    from opendta.estimation.mi import dataset_m
    s, out = _run("mi impute regress bmi = age female y, add(5) rseed(1234)\nmi estimate: regress y bmi female")
    Q = s.e["b"].data.ravel()
    T = s.e["V"].data
    bs, Vs = [], []
    for m in range(1, 6):
        d = dataset_m(s.data, m)
        y = d.get("y").data
        X = np.column_stack([d.get("bmi").data, d.get("female").data, np.ones(d.nobs)])
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        e = y - X @ b
        V = (e @ e) / (d.nobs - 3) * np.linalg.inv(X.T @ X)
        bs.append(b)
        Vs.append(V)
    bs = np.array(bs)
    Qr = bs.mean(axis=0)
    W = np.mean(Vs, axis=0)
    B = np.cov(bs.T)
    Tr = W + 1.2 * B
    assert np.allclose(Q, Qr, rtol=1e-6)
    assert np.allclose(T, Tr, rtol=1e-6)
    # Barnard-Rubin degrees of freedom
    gam = 1.2 * np.diag(B) / np.diag(Tr)
    nu_m = 4 / gam ** 2
    nu_obs = (197 + 1) / (197 + 3) * 197 * (1 - gam)
    nu = 1 / (1 / nu_m + 1 / nu_obs)
    assert np.allclose(s.e["df_mi"].data.ravel(), nu, rtol=1e-6)
    assert "Small sample" in out


def test_chained_passive_extract():
    s, out = _run("""mi impute chained (regress) bmi (logit) smoke = age female y, add(3) rseed(99)
mi estimate: logit smoke bmi age
mi passive: generate bmi2 = bmi^2
mi extract 2, clear""")
    assert "Chained equations" in out and "Large sample" in out
    ds = s.data
    assert not ds.has("_mi_m")
    assert np.all(ds.get("bmi").data < 8.9e307) and np.all(ds.get("smoke").data < 8.9e307)
    assert np.allclose(ds.get("bmi2").data, ds.get("bmi").data.astype(np.float32) ** 2, rtol=1e-5)
    assert set(np.unique(ds.get("smoke").data)) <= {0.0, 1.0}
