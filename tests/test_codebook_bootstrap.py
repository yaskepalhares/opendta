"""Fase 8: codebook, misstable e bootstrap."""

import numpy as np
import pytest

DATA = """clear
set seed 803
set obs 80
gen x = rnormal(10, 2)
gen y = 1 + 2*x + rnormal()
gen g = mod(_n, 4)
label define gl 0 "zero" 1 "um" 2 "dois" 3 "tres"
label values g gl
replace x = . in 1/5
gen str8 s = "a b" if mod(_n, 2)
"""


def test_codebook_numeric_and_labels(run):
    run(DATA)
    out = run("codebook x g s")
    assert "type:  numeric (float)" in out
    assert "missing .:  5/80" in out
    assert "Freq.   Numeric  Label" in out
    assert "zero" in out
    assert 'missing "":  40/80' in out
    assert "variable has embedded blanks" in out


def test_codebook_compact(run):
    run(DATA)
    out = run("codebook x y, compact")
    lines = [ln for ln in out.splitlines() if ln.startswith("x ")]
    assert lines and "75" in lines[0]


def test_misstable_summarize_and_patterns(run):
    run(DATA)
    out = run("misstable summarize")
    assert "x |" in out and "5" in out
    assert "y |" not in out
    out = run("misstable patterns, frequency")
    assert "Missing-value patterns" in out
    assert "(1) x" in out


def test_bootstrap_eclass_matches_resampling_spread(run):
    run(DATA)
    out = run("bootstrap, reps(200) seed(1) nodots: regress y x")
    assert "Observed   Bootstrap" in out
    e = run.session.e
    assert e["N_reps"] == 200
    assert e["vce"] == "bootstrap" and e["vcetype"] == "Bootstrap"
    se_boot = np.sqrt(np.diag(e["V"].data))
    run("quietly regress y x")
    se_ols = np.sqrt(np.diag(run.session.e["V"].data))
    # with homoskedastic errors, the bootstrap stays close to the usual standard error
    assert se_boot[0] == pytest.approx(se_ols[0], rel=0.3)


def test_bootstrap_rclass_and_seed_reproducible(run):
    run(DATA)
    run("bootstrap m=r(mean), reps(50) seed(7) nodots: summarize y")
    v1 = run.session.e["V"].data[0, 0]
    assert run.session.e["cmd"] == "bootstrap"
    run("bootstrap m=r(mean), reps(50) seed(7) nodots: summarize y")
    assert run.session.e["V"].data[0, 0] == v1


def test_vce_bootstrap_option_and_estat(run):
    run(DATA)
    out = run("regress y x, vce(bootstrap, reps(30) seed(2) nodots)")
    assert "Bootstrap" in out
    out = run("estat bootstrap, all")
    assert "(BC)" in out and "(P)" in out and "(N)" in out
    run("predict yh")
    assert run.session.data.has("yh")


def test_bootstrap_requires_explist_for_rclass(run):
    run(DATA)
    run("capture bootstrap, reps(5): summarize y")
    assert run.session.rc == 198
