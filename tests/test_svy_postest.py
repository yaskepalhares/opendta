"""Postestimation after svy: (compat 2004): equation names, adjusted Wald tests, margins and estat."""

import pytest

DATA = """clear
set seed 123
set obs 600
gen psu = ceil(_n/20)
gen strat = mod(psu,3)+1
gen w = 1+round(runiform()*3,.1)
gen ano = 2006 + floor(runiform()*10)
gen sexo = runiform()<.5
gen idade = 18+floor(runiform()*60)
gen faixa = 1+floor(runiform()*3)
gen pc = 1+floor(runiform()*10)+sexo
gen yb = runiform() < .3+.1*sexo
svyset psu [pw=w], strata(strat)
"""


@pytest.fixture
def svyrun(run):
    run(DATA)
    return run


def test_equation_prefix_after_svy_poisson(svyrun):
    svyrun("quietly svy: poisson pc ano i.sexo ib2.faixa c.idade")
    for cmd, line in [("lincom ano", "[pc]ano = 0"), ("lincom [pc]ano", "[pc]ano = 0"),
                      ("lincom [pc]_b[ano]", "[pc]ano = 0"), ("lincom [pc]c.idade", "[pc]idade = 0"),
                      ("lincom ano + [pc]c.idade", "[pc]ano + [pc]idade = 0"),
                      ("test [pc]ano", "[pc]ano = 0"), ("test [pc]: ano", "[pc]ano = 0"),
                      ("test 3.faixa = 1.faixa", "- [pc]1.faixa + [pc]3.faixa = 0")]:
        out = svyrun(cmd)
        assert f" ( 1)  {line}\n" in out, cmd
    assert "-.0042621" in svyrun("lincom [pc]ano")
    assert "_nl_1:  [pc]_b[ano]" in svyrun("nlcom [pc]_b[ano]")
    assert "-.0042621" in svyrun("nlcom [pc]_b[ano]")
    for e in ("_b[pc:ano]", "[pc]_b[ano]", "[pc]ano"):
        assert svyrun(f"display {e}").strip() == "-.0042621", e
    assert svyrun("display _b[2b.faixa]").strip() == "0"
    assert svyrun("display _b[c.idade]").strip() == ".00090785"
    assert svyrun("lincom ib2.faixa").startswith("ib2:  operator invalid")


def test_equation_prefix_after_svy_logit_and_regress(svyrun):
    svyrun("quietly svy: logit yb ano i.sexo i.faixa, or")
    assert " ( 1)  [yb]ano = 0\n" in svyrun("lincom [yb]ano, or")
    assert " ( 1)  [yb]2.faixa - [yb]3.faixa = 0\n" in svyrun("lincom 2.faixa - 3.faixa, or")
    svyrun("quietly svy: poisson pc ano i.sexo##c.idade")
    float(svyrun("display _b[1.sexo#c.idade]"))   # accepted (no Stata value in the log)
    svyrun("quietly svy: regress pc ano i.sexo")
    assert " ( 1)  ano = 0\n" in svyrun("lincom ano")
    assert svyrun("lincom [pc]ano").startswith("equation [pc] not found")


def test_adjusted_wald_after_svy(svyrun):
    svyrun("quietly svy: poisson pc ano i.sexo ib2.faixa c.idade")
    out = svyrun("testparm i.faixa")
    assert out.startswith("\nAdjusted Wald test\n\n ( 1)  [pc]1.faixa = 0\n")
    assert "F(  2,    26) =    0.10" in out and "0.9049" in out
    out = svyrun("test ano")
    assert "Adjusted Wald test" in out and "F(  1,    27) =    0.33" in out and "0.5707" in out
    out = svyrun("testparm i.faixa, nosvyadjust")
    assert "Adjusted" not in out and "F(  2,    27) =    0.10" in out and "0.9015" in out
    out = svyrun("test 1.sexo, accumulate")
    assert " ( 2)  [pc]1.sexo = 0" in out and "F(  2,    26)" in out
    svyrun("quietly svy: regress pc ano i.sexo")
    assert "Adjusted Wald test" in svyrun("test ano")
    svyrun("quietly regress pc ano i.sexo")
    assert "Adjusted" not in svyrun("test ano")