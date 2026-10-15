"""OpenDTA error explanations (core/hints): what went wrong, why, how to fix it and a help link,
looking at the command and at the data in memory."""

import pytest

from opendta.core.errors import StataError
from opendta.core.hints import Explanation, explain, render

SETUP = """clear
set obs 10
gen x = _n
gen y = mod(_n, 2)
gen s = "a"
label define sl 0 "masc" 1 "fem"
label values y sl
replace x = . in 1/9
gen renda = .
gen str5 cod = string(_n)
gen idade = 20 + _n
"""


@pytest.fixture
def hrun(run):
    run.session.settings["hints"] = "on"
    run(SETUP)
    return run


# (command, phrases the explanation must contain)
CASES = [
    ("list if y = 1", ['uses a single "="', "How to fix: list if y == 1", "help operators"]),
    ("keep if y = 1 & idade > 2", ["How to fix: keep if y == 1 & idade > 2"]),
    ("gen z == 1", ['single "="', "How to fix: gen z = 1"]),
    ("regress y x", ["Only 1 of the 10 observations has a value in all of y and x", "x is missing in 9",
                     'misstable summarize y x', "help missing"]),
    ("regress y x renda", ["No observation has a value in all of y, x and renda", "renda is missing in 10"]),
    ("regress y idade if y == 5", ['No observation satisfies "if y == 5"', 'count if y == 5']),
    ("regress y cod", ["cod is a string (text) variable", 'destring cod, replace']),
    ("regress y s", ['encode s, generate(s_n)']),
    ("logit x y", ["x is never 0 (always positive)", "its values: 10 (1)"]),
    ("by y: gen z = _n", ["by y: needs the data sorted by y", "How to fix: bysort y: gen z = _n", "help by"]),
    ("gen 2a = 1", ["it starts with a digit", 'use "_2a" instead']),
    ("gen açúcar = 1", ['use "acucar" instead']),
    ("gen nova var = 1", ['"nova var" has a space', "nova_var"]),
    ("summarize idade, detal", ['did you mean "detail"?']),
    ("summarize idade, foo", ["the options of summarize are:", "detail"]),
    ("xtreg y idade", ['xtset cod timevar', "help xtset"]),
    ("gen m = mean(idade)", ['"mean()" is an egen function', "How to fix: egen m = mean(idade)", "help egen"]),
    ("gen m = nchar(s)", ["in Stata use strlen()", "How to fix: gen m = strlen(s)"]),
    ("list if y == masc", ["value label of code 0 of y", '"masc":sl']),
    ("list if s == abc", ['How to fix: list if s == "abc"']),
    ("replace s = 5", ["s is a string (text) variable, but the new value is a number"]),
    ("generate t = idade + s", ["mixes numeric variables (idade) with string variables (s)"]),
    ("count if s == 1", ['count if s == "1"']),
    ('count if y == "1"', ['y is numeric, but "1" is text in quotes', "How to fix: count if y == 1"]),
    ("x = 1", ["must start with replace", "How to fix: replace x = 1"]),
    ("novavar = 1", ["must start with generate", "How to fix: generate novavar = 1"]),
    ("summary idade", ["Stata's equivalent is summarize"]),
    ("sumarize idade", ['did you mean "summarize"']),
    ("egen q = mean(idade) if", ['"if" is not followed by a condition']),
    ('display "abc', ["a double quote is not closed"]),
    ("rename x", ["rename needs two names"]),
    ("keep y x x2", ['No variable is named "x2"', 'did you mean "x"']),
    ("predict p", ["no model has been fitted yet"]),
    ("assert idade < 25", ["false in 6 observations of 10", "list if !(idade < 25)"]),
    ("help regres", ['did you mean "regress"']),
    ("gen idade = 1", ['replace idade = ...']),
    ("set obs 2", ["set obs can only add observations"]),
    ("stcox idade", ["stset"]),
    ("forvalues i = 1/3 display `i'", ["forvalues i = 1/10 {"]),
    ("}", ["closes a block that was never opened"]),
    ("mata: (1,2) + (1,2,3)", ["do not have compatible sizes", "help mata"]),
    ("mata: zz + 1", ["In Mata, zz does not exist yet"]),
    ("list in 50", ['"in 50" asks for observations that do not exist; the dataset has 10 observations.']),
    ("use nada", ["add the clear option"]),
]


@pytest.mark.parametrize("command,phrases", CASES, ids=[c for c, _ in CASES])
def test_explanations(hrun, command, phrases):
    out = hrun(command)
    flat = " ".join(out.split())
    for phrase in phrases:
        assert " ".join(phrase.split()) in flat, (phrase, out)


@pytest.mark.parametrize("command", [c for c, _ in CASES])
def test_every_case_is_explained_in_portuguese(hrun, command):
    hrun("set language pt")
    out = hrun(command)
    assert "→" in out and "Como corrigir:" in out or "Saiba mais:" in out, out
    assert "How to fix" not in out and "Learn more" not in out


def test_matrix_sizes(hrun):
    out = hrun("matrix A = (1,2\\3,4)\nmatrix B = (1,2,3)\nmatrix C = A*B")
    assert "A 2×2, B 1×3" in out


def test_merge_key_not_unique(hrun, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    hrun("save outra, replace")
    out = hrun("merge 1:1 y using outra")
    assert "repeats across observations" in out and "duplicates report y" in out and "m:1" in out


def test_files(hrun, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "pessoas.dta").write_bytes(b"")
    (tmp_path / "dados.csv").write_text("a,b\n1,2\n")
    hrun("clear")
    out = hrun("use pesoas")
    assert "files with similar names there: pessoas.dta" in out
    out = hrun("use dados")
    assert 'import delimited "dados.csv", clear' in out
    hrun("set obs 1\ngen v = 1\nsave arq")
    out = hrun("save arq")
    assert "save arq, replace" in out


def test_layout_and_styles():
    pieces = render(Explanation("What happened.", "Why it happened.", "do this", "regress"), 79)
    text = "".join(t for t, _ in pieces)
    assert text == ("  → What happened.\n    Why it happened.\n    How to fix: do this\n"
                    "    Learn more: help regress\n")
    assert ("How to fix: ", "hint_label") in pieces and ("help regress", "hint_link") in pieces
    long = render(Explanation("word " * 40), 40)
    lines = "".join(t for t, _ in long).splitlines()
    assert all(len(ln) <= 40 for ln in lines) and all(ln.startswith("    ") for ln in lines[1:])


def test_off_capture_and_error_command(hrun):
    out = hrun("error 198")
    assert out == "invalid syntax\nr(198);\n"          # error raised on purpose
    hrun("set hints off")
    assert hrun("foo") == "command foo is unrecognized\nr(199);\n"
    out = hrun("capture noisily foo")
    assert out == "command foo is unrecognized\n"


def test_unknown_errors_have_no_explanation(hrun):
    assert explain(StataError(9999, "some message nobody explains"), hrun.session) is None


def test_link_only_to_existing_help(hrun):
    expl = explain(StataError(5, "not sorted"), hrun.session)
    assert expl.help == "by"
    err = StataError(198, "invalid syntax")
    err.context = ("nosuchcommand", "nosuchcommand x")
    expl = explain(err, hrun.session)
    assert expl is None or expl.help == ""
