"""Explicações de erro do OpenDTA (core/hints.py)."""

import pytest


@pytest.fixture
def hrun(run):
    run.session.settings["hints"] = "on"
    return run


def lines(out):
    return out.splitlines()


def test_out_of_range(hrun):
    hrun("clear\ngen nota = 9")
    out = hrun("replace nota = 9 in 2")
    assert lines(out) == [
        "Obs. nos. out of range",
        '  → "in 2" asks for observations that do not exist; the dataset has 0 observations.',
        '    Hint: there is no data in memory: create observations with "set obs #" '
        'or load a dataset with "use".',
        "r(198);",
    ]
    hrun("set obs 3")
    out = hrun("list in 5")
    assert "valid numbers go from 1 to 3" in out


def test_variables(hrun):
    hrun("clear\nset obs 1\ngen idade = 1\ngen renda = 2\ngen rendimento = 3")
    out = hrun("list idad3")
    assert 'No variable is named "idad3"' in out and 'did you mean "idade"?' in out
    out = hrun("list rend")
    assert '"rend" is the beginning of more than one variable: renda, rendimento.' in out
    out = hrun("gen idade = 2")
    assert 'replace idade = ...' in out
    out = hrun('gen z = idade + "a"')
    assert "mixed text (string) and numbers" in out


def test_commands_and_files(hrun, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    out = hrun("sumarize x")
    assert '"sumarize" is not a command' in out
    (tmp_path / "pessoas.dta").write_bytes(b"")
    out = hrun("use pesoas")
    assert "files with similar names there: pessoas.dta" in out
    out = hrun('display "abc')
    assert "a double quote is not closed" in out or "unmatched quote" in out
    hrun("clear\nset obs 1\ngen x = 1")
    out = hrun("use pessoas")
    assert 'add the clear option' in out


def test_off_and_error_command(hrun):
    out = hrun("error 198")
    assert out == "invalid syntax\nr(198);\n"          # erro pedido de propósito
    hrun("set hints off")
    assert hrun("foo") == "command foo is unrecognized\nr(199);\n"
    out = hrun("capture noisily foo")
    assert out == "command foo is unrecognized\n"
