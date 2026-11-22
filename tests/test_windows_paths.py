"""File names in messages keep the / that was typed, also on Windows (compat 2006)."""

from pathlib import Path, PureWindowsPath

from opendta.commands._util import shown_path


def test_shown_path_uses_forward_slashes():
    # on Windows str() of this path is "pasta\\a.txt"; Stata shows what was typed
    assert shown_path(PureWindowsPath("pasta/a.txt")) == "pasta/a.txt"
    assert shown_path(Path("pasta/a.txt")) == "pasta/a.txt"
    assert shown_path("as typed.dta") == "as typed.dta"


def test_messages_show_the_typed_path(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run("clear\nset obs 1\ngen x = 1")
    for cmd, msg in (('save "nada_xyz/a.dta", replace', "file nada_xyz/a.dta could not be opened"),
                     ('export delimited "nada_xyz/a.csv", replace',
                      "file nada_xyz/a.csv could not be opened"),
                     ('use "nada_xyz/a.dta", clear', "file nada_xyz/a.dta not found")):
        assert msg in run(f"capture noisily {cmd}")
