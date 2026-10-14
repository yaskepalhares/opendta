"""Expoente sobrescrito na notação científica exibida (set superscript)."""

from opendta.core.formats import format_value, plain_numbers, set_superscript


def test_format_value_superscript():
    set_superscript(True)
    try:
        assert format_value(1e10, "%12.0g").strip() == "1.0000×10¹⁰"
        assert format_value(2.5e-6, "%9.0g").strip() == "2.5×10⁻⁶"
        assert format_value(1234.5, "%10.3e").strip() == "1.235×10³"
        assert format_value(-1e-10, "%10.0g").strip() == "-1.0×10⁻¹⁰"
        # números em notação fixa não mudam; a largura do formato é mantida
        assert format_value(123456789, "%12.0g") == "   123456789"
        assert len(format_value(1e10, "%12.0g")) == 12
        # texto que vira dado (string(), macros) fica na forma do Stata
        assert format_value(1e10, "%12.0g", dp=False).strip() == "1.00000e+10"
        with plain_numbers():
            assert format_value(1e10, "%12.0g").strip() == "1.00000e+10"
        assert format_value(1e10, "%12.0g").strip() == "1.0000×10¹⁰"
    finally:
        set_superscript(False)


def test_set_superscript_command(run):
    run("set superscript on")
    assert run("display 1e10/3").strip() == "3.333×10⁹"
    assert run('local a : display 1e10/3\ndisplay "`a\'"').strip() == "3.333e+09"
    assert run("display string(1e10/3)").strip() == "3.33e+09"
    assert "1.0000×10¹⁰" in run("mata: 1e10")
    assert run("mata: strofreal(1e10)").strip() == "1.00e+10"
    assert run('mata: sprintf("%g", 1e20)').strip() == "1.00e+20"
    run("set superscript off")
    assert run("display 1e10/3").strip() == "3.333e+09"
    assert "invalid syntax" in run("set superscript maybe") and run.rc == 198


def test_file_write_keeps_stata_form(run, tmp_path):
    f = tmp_path / "out.txt"
    run("set superscript on")
    run(f'file open h using "{f}", write text replace\nfile write h (1e10/3) _n\nfile close h')
    assert f.read_text().strip() == "3.333e+09"


def test_new_session_default_on():
    from opendta.session import Session
    from opendta.core.formats import superscript
    s = Session()
    assert s.settings["superscript"] == "on" and superscript()
    set_superscript(False)
