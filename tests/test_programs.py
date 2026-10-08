"""program define, chamada, return/ereturn/sreturn e escopos."""


def test_define_and_call(run):
    run("program define hello\n display `\"oi `1' [`0']\"'\nend")
    assert run('hello Ana "Bia Lima"') == 'oi Ana [Ana "Bia Lima"]\n'
    assert "hello" in run.session.programs
    out = run("program define hello\nend")
    assert run.rc == 110
    run("program drop hello")
    run("hello")
    assert run.rc == 199
    run("capture program drop naoexiste")
    assert run.session.rc == 111


def test_locals_are_scoped(run):
    run('local x "fora"\nprogram p\n local x "dentro"\n display "`x\'"\nend')
    assert run("p") == "dentro\n"
    assert run('display "`x\'"') == "fora\n"


def test_rclass(run):
    run("program media, rclass\n args a b\n return scalar m = (`a' + `b') / 2\n"
        " return local txt \"de `a' e `b'\"\nend")
    run("media 3 4")
    assert run.session.r["m"] == 3.5 and run.session.r["txt"] == "de 3 e 4"
    assert run("display \"`r(txt)'\" r(m)") == "de 3 e 43.5\n"
    run("program falha, rclass\n return scalar m = 99\n error 459\nend\ncapture falha")
    assert run.session.r["m"] == 3.5                  # erro: r() não muda
    run("return scalar z = 1")
    assert run.rc == 151


def test_exit_and_nesting(run):
    run('program sai\n display "a"\n exit\n display "b"\nend')
    assert run("sai") == "a\n"
    run('program erro\n exit 459\nend')
    run("capture erro")
    assert run.session.rc == 459
    run("program interno, rclass\n return scalar v = 7\nend\n"
        "program externo, rclass\n interno\n return scalar w = r(v) * 2\nend\nexterno")
    assert run.session.r == {"w": 14.0}


def test_eclass_sclass(run):
    run('program e1, eclass\n ereturn clear\n ereturn scalar N = 10\n ereturn local cmd "e1"\nend\ne1')
    assert run.session.e == {"N": 10.0, "cmd": "e1"}
    run('program s1, sclass\n sreturn local k "abc"\nend\ns1')
    assert run('display "`s(k)\'"') == "abc\n"
    out = run("ereturn list")
    assert "e(N) =  10" in out and 'e(cmd) : "e1"' in out


def test_program_in_dofile_echo(run, tmp_path):
    f = tmp_path / "d.do"
    f.write_text('program define hi\n    display "hi"\nend\nhi\n')
    out = run(f'do "{f}"')
    assert '. program define hi\n  1.     display "hi"\n  2. end\n' in out
    assert "\n. hi\nhi\n" in out
