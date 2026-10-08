from pathlib import Path


def test_display(run):
    assert run("display 2+2") == "4\n"
    assert run('di "a" "b"') == "ab\n"
    assert run("di 1/3") == ".33333333\n"
    assert run("di %9.2f _pi") == "     3.14\n"
    assert run('di "x" _col(5) "y"') == "x   y\n"
    assert run('di _dup(3) "-"') == "---\n"
    assert run('di "a" _continue') == "a"


def test_locals_and_globals(run):
    assert run('local x = 2+3\ndi `x\'') == "5\n"
    assert run('local s "hello world"\ndi "`s\'"') == "hello world\n"
    assert run('global g 7\ndi $g + ${g}') == "14\n"
    assert run('local i 2\nlocal v2 ok\ndi "`v`i\'\'"') == "ok\n"
    assert run('local n 1\ndi `++n\'\ndi `n\'') == "2\n2\n"
    assert run("di \"`undefined'\"") == "\n"
    assert run('di "`=2*21\'"') == "42\n"


def test_extended_functions(run):
    assert run('local n : word count a "b c" d\ndi `n\'') == "3\n"
    assert run('local w : word 2 of a b c\ndi "`w\'"') == "b\n"
    assert run('local f : display %5.2f 1/3\ndi "`f\'"') == " 0.33\n"
    assert run('local a x y\nlocal b y z\nlocal u : list a | b\ndi "`u\'"') == "x y z\n"


def test_loops(run):
    assert run("forvalues i = 1/3 {\n di `i'\n}") == "1\n2\n3\n"
    assert run("forvalues i = 10(-5)0 {\n di `i'\n}") == "10\n5\n0\n"
    assert run("forvalues i = 1 3 to 7 {\n di `i'\n}") == "1\n3\n5\n7\n"
    assert run('foreach x in a "b c" {\n di "`x\'"\n}') == "a\nb c\n"
    assert run("foreach n of numlist 1/2 5 {\n di `n'\n}") == "1\n2\n5\n"
    assert run("local i 0\nwhile `i' < 3 {\n local ++i\n di `i'\n}") == "1\n2\n3\n"


def test_continue_break(run):
    src = "forvalues i = 1/5 {\n if `i' == 2 continue\n if `i' == 4 continue, break\n di `i'\n}"
    assert run(src) == "1\n3\n"


def test_if_else(run):
    src = "local x 5\nif `x' > 3 {\n di \"maior\"\n}\nelse {\n di \"menor\"\n}"
    assert run(src) == "maior\n"
    src = "local x 1\nif `x' > 3 {\n di 1\n}\nelse if `x' > 0 {\n di 2\n}\nelse {\n di 3\n}"
    assert run(src) == "2\n"
    assert run("if 1 di \"sim\"") == "sim\n"
    assert run("if 0 di \"sim\"") == ""


def test_prefixes(run):
    assert run('quietly di "x"') == ""
    assert run('quietly {\n di "a"\n noisily di "b"\n}') == "b\n"
    assert run("capture error 111\ndi _rc") == "111\n"
    assert run("capture di 1\ndi _rc") == "0\n"
    assert run("capture noisily error 198").endswith("r(198);\n")


def test_errors(run):
    out = run("foo")
    assert out == "command foo is unrecognized\nr(199);\n"
    assert run.rc == 199
    run("di nada")
    assert run.rc == 111
    run('di "a" + 1')
    assert run.rc == 109


def test_abbreviations(run):
    assert run("di 1") == "1\n"
    assert run("disp 1") == "1\n"
    out = run("d")
    assert run.rc == 0 and "Contains data" in out   # 'd' é describe


def test_scalar(run):
    assert run("scalar a = 2\nsca b = a * 3\ndi b") == "6\n"
    assert run('scalar s = "txt"\ndi s') == "txt\n"


def test_do_file(run, tmp_path: Path):
    f = tmp_path / "t.do"
    f.write_text('args a b\ndi "`a\'-`b\'"\nlocal z 1\n')
    out = run(f'do "{f}" um dois')
    assert "um-dois" in out
    assert out.rstrip().endswith("end of do-file")
    # locais do do-file não vazam
    assert run("di \"`z'\"") == "\n"


def test_do_file_error_stops(run, tmp_path: Path):
    f = tmp_path / "e.do"
    f.write_text("di 1\nfoo\ndi 2\n")
    out = run(f'do "{f}"')
    assert "command foo is unrecognized" in out
    assert "\n2\n" not in out
    assert run.rc == 199


def test_tokenize(run):
    assert run('tokenize "a b c"\ndi "`2\'"') == "b\n"
    assert run('tokenize "x,y", parse(",")\ndi "`3\'"') == "y\n"


def test_confirm(run):
    run("confirm number 3")
    assert run.rc == 0
    run("confirm integer number 3.5")
    assert run.rc == 7
