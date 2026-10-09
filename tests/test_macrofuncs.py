"""Funções estendidas de macro (local x : ...)."""


def test_variable_properties(run):
    run('clear\nset obs 2\ngen byte idade = 1\nlabel variable idade "Idade (anos)"\n'
        'label define sx 1 "Masc" 2 "Fem"\nlabel values idade sx\nformat idade %5.1f\n'
        'label data "Base"\nchar idade[fonte] "censo"\nsort idade')
    loc = run.session.macros.get_local
    run("local a : type idade\nlocal b : format idade\nlocal c : value label idade\n"
        "local d : variable label idade\nlocal e : data label\nlocal f : sortedby\n"
        "local g : label sx 2\nlocal h : label (idade) 1\nlocal i : label sx 9\n"
        "local j : label sx 9, strict\nlocal k : char idade[fonte]\nlocal l : char idade[]\n"
        "local m : label sx maxlength")
    assert [loc(x) for x in "abcdefghijklm"] == [
        "byte", "%5.1f", "sx", "Idade (anos)", "Base", "idade", "Fem", "Masc", "9", "",
        "censo", "fonte", "4"]
    run("local z : type nada")
    assert run.rc == 111


def test_strings_and_dirs(run, tmp_path):
    loc = run.session.macros.get_local
    run('local s "a b a c a"\nlocal t : subinstr local s "a" "x", all count(local n)')
    assert loc("t") == "x b x c x" and loc("n") == "3"
    run('local t : subinstr local s "a" "", word count(local n)')
    assert loc("t") == " b a c a"     # o Stata mantém o espaço and loc("n") == "1"
    run("local u : copy local s\nlocal v : strlen local s")
    assert loc("u") == "a b a c a" and loc("v") == "9"
    (tmp_path / "um.do").write_text("")
    (tmp_path / "dois.dta").write_text("")
    (tmp_path / "sub").mkdir()
    run(f'local f : dir "{tmp_path}" files "*.do"\nlocal g : dir "{tmp_path}" dirs "*"')
    assert loc("f") == '"um.do"' and loc("g") == '"sub"'
    run("clear\nset obs 1\ngen x = 1\nlocal p : permname x")
    assert loc("p") == "x1"
    run("local q : foo bar")
    assert run.rc == 198


def test_results_lists(run):
    run('program r1, rclass\n return scalar n = 1\n return local cmd "r1"\nend\nr1')
    run("local a : r(scalars)\nlocal b : r(macros)")
    assert run.session.macros.get_local("a") == "n" and run.session.macros.get_local("b") == "cmd"
