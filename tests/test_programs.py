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


# -- syntax ------------------------------------------------------------------------------

def _data(run):
    run("clear\nset obs 5\ngen idade = _n * 10\ngen renda = _n\nreplace renda = . in 2\n"
        'gen str3 uf = "MG"')


def test_syntax_varlist_if_in_options(run):
    _data(run)
    run("program p\n syntax [varlist] [if] [in] [, Detail BY(varname) "
        "Level(integer 95) noLOG Title(string) *]\n"
        " display \"[`varlist'][`if'][`in'][`detail'][`by'][`level'][`log'][`title'][`options']\"\nend")
    out = run('p idade ren if idade > 10 in 1/4, d by(uf) l(90) nolog title("Um título") xyz(3)')
    assert out == '[idade renda][if idade > 10][in 1/4][detail][uf][90][nolog][Um título][xyz(3)]\n'
    assert run("p") == "[idade renda uf][][][][][95][][][]\n"        # varlist padrão: todas
    run("program pn\n syntax varlist(numeric)\nend")
    run("pn uf")
    assert run.rc == 109
    run("pn")
    assert run.rc == 100
    run("p idade, by(idade renda)")
    assert run.rc == 103


def test_syntax_required_and_errors(run):
    _data(run)
    run("program q\n syntax varname [using/] [fweight] =/exp , Gen(name)\n"
        " display \"`varlist'|`using'|`weight'|`exp'|`gen'\"\nend")
    assert run("q idade = 2*3 using \"a b.dta\", gen(z)") == "idade|a b.dta||2*3|z\n"
    run("q idade = 1")
    assert run.rc == 198                                     # gen() obrigatória
    run("q idade renda = 1, gen(z)")
    assert run.rc == 103
    run("program r\n syntax [anything] [, opt]\n display `\"`anything'\"'\nend")
    assert run('r a "b c" (d)') == 'a "b c" (d)\n'
    run("r, outra")
    assert run.rc == 198


def test_syntax_newvarlist_numlist(run):
    _data(run)
    run("program nv\n syntax newvarlist(max=2) [, Values(numlist ascending)]\n"
        " display \"`varlist'|`typlist'|`values'\"\nend")
    assert run("nv a double b, v(1/3 10)") == "a b|float double|1 2 3 10\n"
    run("nv idade")
    assert run.rc == 110


def test_gettoken(run):
    run('local s `"um "dois tres" (quatro cinco) seis"\'')
    run("gettoken a s : s\ngettoken b s : s\ngettoken c s : s, match(par)")
    loc = run.session.macros.get_local
    assert (loc("a"), loc("b"), loc("c"), loc("par")) == ("um", "dois tres", "quatro cinco", "(")
    assert loc("s") == " seis"
    run('local t "x=1,y=2"\ngettoken k t : t, parse("=,")')
    assert loc("k") == "x" and loc("t") == "=1,y=2"
    run('local u `"  "q" r"\'\ngettoken v : u, quotes')
    assert loc("v") == '"q"' and loc("u") == '  "q" r'


def test_marksample(run):
    _data(run)
    run("program m\n syntax varlist [if] [in]\n marksample touse\n"
        " count if `touse'\n markout `touse' idade\n count if `touse'\nend")
    assert run("m renda if idade < 50") == "  3\n  3\n"
    assert not any(n.startswith("__") for n in run.session.data.names)   # temporária apagada


# -- temporários e preserve ----------------------------------------------------------------

def test_tempvar_tempname_tempfile(run, tmp_path):
    _data(run)
    run("program t\n tempvar a b\n tempname sc\n tempfile f\n gen `a' = idade * 2\n"
        " scalar `sc' = 5\n save `f'\n global keep \"`a' `sc' `f'\"\n display \"`a'\"\nend")
    out = run("t")
    a, sc, f = run.session.macros.get_global("keep").split(" ", 2)
    assert out.splitlines()[-1].startswith("__") and not run.session.data.has(a)
    assert sc not in run.session.scalars
    from pathlib import Path
    assert not Path(f).exists() and not Path(f + ".dta").exists()


def test_preserve_restore(run):
    _data(run)
    run("preserve\ndrop if idade > 20")
    assert run.session.data.nobs == 2
    run("restore")
    assert run.session.data.nobs == 5
    run("restore")
    assert run.rc == 622
    run("preserve\npreserve")
    assert run.rc == 621
    run("restore, not\nkeep in 1")
    assert run.session.data.nobs == 1
    _data(run)
    run("preserve\nkeep in 1/2\nrestore, preserve\nkeep in 1\nrestore")
    assert run.session.data.nobs == 5
    # dentro de programa: restaura sozinho no fim, mesmo com erro
    run("program pr\n preserve\n drop in 1/3\n error 459\nend")
    run("capture pr")
    assert run.session.data.nobs == 5


def test_preserve_in_dofile(run, tmp_path):
    _data(run)
    f = tmp_path / "p.do"
    f.write_text("preserve\nkeep in 1\ntempvar z\ngen `z' = 1\n")
    run(f'quietly do "{f}"')
    assert run.session.data.nobs == 5 and run.session.data.nvars == 3


# -- adopath ---------------------------------------------------------------------------

def test_ado_autoload_which_findfile(run, tmp_path):
    plus = tmp_path / "plus"
    (plus / "o").mkdir(parents=True)
    (plus / "o" / "oi.ado").write_text('*! version 1.0  exemplo\nprogram oi\n display "oi `0\'"\nend\n')
    run(f'sysdir set PLUS "{plus}"')
    assert run("oi mundo") == "oi mundo\n"
    out = run("which oi")
    assert str(plus / "o" / "oi.ado") in out and "*! version 1.0  exemplo" in out
    assert run("which display") == "built-in command:  display\n"
    run("which naoexiste")
    assert run.rc == 111
    run("findfile oi.ado")
    assert run.session.r["fn"] == str(plus / "o" / "oi.ado")
    run("discard")
    assert "oi" not in run.session.programs
    (tmp_path / "ruim.ado").write_text("display 1\n")
    run(f'adopath + "{tmp_path}"')
    run("ruim")
    assert run.rc == 199
    out = run("adopath")
    assert "(PLUS)" in out and str(tmp_path) in out
