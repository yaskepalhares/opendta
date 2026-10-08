"""log, cmdlog, assert, timer, creturn e set trace."""

import re


def test_log_text_and_smcl(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    out = run("log using t, text")
    assert "      name:  <unnamed>" in out and "  log type:  text" in out
    run('display "dentro"')
    run("log close")
    text = (tmp_path / "t.log").read_text()
    assert "opened on:" in text and "dentro\n" in text and "closed on:" in text
    run("log using t, text")
    assert run.rc == 602
    run("log using s")
    run('display "a{b}"')
    run("log close")
    smcl = (tmp_path / "s.smcl").read_text()
    assert smcl.startswith("{smcl}") and "a{c -(}b{c )-}" in smcl
    run("log close")
    assert run.rc == 606
    run("log using a, text name(um)\nlog using b, text name(dois)\nlog off um\n"
        'display "MARCADOR"\nlog close _all')
    assert "MARCADOR" not in (tmp_path / "a.log").read_text()
    assert "MARCADOR\n" in (tmp_path / "b.log").read_text()


def test_cmdlog(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run.session.run_command("cmdlog using c")
    run.session.run_command('display "oi"')
    run.session.run_command("cmdlog close")
    assert (tmp_path / "c.txt").read_text() == 'display "oi"\ncmdlog close\n'


def test_assert(run):
    run("clear\nset obs 5\ngen x = _n")
    assert run("assert x > 0") == ""
    out = run("assert x < 3")
    assert run.rc == 9 and out.startswith("3 contradictions in 5 observations\nassertion is false")
    out = run("assert x < 3, rc0")
    assert run.rc == 0 and "assertion is false" in out
    assert run("assert x == 1 in 1") == ""


def test_timer_creturn_trace(run):
    run("timer clear\ntimer on 1\ntimer off 1\ntimer list 1")
    assert run.session.r["nt1"] == 1.0
    out = run("creturn list")
    assert "c(pi) = 3.141592653589793" in out and "c(stata_version) = 14" in out
    run('program tr\n local a 5\n display `a\'\nend\nset trace on')
    out = run("tr")
    assert re.search(r"-+ begin tr ---", out) and "- display `a'" in out and "= display 5" in out
    run("set trace off")
    assert run("tr") == "5\n"


def test_help_text_and_smcl(run, tmp_path):
    from opendta.core.smcl import render_html, render_text
    out = run("help describe")
    assert "describe --" in out and "Describe the dataset" in out
    out = run("help gen")
    assert "generate --" in out
    run("help nada")
    assert run.rc == 111
    src = "{smcl}\n{title:T}\n{pstd}\nUse {cmd:x} and {help list}.{p_end}\n{synopt:{opt s:hort}}desc{p_end}\n"
    assert "Use x and list." in render_text(src) and "short" in render_text(src)
    html = render_html(src)
    assert '<a href="help:list"' in html and "<b>x</b>" in html and "<u>s</u>" in html
    (tmp_path / "a.smcl").write_text("{smcl}\n{res}resultado{txt}\n")
    assert run(f'view "{tmp_path / "a.smcl"}"') == "resultado\n"
    out = run("search matrix")
    assert "[matrix]" in out
