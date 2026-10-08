from opendta.lang.lexer import split_commands


def texts(src):
    return [ln.text for ln in split_commands(src) if ln.kind != "comment"]


def test_line_comments():
    assert texts("* comentário\ndisplay 1\n") == ["display 1"]
    assert texts("display 1 // fim\n") == ["display 1"]
    assert texts('display "a//b"\n') == ['display "a//b"']


def test_block_comments():
    assert texts("display /* x */ 1\n") == ["display   1"]
    assert texts("di 1 /* várias\nlinhas */\ndi 2\n") == ["di 1", "di 2"]
    assert texts("/* /* aninhado */ ainda */ di 3\n") == ["di 3"]


def test_continuation():
    assert texts("display 1 + ///\n 2\n") == ["display 1 +   2"]


def test_delimit():
    src = "#delimit ;\ndisplay 1\n + 2;\ndisplay 3;\n#delimit cr\ndisplay 4\n"
    assert texts(src) == ["#delimit ;", "display 1  + 2", "display 3", "#delimit cr", "display 4"]


def test_line_numbers():
    lines = split_commands("\n\ndisplay 1\n* c\ndisplay 2\n")
    assert [(ln.lineno, ln.kind) for ln in lines] == [(3, "cmd"), (4, "comment"), (5, "cmd")]


def test_raw_lines_for_echo():
    src = '#delimit ;\ndisplay "a"\n    "b";\n#delimit cr\n    di 1 // x\n'
    lines = split_commands(src)
    assert lines[1].raw == ('display "a"', '    "b";')
    assert lines[3].raw == ("    di 1 // x",)
