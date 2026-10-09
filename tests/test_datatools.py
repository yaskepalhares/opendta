"""duplicates, isid, levelsof, encode/decode, destring/tostring, split,
mvencode/mvdecode e recode."""

import numpy as np
import pytest

MISS = 8.98846567431158e307


def col(run, name):
    return list(np.array(run.session.data.get(name).data))


def raw(run, name):
    return list(run.session.data.get(name).raw)


@pytest.fixture
def people(run):
    run("""clear
input str8 nome x str6 v
"ana" 1 "10"
"bia" 2 "2,5"
"ana" 1 "$30"
"caio" 3 "abc"
"ana" 1 "7"
end
""")
    return run


def test_duplicates_report(people):
    out = people("duplicates report nome x")
    assert out == (
        "\nDuplicates in terms of nome x\n\n"
        "--------------------------------------\n"
        "   copies | observations       surplus\n"
        "----------+---------------------------\n"
        "        1 |            2             0\n"
        "        3 |            3             2\n"
        "--------------------------------------\n")
    assert people.session.r["unique_value"] == 3
    out = people("duplicates report")
    assert "Duplicates in terms of all variables" in out


def test_duplicates_tag_list_drop(people):
    run = people
    run("duplicates tag nome, gen(d)")
    assert col(run, "d") == [2, 0, 2, 0, 2]
    out = run("duplicates list nome x")
    # um só grupo de duplicatas: sem a coluna group: (Stata 14)
    assert "| obs:   nome   x |" in out
    assert out.count("ana") == 3
    run("duplicates drop nome")
    assert run.rc == 198
    out = run("duplicates drop nome, force")
    assert "(2 observations deleted)" in out
    assert raw(run, "nome") == ["ana", "bia", "caio"]
    out = run("duplicates drop")
    assert "(0 observations deleted)" in out


def test_isid(people):
    run = people
    run("isid nome")
    assert run.rc == 459
    out = run("isid nome")
    assert "variable nome does not uniquely identify the observations" in out
    run("isid nome v")
    assert run.rc == 0
    run("isid v nome, sort")
    assert raw(run, "v")[0] == "$30"
    run("replace x = . in 2\nisid x v")
    assert run.rc == 459
    run("isid x v, missok")
    assert run.rc == 0


def test_levelsof(people):
    run = people
    out = run("levelsof nome")
    assert out == "`\"ana\"' `\"bia\"' `\"caio\"'\n"
    out = run("levelsof nome, clean separate(,)")
    assert out == "ana,bia,caio\n"
    run("replace x = 1.5 in 2\nreplace x = . in 4\nlevelsof x, local(xs)")
    assert run.session.macros.get_local("xs") == "1 1.5"
    assert run.session.r["levels"] == "1 1.5"
    out = run("levelsof x, missing")
    assert out == "1 1.5 .\n"


def test_encode_decode(people):
    run = people
    run("encode nome, gen(cod)")
    ds = run.session.data
    assert col(run, "cod") == [1, 2, 1, 3, 1]
    assert ds.value_labels["cod"] == {1: "ana", 2: "bia", 3: "caio"}
    assert ds.get("cod").vtype == "long"
    run("label define nm 1 \"caio\"\nencode nome, gen(c2) label(nm)")
    assert col(run, "c2") == [2, 3, 2, 1, 2]
    run("decode cod, gen(back)")
    assert raw(run, "back") == raw(run, "nome")
    run("decode x, gen(z)")
    assert run.rc == 182
    run("encode x, gen(z)")
    assert run.rc == 108


def test_destring(people):
    run = people
    out = run("destring v, gen(n1)")
    assert out == "v contains nonnumeric characters; no generate\n"
    run("replace v = \"40\" in 4")
    out = run("destring v, gen(n1) ignore(\"$,\")")
    assert out == "v: characters $ , removed; generated as byte\n"
    assert col(run, "n1") == [10, 25, 30, 40, 7]
    out = run("destring v, replace force")
    assert "v contains nonnumeric characters; v replaced as byte" in out
    assert "(2 missing values generated)" in out
    assert run.session.data.names.index("v") == 2
    run("destring v, replace dpcomma")
    assert run.rc == 0


def test_destring_dpcomma_percent(run):
    run('clear\ninput str6 a\n"1,5"\n"20%"\nend\ndestring a, gen(b) dpcomma percent')
    assert col(run, "b") == pytest.approx([1.5, 0.2])


def test_tostring(people):
    run = people
    out = run("tostring x, replace")
    assert out == "x was float now str1\n"
    assert raw(run, "x") == ["1", "2", "1", "3", "1"]
    run('gen y = 1/3\ntostring y, gen(ys)')
    assert "ys" not in run.session.data.names
    run('tostring y, gen(ys) format(%5.2f)')
    assert raw(run, "ys")[0] == "0.33"
    run('tostring y, gen(yf) force')
    assert raw(run, "yf")[0] == ".3333333433"     # float exibido com %12.0g


def test_split(run):
    run('clear\ninput str20 s\n"a b  c"\n"x,y"\n""\nend')
    out = run("split s")
    assert out == "variables created as string: \ns1  s2  s3\n"
    assert raw(run, "s3") == ["c", "", ""]
    run("split s, parse(,) gen(p) limit(1)")
    assert raw(run, "p1") == ["a b  c", "x", ""]
    assert run.session.r["nvars"] == 1
    run('clear\ninput str10 n\n"1-2"\n"3-40"\nend\nsplit n, parse(-) destring')
    assert col(run, "n2") == [2, 40]


def test_mvdecode_mvencode(run):
    run("clear\nset obs 4\ngen a = _n\nreplace a = -9 in 2\nreplace a = -8 in 3")
    out = run("mvdecode a, mv(-9 -8)")
    assert out == "           a: 2 missing values generated\n"
    run("mvencode a, mv(-1)")
    assert col(run, "a") == [1, -1, -1, 4]
    run("replace a = .a in 1\nmvdecode a, mv(-1=.b)")
    out = run("list a, clean noobs")
    assert ".a" in out and ".b" in out
    run("mvencode a, mv(4)")
    assert run.rc == 9
    run("mvencode a, mv(.a=0 \\ .b=-1) override")
    assert col(run, "a") == [0, -1, -1, 4]


def test_recode(run):
    run("clear\nset obs 6\ngen v = _n\nreplace v = . in 6")
    out = run("recode v (1 2 = 1) (3/4 = 2) (missing = 9) (else = 0)")
    assert out == "(v: 5 changes made)\n"     # 1 → 1 não conta
    assert col(run, "v") == [1, 1, 2, 2, 0, 9]


def test_recode_counts_only_changes(run):
    run("clear\nset obs 5\ngen v = _n")
    out = run("recode v (1 2 = 1) (min/3 = 7)")
    assert out == "(v: 2 changes made)\n"
    assert col(run, "v") == [1, 1, 7, 4, 5]


def test_recode_generate_labels_and_if(run):
    run("clear\nset obs 5\ngen v = _n")
    out = run('recode v (1/2 = 1 "low") (3/max = 2 "high") if v != 4, gen(g)')
    assert out == "(3 differences between v and g)\n"     # só na amostra (VERIFICAR)
    assert col(run, "g")[:3] == [1, 1, 2] and col(run, "g")[3] >= MISS
    ds = run.session.data
    assert ds.value_labels["g"] == {1: "low", 2: "high"}
    assert ds.get("g").value_label == "g"
    run("recode v (1 = 0), prefix(r_)")
    assert col(run, "r_v")[0] == 0 and col(run, "r_v")[1] == 2
    run("recode v (1 = 2.5)")
    assert run.session.data.get("v").vtype == "float"
    run("recode v 1=2")
    assert run.rc == 198
