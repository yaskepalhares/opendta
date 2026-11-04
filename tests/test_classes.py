"""Stata's class system ([P] class, [P] class exit, [P] classutil), from the manual's examples."""

import pytest

COORD = """version 14
class coordinate {
    double x
    double y
}
program .new
    if "`0'" != "" {
        .set `0'
    }
end
program .set
    args x y
    .x = `x'
    .y = `y'
end
program .macroexpand
    local tosub : display "(" `.x' "," `.y' ")"
    class exit "`tosub'"
end
"""

LINE = """version 14
class line {
    coordinate c0
    coordinate c1
}
program .set
    args x0 y0 x1 y1
    .c0.set `x0' `y0'
    .c1.set `x1' `y1'
end
program .length
    class exit sqrt((`.c0.y'-`.c1.y')^2 + (`.c0.x'-`.c1.x')^2)
end
program .midpoint
    local cx = (`.c0.x' + `.c1.x')/2
    local cy = (`.c0.y' + `.c1.y')/2
    tempname b
    .`b' = .coordinate.new
    .`b'.set `cx' `cy'
    class exit .`b'
end
"""

COORD2 = """*! version 1.0.1
class coordinate2 {
classwide:
    double x_origin = 0
instancespecific:
    double x = 0
}
"""

TODO = """class todolist {
    double n = 0
    string name = "nobody"
    array list = {"show second syntax", "mark as done"}
    double b = (.b)
    string s3 = ("no" + "body")
    mixed = {1, , "three"}
}
"""

VARC = """class varc {
classwide:
    n = 0
instancespecific:
    string x
}
program .new
    .n = `.n' + 1
    .x = "vc_`.n'"
    generate `.x' = 1
end
program .oncopy
    .n = `.n' + 1
    .x = "vc_`.n'"
    generate `.x' = `.`.oncopy_src'.x' + 1
end
program .destructor
    capture drop `.x'
end
"""

BASE = """class base {
    double size = 1
}
program .new
    .size = 10
end
program .zifl
    display "base `.size'"
end
program .who
    class exit "base"
end
"""

ALT = """class alt {
    double extra = 5
}, inherit(base)
program .zifl
    .Super.zifl
    display "alt `.extra'"
end
"""


@pytest.fixture
def cls(run, tmp_path, monkeypatch):
    for name, text in (("coordinate", COORD), ("line", LINE), ("coordinate2", COORD2),
                       ("todolist", TODO), ("varc", VARC), ("base", BASE), ("alt", ALT)):
        (tmp_path / f"{name}.class").write_text(text)
    monkeypatch.chdir(tmp_path)          # "." is on the adopath
    return run


def show(run, ref):
    return run(f'display `"{ref}"\'').rstrip("\n")


def test_instances_and_member_programs(cls):
    cls(".coord = .coordinate.new")
    assert show(cls, "`.coord.x' `.coord.y'") == ". ."
    cls(".coord.set 1 2")
    assert show(cls, "`.coord.x' `.coord.y' `.coord'") == "1 2 (1,2)"
    cls(".c2 = .coordinate.new 2 4\n.li = .line.new\n.li.set 1 2 2 4")
    assert show(cls, "`.c2' `.li.c1'") == "(2,4) (2,4)"
    cls(".len = .li.length\n.mid = .li.midpoint")
    assert show(cls, "`.len' `.mid' `.mid.objtype' `.len.objtype'") == \
        "2.23606797749979 (1.5,3) coordinate double"
    assert cls.rc == 0
    # the tempname object of .midpoint is gone
    assert sorted(cls.session._classes.globals) == ["c2", "coord", "len", "li", "mid"]


def test_type_mismatch_and_references(cls):
    cls(".coord = .coordinate.new 1 2\n.li = .line.new\n.newthing = 2")
    for line in ('.newthing = "new"', ".coord.x = .li.midpoint", ".li.c0 = .li.length"):
        cls(line)
        assert cls.rc == 109
    cls(".a.ref = .coord.ref\n.a.x = 99")
    assert show(cls, "`.coord.x' `.coord.ref_n'") == "99 2"
    cls(".a.ref = .a\n.a.x = 5")
    assert show(cls, "`.coord.x' `.a.x' `.coord.ref_n'") == "99 5 1"
    cls(".b = .coord\n.b.x = 0")
    assert show(cls, "`.coord.x' `.b.x'") == "99 0"


def test_classwide_and_member_arrays(cls):
    cls(".ac = .coordinate2.new\n.bc = .coordinate2.new\n.ac.x_origin = 2")
    assert show(cls, "`.bc.x_origin' `.ac.classmv.arrnels' `.ac.instancemv.arrnels'") == "2 1 1"
    assert show(cls, "0`.ac.dynamicmv.arrnels'") == "0"
    cls(".ac.classmv[1] = 7")
    assert show(cls, "`.bc.x_origin'") == "7"
    # once the last instance is gone, classwide variables start again
    cls("classutil drop .ac .bc\n.ac = .coordinate2.new")
    assert show(cls, "`.ac.x_origin'") == "0"


def test_initialization_and_arrays(cls):
    cls(".t = .todolist.new")
    assert show(cls, "`.t.n' `.t.name' `.t.list[2]' `.t.b' `.t.s3'") == \
        "0 nobody mark as done .b nobody"
    assert show(cls, "`.t.mixed.arrnels' |`.t.mixed[2]'| `.t.mixed[3].objtype'") == \
        "3 || string"
    cls('.t.list[4] = "four"\n.t.list[1].Arrdropel\n.t.list.Arrpush "five"')
    assert show(cls, "`.t.list.arrnels' `.t.list.arrindexof \"four\"' |`.t.list[1]'|") == "5 4 ||"
    cls(".t.list.Arrpop\n.t.list = {\"a\", , (1+1)}")
    assert show(cls, "`.t.list.arrnels' `.t.list[3]' `.t.list[3].objtype'") == "3 2 double"
    cls('.t.list[3] = "x"')
    assert cls.rc == 109
    cls(".t.list.Arrdropall")
    assert show(cls, "`.t.list.arrnels'") == "0"


def test_declare_keys_and_builtins(cls):
    cls(".t = .todolist.new\n.t.Declare double extra = 3\n.t.Declare more = \"x\"")
    assert show(cls, "`.t.extra' `.t.more' `.t.dynamicmv.arrnels'") == "3 x 2"
    cls(".t.Declare double extra")
    assert cls.rc == 110
    cls("local k `.t.objkey'")
    assert show(cls, "`.`k'.name'") == "nobody"
    assert show(cls, "`.t.isa' `.todolist.isa' `.t.n.isa' `.t.classname' |`.t.n.classname'|") \
        == "class classtype double todolist ||"
    assert show(cls, "|`.t.nothere'|") == "||"          # errors expand to nothing


def test_inheritance_and_super(cls):
    cls(".o = .alt.new")
    assert cls(".o.zifl") == "base 10\nalt 5\n"
    assert show(cls, "`.o.who' `.o.isofclass base' `.o.isofclass line' `.o.superclass.arrnels'") \
        == "base 1 0 1"


def test_oncopy_and_destructor(cls):
    cls("set obs 2\n.v = .varc.new\n.w = .v")
    ds = cls.session.data
    assert ds.has("vc_1") and ds.has("vc_2")
    assert list(ds.get("vc_2").data) == [2.0, 2.0]
    cls("classutil drop .v")
    assert not ds.has("vc_1")
    cls("discard")
    assert not ds.has("vc_2") and not cls.session._classes.globals


def test_global_and_local_prefixes(cls, tmp_path):
    (tmp_path / "ctx.class").write_text(
        "class ctx {\n double q = 1\n}\nprogram .test\n"
        " display \"`.q' `.Global.q' |`.Local.none'|\"\n .Global.q = 42\nend\n")
    cls(".q = 7\n.c = .ctx.new")
    assert cls(".c.test") == "1 7 ||\n"
    assert show(cls, "`.q' `.c.q'") == "42 1"


def test_class_exit_outside_and_unknown_class(cls):
    cls("class exit 1")
    assert cls.rc == 198
    cls(".z = .nosuchclass.new")
    assert cls.rc == 111


def test_classutil(cls):
    cls(".coord = .coordinate.new 1 2\n.n1 = 3")
    assert cls("classutil dir") == "  .coord\n  .n1\n"
    assert cls.session.r["list"] == ".coord .n1"
    cls("classutil describe .coord")
    assert cls.session.r == {"type": "classname", "bitype": "instance"}
    cls("cutil describe coordinate")
    assert cls.session.r["bitype"] == "class"
    cls("classutil cdir coord*")
    assert cls.session.r["list"] == "coordinate coordinate2"
    out = cls("classutil which coordinate2")
    assert out.endswith("coordinate2.class\n*! version 1.0.1\n")
    cls("classutil which nosuch")
    assert cls.rc == 601
    cls("classutil drop .coord.x")
    assert cls.rc == 198
    cls("classutil drop coord n1")
    assert cls("classutil dir") == ""
