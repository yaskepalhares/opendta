* Stata's class system ([P] class, [P] class exit, [P] classutil)
* The .class files are in compat/do/classes; Stata runs the cases from compat/ and OpenDTA from
* compat/out/, so both forms of the folder go on the adopath.
discard
clear
capture quietly adopath ++ "do/classes"
capture quietly adopath ++ "../do/classes"

* -- instances, member programs, return values ([P] class, sections 1, 4.4, 6)
.coord = .oc_coord.new
display "`.coord.x' `.coord.y'"
.coord.set 1 2
display "`.coord.x' `.coord.y'"
display "`.coord'"
.c2 = .oc_coord.new 2 4
display "`.c2'"
.li = .oc_line.new
.li.set 1 2 2 4
.len = .li.length
display `.len'
display "`.len'"
.mid = .li.midpoint
display "`.mid'"
display `.li.length'
.li.length
.coord.x = .li.length
display "`.coord.x'"
* Not tested: in Stata 14 an object returned by class exit (.li.midpoint) assigned to a declared
* member (.li.c0), or to an object that already exists, is left dangling: its members vanish,
* later calls print garbage and Stata segfaults.
.l2 = .oc_line.new
.l2.set 1 2 3 4
.l2.c0 = .coord
display "`.l2.c0'"
capture noisily .l2.c0.x = 7
display _rc
.l3 = .oc_line.new
.l3.c1 = .oc_coord.new 5 6
display "`.l3.c1'"
capture noisily .l3.c1.x = 7
display _rc
.l4 = .oc_line.new
.l4.set 1 2 3 4
capture noisily .l4.c0 = .l4.back
display _rc
display "`.l4.c0'"
capture noisily .l4.c0.x = 7
display _rc

* -- built-in functions
display "`.mid.objtype' `.mid.isa' `.oc_coord.isa' `.mid.classname' `.len.objtype' `.len.isa'"
display "|`.len.classname'| `.mid.isofclass oc_coord' `.mid.isofclass oc_line'"
.c = .li.back
display "`.c' `.c.ref_n'"
.d = .c.copy
.d.x = 0
display "`.c' `.d'"
.e = .c.new 3 3
display "`.e' `.c'"
display ("`.e.uname'" == "`.c.uname'")
.f.ref = .e.ref
display ("`.f.uname'" == "`.e.uname'")

* -- type matching (section 7.1)
.newthing = 2
capture noisily .newthing = "new"
display _rc
capture noisily .coord.x = .li.midpoint
display _rc
capture noisily .li.c0 = .li.length
display _rc

* -- references (section 7.4)
.a.ref = .coord.ref
.a.x = 99
display "`.coord.x' `.coord.ref_n' `.a.ref_n'"
.i.ref = .a.ref
display "`.coord.ref_n'"
.a.ref = .a
.a.x = 5
display "`.coord.x' `.a.x' `.coord.ref_n' `.a.ref_n'"

* -- classwide variables and arrays of member variables (sections 4.6, 12.3)
.ac = .oc_coord2.new
.bc = .oc_coord2.new
.ac.x_origin = 2
display "`.bc.x_origin' `.bc.x'"
display "`.ac.classmv.arrnels' `.ac.instancemv.arrnels' 0`.ac.dynamicmv.arrnels'"
.ac.classmv[1] = 7
display "`.bc.x_origin'"
classutil drop .ac .bc
.ac = .oc_coord2.new
display "`.ac.x_origin'"

* -- initialization and arrays (sections 4.3, 7.2, 8.2)
.t = .oc_todo.new
display `"`.t.n' `.t.name' `.t.list[1]' `.t.list.arrnels' `.t.z' `.t.b' |`.t.s2'| `.t.s3'"'
display "`.t.mixed[3]' `.t.mixed.arrnels' `.t.mixed[4].objtype' `.t.mixed[3].objtype'"
.t.list[4] = "four"
display "`.t.list.arrnels' |`.t.list[3]'| `.t.list.arrindexof "four"' `.t.list.arrindexof "zz"'"
.t.list[1].Arrdropel
display "|`.t.list[1]'| `.t.list.arrnels'"
.t.list.Arrpush "five"
display "`.t.list.arrnels' `.t.list[5]'"
.t.list.Arrpop
.t.list.Arrpop
display "`.t.list.arrnels'"
.t.list = {"a", , (1+1), [.li.length]}
display "`.t.list.arrnels' `.t.list[1]' `.t.list[3]' `.t.list[4]' |`.t.list[2].objtype'|"
.t.list[2] = .li.c1
display "`.t.list[2]' `.t.list[2].objtype'"
capture noisily .t.list[2] = "Clear the coordinate"
display _rc
.t.list[2].Arrdropel
.t.list[2] = "Clear the coordinate"
display "`.t.list[2]'"
.t.list.Arrdropall
display "`.t.list.arrnels'"

* -- adding members dynamically (section 4.7)
.t.Declare double extra = 3
.t.Declare string who
.t.who = "me"
.t.Declare more = "inferred"
display "`.t.extra' `.t.who' `.t.more' `.t.dynamicmv.arrnels' `.t.more.objtype'"
.t.Declare oc_coord where = .oc_coord.new 5 6
display "`.t.where'"

* -- keys (section 12.1)
local k `.t.objkey'
display "`.`k'.name'"
.keep = .t.objkey
display "`.`.keep'.n'"

* -- inheritance (section 5)
.o = .oc_alt.new
.o.zifl
display "`.o.who' `.o.isofclass oc_base' `.o.isofclass oc_alt' `.o.size'"
.q = 1
.m = .oc_multi.new
.m.test
display "`.q' `.m.superclass.arrnels' `.oc_multi.superclass[1].classname'"

* -- .oncopy and .destructor (sections 4.8, 4.9)
set obs 3
generate a = _n
generate b = 2*_n
.v = .oc_varc.new a b
.w = .v
list
.w = .oc_varc.new
describe, short
classutil drop .v
describe, short
discard
describe, short

* -- classutil ([P] classutil)
.coord = .oc_coord.new 1 2
.li = .oc_line.new
.n1 = 3
.s1 = "text"
capture noisily classutil dir
capture noisily classutil dir, detail
capture noisily classutil dir c*
return list
capture noisily classutil describe .coord
return list
capture noisily classutil describe .li, recurse
capture noisily classutil describe oc_coord2, newok
return list
capture noisily classutil describe .n1
return list
classutil cdir oc_c*
return list
quietly classutil which oc_coord2
display strpos(`"`r(fn)'"', "oc_coord2.class") > 0
capture noisily classutil which nosuchclass
display _rc
cutil drop .n1 .s1
capture noisily classutil drop .li.c0
display _rc
capture noisily classutil drop .nothere
display _rc

* -- errors the manual does not describe
capture noisily .z = .nosuchclass.new
display _rc
capture noisily .coord.nomember = 1
display _rc
capture noisily display "|`.coord.nomember'|"
capture noisily class exit 1
display _rc
capture noisily .coord.Arrpush "x"
display _rc
capture noisily .coord.Declare double x
display _rc
capture noisily .coord
display _rc

discard
capture quietly adopath - "do/classes"
capture quietly adopath - "../do/classes"
