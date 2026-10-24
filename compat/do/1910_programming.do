* Programming, second part: makecns with a list of constraints and its displaycns option,
* matrix dispCns, matlist with equation names (showcoleq(), roweqonly, coleqonly, keepcoleq,
* lines(eq)), matrix rownames/colnames with equation names, and matrix list and matlist of a 1x1
* matrix defined by the user. The window dialogs and menus are left out: they need someone at the
* screen. Every command whose behavior the manual does not settle runs under capture noisily, so
* that one Stata run shows the whole case.
clear
set obs 60
gen x1 = mod(_n*7, 11)
gen x2 = mod(_n*5, 13)
gen x3 = mod(_n*3, 7)
gen y = 2 + x1 + x2 - x3 + mod(_n*11, 9)

* constraint without the word define, and constraint dir
constraint 1 x1 = x2
constraint define 2 _cons = 1
constraint 3 2*x1 = 2*x2
constraint 4 zz = 0
constraint 5 x3 + x1 = 0.5
capture noisily constraint dir

* makecns with a list of constraints
quietly regress y x1 x2 x3
capture noisily makecns 1-5, displaycns r
capture noisily return list
capture noisily matrix list e(Cns)
capture noisily matrix dispCns
capture noisily matrix dispCns, r
capture noisily return list
capture noisily makecns 1 4, nocnsnotes r
capture noisily return list
capture noisily makecns 3 1
capture noisily matrix dispCns
capture noisily makecns 99
capture noisily display _rc
capture noisily matrix dispCns
capture noisily makecns 2 5, displaycns
matrix C = (1, -1, 0, 0, 0)
capture noisily makecns C, displaycns r
capture noisily return list

* equations: makecns after ereturn post of a two-equation vector
matrix b = (1, 2, 3, 4)
matrix colnames b = eq1:x1 eq1:_cons eq2:x1 eq2:_cons
capture noisily matrix list b
ereturn post b
constraint 6 [eq1]x1 = [eq2]x1
constraint 7 [eq2]_cons = 2
capture noisily makecns 6 7, displaycns r
capture noisily return list
capture noisily matrix list e(Cns)

* matlist with equation names
matrix A = (1,2,3,4 \ 2,5,6,7 \ 3,6,8,9 \ 4,7,9,10)
matrix rownames A = y:x1 y:_cons z:x1 z:_cons
matrix colnames A = y:x1 y:_cons z:x1 z:_cons
capture noisily matrix list A
capture noisily matlist A
capture noisily matlist A, lines(eq)
capture noisily matlist A, showcoleq(first) nohalf
capture noisily matlist A, showcoleq(each)
capture noisily matlist A, showcoleq(lcombined) lines(eq)
capture noisily matlist A, showcoleq(rcombined) border(all)
capture noisily matlist A, roweqonly
capture noisily matlist A, coleqonly
capture noisily matlist A, keepcoleq linesize(40)
capture noisily matlist A, linesize(40)
capture noisily matlist A, colorcoleq(res) lines(eq)
capture noisily matlist A, showcoleq(sideways)
matrix B = (1, 2, 3 \ 4, 5, 6)
matrix colnames B = a:u a:v b:w
capture noisily matlist B
capture noisily matlist B, lines(eq) showcoleq(each)
matrix colnames B = u v w
capture noisily matrix list B

* a 1x1 matrix defined by the user
matrix S = (7)
capture noisily matrix list S
capture noisily matlist S
matrix S2 = (7, 1 \ 1, 3)
capture noisily matrix list S2
capture noisily matlist S2
