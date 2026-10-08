* Macros locais e globais.
local a hello
local b "hello world"
local c `"com "aspas" dentro"'
display "`a'"
display "`b'"
display `"`c'"'
global g 7
display $g + ${g}
local i 2
local v2 segundo
display "`v`i''"
local n 1
display `++n'
display `n'
local m 5
display `m--'
display `m'
local ++n
display `n'
local --n
display `n'
display "[`naoexiste']"
display "`=2*21'"
local w : word count a "b c" d
display `w'
local w2 : word 2 of a b c
display "`w2'"
local f : display %5.2f 1/3
display "[`f']"
local lista1 x y z
local lista2 y z w
local u : list lista1 | lista2
display "`u'"
local d : list lista1 - lista2
display "`d'"
local s : list sizeof lista1
display `s'
tokenize "um dois tres"
display "`1' `3'"
tokenize "a,b,c", parse(",")
display "`1'|`2'|`3'|`5'"
