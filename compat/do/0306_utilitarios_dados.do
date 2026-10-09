* Fase 3c: duplicates, isid, levelsof, encode/decode, destring/tostring, split,
* mvencode/mvdecode e recode.
clear
input str8 nome x str6 v
"ana" 1 "10"
"bia" 2 "2,5"
"ana" 1 "$30"
"caio" 3 "abc"
"ana" 1 "7"
end
duplicates report
duplicates report nome x
duplicates examples nome
duplicates list nome x
duplicates tag nome, generate(dup)
list
capture noisily isid nome
isid nome v
levelsof nome
levelsof x
levelsof nome, clean separate(", ")
display "`r(levels)'"
encode nome, generate(cod)
label list cod
decode cod, generate(nome2)
describe cod nome2
destring v, generate(n1)
destring v, generate(n1) ignore("$,")
destring v, replace force
tostring x, replace
tostring x, replace
split nome2, parse(n) generate(pp)
recode cod (1=10 "um") (2/3=20) (else=.), generate(rc)
label list rc
recode dup (0=5) (nonmissing=9)
mvdecode dup, mv(5)
mvencode dup, mv(-1)
duplicates drop nome, force
list
