* Fase 2: funções estendidas de macro.
clear
set obs 2
gen byte idade = 1
label variable idade "Idade (anos)"
label define sx 1 "Masc" 2 "Fem"
label values idade sx
format idade %5.1f
label data "Base de teste"
char idade[fonte] "censo"
sort idade
local a : type idade
local b : format idade
local c : value label idade
local d : variable label idade
local e : data label
local f : sortedby
local g : label sx 2
local h : label (idade) 1
local i : label sx 9
local j : char idade[fonte]
display "`a'|`b'|`c'|`d'|`e'|`f'|`g'|`h'|`i'|`j'"
local s "a b a c a"
local t : subinstr local s "a" "x", all count(local n)
display "`t' (`n')"
local t : subinstr local s "a" "", word
display "[`t']"
local u : list uniq s
local v : list sizeof s
local w : word count `s'
local x : word 2 of `s'
display "`u'|`v'|`w'|`x'"
local y : permname idade
display "`y'"
capture noisily local z : type naoexiste
