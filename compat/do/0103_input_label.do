* Fase 1a: input, label, rename, order, format, drop/keep, compress.
clear
input str8 nome idade sexo
"Ana" 30 2
"Bruno" 25 1
"Carla" . 2
"Davi" 41 1
end
label define sx 1 "Masculino" 2 "Feminino"
label values sexo sx
label list sx
list
list, nolabel
capture noisily label define sx 3 "Outro"
label define sx 3 "Outro", add
label list
rename idade age
order sexo
format age %5.1f
describe
list
drop if missing(age)
keep nome age
compress
describe
capture noisily drop naoexiste
capture noisily gen nome = 1
capture noisily gen n2 = nome + 1
