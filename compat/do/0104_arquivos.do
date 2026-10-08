* Fase 1b: use, save, saveold, notes, char, describe com arquivo.
* Cria e apaga arquivos temporários na pasta atual.
clear
input str6 nome idade
"Ana" 30
"Bia" .
"Caio" 41
end
label data "Pessoas"
label variable idade "Idade (anos)"
notes: arquivo criado no teste
notes: segunda nota
notes idade: idade em anos completos
char idade[fonte] "questionario"
char list
notes
sort nome
capture erase odta_pessoas.dta
save odta_pessoas
capture noisily save odta_pessoas
save odta_pessoas, replace
capture erase odta_outro.dta
save odta_outro, replace
gen x = 1
capture noisily use odta_pessoas
use odta_pessoas, clear
describe
notes list idade
list
capture noisily use nome using odta_pessoas if idade < ., clear
list
use odta_pessoas in 2/3, clear
list
capture noisily use odta_naoexiste, clear
capture erase odta_v13.dta
saveold odta_v13
use odta_v13, clear
describe
notes drop _dta
notes
gen double d = _n + 0.5
capture noisily recast int d
recast int d, force
recast str3 nome, force
capture noisily recast byte nome
describe
list
capture erase odta_v12.dta
saveold odta_v12, version(12)
use odta_v12, clear
describe
capture erase odta_v11.dta
saveold odta_v11, version(11)
capture noisily saveold odta_v9, version(9)
erase odta_pessoas.dta
erase odta_outro.dta
erase odta_v13.dta
erase odta_v12.dta
erase odta_v11.dta
capture noisily erase odta_pessoas.dta
