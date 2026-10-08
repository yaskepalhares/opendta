* Fase 1c: export excel e import excel.
* Cria e apaga arquivos temporários na pasta atual.
clear
input str10 nome idade double renda dnasc
"Ana" 30 1500.5 21915
"Bia" . 2000.25 .
"Caio" 41 . 22000
end
format dnasc %td
label define id 30 "trinta"
label values idade id
label variable renda "Renda mensal"
capture erase odta_x.xlsx
export excel using odta_x.xlsx, firstrow(variables)
capture noisily export excel using odta_x.xlsx
export excel nome renda using odta_x.xlsx if idade < ., sheet("Outra") sheetmodify firstrow(varlabels)
import excel odta_x.xlsx, describe
return list
import excel odta_x.xlsx, clear firstrow
describe
list
import excel odta_x.xlsx, clear sheet("Outra") firstrow
describe
list
import excel odta_x.xlsx, clear cellrange(A2:B4)
describe
list
import excel odta_x.xlsx, clear firstrow allstring
describe
list
import excel odta_x.xlsx, clear
describe
erase odta_x.xlsx
