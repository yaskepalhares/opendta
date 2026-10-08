* Fase 1b: export delimited, import delimited, insheet, outsheet, type.
* Cria e apaga arquivos temporários na pasta atual.
clear
input str20 nome idade double renda str2 uf
"Ana" 30 1500.5 "SP"
"Bia, a grande" . 2000 "RJ"
"Caio" 41 0.25 "MG"
end
label define id 30 "trinta"
label values idade id
export delimited odta_p.csv, replace
type odta_p.csv
capture noisily export delimited odta_p.csv
export delimited nome renda using odta_q if idade < ., replace quote nolabel
type odta_q.csv
export delimited using odta_r.txt, replace delimiter(tab) novarnames
type odta_r.txt, showtabs
import delimited odta_p.csv, clear
describe
list
import delimited odta_p.csv, clear varnames(nonames)
describe
list
import delimited odta_p.csv, clear stringcols(2) asdouble
describe
outsheet using odta_o, replace
type odta_o.out
outsheet nome uf using odta_o.csv, comma noquote replace
type odta_o.csv
insheet using odta_o.csv, clear
list
erase odta_p.csv
erase odta_q.csv
erase odta_r.txt
erase odta_o.out
erase odta_o.csv
