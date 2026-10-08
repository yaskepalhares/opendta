* Fase 3a: tabulate (uma e duas vias), tab1, tab2, table.
clear
input sexo resp str5 regiao
1 1 "norte"
1 2 "sul"
2 1 "norte"
2 2 "sul"
1 1 "sul"
2 3 "norte"
1 . "sul"
2 2 "norte"
1 3 "sul"
2 1 "sul"
end
label define sx 1 "masc" 2 "fem"
label values sexo sx
tabulate sexo
tabulate resp, missing
tabulate regiao, sort
tabulate sexo resp
tabulate sexo resp, row col cell
tabulate sexo resp, chi2 lrchi2 V exact
tabulate sexo resp, expected nofreq
tabulate sexo resp, gamma taub
tabulate sexo, generate(d_)
describe d_*
tabulate sexo, matcell(F) matrow(R)
matrix list F
tab1 sexo regiao
tab2 sexo resp regiao
table sexo resp
table sexo, contents(freq mean resp)
