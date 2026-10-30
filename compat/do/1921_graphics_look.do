* Graphics checked against Stata: the kdensity bandwidth with aweight and fweight (n is the
* number of observations with aweight, the sum of the weights with fweight), over() nested twice
* in graph bar/box, over(, sort() descending) in graph dot, blabel(, format()) and the drawing
* order of scatter and lfit. The text output is compared line by line; the graphs are exported
* to PNG (odta_*.png) so that they can be compared by eye with the OpenDTA ones.
clear
set seed 1921
set obs 600
gen sexo = 1 + (runiform() < .6)
label define sx 1 "Masculino" 2 "Feminino"
label values sexo sx
gen faixa = 1 + floor(runiform()*3)
label define fx 1 "18-34" 2 "35-54" 3 "55+"
label values faixa fx
gen cid = 1 + floor(runiform()*8)
gen idade = 18 + floor(runiform()*60)
gen imc = 20 + .05*idade + exp(.2*invnormal(runiform()))*5
gen obeso = imc >= 30
gen w = ceil(runiform()*300)

kdensity imc, nograph
display r(bwidth)
kdensity imc [aweight=w], nograph
display r(bwidth)
display r(n)
kdensity imc [fweight=w], nograph
display r(bwidth)
kdensity imc if sexo==1 [aweight=w], kernel(gaussian) nograph
display r(bwidth)

twoway (kdensity imc if sexo==1 [aweight=w]) (kdensity imc if sexo==2 [aweight=w]), legend(order(1 "Masculino" 2 "Feminino"))
graph export odta_kden.png, replace
graph bar (mean) imc [pweight=w], over(sexo) blabel(bar, format(%4.1f))
graph export odta_bar1.png, replace
graph bar (mean) obeso [pweight=w], over(faixa) over(sexo) blabel(bar, format(%4.2f))
graph export odta_bar2.png, replace
graph hbar (mean) obeso [pweight=w], over(faixa) blabel(bar, format(%4.2f))
graph export odta_hbar.png, replace
graph box imc, over(faixa) over(sexo)
graph export odta_box.png, replace
graph dot (mean) obeso [pweight=w], over(cid, sort(1) descending)
graph export odta_dot.png, replace
twoway (scatter imc idade, msize(vsmall)) (lfit imc idade)
graph export odta_scatter.png, replace
