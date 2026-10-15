* Graphics: histogram messages, kdensity with generate(), graph
* memory (graph dir/drop/rename), graph export/save/use and errors
graph drop _all
clear
set obs 74
gen x = mod(_n*37, 101)/10
gen g = mod(_n, 3)
label variable x "Variável x"

histogram x
histogram x, bin(5)
histogram x, width(2)
histogram x, width(2) start(0)
histogram g, discrete
histogram g, discrete percent
histogram x, frequency by(g)

kdensity x, generate(kx kd) n(10) nograph
display r(bwidth)
display r(n)
display "`r(kernel)'"
list kx kd in 1/10
kdensity x, kernel(gaussian) bwidth(1) at(x) generate(kg) nograph
summarize kg
kdensity x, kernel(epan2) nograph
display r(bwidth)
kdensity x, kernel(triangle) bwidth(.5) n(5) generate(tx td) nograph
list tx td in 1/5

graph drop _all
scatter x g, name(a)
line x g, sort name(b)
graph dir, memory
capture noisily scatter x g, name(a)
scatter x g, name(a, replace)
graph rename b c
graph dir, memory
graph export odta_a.png, name(a)
capture noisily graph export odta_a.png, name(a)
graph export odta_a.png, name(a) replace
graph export odta_a.pdf, name(a)
graph save a odta_a, replace
graph drop _all
graph use odta_a
graph dir, memory
capture noisily graph display nada
capture noisily graph drop nada
capture noisily scatter x, name(z)
capture noisily twoway foo x g
capture noisily histogram x g
twoway (scatter x g) (lfit x g), name(d)
graph bar (mean) x, over(g) name(e)
graph box x, over(g) name(f)
graph pie, over(g) name(p)
graph dir, memory
erase odta_a.png
erase odta_a.pdf
erase odta_a.gph
graph drop _all
