* Fase 3d: dados para reproduzir rnormal() e as outras distribuições do
* Stata. Os runiform() da mesma semente mostram os uniformes disponíveis;
* comparando, dá para descobrir quantos cada função consome e como.
clear
set obs 8
set seed 123
gen double u = runiform()
set seed 123
gen double z = rnormal()
set seed 123
gen double k = runiformint(1, 6)
set seed 123
gen double b = rbinomial(10, .5)
set seed 123
gen double p = rpoisson(3)
set seed 123
gen double c = rchi2(2)
set seed 123
gen double e = rexponential(1)
set seed 123
gen double t = rt(5)
format u-t %20.17f
list u z k, noobs
list b p c, noobs
list e t, noobs
set seed 123
display %20.17f rnormal() _n %20.17f runiform()
set seed 123
display %20.17f rnormal() _n %20.17f rnormal() _n %20.17f runiform()
