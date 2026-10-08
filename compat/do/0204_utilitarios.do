* Fase 2: assert, log com nome, set trace, creturn de valores fixos.
clear
set obs 5
gen x = _n
assert x > 0
capture noisily assert x < 3
assert x < 3, rc0
capture erase odta_log.log
log using odta_log, text name(teste)
display "dentro do log"
log close teste
type odta_log.log
erase odta_log.log
capture program drop tr
program define tr
    local a 5
    display `a'
end
set trace on
tr
set trace off
display c(pi)
display c(maxdouble)
display "`c(dirsep)'"
