* Gera os logs de referência no Stata 14.
*
* Como usar (no Stata 14):
*   1. cd para a pasta compat do repositório, por exemplo:
*        cd "~/opendta/compat"
*   2. do gerar_esperados.do
*
* Cada do-file de compat/do/ é executado e sua saída é gravada em
* compat/expected/<nome>.log.

version 14
set more off
set linesize 255
capture log close

* as referências usam ponto decimal; a configuração original é restaurada no fim
local dp_original = c(dp)
set dp period

local arquivos : dir "do" files "*.do"
foreach f of local arquivos {
    local nome = subinstr("`f'", ".do", "", .)
    log using "expected/`nome'.log", text replace
    capture noisily do "do/`f'"
    log close
}
set dp `dp_original'
display as result "Logs gerados em compat/expected/"
