* Gera os logs de referência no Stata 14.
*
* Como usar (no Stata 14):
*   1. cd para a pasta compat do repositório, por exemplo:
*        cd "C:/Users/voce/opendta/compat"
*   2. do gerar_esperados.do
*
* Cada do-file de compat/do/ é executado e sua saída é gravada em
* compat/expected/<nome>.log. Depois, faça commit desses logs.

version 14
set more off
set linesize 255
capture log close _all

local arquivos : dir "do" files "*.do"
foreach f of local arquivos {
    local nome = subinstr("`f'", ".do", "", .)
    log using "expected/`nome'.log", text replace name(ref)
    capture noisily do "do/`f'"
    log close ref
}
display as result "Logs gerados em compat/expected/"
