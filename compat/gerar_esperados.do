* Produces the reference logs in Stata 14.
*
* How to use it (in Stata 14):
*   1. cd to the compat folder of the repository, for example:
*        cd "~/opendta/compat"
*   2. do gerar_esperados.do
*
* Each do-file in compat/do/ is run and its output is written to
* compat/stata/<name>.log. That folder is not tracked by git, so a new run never
* blocks git pull; the logs reach compat/expected/ (the versioned references)
* through a pull request, after they are reviewed.

version 14
set more off
set linesize 255
capture log close

* the references use a decimal point; the original setting is restored at the end
local dp_original = c(dp)
set dp period

capture mkdir "stata"
local arquivos : dir "do" files "*.do"
foreach f of local arquivos {
    local nome = subinstr("`f'", ".do", "", .)
    * deletes leftovers (odta_*) of an interrupted case before running the next one
    local restos : dir . files "odta_*"
    foreach r of local restos {
        capture erase "`r'"
    }
    * each case starts with an empty Mata (functions of one case do not leak into another)
    capture mata: mata clear
    log using "stata/`nome'.log", text replace
    capture noisily do "do/`f'"
    log close
}
set dp `dp_original'
display as result "Logs written to compat/stata/"
