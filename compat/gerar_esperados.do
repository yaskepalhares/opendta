* Produces the reference logs in Stata 14.
*
* How to use it (in Stata 14):
*   1. cd to the compat folder of the repository, for example:
*        cd "~/opendta/compat"
*   2. do gerar_esperados.do            only the cases that have no log in expected/ yet
*      do gerar_esperados.do 1919       the cases whose name starts with 1919
*      do gerar_esperados.do all        every case
*
* Each selected do-file in compat/do/ is run and its output is written to
* compat/stata/<name>.log. That folder is not tracked by git, so a new run never
* blocks git pull; the logs reach compat/expected/ (the versioned references)
* through a pull request, after they are reviewed.
*
* compat/stata/ holds only the current round: before running, the logs of the
* cases that already have a reference in expected/ are deleted, so what is left
* after the run is exactly the new logs to review.

version 14
args padrao
set more off
set linesize 255
capture log close

* the references use a decimal point; the original setting is restored at the end
local dp_original = c(dp)
set dp period

capture mkdir "stata"
* logs of earlier rounds whose reference is already in expected/ are not needed
local antigos : dir "stata" files "*.log"
foreach f of local antigos {
    if fileexists("expected/`f'") capture erase "stata/`f'"
}
local glob "`padrao'*.do"
if "`padrao'" == "all" local glob "*.do"
local todos : dir "do" files "`glob'"
local arquivos
foreach f of local todos {
    local nome = subinstr("`f'", ".do", "", .)
    * without an argument, a case that already has a reference log is not run again
    if "`padrao'" != "" | !fileexists("expected/`nome'.log") {
        local arquivos `"`arquivos' "`f'""'
    }
}
display as text `: word count `arquivos'' " case(s) to run:" `"`arquivos'"'
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
