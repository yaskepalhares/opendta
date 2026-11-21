* Probe (survey): c() values and S_ macros that do-files use. Each block is independent.
version 14
set more off
clear
display "--- c() values"
foreach n in flavor machine_type username memory max_N_theory max_k_theory max_width_theory sortlist byteorder hostname processors mode console adosize eps epsfloat maxbyte maxint maxlong minbyte minint minlong maxfloat minfloat maxstrvarlen max_macrolen macrolen namelen cmdlen tmpdir sysdir_stata sysdir_plus sysdir_personal sysdir_base sysdir_site pagesize logtype born_date SE MP eqlen niceness noisily trace varabbrev matsize maxvar fastscroll lastcmd httpproxy checksum reventries virtual scrollbufsize rmsg emptycells adoplace tracedepth showbaselevels showemptycells showomitted fvwrap lstretch cformat pformat sformat iterlog level N k width changed filename filedate current_date current_time os osdtl dirsep rc version rng stata_version type linesize more dp graphics scheme {
    capture noisily display "c(`n') = [" c(`n') "]"
}
display "--- unknown c() name"
capture noisily display "[" c(nao_existe_xyz) "]"
display "rc = " _rc
display "--- S_ macros"
foreach m in S_DATE S_TIME S_OS S_OSDTL S_FN S_FNDATE S_level S_ADO S_CONSOLE S_MODE S_StataMP S_StataSE S_StataIC {
    display "`m' = [${`m'}]"
}
display "--- creturn list"
capture noisily creturn list
