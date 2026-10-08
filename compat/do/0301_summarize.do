* Fase 3a: summarize, ci, ttest, prtest, correlate, pwcorr, centile, pctile, xtile, tabstat.
clear
input grupo x y peso
1 10 3.2 1
1 12 4.1 2
1 9 2.8 1
1 15 5.0 3
2 20 6.3 1
2 18 5.9 2
2 25 7.7 1
2 . 6.1 2
2 22 7.0 1
1 11 3.9 1
end
summarize
summarize x, detail
summarize x [fw=peso]
summarize x [aw=peso], detail
summarize x, meanonly
display r(mean) " " r(N)
summarize x y, format separator(1)
by grupo, sort: summarize x
ci x y
ci x, level(90)
ttest x == 15
ttest x, by(grupo)
ttest x, by(grupo) unequal
ttest x == y
gen d = grupo == 1
prtest d == .5
capture noisily prtest grupo == .5
correlate x y
correlate x y, covariance
pwcorr x y peso, sig obs
centile x, centile(25 50 75)
pctile p = x, nq(4)
list p in 1/4
xtile q = x, nq(3)
list x q
tabstat x y, stats(mean sd min max n) by(grupo)
tabstat x y, statistics(p25 p50 p75) columns(statistics)
