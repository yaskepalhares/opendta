* Fase 3b: egen (funções de grupo e de linha).
clear
input g x y str4 s
1 2 5 "a b"
1 4 . "c"
1 . 7 "d e"
2 8 1 "f"
2 6 2 "g h"
2 6 3 "i"
end
egen m = mean(x), by(g)
egen sd = sd(x), by(g)
egen n = count(x), by(g)
egen t = total(x)
egen md = median(x)
egen p25 = pctile(x), p(25)
egen mo = mode(x)
egen r = rank(x)
egen rf = rank(x), field
egen z = std(x)
egen grp = group(g y)
egen tg = tag(g)
egen sq = seq(), by(g)
egen f = fill(10 20)
egen ct = cut(x), at(0 5 10)
egen cg = cut(x), group(2)
egen rm = rowmean(x y)
egen rt = rowtotal(x y)
egen rn = rownonmiss(x y)
egen rmax = rowmax(x y)
egen ac = anycount(x y), values(6 7)
egen cc = concat(g s), punct(-)
egen e1 = ends(s)
egen e2 = ends(s), last
list, abbreviate(4)
describe
bysort g: egen mx = max(y)
list g y mx
capture noisily egen m = mean(x)
capture noisily egen q = naoexiste(x)
