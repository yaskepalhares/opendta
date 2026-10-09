* Fase 1c (continuação de 0107, que parava no label list): infile com
* :rótulo e automatic, x1-x2, infix, dicionário e file read.
capture file close _all
capture erase odta_f.raw
file open h using odta_f.raw, write replace
file write h "Ana 30 1.5" _n `""Bia Lima" . 2"' _n "Caio, 41, abc" _n "Davi 7" _n
file close h
clear
infile str10 nome idade renda:rl using odta_f.raw, automatic
list
label list rl
describe
infile str10 nome x1-x2 using odta_f.raw, clear
describe
capture noisily infile a using odta_f.raw
file open h using odta_x.raw, write replace
file write h "000123Ana       031150" _n "000124Bia Lima  .  2000" _n
file close h
infix long id 1-6 str nome 7-16 idade 17-19 renda 20-23 using odta_x.raw, clear
describe
list
file open h using odta_m.raw, write replace
file write h "1001Ana" _n "  25 1500" _n "1002Bia" _n "  31  200" _n
file close h
infix 2 lines 1: id 1-4 str nome 5-7 2: idade 1-4 renda 5-9 using odta_m.raw, clear
list
file open h using odta_d.dct, write replace
file write h "dictionary using odta_x.raw {" _n
file write h `"  _column(1) long id %6f "Identificador""' _n
file write h `"  str10 nome %10s "Nome""' _n
file write h "  int idade %3f" _n "  renda %4.1f" _n "}" _n
file close h
infile using odta_d.dct, clear
describe
list
file open r using odta_f.raw, read
file read r linha
display `"`linha'"' " eof=" r(eof)
file close r
erase odta_f.raw
erase odta_x.raw
erase odta_m.raw
erase odta_d.dct
