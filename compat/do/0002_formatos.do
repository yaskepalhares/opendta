* Formatos explícitos no display.
display %9.2f _pi
display %9.2f -_pi
display %09.2f _pi
display %-9.2f _pi "|"
display %12.2fc 1234567.891
display %10.3e 1234.5
display %10.0g 1/3
display %5.0g 1/3
display %8.0g 123456789
display %-10s "ab" "|"
display %10s "ab" "|"
display %td 0
display %td mdy(10, 8, 2026)
display %td mdy(2, 29, 2024)
display %9.2f .
display %9.2f .a
display "a" _col(10) "b" _col(5) "c"
display "x" _skip(3) "y"
display _dup(5) "=" " fim"
display "sem quebra " _continue
display "continuação"
display as result "resultado" as text " texto" as error " erro"
display _n "linha nova"
display _char(65) _char(66)
