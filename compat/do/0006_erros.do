* Mensagens de erro e códigos de retorno.
capture noisily foo
display _rc
capture noisily display nada
display _rc
capture noisily display "a" + 1
display _rc
capture noisily display 1 +
display _rc
capture noisily error 111
capture noisily error 198
capture noisily error 459
capture noisily nope(1)
capture noisily display nope(1)
display _rc
capture error 7
display _rc
capture display 1
display _rc
capture noisily confirm number abc
capture noisily confirm integer number 1.5
capture noisily quietly display nada
scalar s1 = 10
scalar s2 = "texto"
display s1 * 2
display s2
scalar list s1
capture noisily scalar drop naoexiste
