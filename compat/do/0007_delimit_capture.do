* #delimit, capture e do-files aninhados.
#delimit ;
display "a"
    "b";
forvalues k = 1/2 {;
    display `k';
};
#delimit cr
display "de volta"
capture noisily error 7
capture {
    display "não aparece"
    error 111
}
display _rc
capture noisily {
    display "aparece"
    error 198
}
display _rc
