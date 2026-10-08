* Funções de string.
display "a" + "b"
display "ab" * 3
display strlen("abc")
display length("abc")
display upper("Abc") lower("ABC") proper("joão da silva")
display substr("abcdef", 2, 3)
display substr("abcdef", -2, .)
display substr("abcdef", 10, 2) "|"
display strpos("abcabc", "c")
display subinstr("aaa", "a", "b", 2)
display subinstr("aaa", "a", "b", .)
display subinword("a b a b", "a", "x", .)
display word("um dois tres", 2)
display word("um dois tres", -1)
display wordcount("  um   dois  ")
display trim("  x  ") "|"
display ltrim("  x  ") "|"
display itrim("a    b")
display reverse("abc")
display real("1.5") + 1
display real("abc")
display strmatch("abc", "a*")
display strmatch("abc", "?b?")
display regexm("abc123", "[0-9]+")
display regexs(0)
display regexr("abc123", "[0-9]+", "X")
display "abc" < "abd"
display "B" < "a"
display abbrev("variavel_muito_longa", 10)
display inlist("b", "a", "b")
