"""replace with time-series operators on the variable being replaced is sequential, as in
Stata: each observation sees the values already replaced in the earlier ones."""


def test_replace_lag_is_sequential(run):
    run("clear\nset obs 5\ngen t = _n\ntsset t\ngen y = 1\nreplace y = 2*L.y in 2/l")
    assert run("list y, noobs clean").split()[1:] == ["1", "2", "4", "8", "16"]
    run("gen w = 1\nreplace w = D.w + 1 in 2/l")
    assert run("list w, noobs clean").split()[1:] == ["1", "1", "1", "1", "1"]
    run("gen v = _n\nreplace v = L2.v + v in 3/l")
    assert run("list v, noobs clean").split()[1:] == ["1", "2", "4", "6", "9"]
