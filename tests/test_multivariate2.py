"""Block 19: clustermat, manovatest, camat, the plots after ca, mca, mds and procrustes, mdslong,
the binary similarity measures of cluster, kmeans keepcenters and procrustes
transform(unrestricted)."""

import numpy as np

from opendta.session import Session


def _run(text, s=None):
    s = s or Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text("set graphics off\n" + text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def _col(s, name):
    return s.data.get(name).data.astype(float)


D5 = """matrix D = (0,1,5,6,8 \\ 1,0,4,5,7 \\ 5,4,0,1.5,3 \\ 6,5,1.5,0,2.5 \\ 8,7,3,2.5,0)
matrix rownames D = a b c d e
matrix colnames D = a b c d e
"""


def test_clustermat_matches_cluster_on_the_same_distances():
    s, _ = _run("clear\n" + D5 + "clustermat averagelinkage D, labelvar(nm)\ncluster generate g = groups(2)\n")
    assert list(s.data.get("nm").raw) == ["a", "b", "c", "d", "e"]
    assert list(_col(s, "g")) == [1, 1, 2, 2, 2]
    h = _col(s, "_clus_1_hgt")
    assert np.isclose(sorted(h[:4])[0], 1.0)
    # points on a line with these distances give the same tree with cluster
    s2, _ = _run("clear\ninput x\n0\n1\n5\n6.5\n8\nend\ncluster averagelinkage x, measure(L1)\n"
                 "cluster generate g = groups(2)\n")
    assert list(_col(s2, "g")) == [1, 1, 2, 2, 2]
    s3 = Session()
    assert s3.run_text("clear\nset obs 3\ngen k = 1\n" + D5 + "clustermat singlelinkage D\n") == 4
    s4, _ = _run("clear\nset obs 5\ngen k = _n\n" + D5 + "clustermat completelinkage D, add\n")
    assert s4.data.nobs == 5 and s4.data.has("_clus_1_id")


def test_manovatest_matches_the_manova_rows():
    base = """clear
set seed 12
set obs 60
gen a = mod(_n,3)
gen b = mod(_n,2)
gen y1 = a + rnormal()
gen y2 = b + .5*a + rnormal()
gen y3 = rnormal()
manova y1 y2 y3 = a b
"""
    s, out = _run(base)
    E, H1 = s.e["E"].data, s.e["H_1"].data
    s, out = _run("manovatest a\n", s)
    wilks = np.prod(1 / (1 + np.linalg.eigvals(np.linalg.solve(E, H1)).real))
    assert abs(s.r["stat"].data[0, 0] - wilks) < 1e-6
    assert s.r["df"] == 2 and s.r["df_r"] == 56
    s, out = _run("matrix M = (1,-1,0 \\ 0,1,-1)\nmanovatest b, ytransform(M)\n", s)
    assert s.r["E"].data.shape == (2, 2)
    assert Session().run_text(base + "manovatest c\n") == 111


def test_camat_matches_ca_on_the_data():
    s, _ = _run("""clear
input r c n
1 1 20
1 2 10
1 3 5
2 1 8
2 2 15
2 3 12
3 1 3
3 2 9
3 3 25
end
quietly ca r c [fw=n]
""")
    sv = s.e["Sv"].data.ravel().copy()
    s, out = _run("matrix F = (20,10,5 \\ 8,15,12 \\ 3,9,25)\ncamat F, rowname(r) colname(c)\n"
                  "cabiplot, origin\ncaprojection\n", s)
    assert np.allclose(s.e["Sv"].data.ravel(), sv)
    assert s.e["ca_data"] == "matrix" and s.e["Rname"] == "r"
    assert "Correspondence analysis" in out


def test_mds_plots_and_mdslong():
    s, _ = _run(D5 + "quietly mdsmat D\n")
    Y = s.e["Y"].data.copy()
    s, out = _run("""clear
input i j d
1 2 1
1 3 5
1 4 6
1 5 8
2 3 4
2 4 5
2 5 7
3 4 1.5
3 5 3
4 5 2.5
end
mdslong d, id(i j)
mdsconfig
mdsshepard
""", s)
    assert np.allclose(np.abs(s.e["Y"].data), np.abs(Y))
    assert s.e["cmd"] == "mdslong"


def test_mca_plots_and_procoverlay():
    _run("""clear
set seed 3
set obs 80
gen a = 1 + mod(_n,3)
gen b = 1 + mod(floor(_n/2),4)
gen c = 1 + (runiform() > .5)
quietly mca a b c
mcaplot
mcaplot a b, origin
mcaprojection
gen x1 = rnormal()
gen x2 = x1 + rnormal()
gen x3 = rnormal()
quietly procrustes (x1 x2) (x2 x3)
procoverlay
""")


def test_binary_measures():
    s, _ = _run("""clear
input b1 b2 b3 b4
1 1 0 0
1 0 0 0
0 1 1 1
0 0 1 1
end
cluster singlelinkage b1-b4, measure(Jaccard)
""")
    # rows 3 and 4: a = 2, b = 1, c = 0 -> Jaccard 2/3, dissimilarity 1/3 (the first merge)
    h = _col(s, "_clus_1_hgt")
    assert np.isclose(np.nanmin(np.where(h < 1e300, h, np.nan)), 1 / 3)
    from opendta.core.cluster import binary_similarity
    X = np.array([[1, 1, 0, 0], [0, 1, 1, 1]], dtype=float)
    a, b, c, d = 1, 1, 2, 0
    assert np.isclose(binary_similarity(X, X, "matching")[0, 1], (a + d) / 4)
    assert np.isclose(binary_similarity(X, X, "Dice")[0, 1], 2 * a / (2 * a + b + c))
    assert np.isclose(binary_similarity(X, X, "Rogers")[0, 1], (a + d) / (a + d + 2 * (b + c)))


def test_keepcenters_and_unrestricted_procrustes():
    s, _ = _run("""clear
set seed 4
set obs 20
gen x = rnormal()
gen y = rnormal()
gen z = x + .5*y + .1*rnormal()
gen w = y - x + .1*rnormal()
cluster kmeans x y, k(3) start(firstk) keepcenters name(km)
""")
    assert s.data.nobs == 23
    x, km = _col(s, "x"), _col(s, "km")
    centers = x[20:]
    assert np.allclose(sorted(centers), sorted([x[:20][km[:20] == g].mean() for g in (1, 2, 3)]))
    s, _ = _run("drop in 21/23\nprocrustes (z w) (x y), transform(unrestricted)\n", s)
    A = s.e["A"].data
    X = np.column_stack([_col(s, "x"), _col(s, "y"), np.ones(20)])
    B = np.linalg.lstsq(X, np.column_stack([_col(s, "z"), _col(s, "w")]), rcond=None)[0]
    assert np.allclose(A, B[:2])
    assert s.e["transform"] == "unrestricted"
