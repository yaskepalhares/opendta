"""Phase 7: graphics (computations, messages, graph memory and export)."""

import math
import os

import numpy as np
import pytest

from opendta.graphics import data as D
from opendta.graphics.store import store
from opendta.graphics.twoway import split_plots

SEQ = """clear
set obs 10
gen x = _n
gen y = 2*x + 1
gen g = mod(_n, 2)
"""


def test_histogram_message_default_bins(run):
    run(SEQ)
    # k = min(sqrt(10), 10*log10(10)) = 3.16 -> 3 bins of width (10-1)/3
    assert run("histogram x").strip() == "(bin=3, start=1, width=3)"


def test_histogram_message_discrete(run):
    run(SEQ)
    assert run("histogram g, discrete").strip() == "(start=0, width=1)"


def test_histogram_heights_density_integrate_to_one():
    x = np.arange(1.0, 11.0)
    h = D.hist_bins(x, None, {})
    assert h["bins"] == 3
    assert sum(h["heights"]) * h["width"] == pytest.approx(1.0)
    f = D.hist_bins(x, None, {"scale": "frequency"})
    assert f["heights"] == [3.0, 3.0, 4.0]      # the maximum goes into the last bin


def test_kdensity_bandwidth_and_generate(run):
    run("clear\nset seed 7\nset obs 300\ngen x = rnormal()")
    run("kdensity x, gen(kx kd) n(40) nograph")
    x = run.session.data.get("x").data
    sd = x.std(ddof=1)
    q1, q3 = np.percentile(x, [25, 75])
    h = run.session.r["bwidth"]
    # same rule; Stata's quartiles may differ slightly from numpy's
    assert h == pytest.approx(0.9 * min(sd, (q3 - q1) / 1.349) / 300 ** 0.2, rel=0.02)
    kd = run.session.data.get("kd").data[:40]
    kx = run.session.data.get("kx").data[:40]
    # the grid goes from minimum - h to maximum + h (compat 0701)
    assert kx[0] == pytest.approx(x.min() - h, rel=1e-12)
    assert kx[-1] == pytest.approx(x.max() + h, rel=1e-12)
    assert np.trapezoid(kd, kx) == pytest.approx(1.0, abs=0.05)
    assert run.session.r["n"] == 40


def test_epanechnikov_kernel_integrates_to_one():
    z = np.linspace(-8, 8, 40001)
    for k in D.KERNELS:
        assert np.trapezoid(D._kernel(k)(z), z) == pytest.approx(1.0, abs=1e-3), k


def test_box_stats_adjacent_values():
    x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 100.0])
    b = D.box_stats(x)
    assert b["med"] == 5.5
    assert b["hi"] == 9.0
    assert b["out"] == [100.0]


def test_split_plots_forms():
    assert split_plots("(scatter y x) (lfit y x), title(a)") == (["scatter y x", "lfit y x"], "title(a)")
    assert split_plots("scatter y x || lfit y x ||, title(a)") == (["scatter y x", "lfit y x"], "title(a)")
    assert split_plots('scatter y x, title("a || b")') == (['scatter y x, title("a || b")'], "")


def test_twoway_layers_and_by(run):
    run(SEQ)
    run("twoway (scatter y x) (lfit y x), by(g)")
    g = store(run.session).get()
    assert [p.title for p in g.panels] == ["0", "1"]
    fit = g.panels[0].layers[1]
    assert fit.kind == "lfit"
    assert fit.y[0] == pytest.approx(5.0) and fit.y[-1] == pytest.approx(21.0)
    assert g.opts["note"]["lines"] == ["Graphs by g"]


def test_function_plot(run):
    run("twoway function y = x^2, range(0 2) n(5)")
    ly = store(run.session).get().panels[0].layers[0]
    assert ly.x == [0.0, 0.5, 1.0, 1.5, 2.0]
    assert ly.y == pytest.approx([0, 0.25, 1, 2.25, 4])


def test_graph_bar_means(run):
    run(SEQ)
    run("graph bar (mean) y (count) x, over(g)")
    ex = store(run.session).get().panels[0].layers[0].extra
    assert ex["cats"] == ["0", "1"]
    assert ex["series"][0]["values"] == pytest.approx([13.0, 11.0])
    assert ex["series"][1]["values"] == [5.0, 5.0]
    assert ex["series"][0]["label"] == "mean of y"


def test_graph_memory_and_names(run):
    run(SEQ)
    run("scatter y x, name(a)")
    run("line y x, name(b)")
    out = run("scatter y x, name(a)")
    assert run.rc == 110 and "already exists" in out
    run("scatter y x, name(a, replace)")
    assert run.rc == 0
    assert run("graph dir, memory").strip() == "a  b"
    run("graph rename b c")
    run("graph drop a")
    assert list(store(run.session).graphs) == ["c"]
    run("graph drop _all")
    out = run("graph export x.png")
    assert run.rc == 111


def test_export_save_use(run, tmp_path):
    run(SEQ)
    run("histogram x, name(h)")
    png = tmp_path / "h.png"
    out = run(f'graph export "{png}"')
    assert out.strip() == f"(file {png} written in PNG format)"
    assert png.stat().st_size > 1000
    out = run(f'graph export "{png}"')
    assert run.rc == 602
    gph = tmp_path / "h"
    out = run(f'graph save h "{gph}"')
    assert out.strip() == f"(file {gph}.gph saved)"
    run("graph drop _all")
    run(f'graph use "{gph}"')
    assert store(run.session).current == "h"
    for ext in ("pdf", "svg"):
        p = tmp_path / f"h.{ext}"
        assert run(f'graph export "{p}"').strip().endswith(f"{ext.upper()} format)")


def test_render_all_families(run, tmp_path):
    run(SEQ)
    for i, cmd in enumerate(["graph box y, over(g)", "graph hbar (sum) y, over(g) blabel(bar)",
                             "graph pie, over(g) plabel(_all percent)", "graph matrix x y g, half",
                             "twoway rcap y x x || scatter y x", "scatter y x, name(a)",
                             "line y x, name(b)", "graph combine a b, cols(1)",
                             'scatter y x, title("t") subtitle("s") note("n") xlabel(0(2)10) ylog']):
        run(cmd)
        assert run.rc == 0 or cmd.endswith("ylog"), cmd
        if run.rc == 0:
            p = tmp_path / f"g{i}.png"
            run(f'graph export "{p}"')
            assert p.exists(), cmd


def test_graphics_off_and_scheme(run):
    shown = []
    run.session.ui_hooks["graph"] = shown.append
    run(SEQ)
    run("scatter y x")
    run("set graphics off")
    run("scatter y x")
    assert len(shown) == 1
    assert run.eval('c(graphics)') == "off"
    assert run.eval('c(scheme)') == "s2color"


def test_scalar_abbreviation_not_taken_by_scatter(run):
    run("sca z = 3")
    assert run.eval("z") == 3


def test_kdensity_bandwidth_counts_observations_with_aweight(run):
    # aweight: n is the number of observations (not the sum of the weights), so a constant
    # weight changes nothing; with fweight n is the sum of the weights
    run("clear\nset seed 7\nset obs 200\ngen x = rnormal()\ngen w = 1000")
    run("kdensity x, nograph")
    h0 = run.session.r["bwidth"]
    run("kdensity x [aweight=w], nograph")
    assert run.session.r["bwidth"] == pytest.approx(h0, rel=1e-9)
    run("kdensity x [fweight=w], nograph")
    assert run.session.r["bwidth"] == pytest.approx(h0 * 1000 ** -0.2, rel=1e-9)


def test_graph_bar_nested_over_and_sort(run):
    run("clear\nset obs 40\ngen a = mod(_n, 4)\ngen b = mod(_n, 2)\ngen y = a*10 + b")
    run("graph bar (mean) y, over(a) over(b)")
    ex = store(run.session).get().panels[0].layers[0].extra
    # first over() is the innermost: its categories repeat inside each group of the second
    assert ex["cats"] == ["0", "2", "1", "3"] and ex["levels"] == [["0", "0", "1", "1"]]
    run("graph dot (mean) y, over(a, sort(1) descending)")
    ex = store(run.session).get().panels[0].layers[0].extra
    assert ex["cats"] == ["3", "2", "1", "0"]


def test_graph_bar_blabel_format(run):
    run(SEQ)
    run("graph bar (mean) y, over(g) blabel(bar, format(%4.1f))")
    assert store(run.session).get().panels[0].layers[0].extra["blabel_fmt"] == "%4.1f"


def test_histogram_x_axis_starts_at_zero(run):
    # Stata starts the x axis of a histogram of positive data at zero (also with normal)
    from opendta.graphics.render import figure
    run("clear\nset seed 7\nset obs 300\ngen x = 20 + 5*runiform()")
    run("histogram x, normal name(h)")
    ax = figure(store(run.session).get()).axes[0]
    assert ax.get_xlim()[0] == 0
    assert ax.get_xlim()[1] > 24
    # data that already include negatives are left alone
    run("histogram x, name(h2)")
    run("replace x = x - 22")
    run("histogram x, name(h3)")
    ax = figure(store(run.session).get()).axes[0]
    assert ax.get_xlim()[0] < 0


def test_kdensity_bandwidth_fweight_equals_replicated_data():
    # an fweight means the observation is repeated w times (variance divisor sum(w) - 1)
    rng = np.random.default_rng(1921)
    x = np.sort(rng.uniform(size=60))   # sd < IQR/1.349
    w = rng.integers(1, 6, size=60).astype(float)
    h_fw = D.default_bandwidth(x, w, True)
    h_rep = D.default_bandwidth(np.repeat(x, w.astype(int)), None)
    assert h_fw == pytest.approx(h_rep, rel=1e-12)


def test_box_value_axis_follows_the_data_but_bars_start_at_zero(run):
    # Stata: the y axis of graph box fits the data (here about 24-32); bars start at zero
    from opendta.graphics.render import figure
    run("clear\nset seed 7\nset obs 200\ngen g = mod(_n, 2)\ngen y = 25 + 3*runiform()")
    run("graph box y, over(g)")
    lo, hi = figure(store(run.session).get()).axes[0].get_ylim()
    assert 24 < lo < 25.1 and 28 < hi < 29
    run("graph bar (mean) y, over(g)")
    assert figure(store(run.session).get()).axes[0].get_ylim()[0] == 0
    run("graph hbox y, over(g)")
    lo, hi = figure(store(run.session).get()).axes[0].get_xlim()
    assert lo > 24
