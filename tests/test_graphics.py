"""Fase 7: gráficos (contas, mensagens, memória de gráficos e exportação)."""

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
    # k = min(sqrt(10), 10*log10(10)) = 3,16 -> 3 barras de largura (10-1)/3
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
    assert f["heights"] == [3.0, 3.0, 4.0]      # o máximo entra na última barra


def test_kdensity_bandwidth_and_generate(run):
    run("clear\nset seed 7\nset obs 300\ngen x = rnormal()")
    run("kdensity x, gen(kx kd) n(40) nograph")
    x = run.session.data.get("x").data
    sd = x.std(ddof=1)
    q1, q3 = np.percentile(x, [25, 75])
    h = run.session.r["bwidth"]
    # mesma regra; os quartis do Stata podem diferir um pouco dos do numpy
    assert h == pytest.approx(0.9 * min(sd, (q3 - q1) / 1.349) / 300 ** 0.2, rel=0.02)
    kd = run.session.data.get("kd").data[:40]
    kx = run.session.data.get("kx").data[:40]
    # a grade vai de mínimo - h a máximo + h (compat 0701)
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
