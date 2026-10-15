"""Phase 9: interface language, menu dialogs and Graph Editor."""

import os

import pytest

from opendta import i18n
from opendta.session import Session

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(autouse=True)
def reset_language(tmp_path, monkeypatch):
    monkeypatch.setenv("OPENDTA_SETTINGS", str(tmp_path / "prefs.ini"))
    yield
    i18n.set_language("en")
    i18n.set_messages("stata")


def test_tr_and_accelerators():
    assert i18n.tr("&File") == "&File"
    i18n.set_language("pt")
    assert i18n.tr("&File") == "&Arquivo"
    assert i18n.tr("Open...") == "Abrir..."
    assert i18n.tr("texto sem tradução") == "texto sem tradução"
    assert i18n.L("a", "b") == "b"


def test_translate_error_only_when_asked():
    msg = "variable zz not found"
    i18n.set_language("pt")
    assert i18n.translate_error(msg) == msg          # Stata's messages stay in English
    i18n.set_messages("translated")
    assert i18n.translate_error(msg) == "variável zz não encontrada"
    assert i18n.translate_error("invalid syntax") == "sintaxe inválida"
    i18n.set_language("en")
    assert i18n.translate_error(msg) == msg


def _session_text(cmds: str) -> str:
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    s.run_text(cmds)
    return "".join(lines)


def test_set_language_hints_in_portuguese():
    out = _session_text("set obs 2\ngen x = 1\nset language pt\nsummarize y")
    assert "Não existe" not in out          # there are variables: the explanation is different
    assert "Nenhuma variável se chama" in out
    assert "variable y not found" in out    # Stata's message stays in English
    out = _session_text("set language pt\nset messages translated\nsummarize y")
    assert "variável y não encontrada" in out
    assert "Como corrigir:" in out


def test_set_language_invalid():
    out = _session_text("capture set language xx\ndisplay _rc")
    assert "198" in out


def test_dialog_commands():
    from opendta.gui import dialogs as D
    sp = D.spec("regress")
    v = D.initial_values(sp)
    v.update(dep="price", indep="mpg weight", robust=True, _if="foreign==0",
             _weight="aweight", _weightvar="w")
    assert sp.command(v) == "regress price mpg weight [aweight=w] if foreign==0, vce(robust)"
    v.update(robust=False, cluster="rep78", level="90")
    assert sp.command(v).endswith(", vce(cluster rep78) level(90)")
    assert D.missing_required(sp, {"dep": ""}) == ["Dependent variable:"]

    sp = D.spec("ttest")
    assert sp.command({"var": "mpg", "value": "20"}) == "ttest mpg == 20"
    assert sp.command({"var": "mpg", "value": "20", "by": "foreign", "unequal": True}) == \
        "ttest mpg, by(foreign) unequal"
    sp = D.spec("svymean")
    assert sp.command({"vars": "bmi", "over": "sexo", "subpop": "adulto"}) == \
        "svy, subpop(adulto): mean bmi, over(sexo)"
    sp = D.spec("svyset")
    assert sp.command({"psu": "upa", "pw": "peso", "strata": "estrato", "single": "centered"}) == \
        "svyset upa [pweight=peso], strata(estrato) singleunit(centered)"
    sp = D.spec("scatter")
    assert sp.command({"y": "y", "x": "x", "lfit": True}) == "twoway (scatter y x) (lfit y x)"
    sp = D.spec("glm")
    v = D.initial_values(sp)
    v.update(dep="y", indep="x", family="poisson", link="log", eform=True)
    assert sp.command(v) == "glm y x, family(poisson) link(log) eform"
    assert D.initial_values(D.spec("xtreg"), "xtreg , re")["est"] == "re"
    assert D.spec_for("tabulate  , chi2").ident == "tab2"


def test_dialog_commands_run():
    """The commands built by the dialogs run in the interpreter."""
    from opendta.gui import dialogs as D
    s = Session()
    s.run_text("set obs 40\ngen x = _n\ngen y = 2*x + mod(_n, 3)\ngen g = mod(_n, 2)")
    cmds = [D.spec("summarize").command({"vars": "x y", "detail": True}),
            D.spec("tab2").command({"row": "g", "col": "g", "chi2": True}),
            D.spec("regress").command({"dep": "y", "indep": "x", "robust": True}),
            D.spec("ttest").command({"var": "y", "by": "g"}),
            D.spec("tabstat").command({"vars": "x y", "stats": "mean sd n", "by": "g"})]
    for c in cmds:
        assert s.run_command(c) == 0, c


def test_graph_editor_edits():
    from opendta.graphics.model import Graph
    from opendta.gui import graph_editor as GE
    g = Graph(opts={"_xtitle": "x", "_ytitle": "y"})
    vals = GE.current_values(g)
    assert vals["xtitle"] == "x" and vals["legend"] and not vals["ylog"]
    vals.update(title="Título", xtitle="", ylog=True, legend=False, xsize=7)
    GE.apply_edits(g, vals)
    o = g.opts
    assert o["title"]["lines"] == ["Título"]
    assert o["xtitle"]["lines"] == []        # axis title cleared on purpose
    assert o["yscale"]["log"] and o["legend"]["off"] and o["xsize"] == 7
    assert "ytitle" not in o                 # what did not change is not written
    vals = GE.current_values(g)
    vals["title"] = ""
    GE.apply_edits(g, vals)
    assert "title" not in g.opts


def test_gui_language_switch():
    QtWidgets = pytest.importorskip("PySide6.QtWidgets")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841
    from opendta.gui.main_window import MainWindow
    w = MainWindow()
    assert [a.text() for a in w.menuBar().actions()][0] == "&File"
    w.prefs.language = "pt"
    w.apply_preferences()
    assert [a.text() for a in w.menuBar().actions()][0] == "&Arquivo"
    assert w.dock_review.windowTitle() == "Histórico"
    assert w.dock_review.objectName() == "Review"      # the saved layout does not depend on the language
    w.run_command("set obs 3")
    data = {w._prop_data.child(k).data(0, 256): w._prop_data.child(k).text(1)
            for k in range(w._prop_data.childCount())}
    assert data["Observations"] == "3"
    stat = [a for a in w.menuBar().actions() if a.text() == "&Estatística"][0].menu()
    first = stat.actions()[0].menu().actions()[0]
    assert first.text() == "Estatísticas descritivas..."      # opens the dialog
    w.prefs.language = "en"
    w.apply_preferences()
    assert w.dock_review.windowTitle() == "Review"
    w.close()


def test_gui_dialog_and_graph_editor():
    QtWidgets = pytest.importorskip("PySide6.QtWidgets")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841
    from opendta.gui import dialogs as D
    from opendta.gui.main_window import MainWindow
    w = MainWindow()
    w.run_command("set obs 20")
    w.run_command("gen x = _n")
    w.run_command("gen y = x^2")
    d = D.CommandDialog(D.spec("summarize"), w.session, w, run=w.run_command,
                        to_command=w.insert_command)
    d.widgets["vars"].set_value("x y")
    d.widgets["detail"].setChecked(True)
    assert d.preview.text() == "summarize x y, detail"
    d._copy()
    assert w.command.text() == "summarize x y, detail"
    d._submit()
    assert w.review.topLevelItem(w.review.topLevelItemCount() - 1).text(0) == "summarize x y, detail"

    w.run_command("scatter y x")
    gw = w.graph_window
    gw.act_edit.setChecked(True)
    vals = gw.editor.values()
    vals["title"] = "Editado"
    gw._apply_edits(vals)
    from opendta.graphics.store import store
    g = store(w.session).graphs["Graph"]
    assert g.opts["title"]["lines"] == ["Editado"]
    gw._revert_edits()
    assert "title" not in g.opts
    w.close()


def test_every_command_has_help():
    from opendta.commands.helpcmd import _read, find_help
    from opendta.commands.registry import REGISTRY
    from opendta.core.smcl import render_text
    s = Session()
    names = sorted({c.name for c in REGISTRY.values()})
    assert [n for n in names if find_help(s, n) is None] == []
    for n in names:
        assert render_text(_read(find_help(s, n)), 79).strip()
