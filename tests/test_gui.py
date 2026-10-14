"""Smoke test of the interface (runs without a screen, with QT_QPA_PLATFORM=offscreen)."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PySide6.QtWidgets")


@pytest.fixture(autouse=True)
def isolated_prefs(tmp_path, monkeypatch):
    """Preferences in a temporary .ini: the tests depend neither on the
    real preferences of whoever runs pytest (nor change them)."""
    monkeypatch.setenv("OPENDTA_SETTINGS", str(tmp_path / "opendta-prefs.ini"))


@pytest.fixture(scope="module")
def app():
    from PySide6.QtCore import QStandardPaths
    QStandardPaths.setTestModeEnabled(True)   # does not touch the real preferences
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_main_window_runs_commands(app):
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.run_command("display 40 + 2")
    w.run_command("foo")
    text = w.results.toPlainText()
    assert ". display 40 + 2" in text
    assert "\n42\n" in text
    assert "r(199);" in text
    assert w.review.topLevelItemCount() == 2
    assert w.review.topLevelItem(1).text(1) == "199"
    titles = [d.windowTitle() for d in (w.dock_review, w.dock_variables, w.dock_properties)]
    assert titles == ["Review", "Variables", "Properties"]
    w.close()


def test_preferences_icon_and_font(app, tmp_path):
    from PySide6.QtCore import QSettings

    from opendta.gui.main_window import MainWindow
    from opendta.gui.settings import Preferences
    from opendta.gui.theme import default_monospace

    store = QSettings(str(tmp_path / "prefs.ini"), QSettings.Format.IniFormat)
    prefs = Preferences(store)
    assert prefs.app_icon == "light"                      # default: light icon
    assert prefs.font_family == default_monospace()[0]

    w = MainWindow()
    w.prefs = prefs
    prefs.app_icon = "dark"
    prefs.font_size = 14
    w.apply_preferences()
    assert not w.windowIcon().isNull()
    assert w.results.font().pointSize() == 14
    w.close()


def test_variables_panel_follows_data(app):
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.run_command("set obs 3")
    w.run_command("gen idade = 20 + _n")
    w.run_command('label variable idade "Idade em anos"')
    assert w.variables.topLevelItemCount() == 1
    item = w.variables.topLevelItem(0)
    assert (item.text(0), item.text(1)) == ("idade", "Idade em anos")
    item.setSelected(True)
    props = {w._prop_vars.child(k).text(0): w._prop_vars.child(k).text(1)
             for k in range(w._prop_vars.childCount())}
    assert props["Type"] == "float" and props["Format"] == "%9.0g"
    data = {w._prop_data.child(k).text(0): w._prop_data.child(k).text(1)
            for k in range(w._prop_data.childCount())}
    assert data["Observations"] == "3" and data["Variables"] == "1"
    w.close()


def test_theme_switch(app):
    from PySide6.QtGui import QPalette

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.run_command("display 1")
    w.apply_theme("dark")
    assert app.palette().color(QPalette.ColorRole.Window).lightness() < 100
    assert "#1e1f22" in w.results.styleSheet()
    assert "\n1\n" in w.results.toPlainText()          # output redrawn
    w.apply_theme("light")
    assert app.palette().color(QPalette.ColorRole.Window).lightness() > 200
    assert "#ffffff" in w.results.styleSheet()
    w.close()


def test_command_pane_resizable_and_multiline(app):
    from PySide6.QtCore import Qt
    from PySide6.QtTest import QTest

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.resize(1000, 700)
    w.show()
    app.processEvents()
    # the splitter lets the Command window grow
    w.splitter.setSizes([300, 300])
    app.processEvents()
    assert w.command.height() > 200

    # Enter executa; Shift+Enter quebra a linha
    w.command.setFocus()
    QTest.keyClicks(w.command, "display 1")
    QTest.keyClick(w.command, Qt.Key.Key_Return, Qt.KeyboardModifier.ShiftModifier)
    QTest.keyClicks(w.command, "display 2")
    assert w.command.text() == "display 1\ndisplay 2"
    QTest.keyClick(w.command, Qt.Key.Key_Return)
    text = w.results.toPlainText()
    assert ". display 1\n1\n" in text and ". display 2\n2\n" in text
    assert w.command.text() == ""

    # history with PgUp
    QTest.keyClick(w.command, Qt.Key.Key_PageUp)
    assert w.command.text() == "display 1\ndisplay 2"

    # a pasted block runs whole, as a piece of a do-file
    w.command.setText('forvalues i = 1/2 {\n    display "i = `i\'"\n}')
    QTest.keyClick(w.command, Qt.Key.Key_Return)
    text = w.results.toPlainText()
    assert "  2.     display" in text and "i = 1\ni = 2\n" in text
    last = w.review.topLevelItem(w.review.topLevelItemCount() - 1)
    assert last.text(0) == "forvalues i = 1/2 { …" and last.text(1) == ""
    w.close()


def test_open_save_and_browser(app, tmp_path, monkeypatch):
    from PySide6.QtCore import Qt

    from opendta.gui.main_window import MainWindow

    monkeypatch.chdir(tmp_path)
    w = MainWindow()
    w.run_commands('clear\ninput str4 nome nota\n"Ana" 1\n"Bia" 2\nend\n'
                   'label define n 1 "um"\nlabel values nota n')
    w.run_command("browse")
    b = w.browser
    assert b is not None and b.isVisible()
    m = b.model
    assert (m.rowCount(), m.columnCount()) == (2, 2)
    assert m.headerData(0, Qt.Orientation.Horizontal) == "nome"
    assert m.data(m.index(0, 1)) == "um"           # value label
    assert m.data(m.index(1, 1)) == "2"            # no label
    w.run_command("browse nota if nota == 2")
    assert (m.rowCount(), m.columnCount()) == (1, 1)
    assert m.headerData(0, Qt.Orientation.Vertical) == "2"

    # Save (no name yet) uses Save as; here the path comes directly
    w.run_command(f'save "{tmp_path / "d.dta"}"')
    w.run_command("replace nota = 3 in 1")
    assert w.save_dataset() is True                 # saves again through the absolute path
    w.run_command("clear")
    w.open_dataset(str(tmp_path / "d.dta"))
    assert w.session.data.nobs == 2
    app.processEvents()
    assert m.rowCount() == 2                        # the browser follows the data
    assert float(w.session.data.get("nota").data[0]) == 3
    w.close()


def test_edit_commands_unit(app):
    import numpy as np
    from opendta.core.dataset import Dataset, Variable
    from opendta.gui.data_browser import edit_commands
    ds = Dataset()
    ds.nobs = 2
    ds.vars = [Variable("x", "byte", np.array([1.0, 2.0])), Variable("s", "str5", ["a", "b"])]
    ds.value_labels = {"sim": {1: "Sim", 2: "Não"}}
    ds.vars[0].value_label = "sim"
    assert edit_commands(ds, ds.vars[0], 1, "7") == ["replace x = 7 in 2"]
    assert edit_commands(ds, ds.vars[0], 0, "Não") == ["replace x = 2 in 1"]
    assert edit_commands(ds, ds.vars[0], 0, "") == ["replace x = . in 1"]
    assert edit_commands(ds, ds.vars[1], 2, 'di"z') == ["set obs 3", 'replace s = `"di"z"\' in 3']
    assert edit_commands(ds, None, 0, "abc") == ['generate var1 = "abc" in 1']
    with pytest.raises(ValueError):
        edit_commands(ds, ds.vars[0], 0, "talvez")


def test_data_editor_and_variables_manager(app, tmp_path, monkeypatch):
    from PySide6.QtCore import Qt

    from opendta.gui.main_window import MainWindow

    monkeypatch.chdir(tmp_path)
    w = MainWindow()
    w.run_commands('clear\ninput str4 nome nota\n"Ana" 1\n"Bia" 2\nend')
    w.run_command("edit")
    b = w.browser
    m = b.model
    assert b.windowTitle() == "Data Editor (Edit)"
    assert (m.rowCount(), m.columnCount()) == (3, 3)          # row and column for new data
    assert m.setData(m.index(1, 1), "9")
    app.processEvents()
    assert float(w.session.data.get("nota").data[1]) == 9
    assert ". replace nota = 9 in 2" in w.results.toPlainText()
    assert m.setData(m.index(2, 0), "Caio")                    # new observation
    app.processEvents()
    assert w.session.data.nobs == 3 and w.session.data.get("nome").data[2] == "Caio"
    assert m.setData(m.index(0, 2), "5")                       # new variable
    app.processEvents()
    assert w.session.data.names == ["nome", "nota", "var1"]
    review = [w.review.topLevelItem(k).text(0) for k in range(w.review.topLevelItemCount())]
    assert "set obs 3" in review and "generate var1 = 5 in 1" in review
    w.run_command("browse")
    assert b.windowTitle() == "Data Editor (Browse)" and not (m.flags(m.index(0, 0)) & Qt.ItemFlag.ItemIsEditable)

    w.run_command("varmanage")
    vm = w.varmanager.model
    assert vm.rowCount() == 3 and vm.data(vm.index(1, 0)) == "nota"
    assert vm.setData(vm.index(1, 1), "Nota final")
    assert vm.setData(vm.index(1, 0), "nf")
    app.processEvents()
    d = w.session.data
    assert d.names[1] == "nf" and d.get("nf").label == "Nota final"
    assert vm.data(vm.index(1, 0)) == "nf"
    w.close()


def test_help_hints_toggle(app):
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    assert w.act_hints.isChecked()
    w.act_hints.trigger()                        # turned off through the Help menu
    assert w.session.settings["hints"] == "off" and not w.prefs.hints
    assert ". set hints off, permanently" in w.results.toPlainText()
    w.run_command("set hints on")                # without permanently: not saved
    assert w.act_hints.isChecked() and not w.prefs.hints
    w.act_hints.trigger()
    w.act_hints.trigger()
    assert w.prefs.hints
    w.close()


def test_help_superscript_toggle(app):
    from opendta.core.formats import superscript
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    assert w.act_superscript.isChecked() and superscript()
    w.run_command("display 1e10/3")
    assert "3.333×10⁹" in w.results.toPlainText()
    w.act_superscript.trigger()                  # turned off through the Help menu
    assert w.session.settings["superscript"] == "off" and not w.prefs.superscript
    assert not superscript()
    w.close()
    w2 = MainWindow()                            # remembered at the next startup
    assert not w2.act_superscript.isChecked() and not superscript()
    w2.act_superscript.trigger()
    assert w2.prefs.superscript and superscript()
    w2.close()


def test_viewer_help_and_log(app, tmp_path, monkeypatch):
    from opendta.gui.main_window import MainWindow

    monkeypatch.chdir(tmp_path)
    w = MainWindow()
    w.run_command("help describe")
    v = w.viewer
    assert v is not None and v.isVisible() and "describe" in v.windowTitle()
    assert "Describe the dataset" in v.browser.toPlainText()
    v.open_topic("hints")
    assert "set hints" in v.browser.toPlainText()
    v.go(-1)
    assert "Describe the dataset" in v.browser.toPlainText()
    w.run_command("log using t, text")
    assert w.session.logs
    w.toggle_log()                                   # fecha o log aberto
    assert not w.session.logs and (tmp_path / "t.log").exists()
    w.run_command("help naoexiste")
    assert w.session.rc == 111
    w.close()


def test_graph_window_tabs(app):
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.run_command("set obs 30")
    w.run_command("gen x = _n")
    w.run_command("gen y = x^2")
    w.run_command("scatter y x, name(a)")
    w.run_command("histogram x")
    gw = w.graph_window
    assert list(gw.canvases) == ["a", "Graph"]
    assert gw.current_name() == "Graph"
    w.run_command("graph drop a")
    assert list(gw.canvases) == ["Graph"]
    w.run_command("set graphics off")
    w.run_command("scatter y x, name(b)")
    assert "b" not in gw.canvases
    w.close()


def test_dofile_editor_runs_selection(app, tmp_path):
    from PySide6.QtGui import QTextCursor

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    f = tmp_path / "prog.do"
    f.write_text("display 1 + 1\ndisplay 2 + 2\n")
    w.run_command(f'doedit "{f}"')
    ed = w.dofile_editor.current()
    assert ed.path == f and ed.toPlainText().startswith("display 1")
    # selection in the middle of the 2nd line: runs the whole line
    cur = ed.textCursor()
    cur.setPosition(len("display 1 + 1\ndisp"))
    cur.setPosition(len("display 1 + 1\ndisplay 2"), QTextCursor.MoveMode.KeepAnchor)
    ed.setTextCursor(cur)
    w.dofile_editor.do_selection()
    text = w.results.toPlainText()
    assert "\n4\n" in text and "\n2\n" not in text
    ed.moveCursor(QTextCursor.MoveOperation.End)
    ed.insertPlainText("* x")
    assert w.dofile_editor.tabs.tabText(w.dofile_editor.tabs.currentIndex()).endswith("*")
    assert w.dofile_editor.save()
    assert "x" in f.read_text()
    assert w.dofile_editor.replace_all("display", "di") == 2
    ed.document().setModified(False)
    w.dofile_editor.cleanup()
    w.close()


def test_menu_templates_fill_command(app):
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.insert_command("histogram ")
    assert w.command.toPlainText() == "histogram "
    w.close()
