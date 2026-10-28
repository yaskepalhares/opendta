"""Smoke test of the interface (runs without a screen, with QT_QPA_PLATFORM=offscreen)."""

import os
import sys

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PySide6.QtWidgets")


@pytest.fixture(autouse=True)
def isolated_prefs(tmp_path, monkeypatch):
    """Preferences in a temporary .ini: the tests depend neither on the
    real preferences of whoever runs pytest (nor change them)."""
    monkeypatch.setenv("OPENDTA_SETTINGS", str(tmp_path / "opendta-prefs.ini"))


@pytest.fixture(autouse=True)
def no_modal_boxes(monkeypatch):
    """A modal "Save changes?" box blocks forever without a screen: answer Discard."""
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Discard)


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


def test_error_explanation_links_open_help(app):
    """The "help topic" at the end of an explanation is a link in Results: clicking it runs
    the help command."""
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtTest import QTest

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.session.settings["hints"] = "on"
    w.run_command("set obs 2")
    w.run_command("generate y = 1")
    w.run_command("list if y = 1")
    assert "Learn more: help operators" in w.results.toPlainText()
    ran = []
    w.run_command = ran.append          # the click runs a command; record it instead
    w.results.linkActivated.disconnect()
    w.results.linkActivated.connect(w.run_command)
    w.results.resize(900, 600)
    w.results.show()
    found = w.results.document().find("help operators")
    assert not found.isNull()
    at = w.results.textCursor()
    at.setPosition(found.selectionStart() + 3)
    w.results.setTextCursor(at)
    w.results.ensureCursorVisible()
    rect = w.results.cursorRect(at)
    QTest.mouseClick(w.results.viewport(), Qt.MouseButton.LeftButton, pos=QPoint(rect.right() + 2, rect.center().y()))
    assert ran == ["help operators"]
    w.close()


def test_results_highlights_deleted_observations(app):
    from opendta.gui.panes import ResultsView

    rv = ResultsView()
    rv.append_styled("(519,227 observations deleted)\n(0 observations deleted)\n", "text")
    doc = rv.document()

    def color_at(block: int) -> str:
        cur = rv.textCursor()
        cur.setPosition(doc.findBlockByNumber(block).position() + 1)
        cur.movePosition(cur.MoveOperation.Right, cur.MoveMode.KeepAnchor)
        return cur.charFormat().foreground().color().name()

    assert color_at(0) == rv._scheme.colors["alert"]
    assert color_at(1) == rv._scheme.colors["text"]       # nothing deleted: no highlight
    assert rv.toPlainText().startswith("(519,227 observations deleted)")


def test_review_lists_do_file_commands(app, tmp_path):
    """Settings > General > Session: with the option on, the commands of the do-files (nested ones
    too, and run) go under the do line in Review, with the failed ones marked; off, one line (Stata)."""
    from PySide6.QtGui import QTextCursor

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    inner = tmp_path / "inner.do"
    inner.write_text("display 5\n")
    f = tmp_path / "prog.do"
    f.write_text(f'* comment\nclear\nset obs 2\nforvalues i = 1/2 {{\n    display `i\'\n}}\ndo "{inner}"\nsummarize nothere\n')
    w.run_command(f'do "{f}"')
    top = w.review.topLevelItem(w.review.topLevelItemCount() - 1)
    assert top.childCount() == 0 and top.text(1) == "111"
    w.prefs.review_do_commands = True
    w.run_command(f'do "{f}"')
    top = w.review.topLevelItem(w.review.topLevelItemCount() - 1)
    kids = [top.child(k).text(0) for k in range(top.childCount())]
    assert kids == ["clear", "set obs 2", "forvalues i = 1/2 { …", f'do "{inner}"', "summarize nothere"]
    assert top.child(4).text(1) == "111" and top.child(0).text(1) == ""
    assert top.child(3).childCount() == 1 and top.child(3).child(0).text(0) == "display 5"
    assert top.child(2).data(0, 256).startswith("forvalues i = 1/2 {\n    display")
    # the filter keeps the do line when a command under it matches
    w._filter_review("display 5")
    assert not top.isHidden() and not top.child(3).isHidden() and top.child(0).isHidden()
    w._filter_review("")
    # Do-file Editor: the temporary file is shown by what it holds
    w.run_command(f'doedit "{f}"')
    ed = w.dofile_editor.current()
    cur = ed.textCursor()
    cur.setPosition(len("* comment\n"))
    cur.setPosition(len("* comment\nclear\nset obs 2"), QTextCursor.MoveMode.KeepAnchor)
    ed.setTextCursor(cur)
    w.dofile_editor.do_selection()
    top = w.review.topLevelItem(w.review.topLevelItemCount() - 1)
    assert top.text(0) == "do (selection: prog.do)" and top.childCount() == 2
    ed.document().setModified(False)
    w.dofile_editor.cleanup()
    w.prefs.review_do_commands = False
    w.close()


def test_dofile_editor_line_tools(app):
    from PySide6.QtGui import QTextCursor

    from opendta.gui.dofile_editor import indent_lines, outdent_lines, toggle_comment
    from opendta.gui.main_window import MainWindow

    assert toggle_comment(["a", "  b", ""]) == ["// a", "  // b", ""]
    assert toggle_comment(["// a", "  // b", ""]) == ["a", "  b", ""]
    assert indent_lines(["a", ""]) == ["\ta", ""]
    assert outdent_lines(["\ta", "      b", "c"]) == ["a", "  b", "c"]

    w = MainWindow()
    w.show_dofile_editor("")
    de = w.dofile_editor
    ed = de.current()
    text = "display 1\ndisplay 2 ///\n  + 3\ndisplay 4\n"
    ed.setPlainText(text)
    cur = ed.textCursor()
    cur.setPosition(len("display 1\ndisp"))
    ed.setTextCursor(cur)
    assert ed.current_line_code() == "display 2 ///\n  + 3"
    de.do_line()
    out = w.results.toPlainText()
    assert "\n5\n" in out and "\n1\n" not in out
    # comment two lines: one undo step, and the selection is kept
    cur = ed.textCursor()
    cur.setPosition(0)
    cur.setPosition(len("display 1\ndisp"), QTextCursor.MoveMode.KeepAnchor)
    ed.setTextCursor(cur)
    de.tool_actions["comment"].trigger()
    assert ed.toPlainText().startswith("// display 1\n// display 2 ///\n  + 3")
    assert ed.textCursor().hasSelection()
    de.tool_actions["comment"].trigger()
    assert ed.toPlainText() == text
    ed.undo()
    ed.undo()
    assert ed.toPlainText() == text
    assert de.tool_actions["redo"].isEnabled()
    ed.document().setModified(False)
    de.cleanup()
    w.close()


def test_toolbar_settings_button_is_last(app):
    from PySide6.QtWidgets import QToolBar

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    act = w.toolbar_actions["settings"]
    assert w.findChild(QToolBar, "Toolbar").actions()[-1] is act
    act.trigger()
    assert w._settings is not None and w._settings.isVisible()
    w._settings.close()
    w.close()


def test_color_schemes_for_results_and_editor(app):
    from opendta.gui import theme
    from opendta.gui.main_window import MainWindow

    assert theme.results_scheme("light", "monokai").background == "#272822"
    assert theme.results_scheme("light", "high_contrast").background != theme.results_scheme("dark", "high_contrast").background
    assert theme.results_scheme("light") is theme.STANDARD and theme.results_scheme("dark", "classic") is theme.DARK_RESULTS
    w = MainWindow()
    assert w.prefs.results_palette == "classic" and w.prefs.editor_palette == "classic"
    w.prefs.results_palette = "dracula"
    w.prefs.editor_palette = "nord"
    w.prefs.results_palette = "no such scheme"          # an invalid key falls back to classic
    assert w.prefs.results_palette == "classic"
    w.prefs.results_palette = "dracula"
    w.show_dofile_editor("")
    w.apply_preferences()
    assert w.results._scheme.background == "#282a36"
    ed = w.dofile_editor.current()
    assert ed.colors["bg"] == "#2e3440" and "#2e3440" in ed.styleSheet()
    w.prefs.editor_palette = "classic"
    w.apply_preferences()
    assert ed.styleSheet() == "" and "bg" not in ed.colors
    w.dofile_editor.cleanup()
    w.close()


def test_editor_reopens_last_session_files(app, tmp_path):
    from opendta.gui.main_window import MainWindow
    from opendta.gui.settings import Preferences

    a, b = tmp_path / "a.do", tmp_path / "b.do"
    a.write_text("display 1\n")
    b.write_text("display 2\n")
    prefs = Preferences()
    prefs.editor_files = [str(a), str(tmp_path / "gone.do")]
    prefs.sync()
    w = MainWindow()
    w.show_dofile_editor("")                       # off by default: nothing is reopened
    assert w.dofile_editor.open_paths() == []
    w.dofile_editor.cleanup()
    w.close()

    prefs.reopen_files = True
    prefs.sync()
    w = MainWindow()
    w.show_dofile_editor("")
    assert w.dofile_editor.open_paths() == [str(a)]      # files that no longer exist are skipped
    w.dofile_editor.open_path(b)
    assert w.prefs.editor_files == [str(a), str(b)]      # the list follows the open tabs
    w.dofile_editor.close()                              # closing the window keeps it for next time
    assert w.prefs.editor_files == [str(a), str(b)]
    w.dofile_editor.cleanup()
    w.close()


def test_initial_folder_and_automatic_log(app, tmp_path, monkeypatch):
    from opendta.gui.main_window import MainWindow
    from opendta.gui.settings import Preferences

    monkeypatch.chdir(tmp_path)
    work = tmp_path / "project"
    work.mkdir()
    w = MainWindow()                                # defaults: no initial folder, no log
    assert os.path.samefile(os.getcwd(), tmp_path) and not (work / "logs").exists()
    w.close()

    prefs = Preferences()
    prefs.start_dir = str(work)
    prefs.auto_log = True
    prefs.sync()
    w = MainWindow()
    assert os.path.samefile(os.getcwd(), work)
    logs = list((work / "logs").glob("opendta-*.log"))
    assert len(logs) == 1
    w.run_command("display 6 * 7")
    w.close()
    assert "42" in logs[0].read_text()
    assert not getattr(w.session, "logs", None)    # closed with the window


def test_confirmation_before_risky_commands(app, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    from opendta import guard
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.run_command("set obs 20000")
    w.run_command("generate year = 2000 + mod(_n, 20)")
    asked = []
    answer = {"run": False}
    with monkeypatch.context() as m:
        m.setattr(MainWindow, "_ask_risk", lambda self, risk, line: asked.append(risk.kind) or answer["run"])
        w.run_command("keep if year >= 2016")       # declined: nothing happens
        assert asked == ["drop_obs"] and w.session.data.nobs == 20000
        assert "not run" in w.results.toPlainText()
        answer["run"] = True
        w.run_command("keep if year >= 2016")
        assert w.session.data.nobs == 4000
        w.run_command("generate z = 1")             # dirty again so that clear is risky
        w.prefs.set_confirm("keep", False)          # each command has its own switch
        w.run_command("keep if year >= 2019")
        assert asked == ["drop_obs", "drop_obs"] and w.session.data.nobs == 1000    # keep: no question
        w.run_command("clear")
        assert asked == ["drop_obs", "drop_obs", "clear"] and w.session.data.nobs == 0
        w.prefs.set_confirm("clear", False)
        w.run_command("set obs 5")
        w.run_command("clear")
        assert asked[-1] == "clear" and len(asked) == 3 and w.session.data.nobs == 0
    w.close()

    # the first confirmation says where to turn it off, with a way to open Settings there
    w = MainWindow()
    seen = []
    monkeypatch.setattr(QMessageBox, "exec", lambda self: seen.append(
        (self.informativeText(), [b.text() for b in self.buttons()])) or 0)
    risk = guard.Risk("clear", "clear", "The data in memory have unsaved changes.")
    assert w._ask_risk(risk, "clear") is False
    assert w._ask_risk(risk, "clear") is False
    assert "Settings > General > Safety" in seen[0][0] and "Cancel and open Settings" in seen[0][1]
    assert "Settings" not in seen[1][0] and "Cancel and open Settings" not in seen[1][1]
    w.close()


def test_settings_new_options_and_reveal(app, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.show_preferences("general", reveal="confirm")
    dlg = w._settings
    assert all(sw.isChecked() for sw in dlg.confirms.values()) and set(dlg.confirms) == {"clear", "drop", "keep", "save"}
    assert not dlg.autolog.isChecked() and not dlg.reopen.isChecked()
    dlg.reveal("confirm_keep")
    assert dlg.stack.currentWidget() is [p for k, p, _ in dlg.pages if k == "general"][0]
    dlg.rpalette.setCurrentIndex(dlg.rpalette.findData("monokai"))
    dlg.epalette.setCurrentIndex(dlg.epalette.findData("gruvbox_dark"))
    assert (w.prefs.results_palette, w.prefs.editor_palette) == ("monokai", "gruvbox_dark")
    assert w.results._scheme.background == "#272822"
    dlg.confirms["keep"].setChecked(False)
    dlg.reopen.setChecked(True)
    dlg.autolog.setChecked(True)
    assert not w.prefs.confirm_enabled("keep") and w.prefs.confirm_enabled("drop")      # one at a time
    assert w.prefs.reopen_files and w.prefs.auto_log
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Yes)
    dlg._reset()                                    # Restore Defaults
    assert (w.prefs.results_palette, w.prefs.editor_palette) == ("classic", "classic")
    assert w.prefs.confirm_enabled("keep") and not w.prefs.reopen_files and not w.prefs.auto_log
    assert dlg.confirms["keep"].isChecked() and not dlg.reopen.isChecked() and not dlg.autolog.isChecked()
    assert dlg.rpalette.currentData() == "classic" and dlg.epalette.currentData() == "classic"
    dlg.close()
    w.close()


def test_qt_own_texts_follow_the_interface_language(app):
    from PySide6.QtCore import QCoreApplication

    from opendta import i18n
    from opendta.gui import qt_i18n

    def select_all() -> str:
        return QCoreApplication.translate("QWidgetTextControl", "Select All")

    try:
        i18n.set_language("pt")
        qt_i18n.apply("pt")
        assert select_all() != "Select All"             # Qt's right-click menu and buttons, in Portuguese
        qt_i18n.apply("en")
        assert select_all() == "Select All"
    finally:
        i18n.set_language("en")
        qt_i18n.apply("en")


def test_side_windows_open_and_return_at_their_default_place(app):
    from PySide6.QtCore import Qt

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.show()
    w.dock_review.setFloating(True)
    w.dock_variables.hide()
    w.dock_properties.hide()
    w.close()                                       # the layout is saved with the windows moved and closed

    w = MainWindow()
    w.show()
    assert all(d.isVisible() and not d.isFloating()
               for d in (w.dock_review, w.dock_variables, w.dock_properties))
    assert w.dockWidgetArea(w.dock_review) == Qt.DockWidgetArea.LeftDockWidgetArea
    assert w.dockWidgetArea(w.dock_variables) == Qt.DockWidgetArea.RightDockWidgetArea
    # closed and moved, they come back to their place when chosen again in the Window menu
    for dock, area in ((w.dock_review, Qt.DockWidgetArea.LeftDockWidgetArea),
                       (w.dock_variables, Qt.DockWidgetArea.RightDockWidgetArea),
                       (w.dock_properties, Qt.DockWidgetArea.RightDockWidgetArea)):
        dock.setFloating(True)
        dock.hide()
        dock.toggleViewAction().trigger()
        assert dock.isVisible() and not dock.isFloating() and w.dockWidgetArea(dock) == area
    w.close()


def test_closed_side_windows_get_a_toolbar_button(app):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QToolBar

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.show()
    tb = w.findChild(QToolBar, "Toolbar")
    panels = {"review": (w.dock_review, Qt.DockWidgetArea.LeftDockWidgetArea),
              "variables": (w.dock_variables, Qt.DockWidgetArea.RightDockWidgetArea),
              "properties": (w.dock_properties, Qt.DockWidgetArea.RightDockWidgetArea)}

    def shown(name: str) -> bool:
        app.processEvents()                         # the toolbar lays its buttons out on the next pass
        act = w.toolbar_actions[name]
        if sys.platform == "win32":                 # the Windows runner scales the fonts up: the last buttons
            return act.isVisible()                  # go to the overflow menu and their widget is hidden
        return tb.widgetForAction(act).isVisible()

    def separators() -> int:
        app.processEvents()
        return sum(1 for a in tb.actions() if a.isSeparator() and a.isVisible())

    # open at the start: no button for them, and the groups are still divided by one separator each
    assert not any(shown(n) for n in panels)
    open_separators = separators()
    # closed by the X of the window (or by the Window menu): the button takes its place in the toolbar
    for name, (dock, _area) in panels.items():
        dock.close()
        assert shown(name)
    assert separators() == open_separators + 1
    # the button brings the window back to its place and leaves the toolbar; Manager is untouched
    for name, (dock, area) in panels.items():
        dock.setFloating(True)
        w.toolbar_actions[name].trigger()
        assert dock.isVisible() and not dock.isFloating() and w.dockWidgetArea(dock) == area
        assert not shown(name)
    assert separators() == open_separators
    dock = w.dock_variables
    dock.toggleViewAction().trigger()               # closed from the Window menu
    assert shown("variables") and not shown("review") and not shown("properties")
    # open windows are not buttons after the main window is hidden and shown again
    dock.toggleViewAction().trigger()
    w.hide()
    w.show()
    assert not any(shown(n) for n in panels)
    w.close()


def test_side_window_buttons_follow_the_language(app):
    from opendta import i18n
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    try:
        act = w.toolbar_actions["variables"]
        assert (act.text(), act.toolTip()) == ("Variables", "Show the Variables window")
        i18n.set_language("pt")
        w.retranslate()
        assert [w.toolbar_actions[n].text() for n in ("review", "variables", "properties", "varmanager")] == \
            ["Histórico", "Variáveis", "Propriedades", "Gerenciador"]
        assert w.toolbar_actions["review"].toolTip() == "Mostrar a janela Histórico"
    finally:
        i18n.set_language("en")
        w.close()


def test_variables_manager_button_opens_the_manager(app):
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.show()
    w.toolbar_actions["varmanager"].trigger()
    assert w.varmanager is not None and w.varmanager.isVisible()
    assert w.dock_variables.isVisible()             # the side window is not the manager
    w.varmanager.close()
    w.close()


def test_dofile_editor_toolbar_layout(app):
    from PySide6.QtWidgets import QToolBar

    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    w.show_dofile_editor("")
    de = w.dofile_editor
    names = [a.data()[0] for a in de.findChild(QToolBar, "DoToolbar").actions() if a.data()]
    assert names[-3:] == ["play-line", "play-quiet", "play"]        # Do at the far right
    assert "floppy" in names and "floppy-saveas" in names
    assert de.tabs.count() == 1
    de.new_tab_button.click()                                       # the + at the end of the tab bar
    assert de.tabs.count() == 2 and de.tabs.currentIndex() == 1
    for i in (1, 0):
        de.tabs.widget(i).document().setModified(False)
    de.cleanup()
    w.close()


def test_window_commands_in_the_interface(app, monkeypatch):
    from PySide6.QtWidgets import QFileDialog, QMessageBox
    from opendta.gui.main_window import MainWindow
    from opendta.gui.window_cmds import qt_filter

    assert qt_filter("Data (*.dta)|*.dta|Text|*.txt;*.csv") == "Data (*.dta);;Text (*.txt *.csv)"
    w = MainWindow()
    w.run_command('window menu append submenu "stUser" "Tools"')
    w.run_command('window menu append item "Tools" "Say hi" "display 123"')
    w.run_command("window menu refresh")
    tools = [a.menu() for a in w.user_menu.actions() if a.menu() and a.text() == "Tools"][0]
    tools.actions()[0].trigger()
    assert "\n123\n" in w.results.toPlainText()
    w.run_command('window push display "pushed"')
    texts = [w.review.topLevelItem(k).text(0) for k in range(w.review.topLevelItemCount())]
    assert texts[-2:] == ['display "pushed"', 'window push display "pushed"']
    monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: ("/tmp/x.dta", ""))
    w.run_command('window fopen fname "Open" "Data (*.dta)|*.dta"')
    w.run_command('display "$fname"')
    assert "/tmp/x.dta" in w.results.toPlainText()
    monkeypatch.setattr(QMessageBox, "question", lambda *a, **k: QMessageBox.StandardButton.Cancel)
    w.run_command('capture window stopbox rusure "Sure?"')
    w.run_command("display _rc")
    assert w.results.toPlainText().rstrip().endswith("1")
    w.run_command('window manage maintitle "My title"')
    assert w.windowTitle() == "My title"


def test_toolbar_buttons_show_their_name_under_the_icon(app):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QToolBar, QToolButton

    from opendta import i18n
    from opendta.gui.main_window import MainWindow

    w = MainWindow()
    try:
        w.show()
        tb = w.findChild(QToolBar, "Toolbar")
        assert tb.toolButtonStyle() == Qt.ToolButtonStyle.ToolButtonTextUnderIcon
        names = ["open", "save", "print", "log", "viewer", "graph", "dofile", "dataeditor", "databrowser",
                 "varmanager", "more", "break", "settings"]
        en = ["Open", "Save", "Print", "Log", "Viewer", "Graph", "Do-file Editor", "Data Editor",
              "Data Browser", "Manager", "More", "Break", "Settings"]
        assert [w.toolbar_actions[n].text() for n in names] == en
        # the long text stays in the tooltip
        assert w.toolbar_actions["log"].toolTip().startswith("Log (begin or close)")
        # a disabled button keeps its label (the style greys it out)
        btn = tb.widgetForAction(w.toolbar_actions["more"])
        assert isinstance(btn, QToolButton) and not btn.isEnabled() and btn.text() == "More"
        # everything fits in the window, in both languages: no overflow chevron
        for lang in ("en", "pt"):
            i18n.set_language(lang)
            w.retranslate()
            app.processEvents()
            widths = {n: tb.widgetForAction(a).sizeHint().width() for n, a in w.toolbar_actions.items()}
            if sys.platform != "win32":   # the Windows runner scales the fonts up (unlike the other two)
                assert tb.sizeHint().width() <= w.width(), (lang, tb.sizeHint().width(), widths)
                # the buttons of the side windows are shown only while their window is closed
                assert all(tb.widgetForAction(a).isVisible() for n, a in w.toolbar_actions.items()
                           if n not in ("review", "variables", "properties"))
                # and with the three side windows closed, when their buttons are all there too
                docks = (w.dock_review, w.dock_variables, w.dock_properties)
                for d in docks:
                    d.close()
                app.processEvents()
                assert tb.sizeHint().width() <= w.width(), (lang, tb.sizeHint().width())
                assert all(tb.widgetForAction(a).isVisible() for a in w.toolbar_actions.values())
                for d in docks:
                    d.toggleViewAction().trigger()
        assert w.toolbar_actions["dofile"].text() == "Editor de do-files"
        assert w.toolbar_actions["databrowser"].text() == "Navegador de dados"
    finally:
        i18n.set_language("en")
        w.close()


def test_data_editor_empty_grid_and_panes(app, tmp_path, monkeypatch):
    from opendta.gui.main_window import MainWindow

    monkeypatch.chdir(tmp_path)
    w = MainWindow()
    w.run_command("edit")
    app.processEvents()
    b = w.browser
    m = b.model
    assert m.rowCount() > 10 and m.columnCount() > 3             # empty cells fill the view
    assert b.cell_name.text() == "var1[1]"
    b.table.setCurrentIndex(m.index(0, 0))
    b.cell_value.setText("4")
    b._commit_cell()
    app.processEvents()
    assert w.session.data.names == ["var1"] and w.session.data.nobs == 1
    assert m.rowCount() == 2 and m.columnCount() == 2              # no padding once there is data
    assert b.vars_list.topLevelItemCount() == 1
    w.close()


def test_do_editor_closes_empty_tab_without_asking(app):
    from opendta.gui.dofile_editor import DoFileEditor
    from opendta.session import Session

    ed = DoFileEditor(Session(), run=lambda c: None)
    ed.new_file()
    ed.current().setPlainText("   ")
    ed.current().document().setModified(True)
    assert ed.close_tab(ed.tabs.currentIndex()) is True            # blank: no confirmation
    assert ed.tab_bar.plus.x() > ed.tab_bar.tabRect(ed.tabs.count() - 1).right() - 1
