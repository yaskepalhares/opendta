"""Teste de fumaça da interface (roda sem tela, com QT_QPA_PLATFORM=offscreen)."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PySide6.QtWidgets")


@pytest.fixture(scope="module")
def app():
    from PySide6.QtCore import QStandardPaths
    QStandardPaths.setTestModeEnabled(True)   # não toca nas preferências reais
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
    assert prefs.app_icon == "light"                      # padrão: ícone claro
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
    assert "\n1\n" in w.results.toPlainText()          # saída redesenhada
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
    # a divisória deixa a Command crescer
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

    # histórico com PgUp
    QTest.keyClick(w.command, Qt.Key.Key_PageUp)
    assert w.command.text() == "display 1\ndisplay 2"

    # um bloco colado roda inteiro, como trecho de do-file
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
    assert m.data(m.index(0, 1)) == "um"           # rótulo de valor
    assert m.data(m.index(1, 1)) == "2"            # sem rótulo
    w.run_command("browse nota if nota == 2")
    assert (m.rowCount(), m.columnCount()) == (1, 1)
    assert m.headerData(0, Qt.Orientation.Vertical) == "2"

    # Save (sem nome ainda) usa Save as; aqui o caminho vem direto
    w.run_command(f'save "{tmp_path / "d.dta"}"')
    w.run_command("replace nota = 3 in 1")
    assert w.save_dataset() is True                 # regrava pelo caminho absoluto
    w.run_command("clear")
    w.open_dataset(str(tmp_path / "d.dta"))
    assert w.session.data.nobs == 2
    app.processEvents()
    assert m.rowCount() == 2                        # o browser acompanha os dados
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
    assert (m.rowCount(), m.columnCount()) == (3, 3)          # linha e coluna para novos
    assert m.setData(m.index(1, 1), "9")
    app.processEvents()
    assert float(w.session.data.get("nota").data[1]) == 9
    assert ". replace nota = 9 in 2" in w.results.toPlainText()
    assert m.setData(m.index(2, 0), "Caio")                    # nova observação
    app.processEvents()
    assert w.session.data.nobs == 3 and w.session.data.get("nome").data[2] == "Caio"
    assert m.setData(m.index(0, 2), "5")                       # nova variável
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
    w.act_hints.trigger()                        # desliga pelo menu Help
    assert w.session.settings["hints"] == "off" and not w.prefs.hints
    assert ". set hints off, permanently" in w.results.toPlainText()
    w.run_command("set hints on")                # sem permanently: não grava
    assert w.act_hints.isChecked() and not w.prefs.hints
    w.act_hints.trigger()
    w.act_hints.trigger()
    assert w.prefs.hints
    w.close()
