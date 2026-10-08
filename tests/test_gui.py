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
