"""Teste de fumaça da interface (roda sem tela, com QT_QPA_PLATFORM=offscreen)."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PySide6.QtWidgets")


@pytest.fixture(scope="module")
def app():
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
