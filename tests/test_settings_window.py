"""The Settings window: immediate apply, search, set commands and reset."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
QtWidgets = pytest.importorskip("PySide6.QtWidgets")


@pytest.fixture
def win(tmp_path, monkeypatch):
    monkeypatch.setenv("OPENDTA_SETTINGS", str(tmp_path / "prefs.ini"))
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841
    from opendta import i18n
    from opendta.gui.main_window import MainWindow
    w = MainWindow()
    yield w
    w.close()
    i18n.set_language("en")
    i18n.set_messages("stata")


def test_theme_icon_font_apply_at_once(win):
    win.show_preferences("appearance")
    d = win._settings
    d.theme.buttons["dark"].click()
    assert win.prefs.theme == "dark" and win.prefs.effective_theme == "dark"
    d.appicon.buttons["dark"].click()
    assert win.prefs.app_icon == "dark"
    d.size.setCurrentIndex(d.size.findData("14"))
    assert win.prefs.font_size == 14 and win.results.font().pointSize() == 14
    d.theme.buttons["system"].click()
    assert win.prefs.theme == "system" and win.prefs.effective_theme in ("light", "dark")
    d.close()


def test_switches_run_set_permanently(win):
    win.show_preferences("general")
    d = win._settings
    d.hints.setChecked(False)
    assert win.session.settings["hints"] == "off" and win.prefs.hints is False
    d.hints.setChecked(True)
    assert win.session.settings["hints"] == "on"
    d.numerics.buttons["precise"].click()
    assert win.prefs.numerics == "precise"
    d.superscript.setChecked(False)
    assert win.prefs.superscript is False
    d.close()


def test_search_filters_pages(win):
    win.show_preferences()
    d = win._settings
    d.search.setText("font size")
    shown = {k for k, _p, it in d.pages if not it.isHidden()}
    assert shown == {"appearance"}
    assert d.stack.currentWidget() is dict((k, p) for k, p, _i in d.pages)["appearance"]
    d.search.setText("")
    assert all(not it.isHidden() for _k, _p, it in d.pages)
    d.close()


def test_reset(win, monkeypatch):
    win.show_preferences("about")
    d = win._settings
    d.theme.buttons["dark"].click()
    monkeypatch.setattr(QtWidgets.QMessageBox, "question",
                        lambda *a, **k: QtWidgets.QMessageBox.StandardButton.Yes)
    d._reset()
    assert win.prefs.theme == "light" and win.prefs.check_updates is True
    d.close()


def test_dropped_dofile_opens_in_the_editor_unless_asked_to_run(win, tmp_path, monkeypatch):
    from PySide6.QtCore import QMimeData, QUrl
    do = tmp_path / "a.do"
    do.write_text("display 1\n")

    class Event:
        def __init__(self):
            self._m = QMimeData()
            self._m.setUrls([QUrl.fromLocalFile(str(do))])

        def mimeData(self):
            return self._m

    path = QUrl.fromLocalFile(str(do)).toLocalFile()     # Qt uses "/" also on Windows
    ran, opened = [], []
    monkeypatch.setattr(win, "run_command", ran.append)
    monkeypatch.setattr(win, "show_dofile_editor", opened.append)
    assert not win.prefs.run_dropped_dofiles
    win.dropEvent(Event())
    assert opened == [path] and ran == []
    win.prefs.run_dropped_dofiles = True
    win.dropEvent(Event())
    assert ran == [f'do "{path}"'] and opened == [path]
