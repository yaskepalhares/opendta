"""Solid table lines in the Results window: which characters become lines (pure logic) and the
Appearance switch."""

import os

import pytest

from opendta.gui.rules import DOWN, LEFT, RIGHT, UP, is_rule, table_cells

HEADER = "    Variable |        Obs        Mean"
RULE = "-------------+---------------------------"
ROW = "         bmi |        250    26.18413"


def test_is_rule():
    assert is_rule(RULE) and is_rule("  +-----+  ") and is_rule("|-----|") and is_rule("---")
    assert not is_rule(ROW) and not is_rule("--") and not is_rule("") and not is_rule("+++")
    assert not is_rule("-5 -3 -2")


def test_summarize_table():
    # header and rows: the bar is a full-height vertical line
    assert table_cells(None, HEADER, RULE) == [(13, UP | DOWN)]
    assert table_cells(RULE, ROW, ROW) == [(13, UP | DOWN)]
    # the rule: dashes, and a crossing where the bars meet it
    cells = dict(table_cells(HEADER, RULE, ROW))
    assert cells[0] == LEFT | RIGHT and cells[13] == UP | DOWN | LEFT | RIGHT
    assert len(cells) == len(RULE)


def test_table_ends_have_tees_and_corners():
    top = "------------+-----------"
    bottom = "------------+-----------"
    row = "        sex |      Freq."
    assert dict(table_cells(None, top, row))[12] == LEFT | RIGHT | DOWN       # ┬
    assert dict(table_cells(row, bottom, None))[12] == LEFT | RIGHT | UP      # ┴
    box_top, box_row, box_rule = "+------+", "| make |", "|------|"
    assert dict(table_cells(None, box_top, box_row))[0] == RIGHT | DOWN        # ┌
    assert dict(table_cells(box_row, box_top, None))[7] == LEFT | UP           # ┘
    assert dict(table_cells(box_row, box_rule, box_row))[0] == UP | DOWN | RIGHT   # ├
    assert dict(table_cells(box_row, box_rule, box_row))[7] == UP | DOWN | LEFT    # ┤


def test_text_that_is_not_a_table_is_left_alone():
    assert table_cells(None, "x | y", None) == []
    assert table_cells("a", "x | y", "b") == []
    assert table_cells(None, "- a negative number: -5", None) == []
    assert table_cells(None, "", None) == []


@pytest.fixture
def win(tmp_path, monkeypatch):
    monkeypatch.setenv("OPENDTA_SETTINGS", str(tmp_path / "prefs.ini"))
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    QtWidgets = pytest.importorskip("PySide6.QtWidgets")
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841
    from opendta import i18n
    from opendta.gui.main_window import MainWindow
    w = MainWindow()
    yield w
    w.close()
    i18n.set_language("en")


def test_appearance_switch_toggles_the_lines(win):
    assert win.prefs.solid_table_lines and win.results._solid_rules
    win.show_preferences("appearance")
    d = win._settings
    d.solid_lines.click()
    assert not win.prefs.solid_table_lines and not win.results._solid_rules
    d.solid_lines.click()
    assert win.prefs.solid_table_lines and win.results._solid_rules
    d.solid_lines.click()
    assert not win.prefs.solid_table_lines
    d.close()


def test_results_view_draws_and_keeps_the_text(win):
    text = f"{HEADER}\n{RULE}\n{ROW}\n"
    win.results.append_styled(text, "text")
    win.show()
    for solid in (True, False):
        win.results.set_solid_rules(solid)
        pix = win.results.grab()                # runs paintEvent
        assert not pix.isNull()
    assert RULE in win.results.toPlainText() and "|" in win.results.toPlainText()

