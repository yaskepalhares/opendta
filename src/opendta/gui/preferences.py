"""Settings window (Edit > Preferences, Ctrl+,).

A sidebar with a search field and the pages (General, Appearance, Results,
Do-file Editor, Updates, About); each page has sections made of rows: a title and a short
description on the left, the control on the right (segmented choices,
switches, menus, buttons). Changes apply at once, as in current desktop
apps, so there is no OK or Cancel; Esc or the close button closes the
window. The colors follow the light or dark theme of the interface.
"""

from __future__ import annotations

import html
import os
import tempfile
import time
from pathlib import Path

from PySide6.QtCore import (Property, QEasingCurve, QPropertyAnimation, QRectF, QSize, Qt, QTimer,
                            Signal)
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter
from PySide6.QtWidgets import (QAbstractButton, QButtonGroup, QComboBox, QDialog, QFileDialog,
                               QFontComboBox, QFrame, QHBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem,
                               QMessageBox, QPushButton, QScrollArea, QSizePolicy, QStackedWidget,
                               QToolButton, QVBoxLayout, QWidget)

from .. import __version__
from ..i18n import LANGUAGES, L, tr
from .dofile_editor import editor_colors
from .icons import app_icon, icon, is_dark
from .palettes import labels as palette_labels
from .settings import Preferences
from .theme import default_monospace, results_scheme

# colors of the window in each theme (the same neutrals as theme.py)
COLORS = {
    "light": {"content": "#ffffff", "sidebar": "#f4f4f6", "text": "#1d1d1f", "muted": "#6e6e73",
              "sep": "#e8e8ec", "control": "#f2f2f5", "border": "#dcdce1", "segsel": "#ffffff",
              "navsel": "#e4e4e9", "hover": "#ebebef", "accent": "#0a84ff", "off": "#d1d1d6",
              "danger": "#d70015"},
    "dark": {"content": "#1e1f22", "sidebar": "#26272b", "text": "#e5e5ea", "muted": "#9a9aa1",
             "sep": "#34353a", "control": "#2b2c30", "border": "#3d3e44", "segsel": "#404148",
             "navsel": "#36373c", "hover": "#303136", "accent": "#0a84ff", "off": "#505157",
             "danger": "#ff6961"},
}

_CHEVRON = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="{c}" '
            'stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M7 10l5 5 5-5"/></svg>')


def _chevron_file(color: str) -> str:
    path = Path(tempfile.gettempdir()) / f"opendta-chevron-{color.strip('#')}.svg"
    if not path.exists():
        path.write_text(_CHEVRON.format(c=color), encoding="utf-8")
    return path.as_posix()


def stylesheet(c: dict) -> str:
    chevron = _chevron_file(c["muted"])
    return f"""
#SettingsRoot {{ background: {c['content']}; }}
#Sidebar {{ background: {c['sidebar']}; border-right: 1px solid {c['sep']}; }}
#Search {{ background: {c['content']}; border: 1px solid {c['border']}; border-radius: 9px;
           padding: 7px 10px; color: {c['text']}; font-size: 13px; }}
#Search:focus {{ border: 1px solid {c['accent']}; }}
QListWidget#Nav {{ background: transparent; border: none; outline: none; font-size: 13.5px; }}
QListWidget#Nav::item {{ padding: 7px 8px; border-radius: 8px; color: {c['text']}; margin: 1px 0; }}
QListWidget#Nav::item:selected {{ background: {c['navsel']}; color: {c['text']}; }}
QListWidget#Nav::item:hover:!selected {{ background: {c['hover']}; }}
QListWidget#Nav::item:disabled {{ color: {c['muted']}; }}
QScrollArea#Page {{ background: {c['content']}; border: none; }}
#PageBody {{ background: {c['content']}; }}
QLabel#SectionTitle {{ font-size: 15px; font-weight: 600; color: {c['text']}; }}
QLabel#RowTitle {{ font-size: 13.5px; color: {c['text']}; }}
QLabel#RowDesc {{ font-size: 12px; color: {c['muted']}; }}
QLabel#Muted {{ font-size: 12px; color: {c['muted']}; }}
QLabel#BigTitle {{ font-size: 20px; font-weight: 600; color: {c['text']}; }}
QFrame#Sep {{ background: {c['sep']}; border: none; min-height: 1px; max-height: 1px; }}
QFrame#Seg {{ background: {c['control']}; border: 1px solid {c['border']}; border-radius: 10px; }}
QToolButton#SegBtn {{ background: transparent; border: 1px solid transparent; border-radius: 8px;
                      padding: 5px 12px; color: {c['muted']}; font-size: 13px; }}
QToolButton#SegBtn:checked {{ background: {c['segsel']}; border: 1px solid {c['border']}; color: {c['text']}; }}
QToolButton#SegBtn:hover:!checked {{ color: {c['text']}; }}
QToolButton#SegBtn:disabled {{ color: {c['off']}; }}
QComboBox#Menu {{ background: {c['control']}; border: 1px solid {c['border']}; border-radius: 8px;
                  padding: 5px 30px 5px 10px; color: {c['text']}; font-size: 13px; min-width: 190px; }}
QComboBox#Menu:hover {{ background: {c['hover']}; }}
QComboBox#Menu::drop-down {{ border: none; width: 26px; }}
QComboBox#Menu::down-arrow {{ image: url({chevron}); width: 14px; height: 14px; }}
QComboBox#Menu:disabled {{ color: {c['off']}; }}
QComboBox#Menu QAbstractItemView {{ background: {c['content']}; color: {c['text']};
                                     border: 1px solid {c['border']}; selection-background-color: {c['accent']};
                                     selection-color: #ffffff; outline: none; }}
QPushButton#Secondary {{ background: {c['control']}; border: 1px solid {c['border']}; border-radius: 8px;
                         padding: 6px 14px; color: {c['text']}; font-size: 13px; }}
QPushButton#Secondary:hover {{ background: {c['hover']}; }}
QPushButton#Danger {{ background: transparent; border: 1px solid {c['border']}; border-radius: 8px;
                      padding: 6px 14px; color: {c['danger']}; font-size: 13px; }}
QPushButton#Danger:hover {{ background: {c['hover']}; }}
QPushButton#Link {{ background: transparent; border: none; color: {c['accent']}; font-size: 13px;
                    padding: 2px 0; text-align: left; }}
QPushButton#Link:hover {{ text-decoration: underline; }}
QFrame#Preview {{ border: 1px solid {c['border']}; border-radius: 10px; }}
QWidget#Row[flash="true"] {{ background: rgba(10, 132, 255, 0.14); border-radius: 8px; }}
"""


# ---------------------------------------------------------------------------
# controls
# ---------------------------------------------------------------------------

class Switch(QAbstractButton):
    """An on/off switch, animated, in the accent color when on."""

    def __init__(self, checked: bool = False, parent=None):
        super().__init__(parent)
        self.setCheckable(True)
        self.setChecked(checked)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self._pos = 1.0 if checked else 0.0
        self._anim = QPropertyAnimation(self, b"knob", self)
        self._anim.setDuration(140)
        self._anim.setEasingCurve(QEasingCurve.Type.InOutCubic)
        self.toggled.connect(self._animate)
        self.colors = COLORS["light"]

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(40, 24)

    def _get(self) -> float:
        return self._pos

    def _set(self, v: float) -> None:
        self._pos = v
        self.update()

    knob = Property(float, _get, _set)

    def _animate(self, on: bool) -> None:
        self._anim.stop()
        self._anim.setStartValue(self._pos)
        self._anim.setEndValue(1.0 if on else 0.0)
        self._anim.start()

    def paintEvent(self, _event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        r = QRectF(1, 2, 38, 20)
        off, on = QColor(self.colors["off"]), QColor(self.colors["accent"])
        t = self._pos
        track = QColor(int(off.red() + (on.red() - off.red()) * t), int(off.green() + (on.green() - off.green()) * t),
                       int(off.blue() + (on.blue() - off.blue()) * t))
        if not self.isEnabled():
            track.setAlpha(110)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(track)
        p.drawRoundedRect(r, 10, 10)
        p.setBrush(QColor("#ffffff"))
        x = r.left() + 2 + t * (r.width() - 20)
        p.drawEllipse(QRectF(x, r.top() + 2, 16, 16))
        p.end()


class Segmented(QFrame):
    """Mutually exclusive choices drawn as a pill (theme, icon, numerics...)."""
    changed = Signal(str)

    def __init__(self, options, current: str, parent=None):
        super().__init__(parent)
        self.setObjectName("Seg")
        lay = QHBoxLayout(self)
        lay.setContentsMargins(3, 3, 3, 3)
        lay.setSpacing(2)
        self.group = QButtonGroup(self)
        self.group.setExclusive(True)
        self.buttons: dict[str, QToolButton] = {}
        for opt in options:
            key, label = opt[0], opt[1]
            ic = opt[2] if len(opt) > 2 else None
            b = QToolButton()
            b.setObjectName("SegBtn")
            b.setCheckable(True)
            b.setCursor(Qt.CursorShape.PointingHandCursor)
            if label:
                b.setText(label)
            if ic is not None:
                b.setIcon(ic)
                b.setIconSize(QSize(18, 18))
                b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon if label
                                     else Qt.ToolButtonStyle.ToolButtonIconOnly)
            b.setChecked(key == current)
            b.clicked.connect(lambda _c=False, k=key: self.changed.emit(k))
            self.group.addButton(b)
            self.buttons[key] = b
            lay.addWidget(b)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

    def set_current(self, key: str) -> None:
        if key in self.buttons:
            self.buttons[key].setChecked(True)


class Row(QWidget):
    """Title and description on the left, control on the right."""

    def __init__(self, title: str, desc: str, control: QWidget | None, keywords: str = ""):
        super().__init__()
        self.setObjectName("Row")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.search_text = f"{title} {desc} {keywords}".lower()
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 14, 0, 14)
        lay.setSpacing(24)
        text = QVBoxLayout()
        text.setSpacing(3)
        t = QLabel(title)
        t.setObjectName("RowTitle")
        text.addWidget(t)
        if desc:
            d = QLabel(desc)
            d.setObjectName("RowDesc")
            d.setWordWrap(True)
            text.addWidget(d)
        lay.addLayout(text, 1)
        if control is not None:
            lay.addWidget(control, 0, Qt.AlignmentFlag.AlignVCenter)


class Page(QScrollArea):
    def __init__(self, title: str):
        super().__init__()
        self.setObjectName("Page")
        self.setWidgetResizable(True)
        self.setFrameShape(QFrame.Shape.NoFrame)
        body = QWidget()
        body.setObjectName("PageBody")
        self.lay = QVBoxLayout(body)
        self.lay.setContentsMargins(36, 18, 36, 28)
        self.lay.setSpacing(0)
        self.setWidget(body)
        self.title = title
        self.sections: list[tuple[QLabel, list[QWidget]]] = []

    def section(self, title: str) -> None:
        lab = QLabel(title)
        lab.setObjectName("SectionTitle")
        lab.setContentsMargins(0, 22 if self.sections else 6, 0, 6)
        self.lay.addWidget(lab)
        self.sections.append((lab, []))

    def row(self, row: Row) -> Row:
        items = self.sections[-1][1]
        if items:
            sep = QFrame()
            sep.setObjectName("Sep")
            self.lay.addWidget(sep)
            row.sep = sep
        else:
            row.sep = None
        self.lay.addWidget(row)
        items.append(row)
        return row

    def widget_(self, w: QWidget) -> None:
        self.lay.addWidget(w)

    def finish(self) -> None:
        self.lay.addStretch(1)

    def filter(self, text: str) -> int:
        """Shows only the rows that contain the text; returns how many."""
        total = 0
        for lab, rows in self.sections:
            shown = 0
            for r in rows:
                ok = not text or text in r.search_text
                r.setVisible(ok)
                if r.sep is not None:
                    r.sep.setVisible(ok and shown > 0)
                shown += ok
            lab.setVisible(shown > 0 or not rows)
            total += shown
        return total


def _menu(items, current) -> QComboBox:
    cb = QComboBox()
    cb.setObjectName("Menu")
    for key, label in items:
        cb.addItem(label, key)
    cb.setCurrentIndex(max(0, cb.findData(current)))
    cb.setCursor(Qt.CursorShape.PointingHandCursor)
    return cb


def _button(text: str, kind: str = "Secondary") -> QPushButton:
    b = QPushButton(text)
    b.setObjectName(kind)
    b.setCursor(Qt.CursorShape.PointingHandCursor)
    return b


# ---------------------------------------------------------------------------
# the window
# ---------------------------------------------------------------------------

_PALETTE_WORDS = "cores colors paleta palette tema scheme solarized monokai dracula nord gruvbox github"

SAMPLE = (". summarize bmi age\n\n"
          "    Variable |        Obs        Mean    Std. Dev.       Min        Max\n"
          "-------------+---------------------------------------------------------\n"
          "         bmi |        250    26.18413    4.120871      17.2       38.9\n"
          "         age |        250      46.312    14.80267        18         80\n")


def _span(text: str, color: str, bold: bool = False, underline: bool = False) -> str:
    style = f"color:{color};" + ("font-weight:600;" if bold else "") + ("text-decoration:underline;" if underline else "")
    return f'<span style="{style}">{html.escape(text)}</span>'


def _results_html(rs) -> str:
    c = rs.colors
    return "<pre style='margin:0'>" + "\n".join([
        _span(". ", c["text"]) + _span("keep if year >= 2016", c["input"], True),
        _span("(519,227 observations deleted)", c["alert"], True),
        _span("variable zz not found", c["error"]),
        _span("    How to fix: ", c["hint_label"], True) + _span("check the spelling", c["hint"]),
        _span("    ", c["hint"]) + _span("help varlist", c["hint_link"], underline=True),
    ]) + "</pre>"


def _editor_html(ec: dict, fg: str) -> str:
    def code(text: str) -> str:
        return _span(text, fg)
    return "<pre style='margin:0'>" + "\n".join([
        _span("* mean BMI by sex", ec["comment"]),
        _span("use", ec["command"], True) + code(" ") + _span('"survey.dta"', ec["string"]) + code(", clear"),
        _span("local", ec["command"], True) + code(" lim = ") + _span("18", ec["number"]),
        _span("keep if", ec["command"], True) + code(" age >= ") + _span("`lim'", ec["macro"]),
        _span("summarize", ec["command"], True) + code(" bmi"),
    ]) + "</pre>"


class PreferencesDialog(QDialog):
    """The Settings window. `applied` is emitted after every change."""
    applied = Signal()

    def __init__(self, prefs: Preferences, parent=None):
        super().__init__(parent)
        self.prefs = prefs
        self.win = parent
        self.setWindowTitle(tr("Settings"))
        self.setObjectName("SettingsRoot")
        self.resize(940, 660)
        self.setMinimumSize(760, 520)
        self.switches: list[Switch] = []
        self.anchors: dict[str, tuple[str, Row]] = {}
        outer = QHBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # sidebar
        side = QFrame()
        side.setObjectName("Sidebar")
        side.setFixedWidth(232)
        sl = QVBoxLayout(side)
        sl.setContentsMargins(14, 16, 12, 14)
        sl.setSpacing(8)
        self.search = QLineEdit()
        self.search.setObjectName("Search")
        self.search.setPlaceholderText(tr("Search"))
        self.search.setClearButtonEnabled(True)
        self._search_action = self.search.addAction(icon("search"), QLineEdit.ActionPosition.LeadingPosition)
        self.search.textChanged.connect(self._filter)
        sl.addWidget(self.search)
        sl.addSpacing(6)
        self.nav = QListWidget()
        self.nav.setObjectName("Nav")
        self.nav.setIconSize(QSize(18, 18))
        self.nav.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        sl.addWidget(self.nav, 1)
        outer.addWidget(side)

        self.stack = QStackedWidget()
        outer.addWidget(self.stack, 1)
        self.pages: list[tuple[str, Page, QListWidgetItem]] = []
        self._caption(tr("Settings"))
        self._build_general()
        self._build_appearance()
        self._build_results()
        self._build_editor()
        self._build_updates()
        self._caption("OpenDTA")
        self._build_about()
        self.nav.currentRowChanged.connect(self._nav_changed)
        self.nav.setCurrentRow(1)
        self.nav.setFocus()             # the search field only gets focus when clicked (or Ctrl+F)
        self._restyle()
        ctl = getattr(parent, "updates", None)
        if ctl is not None and hasattr(ctl, "checked"):
            ctl.checked.connect(self._refresh_update_info)

    # -- structure -------------------------------------------------------------------
    def _caption(self, text: str) -> None:
        it = QListWidgetItem(text)
        it.setFlags(Qt.ItemFlag.NoItemFlags)
        f = QFont()
        f.setPointSizeF(f.pointSizeF() * 0.9)
        it.setFont(f)
        it.setData(Qt.ItemDataRole.UserRole, "caption")
        self.nav.addItem(it)

    def _page(self, key: str, title: str, icon_name: str) -> Page:
        page = Page(title)
        it = QListWidgetItem(title)
        it.setData(Qt.ItemDataRole.UserRole, key)
        it.setData(Qt.ItemDataRole.UserRole + 1, icon_name)
        self.nav.addItem(it)
        self.stack.addWidget(page)
        self.pages.append((key, page, it))
        return page

    def _nav_changed(self, row: int) -> None:
        it = self.nav.item(row)
        if it is None or it.data(Qt.ItemDataRole.UserRole) == "caption":
            return
        for _k, page, item in self.pages:
            if item is it:
                self.stack.setCurrentWidget(page)

    def select(self, key: str) -> None:
        for _k, _p, it in self.pages:
            if _k == key:
                self.nav.setCurrentItem(it)

    def reveal(self, name: str) -> None:
        """Goes to a row of a page, scrolls to it and flashes it for a moment."""
        key, row = self.anchors[name]
        self.select(key)
        for k, page, _it in self.pages:
            if k == key:
                page.ensureWidgetVisible(row, 0, 80)

        def flash(on: bool) -> None:
            try:
                row.setProperty("flash", on)
                row.style().unpolish(row)
                row.style().polish(row)
            except RuntimeError:        # the window was closed meanwhile
                pass
        flash(True)
        QTimer.singleShot(2200, lambda: flash(False))

    def _switch(self, on: bool) -> Switch:
        sw = Switch(on)
        self.switches.append(sw)
        return sw

    # -- pages ----------------------------------------------------------------------
    def _build_general(self) -> None:
        p = self._page("general", tr("General"), "set-general")
        p.section(tr("Language"))
        self.lang = _menu(list(LANGUAGES.items()), self.prefs.language)
        self.lang.currentIndexChanged.connect(self._lang_changed)
        p.row(Row(tr("Interface language"), tr("Menus, windows, dialogs and OpenDTA's own messages."), self.lang,
                  "idioma language english português"))
        self.translate = self._switch(self.prefs.messages == "translated")
        self.translate.toggled.connect(self._messages_changed)
        p.row(Row(tr("Translate Stata-style messages"),
                  tr("Error messages in Portuguese. Keep them in English to compare with Stata logs "
                     "and manuals."), self.translate, "errors erros mensagens"))
        p.section(tr("Errors"))
        self.hints = self._switch(self._session_setting("hints", self.prefs.hints))
        self.hints.toggled.connect(lambda on: self._set_cmd("hints", "on" if on else "off"))
        p.row(Row(tr("Explain errors"),
                  tr("After an error message, explain what went wrong and suggest a fix. Same as "
                     "set hints on."), self.hints, "hints dicas explicações"))
        p.section(tr("Session"))
        self.startdir_label = QLabel()
        self.startdir_label.setObjectName("Muted")
        choose = _button(tr("Choose..."))
        choose.clicked.connect(self._choose_start_dir)
        default = _button(tr("Default"))
        default.clicked.connect(self._clear_start_dir)
        box = QWidget()
        bl = QHBoxLayout(box)
        bl.setContentsMargins(0, 0, 0, 0)
        bl.setSpacing(8)
        for w in (self.startdir_label, choose, default):
            bl.addWidget(w)
        p.row(Row(tr("Initial folder"),
                  tr("The working directory OpenDTA starts in. Takes effect the next time OpenDTA opens."),
                  box, "pasta folder diretório directory cd início start"))
        self._refresh_start_dir()
        self.autolog = self._switch(self.prefs.auto_log)
        self.autolog.toggled.connect(self._autolog_changed)
        p.row(Row(tr("Save a log of each session"),
                  tr("Starts a text log when OpenDTA opens, in a logs folder inside the initial folder (or "
                     "your home folder). Takes effect the next time OpenDTA opens."), self.autolog,
                  "log registro histórico"))
        self.review_do = self._switch(self.prefs.review_do_commands)
        self.review_do.toggled.connect(self._review_do_changed)
        p.row(Row(tr("Show the commands of do-files in Review"),
                  tr("After do, run or the Do button of the Do-file Editor, Review lists the commands the "
                     "do-files ran, under the do line (click the arrow to open them; commands with errors in "
                     "red). Off: one line per do, as in Stata."), self.review_do,
                  "review histórico history do-file comandos commands"))
        p.section(tr("Safety"))
        self.confirms: dict[str, Switch] = {}
        for command, title, desc in (
                ("clear", "Confirm clear", "In the Command window, ask before clear when the data in memory "
                                           "have unsaved changes."),
                ("drop", "Confirm drop", "In the Command window, ask before drop removes many observations or "
                                         "most of the variables."),
                ("keep", "Confirm keep", "In the Command window, ask before keep removes many observations or "
                                         "most of the variables."),
                ("save", "Confirm save, replace", "In the Command window, ask before save, replace overwrites a "
                                                  "file much larger than the data in memory.")):
            sw = self._switch(self.prefs.confirm_enabled(command))
            sw.toggled.connect(lambda on, c=command: self._confirm_changed(c, on))
            self.confirms[command] = sw
            row = p.row(Row(tr(title), tr(desc), sw, f"confirmar confirmação segurança safety {command}"))
            self.anchors[f"confirm_{command}"] = ("general", row)
        p.finish()
        self.translate.setEnabled(self.prefs.language == "pt")

    def _build_appearance(self) -> None:
        p = self._page("appearance", tr("Appearance"), "set-appearance")
        p.section(tr("Appearance"))
        self.theme = Segmented([("system", "", icon("theme-system")), ("light", "", icon("theme-light")),
                                ("dark", "", icon("theme-dark"))], self.prefs.theme)
        for key, tip in (("system", "Follow the system"), ("light", "Light"), ("dark", "Dark")):
            self.theme.buttons[key].setToolTip(tr(tip))
        self.theme.changed.connect(self._theme_changed)
        p.row(Row(tr("Theme"), tr("Light, dark, or the same as the operating system."), self.theme,
                  "tema dark light escuro claro sistema"))
        self.appicon = Segmented([("light", tr("Light"), app_icon("light")),
                                  ("dark", tr("Dark"), app_icon("dark"))], self.prefs.app_icon)
        self.appicon.changed.connect(self._icon_changed)
        p.row(Row(tr("Application icon"), tr("The icon in the Dock or taskbar and in the windows."),
                  self.appicon, "ícone icon"))
        p.section(tr("Tables"))
        self.solid_lines = self._switch(self.prefs.solid_table_lines)
        self.solid_lines.toggled.connect(self._solid_lines_changed)
        p.row(Row(tr("Solid table lines"),
                  tr("Draws the lines of the tables in the Results window as continuous lines. Off: dashes, "
                     "as in Stata and in the logs."), self.solid_lines,
                  "tabelas tables linhas lines contínuas solid tracejadas dashed"))
        p.section(tr("Font"))
        self.family = QFontComboBox()
        self.family.setObjectName("Menu")
        self.family.setFontFilters(QFontComboBox.FontFilter.MonospacedFonts)
        self.family.setCurrentFont(QFont(self.prefs.font_family))
        self.family.currentFontChanged.connect(self._font_changed)
        fam, size = default_monospace()
        p.row(Row(tr("Font family"),
                  tr("Results, Command and Review windows.") + " " + L(f"Default on this system: {fam}.",
                                                                     f"Padrão neste sistema: {fam}."),
                  self.family, "fonte font monospace"))
        self.size = _menu([(str(k), f"{k} pt") for k in range(8, 25)], str(self.prefs.font_size))
        self.size.setMinimumWidth(100)
        self.size.currentIndexChanged.connect(self._font_changed)
        p.row(Row(tr("Font size"), L(f"Default: {size} pt.", f"Padrão: {size} pt."), self.size,
                  "tamanho size"))
        from .panes import ResultsView          # the same widget as the Results window, so the lines show
        self.preview = ResultsView()
        self.preview.append_styled(SAMPLE.rstrip("\n"), "text")
        self.preview.setFixedHeight(150)
        self.preview.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        p.widget_(self._spacer(10))
        p.widget_(self.preview)
        p.finish()

    def _build_results(self) -> None:
        p = self._page("results", tr("Results"), "set-results")
        p.section(tr("Numbers"))
        self.superscript = self._switch(self._session_setting("superscript", self.prefs.superscript))
        self.superscript.toggled.connect(lambda on: self._set_cmd("superscript", "on" if on else "off"))
        p.row(Row(tr("Superscript exponents"),
                  tr("Show scientific notation as 1.0000×10¹⁰ instead of 1.00000e+10. Text that can be "
                     "read back (exports, macros) always uses Stata's form."), self.superscript,
                  "notação científica scientific notation"))
        self.numerics = Segmented([("stata", tr("Stata")), ("precise", tr("Precise"))], self.prefs.numerics)
        self.numerics.changed.connect(lambda v: self._set_cmd("numerics", v))
        p.row(Row(tr("Arithmetic"),
                  tr("Stata reproduces Stata's rounding. Precise gives the same estimates with more exact "
                     "arithmetic (set numerics)."), self.numerics, "numerics precisão precision"))
        p.section(tr("Color scheme"))
        self.rpalette = _menu(palette_labels(), self.prefs.results_palette)
        self.rpalette.currentIndexChanged.connect(self._rpalette_changed)
        p.row(Row(tr("Color scheme"),
                  tr("Colors of the Results window and the Review list: the classic look, high contrast, a "
                     "soft one, or a scheme common in code editors."), self.rpalette, _PALETTE_WORDS))
        self.rpreview = QLabel()
        self.rpreview.setObjectName("PalettePreview")
        self.rpreview.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        self.rpreview.setContentsMargins(16, 12, 16, 12)
        p.widget_(self._spacer(10))
        p.widget_(self.rpreview)
        p.finish()

    def _build_editor(self) -> None:
        p = self._page("editor", tr("Do-file Editor"), "set-editor")
        p.section(tr("Color scheme"))
        self.epalette = _menu(palette_labels(), self.prefs.editor_palette)
        self.epalette.currentIndexChanged.connect(self._epalette_changed)
        p.row(Row(tr("Color scheme"),
                  tr("Syntax colors and background of the Do-file Editor, independent of the Results window."),
                  self.epalette, _PALETTE_WORDS))
        self.epreview = QLabel()
        self.epreview.setObjectName("PalettePreview")
        self.epreview.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        self.epreview.setContentsMargins(16, 12, 16, 12)
        p.widget_(self._spacer(10))
        p.widget_(self.epreview)
        p.section(tr("Session"))
        self.reopen = self._switch(self.prefs.reopen_files)
        self.reopen.toggled.connect(self._reopen_changed)
        p.row(Row(tr("Reopen last session's files"),
                  tr("When the Do-file Editor opens, reopen the files that were open when it was last used."),
                  self.reopen, "reabrir arquivos reopen files sessão session"))
        p.section(tr("Drag and drop"))
        self.run_dropped = self._switch(self.prefs.run_dropped_dofiles)
        self.run_dropped.toggled.connect(self._run_dropped_changed)
        p.row(Row(tr("Run do-files dropped on the window"),
                  tr("Off: a do-file dragged onto OpenDTA opens in the Do-file Editor without running."),
                  self.run_dropped, "arrastar soltar drag drop do-file executar run abrir open"))
        p.finish()

    def _build_updates(self) -> None:
        p = self._page("updates", tr("Updates"), "set-updates")
        p.section(tr("Updates"))
        self.auto = self._switch(self.prefs.check_updates)
        self.auto.toggled.connect(self._auto_changed)
        p.row(Row(tr("Check automatically"),
                  tr("Look for a new version on GitHub when OpenDTA opens, at most once a day."), self.auto,
                  "atualização update github"))
        self.check_btn = _button(tr("Check Now"))
        self.check_btn.clicked.connect(self._check_now)
        self.update_info = Row(tr("Version"), "", self.check_btn, "versão version")
        p.row(self.update_info)
        self._refresh_update_info()
        p.finish()

    def _build_about(self) -> None:
        p = self._page("about", tr("About"), "set-about")
        head = QWidget()
        hl = QHBoxLayout(head)
        hl.setContentsMargins(0, 10, 0, 6)
        hl.setSpacing(16)
        logo = QLabel()
        logo.setPixmap(app_icon(self.prefs.app_icon).pixmap(64, 64))
        self.logo = logo
        hl.addWidget(logo)
        tv = QVBoxLayout()
        t = QLabel(f"OpenDTA {__version__}")
        t.setObjectName("BigTitle")
        tv.addWidget(t)
        d = QLabel(tr("Open-source interpreter for do-files and .dta data, compatible with Stata 14."))
        d.setObjectName("Muted")
        d.setWordWrap(True)
        tv.addWidget(d)
        hl.addLayout(tv, 1)
        p.widget_(head)
        p.section(tr("Project"))
        for label, url, desc in (
                ("GitHub", "https://github.com/yaskepalhares/opendta", tr("Source code, issues and releases.")),
                (tr("Report a problem"), "https://github.com/yaskepalhares/opendta/issues/new/choose",
                 tr("Differences from Stata and bugs.")),
                (tr("Release notes"), "https://github.com/yaskepalhares/opendta/releases", tr("What changed in each version."))):
            b = _button(tr("Open"))
            b.clicked.connect(lambda _c=False, u=url: self._open(u))
            p.row(Row(label, desc, b, url))
        p.section(tr("License"))
        p.row(Row("MIT", tr("Free to use, copy, modify and distribute. Stata is a registered trademark of "
                            "StataCorp LLC; OpenDTA is not affiliated with StataCorp."), None, "licença license"))
        p.section(tr("Reset"))
        reset = _button(tr("Restore Defaults"), "Danger")
        reset.clicked.connect(self._reset)
        p.row(Row(tr("Restore default settings"),
                  tr("Theme, icon, font, language and the other settings on these pages."), reset,
                  "redefinir reset padrão"))
        p.finish()

    @staticmethod
    def _spacer(h: int) -> QWidget:
        w = QWidget()
        w.setFixedHeight(h)
        return w

    # -- applying changes -----------------------------------------------------------
    def _session_setting(self, name: str, default: bool) -> bool:
        s = getattr(self.win, "session", None)
        if s is None:
            return default
        return s.settings.get(name, "on" if default else "off") == "on"

    def _apply(self) -> None:
        self.prefs.sync()
        self.applied.emit()
        self._restyle()

    def _set_cmd(self, name: str, value: str) -> None:
        """hints, superscript and numerics go through `set ..., permanently`, like the menu items."""
        s = getattr(self.win, "session", None)
        if s is not None:
            s.run(f"quietly set {name} {value}, permanently")
            if hasattr(self.win, "refresh_state"):
                self.win.refresh_state()
        else:
            setattr(self.prefs, name, value if name == "numerics" else value == "on")
        self.prefs.sync()

    def _lang_changed(self) -> None:
        self.prefs.language = self.lang.currentData() or "en"
        if self.prefs.language != "pt":
            self.prefs.messages = "stata"
            self.translate.setChecked(False)
        self.translate.setEnabled(self.prefs.language == "pt")
        self._apply()
        # the window's own texts change language: open it again on the same page
        if hasattr(self.win, "show_preferences"):
            QTimer.singleShot(0, lambda: (self.close(), self.win.show_preferences("general")))

    def _messages_changed(self, on: bool) -> None:
        self.prefs.messages = "translated" if on and self.prefs.language == "pt" else "stata"
        self._apply()

    def _theme_changed(self, key: str) -> None:
        self.prefs.theme = key
        self._apply()

    def _icon_changed(self, key: str) -> None:
        self.prefs.app_icon = key
        self.logo.setPixmap(app_icon(key).pixmap(64, 64))
        self._apply()

    def _solid_lines_changed(self, on: bool) -> None:
        self.prefs.solid_table_lines = on
        self._apply()

    def _font_changed(self, *_a) -> None:
        self.prefs.font_family = self.family.currentFont().family()
        try:
            self.prefs.font_size = int(self.size.currentData())
        except (TypeError, ValueError):
            pass
        self._apply()

    def _refresh_start_dir(self) -> None:
        path = self.prefs.start_dir
        shown = QFontMetrics(self.startdir_label.font()).elidedText(path, Qt.TextElideMode.ElideMiddle, 220) \
            if path else tr("Not set")
        self.startdir_label.setText(shown)
        self.startdir_label.setToolTip(path)

    def _choose_start_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, tr("Initial folder"), self.prefs.start_dir or os.getcwd())
        if path:
            self.prefs.start_dir = path
            self.prefs.sync()
            self._refresh_start_dir()

    def _clear_start_dir(self) -> None:
        self.prefs.start_dir = ""
        self.prefs.sync()
        self._refresh_start_dir()

    def _autolog_changed(self, on: bool) -> None:
        self.prefs.auto_log = on
        self.prefs.sync()

    def _review_do_changed(self, on: bool) -> None:
        self.prefs.review_do_commands = on
        self.prefs.sync()

    def _confirm_changed(self, command: str, on: bool) -> None:
        self.prefs.set_confirm(command, on)
        self.prefs.sync()

    def _run_dropped_changed(self, on: bool) -> None:
        self.prefs.run_dropped_dofiles = on
        self.prefs.sync()

    def _reopen_changed(self, on: bool) -> None:
        self.prefs.reopen_files = on
        self.prefs.sync()

    def _rpalette_changed(self, *_a) -> None:
        self.prefs.results_palette = self.rpalette.currentData() or "classic"
        self._apply()

    def _epalette_changed(self, *_a) -> None:
        self.prefs.editor_palette = self.epalette.currentData() or "classic"
        self._apply()

    def _auto_changed(self, on: bool) -> None:
        self.prefs.check_updates = on
        self.prefs.sync()

    def _check_now(self) -> None:
        ctl = getattr(self.win, "updates", None)
        if ctl is not None:
            ctl.check_now()
            QTimer.singleShot(4000, self._refresh_update_info)

    def _refresh_update_info(self, *_a) -> None:
        last = self.prefs.last_update_check
        when = time.strftime("%d/%m/%Y %H:%M", time.localtime(last)) if last else L("never", "nunca")
        desc = L(f"You have OpenDTA {__version__}. Last check: {when}.",
                 f"Você tem o OpenDTA {__version__}. Última verificação: {when}.")
        for lab in self.update_info.findChildren(QLabel):
            if lab.objectName() == "RowDesc":
                lab.setText(desc)
                break
        else:
            d = QLabel(desc)
            d.setObjectName("RowDesc")
            d.setWordWrap(True)
            self.update_info.layout().itemAt(0).layout().addWidget(d)

    def _open(self, url: str) -> None:
        from ..updater import open_externally
        open_externally(url)

    def _reset(self) -> None:
        ans = QMessageBox.question(self, tr("Restore Defaults"),
                                   tr("Restore all settings to their defaults?"))
        if ans != QMessageBox.StandardButton.Yes:
            return
        keep_layout = self.prefs.layout()
        self.prefs.reset()
        if all(v is not None for v in keep_layout):
            self.prefs.save_layout(*keep_layout)       # the window layout is not a setting
        for name, value in (("hints", "on"), ("superscript", "on"), ("numerics", "stata")):
            self._set_cmd(name, value)
        self._apply()
        # rebuild the controls with the default values
        self.lang.setCurrentIndex(max(0, self.lang.findData(self.prefs.language)))
        self.theme.set_current(self.prefs.theme)
        self.appicon.set_current(self.prefs.app_icon)
        self.family.setCurrentFont(QFont(self.prefs.font_family))
        self.size.setCurrentIndex(max(0, self.size.findData(str(self.prefs.font_size))))
        self.numerics.set_current("stata")
        for combo, key in ((self.rpalette, self.prefs.results_palette), (self.epalette, self.prefs.editor_palette)):
            combo.blockSignals(True)
            combo.setCurrentIndex(max(0, combo.findData(key)))
            combo.blockSignals(False)
        self._refresh_start_dir()
        for sw, on in ((self.hints, True), (self.superscript, True), (self.auto, True), (self.translate, False),
                       (self.solid_lines, True), (self.run_dropped, False),
                       (self.reopen, False), (self.autolog, False),
                       *((sw, True) for sw in self.confirms.values())):
            sw.blockSignals(True)
            sw.setChecked(on)
            sw._set(1.0 if on else 0.0)
            sw.blockSignals(False)

    # -- search and colors --------------------------------------------------------------
    def _filter(self, text: str) -> None:
        q = text.strip().lower()
        first = None
        for _k, page, it in self.pages:
            n = page.filter(q)
            it.setHidden(bool(q) and n == 0)
            if n and first is None:
                first = it
        for i in range(self.nav.count()):
            it = self.nav.item(i)
            if it.data(Qt.ItemDataRole.UserRole) == "caption":
                it.setHidden(bool(q))
        cur = self.nav.currentItem()
        if q and first is not None and (cur is None or cur.isHidden()):
            self.nav.setCurrentItem(first)

    def _restyle(self) -> None:
        dark = is_dark()
        c = COLORS["dark" if dark else "light"]
        self.setStyleSheet(stylesheet(c))
        for sw in self.switches:
            sw.colors = c
            sw.update()
        for _k, _p, it in self.pages:
            it.setIcon(icon(it.data(Qt.ItemDataRole.UserRole + 1)))
        self._search_action.setIcon(icon("search"))
        for key in ("system", "light", "dark"):
            self.theme.buttons[key].setIcon(icon(f"theme-{key}"))
        # preview in the Results window's colors and font
        from .panes import monospace_font
        theme = "dark" if dark else "light"
        rs = results_scheme(theme, self.prefs.results_palette)
        self.preview.setFont(monospace_font(self.prefs))
        self.preview.set_scheme(rs)
        self.preview.set_solid_rules(self.prefs.solid_table_lines)
        self.preview.setStyleSheet(f"QPlainTextEdit#Results {{ background: {rs.background}; "
                                   f"color: {rs.colors['text']}; border: 1px solid {c['border']}; "
                                   f"border-radius: 10px; padding: 6px 10px; }}")
        ec = editor_colors(self.prefs.editor_palette, theme)
        efg = ec.get("fg", c["text"])
        for lab, bg, fg, body in ((self.rpreview, rs.background, rs.colors["text"], _results_html(rs)),
                                  (self.epreview, ec.get("bg", c["content"]), efg, _editor_html(ec, efg))):
            lab.setFont(monospace_font(self.prefs))
            lab.setTextFormat(Qt.TextFormat.RichText)
            lab.setText(body)
            lab.setStyleSheet(f"QLabel#PalettePreview {{ background: {bg}; color: {fg}; "
                              f"border: 1px solid {c['border']}; border-radius: 10px; }}")

    def keyPressEvent(self, e) -> None:  # noqa: N802
        if e.key() == Qt.Key.Key_Escape:
            self.close()
            return
        super().keyPressEvent(e)
