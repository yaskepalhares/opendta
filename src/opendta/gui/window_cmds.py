"""The graphical side of `window` ([P] window): file dialogs, message boxes, the User menu, the
Review window and the management of the windows. commands/progutil.py parses the command and calls
these functions through `Session.ui_hooks`; without the interface it does what the manual says for
a window that does not exist.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING, Callable

from PySide6.QtWidgets import QFileDialog, QMenu, QMessageBox

from ..i18n import tr

if TYPE_CHECKING:
    from .main_window import MainWindow

# the submenus of User that the menu commands name (Stata's stUserData, ...)
USER_SUBMENUS = (("stUserData", "Data"), ("stUserGraphics", "Graphics"), ("stUserStatistics", "Statistics"))


def qt_filter(spec: str) -> str:
    """'Data (*.dta)|*.dta|All (*.*)|*.*' -> 'Data (*.dta);;All (*.*)': Stata's pairs of a
    description and its patterns (separated by semicolons) as a filter of a Qt file dialog."""
    parts = spec.split("|")
    out = []
    for k in range(0, len(parts) - 1, 2):
        desc, pats = parts[k].strip(), " ".join(p.strip() for p in parts[k + 1].split(";") if p.strip())
        out.append(desc if f"({pats})" in desc else f"{desc} ({pats})")
    if len(parts) % 2:
        last = parts[-1].strip()
        if last:
            out.append(last)
    return ";;".join(out)


def build_user_menu(menu: QMenu, items: list[list[str]], run: Callable[[str], None]) -> None:
    """Rebuilds the User menu from the `window menu append` calls, in order: each one is
    ['append', 'submenu'|'item'|'separator', parent, ...]. A parent is stUser, one of the
    USER_SUBMENUS or the title of a submenu appended before."""
    menu.clear()
    menus: dict[str, QMenu] = {"stUser": menu}
    for key, title in USER_SUBMENUS:
        menus[key] = menu.addMenu(tr(title))
    for words in items:
        if len(words) < 3 or words[0] not in ("append", "add"):
            continue
        kind, parent = words[1], menus.get(words[2])
        if parent is None:
            continue
        if kind == "submenu" and len(words) > 3:
            menus[words[3]] = parent.addMenu(words[3])
        elif kind == "item" and len(words) > 4:
            command = words[4]
            parent.addAction(words[3]).triggered.connect(lambda _=False, c=command: run(c))
        elif kind == "separator":
            parent.addSeparator()


def install(win: "MainWindow") -> None:
    """Registers the window_* hooks of `win`'s session."""

    def files(sub, title, filt, ext, default):
        start = default or os.getcwd()
        if sub == "fopen":
            path, _ = QFileDialog.getOpenFileName(win, title, start, qt_filter(filt))
        else:
            path, _ = QFileDialog.getSaveFileName(win, title, start, qt_filter(filt))
            if path and ext and not os.path.splitext(path)[1]:
                path += "." + ext.lstrip(".")
        return path

    def stopbox(kind, lines):
        text = "\n".join(lines)
        if kind == "rusure":
            buttons = QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel
            return QMessageBox.question(win, "OpenDTA", text, buttons) == QMessageBox.StandardButton.Ok
        icon = QMessageBox.critical if kind == "stop" else QMessageBox.information
        icon(win, "OpenDTA", text)
        return True

    def manage(rest: str):
        from .main_window import _present
        words = rest.split()
        what = words[0]
        arg = words[1] if len(words) > 1 else ""
        if what == "minimize":
            win.showMinimized()
        elif what == "restore":
            win.showNormal()
        elif what == "maintitle":
            title = rest.split(None, 1)[1].strip().strip('"') if len(words) > 1 else ""
            win.main_title = "" if title in ("", "reset") else title
            win.refresh_state()
        elif what == "forward":
            targets = {"command": win.command, "results": win.results, "review": win.dock_review,
                       "variables": win.dock_variables, "properties": win.dock_properties}
            if arg in targets:
                targets[arg].show()
                targets[arg].raise_()
                targets[arg].setFocus()
            elif arg == "viewer" and win.viewer is not None:
                _present(win.viewer)
            elif arg == "doeditor":
                win.show_dofile_editor("")
            elif arg == "graph":
                win.show_graph_window()
            elif arg == "data":
                win.run_command("browse")
        elif what == "close":
            if arg == "viewer" and win.viewer is not None:
                win.viewer.close()
            elif arg == "graph":
                win.session.run_command("graph close " + (words[2] if len(words) > 2 else ""))
        # prefs, print, associate, docklabel and rename change nothing in OpenDTA's windows

    win.session.ui_hooks.update(
        window_fopen=files, window_fsave=files, window_stopbox=stopbox, window_manage=manage,
        window_push=lambda text: win._add_review(text, 0),
        window_menu=lambda items: _refresh_menu(win, items))


def _refresh_menu(win: "MainWindow", items: list[list[str]]) -> None:
    win._user_items = list(items)      # kept so that a change of language rebuilds the menu
    build_user_menu(win.user_menu, win._user_items, win.run_command)
