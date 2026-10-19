"""pause, sleep, window, plugin and viewsource ([P] pause, [P] sleep, [P] window, [P] plugin,
[P] viewsource).

The commands that need the graphical interface (window fopen, fsave, stopbox, push, manage and
menu) call a hook in `Session.ui_hooks`; without the interface they do what the manual says for a
window that does not exist, or act as if the user had cancelled.
"""

from __future__ import annotations

import inspect
import re
import sys
import time
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..lang.adopath import find_file
from ..lang.words import strip_outer_quotes
from .registry import command, lookup

if TYPE_CHECKING:
    from ..session import Session


# ---------------------------------------------------------------------------
# sleep, pause
# ---------------------------------------------------------------------------

@command("sleep")
def cmd_sleep(s: "Session", args: str) -> None:
    if not args.strip():
        raise StataError(111, " not found")        # Stata evaluates the empty argument
    ms = s.eval(args.strip())
    if isinstance(ms, str) or ms != ms or ms < 0:
        raise StataError(198, "invalid syntax")
    time.sleep(float(ms) / 1000.0)


def _read_line(s: "Session"):
    """A line typed in pause mode; None at the end of the input."""
    if s.input_provider is not None:
        text = s.input_provider()
        if text is None:
            raise StataError(1, "--Break--")
        return str(text)
    line = sys.stdin.readline()
    return None if line == "" else line.rstrip("\r\n")


@command("pause")
def cmd_pause(s: "Session", args: str) -> None:
    t = args.strip()
    if t in ("on", "off"):
        s.pause_on = t == "on"
        return
    if not getattr(s, "pause_on", False):
        return
    out = s.output
    out.write(f"pause: {t}\n" if t else "pause:\n", "text")
    while True:
        out.write("-> . ", "text")
        line = _read_line(s)
        if line is None:
            out.write("\n", "text")
            break
        out.write(line + "\n", "input")
        word = line.strip()
        if word in ("end", "q"):
            break
        if word == "BREAK":
            out.write("sending Break to calling program...\n", "text")
            raise StataError(1, "Break")
        s.run(line)
    out.write("execution resumes...\n", "text")


# ---------------------------------------------------------------------------
# window
# ---------------------------------------------------------------------------

_MANAGE = {"minimize", "restore", "maintitle", "associate", "docklabel", "prefs", "forward",
           "print", "close", "rename"}
_MENU = {"clear", "append", "refresh", "add"}


def _quoted(text: str) -> list[str]:
    """The strings of '"a" "b" c': quoted ones without the quotes, the others as they are."""
    return [m.group(1) if m.group(1) is not None else m.group(2)
            for m in re.finditer(r'"([^"]*)"|(\S+)', text)]


@command("window", "wind")
def cmd_window(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    rest = rest.strip()
    hook = s.ui_hooks.get(f"window_{sub}")
    if sub in ("fopen", "fsave"):
        words = _quoted(rest)
        if len(words) < 3 or not re.fullmatch(r"[A-Za-z_]\w*", words[0]):
            raise StataError(198, "invalid syntax")
        name, title, filt = words[:3]
        ext = words[3] if len(words) > 3 else ""
        path = hook(sub, title, filt, ext, s.macros.get_global(name)) if hook else None
        s.macros.set_global(name, path or "")
        if not path:
            # cancelled (and, without the interface, there is nobody to choose a file)
            raise StataError(601, "")
    elif sub == "stopbox":
        kind, _, text = rest.partition(" ")
        if kind not in ("stop", "note", "rusure"):
            raise StataError(198, "invalid syntax")
        lines = _quoted(text)[:4]
        if not lines:
            raise StataError(198, "invalid syntax")
        answer = hook(kind, lines) if hook else None
        if hook is None:
            for line in lines:
                s.output.write(line + "\n", "text")
            answer = kind == "note"      # VERIFY: nobody to answer without the interface
        if kind == "stop" or (kind == "rusure" and not answer):
            raise StataError(1, "")
    elif sub == "push":
        if not rest:
            raise StataError(198, "invalid syntax")
        if hook:
            hook(rest)
    elif sub == "manage":
        word = rest.split(None, 1)[0] if rest else ""
        if word not in _MANAGE:
            raise StataError(198, "invalid syntax")
        if hook:
            hook(rest)
    elif sub == "menu":
        word = rest.split(None, 1)[0] if rest else ""
        if word not in _MENU:
            raise StataError(198, "invalid syntax")
        if word == "clear":
            s.window_menu = []
        elif word in ("append", "add"):
            s.window_menu = getattr(s, "window_menu", []) + [_quoted(rest)]
        elif hook:
            hook(getattr(s, "window_menu", []))
    else:
        raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# plugin, viewsource
# ---------------------------------------------------------------------------

@command("plugin")
def cmd_plugin(s: "Session", args: str) -> None:
    sub = args.split(None, 1)[0] if args.strip() else ""
    if sub != "call":
        raise StataError(198, "invalid syntax")
    # Stata says this when the plugin was not loaded; OpenDTA never loads one
    raise StataError(199, "unrecognized command")


@command("viewsource")
def cmd_viewsource(s: "Session", args: str) -> None:
    """Shows a file found on the adopath or, for a command that OpenDTA has built in, the source of
    OpenDTA's own implementation (not Stata's)."""
    fname = strip_outer_quotes(args.strip())
    if not fname:
        raise StataError(601, 'file "" not found')
    hits = find_file(s, fname)
    if hits:
        try:
            text = hits[0].read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = hits[0].read_text(encoding="latin-1")
    else:
        spec = lookup(fname[:-4]) if fname.endswith(".ado") else None
        try:
            fn = spec.fn if spec is not None else None
            while fn is not None and "_inner" in inspect.signature(fn).parameters:
                fn = inspect.signature(fn).parameters["_inner"].default    # a wrapper of the command
            text = inspect.getsource(fn) if fn is not None else None
        except (OSError, TypeError):
            text = None
        if text is None:
            raise StataError(601, f'file "{fname}" not found')
    for line in text.splitlines():
        s.output.write(line + "\n", "text")
