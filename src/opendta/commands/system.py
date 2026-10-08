"""set, version, clear, pwd, cd, about, return list, ereturn list, creturn list."""

from __future__ import annotations

import os
import platform
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from .. import __version__
from ..core.errors import StataError
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


@command("set")
def cmd_set(s: "Session", args: str) -> None:
    words = args.split()
    if not words:
        raise StataError(198, "invalid syntax")
    name, value = words[0], " ".join(words[1:])
    if name == "obs":
        from .data import set_obs
        set_obs(s, value)
        s.notify_state()
        return
    if name in ("more", "rmsg", "varabbrev", "trace", "hints") and value.split(",")[0].strip() not in ("on", "off"):
        raise StataError(198, "invalid syntax")
    if name == "dp":
        v = value.split(",")[0].strip()
        if v not in ("comma", "period"):
            raise StataError(198, "invalid syntax")
        from ..core.formats import set_decimal_comma
        set_decimal_comma(v == "comma")
    s.settings[name] = value.split(",")[0].strip()
    opts = value.partition(",")[2].strip()
    if opts and "permanently".startswith(opts.split()[0]) and len(opts.split()[0]) >= 4:
        hook = s.ui_hooks.get("set_permanently")
        if hook is not None:
            hook(name, s.settings[name])


@command("version", "vers")
def cmd_version(s: "Session", args: str) -> None:
    t = args.strip()
    if not t:
        s.output.write(f"version {s.version:g}\n", "result")
        return
    head, sep, tail = t.partition(":")
    try:
        v = float(head.strip())
    except ValueError:
        raise StataError(198, "invalid syntax")
    if sep:
        saved = s.version
        s.version = v
        try:
            s.interp._execute_expanded(tail)
        finally:
            s.version = saved
    else:
        s.version = v


@command("clear")
def cmd_clear(s: "Session", args: str) -> None:
    t = args.strip()
    if t in ("", "all", "*", "programs", "results", "mata", "matrix", "ado", "frames"):
        if t in ("all", "*", "results"):
            s.r.clear()
            s.e.clear()
            s.sret.clear()
        if t in ("all", "*"):
            s.scalars.clear()
        if t in ("", "all", "*"):
            s.data.clear()
        s.notify_state()
        return
    raise StataError(198, "invalid syntax")


@command("pwd")
def cmd_pwd(s: "Session", args: str) -> None:
    s.output.write(os.getcwd() + "\n", "text")


@command("cd")
def cmd_cd(s: "Session", args: str) -> None:
    t = args.strip().strip('"')
    if not t:
        if platform.system() == "Windows":
            s.output.write(os.getcwd() + "\n", "text")
            return
        t = str(Path.home())
    try:
        os.chdir(os.path.expanduser(t))
    except OSError:
        raise StataError(170, f"unable to change to {t}")
    if platform.system() == "Windows":
        s.output.write(os.getcwd() + "\n", "text")
    s.notify_state()


@command("about")
def cmd_about(s: "Session", args: str) -> None:
    out = s.output
    out.write(f"OpenDTA {__version__}\n", "result")
    out.write("Interpretador de do-files compatível com a sintaxe do Stata 14\n", "text")
    out.write(f"Python {sys.version.split()[0]} em {platform.system()}\n", "text")


