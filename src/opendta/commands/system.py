"""set, version, clear, pwd, cd, about, return list, ereturn list, creturn list."""

from __future__ import annotations

import os
import platform
import sys
from pathlib import Path
from typing import TYPE_CHECKING

from .. import __version__
from ..core.errors import StataError
from ..core.formats import format_value
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
    if name in ("more", "rmsg", "varabbrev", "trace") and value.split(",")[0].strip() not in ("on", "off"):
        raise StataError(198, "invalid syntax")
    s.settings[name] = value.split(",")[0].strip()


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


def _results_list(s: "Session", store: dict, prefix: str) -> None:
    out = s.output
    scalars = {k: v for k, v in store.items() if not isinstance(v, str)}
    macros = {k: v for k, v in store.items() if isinstance(v, str)}
    if scalars:
        out.write("\nscalars:\n", "text")
        for k, v in scalars.items():
            out.write(f"{prefix + '(' + k + ')':>20} =  ", "text")
            out.write(format_value(v, "%9.0g", pad=False, sign_outside_width=True) + "\n", "result")
    if macros:
        out.write("\nmacros:\n", "text")
        for k, v in macros.items():
            out.write(f"{prefix + '(' + k + ')':>20} : ", "text")
            out.write(f'"{v}"\n', "result")


@command("return")
def cmd_return(s: "Session", args: str) -> None:
    sub = args.strip()
    if sub in ("list", "li", "l"):
        _results_list(s, s.r, "r")
        return
    if sub in ("clear",):
        s.r.clear()
        return
    raise StataError(198, "return scalar/local chegam com program define (fase 2)")


@command("ereturn")
def cmd_ereturn(s: "Session", args: str) -> None:
    sub = args.strip()
    if sub in ("list", "li", "l"):
        _results_list(s, s.e, "e")
        return
    if sub == "clear":
        s.e.clear()
        return
    raise StataError(198, "invalid syntax")