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
    if name in ("seed", "rng"):
        _set_rng(s, name, value)
        return
    if name in ("more", "rmsg", "varabbrev", "trace", "hints", "superscript", "graphics") and value.split(",")[0].strip() not in ("on", "off"):
        raise StataError(198, "invalid syntax")
    if name == "dp":
        v = value.split(",")[0].strip()
        if v not in ("comma", "period"):
            raise StataError(198, "invalid syntax")
        from ..core.formats import set_decimal_comma
        set_decimal_comma(v == "comma")
    if name == "superscript":
        from ..core.formats import set_superscript
        set_superscript(value.split(",")[0].strip() == "on")
    if name == "numerics":
        v = value.split(",")[0].strip()
        if v not in ("stata", "precise"):
            raise StataError(198, "invalid syntax")
        from ..estimation.numerics import set_precise
        set_precise(v == "precise")
    if name in ("language", "messages"):
        # OpenDTA extension (not in Stata): language of the interface and of the error
        # explanations, and optional translation of Stata's messages
        from .. import i18n
        v = value.split(",")[0].strip().lower()
        if name == "language":
            if v not in ("en", "pt", "english", "portuguese", "portugues", "português"):
                raise StataError(198, "invalid syntax")
            i18n.set_language(v)
            value = i18n.language() + value[len(value.split(",")[0]):]
        else:
            if v not in ("stata", "translated"):
                raise StataError(198, "invalid syntax")
            i18n.set_messages(v)
        hook = s.ui_hooks.get("language")
        if hook is not None:
            hook()
    s.settings[name] = value.split(",")[0].strip()
    opts = value.partition(",")[2].strip()
    if opts and "permanently".startswith(opts.split()[0]) and len(opts.split()[0]) >= 4:
        hook = s.ui_hooks.get("set_permanently")
        if hook is not None:
            hook(name, s.settings[name])


def _set_rng(s: "Session", name: str, value: str) -> None:
    from ..core import rng
    v = value.strip()
    if name == "rng":
        # VERIFY: kiss32 (generator of Stata 13 and earlier) not implemented
        if v not in ("default", "mt64"):
            raise StataError(198, f"rng {v} not available")
        s.settings["rng"] = v
        return
    if v.startswith("X"):
        try:
            rng.RNG.set_state(v)
        except ValueError:
            raise StataError(198, "invalid seed state")   # VERIFY
        return
    try:
        seed = int(v)
    except ValueError:
        raise StataError(198, "invalid syntax")
    # Stata 14 accepts negative seeds without an error (compat 0307); how it turns them into a
    # state is undocumented: here, 64-bit two's complement
    rng.RNG.seed(seed % 2**64)   # VERIFY sementes negativas
    s.settings["seed"] = str(seed)


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
        if t in ("all", "*"):
            from ..graphics.store import store
            store(s).clear()
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
    out.write("Do-file interpreter compatible with Stata 14 syntax\n", "text")
    out.write(f"Python {sys.version.split()[0]} on {platform.system()}\n", "text")


