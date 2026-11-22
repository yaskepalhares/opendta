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
    if name.startswith("httpproxy"):
        # these are always remembered (compat 2008 shows it for httpproxyhost; VERIFY the
        # other httpproxy settings)
        s.output.write(f"(set {name} preference recorded)\n", "text")
    opts = value.partition(",")[2].strip()
    if opts and "permanently".startswith(opts.split()[0]) and len(opts.split()[0]) >= 4:
        hook = s.ui_hooks.get("set_permanently")
        if hook is not None:
            hook(name, s.settings[name])


def _set_rng(s: "Session", name: str, value: str) -> None:
    from ..core import rng
    v = value.strip()
    if name == "rng":
        if v not in ("default", "mt64", "kiss32"):
            raise StataError(198, f"rng {v} not available")   # VERIFY message
        s.settings["rng"] = v
        rng.use("kiss32" if v == "kiss32" else "mt64")
        return
    if v.startswith(("X", "K")):
        try:
            rng._GENERATORS["kiss32" if v.startswith("K") else "mt64"].set_state(v)
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


# ---------------------------------------------------------------------------
# query and set_defaults ([R] query, [R] set_defaults)
# ---------------------------------------------------------------------------

_CATEGORIES = {"memory": "Memory settings", "output": "Output settings",
               "interface": "Interface settings", "graphics": "Graphics settings",
               "efficiency": "Efficiency settings", "network": "Network settings",
               "update": "Update settings", "trace": "Trace (programming debugging) settings",
               "mata": "Mata settings", "unicode": "Unicode", "other": "Other settings"}

# What query shows, category by category, as Stata 14 lists it (compat 2005 and 2008): the
# setting and its note; "-" is a separator line inside the category.
_P, _E3 = "may be empty, off, or on", "may be empty or a numerical format"
_QUERY = {
    "memory": [("maxvar", "2048-32767; max. vars allowed"),
               ("matsize", "10-11000; max. # vars in models"), ("niceness", "0-10"),
               ("min_memory", "0-1600g"), ("max_memory", "32m-1600g or ."),
               ("segmentsize", "1m-32g")],
    "output": [("more", ""), ("rmsg", ""), ("dp", "may be period or comma"),
               ("linesize", "characters"), ("pagesize", "lines"), "-",
               ("level", "percent confidence intervals"),
               ("clevel", "percent credible intervals"), "-",
               ("showbaselevels", "may be empty, off, on, or all"), ("showemptycells", _P),
               ("showomitted", _P), ("fvlabel", ""), ("fvwrap", ""),
               ("fvwrapon", "may be word or width"), "-", ("lstretch", _P), "-",
               ("cformat", _E3), ("pformat", _E3), ("sformat", _E3), "-",
               ("coeftabresults", ""), "-", ("logtype", "may be smcl or text"), "-",
               ("eolchar", "may be mac or unix"), ("notifyuser", ""), ("playsnd", ""),
               ("include_bitmap", "")],
    "interface": [("revkeyboard", ""), ("varkeyboard", ""), ("smoothfonts", ""), "-",
                  ("linegap", "pixels"), ("scrollbufsize", "characters"),
                  ("varlabelpos", "(not relevant)"), ("reventries", "lines"), "-",
                  ("maxdb", "dialog boxes")],
    "graphics": [("graphics", ""), ("scheme", ""),
                 ("printcolor", "may be automatic, asis, gs1, gs2, gs3"),
                 ("copycolor", "may be automatic, asis, gs1, gs2, gs3")],
    "efficiency": [("adosize", "kilobytes")],
    "network": [("checksum", ""), ("timeout1", "seconds"), ("timeout2", "seconds"), "-",
                ("httpproxy", ""), ("httpproxyhost", ""), ("httpproxyport", ""), "-",
                ("httpproxyauth", ""), ("httpproxyuser", ""), ("httpproxypw", "")],
    "update": [("update_query", ""), ("update_interval", ""), ("update_prompt", "")],
    "trace": [("trace", ""), ("tracedepth", ""), ("traceexpand", ""), ("tracesep", ""),
              ("traceindent", ""), ("tracenumber", ""), ("tracehilite", "")],
    "mata": [("matastrict", ""), ("matalnum", ""), ("mataoptimize", ""),
             ("matafavor", "may be space or speed"), ("matacache", "kilobytes"),
             ("matalibs", ""), ("matamofirst", "")],
    "unicode": [("locale_ui", ""), ("locale_functions", "")],
    "other": [("type", "may be float or double"),
              ("maxiter", "max iterations for estimation commands"),
              ("searchdefault", "may be local, net, or all"),
              ("rng", "may be default, mt64, or kiss32"), ("rngstate", ""), ("varabbrev", ""),
              ("emptycells", "may be keep or drop"), ("processors", ""),
              ("odbcmgr", "may be iodbc or unixodbc"), ("odbcdriver", "may be ansi or utf16")],
}
# settings query shows that have no c() value
_NOT_IN_C = {"odbcdriver": "ansi"}

# What set_defaults output resets, in order (compat 2008): linesize and pagesize are left as
# they are, and the settings after logtype are not touched. VERIFY: the other categories, which
# here reset every setting query lists that has a known default.
_RESET = {"output": ["more", "rmsg", "dp", "linesize", "pagesize", "level", "clevel",
                     "showbaselevels", "showemptycells", "showomitted", "fvlabel", "fvwrap",
                     "fvwrapon", "lstretch", "cformat", "pformat", "sformat",
                     "coeftabresults", "logtype"]}
_DEFAULTS = {"more": "on", "rmsg": "off", "dp": "period", "level": "95", "clevel": "95",
             "showbaselevels": "", "showemptycells": "", "showomitted": "", "fvlabel": "on",
             "fvwrap": "1", "fvwrapon": "word", "lstretch": "", "cformat": "", "pformat": "",
             "sformat": "", "coeftabresults": "on", "logtype": "smcl"}


def _settable(category: str) -> list[str]:
    return [r[0] for r in _QUERY[category] if isinstance(r, tuple)]


def _bytes(x: float) -> str:
    for unit, k in (("g", 1024 ** 3), ("m", 1024 ** 2), ("k", 1024)):
        if x and x % k == 0:
            return f"{int(x // k)}{unit}"
    return f"{int(x)}"


def _shown(s: "Session", name: str) -> str:
    from ..core import missing as M
    if name == "varlabelpos":
        return ""
    if name in _NOT_IN_C:
        return s.settings.get(name, _NOT_IN_C[name])
    v = s.creturn(name)
    if isinstance(v, str):
        # VERIFY: the cut of long values (the log shows only rngstate cut, at 40 characters)
        return v[:40] + "..." if name == "rngstate" and len(v) > 40 else v
    if v >= M.SYSMISS:
        return "."
    if name in ("segmentsize", "min_memory", "max_memory"):
        return _bytes(v)
    return f"{v:g}"


@command("query", "q")
def cmd_query(s: "Session", args: str) -> None:
    """query [category], in the layout of Stata 14 (compat 2005 and 2008). VERIFY: query with
    no category, shown here as every category in turn."""
    word = args.strip()
    cats = list(_CATEGORIES)
    if word:
        match = [c for c in cats if c.startswith(word) and len(word) >= 3] or \
            ([word] if word in cats else [])
        if not match:
            raise StataError(198, f"{word} invalid")   # VERIFY
        cats = match[:1]
    out = s.output
    ls = int(s.settings.get("linesize", "80"))
    for cat in cats:
        out.write("-" * ls + "\n", "text")
        out.write(f"    {_CATEGORIES[cat]}\n", "text")
        for row in _QUERY[cat]:
            if row == "-":
                out.write("      " + "-" * (ls - 6) + "\n", "text")
                continue
            name, note = row
            line = f"      set {name:<17}{_shown(s, name):<11}{note}"
            out.write(line.rstrip() + "\n", "text")


@command("set_defaults")
def cmd_set_defaults(s: "Session", args: str) -> None:
    """set_defaults category|_all [, permanently]: back to the default values, listing each
    setting as Stata 14 does (compat 2008)."""
    body, _, opts = args.partition(",")
    word = body.strip()
    if opts.strip() and not "permanently".startswith(opts.strip()):
        raise StataError(198, f"option {opts.strip()} not allowed")
    names = [c for c in _CATEGORIES if c != "graphics"] + ["_all"]
    expected = ", ".join(names[:-1]) + ", or _all expected."
    if word not in names:
        if not word:
            raise StataError(100, "set category required\n    Nothing was found where " +
                             expected)   # compat 2005
        # the stray quote is in Stata's message (compat 2008)
        raise StataError(198, f'    {word} found where" {expected}')
    out = s.output
    for cat in (names[:-1] if word == "_all" else [word]):
        for name in _RESET.get(cat, _settable(cat)):
            if name in ("linesize", "pagesize"):
                out.write(f"-> ({name} left at current value)\n", "text")
            elif name in _DEFAULTS:
                value = _DEFAULTS[name]
                out.write(f"-> set {name} {value}".rstrip() + "\n", "text")
                cmd_set(s, f"{name} {value}")
    out.write("(preferences reset)\n", "text")


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
            from ..lang.classes import clear_classes
            clear_classes(s)
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


@command("memory", "mem")
def cmd_memory(s: "Session", args: str) -> None:
    """memory: the layout of Stata 14's report. VERIFY: OpenDTA does not preallocate, so allocated equals
    used, and what it does not track (ado-files, Mata, stored results, overhead) is 0; compat/do/2001."""
    if args.strip():
        raise StataError(198, "invalid syntax")
    ds = s.data
    mats = getattr(s, "matrices", {}) or {}
    msize = sum(int(getattr(m, "size", 0)) * 8 for m in mats.values())
    names = sum(len(v.name) + len(v.fmt) + len(v.label) + len(v.value_label) for v in ds.vars)
    data = ds.nobs * ds.width()
    rows = {"data": data, "strLs": 0, "names": names, "overhead": 0, "matrices": msize}
    total = sum(rows.values())
    w = s.output.write

    def line(label: str, n: int) -> None:
        w(f"    {label:<30}{n:>14,}{n:>25,}\n", "text")

    rule = "    " + "-" * 69 + "\n"
    w("\n  Memory usage\n", "text")
    w(f"{'used':>48}{'allocated':>25}\n", "text")
    w(rule, "text")
    line("data", data)
    line("strLs", 0)
    w(rule, "text")
    line("data & strLs", data)
    w("\n" + rule, "text")
    line("data & strLs", data)
    line("var. names, %fmts, ...", names)
    line("overhead", 0)
    w("\n", "text")
    line("Stata matrices", msize)
    line("ado-files", 0)
    line("stored results", 0)
    w("\n", "text")
    line("Mata matrices", 0)
    line("Mata functions", 0)
    w("\n", "text")
    line("set maxvar usage", 0)
    w("\n", "text")
    line("other", 0)
    w(rule, "text")
    line("grand total", total)


@command("about")
def cmd_about(s: "Session", args: str) -> None:
    out = s.output
    out.write(f"OpenDTA {__version__}\n", "result")
    out.write("Do-file interpreter compatible with Stata 14 syntax\n", "text")
    out.write(f"Python {sys.version.split()[0]} on {platform.system()}\n", "text")


