"""log, cmdlog, assert, timer and creturn list.

log using writes the output as text (.log) or SMCL (.smcl, the Stata 14
default). Several logs can be open at the same time with name().
cmdlog writes only the typed commands.
"""

from __future__ import annotations

import datetime as _dt
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Callable

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import eval_vector, touse
from .registry import command
from ._util import shown_path

if TYPE_CHECKING:
    from ..session import Session

_SMCL_TAG = {"text": "{txt}", "result": "{res}", "error": "{err}", "command": "{com}",
             "input": "{inp}", "hint": "{txt}"}


def _stamp() -> str:
    now = _dt.datetime.now()
    return f"{now.day:2d} {now.strftime('%b %Y, %H:%M:%S')}"


def smcl_escape(text: str) -> str:
    return text.replace("{", "\x00").replace("}", "{c )-}").replace("\x00", "{c -(}")


@dataclass
class LogFile:
    name: str
    path: Path
    kind: str                 # text | smcl
    fh: object
    on: bool = True
    style: str = ""
    listener: Callable | None = None
    extra: dict = field(default_factory=dict)

    def write(self, text: str, style: str) -> None:
        if not self.on:
            return
        if self.kind == "text":
            self.fh.write(text)
        else:
            tag = _SMCL_TAG.get(style, "{txt}")
            if tag != self.style and text.strip("\n"):
                self.fh.write(tag)
                self.style = tag
            self.fh.write(smcl_escape(text))
        self.fh.flush()


def _logs(s: "Session") -> dict[str, LogFile]:
    if not hasattr(s, "logs"):
        s.logs = {}
    return s.logs


def _info(log: LogFile, verb: str) -> list[tuple[str, str]]:
    return [("      name:  ", log.name), ("       log:  ", str(log.path)),
            ("  log type:  ", log.kind), (f" {verb} on:  ", _stamp())]


def _show_info(s: "Session", log: LogFile, verb: str, *, rule_before: bool) -> None:
    out = s.output
    if rule_before:
        out.write("-" * 79 + "\n", "text")
    for label, value in _info(log, verb):
        out.write(label, "text")
        out.write(value + "\n", "result")
    if not rule_before:
        out.write("-" * 79 + "\n", "text")


@command("log")
def cmd_log(s: "Session", args: str) -> None:
    t = args.strip()
    sub, _, rest = t.partition(" ")
    logs = _logs(s)
    if sub == "using":
        head, comma, opts = rest.partition(",")
        o = match_options(opts, {"append": 1, "replace": 1, "text": 1, "smcl": 1,
                                 "name": 4}) if comma else {}
        name = str(o.get("name", "<unnamed>")).strip()
        # in batch mode the batch log is the unnamed log (compat 2001: r(604) before the path)
        if name in logs or (name == "<unnamed>" and getattr(s, "batch_log", False)):
            # VERIFY: Stata 14 text
            raise StataError(604, "log file already open")
        path = Path(strip_outer_quotes(head.strip())).expanduser()
        if path.suffix == "":
            path = path.with_suffix(".log" if o.get("text") else ".smcl")
        kind = "smcl" if o.get("smcl") or (path.suffix == ".smcl" and not o.get("text")) else "text"
        if path.exists() and not (o.get("replace") or o.get("append")):
            raise StataError(602, f"file {shown_path(path)} already exists")
        if o.get("replace") and not path.exists():
            # Stata shows the full path here (compat/do/2006)
            s.output.write(f"(note: file {shown_path(path.absolute())} not found)\n", "text")
        try:
            fh = open(path, "a" if o.get("append") else "w", encoding="utf-8")
        except OSError:
            raise StataError(603, f"file {shown_path(path.absolute())} could not be opened")
        log = LogFile(name, path.resolve(), kind, fh)
        if kind == "smcl" and not o.get("append"):
            fh.write("{smcl}\n")
        # header: on screen and in the file
        if kind == "smcl":
            fh.write("{txt}{sf}{ul off}{.-}\n")
            for label, value in _info(log, "opened"):
                fh.write("{txt}" + label + "{res}" + smcl_escape(value) + "\n")
            fh.write("{txt}\n")
        else:
            fh.write("-" * 79 + "\n")
            for label, value in _info(log, "opened"):
                fh.write(label + value + "\n")
        _show_info(s, log, "opened", rule_before=True)
        log.listener = log.write
        s.output.add_listener(log.listener)
        logs[name] = log
        return
    if sub in ("close", "off", "on", "query", ""):
        names = rest.split() if rest.strip() else []
        if sub == "close" and names == ["_all"]:
            names = list(logs)
        if not logs:
            if sub == "query" or sub == "":
                s.output.write("no log file open\n", "text")
                return
            raise StataError(606, "no log file open")
        if not names:
            names = ["<unnamed>"] if "<unnamed>" in logs else list(logs)[:1]
        for n in names:
            if n not in logs:
                raise StataError(111, f"log {n} not found")   # VERIFY
            log = logs[n]
            if sub == "close":
                s.output.remove_listener(log.listener)
                if log.kind == "smcl":
                    log.fh.write("{txt}{sf}{ul off}{.-}\n")
                    for label, value in _info(log, "closed"):
                        log.fh.write("{txt}" + label + "{res}" + smcl_escape(value) + "\n")
                    log.fh.write("{txt}{sf}{ul off}{.-}\n{smcl}\n{txt}{sf}{ul off}\n")
                else:
                    for label, value in _info(log, "closed"):
                        log.fh.write(label + value + "\n")
                    log.fh.write("-" * 79 + "\n")
                log.fh.close()
                del logs[n]
                for label, value in _info(log, "closed"):
                    s.output.write(label, "text")
                    s.output.write(value + "\n", "result")
            elif sub == "off":
                log.on = False
            elif sub == "on":
                log.on = True
            else:
                status = "on" if log.on else "off"
                for label, value in [("      name:  ", log.name), ("       log:  ", str(log.path)),
                                     ("  log type:  ", log.kind), ("    status:  ", status)]:
                    s.output.write(label, "text")
                    s.output.write(value + "\n", "result")
        return
    raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# cmdlog
# ---------------------------------------------------------------------------

def record_command(s: "Session", line: str) -> None:
    """Called by the session for each typed command (not for do-files)."""
    cl = getattr(s, "cmdlog", None)
    if cl is not None and cl.get("on", True):
        cl["fh"].write(line.rstrip("\n") + "\n")
        cl["fh"].flush()


@command("cmdlog")
def cmd_cmdlog(s: "Session", args: str) -> None:
    t = args.strip()
    sub, _, rest = t.partition(" ")
    cl = getattr(s, "cmdlog", None)
    if sub == "using":
        head, comma, opts = rest.partition(",")
        o = match_options(opts, {"append": 1, "replace": 1}) if comma else {}
        if cl is not None:
            raise StataError(604, "log file already open")
        path = Path(strip_outer_quotes(head.strip())).expanduser()
        if path.suffix == "":
            path = path.with_suffix(".txt")
        if path.exists() and not (o.get("replace") or o.get("append")):
            raise StataError(602, f"file {shown_path(path)} already exists")
        s.cmdlog = {"path": path.resolve(), "fh": open(path, "a" if o.get("append") else "w",
                                                       encoding="utf-8"), "on": True}
        s.output.write(f"(cmdlog {path.resolve()} opened)\n", "text")
        return
    if cl is None:
        raise StataError(606, "no cmdlog file open")   # VERIFY
    if sub == "close":
        cl["fh"].close()
        s.output.write(f"(cmdlog {cl['path']} closed)\n", "text")
        s.cmdlog = None
    elif sub in ("off", "on"):
        cl["on"] = sub == "on"
    elif sub in ("query", ""):
        s.output.write(f"cmdlog: {cl['path']} ({'on' if cl['on'] else 'off'})\n", "text")
    else:
        raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# assert
# ---------------------------------------------------------------------------

@command("assert", byable=True)
def cmd_assert(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"rc0": 3, "null": 4, "fast": 4}) if p.options.strip() else {}
    expr = (p.varlist + (" = " + p.exp if p.exp is not None else "")).strip()
    if not expr:
        raise StataError(198, "invalid syntax")
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    n = int(mask.sum())
    if n == 0 and not o.get("null"):
        # VERIFY: assert without observations
        raise StataError(2000, "no observations")
    v = eval_vector(s, expr)
    if isinstance(v, str):
        raise StataError(109, "type mismatch")
    vals = np.broadcast_to(np.asarray(v, dtype=np.float64), (s.data.nobs,))[mask]
    bad = int(((vals == 0)).sum())
    if bad:
        if not o.get("rc0"):
            s.output.write(f"{bad:,} contradiction{'s' if bad != 1 else ''} in {n:,} "
                           f"observation{'s' if n != 1 else ''}\n", "error")
            raise StataError(9, "assertion is false")
        s.output.write(f"{bad:,} contradiction{'s' if bad != 1 else ''} in {n:,} "
                       f"observation{'s' if n != 1 else ''}\n", "error")
        s.output.write("assertion is false\n", "error")


# ---------------------------------------------------------------------------
# timer
# ---------------------------------------------------------------------------

def _timers(s: "Session") -> dict[int, dict]:
    if not hasattr(s, "timers"):
        s.timers = {}
    return s.timers


@command("timer")
def cmd_timer(s: "Session", args: str) -> None:
    words = args.split()
    if not words:
        raise StataError(198, "invalid syntax")
    sub = words[0]
    ts = _timers(s)
    nums = [int(w) for w in words[1:]] if len(words) > 1 else []
    for k in nums:
        if not 1 <= k <= 100:
            raise StataError(198, "timer number must be between 1 and 100")
    if sub == "on":
        if len(nums) != 1:
            raise StataError(198, "invalid syntax")
        t = ts.setdefault(nums[0], {"total": 0.0, "count": 0, "start": None})
        if t["start"] is not None:
            raise StataError(198, f"timer {nums[0]} already on")   # VERIFY
        t["start"] = time.perf_counter()
    elif sub == "off":
        if len(nums) != 1:
            raise StataError(198, "invalid syntax")
        t = ts.get(nums[0])
        if t is None or t["start"] is None:
            raise StataError(198, f"timer {nums[0]} not on")   # VERIFY
        t["total"] += time.perf_counter() - t["start"]
        t["count"] += 1
        t["start"] = None
    elif sub == "clear":
        if nums:
            for k in nums:
                ts.pop(k, None)
        else:
            ts.clear()
    elif sub == "list":
        out = s.output
        s.r = {}
        for k in sorted(nums or ts):
            t = ts.get(k)
            if not t:
                continue
            total = t["total"]
            out.write(f"{k:4d}: {total:9.2f} / {t['count']:8d} = {total / max(t['count'], 1):10.4f}\n",
                      "text")   # VERIFY: layout
            s.r[f"t{k}"] = total
            s.r[f"nt{k}"] = float(t["count"])
    else:
        raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# creturn list
# ---------------------------------------------------------------------------

# groups, rows and annotations of Stata 14's creturn list (compat 2003)
_CRETURN_GROUPS = [
    ("System values", [
        "-", "current_date", "current_time", ("rmsg_time", "seconds, from set rmsg"), "-",
        "stata_version", ("version", "version"), ("userversion", "version"), "-", "born_date",
        "flavor", "bit", "SE", "MP", ("processors", "Stata/MP, set processors"),
        "processors_lic", "processors_mach", "processors_max", "mode", "console", "-", "os",
        "osdtl", "hostname", "machine_type", "byteorder", "username", "-",
    ]),
    ("Directories and paths", [
        "-", ("sysdir_stata", "sysdir"), ("sysdir_base", "sysdir"), ("sysdir_site", "sysdir"),
        ("sysdir_plus", "sysdir"), ("sysdir_personal", "sysdir"), ("sysdir_oldplace", "sysdir"),
        "tmpdir", "-", ("adopath", "adopath"), ("pwd", "cd"), "dirsep", "-",
    ]),
    ("System limits", [
        "-", "max_N_theory", ("max_k_theory", "set maxvar"), ("max_width_theory", "set maxvar"),
        "-", "max_matsize", "min_matsize", "-", "max_macrolen", ("macrolen", "set maxvar"),
        "charlen", "max_cmdlen", ("cmdlen", "set maxvar"), "namelenbyte", "namelenchar", "eqlen",
        "-",
    ]),
    ("Numerical and string limits", [
        "-", "mindouble", "maxdouble", "epsdouble", "smallestdouble", "-", "minfloat",
        "maxfloat", "epsfloat", "-", "minlong", "maxlong", "-", "minint", "maxint", "-",
        "minbyte", "maxbyte", "-", "maxstrvarlen", "maxstrlvarlen", "maxvlabellen", "-",
    ]),
    ("Current dataset", [
        "-", "N", "k", "width", "changed", "filename", "filedate", "-",
    ]),
    ("Memory settings", [
        "-", "memory", ("maxvar", "set maxvar"), ("matsize", "set matsize"),
        ("niceness", "set min_memory"), ("min_memory", "set min_memory"),
        ("max_memory", "set max_memory"), ("segmentsize", "set segmentsize"), "-",
    ]),
    ("Output settings", [
        "-", ("more", "set more"), ("rmsg", "set rmsg"), ("dp", "set dp"),
        ("linesize", "set linesize"), ("pagesize", "set pagesize"), ("logtype", "set logtype"),
        "noisily", ("eolchar", "set eolchar"), ("notifyuser", "set notifyuser"),
        ("playsnd", "set playsnd"), ("include_bitmap", "set include_bitmap"), "-",
        ("level", "set level"), ("clevel", "set clevel"), "-",
        ("showbaselevels", "set showbaselevels"), ("showemptycells", "set showemptycells"),
        ("showomitted", "set showomitted"), ("fvlabel", "set fvlabel"), ("fvwrap", "set fvwrap"),
        ("fvwrapon", "set fvwrapon"), ("lstretch", "set lstretch"), "-",
        ("cformat", "set cformat"), ("sformat", "set sformat"), ("pformat", "set pformat"), "-",
        ("coeftabresults", "set coeftabresults"),
    ]),
    ("Interface settings", [
        "-", ("reventries", "set reventries"), ("revkeyboard", "set revkeyboard"),
        ("varkeyboard", "set varkeyboard"), ("smoothfonts", "set smoothfonts"),
        ("linegap", "set linegap"), ("scrollbufsize", "set scrollbufsize"),
        ("maxdb", "set maxdb"), "-",
    ]),
    ("Graphics settings", [
        "-", ("graphics", "set graphics"), ("scheme", "set scheme"),
        ("printcolor", "set printcolor"), ("copycolor", "set copycolor"), "-",
    ]),
    ("Efficiency settings", [
        "-", ("adosize", "set adosize"), "-",
    ]),
    ("Network settings", [
        "-", ("checksum", "set checksum"), ("timeout1", "set timeout1"),
        ("timeout2", "set timeout2"), "-", ("httpproxy", "set httpproxy"),
        ("httpproxyhost", "set httpproxyhost"), ("httpproxyport", "set httpproxyport"), "-",
        ("httpproxyauth", "set httpproxyauth"), ("httpproxyuser", "set httpproxyuser"),
        ("httpproxypw", "set httpproxypw"), "-",
    ]),
    ("Update settings", [
        "-", ("update_query", "set update_query"), ("update_interval", "set update_interval"),
        ("update_prompt", "set update_prompt"), "-",
    ]),
    ("Trace (program debugging) settings", [
        "-", ("trace", "set trace"), ("tracedepth", "set tracedepth"),
        ("tracesep", "set tracesep"), ("traceindent", "set traceindent"),
        ("traceexpand", "set traceexpand"), ("tracenumber", "set tracenumber"),
        ("tracehilite", "set tracehilite"), "-",
    ]),
    ("Mata settings", [
        "-", ("matastrict", "set matastrict"), ("matalnum", "set matalnum"),
        ("mataoptimize", "set mataoptimize"), ("matafavor", "set matafavor"),
        ("matacache", "set matacache"), ("matalibs", "set matalibs"),
        ("matamofirst", "set matamofirst"), "-",
    ]),
    ("Unicode settings", [
        "-", ("locale_ui", "set locale_ui"), ("locale_functions", "set locale_functions"),
        ("locale_icudflt", "unicode locale"), "-",
    ]),
    ("Other settings", [
        "-", ("type", "set type"), ("maxiter", "set maxiter"),
        ("searchdefault", "set searchdefault"), ("rng", "set rng"), "rng_current",
        ("rngstate", "set rngstate"), ("varabbrev", "set varabbrev"),
        ("emptycells", "set emptycells"), ("odbcmgr", "set odbcmgr"), "-",
    ]),
    ("Other", [
        "-", "pi", "alpha", "ALPHA", "Mons", "Months", "Wdays", "Weekdays", ("rc", "capture"),
        "-",
    ]),
]


@command("creturn")
def cmd_creturn(s: "Session", args: str) -> None:
    if args.strip() not in ("list", "l", "li", "lis"):
        raise StataError(198, "invalid syntax")
    out = s.output
    rule = "    " + "-" * (int(s.settings.get("linesize", "80")) - 4)
    for title, rows in _CRETURN_GROUPS:
        out.write(f"\n{title}\n\n", "text")
        for row in rows:
            if row == "-":
                out.write(rule + "\n", "text")
                continue
            n, note = row if isinstance(row, tuple) else (row, "")
            v = s.creturn(n)
            if isinstance(v, str):
                # VERIFY: strings longer than 19 characters are cut to 17 and ".."
                shown = '"' + (v if len(v) <= 19 else v[:17] + "..") + '"'
            elif v >= M.SYSMISS:
                shown = "."
            else:
                shown = format_value(float(v), "%18.0g", pad=False).strip()
            out.write(f"{'c(' + n + ')':>23} = ", "text")
            out.write((f"{shown:<27}({note})" if note else shown) + "\n", "result")
