"""log, cmdlog, assert, timer e creturn list.

log using grava a saída em texto (.log) ou SMCL (.smcl, o padrão do
Stata 14). Vários logs podem ficar abertos ao mesmo tempo com name().
cmdlog grava só os comandos digitados.
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
        if name in logs:
            # VERIFICAR: texto do Stata 14
            raise StataError(604, "log file already open")
        path = Path(strip_outer_quotes(head.strip())).expanduser()
        if path.suffix == "":
            path = path.with_suffix(".log" if o.get("text") else ".smcl")
        kind = "smcl" if o.get("smcl") or (path.suffix == ".smcl" and not o.get("text")) else "text"
        if path.exists() and not (o.get("replace") or o.get("append")):
            raise StataError(602, f"file {path} already exists")
        fh = open(path, "a" if o.get("append") else "w", encoding="utf-8")
        log = LogFile(name, path.resolve(), kind, fh)
        if kind == "smcl" and not o.get("append"):
            fh.write("{smcl}\n")
        # cabeçalho: na tela e no arquivo
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
                raise StataError(111, f"log {n} not found")   # VERIFICAR
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
    """Chamado pela sessão para cada comando digitado (não para do-files)."""
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
            raise StataError(602, f"file {path} already exists")
        s.cmdlog = {"path": path.resolve(), "fh": open(path, "a" if o.get("append") else "w",
                                                       encoding="utf-8"), "on": True}
        s.output.write(f"(cmdlog {path.resolve()} opened)\n", "text")
        return
    if cl is None:
        raise StataError(606, "no cmdlog file open")   # VERIFICAR
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
        # VERIFICAR: assert sem observações
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
            raise StataError(198, f"timer {nums[0]} already on")   # VERIFICAR
        t["start"] = time.perf_counter()
    elif sub == "off":
        if len(nums) != 1:
            raise StataError(198, "invalid syntax")
        t = ts.get(nums[0])
        if t is None or t["start"] is None:
            raise StataError(198, f"timer {nums[0]} not on")   # VERIFICAR
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
                      "text")   # VERIFICAR: layout
            s.r[f"t{k}"] = total
            s.r[f"nt{k}"] = float(t["count"])
    else:
        raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# creturn list
# ---------------------------------------------------------------------------

_CRETURN_GROUPS = [
    ("System values", ["current_date", "current_time", "rmsg_time", "stata_version", "version",
                       "opendta_version", "N", "k", "width", "changed", "filename", "rc"]),
    ("Directories and paths", ["pwd", "dirsep"]),
    ("System limits and constants", ["maxdouble", "mindouble", "epsdouble", "smallestdouble", "pi"]),
    ("Settings", ["os", "level", "more", "dp", "linesize", "type", "seed"]),
    ("Other", ["alpha", "ALPHA", "Mons", "Weekdays"]),
]


@command("creturn")
def cmd_creturn(s: "Session", args: str) -> None:
    if args.strip() not in ("list", "l", "li", "lis"):
        raise StataError(198, "invalid syntax")
    out = s.output
    for title, names in _CRETURN_GROUPS:
        out.write("\n" + "-" * 79 + "\n", "text")
        out.write(f"    {title}\n", "text")
        out.write("-" * 79 + "\n", "text")
        for n in names:
            v = s.creturn(n)
            if isinstance(v, str):
                shown = f'"{v}"'
            elif v >= M.SYSMISS:
                shown = "."
            else:
                shown = format_value(float(v), "%18.0g", pad=False).strip()
            out.write(f"{'c(' + n + ')':>22} = ", "text")
            out.write(shown + "\n", "result")


