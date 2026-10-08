"""Sessão: o estado de uma instância do OpenDTA.

Equivale ao que vive na memória de uma janela do Stata: macros, scalars,
resultados r() e e(), configurações (set), diretório de trabalho e, a partir
da fase 1, o conjunto de dados.
"""

from __future__ import annotations

import datetime as _dt
import math
import os
import platform
from pathlib import Path
from typing import Any, Callable

from . import __version__
from .core import missing as M
from .core.errors import ExitRequest, StataError
from .core.formats import number_to_macro
from .core.output import Output
from .lang import functions as F
from .lang.expr import parse, evaluate
from .lang.lexer import split_commands
from .lang.macros import MacroStore, expand

Value = Any

_MONTH_ABBR = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
               "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

DEFAULT_GLOBALS = {
    "S_ADO": "BASE;SITE;.;PERSONAL;PLUS;OLDPLACE",
    "S_level": "95",
    "F1": "help advice;",
    "F2": "describe;",
    "F7": "save ",
    "F8": "use ",
}

DEFAULT_SETTINGS: dict[str, str] = {
    "more": "off",
    "rmsg": "off",
    "linesize": "80",
    "type": "float",
    "level": "95",
    "varabbrev": "on",
    "trace": "off",
    "seed": "123456789",
}


class EvalContext:
    """Liga o avaliador de expressões ao estado da sessão."""

    def __init__(self, session: "Session"):
        self.s = session

    def resolve_name(self, name: str) -> Value:
        s = self.s
        # fase 1: variáveis têm precedência sobre scalars
        if name in s.scalars:
            return s.scalars[name]
        if name == "_pi":
            return math.pi
        if name == "_rc":
            return float(s.rc)
        if name == "_N":
            return float(s.nobs)
        if name == "_n":
            return 1.0
        raise StataError(111, f"{name} not found")

    def resolve_subscript(self, name: str, index: Any) -> Value:
        if name in ("_b", "_se", "_coef"):
            raise StataError(111, f"[{index}] not found")
        raise StataError(111, f"{name} not found")

    def resolve_result(self, kind: str, raw: str) -> Value:
        s = self.s
        if kind == "c":
            return s.creturn(raw)
        store = {"r": s.r, "e": s.e, "s": s.sret}[kind]
        if raw in store:
            return store[raw]
        return M.SYSMISS

    def call_function(self, name: str, args: list[Value]) -> Value:
        if name == "scalar" and len(args) == 1:
            raise StataError(198, "invalid syntax")
        return F.call(name, args)


class Session:
    def __init__(self) -> None:
        self.output = Output()
        self.macros = MacroStore()
        self.macros.globals.update(DEFAULT_GLOBALS)
        self.scalars: dict[str, Value] = {}
        self.r: dict[str, Value] = {}
        self.e: dict[str, Value] = {}
        self.sret: dict[str, Value] = {}
        self.rc = 0
        self.settings: dict[str, str] = dict(DEFAULT_SETTINGS)
        self.version = 14.0
        self.nobs = 0          # fase 1: dados em memória
        self.context = EvalContext(self)
        self._state_listeners: list[Callable[[], None]] = []
        self.do_depth = 0
        self.current_dofile = ""

        from .lang.interpreter import Interpreter
        from . import commands  # noqa: F401  (registra os comandos)
        self.interp = Interpreter(self)

    # -- notificações para a interface --------------------------------------
    def add_state_listener(self, fn: Callable[[], None]) -> None:
        self._state_listeners.append(fn)

    def notify_state(self) -> None:
        for fn in list(self._state_listeners):
            fn()

    # -- utilidades ---------------------------------------------------------
    def set_rc(self, rc: int) -> None:
        self.rc = int(rc)

    def expand(self, text: str) -> str:
        return expand(text, self.macros,
                      eval_inline=self._inline_eval,
                      extended=self._extended)

    def _inline_eval(self, text: str) -> str:
        v = self.eval(text)
        return v if isinstance(v, str) else number_to_macro(v)

    def _extended(self, text: str) -> str:
        from .commands.macro import extended_function
        return extended_function(self, text)

    def eval(self, text: str) -> Value:
        return evaluate(parse(text), self.context)

    def expand_varlist(self, text: str) -> list[str]:
        raise StataError(111, "variable not found (dados chegam na fase 1)")

    def creturn(self, name: str) -> Value:
        now = _dt.datetime.now()
        values: dict[str, Value] = {
            "pi": math.pi,
            "current_date": f"{now.day:2d} {_MONTH_ABBR[now.month - 1]} {now.year}",
            "current_time": now.strftime("%H:%M:%S"),
            "pwd": os.getcwd(),
            "dirsep": "/",
            "os": {"Windows": "Windows", "Darwin": "MacOSX"}.get(platform.system(), "Unix"),
            "version": self.version,
            "stata_version": 14.0,
            "opendta_version": __version__,
            "maxdouble": M.MAXDOUBLE,
            "mindouble": M.MINDOUBLE,
            "epsdouble": M.EPSDOUBLE,
            "smallestdouble": 2.0 ** -1022,
            "N": float(self.nobs),
            "k": 0.0,
            "rc": float(self.rc),
            "filename": "",
            "level": float(self.settings["level"]),
            "more": self.settings["more"],
            "linesize": float(self.settings["linesize"]),
            "type": self.settings["type"],
            "seed": self.settings["seed"],
            "alpha": "a b c d e f g h i j k l m n o p q r s t u v w x y z",
            "ALPHA": "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z",
            "Mons": "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec",
            "Weekdays": "Sun Mon Tue Wed Thu Fri Sat",
        }
        if name not in values:
            return M.SYSMISS
        return values[name]

    # -- execução -------------------------------------------------------------
    def run(self, source: str, *, echo: bool = False) -> int:
        """Executa texto (uma ou mais linhas) e devolve o código de retorno.
        Erros são impressos como no Stata."""
        try:
            lines = split_commands(source)
            self.interp.run_lines(lines, echo=echo)
        except StataError as e:
            self.output.error(e.message, e.rc)
            self.set_rc(e.rc)
            return e.rc
        return 0

    def run_command(self, line: str) -> int:
        """Comando digitado na janela Command: eco '. linha' e execução."""
        self.output.echo_command(line)
        try:
            rc = self.run(line, echo=False)
        except ExitRequest:
            raise
        finally:
            self.notify_state()
        return rc

    def run_file(self, path: str | Path, args: str = "", *, echo: bool = True) -> int:
        quoted = f'"{path}"' + (f" {args}" if args else "")
        cmd = "do" if echo else "run"
        return self.run(f"{cmd} {quoted}")
