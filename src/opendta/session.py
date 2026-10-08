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
from .core.dataset import Dataset
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
    "dp": "period",
    "hints": "on",       # OpenDTA: explicação depois das mensagens de erro
}


class EvalContext:
    """Liga o avaliador de expressões ao estado da sessão."""

    def __init__(self, session: "Session"):
        self.s = session

    def _variable(self, name: str):
        ds = self.s.data
        if ds.has(name):
            return ds.get(name)
        if self.s.settings.get("varabbrev", "on") == "on" and name not in self.s.scalars:
            matches = [n for n in ds.names if n.startswith(name)]
            if len(matches) == 1:
                return ds.get(matches[0])
            if len(matches) > 1:
                raise StataError(111, f"{name} ambiguous abbreviation")
        return None

    def resolve_name(self, name: str) -> Value:
        """Fora do contexto de observações (display, scalar, if de programação),
        uma variável vale o seu valor na 1ª observação, como x[1]."""
        s = self.s
        if s.data.has(name):
            return self.resolve_subscript(name, 1.0)
        if name in s.scalars:
            return s.scalars[name]
        if name == "_pi":
            return math.pi
        if name == "_rc":
            return float(s.rc)
        if name == "_N":
            return float(s.data.nobs)
        if name == "_n":
            return 1.0
        if self._variable(name) is not None:
            return self.resolve_subscript(name, 1.0)
        raise StataError(111, f"{name} not found")

    def resolve_subscript(self, name: str, index: Any) -> Value:
        if name in ("_b", "_se", "_coef"):
            raise StataError(111, f"[{index}] not found")
        var = self._variable(name)
        if var is None:
            raise StataError(111, f"{name} not found")
        if isinstance(index, str):
            raise StataError(109, "type mismatch")
        k = int(index) if not M.is_missing(index) else 0
        if 1 <= k <= len(var):
            return var.value(k - 1)
        return "" if var.is_string else M.SYSMISS

    def resolve_result(self, kind: str, raw: str) -> Value:
        s = self.s
        if kind == "c":
            return s.creturn(raw)
        store = {"r": s.r, "e": s.e, "s": s.sret}[kind]
        if raw in store:
            v = store[raw]
            if hasattr(v, "rows") and hasattr(v, "data"):
                raise StataError(109, "type mismatch")   # matriz numa expressão escalar
            return v
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
        self.data = Dataset()
        from .core import rng
        rng.reset(int(self.settings["seed"]))
        self.by_groups = None   # Groups ativo durante um prefixo by
        self.context = EvalContext(self)
        self._state_listeners: list[Callable[[], None]] = []
        self.do_depth = 0
        self.current_dofile = ""
        # ganchos da interface gráfica (por exemplo "browse"); sem GUI, ficam vazios
        self.ui_hooks: dict[str, Callable[..., None]] = {}
        # programas do usuário e escopos (programas e do-files em execução)
        from .lang.programs import Program, Scope
        self.programs: dict[str, Program] = {}
        self.scopes: list[Scope] = []

        from .lang.interpreter import Interpreter
        from . import commands  # noqa: F401  (registra os comandos)
        self.interp = Interpreter(self)

    # -- notificações para a interface --------------------------------------
    def add_state_listener(self, fn: Callable[[], None]) -> None:
        self._state_listeners.append(fn)

    def notify_state(self) -> None:
        for fn in list(self._state_listeners):
            fn()

    # -- escopos ------------------------------------------------------------------
    def close_scope(self, scope) -> None:
        """Desfaz o que o programa/do-file deixou: temporários e preserve."""
        ds = self.data
        if scope.tempvars:
            gone = [n for n in scope.tempvars if ds.has(n)]
            if gone:
                changed = ds.changed
                ds.drop_vars(gone)
                ds.changed = changed
        for name in scope.tempnames:
            self.scalars.pop(name, None)
            getattr(self, "matrices", {}).pop(name, None)
        for path in scope.tempfiles:
            for p in (path, path.with_suffix(".dta")):
                try:
                    p.unlink()
                except OSError:
                    pass
        if scope.preserved is not None:
            from .commands.preserve import restore_preserved
            restore_preserved(self, scope)

    # -- utilidades ---------------------------------------------------------
    def set_rc(self, rc: int) -> None:
        self.rc = int(rc)

    def expand(self, text: str) -> str:
        return expand(text, self.macros,
                      eval_inline=self._inline_eval,
                      extended=self._extended,
                      results=self._result_text)

    def _result_text(self, kind: str, name: str) -> str:
        """`r(x)' e afins: texto da macro ou do escalar; vazio se não existir."""
        if kind == "c":
            v = self.creturn(name)
        else:
            store = {"r": self.r, "e": self.e, "s": self.sret}[kind]
            if name not in store:
                return ""
            v = store[name]
        if isinstance(v, str):
            return v
        from .commands.matrix import Matrix
        if isinstance(v, Matrix):
            return "matrix"     # VERIFICAR
        return number_to_macro(float(v))

    def _inline_eval(self, text: str) -> str:
        v = self.eval(text)
        return v if isinstance(v, str) else number_to_macro(v)

    def _extended(self, text: str) -> str:
        from .commands.macro import extended_function
        return extended_function(self, text)

    def eval(self, text: str) -> Value:
        return evaluate(parse(text), self.context)

    def expand_varlist(self, text: str) -> list[str]:
        from .core.varlist import expand as _expand
        return _expand(self.data, text, abbrev=self.settings.get("varabbrev", "on") == "on")

    @property
    def nobs(self) -> int:
        return self.data.nobs

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
            "N": float(self.data.nobs),
            "k": float(self.data.nvars),
            "width": float(self.data.width()),
            "changed": float(self.data.changed),
            "rc": float(self.rc),
            "filename": self.data.filename,
            "level": float(self.settings["level"]),
            "more": self.settings["more"],
            "dp": self.settings.get("dp", "period"),
            "linesize": float(self.settings["linesize"]),
            "type": self.settings["type"],
            "seed": _rng_state(),
            "rngstate": _rng_state(),
            "rng": self.settings.get("rng", "default"),
            "rng_current": "mt64",
            "alpha": "a b c d e f g h i j k l m n o p q r s t u v w x y z",
            "ALPHA": "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z",
            "Mons": "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec",
            "Weekdays": "Sun Mon Tue Wed Thu Fri Sat",
        }
        if name not in values:
            return M.SYSMISS
        return values[name]

    # -- by -----------------------------------------------------------------
    def run_by(self, bp, command_text: str) -> None:
        """Executa `command_text` com o prefixo by/bysort."""
        from .core import sorting
        from .lang.vexpr import Groups

        ds = self.data
        keys = self.expand_varlist(" ".join(bp.keys))
        extra = self.expand_varlist(" ".join(bp.sort_extra)) if bp.sort_extra else []
        if bp.sort:
            sorting.sort(ds, keys + extra)
        elif not sorting.is_sorted(ds, keys + extra):
            raise StataError(5, "not sorted")
        groups = Groups.from_columns([ds.get(k).data for k in keys])
        saved = self.by_groups, getattr(self, "_by_keys", [])
        self.by_groups, self._by_keys = groups, keys
        try:
            self.interp._execute_expanded(command_text)
        finally:
            self.by_groups, self._by_keys = saved
        self.notify_state()

    # -- execução -------------------------------------------------------------
    def run(self, source: str, *, echo: bool = False) -> int:
        """Executa texto (uma ou mais linhas) e devolve o código de retorno.
        Erros são impressos como no Stata."""
        try:
            lines = split_commands(source)
            self.interp.run_lines(lines, echo=echo)
        except StataError as e:
            self.report_error(e)
            self.set_rc(e.rc)
            return e.rc
        return 0

    def report_error(self, e: StataError, *, show_rc: bool = True) -> None:
        """Mensagem de erro (como no Stata), explicação do OpenDTA e r(#)."""
        lines: list[str] = []
        if self.settings.get("hints", "on") == "on":
            from .core.hints import explain
            try:
                lines = explain(e, self)
            except Exception:  # noqa: BLE001 - uma explicação nunca derruba o erro original
                lines = []
        self.output.error(e.message, e.rc, show_rc=show_rc, hints=lines)

    def run_command(self, line: str) -> int:
        """Comando digitado na janela Command: eco '. linha' e execução."""
        from .commands.logcmd import record_command
        record_command(self, line)
        self.output.echo_command(line)
        try:
            rc = self.run(line, echo=False)
        except ExitRequest:
            raise
        finally:
            self.output.end_command()
            self.notify_state()
        return rc

    def run_text(self, text: str) -> int:
        """Várias linhas digitadas ou coladas na janela Command: rodam como um
        trecho de do-file (blocos inteiros, eco linha a linha)."""
        from .commands.logcmd import record_command
        record_command(self, text)
        try:
            self.interp.run_lines(split_commands(text), echo=True)
            rc = 0
        except StataError as e:
            self.report_error(e)
            self.set_rc(e.rc)
            rc = e.rc
        finally:
            self.notify_state()
        return rc

    def run_file(self, path: str | Path, args: str = "", *, echo: bool = True) -> int:
        quoted = f'"{path}"' + (f" {args}" if args else "")
        cmd = "do" if echo else "run"
        return self.run(f"{cmd} {quoted}")


def _rng_state() -> str:
    from .core import rng
    return rng.RNG.state()
