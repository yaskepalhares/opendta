"""Session: the state of one OpenDTA instance.

Equivalent to what lives in the memory of a Stata window: macros, scalars,
r() and e() results, settings (set), the working directory and
the dataset.
"""

from __future__ import annotations

import datetime as _dt
import functools
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
    # set when Stata starts (compat/expected/2003)
    "S_OS": {"Windows": "Windows", "Darwin": "MacOSX"}.get(platform.system(), "Unix"),
    "S_OSDTL": platform.mac_ver()[0] if platform.system() == "Darwin" else platform.release(),
    "S_StataSE": "SE",
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
    "hints": "on",       # OpenDTA: explanation after error messages
    "superscript": "on",  # OpenDTA: 1.0000×10¹⁰ instead of 1.00000e+10 on screen
    "numerics": "stata",  # OpenDTA: precise = contas mais exatas (estimation/numerics.py)
    "graphics": "on",
    "scheme": "s2color",
}


@functools.lru_cache(maxsize=1)
def _locale() -> str:
    """en_ and the region of the system locale (Stata 14 on a Mac set to Brazil: en_BR)."""
    import locale
    import subprocess
    loc = ""
    if platform.system() == "Darwin":
        try:
            loc = subprocess.run(["defaults", "read", "-g", "AppleLocale"], capture_output=True,
                                 text=True, timeout=2).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            loc = ""
    elif platform.system() == "Windows":
        try:   # "en-US"; locale.getlocale() gives "English_United States" there
            import ctypes
            buf = ctypes.create_unicode_buffer(85)
            if ctypes.windll.kernel32.GetUserDefaultLocaleName(buf, 85):
                loc = buf.value.replace("-", "_")
        except (OSError, AttributeError):
            loc = ""
    loc = loc or locale.getlocale()[0] or os.environ.get("LANG", "")
    region = loc.split(".")[0].split("@")[0].split("_")[1] if "_" in loc else "US"
    if not (len(region) == 2 and region.isalpha()):
        region = "US"
    return "en_" + region


class EvalContext:
    """Links the expression evaluator to the session state."""

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
        """Outside the observation context (display, scalar, programming if),
        a variable is worth its value in the 1st observation, like x[1]."""
        s = self.s
        if s.data.has(name):
            return self.resolve_subscript(name, 1.0)
        if "." in name:
            from .lang.vexpr import ts_value
            v = ts_value(s, name)
            return float(v[0]) if len(v) else M.SYSMISS
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
            return self._coef(name, index)
        var = self._variable(name)
        if var is None:
            raise StataError(111, f"{name} not found")
        if isinstance(index, str):
            raise StataError(109, "type mismatch")
        k = int(index) if not M.is_missing(index) else 0
        if 1 <= k <= len(var):
            return var.value(k - 1)
        return "" if var.is_string else M.SYSMISS

    def _coef(self, name: str, index: Any) -> Value:
        """_b[x], _se[x], _b[eq:x], [eq]_b[x] of the last estimation."""
        from .estimation.postest import _coef_index
        from .estimation.results import current
        import numpy as _np
        est = current(self.s)
        text = str(index).strip()
        eq = ""
        if ":" in text:
            eq, text = text.split(":", 1)
        try:
            i = _coef_index(est, self.s, text, eq.strip())
        except StataError:
            raise StataError(111, f"[{text}] not found")
        if name == "_se":
            v = float(est.V[i, i])
            return float(_np.sqrt(v)) if v > 0 else 0.0
        return float(est.b[i])

    def resolve_result(self, kind: str, raw: str) -> Value:
        if kind == "e" and raw.strip() == "sample":
            from .estimation.postest import esample
            m = esample(self.s)
            return 1.0 if len(m) and m[0] else 0.0
        s = self.s
        if kind == "c":
            return s.creturn(raw, strict=True)
        store = {"r": s.r, "e": s.e, "s": s.sret}[kind]
        if raw in store:
            v = store[raw]
            if hasattr(v, "rows") and hasattr(v, "data"):
                raise StataError(109, "type mismatch")   # matrix in a scalar expression
            return v
        return M.SYSMISS

    def call_function(self, name: str, args: list[Value]) -> Value:
        if name == "scalar" and len(args) == 1:
            raise StataError(198, "invalid syntax")
        return F.call(name, args)


class Session:
    def __init__(self) -> None:
        self.output = Output()
        self.output.linesize = lambda: int(float(self.settings.get("linesize", 80)))
        self.macros = MacroStore()
        self.macros.globals.update(DEFAULT_GLOBALS)
        self.scalars: dict[str, Value] = {}
        self.r: dict[str, Value] = {}
        self.e: dict[str, Value] = {}
        self.sret: dict[str, Value] = {}
        self.rc = 0
        self.settings: dict[str, str] = dict(DEFAULT_SETTINGS)
        self.version = 14.0
        # display ..., _request(): who supplies the typed line (the window sets a dialog; the
        # console and batch mode read standard input)
        self.input_provider = None
        self.data = Dataset()
        from .core import rng
        rng.reset(int(self.settings["seed"]))
        from .core.formats import set_superscript
        set_superscript(self.settings["superscript"] == "on")
        from .estimation.numerics import set_precise
        set_precise(self.settings["numerics"] == "precise")
        self.by_groups = None   # Groups active during a by prefix
        self.context = EvalContext(self)
        self._mata = None
        self._state_listeners: list[Callable[[], None]] = []
        self.do_depth = 0
        # Review in the interface: when review_root is a dict, the lines run by do-files
        # (do, run, include) are recorded in it as a tree, {"children": [{"text", "rc",
        # "children"}]}; _review_stack holds the entries being run (the parents of nested files)
        self.review_root: dict | None = None
        self._review_stack: list[dict] = []
        self.current_dofile = ""
        # graphical interface hooks (for example "browse"); without a GUI, they stay empty
        self.ui_hooks: dict[str, Callable[..., None]] = {}
        # user programs and scopes (programs and do-files running)
        from .lang.programs import Program, Scope
        self.programs: dict[str, Program] = {}
        self.scopes: list[Scope] = []

        from .lang.interpreter import Interpreter
        from . import commands  # noqa: F401  (registers the commands)
        self.interp = Interpreter(self)

    # -- notifications to the interface --------------------------------------
    def add_state_listener(self, fn: Callable[[], None]) -> None:
        self._state_listeners.append(fn)

    def default_type(self) -> str:
        """Type of new variables: set type; with set numerics precise, double
        instead of float."""
        vt = self.settings.get("type", "float")
        if vt == "float" and self.settings.get("numerics", "stata") == "precise":
            return "double"
        return vt

    def notify_state(self) -> None:
        for fn in list(self._state_listeners):
            fn()

    # -- escopos ------------------------------------------------------------------
    def close_scope(self, scope) -> None:
        """Undoes what the program/do-file left behind: temporaries and preserve."""
        ds = self.data
        if scope.tempvars:
            gone = [n for n in scope.tempvars if ds.has(n)]
            if gone:
                changed = ds.changed
                ds.drop_vars(gone)
                ds.changed = changed
        if scope.tempnames and getattr(self, "_classes", None) is not None:
            self._classes.drop_names(scope.tempnames)      # class objects named by tempname
        for name in scope.tempnames:
            self.scalars.pop(name, None)
            getattr(self, "matrices", {}).pop(name, None)
            getattr(self, "posts", {}).pop(name, None)
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
                      results=self._result_text,
                      objects=self._class_text)

    def _result_text(self, kind: str, name: str) -> str:
        """`r(x)' and the like: text of the macro or scalar; empty if it does not exist."""
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
            return "matrix"     # VERIFY
        return number_to_macro(float(v))

    def _class_text(self, text: str) -> str:
        from .lang.classeval import macro_text
        return macro_text(self, text)

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

    def creturn(self, name: str, *, strict: bool = False) -> Value:
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
            "flavor": "IC",   # VERIFY: OpenDTA has no Small/IC/SE/MP editions; reports IC
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
            "graphics": self.settings.get("graphics", "on"),
            "scheme": self.settings.get("scheme", "s2color"),
            "seed": _rng_state(),
            "rngstate": _rng_state(),
            "rng": self.settings.get("rng", "default"),
            "rng_current": __import__("opendta.core.rng", fromlist=["current"]).current(),
            "alpha": "a b c d e f g h i j k l m n o p q r s t u v w x y z",
            "ALPHA": "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z",
            "Mons": "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec",
            "Months": "January February March April May June July August September October "
                      "November December",
            "Wdays": "Sun Mon Tue Wed Thu Fri Sat",
            "Weekdays": "Sunday Monday Tuesday Wednesday Thursday Friday Saturday",
            "filedate": self.data.timestamp if self.data.filename else "",
            **self._c_system(),
        }
        if name not in values:
            if strict:
                # c(sortlist), c(eps), c(lastcmd)... do not exist in Stata 14 (compat 2003)
                raise StataError(133, f"c({name}) undefined")
            return M.SYSMISS
        return values[name]

    def _c_system(self) -> dict[str, Value]:
        """c() values that describe the system, the limits and the settings, with the values
        Stata 14 SE reports (compat/expected/2003)."""
        import getpass
        import socket
        import tempfile
        from .lang.adopath import sysdir

        def setting(name: str, default: str) -> str:
            return self.settings.get(name, default)

        def number(name: str, default: str) -> float:
            try:
                return float(self.settings.get(name, default))
            except ValueError:
                return float(default)

        dirs = sysdir(self)
        osdtl = platform.mac_ver()[0] if platform.system() == "Darwin" else platform.release()
        try:
            user = getpass.getuser()
        except (OSError, KeyError):
            user = ""
        # VERIFY: machine_type on Windows and Linux; Stata 14 on any Mac reports Intel
        machine = {"Darwin": "Macintosh (Intel 64-bit)", "Windows": "PC (64-bit x86-64)"}.get(
            platform.system(), "PC (64-bit x86-64)")
        return {
            "machine_type": machine,
            "byteorder": "lohi",
            "username": user,
            "hostname": socket.gethostname(),
            "osdtl": osdtl,
            # ends with / on every system, like c(dirsep) (VERIFY on Windows)
            "tmpdir": tempfile.gettempdir().rstrip("/\\") + "/",
            "sysdir_stata": dirs["STATA"].rstrip("/") + "/",
            "sysdir_base": dirs["BASE"].rstrip("/") + "/",
            "sysdir_site": dirs["SITE"].rstrip("/") + "/",
            "sysdir_plus": dirs["PLUS"].rstrip("/") + "/",
            "sysdir_personal": dirs["PERSONAL"].rstrip("/") + "/",
            "born_date": "02 Apr 2015",   # the Stata 14 build the references come from
            "SE": 1.0, "MP": 0.0, "processors": 1.0,
            "mode": "", "console": "",
            "memory": 33554432.0,   # VERIFY: Stata reports the memory it has allocated
            "max_N_theory": 2147483620.0, "max_k_theory": 5000.0, "max_width_theory": 1048576.0,
            "maxvar": number("maxvar", "5000"), "matsize": number("matsize", "400"),
            "epsfloat": 2.0 ** -23,
            "maxbyte": 100.0, "minbyte": -127.0, "maxint": 32740.0, "minint": -32767.0,
            "maxlong": 2147483620.0, "minlong": -2147483647.0,
            "maxfloat": 1.70141173319e+38, "minfloat": -1.70141173319e+38,
            "maxstrvarlen": 2045.0, "max_macrolen": 4227143.0, "macrolen": 645200.0,
            "namelen": 32.0, "cmdlen": 645216.0, "eqlen": 1337.0,
            "adosize": number("adosize", "1000"), "niceness": number("niceness", "5"),
            "pagesize": number("pagesize", "37"), "logtype": setting("logtype", "smcl"),
            "noisily": 0.0 if self.output.quiet_depth and not self.output.noisy_depth else 1.0,
            "trace": setting("trace", "off"), "tracedepth": number("tracedepth", "32000"),
            "varabbrev": setting("varabbrev", "on"), "rmsg": setting("rmsg", "off"),
            "httpproxy": setting("httpproxy", "off"), "checksum": setting("checksum", "off"),
            "reventries": number("reventries", "5000"), "virtual": setting("virtual", "off"),
            "scrollbufsize": number("scrollbufsize", "204800"),
            "emptycells": setting("emptycells", "keep"),
            "showbaselevels": setting("showbaselevels", ""),
            "showemptycells": setting("showemptycells", ""),
            "showomitted": setting("showomitted", ""),
            "fvwrap": number("fvwrap", "1"), "lstretch": setting("lstretch", ""),
            "cformat": setting("cformat", ""), "pformat": setting("pformat", ""),
            "sformat": setting("sformat", ""),
            # the rest of creturn list, with Stata 14's defaults (compat 2003)
            "rmsg_time": number("rmsg_time", "0"), "userversion": float(self.version),
            "bit": 64.0, "processors_lic": 1.0, "processors_mach": M.SYSMISS,
            "processors_max": 1.0, "sysdir_oldplace": "~/ado/",
            "adopath": "BASE;SITE;.;PERSONAL;PLUS;OLDPLACE",
            "max_matsize": 11000.0, "min_matsize": 10.0, "charlen": 67783.0,
            "max_cmdlen": 4227159.0, "namelenbyte": 128.0, "namelenchar": 32.0,
            "maxstrlvarlen": 2000000000.0, "maxvlabellen": 32000.0,
            "min_memory": number("min_memory", "0"),
            "max_memory": number("max_memory", str(M.SYSMISS)) if setting("max_memory", ".") != "."
            else M.SYSMISS,
            "segmentsize": number("segmentsize", "33554432"),
            "eolchar": setting("eolchar", "unix"), "notifyuser": setting("notifyuser", "on"),
            "playsnd": setting("playsnd", "off"), "include_bitmap": setting("include_bitmap", "on"),
            "clevel": number("clevel", "95"), "fvlabel": setting("fvlabel", "on"),
            "fvwrapon": setting("fvwrapon", "word"),
            "coeftabresults": setting("coeftabresults", "on"),
            "revkeyboard": setting("revkeyboard", "on"), "varkeyboard": setting("varkeyboard", "on"),
            "smoothfonts": setting("smoothfonts", "on"), "linegap": number("linegap", "1"),
            "maxdb": number("maxdb", "50"), "printcolor": setting("printcolor", "automatic"),
            "copycolor": setting("copycolor", "automatic"),
            "timeout1": number("timeout1", "30"), "timeout2": number("timeout2", "180"),
            "httpproxyhost": setting("httpproxyhost", ""),
            "httpproxyport": number("httpproxyport", "80"),
            "httpproxyauth": setting("httpproxyauth", "off"),
            "httpproxyuser": setting("httpproxyuser", ""), "httpproxypw": setting("httpproxypw", ""),
            "update_query": setting("update_query", "on"),
            "update_interval": number("update_interval", "7"),
            "update_prompt": setting("update_prompt", "off"),
            "tracesep": setting("tracesep", "on"), "traceindent": setting("traceindent", "on"),
            "traceexpand": setting("traceexpand", "on"), "tracenumber": setting("tracenumber", "off"),
            "tracehilite": setting("tracehilite", ""),
            "matastrict": setting("matastrict", "off"), "matalnum": setting("matalnum", "off"),
            "mataoptimize": setting("mataoptimize", "on"), "matafavor": setting("matafavor", "space"),
            "matacache": number("matacache", "2000"),
            # the list query shows (compat 2008)
            "matalibs": setting("matalibs", "lmatabase;lmatamcmc;lmatasem;lmatagsem;lmatapss;"
                                            "lmatapostest;lmatapath;lmataopt;lmataado;lmatafc"),
            "matamofirst": setting("matamofirst", "off"),
            # VERIFY: Stata built the locale from the interface language and the region
            "locale_ui": setting("locale_ui", _locale()), "locale_functions": setting(
                "locale_functions", _locale()), "locale_icudflt": _locale(),
            "maxiter": number("maxiter", "16000"), "searchdefault": setting("searchdefault", "all"),
            "odbcmgr": setting("odbcmgr", "iodbc"),
        }

    # -- Mata ---------------------------------------------------------------
    @property
    def mata(self):
        """Mata workspace (created on first use)."""
        if self._mata is None:
            from .mata.interp import MataEngine
            self._mata = MataEngine(self)
        return self._mata

    # -- by -----------------------------------------------------------------
    def run_by(self, bp, command_text: str) -> None:
        """Runs `command_text` with the by/bysort prefix."""
        from .core import sorting
        from .lang.vexpr import Groups

        ds = self.data
        keys = self.expand_varlist(" ".join(bp.keys))
        extra = self.expand_varlist(" ".join(bp.sort_extra)) if bp.sort_extra else []
        if bp.sort:
            sorting.sort(ds, keys + extra)
        elif ds.sortlist[:len(keys + extra)] != keys + extra:
            # Stata goes by the sort marker, not the order of the data (compat 1501)
            raise StataError(5, "not sorted")
        groups = Groups.from_columns([ds.get(k).data for k in keys])
        saved = self.by_groups, getattr(self, "_by_keys", [])
        self.by_groups, self._by_keys = groups, keys
        try:
            self.interp._execute_expanded(command_text)
        finally:
            self.by_groups, self._by_keys = saved
        self.notify_state()

    # -- execution -------------------------------------------------------------
    def run(self, source: str, *, echo: bool = False) -> int:
        """Runs text (one or more lines) and returns the return code.
        Errors are printed as in Stata."""
        try:
            lines = split_commands(source)
            self.interp.run_lines(lines, echo=echo)
        except StataError as e:
            self.report_error(e)
            self.set_rc(e.rc)
            return e.rc
        return 0

    def request_line(self) -> str:
        """A line typed by the user, for display _request()."""
        if self.input_provider is not None:
            text = self.input_provider()
            if text is None:
                raise StataError(1, "--Break--")   # cancelled: like Break
            return str(text)
        try:
            return input()
        except EOFError:
            return ""

    def report_error(self, e: StataError, *, show_rc: bool = True) -> None:
        """Error message (as in Stata), OpenDTA's explanation and r(#)."""
        hint: list[tuple[str, str]] = []
        if self.settings.get("hints", "on") == "on":
            from .core.hints import explain, render
            try:
                expl = explain(e, self)
                if expl is not None:
                    hint = render(expl, int(float(self.settings.get("linesize", 80))) - 1)
            except Exception:  # noqa: BLE001 - an explanation never brings down the original error
                hint = []
        from .i18n import translate_error
        self.output.error(translate_error(e.message), e.rc, show_rc=show_rc, hint=hint)

    def run_command(self, line: str) -> int:
        """Command typed in the Command window: echo '. line' and run."""
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
        """Several lines typed or pasted into the Command window: they run as a
        piece of a do-file (whole blocks, echoed line by line)."""
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
