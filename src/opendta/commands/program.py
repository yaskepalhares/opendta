"""do, run, include, exit, error, continue, args, tokenize, confirm."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import (BreakLoop, ContinueLoop, ExitRequest, StataError,
                           STANDARD_MESSAGES)
from ..lang.lexer import split_commands
from ..lang.words import split_words
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _split_path(args: str) -> tuple[str, str]:
    a = args.strip()
    if a.startswith('`"'):
        j = a.find("\"'")
        return a[2:j], a[j + 2:].strip()
    if a.startswith('"'):
        j = a.find('"', 1)
        if j == -1:
            raise StataError(198, "invalid syntax")
        return a[1:j], a[j + 1:].strip()
    path, _, rest = a.partition(" ")
    return path, rest.strip()


def _resolve_do(path: str) -> Path:
    p = Path(os.path.expanduser(path))
    if p.suffix == "":
        p = p.with_suffix(".do")
    if not p.exists():
        raise StataError(601, f"file {p} not found")
    return p


def _run_file(s: "Session", args: str, *, echo: bool, new_scope: bool) -> None:
    path_text, rest = _split_path(args)
    if "," in rest and rest.strip().startswith(","):
        rest = ""  # opções (nostop) chegam depois
    path = _resolve_do(path_text)
    try:
        source = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        source = path.read_text(encoding="latin-1")
    lines = split_commands(source)

    out = s.output
    scope = None
    if new_scope:
        argv = split_words(rest)
        frame = {str(k): v for k, v in enumerate(argv, start=1)}
        if rest:
            frame["0"] = rest
        s.macros.push_frame(frame)
        from ..lang.programs import Scope
        scope = Scope("do", str(path))
        s.scopes.append(scope)
    if not echo:
        out.quiet_depth += 1
    s.do_depth += 1
    saved_file = s.current_dofile
    s.current_dofile = str(path)
    if echo:
        out.end_command()
    try:
        s.interp.run_lines(lines, echo=echo)
    except ExitRequest as e:
        if e.clear_all:
            raise
        if e.rc:
            raise StataError(e.rc, "")
    except (ContinueLoop, BreakLoop):
        raise StataError(198, "continue outside loop")
    except StataError as e:
        if echo:
            s.report_error(e)
            out.write("\nend of do-file\n", "text", force=True)
            raise StataError(e.rc, "")
        raise
    finally:
        s.do_depth -= 1
        s.current_dofile = saved_file
        if not echo:
            out.quiet_depth -= 1
        if new_scope:
            s.macros.pop_frame()
            s.scopes.remove(scope)
            s.close_scope(scope)      # temporários e preserve do do-file
    if echo:
        out.ensure_line_start()
        out.write(". \nend of do-file\n", "text")


@command("do")
def cmd_do(s: "Session", args: str) -> None:
    _run_file(s, args, echo=True, new_scope=True)


@command("run")
def cmd_run(s: "Session", args: str) -> None:
    _run_file(s, args, echo=False, new_scope=True)


@command("include")
def cmd_include(s: "Session", args: str) -> None:
    _run_file(s, args, echo=True, new_scope=False)


@command("exit")
def cmd_exit(s: "Session", args: str) -> None:
    t = args.strip()
    clear_all = False
    if "," in t:
        t, _, opts = t.partition(",")
        clear_all = "clear" in opts.split() or "STATA" in opts.split()
    t = t.strip()
    rc = int(float(s.eval(t))) if t else 0
    in_program = any(sc.kind == "program" for sc in s.scopes)
    if s.do_depth == 0 and not in_program:
        clear_all = True
    raise ExitRequest(rc, clear_all=clear_all)


@command("error")
def cmd_error(s: "Session", args: str) -> None:
    t = args.strip()
    if not t:
        raise StataError(198, "invalid syntax")
    rc = int(float(s.eval(t)))
    if rc == 0:
        return
    # erro pedido de propósito (error #): sem explicação do OpenDTA
    raise StataError(rc, STANDARD_MESSAGES.get(rc, ""), hint="")


@command("continue")
def cmd_continue(s: "Session", args: str) -> None:
    t = args.strip()
    if t in ("", ","):
        raise ContinueLoop()
    if re.match(r"^,\s*break$", t):
        raise BreakLoop()
    raise StataError(198, "invalid syntax")


@command("args")
def cmd_args(s: "Session", args: str) -> None:
    for k, name in enumerate(args.split(), start=1):
        s.macros.set_local(name, s.macros.get_local(str(k)))


@command("tokenize")
def cmd_tokenize(s: "Session", args: str) -> None:
    text = args
    parse_chars = None
    m = re.search(r",\s*p(?:a|ar|ars|arse)?\(\s*\"?(.*?)\"?\s*\)\s*$", text)
    if m:
        parse_chars = m.group(1)
        text = text[:m.start()]
    text = text.strip()
    if text.startswith('`"') and text.endswith("\"'"):
        text = text[2:-2]
    elif len(text) >= 2 and text.startswith('"') and text.endswith('"'):
        text = text[1:-1]
    loc = s.macros.locals
    for k in [k for k in loc if k.isdigit()]:
        del loc[k]
    if parse_chars:
        pattern = "([" + re.escape(parse_chars) + "])"
        pieces = [p.strip() for p in re.split(pattern, text) if p.strip()]
    else:
        pieces = split_words(text)
    for k, w in enumerate(pieces, start=1):
        s.macros.set_local(str(k), w)


@command("confirm", "conf")
def cmd_confirm(s: "Session", args: str) -> None:
    t = args.strip()
    m = re.match(r"^(?:(integer)\s+)?n(?:u|um|umb|umbe|umber)?\s+(.*)$", t)
    if m:
        val = m.group(2).strip()
        try:
            x = float(val)
        except ValueError:
            raise StataError(7, f"'{val}' found where number expected" if not m.group(1)
                             else f"'{val}' found where integer expected")
        if m.group(1) and x != int(x):
            raise StataError(7, f"'{val}' found where integer expected")
        return
    m = re.match(r"^(new\s+)?f(?:i|il|ile)?\s+(.*)$", t)
    if m:
        path = m.group(2).strip().strip('"')
        exists = Path(path).exists()
        if m.group(1) and exists:
            raise StataError(602, f"file {path} already exists")
        if not m.group(1) and not exists:
            raise StataError(601, f"file {path} not found")
        return
    m = re.match(r"^names?\s+(.*)$", t)
    if m:
        for name in m.group(1).split():
            if not re.match(r"^[A-Za-z_][A-Za-z0-9_]{0,31}$", name):
                raise StataError(7, f"{name} invalid name")
        return
    raise StataError(198, "invalid syntax")
