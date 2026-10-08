"""Linha de comando do OpenDTA.

    opendta                 abre a interface gráfica
    opendta --console       console interativo (prompt ". ")
    opendta -b do arq.do    modo batch: executa e grava arq.log, como `stata -b`
    opendta -e do arq.do    igual a -b (aceito por compatibilidade)
"""

from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

from . import __version__
from .core.errors import ExitRequest
from .session import Session


def _console_listener(text: str, style: str) -> None:
    sys.stdout.write(text)
    sys.stdout.flush()


def run_console() -> int:
    s = Session()
    s.output.add_listener(_console_listener)
    print(f"OpenDTA {__version__} — console. Digite 'exit' para sair.")
    while True:
        try:
            line = input("\n. ")
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        if not line.strip():
            continue
        try:
            s.run(line)
        except ExitRequest:
            return 0


def _log_header(path: Path) -> str:
    now = _dt.datetime.now()
    stamp = now.strftime("%d %b %Y, %H:%M:%S").lstrip("0")
    return (
        "-" * 79 + "\n"
        "      name:  <unnamed>\n"
        f"       log:  {path}\n"
        "  log type:  text\n"
        f" opened on:  {stamp}\n"
    )


def run_batch(dofile: str, args: list[str]) -> int:
    """Executa um do-file e grava <nome>.log no diretório atual."""
    do_path = Path(dofile)
    if do_path.suffix == "":
        do_path = do_path.with_suffix(".do")
    log_path = Path.cwd() / (do_path.stem + ".log")
    s = Session()
    with log_path.open("w", encoding="utf-8") as log:
        log.write(_log_header(log_path))
        s.output.add_listener(lambda text, style: log.write(text))
        cmd = f'do "{do_path}"' + ("" if not args else " " + " ".join(args))
        s.output.echo_command(cmd)
        try:
            rc = s.run(cmd)
        except ExitRequest as e:
            rc = e.rc
        s.output.end_command()
    return rc


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="opendta", description="OpenDTA")
    p.add_argument("--console", action="store_true", help="console interativo, sem janela")
    p.add_argument("-b", "-e", dest="batch", action="store_true", help="modo batch")
    p.add_argument("--version", action="version", version=f"OpenDTA {__version__}")
    p.add_argument("rest", nargs=argparse.REMAINDER)
    ns = p.parse_args(argv)

    if ns.batch:
        rest = list(ns.rest)
        if rest and rest[0] in ("do", "run"):
            rest = rest[1:]
        if not rest:
            p.error("informe o do-file: opendta -b do arquivo.do")
        return run_batch(rest[0], rest[1:])
    if ns.console:
        return run_console()

    from .gui.app import run_gui
    return run_gui(ns.rest)


if __name__ == "__main__":
    raise SystemExit(main())
