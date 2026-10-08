"""Compara a saída do OpenDTA com logs de referência gerados no Stata 14.

Uso:
    python tools/compare.py              # roda todos os compat/do/*.do
    python tools/compare.py 0001         # só os que começam com 0001
    python tools/compare.py --show-diff  # mostra as diferenças linha a linha

Para cada compat/do/NOME.do:
  1. executa no OpenDTA em modo batch, gravando compat/out/NOME.log;
  2. se existir compat/expected/NOME.log (gerado no Stata), compara os dois.

Antes de comparar, cabeçalho e rodapé do log, linhas `. do ...`/`. log close`
e espaços no fim das linhas são descartados; linhas em branco são ignoradas.
"""

from __future__ import annotations

import argparse
import difflib
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

COMPAT = ROOT / "compat"
_HEADER = re.compile(r"^\s*(name|log|log type|opened on|closed on):")
_SKIP = re.compile(r"^\.\s+((capture\s+)?(noisily\s+)?(do|run)|log close|log using|set linesize)\b")


def normalize(text: str) -> list[str]:
    out: list[str] = []
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line or set(line) <= {"-"} or _HEADER.match(line) or _SKIP.match(line):
            continue
        if line in (".",):
            continue
        out.append(line)
    return out


def run_opendta(dofile: Path, outdir: Path) -> Path:
    from opendta.cli import run_batch

    outdir.mkdir(parents=True, exist_ok=True)
    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        run_batch(str(dofile), [])
    finally:
        os.chdir(cwd)
    return outdir / (dofile.stem + ".log")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("prefix", nargs="?", default="")
    p.add_argument("--show-diff", action="store_true")
    ns = p.parse_args(argv)

    files = sorted(f for f in (COMPAT / "do").glob("*.do") if f.name.startswith(ns.prefix))
    passed = failed = missing = 0
    for f in files:
        got = run_opendta(f, COMPAT / "out")
        exp = COMPAT / "expected" / (f.stem + ".log")
        if not exp.exists():
            print(f"  ?  {f.stem:<32} sem log de referência do Stata")
            missing += 1
            continue
        a = normalize(exp.read_text(encoding="utf-8", errors="replace"))
        b = normalize(got.read_text(encoding="utf-8", errors="replace"))
        if a == b:
            print(f"  ok {f.stem}")
            passed += 1
        else:
            print(f"  XX {f.stem}")
            failed += 1
            if ns.show_diff:
                for line in difflib.unified_diff(a, b, "Stata 14", "OpenDTA", lineterm="", n=1):
                    print("     " + line)
    print(f"\n{passed} iguais, {failed} diferentes, {missing} sem referência")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
