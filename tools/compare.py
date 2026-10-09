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
_SKIP = re.compile(r"^\.\s+((capture\s+)?(noisily\s+)?(do|run)|log close|log using|set linesize|quietly set (dp|hints|linesize))\b")


# carimbo de data do .dta (describe): muda a cada execução
_STAMP = re.compile(r"\b\d{1,2} [A-Z][a-z]{2} \d{4} \d{2}:\d{2}\b")


def normalize(text: str) -> list[str]:
    out: list[str] = []
    for raw in text.splitlines():
        line = _STAMP.sub("<data>", raw.rstrip())
        if not line or set(line) <= {"-"} or _HEADER.match(line) or _SKIP.match(line):
            continue
        if line in (".",):
            continue
        out.append(line)
    # o modo batch termina um do-file interrompido com "end of do-file" e r(#);
    # no Stata as referências rodam sob capture noisily, sem o r(#) final
    if len(out) >= 2 and out[-2] == "end of do-file" and re.fullmatch(r"r\(\d+\);", out[-1]):
        out.pop()
    return out


def run_opendta(dofile: Path, outdir: Path, setup: str = "") -> Path:
    from opendta.cli import run_batch

    outdir.mkdir(parents=True, exist_ok=True)
    # explicações de erro do OpenDTA desligadas: o Stata não as tem.
    # `setup` (ex.: set dp comma) roda antes, sem aparecer no log comparado
    target = outdir / dofile.name
    # gerar_esperados.do roda os casos com set linesize 255
    pre = ("quietly set hints off\nquietly set linesize 255\n"
           + (f"quietly {setup}\n" if setup else ""))
    target.write_text(pre + dofile.read_text(encoding="utf-8"), encoding="utf-8")
    # restos de um caso interrompido (odta_*) não podem afetar o seguinte;
    # gerar_esperados.do faz a mesma limpeza no Stata
    for leftover in outdir.glob("odta_*"):
        if leftover.is_file():
            leftover.unlink()
    cwd = os.getcwd()
    os.chdir(outdir)
    try:
        run_batch(str(target), [])
    finally:
        os.chdir(cwd)
    return outdir / (dofile.stem + ".log")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("prefix", nargs="?", default="")
    p.add_argument("--show-diff", action="store_true")
    p.add_argument("--expected", default=str(COMPAT / "expected"),
                   help="pasta com os logs de referência (padrão: compat/expected)")
    p.add_argument("--setup", default="",
                   help="comando rodado antes de cada do-file, ex.: 'set dp comma'")
    ns = p.parse_args(argv)

    files = sorted(f for f in (COMPAT / "do").glob("*.do") if f.name.startswith(ns.prefix))
    passed = failed = missing = 0
    for f in files:
        got = run_opendta(f, COMPAT / "out", ns.setup)
        exp = Path(ns.expected) / (f.stem + ".log")
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
