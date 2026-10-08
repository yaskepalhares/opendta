"""Pastas do sistema (sysdir) e caminho de busca de .ado (adopath).

Ao encontrar um comando que não é interno nem um programa já definido, o
OpenDTA procura `nome.ado` nas pastas do adopath (em cada pasta e na
subpasta com a primeira letra do nome, como plus/o/outreg2.ado), roda o
arquivo e chama o programa que ele define. Isso permite usar pacotes de
terceiros (SSC, pacotes de artigos) instalados nessas pastas.

Pastas padrão (podem ser mudadas com `sysdir set`):

    BASE      .ado próprios do OpenDTA (src/opendta/ado)
    SITE      <pasta do OpenDTA>/site
    PERSONAL  ~/ado/personal
    PLUS      ~/ado/plus
    OLDPLACE  ~/ado

Os .ado oficiais da StataCorp nunca entram no caminho (política
clean-room); os de terceiros, sim.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..session import Session

_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PATH = ["BASE", "SITE", ".", "PERSONAL", "PLUS", "OLDPLACE"]


def sysdir(s: "Session") -> dict[str, str]:
    if not hasattr(s, "sysdirs"):
        home = Path.home()
        s.sysdirs = {
            "STATA": str(_ROOT),
            "BASE": str(_ROOT / "ado"),
            "SITE": str(_ROOT / "site"),
            "PLUS": str(home / "ado" / "plus"),
            "PERSONAL": str(home / "ado" / "personal"),
            "OLDPLACE": str(home / "ado"),
        }
    return s.sysdirs


def adopath(s: "Session") -> list[str]:
    if not hasattr(s, "adopaths"):
        s.adopaths = list(DEFAULT_PATH)
    return s.adopaths


def resolve_entry(s: "Session", entry: str) -> Path:
    dirs = sysdir(s)
    if entry.upper() in dirs:
        return Path(dirs[entry.upper()])
    return Path(entry).expanduser()


def find_file(s: "Session", filename: str, *, all_: bool = False,
              path: list[str] | None = None, descend: bool = True) -> list[Path]:
    """Procura filename nas pastas do adopath (ou em `path`)."""
    found: list[Path] = []
    first = filename[0].lower() if filename else "_"
    for entry in (path if path is not None else adopath(s)):
        folder = resolve_entry(s, entry)
        candidates = [folder / filename]
        if descend:
            candidates.append(folder / first / filename)
        for c in candidates:
            if c.is_file() and c not in found:
                found.append(c)
                if not all_:
                    return found
    return found


def find_ado(s: "Session", name: str) -> Path | None:
    hits = find_file(s, f"{name}.ado")
    return hits[0] if hits else None
