"""System folders (sysdir) and the search path for .ado files (adopath).

When it meets a command that is neither built in nor an already defined program,
OpenDTA looks for `name.ado` in the adopath folders (in each folder and in the
subfolder named after the first letter of the name, as in plus/o/outreg2.ado), runs the
file and calls the program it defines. This allows using third-party
packages (SSC, packages from papers) installed in these folders.

Default folders (can be changed with `sysdir set`):

    BASE      OpenDTA's own .ado files (src/opendta/ado)
    SITE      <OpenDTA folder>/site
    PERSONAL  ~/ado/personal
    PLUS      ~/ado/plus
    OLDPLACE  ~/ado

StataCorp's official .ado files never go on the path (clean-room
policy); third-party ones do.
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
    """Looks for filename in the adopath folders (or in `path`)."""
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
