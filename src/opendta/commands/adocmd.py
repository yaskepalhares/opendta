"""sysdir, adopath, which, findfile, discard ([P] sysdir, [R] which, [P] findfile)."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..lang.adopath import adopath, find_file, resolve_entry, sysdir
from ..lang.syntax import match_options
from ..lang.words import strip_outer_quotes
from .registry import REGISTRY, command

if TYPE_CHECKING:
    from ..session import Session


@command("sysdir")
def cmd_sysdir(s: "Session", args: str) -> None:
    words = args.split(None, 2)
    dirs = sysdir(s)
    if not words or words[0] in ("list", "l"):
        out = s.output
        for key in ("STATA", "BASE", "SITE", "PLUS", "PERSONAL", "OLDPLACE"):
            out.write(f"{key:>10}:  {dirs[key]}{'/' if not dirs[key].endswith('/') else ''}\n", "text")
        return
    if words[0] == "set" and len(words) == 3:
        key = words[1].upper()
        if key not in dirs:
            raise StataError(198, f"{words[1]} invalid codeword")   # VERIFICAR
        dirs[key] = str(Path(strip_outer_quotes(words[2])).expanduser())
        return
    raise StataError(198, "invalid syntax")


@command("adopath")
def cmd_adopath(s: "Session", args: str) -> None:
    t = args.strip()
    path = adopath(s)
    if t.startswith("++"):
        entry = strip_outer_quotes(t[2:].strip())
        if entry in path:
            path.remove(entry)
        path.insert(0, entry)
    elif t.startswith("+"):
        entry = strip_outer_quotes(t[1:].strip())
        if entry in path:
            path.remove(entry)
        path.append(entry)
    elif t.startswith("-"):
        entry = strip_outer_quotes(t[1:].strip())
        if entry.isdigit():
            k = int(entry) - 1
            if not 0 <= k < len(path):
                raise StataError(198, f"adopath element {entry} not found")
            path.pop(k)
        elif entry in path or entry.upper() in path:
            path.remove(entry if entry in path else entry.upper())
        else:
            raise StataError(198, f"{entry} not in adopath")   # VERIFICAR
    elif t:
        raise StataError(198, "invalid syntax")
    out = s.output
    for k, entry in enumerate(path, start=1):
        shown = f'({entry})' if entry.upper() in sysdir(s) else ""
        target = resolve_entry(s, entry) if entry != "." else Path(".")
        out.write(f"  [{k}]  {shown:<11} \"{str(target).rstrip('/') + '/' if entry != '.' else '.'}\"\n",
                  "text")


@command("which")
def cmd_which(s: "Session", args: str) -> None:
    head, comma, opts = args.partition(",")
    name = head.strip()
    if not name:
        raise StataError(100, "name required")   # VERIFICAR
    o = match_options(opts, {"all": 3}) if comma else {}
    out = s.output
    if name in REGISTRY and not o.get("all"):
        out.write("built-in command:  " + name + "\n", "text")
        return
    hits = find_file(s, name if "." in name else f"{name}.ado", all_=bool(o.get("all")))
    if not hits:
        if name in REGISTRY:
            out.write("built-in command:  " + name + "\n", "text")
            return
        raise StataError(111, f"command {name} not found as either built-in or ado-file")
    if name in REGISTRY:
        out.write("built-in command:  " + name + "\n", "text")
    for h in hits:
        out.write(f"{h}\n", "text")
        try:
            for line in h.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("*!"):
                    out.write(line + "\n", "text")
        except OSError:
            pass


@command("findfile")
def cmd_findfile(s: "Session", args: str) -> None:
    head, comma, opts = args.partition(",")
    fname = strip_outer_quotes(head.strip())
    if not fname:
        raise StataError(198, "invalid syntax")
    o = match_options(opts, {"path": 4, "nodescend": 6, "all": 3}) if comma else {}
    path = None
    if o.get("path"):
        path = [p.strip() for p in strip_outer_quotes(str(o["path"])).split(";") if p.strip()]
    hits = find_file(s, fname, all_=bool(o.get("all")), path=path, descend=not o.get("nodescend"))
    if not hits:
        raise StataError(601, f"file \"{fname}\" not found")
    s.r = {"fn": " ".join(f'"{h}"' for h in hits) if o.get("all") else str(hits[0])}
    for h in hits:
        s.output.write(f"{h}\n", "text")


@command("discard")
def cmd_discard(s: "Session", args: str) -> None:
    """Esquece os programas carregados de .ado (são relidos no próximo uso)."""
    for name in [n for n, p in s.programs.items() if p.source.endswith(".ado")]:
        del s.programs[name]



