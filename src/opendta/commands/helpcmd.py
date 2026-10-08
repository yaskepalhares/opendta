"""help, view e search.

help procura o tópico primeiro nas pastas do adopath (arquivos .sthlp e
.hlp de pacotes de terceiros) e depois nas páginas próprias do OpenDTA
(src/opendta/help). Com a interface gráfica, a página abre no Viewer; sem
ela, o texto sai em Results.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.smcl import render_text
from ..lang.adopath import find_file
from ..lang.words import strip_outer_quotes
from .registry import REGISTRY, command

if TYPE_CHECKING:
    from ..session import Session

HELP_DIR = Path(__file__).resolve().parents[1] / "help"

ALIASES = {
    "d": "describe", "de": "describe", "des": "describe", "g": "generate", "gen": "generate",
    "replace": "generate", "l": "list", "li": "list", "saveold": "save", "export": "import",
    "insheet": "import", "outsheet": "import", "import_delimited": "import",
    "import_excel": "import", "export_excel": "import", "export_delimited": "import",
    "infix": "infile", "restore": "preserve", "tempvar": "preserve", "tempfile": "preserve",
    "tempname": "preserve", "return": "program", "ereturn": "program", "sreturn": "program",
    "args": "syntax", "gettoken": "syntax", "marksample": "syntax", "markout": "syntax",
    "mat": "matrix", "cmdlog": "log", "set_hints": "hints", "hint": "hints",
}


def find_help(s: "Session", topic: str) -> Path | None:
    t = "_".join(topic.split()) or "opendta"
    t = ALIASES.get(t, t)
    for ext in (".sthlp", ".hlp"):
        hits = find_file(s, t + ext)
        if hits:
            return hits[0]
    own = HELP_DIR / f"{t}.sthlp"
    if own.exists():
        return own
    full = REGISTRY.get(t)
    if full is not None and (HELP_DIR / f"{full.name}.sthlp").exists():
        return HELP_DIR / f"{full.name}.sthlp"
    return None


def available_topics() -> list[str]:
    return sorted(p.stem for p in HELP_DIR.glob("*.sthlp"))


def _read(path: Path) -> str:
    raw = path.read_bytes()
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


@command("help", "h")
def cmd_help(s: "Session", args: str) -> None:
    topic = args.strip().split(",")[0].strip()
    path = find_help(s, topic)
    if path is None:
        # VERIFICAR: mensagem do Stata 14 para tópico inexistente
        raise StataError(111, f"help for {topic} not found",
                         hint=f'"{topic}" has no help page yet.\nHint: type "help" to see the '
                         "topics available.")
    hook = s.ui_hooks.get("help")
    if hook is not None:
        hook(path, topic or "opendta")
        return
    s.output.write(render_text(_read(path), int(float(s.settings.get("linesize", 80))) - 1), "text")


@command("view")
def cmd_view(s: "Session", args: str) -> None:
    t = args.strip()
    if t.startswith("help "):
        cmd_help(s, t[5:])
        return
    if t.startswith("file "):
        t = t[5:]
    path = Path(strip_outer_quotes(t)).expanduser()
    if not path.exists():
        raise StataError(601, f"file {path} not found")
    hook = s.ui_hooks.get("view")
    if hook is not None:
        hook(path, path.name)
        return
    text = _read(path)
    if path.suffix in (".smcl", ".sthlp", ".hlp") or text.lstrip().startswith("{smcl}"):
        text = render_text(text)
    s.output.write(text if text.endswith("\n") else text + "\n", "text")


@command("search")
def cmd_search(s: "Session", args: str) -> None:
    """Procura palavras nos títulos e no texto das páginas de help."""
    words = [w.lower() for w in args.split() if not w.startswith(",")]
    if not words:
        raise StataError(198, "nothing to search for")   # VERIFICAR
    out = s.output
    hits = []
    for path in sorted(HELP_DIR.glob("*.sthlp")):
        text = render_text(_read(path)).lower()
        if all(w in text for w in words):
            title = next((ln.strip() for ln in render_text(_read(path)).splitlines()
                          if "--" in ln), path.stem)
            hits.append((path.stem, title))
    if not hits:
        out.write(f"no entries found for search on \"{' '.join(words)}\"\n", "text")
        return
    for topic, title in hits:
        out.write(f"  [{topic}]  ", "result")
        out.write(title + "\n", "text")
