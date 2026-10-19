"""Draws the command-coverage heatmap shown in ROADMAP.md.

Implemented commands come from the command registry (plus the few that the
interpreter handles itself); missing ones come from the "After 1.0" section of
ROADMAP.md, where an item counts only when its first word is not an implemented
command (so options of existing commands, such as `ivregress liml`, are left
out). Each tile is an area, sized by its number of commands and coloured by the
share already implemented.

    python tools/roadmap_coverage.py [--out DIR]

writes, in DIR (default _site): roadmap-coverage.svg and roadmap-coverage.pt-BR.svg
(the images in ROADMAP.md) and index.html, an interactive page (hover, click an
area to see its commands, search). The coverage workflow runs it on every push
to main and publishes DIR on GitHub Pages.
"""

from __future__ import annotations

import argparse
import datetime
import html
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

# area of each implemented command, by the module that registers it; AREA_OF overrides it
MODULE_AREA = {
    "data": "Data management", "combine": "Data management", "collapse": "Data management",
    "datatools": "Data management", "files": "Data management", "inspect": "Data management",
    "delimited": "Data management", "infile": "Data management", "filecmd": "Data management",
    "reshape": "Data management", "label": "Data management", "codebook": "Data management",
    "putexcel": "Data management", "egen": "Data management",
    "program": "Programming", "progcmd": "Programming", "syntaxcmd": "Programming",
    "macro": "Programming", "scalar": "Programming", "matrix": "Programming",
    "adocmd": "Programming", "logcmd": "Programming", "system": "Programming",
    "helpcmd": "Programming", "uicmd": "Programming", "display": "Programming",
    "matacmd": "Programming", "preserve": "Programming",
    "summarize": "Base statistics", "tabulate": "Base statistics", "immediate": "Base statistics",
    "nptests": "Base statistics",
    "bootstrap": "Models",
    "graphcmd": "Graphics",
    "power": "Other suites",
    "epitab": "Epidemiology", "epiodds": "Epidemiology", "stdize": "Epidemiology",
    "kappa": "Epidemiology", "roc": "Epidemiology", "rocreg": "Epidemiology",
    "datautil": "Data management", "datacmds": "Data management", "datasets": "Data management",
}
AREA_OF = {
    "preserve": "Data management", "restore": "Data management",
    "bootstrap": "Base statistics", "anova": "Base statistics", "oneway": "Base statistics",
    "mean": "Base statistics", "proportion": "Base statistics", "ratio": "Base statistics",
    "total": "Base statistics", "estat": "Base statistics", "estimates": "Base statistics",
    "lincom": "Base statistics", "nlcom": "Base statistics", "test": "Base statistics",
    "testparm": "Base statistics", "predict": "Base statistics", "margins": "Base statistics",
    "lrtest": "Base statistics", "hausman": "Base statistics",
    "ivregress": "Models", "qreg": "Models",
    "xtset": "Panel data", "xtreg": "Panel data", "xtlogit": "Panel data",
    "xtpoisson": "Panel data", "xtgee": "Panel data", "xtsum": "Panel data",
    "xtdescribe": "Panel data",
    "st": "Survival analysis", "stset": "Survival analysis", "stcox": "Survival analysis",
    "streg": "Survival analysis", "streset": "Survival analysis", "sts": "Survival analysis",
    "stsum": "Survival analysis", "stdescribe": "Survival analysis",
    "tsset": "Time series", "arima": "Time series", "var": "Time series",
    "varsoc": "Time series", "vargranger": "Time series", "dfuller": "Time series",
    "corrgram": "Time series", "newey": "Time series", "irf": "Time series",
    "mixed": "Multilevel models", "xtmixed": "Multilevel models", "melogit": "Multilevel models",
    "meprobit": "Multilevel models", "mepoisson": "Multilevel models",
    "svy": "Other suites", "svyset": "Other suites", "svydescribe": "Other suites",
    "mi": "Other suites", "sem": "Other suites",
}
# handled by the interpreter, not the registry
INTERPRETER = {"by": "Programming", "bysort": "Programming", "capture": "Programming",
               "noisily": "Programming", "quietly": "Programming", "foreach": "Programming",
               "forvalues": "Programming", "while": "Programming", "if": "Programming",
               "input": "Data management"}

SHORT = {"Multivariate statistics": "Multivariate"}     # labels that do not fit their tile

PT = {
    "Data management": "Manipulação de dados", "Programming": "Programação",
    "Base statistics": "Estatística básica", "Models": "Modelos", "Graphics": "Gráficos",
    "Panel data": "Dados em painel", "Survival analysis": "Sobrevivência",
    "Time series": "Séries temporais", "Multilevel models": "Multinível",
    "Multivariate statistics": "Multivariada", "Epidemiology": "Epidemiologia",
    "Other suites": "Outras suítes",
}


def implemented() -> dict[str, str]:
    from opendta.session import Session
    Session()                   # a session loads every command module
    from opendta.commands.registry import REGISTRY
    out = dict(INTERPRETER)
    for name, spec in REGISTRY.items():
        module = spec.fn.__module__.rsplit(".", 1)[-1]
        area = AREA_OF.get(name) or MODULE_AREA.get(module)
        if area is None:
            # a new module: counted under Other suites until MODULE_AREA names its area
            print(f"warning: no area for module {module} (command {name})", file=sys.stderr)
            area = "Other suites"
        out[name] = area
    return out


# items that the "After 1.0" text names but that are options, not commands
NOT_COMMANDS = {"exactm", "exactp"}     # ties options of stcox


def missing(done: dict[str, str]) -> dict[str, set[str]]:
    """Commands listed in "After 1.0". Each bullet is a list of commands followed by clauses
    (after ";") about options of existing commands: a clause whose first item is an
    implemented command (`arima` with `sarima()`...) is about options and is skipped."""
    text = (ROOT / "ROADMAP.md").read_text(encoding="utf-8").split("## After 1.0", 1)[1]
    head = lambda item: item.split()[0].rstrip(",:")
    out: dict[str, set[str]] = {}
    for area, items in re.findall(r"^- \*\*(.+?)\*\*: (.+)$", text, re.M):
        for clause in items.split(";"):
            names = [n.strip() for item in re.findall(r"`([^`]+)`", clause) for n in item.split("/")]
            names = [n for n in names if n]
            if not names or head(names[0]) in done:
                continue
            for name in names:
                if head(name) not in done and "(" not in name and name not in NOT_COMMANDS:
                    out.setdefault(area, set()).add(name)
    return out


def squarify(items: list[tuple[str, float]], x: float, y: float, w: float, h: float):
    """Squarified treemap (Bruls et al.): yields (key, x, y, w, h)."""
    items = sorted(items, key=lambda t: -t[1])
    total = sum(v for _, v in items)
    scaled = [(k, v * w * h / total) for k, v in items]
    while scaled:
        short = min(w, h)
        row, rest = [scaled[0]], scaled[1:]

        def worst(r):
            s = sum(a for _, a in r)
            return max(max(short * short * a / (s * s), (s * s) / (short * short * a)) for _, a in r)
        while rest and worst(row + [rest[0]]) <= worst(row):
            row.append(rest.pop(0))
        s = sum(a for _, a in row)
        if w >= h:
            cw, cy = s / h, y
            for k, a in row:
                yield k, x, cy, cw, a / cw
                cy += a / cw
            x, w = x + cw, w - cw
        else:
            ch, cx = s / w, x
            for k, a in row:
                yield k, cx, y, a / ch, ch
                cx += a / ch
            y, h = y + ch, h - ch
        scaled = rest


def colour(share: float) -> str:
    """Red (nothing) through neutral grey (half) to green (all), as in a stock heatmap."""
    stops = [(0.0, (0xB8, 0x3B, 0x34)), (0.5, (0x5F, 0x63, 0x68)), (1.0, (0x2E, 0x8B, 0x4E))]
    for (p0, c0), (p1, c1) in zip(stops, stops[1:]):
        if share <= p1:
            t = (share - p0) / (p1 - p0)
            return "#" + "".join(f"{round(a + (b - a) * t):02x}" for a, b in zip(c0, c1))
    return "#2e8b4e"


def svg(rows: list[tuple[str, int, int]], lang: str) -> str:
    W, H, top, gap = 960, 520, 64, 2
    done_all = sum(d for _, d, _ in rows)
    total_all = sum(t for _, _, t in rows)
    pt = lang == "pt"
    title = (f"Comandos implementados: {done_all} de {total_all} ({100 * done_all / total_all:.0f}%)"
             if pt else
             f"Commands implemented: {done_all} of {total_all} ({100 * done_all / total_all:.0f}%)")
    sub = ("Área do bloco = número de comandos; cor = parcela implementada"
           if pt else "Tile area = number of commands; colour = share implemented")
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
           f'font-family="-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif">',
           "<style>.bg{fill:#fcfcfb}.ink{fill:#1a1a19}.muted{fill:#5f5f5b}"
           "@media (prefers-color-scheme: dark){.bg{fill:#1a1a19}.ink{fill:#f0efec}.muted{fill:#a8a8a3}}"
           "</style>",
           f'<rect class="bg" width="{W}" height="{H}" rx="8"/>',
           f'<text class="ink" x="12" y="28" font-size="18" font-weight="600">{html.escape(title)}</text>',
           f'<text class="muted" x="12" y="50" font-size="13">{html.escape(sub)}</text>']
    tiles = squarify([(a, t) for a, _, t in rows], 8, top, W - 16, H - top - 8)
    info = {a: (d, t) for a, d, t in rows}
    for area, x, y, w, h in tiles:
        d, t = info[area]
        share = d / t
        name = PT[area] if pt else SHORT.get(area, area)
        x, y, w, h = x + gap / 2, y + gap / 2, w - gap, h - gap
        out.append(f'<g><title>{html.escape(name)}: {d}/{t} ({100 * share:.0f}%)</title>'
                   f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="4" '
                   f'fill="{colour(share)}"/>')
        # font sizes from the tile, shrunk so the name fits its width (~0.6 em per character)
        big = max(14.0, min(22.0, min(w, h) / 5))
        small = max(8.0, min(big * 0.62, (w - 12) / (0.68 * len(name))))
        if w > 70 and h > 44:
            cx, cy = x + w / 2, y + h / 2
            out.append(f'<text x="{cx:.1f}" y="{cy - big * 0.15:.1f}" fill="#ffffff" font-size="{small:.1f}" '
                       f'font-weight="600" text-anchor="middle">{html.escape(name)}</text>')
            out.append(f'<text x="{cx:.1f}" y="{cy + big * 0.95:.1f}" fill="#ffffff" font-size="{big:.1f}" '
                       f'font-weight="700" text-anchor="middle">{100 * share:.0f}%</text>')
            out.append(f'<text x="{cx:.1f}" y="{cy + big * 1.75:.1f}" fill="#ffffff" fill-opacity="0.85" '
                       f'font-size="{max(10.0, big * 0.5):.1f}" text-anchor="middle">{d}/{t}</text>')
        out.append("</g>")
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="_site")
    out = Path(ap.parse_args().out)
    done = implemented()
    miss = missing(done)
    areas = sorted(set(done.values()) | set(miss))
    unknown = set(areas) - set(PT)
    if unknown:
        raise SystemExit(f"areas without a Portuguese name: {sorted(unknown)}")
    rows, data = [], []
    for a in areas:
        mine = sorted(c for c, v in done.items() if v == a)
        todo = sorted(miss.get(a, ()))
        rows.append((a, len(mine), len(mine) + len(todo)))
        data.append({"en": a, "pt": PT[a], "done": mine, "todo": todo})
    out.mkdir(parents=True, exist_ok=True)
    for lang, name in (("en", "roadmap-coverage.svg"), ("pt", "roadmap-coverage.pt-BR.svg")):
        (out / name).write_text(svg(rows, lang), encoding="utf-8")
    commit = os.environ.get("GITHUB_SHA") or subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    payload = json.dumps({"areas": data, "commit": commit or "unknown",
                          "date": datetime.date.today().isoformat()}, ensure_ascii=False)
    page = (ROOT / "tools" / "coverage_page.html").read_text(encoding="utf-8")
    (out / "index.html").write_text(page.replace("/*DATA*/", payload.replace("</", "<\\/")), encoding="utf-8")
    for a, d, t in sorted(rows, key=lambda r: -r[2]):
        print(f"{a:25} {d:4} / {t:4}  {100 * d / t:3.0f}%")


if __name__ == "__main__":
    main()
