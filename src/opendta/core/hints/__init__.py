"""Error explanations for people learning Stata: what went wrong in this command and why,
how to fix it and where to read more.

OpenDTA improvement over Stata: after the error message (kept as in Stata, as is the r(#)
code) comes an explanation that looks at the command and at the data in memory:

    . regress peso altura renda
    no observations
      → No observation has values for all of peso, altura and renda: renda is missing in
        all 120 observations.
        Estimation commands use only the observations with no missing value (.) in any
        variable of the model.
        How to fix: see where the data are missing with "misstable summarize peso altura
        renda" ...
        Learn more: help missing
    r(2000);

The rules live in topic modules (syntax, data, estimation, files, programming); each one
registers itself with base.rule for the messages it understands, and the first rule that
returns an explanation wins. The text follows the interface language (i18n.L) and
`set hints off` turns it off (compare.py does, to compare with Stata's logs).
"""

from __future__ import annotations

import textwrap
from typing import TYPE_CHECKING

from ...i18n import L
# the modules register their rules on import, and the first rule that answers wins: syntax comes
# last because its catch-all for "invalid syntax" must not hide the specific rules
from . import data, estimation, files, programming  # noqa: F401, I001
from . import syntax  # noqa: F401, E402
from .base import RULES, Ctx, Explanation

if TYPE_CHECKING:
    from ...session import Session
    from ..errors import StataError

__all__ = ["Explanation", "explain", "render"]


def explain(err: "StataError", s: "Session") -> Explanation | None:
    """The explanation of the error, or None when there is nothing useful to add."""
    if not err.explain or not err.message:
        return None
    command, text = err.context or ("", "")
    ctx = Ctx(err, s, command, text)
    for rx, rc, fn in RULES:
        if rc is not None and rc != err.rc:
            continue
        m = None
        if rx is not None:
            m = rx.fullmatch(err.message)
            if m is None:
                continue
        expl = fn(ctx, m)
        if expl is not None:
            from ...commands.helpcmd import find_help
            if expl.help and find_help(s, expl.help) is None:
                expl.help = ""
            if not expl.help and command and find_help(s, command) is not None:
                expl.help = command
            return expl
    return None


def render(expl: Explanation, width: int = 79) -> list[tuple[str, str]]:
    """Styled pieces (text, style) for the Results window and the logs: the explanation after
    an arrow, the reason, the fix and a link to the help page."""
    out: list[tuple[str, str]] = []
    first, indent = "  → ", "    "

    def para(text: str, lead: str, label: str = "") -> None:
        lines = textwrap.wrap(label + text, width=max(width, 40) - len(indent),
                              break_long_words=False, break_on_hyphens=False)
        for k, line in enumerate(lines):
            out.append((lead if k == 0 else indent, "hint"))
            if k == 0 and label:
                out.append((label, "hint_label"))
                line = line[len(label):]
            out.append((line + "\n", "hint"))

    para(expl.what, first)
    if expl.why:
        para(expl.why, indent)
    if expl.fix:
        para(expl.fix, indent, L("How to fix: ", "Como corrigir: "))
    if expl.help:
        out += [(indent, "hint"), (L("Learn more: ", "Saiba mais: "), "hint_label"),
                (f"help {expl.help}", "hint_link"), ("\n", "hint")]
    return out
