"""Checks before risky commands typed in the Command window (Settings > General > Safety).

check() looks at one command line before it runs and says whether it would throw away a lot of work:
`clear` with unsaved changes, `drop` or `keep` removing many observations or variables, and
`save, replace` over a file much larger than the data in memory. The limits are conservative: both
the absolute number and the share must be exceeded, so small data sets and small trims never ask.

Do-files are not checked (they run as written); neither is anything after a `by` prefix. The command
itself still reports its own errors: when anything here fails, check() answers "no risk".
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

from .i18n import L

if TYPE_CHECKING:
    from .session import Session

MIN_OBS_REMOVED = 10_000        # observations, and...
MIN_OBS_SHARE = 0.25            # ...a quarter of the data
MIN_VARS_REMOVED = 10           # variables, and...
MIN_VARS_SHARE = 0.5            # ...half of them
MIN_FILE_BYTES = 5 * 1024 * 1024    # an existing file of at least 5 MB, replaced by data...
MAX_SIZE_RATIO = 0.5                # ...under half its size

_PREFIX = re.compile(r"\s*(?:quietly|qui|noisily|capture|cap)\s*:?\s+", re.I)
_WORD = re.compile(r"\s*([A-Za-z_]\w*)(.*)$", re.S)


COMMANDS = ("clear", "drop", "keep", "save")      # each one has its own switch in Settings


@dataclass(frozen=True)
class Risk:
    command: str        # one of COMMANDS
    kind: str           # clear | drop_obs | drop_vars | save
    message: str        # what the command is about to do, in the language of the interface


def check(s: "Session", line: str, skip: frozenset[str] = frozenset()) -> Risk | None:
    """The risk of running the line, or None. Commands in `skip` (switched off) are not even examined."""
    try:
        return _check(s, line, skip)
    except Exception:  # noqa: BLE001 - the command reports its own errors
        return None


def _check(s: "Session", line: str, skip: frozenset[str]) -> Risk | None:
    from .commands.registry import lookup

    text = line.strip()
    while True:
        m = _PREFIX.match(text)
        if not m:
            break
        text = text[m.end():]
    m = _WORD.match(text)
    if not m:
        return None
    spec = lookup(m.group(1))
    if spec is None or spec.name not in COMMANDS or spec.name in skip:
        return None
    rest = s.expand(m.group(2))
    ds = s.data
    if spec.name == "clear":
        if rest.strip() in ("", "all", "*") and ds.changed and ds.nobs > 0:
            return Risk("clear", "clear", L(
                f"The data in memory have unsaved changes ({ds.nobs:,} observations, {ds.nvars:,} variables). "
                "clear will discard them.",
                f"Os dados na memória têm alterações não salvas ({ds.nobs:,} observações, {ds.nvars:,} variáveis). "
                "O clear vai descartá-las."))
        return None
    if spec.name == "save":
        return _check_save(s, rest)
    return _check_drop_keep(s, rest, keep=spec.name == "keep")


def _check_drop_keep(s: "Session", rest: str, keep: bool) -> Risk | None:
    from .commands._util import touse
    from .core.varlist import expand, unique
    from .lang.syntax import parse_standard

    p = parse_standard(rest.strip())
    ds = s.data
    n, k = ds.nobs, ds.nvars
    verb = "keep" if keep else "drop"
    if p.varlist and (p.if_ or p.in_):
        return None
    if p.varlist:
        names = unique(expand(ds, p.varlist))
        removed = k - len(names) if keep else len(names)
        if removed >= MIN_VARS_REMOVED and removed >= MIN_VARS_SHARE * k:
            return Risk(verb, "drop_vars", L(
                f"{verb} will remove {removed:,} of the {k:,} variables.",
                f"O {verb} vai remover {removed:,} das {k:,} variáveis."))
        return None
    if not (p.if_ or p.in_) or n == 0:
        return None
    selected = int(touse(s, p).sum())
    removed = n - selected if keep else selected
    if removed >= MIN_OBS_REMOVED and removed >= MIN_OBS_SHARE * n:
        return Risk(verb, "drop_obs", L(
            f"{verb} will delete {removed:,} of the {n:,} observations ({removed / n:.0%}).",
            f"O {verb} vai apagar {removed:,} das {n:,} observações ({removed / n:.0%})."))
    return None


def _check_save(s: "Session", rest: str) -> Risk | None:
    from .commands.files import dta_path
    from .io.dta import _widths
    from .lang.syntax import match_options

    head, comma, opts = rest.partition(",")
    if not comma or not match_options(opts, {"replace": 3, "emptyok": 5}).get("replace"):
        return None
    ds = s.data
    path = dta_path(head) if head.strip() else (ds.fullpath or ds.filename)
    if not path or not os.path.isfile(path):
        return None
    old = os.path.getsize(path)
    new = sum(_widths(ds)) * ds.nobs        # about: the records only
    if old >= MIN_FILE_BYTES and new < MAX_SIZE_RATIO * old:
        return Risk("save", "save", L(
            f"save, replace will overwrite {os.path.basename(str(path))} ({old / 1048576:,.1f} MB) with "
            f"data of about {new / 1048576:,.1f} MB.",
            f"O save, replace vai sobrescrever {os.path.basename(str(path))} ({old / 1048576:,.1f} MB) com "
            f"dados de cerca de {new / 1048576:,.1f} MB."))
    return None
