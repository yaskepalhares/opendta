"""list and describe (manual [D] list, [D] describe), in the Stata 14 layout.

VERIFY (cases in compat/do/0101_list_describe.do): list column widths,
separators, alignment of headers and the describe header.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from ..core.formats import format_value, parse_format
from ..core.varlist import expand, unique
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .registry import command

if TYPE_CHECKING:
    from ..core.dataset import Dataset, Variable
    from ..session import Session


def cell_text(ds: "Dataset", var: "Variable", i: int, *, use_labels: bool = True) -> str:
    v = var.value(i)
    if var.is_string:
        return str(v)
    x = float(v)
    if use_labels and var.value_label and x < 1e300 and x == int(x):
        lab = ds.value_labels.get(var.value_label, {})
        if int(x) in lab:
            return lab[int(x)]
    return format_value(x, var.fmt, pad=False).strip()


def _left_aligned(var: "Variable", use_labels: bool) -> bool:
    """Strings follow the justification of the format (%-18s left, %18s
    right); numbers, with or without a value label, are right-aligned."""
    if var.is_string:
        return parse_format(var.fmt).left
    return False


@command("list", "l", byable=True)
def cmd_list(s: "Session", args: str) -> None:
    p = parse_standard(args)
    ds = s.data
    names = unique(expand(ds, p.varlist)) if p.varlist else ds.names
    opts = match_options(p.options, {
        "noobs": 5, "clean": 5, "separator": 3, "nolabel": 5, "noheader": 6,
        "abbreviate": 2, "table": 3, "divider": 3, "string": 3, "compress": 3,
        "fast": 4, "constant": 4, "sepby": 5,
    })
    mask = touse(s, p)
    rows = np.flatnonzero(mask)
    if not names:
        return
    out = s.output
    use_labels = not opts.get("nolabel")
    vars_ = [ds.get(n) for n in names]
    show_obs = not opts.get("noobs")
    clean = bool(opts.get("clean"))
    sep = int(opts["separator"]) if isinstance(opts.get("separator"), str) else (0 if clean else 5)
    divider = bool(opts.get("divider"))
    sepby = expand(ds, str(opts["sepby"])) if isinstance(opts.get("sepby"), str) else []
    if sepby:
        sep = 0
        keys = [ds.get(n).raw for n in sepby]

    cells = [[cell_text(ds, v, int(i), use_labels=use_labels) for v in vars_] for i in rows]
    # names longer than the column are abbreviated to max(abbreviate(#), width of the values); the
    # default is ab(8) ("Rendamensal" → "Rendam~l")
    ab = int(opts["abbreviate"]) if isinstance(opts.get("abbreviate"), str) else 8
    from ..lang.functions import _abbrev
    titles = []
    widths = []
    for j, v in enumerate(vars_):
        data_w = max([len(r[j]) for r in cells], default=0)
        title = v.name if len(v.name) <= max(ab, data_w) else _abbrev(v.name, float(max(ab, data_w)))
        titles.append(title)
        widths.append(max([len(title), data_w]))
    left = [_left_aligned(v, use_labels) for v in vars_]

    def fmt_row(texts: list[str]) -> str:
        parts = []
        for t, w, lft in zip(texts, widths, left):
            parts.append(t.ljust(w) if lft else t.rjust(w))
        return (" | " if divider else "   ").join(parts)

    # layout observed in the Stata logs: 3 spaces between columns, observation label "  1."
    # followed by " | " (table) or by 3 spaces (clean)
    obs_w = max(3, len(str(int(rows[-1]) + 1))) if len(rows) else 3
    header = fmt_row(titles)
    inner = len(header) + 2

    out.ensure_line_start()
    out.write("\n", "text")
    if clean:
        lead = " " * (obs_w + 1 + 3) if show_obs else ""
        if not opts.get("noheader"):
            out.write(lead + header + "  \n", "text")
        for i, r in zip(rows, cells):
            label = f"{int(i) + 1:>{obs_w}}.   " if show_obs else ""
            out.write(label, "text")
            out.write(fmt_row(r) + "  \n", "result")
        return

    pad = " " * (obs_w + 2) if show_obs else "  "
    out.write(pad + "+" + "-" * inner + "+\n", "text")
    if not opts.get("noheader"):
        out.write(pad + "| " + header + " |\n", "text")
        out.write(pad + "|" + "-" * inner + "|\n", "text")
    for k, (i, r) in enumerate(zip(rows, cells)):
        if sep and k and k % sep == 0:
            out.write(pad + "|" + "-" * inner + "|\n", "text")
        elif sepby and k and any(c[i] != c[rows[k - 1]] for c in keys):
            out.write(pad + "|" + "-" * inner + "|\n", "text")
        label = f"{int(i) + 1:>{obs_w}}. " if show_obs else "  "
        out.write(label + "| ", "text")
        out.write(fmt_row(r), "result")
        out.write(" |\n", "text")
    out.write(pad + "+" + "-" * inner + "+\n", "text")


# ---------------------------------------------------------------------------
# describe
# ---------------------------------------------------------------------------

def _hline(s: "Session") -> str:
    """Full describe line: takes c(linesize) (255 dashes with set linesize
    255, checked in compat; VERIFY with the default 80)."""
    return "-" * int(s.settings.get("linesize", 80))


@command("describe", "d")
def cmd_describe(s: "Session", args: str) -> None:
    p = parse_standard(args)
    ds = s.data
    opts = match_options(p.options, {"short": 2, "simple": 2, "fullnames": 1, "numbers": 1})
    names = unique(expand(ds, p.varlist)) if p.varlist else ds.names
    out = s.output
    _LINE = _hline(s)
    s.r = {"N": float(ds.nobs), "k": float(ds.nvars), "width": float(ds.width()),
           "changed": float(ds.changed)}

    if opts.get("simple"):
        out.write("  ".join(names) + "\n", "text")
        return

    out.ensure_line_start()
    out.write("\n", "text")
    # with a varlist, Stata 14 shows only the variables table
    only_vars = bool(p.varlist.strip())
    if not only_vars:
        src = f"Contains data from {ds.filename}" if ds.filename else "Contains data"
        out.write(src + "\n", "text")
        label = ds.label
        out.write(f"  obs:{ds.nobs:>14,}" + (" " * 26 + label if label else "") + "\n", "text")
        ts = getattr(ds, "timestamp", "") if ds.filename else ""
        out.write(f" vars:{ds.nvars:>14,}" + (" " * 26 + ts.strip() if ts else "") + "\n", "text")
        dta_notes = " " * 26 + "(_dta has notes)" if ds.chars.get("_dta", {}).get("note0") else ""
        out.write(f" size:{ds.width() * ds.nobs:>14,}{dta_notes}\n", "text")
    if not opts.get("short"):
        if not only_vars:
            out.write(_LINE + "\n", "text")
        out.write("              storage   display    value\n", "text")
        out.write("variable name   type    format     label      variable label\n", "text")
        out.write(_LINE + "\n", "text")
        for n in names:
            v = ds.get(n)
            shown = n if len(n) <= 15 or opts.get("fullnames") else n[:14] + "~"
            star = "*" if ds.chars.get(n, {}).get("note0") else " "
            # formatos longos aparecem cortados ("%tdnn/dd/CCYY" → "%td..")
            fmt = v.fmt if len(v.fmt) <= 9 else v.fmt[:3] + ".."   # VERIFY truncation rule
            line = f"{shown:<16}{v.vtype:<8}{fmt:<11}{v.value_label:<9}{star} {v.label}"
            out.write(line.rstrip() + "\n", "text")
    if only_vars:
        # VERIFY: notes footer with a varlist
        return
    # the note about variables with notes comes right below the table; describe, short has no line
    # before "Sorted by:" (compat 0104 and 0201)
    if not opts.get("short"):
        if any(ds.chars.get(n, {}).get("note0") for n in names):
            out.write(" " * 44 + "* indicated variables have notes\n", "text")
        out.write(_LINE + "\n", "text")
    out.write("Sorted by: " + "  ".join(ds.sortlist) + "\n", "text")   # 2 spaces (Stata 14)
    if ds.changed and ds.nvars:
        out.write("     Note: Dataset has changed since last saved.\n", "text")
