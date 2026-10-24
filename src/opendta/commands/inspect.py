"""list and describe (manual [D] list, [D] describe), in the Stata 14 layout.

VERIFY (cases in compat/do/0101_list_describe.do): list column widths,
separators, alignment of headers and the describe header.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core.errors import StataError
from ..core.formats import format_value, parse_format
from ..core.missing import MISSING_THRESHOLD
from ..core.varlist import expand
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
    names = expand(ds, p.varlist) if p.varlist else ds.names     # duplicates are listed twice (compat 1801)
    opts = match_options(p.options, {
        "noobs": 5, "clean": 5, "separator": 3, "nolabel": 5, "noheader": 6,
        "abbreviate": 2, "table": 3, "divider": 3, "string": 3, "compress": 3,
        "fast": 4, "constant": 4, "sepby": 5, "header": 6, "nodotz": 6, "nocompress": 3,
        "mean": 4, "sum": 3, "N": 1, "labvar": 6, "display": 4,
    })
    if opts.get("noheader") and opts.get("header") is True:
        raise StataError(198, "option noheader not allowed")
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
    if isinstance(opts.get("string"), str):
        # string(#): string variables cut to # characters (VERIFY: cut, not wrapped)
        try:
            cut = int(opts["string"])
        except ValueError:
            raise StataError(198, "option string() incorrectly specified")
        if cut < 1:
            raise StataError(198, "option string() incorrectly specified")
        cells = [[c[:cut] if v.is_string else c for c, v in zip(r, vars_)] for r in cells]
    if opts.get("nodotz"):
        # .z is shown as a blank field (the manual's words; VERIFY with .z in the data)
        cells = [[("" if (c == ".z" and not v.is_string) else c) for c, v in zip(r, vars_)] for r in cells]
    # mean/sum/N: one summary row each, below the observations; the label goes in the observation-number
    # position, or in the column of labvar(); VERIFY: layout of the summary rows (compat 1701)
    labvar = str(opts["labvar"]).strip() if isinstance(opts.get("labvar"), str) else ""
    if labvar and labvar not in names:
        out.write(f"variable {labvar} not found\n", "error")
        raise StataError(198, "option labvar() incorrectly specified")
    summaries: list[tuple[str, list[str]]] = []
    for key, label in (("mean", "Mean"), ("sum", "Sum"), ("N", "N")):
        if not opts.get(key):
            continue
        pick = expand(ds, str(opts[key])) if isinstance(opts[key], str) else \
            [n for n in names if not ds.get(n).is_string]
        row = []
        for n, v in zip(names, vars_):
            if n not in pick or v.is_string:
                row.append("")
                continue
            x = np.asarray(v.data, dtype=float)[rows]
            x = x[x < MISSING_THRESHOLD]
            val = float(len(x)) if key == "N" else (float(x.sum()) if key == "sum" else
                                                    float(x.mean()) if len(x) else MISSING_THRESHOLD)
            row.append(f"{int(val)}" if key == "N" else format_value(val, v.fmt, pad=False).strip())
        if labvar:
            row[names.index(labvar)] = label
        summaries.append((label, row))
    every = int(opts["header"]) if isinstance(opts.get("header"), str) and opts["header"].strip().isdigit() else 0
    # names longer than the column are abbreviated to max(abbreviate(#), width of the values); the
    # default is ab(8) ("Rendamensal" → "Rendam~l")
    ab = int(opts["abbreviate"]) if isinstance(opts.get("abbreviate"), str) else 8
    from ..lang.functions import _abbrev
    titles = []
    widths = []
    for j, v in enumerate(vars_):
        data_w = max([len(r[j]) for r in cells + [sm[1] for sm in summaries]], default=0)
        if opts.get("nocompress") or opts.get("fast"):
            data_w = max(data_w, parse_format(v.fmt).width)     # widths from the display formats (compat 1701)
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
    obs_col = show_obs or bool(summaries and not labvar)      # the label of a summary row needs the column
    if summaries and not labvar:
        obs_w = max(obs_w, max(len(lb) for lb, _ in summaries))
    header = fmt_row(titles)
    inner = len(header) + 2
    linesize = int(float(s.settings.get("linesize", 80)))
    pad = " " * (obs_w + 2) if show_obs else "  "
    if not clean and not summaries and (opts.get("display") or (
            not opts.get("table") and len(pad) + inner + 2 > linesize)):
        _list_display(s, rows, cells, titles, widths, left, pad, obs_w, show_obs, linesize,
                      noheader=bool(opts.get("noheader")))
        return

    out.ensure_line_start()
    out.write("\n", "text")
    if clean:
        lead = " " * (obs_w + 1 + 3) if obs_col else " " * 4     # noobs: 4 spaces (compat 1801)
        if not opts.get("noheader"):
            out.write(lead + header + "  \n", "text")
        for k, (i, r) in enumerate(zip(rows, cells)):
            if every and k and k % every == 0 and not opts.get("noheader"):
                out.write("\n" + lead + header + "  \n", "text")      # VERIFY: layout of header(#) (compat 1701)
            label = f"{int(i) + 1:>{obs_w}}.   " if show_obs else lead
            out.write(label, "text")
            out.write(fmt_row(r) + "  \n", "result")
        for label, r in summaries:
            out.write((lead if labvar else label.rjust(obs_w + 1) + "   ") + fmt_row(r) + "  \n", "result")
        return

    out.write(pad + "+" + "-" * inner + "+\n", "text")
    if not opts.get("noheader"):
        out.write(pad + "| " + header + " |\n", "text")
        out.write(pad + "|" + "-" * inner + "|\n", "text")
    for k, (i, r) in enumerate(zip(rows, cells)):
        if every and k and k % every == 0 and not opts.get("noheader"):
            # a new table with its own header, after a blank line (compat 1701)
            out.write(pad + "+" + "-" * inner + "+\n\n", "text")
            out.write(pad + "+" + "-" * inner + "+\n", "text")
            out.write(pad + "| " + header + " |\n", "text")
            out.write(pad + "|" + "-" * inner + "|\n", "text")
        elif sep and k and k % sep == 0:
            out.write(pad + "|" + "-" * inner + "|\n", "text")
        elif sepby and k and any(c[i] != c[rows[k - 1]] for c in keys):
            out.write(pad + "|" + "-" * inner + "|\n", "text")
        label = f"{int(i) + 1:>{obs_w}}. " if show_obs else "  "
        out.write(label + "| ", "text")
        out.write(fmt_row(r), "result")
        out.write(" |\n", "text")
    for k, (label, r) in enumerate(summaries):
        if not k:
            out.write(pad + "|" + "-" * inner + "|\n", "text")
        lead = (" " * obs_w + "  ") if labvar else (label.rjust(obs_w + 1) + " ")
        out.write((lead if obs_col else "  ") + "| ", "text")
        out.write(fmt_row(r), "result")
        out.write(" |\n", "text")
    out.write(pad + "+" + "-" * inner + "+\n", "text")


def _list_display(s: "Session", rows, cells, titles, widths, left, pad: str, obs_w: int,
                  show_obs: bool, linesize: int, *, noheader: bool) -> None:
    """list in display format: one box per observation, with the variables in groups that fit the
    line, each group a line of names and a line of values. Also used when the table does not fit.
    VERIFY: the whole layout (blank line between boxes, lines not padded to the box width)."""
    out = s.output
    groups: list[list[int]] = [[]]
    room = linesize - len(pad) - 4
    used = -3
    for j, w in enumerate(widths):
        if groups[-1] and used + 3 + w > room:
            groups.append([])
            used = -3
        groups[-1].append(j)
        used += 3 + w

    def line(texts: list[str], idx: list[int]) -> str:
        return "| " + " | ".join(t.ljust(widths[j]) if left[j] else t.rjust(widths[j])
                                  for t, j in zip(texts, idx)) + " |"

    box = max(len(line([""] * len(g), g)) for g in groups)
    out.ensure_line_start()
    for k, (i, r) in enumerate(zip(rows, cells)):
        out.write("\n" + pad + "+" + "-" * (box - 2) + "+\n", "text")
        for gk, g in enumerate(groups):
            if gk:
                out.write(pad + "|" + "-" * (box - 2) + "|\n", "text")
            label = f"{int(i) + 1:>{obs_w}}. " if show_obs and gk == 0 else pad
            if not noheader:
                out.write(label + line([titles[j] for j in g], g) + "\n", "text")
                label = pad
            out.write(label, "text")
            out.write(line([r[j] for j in g], g) + "\n", "result")
        out.write(pad + "+" + "-" * (box - 2) + "+\n", "text")


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
    opts = match_options(p.options, {"short": 2, "simple": 2, "fullnames": 1, "numbers": 1,
                                     "replace": 3, "clear": 5, "varlist": 7})
    if opts.get("clear") and not opts.get("replace"):
        raise StataError(198, "option clear not allowed")
    if p.using is not None:
        if opts.get("fullnames") or opts.get("numbers") or opts.get("replace"):
            raise StataError(198, "option not allowed with using")   # VERIFY
        _describe_using(s, p, opts)
        return
    ds = s.data
    names = expand(ds, p.varlist) if p.varlist else ds.names     # duplicates are listed twice (compat 1801)
    if opts.get("replace"):
        _describe_replace(s, names, bool(opts.get("clear")))
        return
    out = s.output
    _LINE = _hline(s)
    # describe using a file (compat 1701; VERIFY each difference): no r(changed), r(varlist) lists all the
    # variables of the file, the header has no file name and the date stamp on the obs line, size
    # counts 4 more bytes per observation, and the lines around the table are 79 wide
    using_file = getattr(ds, "describe_using", False)
    s.r = {"N": float(ds.nobs), "k": float(ds.nvars), "width": float(ds.width())}
    if not using_file:
        s.r["changed"] = float(ds.changed)
    if opts.get("varlist"):
        s.r["varlist"] = " ".join(ds.names if using_file else names)
        s.r["sortlist"] = " ".join(ds.sortlist)

    if opts.get("simple"):
        w = max((len(n) for n in names), default=0) + 2
        out.write("".join(n.ljust(w) for n in names).rstrip() + "\n", "text")
        return

    out.ensure_line_start()
    out.write("\n", "text")
    # with a varlist, Stata 14 shows only the variables table
    only_vars = bool(p.varlist.strip())
    if not only_vars:
        src = f"Contains data from {ds.filename}" if ds.filename and not using_file else "Contains data"
        out.write(src + "\n", "text")
        label = ds.label
        ts = getattr(ds, "timestamp", "").strip() if ds.filename else ""
        out.write(f"  obs:{ds.nobs:>14,}" + (" " * 26 + (ts if using_file else label) if (ts if using_file else label) else "") + "\n", "text")
        out.write(f" vars:{ds.nvars:>14,}" + (" " * 26 + ts if ts and not using_file else "") + "\n", "text")
        dta_notes = " " * 26 + "(_dta has notes)" if ds.chars.get("_dta", {}).get("note0") else ""
        size = ds.nobs * (ds.width() + (4 if using_file else 0))
        out.write(f" size:{size:>14,}{dta_notes}\n", "text")
    if not opts.get("short"):
        if not only_vars:
            out.write(("-" * 79 if using_file else _LINE) + "\n", "text")
        numbers = bool(opts.get("numbers"))
        if numbers:
            # VERIFY: with numbers the name column is 10 wide (compat 1701)
            out.write("    variable  storage   display    value\n", "text")
            out.write("      name      type    format     label      variable label\n", "text")
        else:
            out.write("              storage   display    value\n", "text")
            out.write("variable name   type    format     label      variable label\n", "text")
        out.write(_LINE + "\n", "text")
        for k, n in enumerate(names, start=1):
            v = ds.get(n)
            shown = n if len(n) <= 15 or opts.get("fullnames") else n[:14] + "~"
            if numbers:
                shown = f"{ds.names.index(n) + 1:>5}. " + (n if len(n) <= 9 or opts.get("fullnames")
                                                             else n[:8] + "~")
            star = "*" if ds.chars.get(n, {}).get("note0") else " "
            # formatos longos aparecem cortados ("%tdnn/dd/CCYY" → "%td..")
            fmt = v.fmt if len(v.fmt) <= 9 else v.fmt[:3] + ".."   # VERIFY truncation rule
            line = f"{shown:<{17 if numbers else 16}}{v.vtype:<8}{fmt:<11}{v.value_label:<9}{star} {v.label}"
            out.write(line.rstrip() + "\n", "text")
    if only_vars:
        # VERIFY: notes footer with a varlist
        return
    # the note about variables with notes comes right below the table; describe, short has no line
    # before "Sorted by:" (compat 0104 and 0201)
    if not opts.get("short"):
        if any(ds.chars.get(n, {}).get("note0") for n in names):
            out.write(" " * 44 + "* indicated variables have notes\n", "text")
        out.write(("-" * 79 if using_file else _LINE) + "\n", "text")
    out.write("Sorted by: " + "  ".join(ds.sortlist) + "\n", "text")   # 2 spaces (Stata 14)
    if ds.changed and ds.nvars:
        out.write("     Note: Dataset has changed since last saved.\n", "text")


def _describe_using(s: "Session", p, opts: dict) -> None:
    """describe [varlist] using filename: the same report for the data in a .dta file, which is
    not loaded (the dataset in memory stays as it is)."""
    from pathlib import Path
    from ..io.dta import read_dta
    from ..lang.words import strip_outer_quotes
    path = Path(strip_outer_quotes(p.using.strip())).expanduser()
    if path.suffix == "":
        path = path.with_suffix(".dta")
    if not path.exists():
        raise StataError(601, f"file {path} not found")
    saved, saved_r = s.data, s.r
    try:
        new = read_dta(path)
        new.filename = path.name
        new.changed = False
        new.describe_using = True
        s.data = new
        rest = []
        for name in ("short", "simple", "varlist"):
            if opts.get(name):
                rest.append(name)
        cmd_describe(s, (p.varlist + (", " + " ".join(rest) if rest else "")).strip())
        saved_r = s.r
    finally:
        s.data = saved
        s.r = saved_r


def _describe_replace(s: "Session", names: list[str], clear: bool) -> None:
    """describe, replace: the data in memory become one observation per variable (position, name,
    type, isnumeric, format, vallab, varlab); nothing is shown or stored in r()."""
    from ..core.dataset import Dataset, Variable, str_type_for
    ds = s.data
    if ds.changed and ds.nvars and not clear:
        s.output.write("no; data in memory would be lost\n", "error")
        s.output.write("    Specify describe ..., replace clear to perform the command anyway.\n", "text")
        raise StataError(4, "")
    vs = [ds.get(n) for n in names]

    def text(attr: str, width_min: int = 1):
        col = [getattr(v, attr) for v in vs]
        return str_type_for(max([len(x.encode("utf-8")) for x in col] + [width_min])), col
    new = Dataset()
    new.nobs = len(vs)
    pos = np.array([float(ds.names.index(v.name) + 1) for v in vs])
    new.vars.append(Variable("position", "int", pos, label="Position of variable in dataset"))
    for name, attr, label in (("name", "name", "Variable name"), ("type", "vtype", "Storage type")):
        vtype, col = text(attr)
        new.vars.append(Variable(name, vtype, col, label=label))
    new.vars.append(Variable("isnumeric", "byte", np.array([0.0 if v.is_string else 1.0 for v in vs]),
                             label="1 if numeric, 0 if string"))
    for name, attr, label in (("format", "fmt", "Display format"), ("vallab", "value_label", "Value label name"),
                              ("varlab", "label", "Variable label")):
        vtype, col = text(attr)
        new.vars.append(Variable(name, vtype, col, label=label))
    # VERIFY: names and values of the characteristics
    # (char list shows them in reverse order of creation: d_filename, d_filedate, d_N, d_sortedby)
    new.chars["_dta"] = {"d_sortedby": " ".join(ds.sortlist), "d_N": str(ds.nobs),
                         "d_filedate": getattr(ds, "timestamp", "").strip(), "d_filename": ds.filename or ""}
    s.r = {}
    new.changed = True
    s.data = new
    s.notify_state()
