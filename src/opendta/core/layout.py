"""Text layout of tables: matlist and tabdisp ([P] matlist, [R] tabdisp).

Plain functions: they take strings and numbers and return lines of text. The commands
(commands/matlist.py, commands/tabdisp.py) read the matrix or the variables, format the cells and
write the lines.

A table is a list of rows of cells, a `Field` (width and alignment) for each column and a `Sep`
for each boundary between columns, the borders included: a separator may carry a vertical line
and spaces before and after it. `rules` says where a horizontal line is drawn: before the first
row, between rows and after the last one.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Sep:
    line: bool = False
    before: int = 0
    after: int = 0

    @property
    def width(self) -> int:
        return self.before + (1 if self.line else 0) + self.after


@dataclass
class Field:
    width: int
    align: str = "R"          # L, R or C


def align(text: str, width: int, how: str) -> str:
    text = text[:width]
    if how == "L":
        return text.ljust(width)
    if how == "C":
        pad = width - len(text)
        return " " * (pad // 2) + text + " " * (pad - pad // 2)
    return text.rjust(width)


def _row_text(cells: list[str], fields: list[Field], seps: list[Sep]) -> str:
    out = []
    for k, f in enumerate(fields):
        sp = seps[k]
        out.append(" " * sp.before + ("|" if sp.line else "") + " " * sp.after)
        out.append(align(cells[k], f.width, f.align))
    sp = seps[len(fields)]
    out.append(" " * sp.before + ("|" if sp.line else "") + " " * sp.after)
    return "".join(out).rstrip() if not seps[len(fields)].line else "".join(out)


def _rule_text(fields: list[Field], seps: list[Sep], outer: bool = False) -> str:
    """A horizontal line. The lines at the top and the bottom of the table have no '+' where the
    vertical lines of the columns would cross them (only at the corners; compat 1801)."""
    last = len(fields)
    out = []
    for k, f in enumerate(fields):
        sp = seps[k]
        cross = "+" if sp.line and (not outer or k == 0) else "-" if sp.line else ""
        out.append((" " if k == 0 else "-") * sp.before + cross + "-" * sp.after)
        out.append("-" * f.width)
    sp = seps[last]
    cross = "+" if sp.line else ""
    out.append("-" * sp.before + cross + "-" * sp.after)
    return "".join(out)


def render_table(rows: list[list[str]], fields: list[Field], seps: list[Sep],
                 rules: list[bool], *, left: int = 0, limit: int | None = None,
                 repeat: int = 1) -> list[str]:
    """The rows as lines. `rules[i]` draws a horizontal line before row i (the last one, after the
    last row). When the table is wider than `limit`, the columns after the first `repeat` ones are
    split into panels that repeat those first columns."""
    n = len(fields)
    chunks: list[list[int]] = []
    if limit is None:
        chunks = [list(range(n))]
    else:
        cur: list[int] = list(range(repeat))
        width = left + seps[0].width + sum(fields[k].width + seps[k + 1].width for k in cur)
        for k in range(repeat, n):
            w = fields[k].width + seps[k + 1].width
            if len(cur) > repeat and width + w > limit:
                chunks.append(cur)
                cur = list(range(repeat))
                width = left + seps[0].width + sum(fields[j].width + seps[j + 1].width for j in cur)
            cur.append(k)
            width += w
        chunks.append(cur)
    lines: list[str] = []
    for cols in chunks:
        f = [fields[k] for k in cols]
        # separators: the left border, the ones after each column and, at the end of a panel, the
        # right border of the table
        sp = [seps[0]] + [seps[k + 1] if k != cols[-1] or cols[-1] == n - 1 else seps[n]
                          for k in cols]
        for i, row in enumerate(rows):
            if rules[i]:
                lines.append(" " * left + _rule_text(f, sp, outer=i == 0))
            lines.append(" " * left + _row_text([row[k] for k in cols], f, sp))
        if rules[len(rows)]:
            lines.append(" " * left + _rule_text(f, sp, outer=True))
    return lines


# ---------------------------------------------------------------------------
# tabdisp
# ---------------------------------------------------------------------------

def wrap_words(text: str, width: int) -> list[str]:
    """Greedy filling of lines of at most `width` characters; a word longer than the line is cut
    (the stub title of tabdisp: "Party Affiliation and Sex" in 10 columns)."""
    lines: list[str] = []
    cur = ""
    for word in text.split():
        while len(word) > width:
            if cur:
                lines.append(cur)
                cur = ""
            lines.append(word[:width])
            word = word[width:]
        if not cur:
            cur = word
        elif len(cur) + 1 + len(word) <= width:
            cur += " " + word
        else:
            lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines or [""]


def center_text(text: str, width: int) -> str:
    """Centered with the extra blank on the left (VERIFY with an odd difference)."""
    pad = max(0, width - len(text))
    left = pad - pad // 2
    return " " * left + text + " " * (pad - left)


def _dashed(label: str, width: int) -> str:
    """'------- Democrat -------': the label between blanks, centered in dashes."""
    inner = f" {label} "
    if len(inner) >= width:
        return label[:width].ljust(width)
    dashes = width - len(inner)
    left = dashes // 2            # VERIFY: which side gets the extra dash
    return "-" * left + inner + "-" * (dashes - left)


def render_tabdisp(*, stub_title: str, body: list, heads: list[str], groups: list[str] | None,
                   group_title: str, super_title: str, width: int, csep: int, scsep: int,
                   stub: int, center: bool, left: bool) -> list[str]:
    """The lines of a tabdisp table.

    `body` has, in order, ("group", title) for the heading of a superrow combination and
    ("row", label, lines) for a row, where `lines` is a list of lists with the text of each column
    (all supercolumns, left to right); several lines are the several cell variables.
    `heads` are the column headings of one supercolumn; `groups` the labels of the supercolumns
    (None without a supercolumn variable); `group_title` is the title over the columns (the
    column variable) and `super_title` the one over the supercolumns."""
    per = len(heads)
    ng = len(groups) if groups else 1
    ncols = per * ng

    def gap(k: int) -> int:
        if k == 0:
            return 1
        return csep + (scsep if k % per == 0 else 0)

    def cells(texts: list[str], heading: bool = False) -> str:
        out = []
        for k, t in enumerate(texts):
            out.append(" " * gap(k) + (center_text(t, width) if center else t.rjust(width)))
        return "".join(out)

    block = len(cells([""] * ncols))
    # heading rows to the right of the stub
    header: list[str] = []
    if groups:
        header.append(center_text(super_title, block).rstrip())
        line, pos = "", 0
        for gi, label in enumerate(groups):
            start = pos + gap(gi * per)
            span = per * width + (per - 1) * csep
            line = line.ljust(start) + _dashed(label, span)
            pos = start + span
        header.append(line)
    elif group_title:
        header.append(center_text(group_title, block).rstrip())
    heads_text = [h for _ in range(ng) for h in heads]
    header.append(cells(heads_text))
    stub_lines = wrap_words(stub_title, stub)
    single = len(stub_lines) == 1
    n = max(len(header), len(stub_lines))
    header = [""] * (n - len(header)) + header
    stub_lines = [""] * (n - len(stub_lines)) + stub_lines
    top = "-" * (stub + 2 + block)
    rule = "-" * (stub + 1) + "+" + "-" * block
    lines = [top]
    for k in range(n):
        st = stub_lines[k]
        st = st.rjust(stub) if single else st.ljust(stub)
        lines.append(f"{st} |{header[k]}")
    lines.append(rule)
    first = True
    for item in body:
        if item[0] == "blank":
            lines.append(" " * (stub + 1) + "|")
            continue
        if item[0] == "group":
            if not first:
                lines.append(rule)
            lines.append(f"{item[1]:<{stub}} |")
            first = False
            continue
        _, label, rows = item
        for j, texts in enumerate(rows):
            lab = label if j == 0 else ""
            lab = lab.ljust(stub) if left else lab.rjust(stub)
            lines.append(f"{lab} |{cells(texts)}")
    lines.append(top)
    return lines
