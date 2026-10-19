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


def _rule_text(fields: list[Field], seps: list[Sep]) -> str:
    out = []
    for k, f in enumerate(fields):
        sp = seps[k]
        out.append("-" * sp.before + ("+" if sp.line else "") + "-" * sp.after)
        out.append("-" * f.width)
    sp = seps[len(fields)]
    out.append("-" * sp.before + ("+" if sp.line else "") + "-" * sp.after)
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
                lines.append(" " * left + _rule_text(f, sp))
            lines.append(" " * left + _row_text([row[k] for k in cols], f, sp))
        if rules[len(rows)]:
            lines.append(" " * left + _rule_text(f, sp))
    return lines
