"""Solid table lines for the Results window.

Commands write their tables as text, as Stata's logs do: dashes for the horizontal lines, "|"
for the vertical ones and "+" where they cross. That text stays as it is (copying and logs give
the same characters as Stata); the Results window draws continuous lines over it. This module
decides, line by line, which character cells are part of a line and which sides of the cell the
line touches; the drawing is in panes.ResultsView.
"""

from __future__ import annotations

LEFT, RIGHT, UP, DOWN = 1, 2, 4, 8

_RULE_CHARS = frozenset("-+|")


def is_rule(line: str) -> bool:
    """A horizontal rule: only "-", "+" and "|" (at least three characters, one of them "-")."""
    s = line.strip()
    return len(s) >= 3 and "-" in s and _RULE_CHARS.issuperset(s)


def _continues_up_or_down(line: str | None, col: int) -> bool:
    """Does this neighbor line have a vertical line at `col`?"""
    if not line or col >= len(line):
        return False
    ch = line[col]
    return ch == "|" or (ch == "+" and is_rule(line))


def table_cells(prev: str | None, line: str, nxt: str | None) -> list[tuple[int, int]]:
    """The cells of `line` drawn as lines: (column, sides). Sides is a combination of LEFT, RIGHT,
    UP and DOWN: the parts of the cell, from its center, that carry the line. `prev` and `nxt` are
    the lines around it (None at the ends)."""
    cells = []
    if is_rule(line):
        for col, ch in enumerate(line):
            if ch == "-":
                cells.append((col, LEFT | RIGHT))
            elif ch in "+|":
                sides = 0
                if _continues_up_or_down(prev, col):
                    sides |= UP
                if _continues_up_or_down(nxt, col):
                    sides |= DOWN
                if col > 0 and line[col - 1] in "-+":
                    sides |= LEFT
                if col + 1 < len(line) and line[col + 1] in "-+":
                    sides |= RIGHT
                if sides:
                    cells.append((col, sides))
    else:
        # a "|" is part of a table when the line above or below has a vertical line in the same column
        for col, ch in enumerate(line):
            if ch == "|" and (_continues_up_or_down(prev, col) or _continues_up_or_down(nxt, col)):
                cells.append((col, UP | DOWN))
    return cells
