"""Splits the text of a do-file into logical command lines.

Rules implemented (manual [U] 16 and [P] comments, #delimit):

* a line whose first non-blank character is `*` is a comment;
* `//` preceded by a space (or at the start of the line) comments to the end of the line;
* `///` preceded by a space joins the next line (continuation);
* `/* ... */` may span lines and may be nested;
* `#delimit ;` switches to `;` as the end of command; `#delimit cr` goes back to
  end of line (the abbreviation `#d` is also accepted). The directive stays in the
  list as the line "#delimit ;" or "#delimit cr", for the echo;
* none of this applies inside "..." or `"..."' strings (compound quotes).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LogicalLine:
    text: str
    lineno: int                  # line (1-based) where the command starts
    raw: tuple[str, ...] = ()    # original physical lines, for the echo
    kind: str = "cmd"            # "cmd", "comment" or "delimit"

    @property
    def echo_lines(self) -> tuple[str, ...]:
        return self.raw or (self.text,)


def _is_delimit(cmd: str) -> str | None:
    """If `cmd` is a #delimit directive, returns 'cr' or ';'."""
    s = cmd.strip()
    if not s.startswith("#d"):
        return None
    head, _, rest = s.partition(" ")
    if not "#delimit".startswith(head) or len(head) < 2:
        return None
    arg = rest.strip()
    if arg == ";":
        return ";"
    if arg == "cr":
        return "cr"
    return None


def split_commands(source: str) -> list[LogicalLine]:
    """Whole-line comments (`*` and `//`) also go into the list, with
    kind="comment": Stata echoes them when running a do-file."""
    physical = source.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    out: list[LogicalLine] = []
    mode = "cr"
    buf: list[str] = []
    start_line = 0
    line = 1
    i = 0
    n = len(source)
    block = 0        # /* */ depth
    in_str = False   # inside "..."
    compound = 0     # `" "' depth

    def at_word_start(pos: int) -> bool:
        return pos == 0 or source[pos - 1] in " \t\n\r"

    def raw_lines(first: int, last: int) -> tuple[str, ...]:
        return tuple(ln.rstrip() for ln in physical[first - 1:last])

    def flush() -> None:
        nonlocal buf, start_line, mode
        text = "".join(buf).strip()
        buf = []
        if not text:
            return
        raw = raw_lines(start_line, line)
        if text.startswith("*"):
            out.append(LogicalLine(text, start_line, raw, "comment"))
            return
        new_mode = _is_delimit(text)
        if new_mode is not None:
            mode = new_mode
            out.append(LogicalLine("#delimit " + new_mode, start_line, raw, "delimit"))
            return
        out.append(LogicalLine(text, start_line, raw))

    def note_start() -> None:
        nonlocal start_line
        if not "".join(buf).strip():
            start_line = line

    while i < n:
        c = source[i]

        if block:
            if source.startswith("*/", i):
                block -= 1
                i += 2
                if block == 0:
                    buf.append(" ")
                continue
            if source.startswith("/*", i):
                block += 1
                i += 2
                continue
            if c == "\n":
                line += 1
            i += 1
            continue

        if in_str:
            if c == "\n":
                in_str = False  # unclosed string: the command reports the error
                continue
            buf.append(c)
            if c == '"':
                in_str = False
            i += 1
            continue

        if compound:
            if source.startswith('`"', i):
                compound += 1
                buf.append('`"')
                i += 2
                continue
            if source.startswith("\"'", i):
                compound -= 1
                buf.append("\"'")
                i += 2
                continue
            if c == "\n":
                compound = 0
                continue
            buf.append(c)
            i += 1
            continue

        # outside strings and comments
        if c == "\r":
            i += 1
            continue
        if source.startswith("/*", i):
            block = 1
            i += 2
            continue
        if source.startswith("///", i) and at_word_start(i):
            j = source.find("\n", i)
            if j == -1:
                i = n
            else:
                i = j + 1
                line += 1
            buf.append(" ")
            continue
        if source.startswith("//", i) and at_word_start(i):
            j = source.find("\n", i)
            if not "".join(buf).strip():
                # whole comment line: only for the echo
                out.append(LogicalLine("", line, raw_lines(line, line), "comment"))
            i = n if j == -1 else j
            continue
        if c == "\n":
            current = "".join(buf).strip()
            if mode == "cr" or _is_delimit(current) is not None:
                flush()
            elif current.startswith("*"):
                # * comment in ; mode runs until the next ;
                buf.append(" ")
            else:
                buf.append(" ")
            line += 1
            i += 1
            continue
        if c == ";" and mode == ";":
            flush()
            i += 1
            continue
        note_start()
        if source.startswith('`"', i):
            compound = 1
            buf.append('`"')
            i += 2
            continue
        if c == '"':
            in_str = True
        buf.append(c)
        i += 1

    flush()
    return out
