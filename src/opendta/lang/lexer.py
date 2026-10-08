"""Divide o texto de um do-file em linhas lógicas de comando.

Regras implementadas (manual [U] 16 e [P] comments, #delimit):

* linha cujo primeiro caractere não branco é `*` é comentário;
* `//` precedido de espaço (ou no início da linha) comenta até o fim da linha;
* `///` precedido de espaço junta a linha seguinte (continuação);
* `/* ... */` pode atravessar linhas e pode ser aninhado;
* `#delimit ;` passa a usar `;` como fim de comando; `#delimit cr` volta ao
  fim de linha (a abreviação `#d` também é aceita). A diretiva continua na
  lista como a linha "#delimit ;" ou "#delimit cr", para o eco;
* nada disso vale dentro de strings "..." ou `"..."' (aspas compostas).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class LogicalLine:
    text: str
    lineno: int                  # linha (1-based) em que o comando começa
    raw: tuple[str, ...] = ()    # linhas físicas originais, para o eco
    kind: str = "cmd"            # "cmd", "comment" ou "delimit"

    @property
    def echo_lines(self) -> tuple[str, ...]:
        return self.raw or (self.text,)


def _is_delimit(cmd: str) -> str | None:
    """Se `cmd` é uma diretiva #delimit, devolve 'cr' ou ';'."""
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
    """Comentários de linha inteira (`*` e `//`) também entram na lista, com
    kind="comment": o Stata os ecoa ao executar um do-file."""
    physical = source.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    out: list[LogicalLine] = []
    mode = "cr"
    buf: list[str] = []
    start_line = 0
    line = 1
    i = 0
    n = len(source)
    block = 0        # profundidade de /* */
    in_str = False   # dentro de "..."
    compound = 0     # profundidade de `" "'

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
                in_str = False  # string não fechada: o comando acusa o erro
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

        # fora de strings e comentários
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
                # linha inteira de comentário: só para o eco
                out.append(LogicalLine("", line, raw_lines(line, line), "comment"))
            i = n if j == -1 else j
            continue
        if c == "\n":
            current = "".join(buf).strip()
            if mode == "cr" or _is_delimit(current) is not None:
                flush()
            elif current.startswith("*"):
                # comentário * em modo ; vai até o próximo ;
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
