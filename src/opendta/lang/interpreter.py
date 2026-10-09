"""Execução de linhas lógicas: blocos, laços, condicionais e prefixos.

Blocos seguem as regras do manual [P]: a chave de abertura `{` é o último
caractere da linha do comando; a chave de fechamento `}` fica sozinha na
sua linha. `else` vem na linha seguinte ao `}` do `if`.

Ao executar um do-file com eco, um bloco de nível superior é ecoado por
inteiro antes da execução, com as linhas internas numeradas ("  2. ..."),
como o Stata faz.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

from ..core.errors import (BreakLoop, ContinueLoop, StataError, syntax_error)
from .expr import Parser, evaluate
from .lexer import LogicalLine
from .words import num_str, parse_numlist, split_words

if TYPE_CHECKING:
    from ..session import Session

_CMD_WORD = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)")


def split_command(text: str) -> tuple[str, str]:
    """'display 2+2' -> ('display', '2+2'); 'di"x"' -> ('di', '"x"')."""
    m = _CMD_WORD.match(text)
    if not m:
        return "", text.strip()
    return m.group(1), text[m.end():].strip()


def _first_word(text: str) -> str:
    return split_command(text)[0]


def _is_abbrev(word: str, full: str, minimum: int) -> bool:
    return len(word) >= minimum and full.startswith(word)


def _is_prefix_word(word: str) -> str | None:
    for full, minimum in (("quietly", 3), ("noisily", 1), ("capture", 3)):
        if _is_abbrev(word, full, minimum):
            return full
    return None


class Interpreter:
    def __init__(self, session: "Session"):
        self.s = session

    # ------------------------------------------------------------------
    # Execução de uma sequência de linhas
    # ------------------------------------------------------------------
    def run_lines(self, lines: list[LogicalLine], *, echo: bool = False) -> None:
        self._run_range(lines, 0, len(lines), echo="top" if echo else None)

    def _run_range(self, lines: list[LogicalLine], start: int, end: int, *, echo: str | None) -> None:
        i = start
        while i < end:
            i = self._run_one(lines, i, end, echo=echo)

    # Modos de eco (observados nos logs do Stata):
    #   "top"    comando de nível superior do do-file: '. linha', saída, linha em branco
    #   "inline" linha dentro de if/else ou de bloco quietly/capture: '. linha', sem
    #            linha em branco depois
    #   None     sem eco (corpo de laços, que são ecoados inteiros antes de rodar)
    def _echo(self, ln: LogicalLine, echo: str | None) -> None:
        if echo:
            self.s.output.echo_command(ln.echo_lines)

    def _after(self, echo: str | None) -> None:
        if echo == "top":
            self.s.output.end_command()

    def _run_one(self, lines: list[LogicalLine], i: int, end: int, *, echo: str | None) -> int:
        ln = lines[i]
        out = self.s.output
        if ln.kind == "comment":
            self._echo(ln, echo)
            return i + 1
        if ln.kind == "delimit":
            self._echo(ln, echo)
            out.write(f"delimiter now {ln.text.split()[-1]}\n", "text")
            return i + 1

        raw = ln.text.strip()
        word = _first_word(raw)
        if raw == "}":
            raise StataError(198, "unexpected }")
        if len(word) >= 2 and "program".startswith(word) and self._is_definition(raw):
            return self._define_program(lines, i, end, echo=echo)
        if raw.endswith("{"):
            j = self._block_end(lines, i, end)
            return self._run_block(lines, i, j, end, echo=echo)
        if word == "input":
            return self._run_input(lines, i, end, echo=echo)
        if word == "if":
            return self._run_if_chain(lines, i, end, echo=echo)
        if word == "else":
            raise StataError(198, "else without if")

        self._echo(ln, echo)
        self.execute(raw)
        self._after(echo)
        return i + 1

    # -- programas --------------------------------------------------------------
    @staticmethod
    def _is_definition(raw: str) -> bool:
        rest = split_command(raw)[1].strip()
        first = rest.split(",")[0].split()
        if not first:
            return False
        return first[0] not in ("drop", "dir", "list", "li", "l")

    def _define_program(self, lines: list[LogicalLine], i: int, end: int, *, echo: str | None) -> int:
        """program [define] nome ... end: guarda o corpo sem executar."""
        from .programs import make_program, parse_definition

        header = self.s.expand(lines[i].text.strip())
        name, options = parse_definition(split_command(header)[1])
        j = i + 1
        while j < end and not (lines[j].kind == "cmd" and lines[j].text.strip() == "end"):
            j += 1
        if j >= end:
            raise StataError(198, "program define: end not found")   # VERIFICAR
        if echo:
            # VERIFICAR: numeração do corpo no eco de program define
            out = self.s.output
            out.echo_command(lines[i].echo_lines)
            for k, ln in enumerate(lines[i + 1:j + 1], start=1):
                first, *rest = ln.echo_lines
                out.write(f"{k:>3}. {first}\n", "command")
                for cont in rest:
                    out.write(f"> {cont}\n", "command")
        if name in self.s.programs:
            raise StataError(110, f"{name} already defined")
        self.s.programs[name] = make_program(name, options, list(lines[i + 1:j]),
                                             self.s.current_dofile)
        self._after(echo)
        return j + 1

    def _autoload(self, word: str):
        """Comando desconhecido: procura word.ado no adopath e roda o arquivo
        (em silêncio, como o Stata), que deve definir o programa word."""
        from .adopath import find_ado
        path = find_ado(self.s, word)
        if path is None:
            return None
        from ..commands.program import _run_file
        _run_file(self.s, f'"{path}"', echo=False, new_scope=True)
        prog = self.s.programs.get(word)
        if prog is None:
            # VERIFICAR: mensagem do Stata quando o .ado não define o programa
            raise StataError(199, f"{path.name} found but program {word} not defined")
        prog.source = str(path)
        return prog

    def run_program(self, prog) -> None:
        """Corpo de um programa: sem eco (set trace mostra as linhas)."""
        tracing = self._tracing()
        # linhas begin/end ocupam a largura da tela (set linesize; observado
        # no Stata 14 com linesize 255)
        width = int(float(self.s.settings.get("linesize", 80)))
        pad = "  " * self._trace_depth()
        if tracing:
            label = f" begin {prog.name} ---"
            self.s.output.write(pad + "-" * max(4, width - len(pad) - len(label)) + label + "\n",
                                "text", force=True)
        try:
            self._run_range(prog.lines, 0, len(prog.lines), echo=None)
        finally:
            if tracing:
                label = f" end {prog.name} ---"
                self.s.output.write(pad + "-" * max(4, width - len(pad) - len(label)) + label + "\n",
                                    "text", force=True)

    def _run_input(self, lines: list[LogicalLine], i: int, end: int, *, echo: str | None) -> int:
        """input var1 var2 ... seguido de linhas de dados até `end`."""
        from ..commands.data import run_input

        raw = self.s.expand(lines[i].text.strip())
        spec = split_command(raw)[1]
        j = i + 1
        rows: list[str] = []
        while j < end and lines[j].text.strip() != "end":
            rows.append(self.s.expand(lines[j].text.strip()))
            j += 1
        if j >= end:
            raise StataError(198, "input: end not found")
        out = self.s.output
        if echo:
            self._echo(lines[i], echo)
            # cabeçalho: 3 espaços e cada nome alinhado à direita numa coluna
            # com a largura do formato do tipo + 2 (float %9.0g → 11,
            # double %10.0g → 12, str20 %20s → 22; observado no Stata 14)
            from ..core.dataset import default_format
            cells = []
            pending_type = None
            for n in spec.split():
                if re.match(r"^(byte|int|long|float|double|str\d*|strL)$", n):
                    pending_type = n
                    continue
                if self.s.data.has(n):
                    fmt = self.s.data.get(n).fmt
                else:
                    fmt = default_format(pending_type or self.s.settings.get("type", "float"))
                pending_type = None
                m = re.match(r"^%-?(\d+)", fmt)
                width = (int(m.group(1)) if m else 9) + 2
                cells.append(f"{n:>{width}}")
            out.write("\n   " + "".join(cells) + "\n", "text")
            for k, ln in enumerate(lines[i + 1:j + 1], start=1):
                out.write(f"{k:>3}. {ln.echo_lines[0].strip()}\n", "command")
        run_input(self.s, spec, rows)
        self._after(echo)
        return j + 1

    def _echo_numbered(self, lines: list[LogicalLine], i: int, j: int) -> None:
        """Eco de laço: cabeçalho e corpo numerado, antes da execução."""
        out = self.s.output
        out.echo_command(lines[i].echo_lines)
        k = 2
        for ln in lines[i + 1:j + 1]:
            first, *rest = ln.echo_lines
            out.write(f"{k:>3}. {first}\n", "command")
            for cont in rest:
                out.write(f"> {cont}\n", "command")
            k += 1

    # ------------------------------------------------------------------
    # Blocos
    # ------------------------------------------------------------------
    @staticmethod
    def _block_end(lines: list[LogicalLine], i: int, end: int) -> int:
        depth = 0
        for k in range(i, end):
            if lines[k].kind != "cmd":
                continue
            t = lines[k].text.strip()
            if t == "}":
                depth -= 1
                if depth == 0:
                    return k
            elif t.startswith("}"):
                raise StataError(198, "} must appear on a line by itself")
            if t.endswith("{"):
                depth += 1
        raise StataError(198, "unexpected end of file")

    def _run_block(self, lines: list[LogicalLine], i: int, j: int, end: int, *, echo: str | None) -> int:
        raw = lines[i].text.strip()
        header = raw[:-1].rstrip()
        word, rest = split_command(header)

        if word == "if":
            return self._run_if_chain(lines, i, end, echo=echo)
        if word == "else":
            raise StataError(198, "else without if")
        loop = None
        if _is_abbrev(word, "forvalues", 4):
            loop = self._forvalues
        elif word == "foreach":
            loop = self._foreach
        elif word == "while":
            loop = self._while
        if loop is not None:
            if echo:
                self._echo_numbered(lines, i, j)
            loop(rest, lines, i + 1, j)
            self._after(echo)
            return j + 1
        prefix = _is_prefix_word(word)
        if prefix is not None or header == "":
            self._prefixed_block(header, lines, i, j, echo=echo)
            return j + 1
        raise StataError(198, "invalid syntax")

    def _prefixed_block(self, header: str, lines: list[LogicalLine], i: int, j: int,
                        *, echo: str | None) -> None:
        """quietly { ... }, capture { ... }, capture noisily { ... }

        O cabeçalho é ecoado; as linhas internas são ecoadas uma a uma ao
        rodar (e somem sob quietly/capture, como no Stata)."""
        words = [_is_prefix_word(w) or w for w in header.replace(":", " ").split()]
        self._echo(lines[i], echo)
        inner = "inline" if echo else None
        self._with_prefixes(words, lambda: self._run_range(lines, i + 1, j, echo=inner), block=True)
        visible = "noisily" in words or not ({"quietly", "capture"} & set(words))
        if echo and visible:
            self._echo(lines[j], inner)
        self._after(echo)

    # -- if / else ------------------------------------------------------------
    def _run_if_chain(self, lines: list[LogicalLine], i: int, end: int, *, echo: str | None) -> int:
        """if / else if / else. Com eco, cada linha aparece ao ser alcançada,
        inclusive as do ramo não executado (que não produzem saída)."""
        taken = False
        idx = i
        first = True
        inner = "inline" if echo else None
        while idx < end:
            ln = lines[idx]
            if ln.kind != "cmd":
                break
            raw = ln.text.strip()
            word, rest = split_command(raw)
            if first:
                if word != "if":
                    raise syntax_error()
                kind, cond_text = "if", rest
            else:
                if word != "else":
                    break
                w2, rest2 = split_command(rest)
                if w2 == "if":
                    kind, cond_text = "if", rest2
                else:
                    kind, cond_text = "else", rest
            first = False
            self._echo(ln, echo)

            is_block = raw.endswith("{")
            if is_block:
                j = self._block_end(lines, idx, end)
                if kind == "else" and cond_text[:-1].strip():
                    raise syntax_error()
                run_it = False
                if not taken:
                    run_it = True if kind == "else" else self._truth(self.s.expand(cond_text[:-1].strip()))
                if run_it:
                    taken = True
                    self._run_range(lines, idx + 1, j, echo=inner)
                elif echo:
                    for ln2 in lines[idx + 1:j]:
                        self._echo(ln2, inner)
                self._echo(lines[j], inner)
                nxt = j + 1
            else:
                if not taken:
                    if kind == "if":
                        expanded = self.s.expand(cond_text)
                        p = Parser(expanded, stop_on_unknown=True)
                        node = p.parse_expr()
                        ok = self._truth_node(node)
                        cmd = expanded[p.offset:].strip()
                        if not cmd:
                            raise StataError(198, "invalid syntax")
                    else:
                        ok, cmd = True, self.s.expand(cond_text)
                    if ok:
                        taken = True
                        self._execute_expanded(cmd)
                nxt = idx + 1
            self._after(echo)
            idx = nxt
            if kind == "else":
                break
        return idx

    def _truth(self, expr_text: str) -> bool:
        p = Parser(expr_text)
        return self._truth_node(p.parse_full())

    def _truth_node(self, node) -> bool:
        v = evaluate(node, self.s.context)
        if isinstance(v, str):
            raise StataError(109, "type mismatch")
        return v != 0

    # -- laços ----------------------------------------------------------------
    def _loop_body(self, lines: list[LogicalLine], a: int, b: int) -> bool:
        """Executa o corpo; devolve False se houve `continue, break`."""
        try:
            self._run_range(lines, a, b, echo=None)
        except ContinueLoop:
            return True
        except BreakLoop:
            return False
        return True

    def _forvalues(self, spec: str, lines: list[LogicalLine], a: int, b: int) -> None:
        text = self.s.expand(spec)
        m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+?)\s*$", text)
        if not m:
            raise StataError(198, "invalid syntax")
        name, rng = m.group(1), m.group(2)
        values = _forvalues_range(rng)
        for v in values:
            self.s.macros.set_local(name, num_str(v))
            if not self._loop_body(lines, a, b):
                break

    def _foreach(self, spec: str, lines: list[LogicalLine], a: int, b: int) -> None:
        text = self.s.expand(spec)
        m = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s+(in|of)\s*(.*)$", text, re.S)
        if not m:
            raise StataError(198, "invalid syntax")
        name, kw, rest = m.group(1), m.group(2), m.group(3).strip()
        if kw == "in":
            items = split_words(rest)
        else:
            ltype, _, arg = rest.partition(" ")
            arg = arg.strip()
            if ltype == "local":
                items = split_words(self.s.macros.get_local(arg))
            elif ltype == "global":
                items = split_words(self.s.macros.get_global(arg))
            elif ltype == "numlist":
                items = [num_str(x) for x in parse_numlist(arg)]
            elif ltype == "newlist":
                items = arg.split()
            elif ltype == "varlist":
                items = self.s.expand_varlist(arg)
            else:
                raise StataError(198, f"invalid syntax")
        for it in items:
            self.s.macros.set_local(name, it)
            if not self._loop_body(lines, a, b):
                break

    def _while(self, cond: str, lines: list[LogicalLine], a: int, b: int) -> None:
        while self._truth(self.s.expand(cond)):
            if not self._loop_body(lines, a, b):
                break

    # ------------------------------------------------------------------
    # Comando individual
    # ------------------------------------------------------------------
    def execute(self, raw: str) -> None:
        """Expande macros e executa um comando de uma linha."""
        text = self.s.expand(raw)
        if self._tracing():
            # set trace on: linha original (- ) e expandida (= ) dentro de programas
            out = self.s.output
            pad = "  " * self._trace_depth()
            out.write(f"{pad}- {raw.strip()}\n", "text", force=True)
            if text.strip() != raw.strip():
                out.write(f"{pad}= {text.strip()}\n", "text", force=True)
        self._execute_expanded(text)

    def _tracing(self) -> bool:
        return (self.s.settings.get("trace", "off") == "on"
                and any(sc.kind == "program" for sc in self.s.scopes))

    def _trace_depth(self) -> int:
        # um programa chamado de um do-file aparece com 4 espaços (Stata 14);
        # VERIFICAR o recuo de programas aninhados
        return sum(1 for sc in self.s.scopes if sc.kind == "program") + 1

    def _execute_expanded(self, text: str) -> None:
        from ..commands.registry import lookup

        text = text.strip()
        if not text or text.startswith("*"):
            return
        word, rest = split_command(text)
        if not word:
            raise StataError(198, "invalid syntax")

        if word == "by" or _is_abbrev(word, "bysort", 3):
            from .syntax import find_top, parse_by
            colon = find_top(rest, ":")
            if colon == -1:
                raise StataError(198, "invalid syntax")
            bp = parse_by(rest[:colon])
            if word != "by":
                bp.sort = True
            self.s.run_by(bp, rest[colon + 1:].strip())
            return

        prefix = _is_prefix_word(word)
        if prefix is not None:
            words = [prefix]
            body = rest
            while True:
                body = body.lstrip()
                if body.startswith(":"):
                    body = body[1:].lstrip()
                w, r = split_command(body)
                p2 = _is_prefix_word(w) if w else None
                if p2 is None:
                    break
                words.append(p2)
                body = r
            self._with_prefixes(words, lambda: self._execute_expanded(body))
            return

        # comandos internos (inclusive abreviados) vêm antes dos programas:
        # com `program define pr`, digitar `pr` ainda chama `program`
        # (observado no Stata 14, compat/expected/0201_programas.log)
        prog = None if lookup(word) is not None else self.s.programs.get(word)
        if prog is None and lookup(word) is None and _plain_name(word):
            prog = self._autoload(word)
        if prog is not None:
            from .programs import call_program
            if self.s.by_groups is not None and not prog.byable:
                raise StataError(190, f"{word} may not be combined with by")
            try:
                call_program(self.s, prog, rest)
            except StataError as err:
                if err.context is None:
                    err.context = (word, text)
                raise
            return
        spec = lookup(word)
        if spec is None:
            err = StataError(199, f"command {word} is unrecognized")
            err.context = (word, text)
            raise err
        if self.s.by_groups is not None and not spec.byable and not spec.prefix:
            raise StataError(190, f"{spec.name} may not be combined with by")
        try:
            spec.fn(self.s, rest)
        except StataError as err:
            if err.context is None:
                err.context = (spec.name, text)   # para a explicação do erro
            raise

    def _with_prefixes(self, words: list[str], action, block: bool = False) -> None:
        out = self.s.output
        full = [_is_prefix_word(w) or w for w in words]
        for w in full:
            if w not in ("quietly", "noisily", "capture"):
                raise StataError(198, "invalid syntax")

        if full and full[0] == "capture":
            inner = full[1:]
            noisy = "noisily" in inner
            if not noisy:
                out.quiet_depth += 1
                out.capture_depth += 1
            try:
                if inner:
                    self._with_prefixes(inner, action, block)
                else:
                    action()
                self.s.set_rc(0)
            except StataError as e:
                if noisy:
                    # capture noisily cmd: só a mensagem; em bloco, também r(#)
                    self.s.report_error(e, show_rc=block)
                self.s.set_rc(e.rc)
            finally:
                if not noisy:
                    out.quiet_depth -= 1
                    out.capture_depth -= 1
            return

        head, tail = full[0], full[1:]
        if head == "quietly":
            out.quiet_depth += 1
            saved = out.noisy_depth
            out.noisy_depth = 0
            try:
                (lambda: self._with_prefixes(tail, action) if tail else action())()
            finally:
                out.quiet_depth -= 1
                out.noisy_depth = saved
        else:  # noisily
            out.noisy_depth += 1
            try:
                (lambda: self._with_prefixes(tail, action) if tail else action())()
            finally:
                out.noisy_depth -= 1


def _forvalues_range(rng: str) -> list[float]:
    r = rng.strip()
    num = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
    m = re.match(rf"^({num})\s*/\s*({num})$", r)
    if m:
        a, b = float(m.group(1)), float(m.group(2))
        return _seq(a, 1.0, b)
    m = re.match(rf"^({num})\s*[\(\[]\s*({num})\s*[\)\]]\s*({num})$", r)
    if m:
        return _seq(float(m.group(1)), float(m.group(2)), float(m.group(3)))
    m = re.match(rf"^({num})\s+({num})\s*(?:to|:)\s*({num})$", r)
    if m:
        a, t, b = (float(g) for g in m.groups())
        return _seq(a, t - a, b)
    raise StataError(198, "invalid syntax")


def _seq(a: float, step: float, b: float) -> list[float]:
    if step == 0:
        raise StataError(198, "invalid syntax")
    out: list[float] = []
    k = 0
    eps = abs(step) * 1e-9
    while True:
        v = a + k * step
        if (step > 0 and v > b + eps) or (step < 0 and v < b - eps):
            return out
        out.append(round(v, 12))
        k += 1


def _plain_name(word: str) -> bool:
    import re as _re
    return _re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]{0,31}", word) is not None
