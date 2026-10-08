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
        self._run_range(lines, 0, len(lines), echo=echo)

    def _run_range(self, lines: list[LogicalLine], start: int, end: int, *, echo: bool) -> None:
        i = start
        while i < end:
            i = self._run_one(lines, i, end, echo=echo)

    def _run_one(self, lines: list[LogicalLine], i: int, end: int, *, echo: bool) -> int:
        raw = lines[i].text.strip()
        word = _first_word(raw)

        if raw == "}":
            raise StataError(198, "unexpected }")

        if raw.endswith("{"):
            j = self._block_end(lines, i, end)
            if echo:
                self._echo_block(lines, i, j)
            return self._run_block(lines, i, j, end)

        if echo:
            self.s.output.echo_command(raw)

        if raw.startswith("#delimit"):
            # VERIFICAR: mensagem exibida pelo Stata após #delimit
            self.s.output.write(f"delimiter now {raw.split()[-1]}\n", "text")
            return i + 1
        if word == "if":
            return self._run_if_chain(lines, i, end)
        if word == "else":
            raise StataError(198, "else without if")

        self.execute(raw)
        return i + 1

    def _echo_block(self, lines: list[LogicalLine], i: int, j: int) -> None:
        out = self.s.output
        out.echo_command(lines[i].text.strip())
        k = 2
        for ln in lines[i + 1:j + 1]:
            out.write(f"{k:>3}. {ln.text.strip()}\n", "command")
            k += 1
        # else encadeados também são ecoados
        nxt = j + 1
        while nxt < len(lines) and _first_word(lines[nxt].text) == "else":
            if lines[nxt].text.strip().endswith("{"):
                jj = self._block_end(lines, nxt, len(lines))
            else:
                jj = nxt
            for ln in lines[nxt:jj + 1]:
                out.write(f"{k:>3}. {ln.text.strip()}\n", "command")
                k += 1
            nxt = jj + 1

    # ------------------------------------------------------------------
    # Blocos
    # ------------------------------------------------------------------
    @staticmethod
    def _block_end(lines: list[LogicalLine], i: int, end: int) -> int:
        depth = 0
        for k in range(i, end):
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

    def _run_block(self, lines: list[LogicalLine], i: int, j: int, end: int) -> int:
        raw = lines[i].text.strip()
        header = raw[:-1].rstrip()
        word, rest = split_command(header)

        if word == "if":
            return self._run_if_chain(lines, i, end)
        if word == "else":
            raise StataError(198, "else without if")
        if _is_abbrev(word, "forvalues", 4):
            self._forvalues(rest, lines, i + 1, j)
            return j + 1
        if word == "foreach":
            self._foreach(rest, lines, i + 1, j)
            return j + 1
        if word == "while":
            self._while(rest, lines, i + 1, j)
            return j + 1
        prefix = _is_prefix_word(word)
        if prefix is not None or header == "":
            self._prefixed_block(header, lines, i + 1, j)
            return j + 1
        raise StataError(198, "invalid syntax")

    def _prefixed_block(self, header: str, lines: list[LogicalLine], a: int, b: int) -> None:
        """quietly { ... }, capture { ... }, capture noisily { ... }"""
        words = header.replace(":", " ").split()
        self._with_prefixes(words, lambda: self._run_range(lines, a, b, echo=False))

    # -- if / else ------------------------------------------------------------
    def _run_if_chain(self, lines: list[LogicalLine], i: int, end: int) -> int:
        taken = False
        idx = i
        first = True
        while idx < end:
            raw = lines[idx].text.strip()
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

            is_block = raw.endswith("{")
            if is_block:
                j = self._block_end(lines, idx, end)
                body = (idx + 1, j)
                nxt = j + 1
                cond_src = cond_text[:-1].strip() if kind == "if" else ""
                if kind == "else" and cond_text[:-1].strip():
                    raise syntax_error()
                single_cmd = None
            else:
                body = None
                nxt = idx + 1
                cond_src = cond_text

            if not taken:
                if kind == "if":
                    expanded = self.s.expand(cond_src)
                    if is_block:
                        ok = self._truth(expanded)
                        single_cmd = None
                    else:
                        p = Parser(expanded, stop_on_unknown=True)
                        node = p.parse_expr()
                        ok = self._truth_node(node)
                        single_cmd = expanded[p.offset:].strip()
                        if not single_cmd:
                            raise StataError(198, "invalid syntax")
                else:
                    ok = True
                    single_cmd = None if is_block else self.s.expand(cond_text)
                if ok:
                    taken = True
                    if body is not None:
                        self._run_range(lines, body[0], body[1], echo=False)
                    elif single_cmd:
                        self._execute_expanded(single_cmd)
            idx = nxt
            if kind == "else" and not raw.startswith("else if") and _first_word(rest) != "if":
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
            self._run_range(lines, a, b, echo=False)
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
        self._execute_expanded(text)

    def _execute_expanded(self, text: str) -> None:
        from ..commands.registry import lookup

        text = text.strip()
        if not text or text.startswith("*"):
            return
        word, rest = split_command(text)
        if not word:
            raise StataError(198, "invalid syntax")

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

        spec = lookup(word)
        if spec is None:
            raise StataError(199, f"command {word} is unrecognized")
        spec.fn(self.s, rest)

    def _with_prefixes(self, words: list[str], action) -> None:
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
                    self._with_prefixes(inner, action)
                else:
                    action()
                self.s.set_rc(0)
            except StataError as e:
                if noisy:
                    out.error(e.message, e.rc)
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
