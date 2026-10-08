"""Renderização de SMCL (o formato de help e logs) para HTML e para texto.

SMCL mistura texto com diretivas entre chaves. Fora de parágrafos, cada
linha sai como está (tabelas alinhadas por espaço); dentro de {p}...{p_end}
(ou {pstd}, {phang}, {pmore}...), as linhas formam um parágrafo corrido.

Diretivas tratadas: {title:} {cmd:} {cmdab:a:b} {opt} {opth} {it:} {bf:}
{ul:} {hi:} {res:} {txt:} {err:} {com:} {inp:} {help} {helpb} {browse}
{view} {manhelp}/{manlink} (só o texto), {p ...} {pstd} {phang} {pmore}
{p2colset}/{p2col} {synoptset}/{synopt}/{p2coldent} {synoptline}
{synopthdr} {hline} {.-} {c X} {col #} {space #} {tab} {break} {...}
{marker} {dlgtab} {right:} {center:} {bind:} {varlist} {newvar} {ifin}
{weight} {depvar} {indepvars} e {* comentário}. Diretivas desconhecidas
mostram o texto depois dos dois-pontos, se houver.
"""

from __future__ import annotations

import html
import re
import textwrap
from dataclasses import dataclass, field

_CHARS = {"|": "|", "-": "-", "+": "+", "TT": "-", "BT": "-", "LT": "+", "RT": "+",
          "TLC": "+", "TRC": "+", "BLC": "+", "BRC": "+", "-(": "{", ")-": "}", "S|": "$",
          "'g": "`", "a'": "á", "e'": "é", "i'": "í", "o'": "ó", "u'": "ú", "c,": "ç",
          "a~": "ã", "o~": "õ", "a^": "â", "e^": "ê", "o^": "ô", "a`": "à"}


@dataclass
class Span:
    text: str
    style: str = ""             # "", cmd, it, bf, ul, title, link, res, err, hi
    href: str = ""


@dataclass
class Block:
    kind: str                   # line | para | rule | title
    spans: list[Span] = field(default_factory=list)
    indent: int = 0             # parágrafos: recuo da 1ª linha e das demais
    subindent: int = 0


def _find_close(text: str, i: int) -> int:
    """Índice da '}' que fecha a '{' em i (com aninhamento)."""
    depth = 0
    j = i
    while j < len(text):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return j
        j += 1
    return -1


def _split_directive(body: str) -> tuple[str, str, str | None]:
    """'cmd:texto' -> ('cmd', '', 'texto'); 'help regress:texto' -> ('help', 'regress', 'texto')."""
    depth = 0
    quoted = False
    for k, ch in enumerate(body):
        if ch == '"':
            quoted = not quoted
        elif quoted:
            continue
        elif ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch == ":" and depth == 0:
            head, text = body[:k], body[k + 1:]
            break
    else:
        head, text = body, None
    parts = head.split(None, 1)
    name = parts[0] if parts else ""
    args = parts[1] if len(parts) > 1 else ""
    return name, args, text


class _Renderer:
    def __init__(self):
        self.blocks: list[Block] = []
        self.para: Block | None = None
        self.p2col = (4, 20)        # recuo e coluna das tabelas p2col/synopt

    def spans(self, text: str, style: str = "") -> list[Span]:
        out: list[Span] = []
        i = 0
        buf = ""
        while i < len(text):
            ch = text[i]
            if ch == "{":
                j = _find_close(text, i)
                if j == -1:
                    buf += text[i:]
                    break
                if buf:
                    out.append(Span(buf, style))
                    buf = ""
                out.extend(self.directive(text[i + 1:j], style))
                i = j + 1
                continue
            buf += ch
            i += 1
        if buf:
            out.append(Span(buf, style))
        return out

    def directive(self, body: str, style: str) -> list[Span]:
        if body.startswith("*"):
            return []
        name, args, text = _split_directive(body)
        n = name.lower()
        styles = {"cmd": "cmd", "bf": "bf", "it": "it", "ul": "ul", "hi": "hi", "res": "res",
                  "result": "res", "err": "err", "error": "err", "txt": "", "text": "",
                  "com": "cmd", "inp": "cmd", "input": "cmd", "sf": "", "title": "title",
                  "hilite": "hi", "bind": style, "right": style, "center": style, "lalign": style,
                  "ralign": style}
        if n in styles:
            return self.spans(text or args or "", styles[n] or style)
        if n in ("cmdab", "opt", "opth", "opt2") and args and text is not None:
            # {opt s:hort}: os dois-pontos marcam a abreviação, não o texto
            args, text = f"{args}:{text}", None
        if n == "cmdab":
            first, _, rest = (text or args).partition(":")
            return [Span(first, "cmdab"), Span(rest, "cmd")]
        if n in ("opt", "opth", "opt2"):
            t = text if text is not None else args
            m = re.fullmatch(r"([^(\s]+)(\(.*\))?", t.strip())
            if not m:
                return self.spans(t, "cmd")
            nm, arg = m.group(1), m.group(2) or ""
            ab, _, rest = nm.partition(":")
            return ([Span(ab, "cmdab"), Span(rest, "cmd")] if rest else [Span(nm, "cmd")]) + \
                ([Span(arg, "it")] if arg else [])
        if n in ("help", "helpb", "manhelp", "manhelpi", "manlink", "manlinki", "view",
                 "browse", "search", "net", "stata", "dialog"):
            target = args.strip().strip('"')
            label = text if text is not None else target
            if n in ("manhelp", "manhelpi"):
                parts = target.split()
                target = parts[0] if parts else target
                label = text if text is not None else target
            if n in ("manlink", "manlinki"):
                return self.spans(label, "bf")
            href = {"browse": target, "view": f"view:{target}", "search": f"help:{target}",
                    "stata": f"stata:{target}"}.get(n, f"help:{target}")
            return [Span(sp.text, "link", href) for sp in self.spans(label, style)]
        if n in ("varlist", "varname", "newvar", "depvar", "indepvars", "ifin", "weight",
                 "dtype", "exp", "newvarlist", "vars"):
            words = {"ifin": "[if] [in]", "weight": "[weight]", "dtype": "[type]",
                     "exp": "exp"}.get(n, n)
            return [Span(words, "it")]
        if n == "c":
            return [Span(_CHARS.get(args.strip(), args.strip()), style)]
        if n in ("hline", ".-"):
            k = args.strip()
            if k.isdigit():
                return [Span("-" * int(k), style)]
            return [Span("\x00RULE", style)]
        if n in ("space",):
            return [Span(" " * int(args or 1), style)]
        if n in ("col", "column"):
            return [Span(f"\x00COL{int(args or 1)}", style)]
        if n == "tab":
            return [Span("    ", style)]
        if n in ("break", "p_break"):
            return [Span("\x00BR", style)]
        if n == "dup":
            k = int(args or 1)
            return self.spans((text or "") * k, style)
        if n == "char":
            return [Span(chr(int(args)) if args.strip().isdigit() else args, style)]
        if n in ("marker", "viewerjumpto", "vieweralsosee", "viewerdialog", "findalias",
                 "smcl", "reset", "ul", "bind"):
            return []
        if n == "dlgtab":
            return [Span(f"\x00DLGTAB{text or args}", style)]
        return self.spans(text, style) if text is not None else []

    # -- linhas e parágrafos -------------------------------------------------------------------
    _SETTINGS_ONLY = re.compile(r"^(\{(p2colset|synoptset|p2colreset|marker|viewerjumpto|"
                                r"vieweralsosee|viewerdialog|findalias|\*)[^{}]*(\{[^{}]*\})?[^{}]*\})+$")

    def feed(self, source: str) -> list[Block]:
        source = source.replace("\r\n", "\n")
        lines = source.split("\n")
        pending = ""
        for raw in lines:
            if raw.endswith("{...}"):
                body = raw[:-5]
                if self._SETTINGS_ONLY.match(body.strip() or "{*}") and not pending:
                    self.line(body)          # só ajustes: não gera linha nem se junta
                    continue
                pending += body              # {...}: o conteúdo continua na linha seguinte
                continue
            self.line(pending + raw)
            pending = ""
        if pending:
            self.line(pending)
        self.end_para()
        return self.blocks

    def end_para(self) -> None:
        if self.para is not None:
            self.blocks.append(self.para)
            self.para = None

    def line(self, raw: str) -> None:
        s = raw
        m = re.match(r"^\{(p|pstd|phang|phang2|phang3|pmore|pmore2|pin|pin2|p2col|p2coldent|"
                     r"synopt|syntab|p2line|synoptline|synopthdr|p_end|pstd|title|marker|"
                     r"p2colset|synoptset|p2colreset|psee|pmore3)\b([^}]*)\}", s)
        if m:
            tag, args = m.group(1), m.group(2)
            rest = s[m.end():]
            if tag == "p_end":
                self.end_para()
                if rest.strip():
                    self.line(rest)
                return
            if tag in ("p2colset", "synoptset"):
                nums = [int(x) for x in re.findall(r"\d+", args)]
                if tag == "synoptset":
                    self.p2col = (2, nums[0] + 2 if nums else 20)
                elif len(nums) >= 2:
                    self.p2col = (nums[0], nums[1])
                if rest.strip():
                    self.line(rest)
                return
            if tag == "p2colreset":
                self.p2col = (4, 20)
                if rest.strip():
                    self.line(rest)
                return
            if tag in ("synoptline", "p2line"):
                self.end_para()
                self.blocks.append(Block("rule", indent=self.p2col[0]))
                return
            if tag == "synopthdr":
                self.end_para()
                label = args.strip(": ") or "options"
                width = self.p2col[1] - self.p2col[0]
                self.blocks.append(Block("line", [Span(" " * self.p2col[0] + f"{label:<{width}}"
                                                        "Description", "bf")]))
                return
            if tag in ("synopt", "p2col", "p2coldent"):
                # {synopt:{opt x}}descrição -> coluna esquerda e texto
                self.end_para()
                left_body = s[m.start():]
                j = _find_close(left_body, 0)
                inner = left_body[1:j]
                _, _, left = _split_directive(inner)
                desc = left_body[j + 1:]
                lead = 2 if tag == "p2coldent" else self.p2col[0]
                self.para = Block("para", indent=lead, subindent=self.p2col[1])
                self.para.spans = self.spans(left or "") + [Span("\x00COL" + str(self.p2col[1] + 1))]
                if desc.strip():
                    self.para.spans += self.spans(desc)
                return
            if tag == "title":
                self.end_para()
                self.blocks.append(Block("title", self.spans(s)))
                return
            if tag == "marker":
                if rest.strip():
                    self.line(rest)
                return
            self.end_para()
            ind = {"pstd": (4, 4), "phang": (4, 8), "phang2": (8, 12), "phang3": (12, 16),
                   "pmore": (8, 8), "pmore2": (12, 12), "pmore3": (16, 16), "pin": (8, 8),
                   "pin2": (12, 12), "psee": (4, 13)}.get(tag)
            if ind is None:
                nums = [int(x) for x in re.findall(r"\d+", args)]
                ind = (nums[0] if nums else 0, nums[1] if len(nums) > 1 else (nums[0] if nums else 0))
            self.para = Block("para", indent=ind[0], subindent=ind[1])
            if rest.strip():
                self.para.spans += self.spans(rest)
            return
        if self.para is not None:
            if not s.strip():
                self.end_para()
                self.blocks.append(Block("line", []))
                return
            if self.para.spans:
                self.para.spans.append(Span(" "))
            self.para.spans += self.spans(s.strip())
            return
        if s.strip().startswith("{smcl}") and not s.strip()[6:].strip():
            return
        spans = self.spans(s)
        if len(spans) == 1 and spans[0].text == "\x00RULE":
            self.blocks.append(Block("rule"))
            return
        self.blocks.append(Block("line", spans))


def parse(source: str) -> list[Block]:
    return _Renderer().feed(source)


# ---------------------------------------------------------------------------
# saída em texto (console, modo batch)
# ---------------------------------------------------------------------------

def _plain(spans: list[Span], col0: int = 0) -> str:
    out = ""
    for sp in spans:
        t = sp.text
        if t.startswith("\x00COL"):
            target = int(t[4:]) - 1
            line_len = len(out.rsplit("\n", 1)[-1]) + col0
            out += " " * max(1, target - line_len)
            continue
        if t.startswith("\x00DLGTAB"):
            out += f"+--- {t[7:]} " + "-" * max(4, 60 - len(t))
            continue
        if t == "\x00RULE":
            out += "-" * 79
            continue
        if t == "\x00BR":
            out += "\n"
            continue
        out += t
    return out


def render_text(source: str, width: int = 79) -> str:
    lines: list[str] = []
    for b in parse(source):
        if b.kind == "rule":
            lines.append(" " * b.indent + "-" * (width - b.indent))
        elif b.kind == "para":
            pieces = _plain(b.spans, b.indent).split("\n")
            for k, piece in enumerate(pieces):
                lead = b.indent if k == 0 else b.subindent
                lines.append(textwrap.fill(" " * lead + piece.strip(), width=width,
                                           subsequent_indent=" " * b.subindent,
                                           break_long_words=False, break_on_hyphens=False))
        else:
            lines.append(_plain(b.spans).rstrip())
    return "\n".join(lines).rstrip() + "\n"


# ---------------------------------------------------------------------------
# saída em HTML (Viewer)
# ---------------------------------------------------------------------------

def _html_spans(spans: list[Span], colors: dict[str, str]) -> str:
    out = []
    col = 0
    for sp in spans:
        t = sp.text
        if t.startswith("\x00COL"):
            target = int(t[4:]) - 1
            pad = max(1, target - col)
            out.append("&nbsp;" * pad)
            col += pad
            continue
        if t.startswith("\x00DLGTAB"):
            label = html.escape(t[7:])
            out.append(f'<span style="color:{colors["dim"]}">+--- </span><b>{label}</b>'
                       f'<span style="color:{colors["dim"]}"> {"-" * 50}</span>')
            continue
        if t == "\x00RULE":
            out.append(f'<span style="color:{colors["dim"]}">{"-" * 79}</span>')
            continue
        if t == "\x00BR":
            out.append("<br>")
            col = 0
            continue
        col += len(t)
        e = html.escape(t).replace("  ", "&nbsp; ")
        if sp.style == "link":
            out.append(f'<a href="{html.escape(sp.href)}" style="color:{colors["link"]}">{e}</a>')
        elif sp.style in ("cmd", "bf", "title"):
            out.append(f"<b>{e}</b>")
        elif sp.style == "cmdab":
            out.append(f"<b><u>{e}</u></b>")
        elif sp.style == "it":
            out.append(f"<i>{e}</i>")
        elif sp.style == "ul":
            out.append(f"<u>{e}</u>")
        elif sp.style == "res":
            out.append(f'<span style="color:{colors["res"]}">{e}</span>')
        elif sp.style == "err":
            out.append(f'<span style="color:{colors["err"]}">{e}</span>')
        elif sp.style == "hi":
            out.append(f'<b style="color:{colors["res"]}">{e}</b>')
        else:
            out.append(e)
    return "".join(out)


LIGHT = {"text": "#000000", "res": "#000000", "err": "#cc0000", "link": "#0a5fd1", "dim": "#8e8e93",
         "bg": "#ffffff"}
DARK = {"text": "#d4d4d6", "res": "#ffffff", "err": "#ff6b68", "link": "#6cb6ff", "dim": "#8e8e93",
        "bg": "#1e1f22"}


def render_html(source: str, *, dark: bool = False, font: str = "Menlo", size: int = 12) -> str:
    colors = DARK if dark else LIGHT
    parts = [f'<html><body style="background:{colors["bg"]}; color:{colors["text"]}; '
             f'font-family:\'{font}\', monospace; font-size:{size}pt;">']
    for b in parse(source):
        if b.kind == "rule":
            parts.append(f'<div style="white-space:pre; color:{colors["dim"]}">'
                         f'{"&nbsp;" * b.indent}{"-" * (79 - b.indent)}</div>')
        elif b.kind == "para":
            # o QTextBrowser não entende a unidade ch: recuo em pixels
            cw = size * 0.8
            parts.append(f'<div style="margin-left:{b.subindent * cw:.0f}px; '
                         f'text-indent:{(b.indent - b.subindent) * cw:.0f}px; '
                         f'white-space:normal; margin-top:2px; margin-bottom:2px;">'
                         f"{_html_spans(b.spans, colors)}</div>")
        elif b.kind == "title":
            parts.append(f'<div style="white-space:pre; font-weight:bold;">{_html_spans(b.spans, colors)}</div>')
        else:
            body = _html_spans(b.spans, colors) or "&nbsp;"
            parts.append(f'<div style="white-space:pre;">{body}</div>')
    parts.append("</body></html>")
    return "\n".join(parts)
