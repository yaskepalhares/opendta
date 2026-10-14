"""syntax, gettoken, marksample e markout ([P] syntax, [P] gettoken, [P] mark).

syntax lê `0' (o que veio depois do nome do programa) e o confere contra a
descrição dada, criando as locais varlist, if, in, using, exp, weight e uma
local por opção. Letras maiúsculas no nome de uma opção marcam a abreviação
mínima (Detail aceita d, de, ..., detail).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import NAME_RE, Variable
from ..core.errors import StataError
from ..core.formats import number_to_macro
from ..lang.syntax import Parsed, parse_options, parse_standard
from ..lang.words import parse_numlist, strip_outer_quotes
from ._util import touse
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


# ---------------------------------------------------------------------------
# descrição
# ---------------------------------------------------------------------------

@dataclass
class OptSpec:
    name: str                # nome completo, minúsculo
    minlen: int              # abreviação mínima
    kind: str                # flag | string | integer | real | numlist | varlist | varname
    #                          | name | namelist | passthru | asis | star
    default: str = ""
    required: bool = False
    no: bool = False         # noOPTION: local sem o "no"
    asis: bool = False
    extra: str = ""          # especificadores (min=, max=, ascending, integer...)


@dataclass
class Spec:
    main: str = ""           # varlist | newvarlist | varname | newvarname | namelist | name | anything
    main_optional: bool = True
    main_args: str = ""
    if_: str = ""            # "" (não aceito) | "opt" | "req"
    if_slash: bool = False
    in_: str = ""
    in_slash: bool = False
    using: str = ""
    using_slash: bool = False
    exp: str = ""
    exp_slash: bool = False
    weights: list[str] = field(default_factory=list)
    weight_slash: bool = False
    options: list[OptSpec] = field(default_factory=list)


_MAIN = {"varlist", "newvarlist", "varname", "newvarname", "namelist", "name", "anything"}
_WEIGHTS = {"fweight", "aweight", "pweight", "iweight"}


def _tokens(text: str) -> list[str]:
    """Palavras, grupos [ ... ] e nome(...) juntos."""
    out: list[str] = []
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch.isspace():
            i += 1
            continue
        if ch == ",":
            out.append(",")
            i += 1
            continue
        if ch == "[":
            depth, j = 0, i
            while j < n:
                if text[j] == "[":
                    depth += 1
                elif text[j] == "]":
                    depth -= 1
                    if depth == 0:
                        break
                j += 1
            out.append(text[i:j + 1])
            i = j + 1
            continue
        j = i
        depth = 0
        while j < n:
            c = text[j]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
            elif c == '"':
                k = text.find('"', j + 1)
                j = n - 1 if k == -1 else k
            elif depth == 0 and (c.isspace() or c in ",["):
                break
            j += 1
        out.append(text[i:j])
        i = j
    return out


def _option_spec(tok: str, required: bool) -> OptSpec:
    if tok == "*":
        return OptSpec("*", 1, "star")
    m = re.fullmatch(r"([A-Za-z_][A-Za-z0-9_]*)(?:\((.*)\))?", tok, re.S)
    if not m:
        raise StataError(197, f"invalid syntax: {tok}")   # VERIFICAR
    raw, arg = m.group(1), m.group(2)
    no = raw.startswith("no") and len(raw) > 2 and raw[2].isupper()
    caps = len(raw) - len(raw.lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
    if no:
        caps = 2 + len(raw[2:]) - len(raw[2:].lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))
    name = raw.lower()
    minlen = caps if caps else len(name)
    if arg is None:
        return OptSpec(name, minlen, "flag", required=required, no=no)
    words = arg.split()
    kind = words[0].lower() if words else "string"
    rest = " ".join(words[1:])
    if kind not in ("string", "integer", "real", "numlist", "varlist", "varname", "name",
                    "namelist", "passthru", "asis"):
        raise StataError(197, f"invalid syntax: {tok}")   # VERIFICAR
    spec = OptSpec(name, minlen, kind, required=required, extra=rest)
    if kind in ("integer", "real") and words[1:]:
        spec.default = words[1]
    if "asis" in words[1:] or kind == "asis":
        spec.asis = True
    return spec


def parse_spec(text: str) -> Spec:
    spec = Spec()
    in_options = False
    for tok in _tokens(text):
        if tok == ",":
            in_options = True
            continue
        optional = tok.startswith("[")
        body = tok[1:-1].strip() if optional else tok
        if in_options:
            for t in _tokens(body):
                if t == ",":
                    continue
                spec.options.append(_option_spec(t, required=not optional))
            continue
        if optional and body.startswith(","):
            in_options_inner = _tokens(body[1:])
            for t in in_options_inner:
                if t != ",":
                    spec.options.append(_option_spec(t, required=False))
            in_options = True
            continue
        parts = _tokens(body)
        for t in parts:
            word = re.match(r"[A-Za-z=/]+", t)
            w = word.group(0) if word else t
            args = t[len(w):].strip()
            if args.startswith("(") and args.endswith(")"):
                args = args[1:-1]
            wl = w.rstrip("/")
            slash = w.endswith("/")
            if wl in _MAIN:
                spec.main, spec.main_optional, spec.main_args = wl, optional, args
            elif wl == "if":
                spec.if_, spec.if_slash = ("opt" if optional else "req"), slash
            elif wl == "in":
                spec.in_, spec.in_slash = ("opt" if optional else "req"), slash
            elif wl == "using":
                spec.using, spec.using_slash = ("opt" if optional else "req"), slash
            elif t.startswith("="):
                spec.exp = "opt" if optional else "req"
                spec.exp_slash = t.startswith("=/")
            elif wl in _WEIGHTS or t == "/":
                if t == "/":
                    spec.weight_slash = True
                else:
                    spec.weights.append(wl)
            else:
                raise StataError(197, f"invalid syntax: {t}")   # VERIFICAR
    return spec


# ---------------------------------------------------------------------------
# conferência
# ---------------------------------------------------------------------------

def _specifiers(args: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for m in re.finditer(r"(\w+)(?:=(\"[^\"]*\"|\S+))?", args):
        out[m.group(1).lower()] = (m.group(2) or "").strip('"')
    return out


def _main(s: "Session", spec: Spec, text: str, locs: dict[str, str]) -> None:
    kind = spec.main
    sp = _specifiers(spec.main_args)
    lname = sp.get("name", "varlist" if "var" in kind else
                   ("namelist" if kind in ("namelist", "name") else "anything"))
    text = text.strip()
    if not text:
        if not spec.main_optional:
            raise StataError(100, "varlist required" if "var" in kind else "invalid syntax")
        if kind == "varlist" and sp.get("default", "all") != "none":
            locs[lname] = " ".join(s.data.names)
        return
    if kind == "anything":
        locs[lname] = text
        return
    if kind in ("namelist", "name"):
        names = text.split()
        for n in names:
            if not NAME_RE.match(n):
                raise StataError(198, f"{n} invalid name")
        if kind == "name" and len(names) > 1:
            raise StataError(103, "too many names specified")   # VERIFICAR
        _count(names, sp, kind)
        locs[lname] = " ".join(names)
        return
    if kind in ("newvarlist", "newvarname"):
        types, names = [], []
        vtype = s.default_type()
        pending = None
        for w in text.split():
            if w in ("byte", "int", "long", "float", "double") or re.fullmatch(r"str\d+|strL", w):
                pending = w
                continue
            if not NAME_RE.match(w):
                raise StataError(198, f"{w} invalid name")
            if s.data.has(w):
                raise StataError(110, f"{w} already defined")
            names.append(w)
            types.append(pending or vtype)
            pending = None
        if kind == "newvarname" and len(names) > 1:
            raise StataError(103, "too many variables specified")
        _count(names, sp, kind)
        locs[lname] = " ".join(names)
        locs["typlist"] = " ".join(types)
        return
    names = s.expand_varlist(text)
    if kind == "varname" and len(names) > 1:
        raise StataError(103, "too many variables specified")
    _count(names, sp, kind)
    if "numeric" in sp:
        for n in names:
            if s.data.get(n).is_string:
                raise StataError(109, f"string variables not allowed in varlist;\n{n} is a string variable")
    if "string" in sp:
        for n in names:
            if not s.data.get(n).is_string:
                raise StataError(109, f"numeric variables not allowed in varlist;\n{n} is not a string variable")
    locs[lname] = " ".join(names)


def _count(items: list[str], sp: dict[str, str], kind: str) -> None:
    lo = int(sp.get("min", "1") or 1)
    hi = sp.get("max")
    what = "variables" if "var" in kind else "names"
    if len(items) < lo:
        raise StataError(102, f"too few {what} specified")
    if hi and len(items) > int(hi):
        raise StataError(103, f"too many {what} specified")


def _option_value(s: "Session", o: OptSpec, arg: str | None, typed: str) -> str:
    if o.kind == "flag":
        if arg is not None:
            raise StataError(198, f"option {typed} incorrectly specified")   # VERIFICAR
        return o.name
    if arg is None:
        raise StataError(198, f"option {o.name}() incorrectly specified")    # VERIFICAR
    if o.kind == "passthru":
        return f"{o.name}({arg})"
    if o.kind in ("string", "asis"):
        return arg if o.asis else strip_outer_quotes(arg) if arg.strip().startswith('"') else arg
    if o.kind in ("integer", "real"):
        try:
            x = float(s.eval(arg))
        except StataError:
            raise StataError(198, f"option {o.name}() incorrectly specified")
        if o.kind == "integer" and x != int(x):
            raise StataError(198, f"option {o.name}() incorrectly specified")
        return number_to_macro(x)
    if o.kind == "numlist":
        nums = parse_numlist(arg)
        if "integer" in o.extra and any(x != int(x) for x in nums):
            raise StataError(126, "noninteger value not allowed")   # VERIFICAR
        return " ".join(number_to_macro(x) for x in nums)
    if o.kind in ("varlist", "varname"):
        names = s.expand_varlist(arg)
        if o.kind == "varname" and len(names) != 1:
            raise StataError(103 if len(names) > 1 else 102,
                             "too many variables specified" if len(names) > 1 else "too few variables specified")
        return " ".join(names)
    if o.kind in ("name", "namelist"):
        names = arg.split()
        for n in names:
            if not NAME_RE.match(n):
                raise StataError(198, f"{n} invalid name")
        if o.kind == "name" and len(names) != 1:
            raise StataError(198, f"option {o.name}() incorrectly specified")
        return " ".join(names)
    return arg


def run_syntax(s: "Session", description: str, text: str) -> dict[str, str]:
    spec = parse_spec(description)
    p: Parsed = parse_standard(text)
    locs: dict[str, str] = {}

    # parse_standard põe "= exp" em p.exp e a lista antes dele em p.varlist
    if p.exp is not None:
        if not spec.exp:
            raise StataError(101, "=exp not allowed")   # VERIFICAR
        locs["exp"] = p.exp if spec.exp_slash else f"= {p.exp}"
    elif spec.exp == "req":
        raise StataError(100, "=exp required")

    if spec.main:
        _main(s, spec, p.varlist, locs)
    elif p.varlist.strip():
        raise StataError(101, "varlist not allowed")

    if p.if_:
        if not spec.if_:
            raise StataError(101, "if not allowed")
        locs["if"] = p.if_ if spec.if_slash else f"if {p.if_}"
    elif spec.if_ == "req":
        raise StataError(100, "if required")
    if p.in_:
        if not spec.in_:
            raise StataError(101, "in range not allowed")
        locs["in"] = p.in_ if spec.in_slash else f"in {p.in_}"
    elif spec.in_ == "req":
        raise StataError(100, "in range required")
    if p.using is not None:
        if not spec.using:
            raise StataError(101, "using not allowed")
        fname = strip_outer_quotes(p.using.strip())
        locs["using"] = fname if spec.using_slash else f'using "{fname}"'
    elif spec.using == "req":
        raise StataError(100, "using required")
    if p.weight is not None:
        wtype, wexp = p.weight
        if not spec.weights:
            raise StataError(101, "weights not allowed")
        if wtype not in spec.weights:
            raise StataError(101, f"{wtype} not allowed")
        locs["weight"] = wtype
        locs["exp"] = wexp if spec.weight_slash else f"= {wexp}"

    # opções
    given = parse_options(p.options) if p.options.strip() else []
    star = any(o.kind == "star" for o in spec.options)
    leftover: list[str] = []
    seen: set[str] = set()
    for typed, arg in given:
        low = typed.lower()
        match = None
        for o in spec.options:
            if o.kind == "star":
                continue
            if len(low) >= o.minlen and o.name.startswith(low):
                match = o
                break
        if match is None:
            if star:
                leftover.append(typed + (f"({arg})" if arg is not None else ""))
                continue
            raise StataError(198, f"option {typed} not allowed")
        key = match.name[2:] if match.no else match.name
        locs[key] = _option_value(s, match, arg, typed)
        seen.add(match.name)
    for o in spec.options:
        if o.kind == "star" or o.name in seen:
            continue
        if o.required:
            raise StataError(198, f"option {o.name}{'()' if o.kind != 'flag' else ''} required")
        key = o.name[2:] if o.no else o.name
        if o.default:
            locs[key] = o.default
    if star and leftover:
        locs["options"] = " ".join(leftover)
    return locs


@command("syntax")
def cmd_syntax(s: "Session", args: str) -> None:
    text = s.macros.get_local("0")
    locs = run_syntax(s, args, text)
    # todas as locais que a descrição pode criar começam vazias
    spec = parse_spec(args)
    names = {"varlist", "typlist", "namelist", "anything", "if", "in", "using", "exp", "weight",
             "options"} | {o.name[2:] if o.no else o.name for o in spec.options if o.kind != "star"}
    if spec.main_args:
        nm = _specifiers(spec.main_args).get("name")
        if nm:
            names.add(nm)
    for n in names:
        s.macros.set_local(n, "")
    for k, v in locs.items():
        s.macros.set_local(k, v)


# ---------------------------------------------------------------------------
# gettoken
# ---------------------------------------------------------------------------

@command("gettoken")
def cmd_gettoken(s: "Session", args: str) -> None:
    head, _, opts = args.partition(",")
    m = re.match(r"^\s*([A-Za-z_]\w*)(?:\s+([A-Za-z_]\w*))?\s*:\s*([A-Za-z_]\w*)\s*$", head)
    if not m:
        raise StataError(198, "invalid syntax")
    first, rest_name, source = m.group(1), m.group(2), m.group(3)
    o = {}
    for name, arg in (parse_options(opts) if opts.strip() else []):
        o[name.lower()] = arg if arg is not None else True
    parse_chars = " "
    if "parse" in o or "p" in o:
        parse_chars = strip_outer_quotes(str(o.get("parse", o.get("p"))))
    keep_quotes = "quotes" in o or "q" in o
    match_paren = "match" in o or "m" in o
    text = s.macros.get_local(source)
    token, remainder = _gettoken(text, parse_chars, keep_quotes, match_paren)
    if match_paren:
        target = o.get("match", o.get("m"))
        if isinstance(target, str):
            s.macros.set_local(target, "(" if token.startswith("(") else "")
        if token.startswith("(") and token.endswith(")"):
            token = token[1:-1]
    s.macros.set_local(first, token)
    if rest_name:                 # sem o 2º nome, a fonte não muda
        s.macros.set_local(rest_name, remainder)


def _gettoken(text: str, parse_chars: str, keep_quotes: bool, match_paren: bool) -> tuple[str, str]:
    i, n = 0, len(text)
    if " " in parse_chars:
        while i < n and text[i] == " ":
            i += 1
    if i >= n:
        return "", ""
    if text.startswith('`"', i):
        depth, j = 1, i + 2
        while j < n and depth:
            if text.startswith('`"', j):
                depth += 1
                j += 2
            elif text.startswith("\"'", j):
                depth -= 1
                j += 2
            else:
                j += 1
        tok = text[i:j]
        return (tok if keep_quotes else tok[2:-2]), text[j:]
    if text[i] == '"':
        j = text.find('"', i + 1)
        j = n if j == -1 else j + 1
        tok = text[i:j]
        return (tok if keep_quotes else tok[1:-1]), text[j:]
    if match_paren and text[i] == "(":
        depth, j = 0, i
        while j < n:
            if text[j] == "(":
                depth += 1
            elif text[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        return text[i:j + 1], text[j + 1:]
    if text[i] in parse_chars and text[i] != " ":
        return text[i], text[i + 1:]
    j = i
    while j < n and text[j] not in parse_chars and text[j] != '"':
        j += 1
    return text[i:j], text[j:]


# ---------------------------------------------------------------------------
# marksample / markout
# ---------------------------------------------------------------------------

def _markout(s: "Session", mask: np.ndarray, names: list[str], strok: bool) -> np.ndarray:
    ds = s.data
    for n in names:
        v = ds.get(n)
        if v.is_string:
            if strok:
                mask &= np.array([x != "" for x in v.raw], dtype=bool)
            # VERIFICAR: sem strok, variáveis string são ignoradas
        else:
            mask &= v.data < M.SYSMISS
    return mask


@command("marksample")
def cmd_marksample(s: "Session", args: str) -> None:
    head, _, opts = args.partition(",")
    name = head.strip()
    if not NAME_RE.match(name):
        raise StataError(198, "invalid syntax")
    o = {n.lower() for n, _ in (parse_options(opts) if opts.strip() else [])}
    from ..lang.programs import current_program_scope
    loc = s.macros.get_local
    if_text = loc("if").strip()
    in_text = loc("in").strip()
    p = Parsed(if_=if_text[3:].strip() if if_text.startswith("if ") else (if_text or None),
               in_=in_text[3:].strip() if in_text.startswith("in ") else (in_text or None))
    mask = touse(s, p)
    if "novarlist" not in o and "novar" not in o:
        names = loc("varlist").split()
        mask = _markout(s, mask, names, "strok" in o)
    wexp = loc("exp")
    if loc("weight") and wexp:
        w = s.eval(wexp.lstrip("= ").strip()) if wexp else None
        if isinstance(w, np.ndarray):
            mask &= (w < M.SYSMISS) & ((w > 0) if "zeroweight" not in o else (w >= 0))
    tmp = _temp_name(s)
    s.data.add(Variable(tmp, "byte", mask.astype(np.float64)))
    sc = current_program_scope(s)
    if sc is not None:
        sc.tempvars.append(tmp)
    s.macros.set_local(name, tmp)


@command("markout")
def cmd_markout(s: "Session", args: str) -> None:
    head, _, opts = args.partition(",")
    words = head.split()
    if not words:
        raise StataError(198, "invalid syntax")
    marker = s.data.get(words[0])
    names = s.expand_varlist(" ".join(words[1:])) if len(words) > 1 else []
    o = {n.lower() for n, _ in (parse_options(opts) if opts.strip() else [])}
    mask = marker.data != 0
    mask &= marker.data < M.SYSMISS
    mask = _markout(s, mask, names, "strok" in o)
    s.data.set_numeric(marker, mask.astype(np.float64), promote=False)


def _temp_name(s: "Session") -> str:
    """__000000, __000001, ... (nomes livres)."""
    k = getattr(s, "_temp_counter", 0)
    while True:
        name = f"__{k:06d}"
        k += 1
        if not s.data.has(name) and name not in s.scalars:
            s._temp_counter = k
            return name


