"""local, global, macro e funções estendidas de macro (manual [P] macro)."""

from __future__ import annotations

import fnmatch
import re
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.formats import number_to_macro
from ..lang.words import split_words, strip_outer_quotes
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_NAME = re.compile(r"^\s*([A-Za-z_0-9][A-Za-z0-9_]*)(.*)$", re.S)


def _assign(s: "Session", args: str, setter) -> None:
    inc = re.match(r"^\s*(\+\+|--)([A-Za-z_][A-Za-z0-9_]*)\s*$", args)
    if inc and setter == s.macros.set_local:
        # local ++i / local --i
        name = inc.group(2)
        current = s.macros.get_local(name) or "0"
        try:
            x = float(current)
        except ValueError:
            raise StataError(198, "invalid syntax")
        x += 1 if inc.group(1) == "++" else -1
        setter(name, number_to_macro(x))
        return
    m = _NAME.match(args)
    if not m:
        raise StataError(198, "invalid syntax")
    name, rest = m.group(1), m.group(2)
    if len(name) > 31:
        raise StataError(198, f"{name} invalid name")
    rest_s = rest.strip()
    if rest_s.startswith("="):
        v = s.eval(rest_s[1:])
        value = v if isinstance(v, str) else number_to_macro(v)
    elif rest_s.startswith(":"):
        value = extended_function(s, rest_s[1:].strip())
    else:
        value = strip_outer_quotes(rest_s) if rest_s else ""
    setter(name, value)


@command("local", "loc")
def cmd_local(s: "Session", args: str) -> None:
    _assign(s, args, s.macros.set_local)


@command("global", "gl")
def cmd_global(s: "Session", args: str) -> None:
    _assign(s, args, s.macros.set_global)


@command("macro", "ma")
def cmd_macro(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    rest = rest.strip()
    out = s.output
    if sub in ("list", "l", "li", "lis", "dir", "di"):
        names = rest.split()
        for name, value in sorted(s.macros.globals.items()):
            if names and name not in names:
                continue
            out.write(f"{name + ':':<16}{value}\n", "text")
        for name, value in sorted(s.macros.locals.items()):
            if names and "_" + name not in names:
                continue
            out.write(f"{'_' + name + ':':<16}{value}\n", "text")
        return
    if sub == "drop":
        for pat in rest.split():
            if pat == "_all":
                s.macros.clear_globals()
                continue
            target = s.macros.locals if pat.startswith("_") else s.macros.globals
            key = pat[1:] if pat.startswith("_") else pat
            for name in [k for k in target if fnmatch.fnmatchcase(k, key)]:
                del target[name]
        return
    if sub == "shift":
        loc = s.macros.locals
        k = 1
        while str(k + 1) in loc:
            loc[str(k)] = loc[str(k + 1)]
            k += 1
        loc.pop(str(k), None)
        return
    raise StataError(198, "invalid syntax")


# ---------------------------------------------------------------------------
# Funções estendidas: local x : <função>
# ---------------------------------------------------------------------------

def extended_function(s: "Session", text: str) -> str:
    t = text.strip()

    m = re.match(r"^word\s+count\s+(.*)$", t, re.S)
    if m:
        return str(len(split_words(m.group(1))))
    m = re.match(r"^word\s+(\S+)\s+of\s+(.*)$", t, re.S)
    if m:
        k = int(float(s.eval(m.group(1))))
        words = split_words(m.group(2))
        return words[k - 1] if 1 <= k <= len(words) else ""
    m = re.match(r"^(?:length|strlen)\s+(local|global)\s+(\S+)$", t)
    if m:
        store = s.macros.get_local if m.group(1) == "local" else s.macros.get_global
        return str(len(store(m.group(2))))
    m = re.match(r"^di(?:s|sp|spl|spla|splay)?\s+(.*)$", t, re.S)
    if m:
        return _capture_display(s, m.group(1))
    m = re.match(r"^list\s+(.*)$", t, re.S)
    if m:
        return _list_function(s, m.group(1).strip())
    m = re.match(r"^env(?:ironment)?\s+(\S+)$", t)
    if m:
        import os
        return os.environ.get(m.group(1), "")
    m = re.match(r"^piece\s+(\S+)\s+(\S+)\s+of\s+(.*)$", t, re.S)
    if m:
        return _piece(int(float(m.group(1))), int(float(m.group(2))), strip_outer_quotes(m.group(3)))
    value = _data_function(s, t)
    if value is not None:
        return value
    word = t.split()[0] if t else ""
    # VERIFICAR: mensagem do Stata para função estendida desconhecida
    raise StataError(198, f"{word} not allowed")


def _var(s: "Session", name: str):
    from ..core.varlist import resolve_name
    return s.data.get(resolve_name(s.data, name.strip()))


def _data_function(s: "Session", t: str) -> str | None:
    """Funções que consultam os dados, rótulos, características, arquivos e
    resultados guardados."""
    import fnmatch
    from pathlib import Path

    m = re.fullmatch(r"(type|format|f|fo|for|form|forma)\s+(\S+)", t)
    if m:
        v = _var(s, m.group(2))
        return v.vtype if m.group(1) == "type" else v.fmt
    m = re.fullmatch(r"value\s+l(?:a|ab|abe|abel)?\s+(\S+)", t)
    if m:
        return _var(s, m.group(1)).value_label
    m = re.fullmatch(r"var(?:iable)?\s+l(?:a|ab|abe|abel)?\s+(\S+)", t)
    if m:
        return _var(s, m.group(1)).label
    if re.fullmatch(r"data\s+l(?:a|ab|abe|abel)?", t):
        return s.data.label
    if t == "sortedby":
        return " ".join(s.data.sortlist)
    m = re.fullmatch(r"label\s+\((\S+)\)\s+(\S+)(?:\s+(\d+))?(\s*,\s*strict)?", t)
    if m:
        v = _var(s, m.group(1))
        return _label_text(s, v.value_label, m.group(2), m.group(3), bool(m.group(4)))
    m = re.fullmatch(r"label\s+(\S+)\s+(\S+)(?:\s+(\d+))?(\s*,\s*strict)?", t)
    if m:
        if m.group(1) not in s.data.value_labels and m.group(2) != "maxlength":
            raise StataError(111, f"value label {m.group(1)} not found")
        return _label_text(s, m.group(1), m.group(2), m.group(3), bool(m.group(4)))
    m = re.fullmatch(r"char\s+(\S+)\[([^\]]*)\]", t)
    if m:
        owner = m.group(1) if m.group(1) == "_dta" else _var(s, m.group(1)).name
        chars = s.data.chars.get(owner, {})
        if m.group(2) == "":
            return " ".join(chars)
        return chars.get(m.group(2), "")
    m = re.fullmatch(r'dir\s+(\S+|"[^"]*")\s+(files|dirs|other)\s+(\S+|"[^"]*")(\s*,.*)?', t)
    if m:
        folder = Path(strip_outer_quotes(m.group(1))).expanduser()
        pattern = strip_outer_quotes(m.group(3))
        opts = (m.group(4) or "").replace(",", " ").split()
        if not folder.is_dir():
            if "nofail" in opts:
                return ""
            raise StataError(601, f"directory {folder} not found")   # VERIFICAR
        kind = m.group(2)
        respect = "respectcase" in opts
        out = []
        for p in sorted(folder.iterdir(), key=lambda x: x.name):
            if kind == "files":
                ok = p.is_file()
            elif kind == "dirs":
                ok = p.is_dir()
            else:
                ok = not (p.is_file() or p.is_dir())
            name = p.name
            if name.startswith(".") and not pattern.startswith("."):
                continue                       # ocultos só com padrão explícito
            hit = (fnmatch.fnmatchcase(name, pattern) if respect
                   else fnmatch.fnmatchcase(name.lower(), pattern.lower()))
            if ok and hit:
                out.append(f'"{name}"')
        return " ".join(out)
    m = re.fullmatch(r"sysdir\s+(\S+)", t)
    if m:
        from ..lang.adopath import sysdir
        return sysdir(s).get(m.group(1).upper(), m.group(1))
    m = re.fullmatch(r"permname\s+(\S+)(?:\s*,\s*length\((\d+)\))?", t)
    if m:
        base = re.sub(r"[^A-Za-z0-9_]", "_", m.group(1))[: int(m.group(2) or 32)]
        name, k = base, 0
        while s.data.has(name):
            k += 1
            name = f"{base[: 32 - len(str(k))]}{k}"
        return name
    m = re.fullmatch(r"(copy|strlen|length|ustrlen|udstrlen)\s+(local|global)\s+(\S+)", t)
    if m:
        store = s.macros.get_local if m.group(2) == "local" else s.macros.get_global
        text = store(m.group(3))
        if m.group(1) == "copy":
            return text
        return str(len(text.encode("utf-8")) if m.group(1) in ("strlen", "length") else len(text))
    m = re.fullmatch(r'subinstr\s+(local|global)\s+(\S+)\s+("[^"]*"|`"[^`]*"\'|\S+)\s+'
                     r'("[^"]*"|`"[^`]*"\'|\S+)(\s*,.*)?', t, re.S)
    if m:
        store = s.macros.get_local if m.group(1) == "local" else s.macros.get_global
        text = store(m.group(2))
        old, new = strip_outer_quotes(m.group(3)), strip_outer_quotes(m.group(4))
        opts = m.group(5) or ""
        count_m = re.search(r"count\((local|global)\s+(\w+)\)", opts)
        all_ = re.search(r"\ball\b", opts) is not None
        word = re.search(r"\bword\b", opts) is not None
        if word:
            words = split_words(text, keep_quotes=True)
            n = 0
            outw = []
            for w in words:
                if w == old and (all_ or n == 0):
                    n += 1
                    if new:
                        outw.append(new)
                else:
                    outw.append(w)
            result = " ".join(outw)
        else:
            n = text.count(old) if all_ else min(1, text.count(old))
            result = text.replace(old, new) if all_ else text.replace(old, new, 1)
        if count_m:
            setter = s.macros.set_local if count_m.group(1) == "local" else s.macros.set_global
            setter(count_m.group(2), str(n))
        return result
    m = re.fullmatch(r"([rse])\((scalars|macros|matrices|functions)\)", t)
    if m:
        store = {"r": s.r, "e": s.e, "s": s.sret}[m.group(1)]
        kind = m.group(2)
        if kind == "scalars":
            names = [k for k, v in store.items() if not isinstance(v, str) and not hasattr(v, "shape")]
        elif kind == "macros":
            names = [k for k, v in store.items() if isinstance(v, str)]
        elif kind == "matrices":
            names = [k for k, v in store.items() if hasattr(v, "rows")]
        else:
            names = []
        return " ".join(sorted(names, reverse=False))
    m = re.fullmatch(r"(rownames|colnames|rowfullnames|colfullnames|roweq|coleq)\s+(\S+)", t)
    if m:
        from .matrix import matrix_names
        return matrix_names(s, m.group(2), m.group(1))
    m = re.fullmatch(r"(rowsof|colsof)\s+(\S+)", t)
    if m:
        from .matrix import get_matrix
        mat = get_matrix(s, m.group(2))
        return str(mat.rows if m.group(1) == "rowsof" else mat.cols)
    return None


def _label_text(s: "Session", lbl: str, value: str, length: str | None, strict: bool) -> str:
    labels = s.data.value_labels.get(lbl, {})
    if value == "maxlength":
        return str(max((len(x) for x in labels.values()), default=0))
    try:
        x = float(value)
    except ValueError:
        raise StataError(198, "invalid syntax")
    text = labels.get(int(x)) if x == int(x) else None
    if text is None:
        text = "" if strict else number_to_macro(x)
    if length:
        text = text[: int(length)]
    return text


def _capture_display(s: "Session", args: str) -> str:
    from ..core.output import Capture
    from .display import cmd_display
    cap = Capture()
    out = s.output
    saved = out._listeners, out.column, out.quiet_depth, out.noisy_depth
    out._listeners = [cap]
    out.quiet_depth = out.noisy_depth = 0
    try:
        cmd_display(s, args + " _continue")
    finally:
        out._listeners, out.column, out.quiet_depth, out.noisy_depth = saved
    return cap.text


def _piece(n: int, length: int, text: str) -> str:
    words = text.split()
    pieces: list[str] = []
    cur = ""
    for w in words:
        cand = (cur + " " + w).strip()
        if len(cand) <= length or not cur:
            cur = cand
        else:
            pieces.append(cur)
            cur = w
    if cur:
        pieces.append(cur)
    return pieces[n - 1] if 1 <= n <= len(pieces) else ""


def _list_function(s: "Session", spec: str) -> str:
    get = s.macros.get_local
    words = spec.split()
    if not words:
        raise StataError(198, "invalid syntax")
    op = words[0]
    if op in ("uniq", "dups", "sort", "clean", "retokenize", "sizeof"):
        items = split_words(get(words[1]), keep_quotes=True)
        if op == "uniq":
            seen: list[str] = []
            for it in items:
                if it not in seen:
                    seen.append(it)
            return " ".join(seen)
        if op == "dups":
            seen, dups = [], []
            for it in items:
                if it in seen:
                    dups.append(it)
                seen.append(it)
            return " ".join(dups)
        if op == "sort":
            return " ".join(sorted(items))
        if op in ("clean", "retokenize"):
            return " ".join(items)
        return str(len(items))
    if op == "posof":
        m = re.match(r'^posof\s+(.+?)\s+in\s+(\S+)$', spec)
        if not m:
            raise StataError(198, "invalid syntax")
        target = strip_outer_quotes(m.group(1))
        items = split_words(get(m.group(2)))
        return str(items.index(target) + 1) if target in items else "0"
    if len(words) == 3 and words[1] in ("|", "&", "-", "==", "===", "in"):
        a = split_words(get(words[0]), keep_quotes=True)
        b = split_words(get(words[2]), keep_quotes=True)
        o = words[1]
        if o == "|":
            return " ".join(a + [x for x in b if x not in a])
        if o == "&":
            return " ".join(x for x in a if x in b)
        if o == "-":
            return " ".join(x for x in a if x not in b)
        if o == "==":
            return "1" if sorted(a) == sorted(b) else "0"
        if o == "===":
            return "1" if a == b else "0"
        return "1" if all(x in b for x in a) else "0"
    raise StataError(198, "invalid syntax")
