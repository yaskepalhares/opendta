"""Utilitários de dados ([D]): duplicates, isid, levelsof, encode, decode,
destring, tostring, split, mvencode, mvdecode e recode.

Mensagens seguem o Stata 14; as marcadas VERIFICAR têm caso em
compat/do/0306_utilitarios_dados.do.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import sorting
from ..core.dataset import Variable, check_name, smallest_type_for, str_len, str_type_for
from ..core.errors import StataError
from ..core.grouping import group_ids
from ..lang.syntax import Parsed, match_options, parse_standard
from ..lang.words import num_str, parse_numlist, split_words, strip_outer_quotes
from ._util import plural, touse
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _opt_text(opts: dict, key: str) -> str | None:
    v = opts.get(key)
    if v in (None, True):
        return None
    return strip_outer_quotes(str(v).strip())


def _missing_mask(v: Variable) -> np.ndarray:
    if v.is_string:
        return np.array([t == "" for t in v.raw], dtype=bool)
    return v.data >= M.SYSMISS


# ---------------------------------------------------------------------------
# duplicates
# ---------------------------------------------------------------------------

def _dup_header(s: "Session", names: list[str], all_vars: bool) -> None:
    what = "all variables" if all_vars else " ".join(names)
    s.output.write(f"\nDuplicates in terms of {what}\n\n", "text")


@command("duplicates")
def cmd_duplicates(s: "Session", args: str) -> None:
    parts = args.strip().split(None, 1)
    if not parts:
        raise StataError(198, "invalid syntax")
    sub = parts[0].lower()
    rest = parts[1] if len(parts) > 1 else ""
    p = parse_standard(rest)
    ds = s.data
    all_vars = not p.varlist.strip()
    names = list(ds.names) if all_vars else s.expand_varlist(p.varlist)
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    ids, k = group_ids(ds, names, mask)
    sizes = np.bincount(ids[ids >= 0], minlength=k) if k else np.zeros(0, dtype=np.int64)
    copies = np.where(ids >= 0, sizes[np.maximum(ids, 0)] if k else 0, 0)

    if sub in ("report", "r", "rep", "repo", "repor"):
        match_options(p.options, {})
        _dup_header(s, names, all_vars)
        w = s.output.write
        w("-" * 38 + "\n", "text")
        w("   copies | observations       surplus\n", "text")
        w("----------+" + "-" * 27 + "\n", "text")
        used = copies[mask]
        for c in sorted(set(used.tolist())):
            nobs = int((used == c).sum())
            w(f"{c:>9} | ", "text")
            w(f"{nobs:>12,}  {nobs - nobs // c:>12,}\n", "result")
        w("-" * 38 + "\n", "text")
        s.r = {"unique_value": float(k), "N": float(mask.sum())}
    elif sub in ("tag",):
        opts = match_options(p.options, {"generate": 3})
        gen = _opt_text(opts, "generate")
        if not gen:
            raise StataError(198, "option generate() required")
        if ds.has(gen):
            raise StataError(110, f"variable {gen} already defined")
        check_name(gen)
        _dup_header(s, names, all_vars)
        tag = np.where(mask, copies - 1, 0).astype(np.float64)
        # VERIFICAR: tipo da variável criada
        ds.add(Variable(gen, smallest_type_for(tag) if len(tag) else "byte", tag))
    elif sub in ("drop",):
        opts = match_options(p.options, {"force": 5})
        if not all_vars and not opts.get("force"):
            raise StataError(198, "force option required with duplicates drop varlist")  # VERIFICAR
        _dup_header(s, names, all_vars)
        keep = np.ones(ds.nobs, dtype=bool)
        seen = np.zeros(k, dtype=bool)
        for i in np.flatnonzero(ids >= 0):
            g = ids[i]
            if seen[g]:
                keep[i] = False
            seen[g] = True
        removed = ds.keep_obs(keep)
        s.output.write(f"({plural(removed, 'observation')} deleted)\n", "text")
    elif sub in ("list", "l", "li", "lis", "examples", "e", "ex", "exa", "exam", "examp", "exampl", "example"):
        match_options(p.options, {})
        _dup_header(s, names, all_vars)
        _dup_list(s, names, ids, k, sizes, examples=sub.startswith("e"))
    else:
        raise StataError(198, f"invalid subcommand {sub}")   # VERIFICAR
    s.notify_state()


def _dup_list(s: "Session", names: list[str], ids: np.ndarray, k: int, sizes: np.ndarray,
              *, examples: bool) -> None:
    """Tabela no estilo de list com as colunas group:, obs: (ou e.g. obs:)."""
    from .inspect import cell_text
    ds = s.data
    dup_groups = [g for g in range(k) if sizes[g] > 1]
    if not dup_groups:
        s.output.write("(0 observations are duplicates)\n", "text")   # VERIFICAR
        return
    order = np.lexsort((np.arange(ds.nobs), ids))
    order = order[ids[order] >= 0]
    rank = {g: j + 1 for j, g in enumerate(dup_groups)}
    rows: list[list[str]] = []
    group_of_row: list[int] = []
    shown: set[int] = set()
    for i in order:
        g = int(ids[i])
        if g not in rank:
            continue
        if examples and g in shown:
            continue
        shown.add(g)
        head = [str(rank[g])]
        if examples:
            head.append(str(int(sizes[g])))
        head.append(str(i + 1))
        rows.append(head + [cell_text(ds, ds.get(n), int(i)) for n in names])
        group_of_row.append(g)
    # com um só grupo de duplicatas o Stata 14 omite a coluna group: (VERIFICAR
    # com vários grupos); em examples o título é "e.g. obs" sem dois-pontos
    show_group = len(dup_groups) > 1
    if not show_group:
        rows = [r[1:] for r in rows]
    titles = (["group:"] if show_group else []) + (["#"] if examples else []) + \
        (["e.g. obs"] if examples else ["obs:"]) + names
    widths = [max(len(t), *(len(r[j]) for r in rows)) for j, t in enumerate(titles)]
    from .inspect import _left_aligned
    left = [False] * (len(titles) - len(names)) + [_left_aligned(ds.get(n), True) for n in names]

    def fmt(cells: list[str]) -> str:
        return "   ".join(c.ljust(w) if lft else c.rjust(w) for c, w, lft in zip(cells, widths, left))
    inner = len(fmt(titles)) + 2
    w = s.output.write
    w("  +" + "-" * inner + "+\n", "text")
    w("  | " + fmt(titles) + " |\n", "text")
    w("  |" + "-" * inner + "|\n", "text")
    for j, r in enumerate(rows):
        if j and not examples and group_of_row[j] != group_of_row[j - 1]:
            w("  |" + "-" * inner + "|\n", "text")
        w("  | " + fmt(r) + " |\n", "text")
    w("  +" + "-" * inner + "+\n", "text")


# ---------------------------------------------------------------------------
# isid / levelsof
# ---------------------------------------------------------------------------

@command("isid")
def cmd_isid(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"sort": 1, "missok": 5})
    ds = s.data
    if p.using:
        from .combine import load_using
        ds = load_using(p.using)
        from ..core.varlist import expand
        names = expand(ds, p.varlist)
    else:
        names = s.expand_varlist(p.varlist)
    if not names:
        raise StataError(100, "varlist required")
    if not opts.get("missok"):
        for n in names:
            if _missing_mask(ds.get(n)).any():
                raise StataError(459, f"variable {n} should never be missing")
    _, k = group_ids(ds, names)
    if k != ds.nobs:
        if len(names) == 1:
            raise StataError(459, f"variable {names[0]} does not uniquely identify the observations")
        raise StataError(459, f"variables {' '.join(names)} do not uniquely identify the observations")
    if opts.get("sort") and not p.using:
        sorting.sort(ds, names)
        s.notify_state()


@command("levelsof")
def cmd_levelsof(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"clean": 1, "local": 1, "missing": 1, "separate": 1,
                                     "matcell": 4, "matrow": 4, "hexadecimal": 3})
    ds = s.data
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103 if len(names) > 1 else 100,
                         "too many variables specified" if len(names) > 1 else "varlist required")
    v = ds.get(names[0])
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    if not opts.get("missing"):
        mask &= ~_missing_mask(v)
    sep = _opt_text(opts, "separate")
    sep = " " if sep is None else sep
    if v.is_string:
        vals = sorted(set(v.raw[mask].tolist()))
        if opts.get("clean"):
            items = [str(t) for t in vals]
        else:
            items = [f'`"{t}"\'' for t in vals]
    else:
        vals = sorted(set(np.asarray(v.data)[mask].tolist()))
        items = [M.missing_name(x) if x >= M.SYSMISS else num_str(float(x)) for x in vals]
    text = sep.join(items)
    s.output.write(text + "\n", "result")
    s.r = {"levels": text}
    if _opt_text(opts, "local"):
        s.macros.set_local(_opt_text(opts, "local"), text)


# ---------------------------------------------------------------------------
# encode / decode
# ---------------------------------------------------------------------------

@command("encode")
def cmd_encode(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"generate": 3, "label": 5, "noextend": 5})
    ds = s.data
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103, "too many variables specified")
    v = ds.get(names[0])
    if not v.is_string:
        raise StataError(108, "not possible with numeric variable")
    gen = _opt_text(opts, "generate")
    if not gen:
        raise StataError(100, "generate() required")   # VERIFICAR
    if ds.has(gen):
        raise StataError(110, f"variable {gen} already defined")
    check_name(gen)
    lname = _opt_text(opts, "label") or gen
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_)) & ~_missing_mask(v)
    lab = dict(ds.value_labels.get(lname, {}))
    rev = {t: c for c, t in lab.items()}
    new_texts = sorted(set(v.raw[mask].tolist()) - set(rev))
    if new_texts and opts.get("noextend") and lname in ds.value_labels:
        raise StataError(459, f"there are values of {names[0]} not in label {lname}")   # VERIFICAR
    nxt = max(lab, default=0) + 1
    for t in new_texts:
        while nxt in lab:
            nxt += 1
        lab[nxt] = t
        rev[t] = nxt
        nxt += 1
    codes = np.full(ds.nobs, M.SYSMISS)
    for i in np.flatnonzero(mask):
        codes[i] = rev[v.raw[i]]
    ds.value_labels[lname] = lab
    ds.add(Variable(gen, "long", codes, fmt="%8.0g", value_label=lname, label=v.label))  # %8.0g: Stata 14
    s.notify_state()


@command("decode")
def cmd_decode(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"generate": 3, "maxlength": 4})
    ds = s.data
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103, "too many variables specified")
    v = ds.get(names[0])
    if v.is_string:
        raise StataError(108, "not possible with string variable")   # VERIFICAR
    if not v.value_label:
        raise StataError(182, f"{names[0]} not labeled")
    gen = _opt_text(opts, "generate")
    if not gen:
        raise StataError(100, "generate() required")
    if ds.has(gen):
        raise StataError(110, f"variable {gen} already defined")
    check_name(gen)
    lab = ds.value_labels.get(v.value_label, {})
    maxlen = int(_opt_text(opts, "maxlength") or 0) or None
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    data = v.data
    texts = []
    for i in range(ds.nobs):
        x = data[i]
        t = ""
        if mask[i] and x < M.SYSMISS and x == int(x):
            t = lab.get(int(x), "")
        texts.append(t[:maxlen] if maxlen else t)
    longest = max((str_len(t) for t in texts), default=1) or 1
    ds.add(Variable(gen, str_type_for(longest), np.array(texts, dtype=object), label=v.label))
    nmiss = sum(1 for t in texts if t == "")
    if nmiss:
        s.output.write(f"({plural(nmiss, 'missing value')} generated)\n", "text")   # VERIFICAR
    s.notify_state()


# ---------------------------------------------------------------------------
# destring / tostring
# ---------------------------------------------------------------------------

def _gen_targets(s: "Session", names: list[str], opts: dict) -> list[str]:
    ds = s.data
    if opts.get("replace") and opts.get("generate"):
        raise StataError(198, "options generate and replace are mutually exclusive")
    if not opts.get("replace") and not opts.get("generate"):
        raise StataError(198, "must specify either generate or replace option")
    if opts.get("replace"):
        return list(names)
    targets = split_words(_opt_text(opts, "generate") or "")
    if len(targets) != len(names):
        raise StataError(198, "number of variables in generate() must equal number in varlist")  # VERIFICAR
    for t in targets:
        if ds.has(t):
            raise StataError(110, f"variable {t} already defined")
        check_name(t)
    return targets


_NUM_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?$")


def _to_number(t: str, *, dpcomma: bool) -> float | None:
    """Texto → número; '' e '.'/.a-.z → missing; None se não numérico."""
    u = t.strip()
    if u == "":
        return M.SYSMISS
    if re.fullmatch(r"\.[a-z]?", u):
        return M.missing_code(u) if hasattr(M, "from_str") else M.SYSMISS
    if dpcomma:
        u = u.replace(".", "").replace(",", ".")
    if not _NUM_RE.match(u):
        return None
    return float(u.replace("d", "e").replace("D", "e"))


@command("destring")
def cmd_destring(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"generate": 3, "replace": 7, "ignore": 3, "force": 5,
                                     "float": 5, "percent": 3, "dpcomma": 7})
    ds = s.data
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else list(ds.names)
    targets = _gen_targets(s, names, opts)
    ignore = _opt_text(opts, "ignore") or ""
    verb = "replaced" if opts.get("replace") else "generated"
    w = s.output.write
    for name, target in zip(names, targets):
        v = ds.get(name)
        if not v.is_string:
            # VERIFICAR texto
            w(f"{name} already numeric; no {'replace' if opts.get('replace') else 'generate'}\n", "text")
            continue
        removed = set()
        vals = np.full(ds.nobs, M.SYSMISS)
        bad = 0
        for i, t in enumerate(v.raw):
            u = t
            for ch in ignore:
                if ch in u:
                    removed.add(ch)
                    u = u.replace(ch, "")
            if opts.get("percent") and u.strip().endswith("%"):
                u = u.strip()[:-1]
                removed.add("%")
            x = _to_number(u, dpcomma=bool(opts.get("dpcomma")))
            if x is None:
                bad += 1
                continue
            if opts.get("percent") and "%" in t and x < M.SYSMISS:
                x /= 100
            vals[i] = x
        # mensagens observadas no Stata 14 (as marcadas VERIFICAR são deduzidas)
        problem = ("contains characters not specified in ignore()" if ignore
                   else "contains nonnumeric characters")
        if bad and not opts.get("force"):
            w(f"{name} {problem}; no {'replace' if opts.get('replace') else 'generate'}\n", "text")
            continue
        vtype = smallest_type_for(vals)
        if vtype == "double" and opts.get("float"):
            vtype = "float"
        if bad:
            w(f"{name} {problem}; {target} {verb} as {vtype}\n", "text")
        elif removed:
            chars = sorted(removed)
            word = "character" if len(chars) == 1 else "characters"
            w(f"{name}: {word} {' '.join(chars)} removed; {verb} as {vtype}\n", "text")   # VERIFICAR
        else:
            w(f"{name} has all characters numeric; {verb} as {vtype}\n", "text")   # VERIFICAR
        nv = Variable(target, vtype, vals, label=v.label)
        if opts.get("replace"):
            pos = ds.index(name)
            ds.drop_vars([name])
            ds.add(nv, position=pos)
        else:
            ds.add(nv)
        if bad:
            w(f"({plural(bad, 'missing value')} generated)\n", "text")
    s.notify_state()


@command("tostring")
def cmd_tostring(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"generate": 3, "replace": 7, "force": 5, "format": 3,
                                     "usedisplayformat": 4})
    ds = s.data
    names = s.expand_varlist(p.varlist)
    targets = _gen_targets(s, names, opts)
    fmt = _opt_text(opts, "format")
    from ..lang.functions import _string
    w = s.output.write
    for name, target in zip(names, targets):
        v = ds.get(name)
        if v.is_string:
            w(f"{name} already string; no {'replace' if opts.get('replace') else 'generate'}\n", "text")
            continue
        use_fmt = fmt or (v.fmt if opts.get("usedisplayformat") else "%12.0g")   # VERIFICAR padrão
        data = v.data
        texts = [("" if x == M.SYSMISS else M.missing_name(x)) if x >= M.SYSMISS else _string(float(x), use_fmt)
                 for x in data]
        if not opts.get("force") and not fmt:
            back = [(_to_number(t, dpcomma=False) if t else M.SYSMISS) for t in texts]
            lossy = any(b is None or (x < M.SYSMISS and b != x and abs(b - x) > 0)
                        for b, x in zip(back, data))
            if lossy:
                w(f"{name} cannot be converted reversibly; no "
                  f"{'replace' if opts.get('replace') else 'generate'}\n", "text")
                continue
        longest = max((str_len(t) for t in texts), default=1) or 1
        vtype = str_type_for(longest)
        nv = Variable(target, vtype, np.array(texts, dtype=object), label=v.label)
        if opts.get("replace"):
            w(f"{name} was {v.vtype} now {vtype}\n", "text")
            pos = ds.index(name)
            ds.drop_vars([name])
            ds.add(nv, position=pos)
        else:
            w(f"{target} generated as {vtype}\n", "text")
            ds.add(nv)
    s.notify_state()


# ---------------------------------------------------------------------------
# split
# ---------------------------------------------------------------------------

@command("split")
def cmd_split(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"generate": 3, "parse": 1, "limit": 1, "notrim": 5,
                                     "destring": 4, "ignore": 3, "force": 5, "float": 5, "percent": 3})
    ds = s.data
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103, "too many variables specified")
    v = ds.get(names[0])
    if not v.is_string:
        raise StataError(108, "not possible with numeric variable")
    stub = _opt_text(opts, "generate") or names[0]
    parse_raw = opts.get("parse")
    seps = split_words(str(parse_raw), keep_quotes=False) if parse_raw not in (None, True) else [" "]
    seps = [x for x in seps] or [" "]
    limit = int(_opt_text(opts, "limit") or 0) or None
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    rx = re.compile("|".join(re.escape(x) for x in sorted(seps, key=len, reverse=True)))
    pieces: list[list[str]] = []
    for i, t in enumerate(v.raw):
        if not mask[i]:
            pieces.append([])
            continue
        u = t if opts.get("notrim") else t.strip()
        parts = rx.split(u) if u != "" else []
        if seps == [" "]:
            parts = [x for x in parts if x != ""]
        if limit:
            parts = parts[:limit]
        pieces.append(parts)
    n_new = max((len(x) for x in pieces), default=0)
    created = []
    for j in range(n_new):
        nm = f"{stub}{j + 1}"
        if ds.has(nm):
            raise StataError(110, f"variable {nm} already defined")
        check_name(nm)
        created.append(nm)
    for j, nm in enumerate(created):
        texts = [x[j] if j < len(x) else "" for x in pieces]
        longest = max((str_len(t) for t in texts), default=1) or 1
        ds.add(Variable(nm, str_type_for(longest), np.array(texts, dtype=object)))
    w = s.output.write
    if created:
        w("variables created as string: \n", "text")   # VERIFICAR
        w("  ".join(created) + "\n", "result")
    s.r = {"nvars": float(len(created)), "varlist": " ".join(created)}
    if opts.get("destring") and created:
        extra = ", replace"
        for key in ("ignore", "force", "float", "percent"):
            if opts.get(key):
                extra += f" {key}" + (f"({opts[key]})" if opts[key] is not True else "")
        cmd_destring(s, " ".join(created) + extra)
    s.notify_state()


# ---------------------------------------------------------------------------
# mvencode / mvdecode
# ---------------------------------------------------------------------------

def _mv_rules(text: str, *, decode: bool) -> list[tuple[list[float], float]]:
    """mv(# | mvc=# [\\ mvc=# ...]) em mvencode; mv(numlist | numlist=mvc ...) em mvdecode.
    Devolve (valores de origem, valor de destino)."""
    rules = []
    for part in text.split("\\"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            lhs, rhs = (x.strip() for x in part.split("=", 1))
        else:
            lhs, rhs = (part, ".") if decode else ("*", part)
        if decode:
            src = parse_numlist(lhs)
            rules.append((src, M.missing_code(rhs) if rhs.startswith(".") else float(rhs)))
        else:
            src = [M.SYSMISS] if lhs == "." else ([] if lhs == "*" else [M.missing_code(lhs)])
            rules.append((src, float(rhs)))
    return rules


@command("mvdecode")
def cmd_mvdecode(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"mv": 2})
    if opts.get("mv") in (None, True):
        raise StataError(198, "option mv() required")
    rules = _mv_rules(str(opts["mv"]), decode=True)
    ds = s.data
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else list(ds.names)
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    for name in names:
        v = ds.get(name)
        if v.is_string:
            continue
        x = np.array(v.data, dtype=np.float64)
        new = x.copy()
        for src, dst in rules:
            hit = mask & np.isin(x, src)
            new[hit] = dst
        n = int((new != x).sum())
        if n:
            ds.set_numeric(v, new, promote=False)
            s.output.write(f"{name:>12}: {plural(n, 'missing value')} generated\n", "text")
    s.notify_state()


@command("mvencode")
def cmd_mvencode(s: "Session", args: str) -> None:
    p = parse_standard(args)
    opts = match_options(p.options, {"mv": 2, "override": 1})
    if opts.get("mv") in (None, True):
        raise StataError(198, "option mv() required")
    rules = _mv_rules(str(opts["mv"]), decode=False)
    ds = s.data
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else list(ds.names)
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    for name in names:
        v = ds.get(name)
        if v.is_string:
            continue
        x = np.array(v.data, dtype=np.float64)
        targets = [dst for _, dst in rules]
        if not opts.get("override") and np.isin(x[x < M.SYSMISS], targets).any():
            raise StataError(9, f"{name}: already {num_str(targets[0])} in "
                                f"{int(np.isin(x, targets).sum())} observations")   # VERIFICAR
        new = x.copy()
        for src, dst in rules:
            hit = mask & ((x >= M.SYSMISS) if not src else np.isin(x, src))
            hit &= new >= M.SYSMISS
            new[hit] = dst
        n = int((new != x).sum())
        if n:
            ds.set_numeric(v, new)
            s.output.write(f"{name:>12}: {plural(n, 'missing value')} recoded\n", "text")
    s.notify_state()


# ---------------------------------------------------------------------------
# recode
# ---------------------------------------------------------------------------

class _Rule:
    def __init__(self, text: str):
        if "=" not in text:
            raise StataError(198, f"invalid rule ({text})")   # VERIFICAR
        lhs, rhs = text.split("=", 1)
        m = re.match(r'^\s*(\S+)\s*(?:"(.*)"|`"(.*)"\')?\s*$', rhs, re.S)
        if not m:
            raise StataError(198, f"invalid rule ({text})")
        self.target = m.group(1)
        self.label = m.group(2) if m.group(2) is not None else m.group(3)
        self.items = split_words(lhs)
        if not self.items:
            raise StataError(198, f"invalid rule ({text})")

    def value(self) -> float:
        t = self.target
        return M.missing_code(t) if t.startswith(".") else float(t)

    def match(self, x: np.ndarray, lo: float, hi: float) -> tuple[np.ndarray, bool]:
        """(observações que casam, é regra else)."""
        hit = np.zeros(len(x), dtype=bool)
        for it in self.items:
            il = it.lower()
            if il in ("else", "*"):
                return np.ones(len(x), dtype=bool), True
            if il in ("nonmissing", "nonm", "nonmi", "nonmis", "nonmiss", "nonmissi", "nonmissin"):
                hit |= x < M.SYSMISS
            elif il in ("missing", "miss", "missi", "missin"):
                hit |= x >= M.SYSMISS
            elif "/" in it:
                a, b = it.split("/", 1)
                av = lo if a.lower() == "min" else float(a)
                bv = hi if b.lower() == "max" else float(b)
                hit |= (x >= av) & (x <= bv)
            elif il == "min":
                hit |= x == lo
            elif il == "max":
                hit |= x == hi
            elif it.startswith("."):
                hit |= x == M.missing_code(it)
            else:
                try:
                    hit |= x == float(it)
                except ValueError:
                    raise StataError(198, "invalid rule")
        return hit, False


@command("recode")
def cmd_recode(s: "Session", args: str) -> None:
    text = args.strip()
    # varlist (regra) [(regra) ...] [if] [in] [, opções]
    start = text.find("(")
    if start == -1:
        raise StataError(198, "invalid syntax")
    varlist_text = text[:start]
    rule_texts: list[str] = []
    i = start
    while i < len(text) and text[i] == "(":
        depth, j = 0, i
        while j < len(text):
            if text[j] == '"':
                q = text.find('"', j + 1)
                j = len(text) - 1 if q == -1 else q
            elif text[j] == "(":
                depth += 1
            elif text[j] == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth:
            raise StataError(198, "invalid syntax")
        rule_texts.append(text[i + 1:j])
        i = j + 1
        while i < len(text) and text[i] in " \t":
            i += 1
    tail = text[i:]
    rules = [_Rule(r) for r in rule_texts]
    p = parse_standard(tail)
    if p.varlist.strip():
        raise StataError(198, "invalid syntax")
    opts = match_options(p.options, {"generate": 3, "prefix": 3, "label": 3, "copyrest": 4, "test": 4})
    ds = s.data
    names = s.expand_varlist(varlist_text)
    if not names:
        raise StataError(100, "varlist required")
    gens = split_words(_opt_text(opts, "generate") or "")
    prefix = _opt_text(opts, "prefix")
    if gens and len(gens) != len(names):
        raise StataError(198, "generate() must contain as many names as varlist")   # VERIFICAR
    targets = gens or ([prefix + n for n in names] if prefix else names)
    new_var = bool(gens or prefix)
    for t in (targets if new_var else []):
        if ds.has(t):
            raise StataError(110, f"variable {t} already defined")
        check_name(t)
    mask = touse(s, Parsed(if_=p.if_, in_=p.in_))
    labelled = [r for r in rules if r.label is not None]
    w = s.output.write
    for name, target in zip(names, targets):
        v = ds.get(name)
        if v.is_string:
            raise StataError(109, "type mismatch")   # VERIFICAR: recode só aceita numéricas
        x = np.array(v.data, dtype=np.float64)
        ok = x < M.SYSMISS
        lo = float(x[ok & mask].min()) if (ok & mask).any() else M.SYSMISS
        hi = float(x[ok & mask].max()) if (ok & mask).any() else M.SYSMISS
        new = x.copy()
        done = np.zeros(len(x), dtype=bool)
        for r in rules:
            hit, _ = r.match(x, lo, hi)
            hit &= mask & ~done
            new[hit] = r.value()
            done |= hit
        if new_var and not opts.get("copyrest"):
            new[~mask] = M.SYSMISS
        vtype = v.vtype
        need = smallest_type_for(new)
        from ..core.dataset import TYPE_ORDER
        if TYPE_ORDER[need] > TYPE_ORDER[vtype]:
            vtype = need if not (vtype == "float" and need == "double") else "float"
        # VERIFICAR: com generate(), conta só as diferenças dentro da amostra
        changes = int((new != x).sum()) if not new_var else int(((new != x) & mask).sum())
        lname = None
        if labelled:
            lname = _opt_text(opts, "label") or target
            lab = dict(ds.value_labels.get(lname, {}))
            for r in labelled:
                val = r.value()
                if val < M.SYSMISS:
                    lab[int(val)] = r.label
            ds.value_labels[lname] = lab
        if new_var:
            nv = Variable(target, vtype, new, label=f"RECODE of {name}"
                          + (f" ({v.label})" if v.label else ""),   # VERIFICAR rótulo
                          value_label=lname or "")
            ds.add(nv)
            w(f"({plural(changes, 'difference')} between {name} and {target})\n", "text")
        else:
            if vtype != v.vtype:
                v.vtype = vtype
            ds.set_numeric(v, new, promote=False)
            if lname:
                v.value_label = lname
            w(f"({name}: {plural(changes, 'change')} made)\n", "text")
    s.notify_state()

