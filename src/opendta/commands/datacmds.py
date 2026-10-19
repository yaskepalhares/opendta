"""Data-management commands of [D]: clonevar, ds, lookfor, inspect, labelbook,
compare, cf, datasignature, insobs, range, ipolate, separate, stack, xpose,
expandcl, drawnorm, corr2data, outfile, putmata and statsby. Written from the
manual's descriptions; layouts not shown there are marked VERIFY and checked in
compat/do/1501_datamgmt.do."""

from __future__ import annotations

import datetime as _dt
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Dataset, Variable
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.varlist import expand, unique
from ..lang.syntax import find_top, match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import eval_vector, touse
from .data import _note
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _missing_text(n: int) -> str:
    return f"({n:,} missing value{'s' if n != 1 else ''} generated)"


def _need_clear(s: "Session", o: dict) -> None:
    ds = s.data
    if ds.changed and ds.nvars and not o.get("clear"):
        raise StataError(4, "no; data in memory would be lost")


def _new_numeric(s: "Session", name: str, values: np.ndarray, vtype: str | None = None) -> Variable:
    return Variable(name, vtype or s.default_type(), np.asarray(values, dtype=float))


# ---------------------------------------------------------------------------
# clonevar
# ---------------------------------------------------------------------------

@command("clonevar", "clonevar")
def cmd_clonevar(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if not re.fullmatch(r"\s*\w+\s*", p.varlist) or not p.exp or not re.fullmatch(r"\s*\w+\s*", p.exp):
        raise StataError(198, "invalid syntax")
    new, old = p.varlist.strip(), s.expand_varlist(p.exp.strip())[0]
    ds = s.data
    src = ds.get(old)
    mask = touse(s, p)
    if src.is_string:
        vals = np.where(mask, np.asarray(src.raw, dtype=object), "")
    else:
        vals = np.where(mask, src.data, SYS)
    v = Variable(new, src.vtype, vals, fmt=src.fmt, label=src.label, value_label=src.value_label,
                 notes=list(src.notes))
    ds.add(v, position=None)
    if old in ds.chars:
        ds.chars[new] = dict(ds.chars[old])
    nmiss = int((~mask).sum()) if not src.is_string else 0
    if p.if_ or p.in_:
        if nmiss:
            _note(s, _missing_text(nmiss))
    s.notify_state()


# ---------------------------------------------------------------------------
# ds, lookfor
# ---------------------------------------------------------------------------

def _spec_match(s: "Session", v: Variable, spec: str) -> bool:
    words = spec.split()
    if not words:
        raise StataError(198, "invalid syntax")
    kind, pats = words[0], words[1:]
    ds = s.data

    def pat_ok(text: str) -> bool:
        if not pats:
            return bool(text)
        from fnmatch import fnmatchcase
        return any(fnmatchcase(text, p) for p in pats)
    if kind == "type":
        t = v.vtype
        for p in pats:
            if p == "numeric" and not v.is_string or p == "string" and v.is_string or p == t:
                return True
            if p.startswith("str") and "-" in p and t.startswith("str") and t != "strL":
                lo, hi = p[3:].split("-")
                if int(lo) <= int(t[3:]) <= int(hi):
                    return True
        return False
    if kind == "varlabel":
        return pat_ok(v.label)
    if kind == "vallabel":
        return pat_ok(v.value_label)
    if kind == "format":
        return pat_ok(v.fmt)
    if kind == "char":
        names = list(ds.chars.get(v.name, {}))
        return bool(names) if not pats else any(pat_ok(n) for n in names)
    raise StataError(198, f"{kind} invalid")   # VERIFY


@command("ds", "ds")
def cmd_ds(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"not": 3, "alpha": 1, "detail": 1, "varwidth": 4, "skip": 2,
                                  "has": 3, "insensitive": 3, "indent": 3})
    ds = s.data
    names = unique(expand(ds, p.varlist)) if p.varlist.strip() else ds.names
    has, nothas = None, None
    if isinstance(o.get("has"), str):
        has = str(o["has"])
    if isinstance(o.get("not"), str):
        nothas = str(o["not"])
    elif o.get("not") is True:
        if not p.varlist.strip():
            raise StataError(198, "varlist required with not")   # VERIFY
        drop = set(names)
        names = [n for n in ds.names if n not in drop]
    if has is not None:
        names = [n for n in names if _spec_match(s, ds.get(n), has)]
    if nothas is not None:
        names = [n for n in names if not _spec_match(s, ds.get(n), nothas)]
    if o.get("alpha"):
        names = sorted(names, key=lambda n: (n.lower(), n) if o.get("insensitive") else n)
    if o.get("detail"):
        if names:
            from .inspect import cmd_describe
            cmd_describe(s, " ".join(names))
    elif names:
        from .summarize import _abbrev
        w = int(o["varwidth"]) if isinstance(o.get("varwidth"), str) else 12   # VERIFY default
        skip = int(o["skip"]) if isinstance(o.get("skip"), str) else 2
        indent = int(o["indent"]) if isinstance(o.get("indent"), str) else 0
        col = w + skip
        per = max(1, (79 - indent) // col)
        shown = [_abbrev(n, w) for n in names]
        out = s.output
        for i in range(0, len(shown), per):
            line = " " * indent + "".join(f"{n:<{col}}" for n in shown[i:i + per])
            out.write(line.rstrip() + "\n", "text")
    s.r = {"varlist": " ".join(names)}


@command("lookfor", "lookfor")
def cmd_lookfor(s: "Session", args: str) -> None:
    words = [strip_outer_quotes(w) for w in re.findall(r'"[^"]*"|\S+', args)]
    if not words:
        raise StataError(198, "invalid syntax")
    ds = s.data
    found = [v.name for v in ds.vars
             if any(w.lower() in v.name.lower() or w.lower() in v.label.lower() for w in words)]
    if found:
        from .inspect import cmd_describe
        cmd_describe(s, " ".join(found))
    s.r = {"varlist": " ".join(found)}


# ---------------------------------------------------------------------------
# inspect
# ---------------------------------------------------------------------------

@command("inspect", "inspect", byable=True)
def cmd_inspect(s: "Session", args: str) -> None:
    p = parse_standard(args)
    ds = s.data
    names = unique(expand(ds, p.varlist)) if p.varlist.strip() else ds.names
    mask = touse(s, p)
    out = s.output
    for n in names:
        v = ds.get(n)
        out.write("\n", "text")
        if v.is_string:
            out.write(f"{n}:  {v.label}\n", "text")
            out.write("-" * min(len(n) + 1, 30) + "\n", "text")
            out.write("    (string variable)\n", "text")   # VERIFY
            continue
        x = v.data[mask]
        ok = x < SYS
        xv = x[ok]
        is_int = xv == np.floor(xv)
        rows = [("Negative", xv < 0), ("Zero", xv == 0), ("Positive", xv > 0)]
        hist = _hist_rows(xv)
        out.write(f"{n + ':':<{max(len(n) + 3, 8)}}{v.label:<{max(0, 44 - max(len(n) + 3, 8))}}Number of Observations\n", "text")
        out.write(f"{'-' * min(len(n) + 1, 30):<44}------------------------------\n", "text")
        out.write(" " * 47 + "Total      Integers   Nonintegers\n", "text")
        for i, (lab, sel) in enumerate(rows):
            t, ints = int(sel.sum()), int((sel & is_int).sum())
            out.write(f"|  {hist[i]:<24}{lab:<14}", "text")
            out.write(f"{_dash(t):>8}  {_dash(ints):>12}  {_dash(t - ints):>12}\n", "result")
        out.write(f"|  {hist[3]:<38}{'-' * 8:>8}  {'-' * 12:>12}  {'-' * 12:>12}\n", "text")
        t, ints = len(xv), int(is_int.sum())
        out.write(f"|  {hist[4]:<24}{'Total':<14}", "text")
        out.write(f"{_dash(t):>8}  {_dash(ints):>12}  {_dash(t - ints):>12}\n", "result")
        out.write(f"|  {hist[5]:<24}{'Missing':<14}", "text")
        out.write(f"{_dash(int((~ok).sum())):>8}\n", "result")
        out.write("+" + "-" * 22 + " " * 17 + "-" * 10 + "\n", "text")
        lo = format_value(float(xv.min()), "%9.0g", pad=False).strip() if len(xv) else ""
        hi = format_value(float(xv.max()), "%9.0g", pad=False).strip() if len(xv) else ""
        out.write(f"{lo:<12}{hi:>11}{len(x):>23}\n", "result")
        u = len(np.unique(xv))
        out.write(f"   ({u if u <= 99 else 'More than 99'} unique value{'s' if u != 1 else ''})\n", "text")
        if v.value_label and v.value_label in ds.value_labels:
            labs = ds.value_labels[v.value_label]
            nlab = sum(1 for val in np.unique(xv) if val == int(val) and int(val) in labs)
            out.write(f"\n   {n} is labeled and {'all' if nlab == u else nlab} values are documented in the "
                      f"label.\n", "text")   # VERIFY wording
        s.r = {"N": float(len(x)), "N_neg": float((xv < 0).sum()), "N_0": float((xv == 0).sum()),
               "N_pos": float((xv > 0).sum()), "N_negint": float(((xv < 0) & is_int).sum()),
               "N_posint": float(((xv > 0) & is_int).sum()), "N_unique": float(u)}


def _dash(n: int) -> str:
    return "-" if n == 0 else f"{n:,}"


def _hist_rows(x: np.ndarray) -> list[str]:
    """Small text histogram: up to 5 columns of '#', 6 rows high (VERIFY: Stata's layout)."""
    rows = [""] * 6
    if not len(x):
        return rows
    k = min(5, len(np.unique(x)))
    counts, _ = np.histogram(x, bins=k)
    top = counts.max()
    for r in range(6):
        level = (5 - r) / 6 * top
        rows[r] = "   ".join("#" if c > level else ("." if c > 0 and r == 5 else " ") for c in counts).rstrip()
    return rows


# ---------------------------------------------------------------------------
# labelbook
# ---------------------------------------------------------------------------

@command("labelbook", "labelbook")
def cmd_labelbook(s: "Session", args: str) -> None:
    words, opt = _split_opts(args)
    o = match_options(opt, {"alpha": 1, "length": 3, "list": 2, "problems": 1, "detail": 1})
    ds = s.data
    names = words or sorted(ds.value_labels)
    if o.get("alpha"):
        names = sorted(names)
    L = int(o["length"]) if isinstance(o.get("length"), str) else 12
    out = s.output
    problems = []
    for name in names:
        if name not in ds.value_labels:
            raise StataError(111, f"value label {name} not found")
        lab = ds.value_labels[name]
        vals = sorted(lab)
        texts = [lab[v] for v in vals]
        users = [v.name for v in ds.vars if v.value_label == name]
        real = [v for v in vals if v < 2 ** 31 - 1 - 26]     # extended missing codes are stored high
        gaps = bool(real) and (real[-1] - real[0] + 1 != len(real))
        lens = [len(t) for t in texts] or [0]
        uniq_full = len(set(texts)) == len(texts)
        uniq_L = len({t[:L] for t in texts}) == len(texts)
        null = any(t == "" for t in texts)
        blanks = any(t != t.strip() for t in texts)
        numnum = any(re.fullmatch(r"\s*-?[\d.]+\s*", t) and float(t) != v for v, t in zip(vals, texts)
                     if re.fullmatch(r"\s*-?\d+(\.\d+)?\s*", t))
        if not users or not uniq_L or null or blanks or numnum or gaps and o.get("problems"):
            problems.append(name)
        if o.get("problems"):
            continue
        out.write("\n" + "-" * 79 + "\n", "text")
        out.write(f"value label {name}\n", "text")
        out.write("-" * 79 + "\n\n", "text")
        rng = f"[{real[0]},{real[-1]}]" if real else "[.]"
        rows = [("range:", rng, "string length:", f"[{min(lens)},{max(lens)}]"),
                ("N:", str(len(vals)), "unique at full length:", "yes" if uniq_full else "no"),
                ("gaps:", "yes" if gaps else "no", f"unique at length {L}:", "yes" if uniq_L else "no"),
                ("missing .*:", str(len(vals) - len(real)), "null string:", "yes" if null else "no"),
                ("", "", "leading/trailing blanks:", "yes" if blanks else "no"),
                ("", "", "numeric -> numeric:", "yes" if numnum else "no")]
        out.write(f"{'values':>12}{'labels':>42}\n", "text")
        for a, b, c, d in rows:
            out.write(f"{a:>12}  {b:<16}{c:>25}  {d}\n", "text")
        limit = int(o["list"]) if isinstance(o.get("list"), str) else 32000
        out.write("\n  definition\n", "text")
        for v, t in list(zip(vals, texts))[:limit]:
            out.write(f"{v:>12}   {t}\n", "text")
        out.write("\n" + f"{'variables:':>12}  {' '.join(users)}\n", "text")
    if o.get("problems") or problems:
        if problems and o.get("problems"):
            out.write("\nValue labels with potential problems: " + " ".join(problems) + "\n", "text")   # VERIFY
    s.r = {"names": " ".join(names), "notused": " ".join(n for n in names if not any(v.value_label == n for v in ds.vars))}


def _split_opts(args: str) -> tuple[list[str], str]:
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    return body.split(), opt


# ---------------------------------------------------------------------------
# compare
# ---------------------------------------------------------------------------

@command("compare", "compare", byable=True)
def cmd_compare(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist)
    if len(names) != 2:
        raise StataError(198, "invalid syntax")
    ds = s.data
    a, b = ds.get(names[0]), ds.get(names[1])
    if a.is_string or b.is_string:
        raise StataError(109, "type mismatch")
    mask = touse(s, p)
    x, y = a.data[mask], b.data[mask]
    okx, oky = x < SYS, y < SYS
    both = okx & oky
    d = x - y
    out = s.output
    out.write("\n" + " " * 40 + "---------- difference ----------\n", "text")
    out.write(" " * 28 + "count       minimum      average     maximum\n", "text")
    out.write("-" * 72 + "\n", "text")

    def line(label, sel, diff=True):
        n = int(sel.sum())
        if not n:
            return
        out.write(f"{label:<23}", "text")
        if diff:
            dd = d[sel]
            out.write(f"{n:>10}{_g(dd.min()):>14}{_g(dd.mean()):>13}{_g(dd.max()):>12}\n", "result")
        else:
            out.write(f"{n:>10}\n", "result")
    an, bn = a.name, b.name
    line(f"{an}<{bn}", both & (x < y))
    line(f"{an}={bn}", both & (x == y))
    line(f"{an}>{bn}", both & (x > y))
    out.write(" " * 23 + "-" * 10 + "\n", "text")
    line("jointly defined", both)
    line(f"{bn} missing only", okx & ~oky, False)
    line(f"{an} missing only", ~okx & oky, False)
    line("jointly missing", ~okx & ~oky, False)
    out.write(" " * 23 + "-" * 10 + "\n", "text")
    out.write(f"{'total':<23}", "text")
    out.write(f"{len(x):>10}\n", "result")


def _g(x: float) -> str:
    return format_value(float(x), "%9.0g", pad=False).strip()


# ---------------------------------------------------------------------------
# cf
# ---------------------------------------------------------------------------

@command("cf", "cf")
def cmd_cf(s: "Session", args: str) -> None:
    from .combine import load_using
    p = parse_standard(args)
    if p.using is None:
        raise StataError(100, "using required")
    o = match_options(p.options, {"all": 3, "verbose": 1})
    ds = s.data
    U = load_using(p.using)
    names = unique(expand(ds, p.varlist)) if p.varlist.strip() else ds.names
    out = s.output
    bad = 0
    if ds.nobs != U.nobs:
        out.write(f"{'':>13}master has {ds.nobs:,} obs.; using {U.nobs:,}\n", "text")   # VERIFY
        raise StataError(9, "")
    for n in names:
        a = ds.get(n)
        if not U.has(n):
            out.write(f"{n:>13}:  does not exist in using\n", "text")
            bad += 1
            continue
        b = U.get(n)
        if a.is_string != b.is_string:
            diff = np.ones(ds.nobs, dtype=bool)
        elif a.is_string:
            diff = np.array([x != y for x, y in zip(a.raw, b.raw)], dtype=bool)
        else:
            diff = ~((a.data == b.data) | ((a.data >= SYS) & (b.data >= SYS) & (a.data == b.data)))
        k = int(diff.sum())
        if k:
            bad += 1
            out.write(f"{n:>13}:  {k} mismatch{'es' if k != 1 else ''}\n", "text")
            if o.get("verbose"):
                for i in np.flatnonzero(diff):
                    out.write(f"{'':>15}obs {i + 1:,}. {a.value(i)} in master; {b.value(i)} in using\n", "text")
        elif o.get("all") or o.get("verbose"):
            out.write(f"{n:>13}:  match\n", "text")
    s.r = {"Nsum": float(bad)}
    if bad:
        raise StataError(9, "")


# ---------------------------------------------------------------------------
# datasignature
# ---------------------------------------------------------------------------

def _signature(ds: Dataset) -> str:
    import zlib
    c1 = c2 = 0
    for v in ds.vars:
        raw = "\x00".join(v.raw.tolist()).encode() if v.is_string else \
            np.ascontiguousarray(v.data, dtype=np.float64).tobytes()
        c1 = zlib.crc32(v.name.encode() + raw, c1)
        c2 = zlib.adler32(raw, c2)
    # VERIFY: Stata's checksum algorithm is undocumented; N:k(width) follow the manual
    return f"{ds.nobs}:{ds.nvars}({ds.width()}):{c1}:{c2}"


@command("datasignature", "datasignature")
def cmd_datasignature(s: "Session", args: str) -> None:
    words, opt = _split_opts(args)
    ds = s.data
    out = s.output
    sub = words[0] if words else ""
    sig = _signature(ds)
    chars = ds.chars.setdefault("_dta", {})
    if sub == "":
        match_options(opt, {"fast": 1, "esample": 2, "nodefault": 5, "nonames": 5})
        out.write(f"  {sig}\n", "result")
        s.r = {"datasignature": sig}
        return
    if sub == "set":
        o = match_options(opt, {"reset": 5, "saving": 6})
        if chars.get("datasignature_si") and not o.get("reset"):
            raise StataError(110, "data signature already set -- specify option reset")   # VERIFY
        stamp = _dt.datetime.now().strftime("%d%b%Y %H:%M").lower()
        chars["datasignature_si"] = sig
        chars["datasignature_dt"] = stamp
        out.write(f"  {sig}       (data signature set)\n", "text")
        s.r = {"datasignature": sig}
        return
    if sub == "confirm":
        match_options(opt, {"strict": 6})
        old = chars.get("datasignature_si")
        if not old:
            raise StataError(9, "data signature not set")   # VERIFY
        if old != sig:
            raise StataError(9, "data have changed since signature was set")   # VERIFY
        out.write(f"  (data unchanged since {chars.get('datasignature_dt', '')})\n", "text")
        return
    if sub == "report":
        old = chars.get("datasignature_si")
        if not old:
            out.write("  (data signature not set)\n", "text")
            return
        out.write(f"  data signature set on {chars.get('datasignature_dt', '')}\n", "text")
        out.write(f"  {'set' :>9}: {old}\n  {'current':>9}: {sig}\n", "text")
        out.write("  (data unchanged)\n" if old == sig else "  (data have changed)\n", "text")
        return
    if sub == "clear":
        chars.pop("datasignature_si", None)
        chars.pop("datasignature_dt", None)
        return
    raise StataError(198, f"datasignature {sub}: invalid subcommand")


# ---------------------------------------------------------------------------
# insobs, range
# ---------------------------------------------------------------------------

@command("insobs", "insobs")
def cmd_insobs(s: "Session", args: str) -> None:
    words, opt = _split_opts(args)
    o = match_options(opt, {"before": 1, "after": 1})
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    k = int(float(s.eval(words[0])))
    if k < 1:
        raise StataError(198, "number of observations must be positive")   # VERIFY
    ds = s.data
    n = ds.nobs
    if isinstance(o.get("before"), str):
        at = int(o["before"]) - 1
    elif isinstance(o.get("after"), str):
        at = int(o["after"])
    else:
        at = n
    if not 0 <= at <= n:
        raise StataError(198, "observation number out of range")   # VERIFY
    ds.set_obs(n + k)
    order = np.concatenate([np.arange(at), np.arange(n, n + k), np.arange(at, n)]).astype(int)
    ds.reorder_obs(order)
    ds.sortlist = []
    s.notify_state()


@command("range", "range")
def cmd_range(s: "Session", args: str) -> None:
    words = args.split()
    if len(words) not in (3, 4):
        raise StataError(198, "invalid syntax")
    name = words[0]
    first, last = float(s.eval(words[1])), float(s.eval(words[2]))
    ds = s.data
    n = int(float(s.eval(words[3]))) if len(words) == 4 else ds.nobs
    if n > ds.nobs:
        _note(s, f"number of observations (_N) was {ds.nobs:,}, now {n:,}")
        ds.set_obs(n)
    vals = np.full(ds.nobs, SYS)
    if n == 1:
        vals[0] = first
    elif n > 1:
        vals[:n] = first + np.arange(n) * (last - first) / (n - 1)
    ds.add(_new_numeric(s, name, vals))
    s.notify_state()


# ---------------------------------------------------------------------------
# ipolate
# ---------------------------------------------------------------------------

@command("ipolate", "ipolate", byable=True)
def cmd_ipolate(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist)
    if len(names) != 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"generate": 3, "epolate": 1})
    if not isinstance(o.get("generate"), str):
        raise StataError(198, "option generate() required")
    ds = s.data
    y, x = ds.get(names[0]).data, ds.get(names[1]).data
    mask = touse(s, p)
    groups = s.by_groups.ids if s.by_groups is not None else np.zeros(ds.nobs, dtype=int)
    out = np.full(ds.nobs, SYS)
    for g in np.unique(groups):
        sel = mask & (groups == g)
        known = sel & (x < SYS) & (y < SYS)
        if not known.any():
            continue
        xs, inv = np.unique(x[known], return_inverse=True)
        ys = np.bincount(inv, weights=y[known]) / np.bincount(inv)     # ties: mean of y
        target = sel & (x < SYS)
        xt = x[target]
        val = np.interp(xt, xs, ys)
        outside = (xt < xs[0]) | (xt > xs[-1])
        if o.get("epolate") and len(xs) >= 2:
            lo = (xt < xs[0])
            hi = (xt > xs[-1])
            val[lo] = ys[0] + (xt[lo] - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0])
            val[hi] = ys[-1] + (xt[hi] - xs[-1]) * (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
        else:
            val[outside] = SYS
        out[target] = val
    new = str(o["generate"])
    ds.add(_new_numeric(s, new, out))
    nmiss = int((out >= SYS).sum())
    if nmiss:
        _note(s, _missing_text(nmiss))
    s.notify_state()


# ---------------------------------------------------------------------------
# separate
# ---------------------------------------------------------------------------

@command("separate", "separate")
def cmd_separate(s: "Session", args: str) -> None:
    from .inspect import cell_text
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"by": 2, "generate": 3, "sequential": 3, "missing": 4, "shortlabel": 5})
    if not isinstance(o.get("by"), str):
        raise StataError(198, "option by() required")
    ds = s.data
    v = ds.get(names[0])
    mask = touse(s, p)
    by = str(o["by"]).strip()
    if ds.has(by) or (re.fullmatch(r"\w+", by) and s.expand_varlist(by)):
        g = ds.get(s.expand_varlist(by)[0])
        gvals = np.asarray(g.raw, dtype=object) if g.is_string else g.data
        valid = (np.array([t != "" for t in gvals]) if g.is_string else gvals < SYS)
        bylabel = g.name
        is_var = True
    else:
        e = np.broadcast_to(np.asarray(eval_vector(s, by), dtype=float), (ds.nobs,))
        gvals, valid, bylabel, g, is_var = e, e < SYS, by, None, False
    if not o.get("missing"):
        mask = mask & valid
    levels = sorted(set(gvals[mask].tolist()), key=lambda t: (isinstance(t, str), t))
    stub = str(o["generate"]) if isinstance(o.get("generate"), str) else v.name
    created = []
    for k, lev in enumerate(levels, start=1):
        sel = mask & (gvals == lev)
        if o.get("sequential") or g is None or g.is_string:
            suffix = str(k) if o.get("sequential") or g is not None else str(int(lev)) if lev == int(lev) else str(k)
        else:
            suffix = str(int(lev)) if lev == int(lev) and lev >= 0 else str(k)
        new = stub + suffix
        if v.is_string:
            vals = np.where(sel, np.asarray(v.raw, dtype=object), "")
        else:
            vals = np.where(sel, v.data, SYS)
        i = int(np.flatnonzero(sel)[0])
        shown = cell_text(ds, g, i) if g is not None else format_value(float(lev), "%9.0g", pad=False).strip()
        label = shown if o.get("shortlabel") else f"{v.name}, {bylabel} == {shown}"
        ds.add(Variable(new, v.vtype, vals, fmt=v.fmt, label=label, value_label=v.value_label))
        created.append(new)
    from .inspect import cmd_describe
    cmd_describe(s, " ".join(created))
    s.r = {"varlist": " ".join(created)}
    s.notify_state()


# ---------------------------------------------------------------------------
# stack, xpose, expandcl
# ---------------------------------------------------------------------------

@command("stack", "stack")
def cmd_stack(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"into": 4, "group": 1, "clear": 5, "wide": 4, "keep": 4})
    ds = s.data
    names = s.expand_varlist(p.varlist)
    if isinstance(o.get("into"), str):
        into = str(o["into"]).split()
        k = len(into)
    elif isinstance(o.get("group"), str):
        G = int(o["group"])
        if len(names) % G:
            raise StataError(198, "number of variables not divisible by group()")   # VERIFY
        k = len(names) // G
        into = names[:k]
    else:
        raise StataError(198, "option into() or group() required")
    if len(names) % k:
        raise StataError(198, "number of variables in varlist not a multiple of into()")   # VERIFY
    _need_clear(s, o)
    keep = s.expand_varlist(str(o["keep"])) if isinstance(o.get("keep"), str) else []
    mask = touse(s, p)
    blocks = [names[i:i + k] for i in range(0, len(names), k)]
    n = int(mask.sum())
    new = Dataset()
    new.set_obs(n * len(blocks))
    stack_id = np.repeat(np.arange(1, len(blocks) + 1), n).astype(float)
    new.add(Variable("_stack", "int", stack_id))
    for j, target in enumerate(into):
        cols = [ds.get(b[j]) for b in blocks]
        if any(c.is_string for c in cols) and not all(c.is_string for c in cols):
            raise StataError(109, "type mismatch")
        if cols[0].is_string:
            vals = np.concatenate([np.asarray(c.raw, dtype=object)[mask] for c in cols])
            width = max(c.max_strlen() for c in cols)
            new.add(Variable(target, f"str{min(max(width, 1), 2045)}", vals))
        else:
            vals = np.concatenate([c.data[mask] for c in cols])
            vt = "double" if any(c.vtype == "double" for c in cols) else s.default_type()
            new.add(Variable(target, vt, vals, label=cols[0].label if not o.get("into") else ""))
    extra = []
    if o.get("wide"):
        extra = [nm for nm in names if nm not in into]
    for nm in keep + extra:
        c = ds.get(nm)
        base = np.asarray(c.raw, dtype=object)[mask] if c.is_string else c.data[mask]
        vals = np.concatenate([base] * len(blocks))
        if o.get("wide") and nm in names:
            # wide: each original variable keeps its values only in its own block
            pos = next(b for b, blk in enumerate(blocks) if nm in blk)
            empty = "" if c.is_string else SYS
            vals = np.concatenate([base if b == pos else np.full(n, empty, dtype=object if c.is_string else float)
                                   for b in range(len(blocks))])
        new.add(Variable(nm, c.vtype, vals, fmt=c.fmt, label=c.label, value_label=c.value_label))
    new.value_labels = dict(ds.value_labels)
    new.changed = True
    s.data = new
    s.notify_state()


@command("xpose", "xpose")
def cmd_xpose(s: "Session", args: str) -> None:
    c = find_top(args, ",")
    opt = "" if c == -1 else args[c + 1:]
    if c == -1 or args[:c].strip():
        raise StataError(198, "invalid syntax" if c != -1 else "option clear required")
    o = match_options(opt, {"clear": 5, "format": 3, "varname": 4, "promote": 3})
    if not o.get("clear"):
        raise StataError(198, "option clear required")
    ds = s.data
    nums = [v for v in ds.vars if not v.is_string and v.name != "_varname"]
    X = np.column_stack([v.data for v in nums]) if nums else np.zeros((ds.nobs, 0))
    new = Dataset()
    new.set_obs(len(nums))
    vt = "double" if o.get("promote") or any(v.vtype == "double" for v in nums) else "float"
    for i in range(ds.nobs):
        fmt = str(o["format"]) if isinstance(o.get("format"), str) else ""
        new.add(Variable(f"v{i + 1}", vt, X[i, :] if nums else np.zeros(0), fmt=fmt))
    if o.get("varname"):
        new.add(Variable("_varname", f"str{max([len(v.name) for v in nums] + [1])}",
                         np.array([v.name for v in nums], dtype=object)))
    new.changed = True
    s.data = new
    s.notify_state()


@command("expandcl", "expandcl")
def cmd_expandcl(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"cluster": 2, "generate": 3})
    if not isinstance(o.get("cluster"), str) or not isinstance(o.get("generate"), str):
        raise StataError(198, "options cluster() and generate() required")
    ds = s.data
    exp = (p.exp or p.varlist.lstrip("=")).strip()
    if not exp:
        raise StataError(198, "invalid syntax")
    times = np.broadcast_to(np.asarray(eval_vector(s, exp), dtype=float), (ds.nobs,)).copy()
    mask = touse(s, p)
    cl = s.expand_varlist(str(o["cluster"]))
    keys = [tuple((ds.get(c).raw[i] if ds.get(c).is_string else float(ds.get(c).data[i])) for c in cl)
            for i in range(ds.nobs)]
    # the expansion count must be constant within cluster (VERIFY: message)
    reps = np.where(mask & (times < SYS), np.floor(times), 1).astype(int)
    reps = np.maximum(reps, 0)
    order = []
    newid = []
    nid = 0
    rows_by_cluster: dict[tuple, list[int]] = {}
    for i in range(ds.nobs):
        rows_by_cluster.setdefault(keys[i], []).append(i)
    for key, rows in rows_by_cluster.items():
        k = reps[rows[0]]
        if any(reps[r] != k for r in rows):
            raise StataError(498, "expansion count must be constant within cluster")   # VERIFY
        for rep in range(max(k, 0)):
            nid += 1
            for r in rows:
                order.append(r)
                newid.append(nid)
    created = len(order) - ds.nobs
    for v in ds.vars:
        v.raw = v.raw[np.array(order, dtype=int)] if order else v.raw[:0]
    ds.nobs = len(order)
    ds.esample = None
    ds.sortlist = []
    ds.add(Variable(str(o["generate"]), "long", np.array(newid, dtype=float)))
    if created >= 0:
        _note(s, f"({created:,} observation{'s' if created != 1 else ''} created)")
    else:
        _note(s, f"({-created:,} observation{'s' if created != -1 else ''} deleted)")
    s.notify_state()


# ---------------------------------------------------------------------------
# drawnorm, corr2data
# ---------------------------------------------------------------------------

def _moments(s: "Session", o: dict, k: int):
    """means, standard deviations and correlation (or covariance) from the options."""
    def vec(name):
        v = o.get(name)
        if not isinstance(v, str):
            return None
        return _matrix_arg(s, v)
    means = vec("means")
    means = np.zeros(k) if means is None else np.ravel(means)
    if isinstance(o.get("cov"), str):
        cov = _square(s, str(o["cov"]), k, o)
    else:
        sds = vec("sds")
        sds = np.ones(k) if sds is None else np.ravel(sds)
        corr = _square(s, str(o["corr"]), k, o) if isinstance(o.get("corr"), str) else np.eye(k)
        cov = corr * np.outer(sds, sds)
    if len(means) != k or cov.shape != (k, k):
        raise StataError(503, "conformability error")
    return means, cov


def _matrix_arg(s: "Session", text: str) -> np.ndarray:
    t = text.strip()
    mats = getattr(s, "matrices", {})
    if t in mats:
        return np.asarray(mats[t].data, float)
    words = t.replace(",", " ").split()
    try:
        return np.array([[float(w) for w in words]], float)      # a list of numbers
    except ValueError:
        pass
    if t.startswith("("):
        rows = [r.split(",") if "," in r else r.split() for r in t.strip("()").split("\\")]
        return np.array([[float(s.eval(x)) for x in r if x.strip()] for r in rows], float)
    raise StataError(111, f"matrix {t} not found")


def _square(s: "Session", text: str, k: int, o: dict) -> np.ndarray:
    m = _matrix_arg(s, text)
    store = str(o.get("cstorage", "full")) if isinstance(o.get("cstorage"), str) else "full"
    if store != "full" or m.shape != (k, k):
        vals = np.ravel(m)
        full = np.zeros((k, k))
        idx = np.tril_indices(k) if store in ("lower", "full") else np.triu_indices(k)
        if len(vals) != len(idx[0]):
            raise StataError(503, "conformability error")
        full[idx] = vals
        m = np.where(full != 0, full, full.T) if store == "lower" else full + full.T - np.diag(np.diag(full))
        if store != "upper":
            m = full + full.T - np.diag(np.diag(full))
    return m


def _chol(cov: np.ndarray, force: bool) -> np.ndarray:
    try:
        return np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        if not force:
            raise StataError(506, "matrix not positive definite")
        w, V = np.linalg.eigh(cov)
        return V @ np.diag(np.sqrt(np.maximum(w, 0)))


def _prepare_obs(s: "Session", o: dict) -> int:
    ds = s.data
    if isinstance(o.get("n"), str):
        n = int(float(s.eval(str(o["n"]))))
        if ds.nobs and n != ds.nobs:
            if not o.get("clear"):
                raise StataError(4, "no; data in memory would be lost")   # VERIFY
            s.data = Dataset()
            ds = s.data
        if n > ds.nobs:
            ds.set_obs(n)
        return n
    if o.get("clear"):
        s.data = Dataset()
    if not s.data.nobs:
        raise StataError(2000, "no observations")   # VERIFY: n() required with empty data
    return s.data.nobs


_GEN_OPTS = {"clear": 5, "double": 3, "n": 1, "sds": 3, "corr": 4, "cov": 3, "cstorage": 2,
             "forcepsd": 5, "means": 2, "seed": 4}


@command("drawnorm", "drawnorm")
def cmd_drawnorm(s: "Session", args: str) -> None:
    from ..core import rng
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    names = body.split()
    if not names:
        raise StataError(198, "invalid syntax")
    o = match_options(opt, _GEN_OPTS)
    if isinstance(o.get("seed"), str):
        s.run(f"quietly set seed {o['seed']}")
    k = len(names)
    means, cov = _moments(s, o, k)
    L = _chol(cov, bool(o.get("forcepsd")))
    n = _prepare_obs(s, o)
    # standard normal draws, one variable after the other (VERIFY: Stata's order of draws)
    Z = np.column_stack([rng.draw("rnormal", [0.0, 1.0], n) for _ in range(k)])
    X = Z @ L.T + means
    ds = s.data
    vt = "double" if o.get("double") else s.default_type()
    for j, nm in enumerate(names):
        col = np.full(ds.nobs, SYS)
        col[:n] = X[:, j]
        ds.add(Variable(nm, vt, col))
    s.notify_state()


@command("corr2data", "corr2data")
def cmd_corr2data(s: "Session", args: str) -> None:
    from ..core import rng
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    names = body.split()
    if not names:
        raise StataError(198, "invalid syntax")
    o = match_options(opt, _GEN_OPTS)
    k = len(names)
    means, cov = _moments(s, o, k)
    L = _chol(cov, bool(o.get("forcepsd")))
    n = _prepare_obs(s, o)
    if n <= k:
        raise StataError(2001, "insufficient observations")   # VERIFY
    Z = np.column_stack([rng.draw("rnormal", [0.0, 1.0], n) for _ in range(k)])
    Z = Z - Z.mean(0)
    S = np.cov(Z, rowvar=False).reshape(k, k)
    Z = Z @ np.linalg.inv(np.linalg.cholesky(S)).T      # exactly identity sample covariance
    X = Z @ L.T + means
    ds = s.data
    vt = "double" if o.get("double") else s.default_type()
    for j, nm in enumerate(names):
        col = np.full(ds.nobs, SYS)
        col[:n] = X[:, j]
        ds.add(Variable(nm, vt, col))
    s.notify_state()


# ---------------------------------------------------------------------------
# outfile
# ---------------------------------------------------------------------------

@command("outfile", "outfile")
def cmd_outfile(s: "Session", args: str) -> None:
    from pathlib import Path
    import os
    p = parse_standard(args)
    if p.using is None:
        raise StataError(100, "using required")
    o = match_options(p.options, {"dictionary": 3, "nolabel": 5, "noquote": 5, "replace": 7, "wide": 1,
                                  "comma": 1, "rjs": 3, "fjs": 3, "runtogether": 3, "missing": 2})
    ds = s.data
    names = unique(expand(ds, p.varlist)) if p.varlist.strip() else ds.names
    path = Path(os.path.expanduser(strip_outer_quotes(p.using)))
    if path.suffix == "":
        path = path.with_suffix(".dct" if o.get("dictionary") else ".raw")
    if path.exists() and not o.get("replace"):
        raise StataError(602, f"file {path} already exists")
    mask = touse(s, p)
    from .inspect import cell_text
    sep = "," if o.get("comma") else " "
    quote = not o.get("noquote")
    lines = []
    if o.get("dictionary"):
        lines.append("dictionary {")
        for n in names:
            v = ds.get(n)
            lab = f' "{v.label}"' if v.label else ""
            lines.append(f"\t{v.vtype:<8} {n:<12}{lab}")
        lines.append("}")
    for i in np.flatnonzero(mask):
        fields = []
        for n in names:
            v = ds.get(n)
            if v.is_string:
                t = str(v.raw[i])
                fields.append(f'"{t}"' if quote and not o.get("runtogether") else t)
                continue
            x = float(v.data[i])
            if x >= SYS:
                fields.append("." if not o.get("missing") or x == SYS else M.missing_name(x))
                continue
            if v.value_label and not o.get("nolabel") and x == int(x) and \
                    int(x) in ds.value_labels.get(v.value_label, {}):
                t = cell_text(ds, v, i)
                fields.append(f'"{t}"' if quote else t)
            elif o.get("comma") or o.get("runtogether"):
                fields.append(format_value(x, v.fmt, pad=False).strip().replace(",", ""))
            else:
                fields.append(format_value(x, v.fmt, pad=True))
        line = ("" if o.get("runtogether") else sep).join(fields)
        if not o.get("wide") and not o.get("comma") and len(line) > 80:
            # VERIFY: Stata breaks long records over several lines
            parts, cur = [], ""
            for f in fields:
                if cur and len(cur) + 1 + len(f) > 80:
                    parts.append(cur)
                    cur = f
                else:
                    cur = f if not cur else cur + sep + f
            parts.append(cur)
            line = "\n".join(parts)
        lines.append(line)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# putmata
# ---------------------------------------------------------------------------

@command("putmata", "putmata")
def cmd_putmata(s: "Session", args: str) -> None:
    from ..mata.values import MV
    p = parse_standard(args)
    o = match_options(p.options, {"omitmissing": 4, "replace": 7, "view": 1})
    ds = s.data
    mask = touse(s, p)
    specs = []          # (Mata name, [variables])
    text = p.varlist + (" = " + p.exp if p.exp else "")
    for m in re.finditer(r"(\w+)\s*=\s*\(([^)]*)\)|(\w+)\s*=\s*(\w+)|(\S+)", text):
        if m.group(1):
            specs.append((m.group(1), s.expand_varlist(m.group(2))))
        elif m.group(3):
            specs.append((m.group(3), s.expand_varlist(m.group(4))))
        else:
            for nm in s.expand_varlist(m.group(5)):
                specs.append((nm, [nm]))
    used = [v for _, vs in specs for v in vs]
    if o.get("omitmissing"):
        for nm in used:
            v = ds.get(nm)
            mask &= (np.array([t != "" for t in v.raw]) if v.is_string else v.data < SYS)
    g = s.mata.globals.vars
    for name, _ in specs:
        if name in g and not o.get("replace"):
            raise StataError(110, f"Mata matrix {name} already exists")   # VERIFY
    nvec = nmat = 0
    for name, vs in specs:
        cols = [ds.get(v) for v in vs]
        if any(c.is_string for c in cols):
            if not all(c.is_string for c in cols):
                raise StataError(109, "type mismatch")
            arr = np.column_stack([np.asarray(c.raw, dtype=object)[mask] for c in cols])
            g[name] = MV(arr, "string")
        else:
            arr = np.column_stack([c.data[mask] for c in cols]) if cols else np.zeros((int(mask.sum()), 0))
            g[name] = MV(arr, "real")
        if len(vs) == 1:
            nvec += 1
        else:
            nmat += 1
    parts = []
    if nvec:
        parts.append(f"{nvec} vector{'s' if nvec != 1 else ''}")
    if nmat:
        parts.append(f"{nmat} matri{'ces' if nmat != 1 else 'x'}")
    _note(s, f"({', '.join(parts)} posted)")


# ---------------------------------------------------------------------------
# statsby
# ---------------------------------------------------------------------------

@command("statsby", "statsby")
def cmd_statsby(s: "Session", args: str) -> None:
    colon = find_top(args, ":")
    if colon == -1:
        raise StataError(198, "invalid syntax")
    head, cmd = args[:colon], args[colon + 1:].strip()
    c = find_top(head, ",")
    explist, opt = (head, "") if c == -1 else (head[:c], head[c + 1:])
    o = match_options(opt, {"by": 2, "clear": 5, "saving": 3, "total": 2, "subsets": 3, "nodots": 5,
                            "noisily": 1, "basepop": 4, "nolegend": 5, "verbose": 1, "force": 5})
    if not cmd:
        raise StataError(198, "command required")
    ds0 = s.data
    if not isinstance(o.get("saving"), str):
        _need_clear(s, o)
    byvars = s.expand_varlist(str(o["by"]).split(",")[0]) if isinstance(o.get("by"), str) else []
    keep_missing = isinstance(o.get("by"), str) and "missing" in str(o["by"]).split(",", 1)[-1]
    exps = _parse_explist(explist)
    full = ds0.copy()
    groups: list[tuple[tuple, np.ndarray]] = []
    if byvars:
        cols = [full.get(b) for b in byvars]
        keys: dict[tuple, list[int]] = {}
        for i in range(full.nobs):
            key = tuple((c.raw[i] if c.is_string else float(c.data[i])) for c in cols)
            if not keep_missing and any((k == "" if c.is_string else k >= SYS) for k, c in zip(key, cols)):
                continue
            keys.setdefault(key, []).append(i)
        for key in sorted(keys, key=lambda t: tuple((isinstance(x, str), x) for x in t)):
            m = np.zeros(full.nobs, dtype=bool)
            m[keys[key]] = True
            groups.append((key, m))
    if o.get("total") or not byvars:
        groups.append((None, np.ones(full.nobs, dtype=bool)))
    out = s.output
    results = []
    for gi, (key, m) in enumerate(groups):
        sub = full.copy()
        sub.keep_obs(m)
        s.data = sub
        s.r, s.e = {}, {}
        ok = True
        if not o.get("noisily"):
            out.quiet_depth += 1
        try:
            s.interp._execute_expanded(cmd)
        except StataError:
            ok = False
        finally:
            if not o.get("noisily"):
                out.quiet_depth -= 1
        if gi == 0:
            if not exps:
                exps = _default_explist(s)
            elif any(e in ("_b", "_se") for _, e in exps):
                exps = [x for name, e in exps
                        for x in ([(nm.replace("_b_", e + "_", 1), ex.replace("_b[", e + "[", 1))
                                   for nm, ex in _default_explist(s, coefs=True)] if e in ("_b", "_se")
                                  else [(name, e)])]
        row = []
        for name, e in exps:
            try:
                val = float(s.eval(e)) if ok else SYS
            except (StataError, TypeError, ValueError):
                val = SYS
            row.append(val)
        results.append((key, row))
    s.data = ds0
    # legend (VERIFY: layout)
    if not o.get("nolegend"):
        out.write(f"\n{'command:':>14}  {cmd}\n", "text")
        for name, e in exps:
            out.write(f"{name + ':':>14}  {e}\n", "text")
        if byvars:
            out.write(f"{'by:':>14}  {' '.join(byvars)}\n", "text")
    if not o.get("nodots"):
        out.write("\nStatsby groups\n", "text")
        out.write("----+--- 1 ---+--- 2 ---+--- 3 ---+--- 4 ---+--- 5 \n", "text")
        out.write("." * len(groups) + "\n", "text")
    new = Dataset()
    new.set_obs(len(results))
    for j, b in enumerate(byvars):
        c = full.get(b)
        vals = [(k[j] if k is not None else ("" if c.is_string else SYS)) for k, _ in results]
        new.add(Variable(b, c.vtype, np.array(vals, dtype=object if c.is_string else float), fmt=c.fmt,
                         label=c.label, value_label=c.value_label))
    for j, (name, e) in enumerate(exps):
        new.add(Variable(name, "float", np.array([r[j] for _, r in results], float), label=e))
    new.value_labels = dict(full.value_labels)
    new.changed = True
    if isinstance(o.get("saving"), str):
        from pathlib import Path
        import os
        from ..io.dta import write_dta
        sv = str(o["saving"])
        cc = find_top(sv, ",")
        fname, rep = (sv, "") if cc == -1 else (sv[:cc], sv[cc + 1:])
        path = Path(os.path.expanduser(strip_outer_quotes(fname.strip())))
        if path.suffix == "":
            path = path.with_suffix(".dta")
        if path.exists() and "replace" not in rep:
            raise StataError(602, f"file {path} already exists")
        write_dta(new, path)
        return
    s.data = new
    s.notify_state()


def _parse_explist(text: str) -> list[tuple[str, str]]:
    out = []
    k = 0
    for m in re.finditer(r"(\w+)\s*=\s*(\([^()]*(?:\([^()]*\)[^()]*)*\)|\S+)|(\S+)", text.strip()):
        if m.group(1):
            out.append((m.group(1), m.group(2).strip()))
        else:
            k += 1
            e = m.group(3)
            if e in ("_b", "_se"):
                out.append((e, e))      # expanded later from e(b)
                continue
            out.append((f"_stat_{k}", e))
    return out


def _default_explist(s: "Session", coefs: bool = False) -> list[tuple[str, str]]:
    """No exp_list: _b after an estimation command, otherwise the r() scalars."""
    b = s.e.get("b") if isinstance(s.e, dict) else None
    if b is not None or coefs:
        if b is None:
            return []
        names = list(getattr(b, "colnames", []))
        return [(f"_b_{n.strip('_') or 'cons'}" if n != "_cons" else "_b_cons", f"_b[{n}]") for n in names]
    return [(k, f"r({k})") for k, v in (s.r or {}).items() if isinstance(v, float)]
