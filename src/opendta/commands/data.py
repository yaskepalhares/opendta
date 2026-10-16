"""Phase 1a data commands: set obs, generate, replace, drop, keep,
count, rename, order, sort, gsort, format, compress, input.

Messages follow Stata 14 ([D] manuals); those that depend on undocumented
details are marked VERIFY and have cases in compat/.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import sorting
from ..core.dataset import (NUMERIC_TYPES, STR_MAX, TYPE_ORDER, Variable,
                            check_name, default_format, fit_numeric,
                            is_string_type, smallest_type_for, str_len,
                            str_type_for)
from ..core.errors import StataError
from ..core.formats import is_format, parse_format
from ..core.varlist import expand, resolve_name, unique
from ..lang.expr import evaluate, parse
from ..lang.syntax import match_options, parse_standard
from ..lang.vexpr import broadcast, is_str
from ._util import (ObsContext, check_kind, eval_sample, groups, plural,
                    references_subscript, touse)
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _note(s: "Session", text: str) -> None:
    s.output.write(text + "\n", "text")


# ---------------------------------------------------------------------------
# set obs
# ---------------------------------------------------------------------------

def set_obs(s: "Session", value: str) -> None:
    try:
        n = int(float(s.expand(value).replace(",", "")))
    except ValueError:
        raise StataError(198, "invalid syntax")
    old = s.data.nobs
    if n < old:
        raise StataError(198, f"obs must be at least {old}" if old else "invalid syntax")
    s.data.set_obs(n)
    _note(s, f"number of observations (_N) was {old}, now {n}")


# ---------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------

_TYPE_RE = re.compile(r"^(byte|int|long|float|double|str\d+|strL|str)$")


@command("generate", "g", byable=True)
def cmd_generate(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.exp is None:
        raise StataError(198, "=exp required")
    words = p.varlist.split()
    vtype = None
    if len(words) == 2 and _TYPE_RE.match(words[0]):
        vtype, target = words
    elif len(words) == 1:
        target = words[0]
    else:
        raise StataError(198, "invalid syntax")
    vlabel = ""
    if ":" in target:
        target, vlabel = target.split(":", 1)
    check_name(target)
    ds = s.data
    if ds.has(target):
        raise StataError(110, f"variable {target} already defined")
    opts = match_options(p.options, {"before": 3, "after": 2})

    mask = touse(s, p)
    value = eval_sample(s, p.exp, mask)
    n = ds.nobs
    want_string = vtype is not None and vtype.startswith("str")
    if vtype is not None:
        check_kind(value, want_string)
    string = is_str(value)

    if string:
        col = broadcast(value, n).astype(object).copy()
        col[~mask] = ""
        longest = max((str_len(x) for x in col), default=1) or 1
        if vtype in (None, "str"):
            vtype = str_type_for(longest)
        elif vtype != "strL":
            cap = int(vtype[3:])
            col = np.array([x.encode("utf-8")[:cap].decode("utf-8", "ignore") for x in col], dtype=object)
        nmiss = int(sum(1 for x in col if x == ""))
    else:
        if vtype is None:
            vtype = s.default_type()
        col = broadcast(value, n).astype(np.float64).copy()
        col[~mask] = M.SYSMISS
        # counts only missing values from the expression and obs excluded by if/in; values out of
        # the type's range also become missing but are not counted (gen byte b = 200 in 1/3 reports
        # 3, observed in Stata)
        nmiss = int(np.sum(col >= M.SYSMISS))
        col, _ = fit_numeric(col, vtype)

    var = Variable(target, vtype, col, value_label=vlabel)
    position = None
    if opts.get("before"):
        position = ds.index(resolve_name(ds, str(opts["before"])))
    elif opts.get("after"):
        position = ds.index(resolve_name(ds, str(opts["after"]))) + 1
    ds.add(var, position=position)
    if nmiss:
        _note(s, f"({plural(nmiss, 'missing value')} generated)")
    s.notify_state()


# ---------------------------------------------------------------------------
# replace
# ---------------------------------------------------------------------------

@command("replace", "replace", byable=True)
def cmd_replace(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.exp is None:
        raise StataError(198, "=exp required")
    ds = s.data
    names = p.varlist.split()
    if len(names) != 1:
        raise StataError(198, "invalid syntax")
    # VERIFY: does Stata require the full name in replace?
    var = ds.get(resolve_name(ds, names[0]))
    opts = match_options(p.options, {"nopromote": 6})
    mask = touse(s, p)
    rows = np.flatnonzero(mask)
    node = parse(p.exp)

    if references_subscript(node, var.name):
        new_values = _replace_sequential(s, var, node, rows)
    else:
        value = eval_sample(s, p.exp, mask)
        check_kind(value, var.is_string)
        new_values = broadcast(value, ds.nobs)[rows]

    old = var.data.copy()
    if var.is_string:
        changed, note = ds.set_string(var, np.asarray(new_values, dtype=object), rows)
        to_missing = int(sum(1 for i in rows if var.raw[i] == "" and old[i] != ""))
    else:
        changed, note = ds.set_numeric(var, np.asarray(new_values, dtype=np.float64), rows,
                                       promote=not opts.get("nopromote"))
        to_missing = int(np.sum((var.data[rows] >= M.SYSMISS) & (old[rows] < M.SYSMISS)))
    if note:
        _note(s, note)
    msg = f"({plural(changed, 'real change')} made"
    if to_missing:
        msg += f", {to_missing:,} to missing"
    _note(s, msg + ")")
    s.notify_state()


def _replace_sequential(s: "Session", var: Variable, node, rows: np.ndarray) -> np.ndarray:
    """replace x = x[_n-1] ... uses values already replaced in earlier obs."""
    ctx = ObsContext(s, groups(s))
    original = var.raw.copy()
    out = []
    try:
        for i in rows:
            ctx.i = int(i)
            v = evaluate(node, ctx)
            check_kind(v, var.is_string)
            if not var.is_string:
                v, _ = fit_numeric(np.array([v]), var.vtype)
                v = float(v[0])
            var.set_value(i, v)
            out.append(v)
    finally:
        var.raw = original
    return np.array(out, dtype=object if var.is_string else np.float64)


# ---------------------------------------------------------------------------
# drop / keep
# ---------------------------------------------------------------------------

def _drop_keep(s: "Session", args: str, keep: bool) -> None:
    p = parse_standard(args)
    ds = s.data
    if p.varlist and (p.if_ or p.in_):
        raise StataError(198, "invalid syntax")
    if p.varlist:
        names = unique(expand(ds, p.varlist))
        if keep:
            names = [n for n in ds.names if n not in set(names)]
        ds.drop_vars(names)
        s.notify_state()
        return
    if not (p.if_ or p.in_):
        raise StataError(100, "varlist required")
    mask = touse(s, p)
    keep_mask = mask if keep else ~mask
    removed = ds.keep_obs(keep_mask)
    _note(s, f"({plural(removed, 'observation')} deleted)")
    s.notify_state()


@command("drop", byable=True)
def cmd_drop(s: "Session", args: str) -> None:
    _drop_keep(s, args, keep=False)


@command("keep", byable=True)
def cmd_keep(s: "Session", args: str) -> None:
    _drop_keep(s, args, keep=True)


# ---------------------------------------------------------------------------
# sample
# ---------------------------------------------------------------------------

@command("sample", byable=True)
def cmd_sample(s: "Session", args: str) -> None:
    """sample # [if] [in] [, count by(groupvars)] ([D] sample): simple random sampling without
    replacement. # is the percentage of observations to keep or, with count, their number.
    Observations outside if/in are kept; with by() or the by prefix, each group is sampled
    separately.

    As observed in Stata 14 (compat 0310): two uniform draws per observation, in the current
    order (the first set orders the observations, the second breaks ties; VERIFY), also with
    by (VERIFY); the data end up sorted by group, then the observations in if/in, then the
    draws, with no sort order recorded; in each group the first n candidates are kept, n being
    # or round(# x size/100). sample 100 draws nothing."""
    from ..core import rng
    from ..core.grouping import group_ids
    p = parse_standard(args)
    if p.weight is not None or p.using is not None:
        raise StataError(101, "weights not allowed" if p.weight is not None else "using not allowed")
    text = p.varlist.strip()
    if not text:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"count": 5, "by": 2}) if p.options.strip() else {}
    count = bool(o.get("count"))
    try:
        num = float(text)
    except ValueError:
        try:
            num = float(s.eval(text))
        except StataError:
            raise StataError(198, "invalid syntax") from None
    if count and num != int(num):
        raise StataError(126, f"'{text}' found where integer expected")   # VERIFY: rc
    if count and num < 0:
        # Stata 14 explains and carries on, keeping no observation (compat 0310)
        s.output.write(f"{text}:  # must be 0 or positive; with the count option, it represents "
                       "the # of observations to be drawn\n", "text")
        num = 0.0
    elif not count and not 0 <= num <= 100:
        s.output.write(f"{text}:  # must be between 0 and 100 and is interpreted as a percent unless "
                       "option count is also specified, in which case # is interpreted as the number "
                       "of observations to be drawn\n", "text")
        raise StataError(198, "invalid syntax")
    ds = s.data
    n = ds.nobs
    mask = touse(s, p)
    if isinstance(o.get("by"), str):
        if s.by_groups is not None:
            raise StataError(190, "option by() may not be combined with by prefix")   # VERIFY
        names = unique(expand(ds, str(o["by"])))
        ids = group_ids(ds, names)[0]
    elif s.by_groups is not None:
        ids = np.asarray(s.by_groups.ids, dtype=np.int64)
    else:
        ids = np.zeros(n, dtype=np.int64)
    if not count and num == 100:
        # nothing to drop: Stata does not draw (compat 0310: the stream advances by 2N for
        # sample 0 but not for sample 100; VERIFY count >= N)
        return
    u1 = rng.RNG.uniform(n) if n else np.zeros(0)
    u2 = rng.RNG.uniform(n) if n else np.zeros(0)
    cand = mask & (ids >= 0)
    # group, then the candidates first, then the draws
    order = np.lexsort((u2, u1, ~cand, ids))
    keep = np.ones(n, dtype=bool)
    oc = cand[order]
    og = ids[order]
    sel = order[oc]
    gs = og[oc]
    if len(sel):
        starts = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
        sizes = np.diff(np.r_[starts, len(sel)])
        for a_, size in zip(starts, sizes):
            want = min(int(num), size) if count else int(np.floor(num * size / 100 + 0.5))
            keep[sel[a_ + want:a_ + size]] = False
    ds.reorder_obs(order)
    ds.sortlist = []
    removed = ds.keep_obs(keep[order])
    if removed:
        _note(s, f"({plural(removed, 'observation')} deleted)")
    s.notify_state()


# ---------------------------------------------------------------------------
# count
# ---------------------------------------------------------------------------

@command("count", "count", byable=True)
def cmd_count(s: "Session", args: str) -> None:
    p = parse_standard(args)
    if p.varlist:
        raise StataError(101, "varlist not allowed")
    mask = touse(s, p)
    if s.by_groups is not None and s.data.nobs:
        g = s.by_groups
        ids = g.ids
        for gid in range(int(ids.max()) + 1 if len(ids) else 0):
            sel = ids == gid
            _by_header(s, int(np.flatnonzero(sel)[0]))
            n = int(mask[sel].sum())
            s.output.write(f"  {n:,}\n", "result")
        s.r = {"N": float(n)} if len(ids) else {"N": 0.0}
        return
    n = int(mask.sum())
    # VERIFY: indentation and thousands separator in count's output
    s.output.write(f"  {n:,}\n", "result")
    s.r = {"N": float(n)}


def _by_header(s: "Session", first_obs: int) -> None:
    """'-> g = 1' header before each group (commands that repeat output)."""
    keys = getattr(s, "_by_keys", [])
    ds = s.data
    parts = []
    for k in keys:
        var = ds.get(k)
        v = var.value(first_obs)
        from ..core.formats import format_value
        shown = v if var.is_string else format_value(float(v), var.fmt, pad=False).strip()
        if not var.is_string and var.value_label:
            lab = ds.value_labels.get(var.value_label, {})
            if float(v) == int(float(v)) and int(float(v)) in lab:
                shown = lab[int(float(v))]
        parts.append(f"{k} = {shown}")
    # the line takes c(linesize) (255 dashes in compat; VERIFY with the default 80)
    width = int(s.settings.get("linesize", 80))
    s.output.write("\n" + "-" * width + "\n", "text")
    s.output.write("-> " + ", ".join(parts) + "\n", "text")


# ---------------------------------------------------------------------------
# rename / order
# ---------------------------------------------------------------------------

@command("rename", "ren")
def cmd_rename(s: "Session", args: str) -> None:
    words = args.split()
    if len(words) != 2:
        raise StataError(198, "invalid syntax")
    ds = s.data
    old = resolve_name(ds, words[0])
    ds.rename(old, words[1])
    s.notify_state()


@command("order")
def cmd_order(s: "Session", args: str) -> None:
    p = parse_standard(args)
    ds = s.data
    names = unique(expand(ds, p.varlist, allow_empty=False))
    opts = match_options(p.options, {"first": 5, "last": 4, "before": 3, "after": 2,
                                     "alphabetic": 5, "sequential": 3})
    if opts.get("alphabetic"):
        names = sorted(names, key=str.lower)
    before = resolve_name(ds, str(opts["before"])) if opts.get("before") else None
    after = resolve_name(ds, str(opts["after"])) if opts.get("after") else None
    ds.order(names, last=bool(opts.get("last")), before=before, after=after)
    s.notify_state()


# ---------------------------------------------------------------------------
# sort / gsort
# ---------------------------------------------------------------------------

@command("sort", "so")
def cmd_sort(s: "Session", args: str) -> None:
    p = parse_standard(args)
    ds = s.data
    names = unique(expand(ds, p.varlist, allow_empty=False))
    match_options(p.options, {"stable": 6})
    if p.in_:
        # [D] sort: only the observations in the range are sorted; the dataset is then not
        # marked as sorted (VERIFY: c(sortlist) after sort ... in)
        from ..lang.syntax import parse_in
        lo, hi = parse_in(p.in_, ds.nobs)
        order = sorting.sort_order(ds, names)
        sub = np.arange(lo, hi)
        rank = np.empty(ds.nobs, dtype=np.int64)
        rank[order] = np.arange(ds.nobs)
        full = np.arange(ds.nobs)
        full[sub] = sub[np.argsort(rank[sub], kind="stable")]
        if not np.all(full == np.arange(ds.nobs)):
            ds.reorder_obs(full)
        ds.sortlist = []
        s.notify_state()
        return
    sorting.sort(ds, names)
    s.notify_state()


@command("gsort")
def cmd_gsort(s: "Session", args: str) -> None:
    p = parse_standard(args)
    ds = s.data
    opts = match_options(p.options, {"mfirst": 2, "generate": 3})
    keys: list[str] = []
    desc: list[bool] = []
    for tok in p.varlist.replace("+ ", "+").replace("- ", "-").split():
        d = tok.startswith("-")
        name = tok.lstrip("+-")
        for n in expand(ds, name):
            keys.append(n)
            desc.append(d)
    if not keys:
        raise StataError(100, "varlist required")
    sorting.sort(ds, keys, desc, missing_first=bool(opts.get("mfirst")))
    if opts.get("generate"):
        newname = str(opts["generate"])
        ds.add(Variable(newname, "long", np.arange(1, ds.nobs + 1, dtype=np.float64)))
    s.notify_state()


# ---------------------------------------------------------------------------
# format
# ---------------------------------------------------------------------------

@command("format")
def cmd_format(s: "Session", args: str) -> None:
    ds = s.data
    words = args.split()
    fmt_pos = [i for i, w in enumerate(words) if w.startswith("%")]
    if not fmt_pos:
        # `format varlist` lists the formats
        for n in expand(ds, args):
            s.output.write(f"{n:<16}{ds.get(n).fmt}\n", "text")
        return
    i = fmt_pos[0]
    fmt_text = words[i]
    names = words[:i] if i > 0 else words[i + 1:]
    if not is_format(fmt_text):
        raise StataError(120, f"{fmt_text} invalid %format")
    fmt = parse_format(fmt_text)
    for n in unique(expand(ds, " ".join(names), allow_empty=False)):
        var = ds.get(n)
        if var.is_string != (fmt.kind == "s"):
            raise StataError(120, f"{fmt_text} invalid %format")
        var.fmt = fmt_text
    s.notify_state()


# ---------------------------------------------------------------------------
# compress
# ---------------------------------------------------------------------------

_BYTES = {"byte": 1, "int": 2, "long": 4, "float": 4, "double": 8}


@command("compress")
def cmd_compress(s: "Session", args: str) -> None:
    ds = s.data
    names = expand(ds, args) if args.strip() else ds.names
    # numeric first, then strings (order observed in Stata)
    names = [n for n in names if not ds.get(n).is_string] + [n for n in names if ds.get(n).is_string]
    saved = 0
    for n in names:
        var = ds.get(n)
        if var.is_string:
            if var.vtype == "strL":
                continue
            cur = int(var.vtype[3:])
            need = max(1, var.max_strlen())
            if need < cur:
                newt = f"str{need}"
                _note(s, f"  variable {n} was {var.vtype} now {newt}")
                saved += (cur - need) * ds.nobs
                if var.fmt == default_format(var.vtype):
                    var.fmt = default_format(newt)
                var.vtype = newt
            continue
        need = smallest_type_for(var.data)
        if var.vtype == "float" and need == "double":
            need = "float"
        if TYPE_ORDER[need] < TYPE_ORDER[var.vtype]:
            # a float only goes down to an integer type if all values are integers
            _note(s, f"  variable {n} was {var.vtype} now {need}")
            saved += (_BYTES[var.vtype] - _BYTES[need]) * ds.nobs
            if var.fmt == default_format(var.vtype):
                var.fmt = default_format(need)
            var.vtype = need
    # VERIFY: exact format of the summary
    _note(s, f"  ({saved:,} bytes saved)")
    if saved:
        ds.changed = True
    s.notify_state()


# ---------------------------------------------------------------------------
# recast
# ---------------------------------------------------------------------------

def _truncate_bytes(text: str, n: int) -> str:
    b = text.encode("utf-8")[:n]
    return b.decode("utf-8", errors="ignore")


@command("recast")
def cmd_recast(s: "Session", args: str) -> None:
    head, comma, opts_text = args.partition(",")
    opts = match_options(opts_text, {"force": 5}) if comma else {}
    words = head.split(None, 1)
    if len(words) < 2:
        raise StataError(100, "varlist required")
    newt = words[0]
    if newt == "str":
        raise StataError(198, "invalid syntax")
    if not (newt in NUMERIC_TYPES or re.match(r"^str(\d+|L)$", newt)):
        raise StataError(198, f"{newt} invalid type")   # VERIFY
    if newt.startswith("str") and newt != "strL" and not 1 <= int(newt[3:]) <= STR_MAX:
        raise StataError(198, f"{newt} invalid type")
    ds = s.data
    for n in unique(expand(ds, words[1])):
        var = ds.get(n)
        if var.vtype == newt:
            continue
        if var.is_string != is_string_type(newt):
            raise StataError(109, f"{n}:  {newt} invalid")   # VERIFY return code
        if var.is_string:
            if newt == "strL":
                new = var.data
                changed = 0
            else:
                width = int(newt[3:])
                new = np.array([_truncate_bytes(x, width) for x in var.data], dtype=object)
                changed = int(sum(a != b for a, b in zip(var.data, new)))
        else:
            new, _ = fit_numeric(var.data, newt)
            same = (new == var.data) | ((new >= M.SYSMISS) & (var.data >= M.SYSMISS)
                                         & (new.view(np.int64) == var.data.view(np.int64)))
            changed = int((~same).sum())
        if changed and not opts.get("force"):
            # VERIFY: text and return code when there is loss without force
            _note(s, f"{n}:  {plural(changed, 'value')} would be changed; not changed")
            continue
        if changed:
            _note(s, f"{n}:  {plural(changed, 'value')} changed")
        # recast keeps the format and the sort order (observed in Stata 14: double %10.0g becomes
        # int %10.0g; "Sorted by" remains)
        var.vtype = newt
        sortlist = list(ds.sortlist)
        if var.is_string:
            var.data = new
        else:
            ds.set_numeric(var, new, promote=False)
        ds.sortlist = sortlist
        ds.changed = True
    s.notify_state()


# ---------------------------------------------------------------------------
# input (the data come in the following lines, up to `end`)
# ---------------------------------------------------------------------------

def run_input(s: "Session", spec: str, rows: list[str]) -> None:
    from ..lang.words import split_words

    ds = s.data
    tokens = spec.split()
    decl: list[tuple[str, str]] = []
    pending_type = None
    for t in tokens:
        if _TYPE_RE.match(t) and pending_type is None and not ds.has(t):
            pending_type = t
            continue
        decl.append((t, pending_type or ""))
        pending_type = None
    if not decl:
        # no varlist: enters data into the existing variables
        decl = [(n, "") for n in ds.names]
    created = []
    for name, t in decl:
        if ds.has(name):
            continue
        vtype = t or s.default_type()
        if vtype == "str":
            vtype = "str1"
        if vtype.startswith("str"):
            col = np.array([""] * ds.nobs, dtype=object)
        else:
            col = np.full(ds.nobs, M.SYSMISS)
        ds.add(Variable(name, vtype, col))
        created.append(name)
    vars_ = [ds.get(n) for n, _ in decl]
    start = ds.nobs
    ds.set_obs(start + len(rows))
    for k, line in enumerate(rows):
        values = split_words(line)
        for j, var in enumerate(vars_):
            raw = values[j] if j < len(values) else ""
            i = start + k
            if var.is_string:
                text = raw
                if var.vtype != "strL" and str_len(text) > int(var.vtype[3:]):
                    text = text.encode("utf-8")[:int(var.vtype[3:])].decode("utf-8", "ignore")
                var.set_value(i, text)
            else:
                if raw in ("", "."):
                    v = M.SYSMISS
                else:
                    try:
                        v = M.missing_code(raw) if raw.startswith(".") else float(raw)
                    except ValueError:
                        raise StataError(198, f"'{raw}' cannot be read as a number")
                fitted, _ = fit_numeric(np.array([v]), var.vtype)
                var.set_value(i, fitted[0])
    s.notify_state()
