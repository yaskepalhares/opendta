"""Time-series data tools: tsappend, tsfill, tsreport and the rolling prefix ([TS] tsappend,
[TS] tsfill, [TS] tsreport, [TS] rolling).

The layouts of tsreport and rolling follow the manual's examples; the alignment of the lines the
manual shows only as text is marked VERIFY (compat/do/1904_timeseries.do)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..estimation.tsops import _fmt_time, require_ts
from ..lang.syntax import find_top, match_options
from .registry import command
from .tshelp import append_periods, kv, panels, to_nan

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# tsappend, tsfill
# ---------------------------------------------------------------------------

_MONTHS = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]


def date_literal(text: str, fmt: str) -> float:
    """The integer of a date literal of the given %t function (td, tw, tm, tq, th, ty)."""
    import datetime as dt
    t = text.strip().lower()
    f = fmt.strip().lower().lstrip("%")
    try:
        if f == "td":
            m = re.fullmatch(r"(\d{1,2})([a-z]{3})(\d{4})", t)
            d = dt.date(int(m.group(3)), _MONTHS.index(m.group(2)) + 1, int(m.group(1)))
            return float((d - dt.date(1960, 1, 1)).days)
        if f in ("tw", "tm", "tq", "th"):
            letter = f[1]
            y, n = re.fullmatch(r"(\d{4})" + letter + r"(\d+)", t).groups()
            per = {"w": 52, "m": 12, "q": 4, "h": 2}[letter]
            return float((int(y) - 1960) * per + int(n) - 1)
        if f == "ty":
            return float(int(t) - 1960)
    except (AttributeError, ValueError):
        raise StataError(198, f"{text} invalid date for {fmt}")
    raise StataError(198, f"tsfmt({fmt}) is not available in OpenDTA")


@command("tsappend", "tsa")
def cmd_tsappend(s: "Session", args: str) -> None:
    info = require_ts(s.data)
    c = find_top(args, ",")
    opts = args[c + 1:] if c != -1 else ""
    o = match_options(opts, {"add": 3, "last": 4, "tsfmt": 3, "panel": 5}) if opts.strip() else {}
    if isinstance(o.get("add"), str) == isinstance(o.get("last"), str):
        raise StataError(198, "option add() or last() required")
    delta = info.delta
    only = None
    if isinstance(o.get("panel"), str):
        only = float(s.eval(o["panel"]))
    entries: list[tuple[float | None, float]] = []
    for pv, rows, t in panels(s):
        if only is not None and pv != only:
            continue
        if not len(t):
            continue
        tl = float(t.max())
        if isinstance(o.get("add"), str):
            k = int(float(s.eval(o["add"])))
            if k < 1:
                raise StataError(198, "add() must be positive")
            entries += [(pv, tl + delta * i) for i in range(1, k + 1)]
        else:
            if not isinstance(o.get("tsfmt"), str):
                raise StataError(198, "option tsfmt() required with last()")
            last = date_literal(str(o["last"]), o["tsfmt"])
            entries += [(pv, x) for x in np.arange(tl + delta, last + delta / 2, delta)]
    append_periods(s, entries)
    kv(s, [("add", float(len(entries)))])
    s.notify_state()


@command("tsfill", "tsf")
def cmd_tsfill(s: "Session", args: str) -> None:
    info = require_ts(s.data)
    o = match_options(args.lstrip(", "), {"full": 1}) if args.strip() else {}
    delta = info.delta
    groups = panels(s)
    if o.get("full") and info.panel:
        gmin = min(float(t.min()) for _p, _r, t in groups if len(t))
        gmax = max(float(t.max()) for _p, _r, t in groups if len(t))
    entries = []
    for pv, rows, t in groups:
        if not len(t):
            continue
        lo, hi = (gmin, gmax) if (o.get("full") and info.panel) else (float(t.min()), float(t.max()))
        have = set(np.round(t / delta).astype(np.int64).tolist())
        for x in np.arange(lo, hi + delta / 2, delta):
            if int(round(x / delta)) not in have:
                entries.append((pv, float(x)))
    append_periods(s, entries)
    s.notify_state()


# ---------------------------------------------------------------------------
# tsreport
# ---------------------------------------------------------------------------

def _gaps(ds, info, rows_by_panel, include_panel_changes: bool):
    """[(row_before, row_after, panel, start, end, nobs)] for consecutive rows that are more
    than one period apart; a panel change counts as a gap (with the period after the end of the
    panel and no end)."""
    t = ds.get(info.tvar).data.astype(float)
    out = []
    prev = None
    for pv, rows in rows_by_panel:
        if prev is not None and include_panel_changes and len(rows):
            out.append((prev[1], int(rows[0]), prev[0], prev[2] + info.delta, None, None))
        for a, b in zip(rows[:-1], rows[1:]):
            if t[b] - t[a] > info.delta + 1e-9:
                n = int(round((t[b] - t[a]) / info.delta)) - 1
                out.append((int(a), int(b), pv, t[a] + info.delta, t[b] - info.delta, n))
        if len(rows):
            prev = (pv, int(rows[-1]), float(t[rows[-1]]))
    return out


def _block(s: "Session", lines: list[tuple[str, str]], width: int = 15) -> None:
    for lab, val in lines:
        s.output.write(f"{lab:>{width}}{val}\n", "text")


def _report_block(s, info, head, lo, hi, n, ng, fmt):
    """Stata's layout of one tsreport block (no panel variable)."""
    w = max(7, len(fmt(lo)), len(fmt(hi)), len(str(n)), len(str(ng)))
    out = s.output
    out.write("\n", "text")
    if head:
        out.write(f"{head[0]:<18}{head[1]}\n", "text")
    out.write(f"{'Time variable:':<18}{info.tvar}\n", "text")
    out.write("-" * (18 + w) + "\n", "text")
    for lab, val in (("Starting period", fmt(lo)), ("Ending period", fmt(hi)),
                     ("Observations", str(n)), ("Number of gaps", str(ng))):
        out.write(f"{lab:<15} = {val:>{w}}\n", "text")


def _report_gaps(s, gl, fmt):
    out = s.output
    out.write("\nGap report\n" + "-" * 50 + "\n", "text")
    out.write(f"{'Obs.':<13}|{'Start':>11}{'End':>12}{'N. Obs.':>13}\n", "text")
    out.write("-" * 13 + "+" + "-" * 36 + "\n", "text")
    for a, b, _pv, st_, en, n in gl:
        out.write(f"{f'{a + 1}--{b + 1}':>13}|{fmt(st_):>11}{(fmt(en) if en is not None else '.'):>12}"
                  f"{(str(n) if n is not None else '.'):>13}\n", "result")
    out.write("-" * 50 + "\n", "text")


@command("tsreport", "tsrep")
def cmd_tsreport(s: "Session", args: str) -> None:
    from ._util import touse
    from ..lang.syntax import parse_standard
    ds = s.data
    info = require_ts(ds)
    p = parse_standard(args)
    o = match_options(p.options, {"detail": 1, "casewise": 1, "panel": 5}) if p.options.strip() else {}
    names = []
    for tok in p.varlist.split():
        names += s.expand_varlist(tok) if ("*" in tok or "?" in tok or ds.has(tok)) else [tok]
    sample = touse(s, p)
    detail, casewise, nopanel = bool(o.get("detail")), bool(o.get("casewise")), bool(o.get("panel"))
    fmt = lambda x: _fmt_time(ds, info.tvar, x)
    t = ds.get(info.tvar).data.astype(float)
    out = s.output

    def groups_for(mask):
        return [(pv, rows) for pv, rows, _t in panels(s, mask)]

    def head(mask, label_lines):
        rows_all = np.flatnonzero(mask & (t < SYS))
        tt = t[rows_all]
        gl = _gaps(ds, info, groups_for(mask), info.panel != "" and not nopanel)
        lines = list(label_lines) + [
            ("Time variable: ", info.tvar) if not info.panel else (" Time variable: ", info.tvar),
            ("Starting period = ", fmt(tt.min())), ("  Ending period = ", fmt(tt.max())),
            ("   Observations = ", f"{len(rows_all)}"), (" Number of gaps = ", f"{len(gl)}")]
        return gl, lines, rows_all, tt

    def gap_table(gl):
        out.write("Gap report\n", "text")
        cols = ["Obs.", info.panel, "Start", "End", "N. Obs."] if info.panel else \
            ["Obs.", "Start", "End", "N. Obs."]
        head_ = f"{'Obs.':>12}" + (f"{info.panel[:12]:>13}" if info.panel else "") + \
            f"{'Start':>22}{'End':>22}{'N. Obs.':>9}"
        del cols
        out.write(head_ + "\n", "text")
        for a, b, pv, st_, en, n in gl:
            line = f"{a + 1:>6}{b + 1:>6}"
            if info.panel:
                line += f"{pv:>13g}"
            line += f"{fmt(st_):>22}{(fmt(en) if en is not None else '.'):>22}{(str(n) if n is not None else '.'):>9}"
            out.write(line + "\n", "result")

    if not names or casewise:
        mask = sample.copy()
        for nm in names:
            mask &= ~np.isnan(to_nan(ds.get(nm).data))
        if names:
            lab = "Variables: " + (" and ".join(names) if len(names) <= 2 else ", ".join(names))
        pre = []
        if info.panel:
            pre.append(("Panel variable: ", info.panel))
        if names:
            out.write(lab + "\n", "text")
        gl, lines, rows_all, tt = head(mask, pre)
        if info.panel:
            _block(s, lines)
        else:
            _report_block(s, info, None, tt.min(), tt.max(), len(rows_all), len(gl), fmt)
        if info.panel and not nopanel:
            out.write("(Gap count includes panel changes)\n", "text")
        if detail and gl:
            gap_table(gl) if info.panel else _report_gaps(s, gl, fmt)
        from ..commands.matrix import Matrix
        pairs = [("N_gaps", float(len(gl))), ("N_obs", float(len(rows_all))),
                 ("start", float(tt.min())), ("end", float(tt.max())),
                 ("tsfmt", ds.get(info.tvar).fmt)]
        if detail and gl:
            tab = np.array([[g[3], g[4] if g[4] is not None else np.nan] for g in gl])
            pairs.append(("table", Matrix(np.where(np.isnan(tab), SYS, tab), [str(i + 1) for i in range(len(gl))],
                                          ["start", "end"])))
        kv(s, pairs)
        return
    # one report per variable
    pairs = []
    summary = []
    for k, nm in enumerate(names, start=1):
        v = to_nan(ds.get(nm).data)
        mask = sample & ~np.isnan(v)
        gl = _gaps(ds, info, groups_for(mask), info.panel != "" and not nopanel)
        rows_all = np.flatnonzero(mask & (t < SYS))
        tt = t[rows_all]
        if not len(tt):
            raise StataError(2000, "no observations")
        summary.append((nm, tt.min(), tt.max(), len(rows_all), len(gl), gl))
        pairs += [(f"N_gaps{k}", float(len(gl))), (f"N_obs{k}", float(len(rows_all))),
                  (f"start{k}", float(tt.min())), (f"end{k}", float(tt.max())), (f"var{k}", nm)]
    if detail:
        for nm, lo, hi, n, ng, gl in summary:
            if not info.panel:
                _report_block(s, info, ("Variable:", nm), lo, hi, n, ng, fmt)
                if gl:
                    _report_gaps(s, gl, fmt)
                continue
            out.write(f"Variable: {nm}\n", "text")
            lines = ([("Panel variable: ", info.panel)] if info.panel else []) + \
                [(" Time variable: ", info.tvar), ("Starting period = ", fmt(lo)),
                 ("  Ending period = ", fmt(hi)), ("   Observations = ", f"{n}"),
                 (" Number of gaps = ", f"{ng}")]
            _block(s, lines)
            if gl:
                gap_table(gl)
    else:
        out.write("Gap summary report\n", "text")
        w = max(8, max(len(x[0]) for x in summary) + 1)
        out.write(f"{'':<{w}}{'':>22}{'':>22}{'Number of':>11}\n", "text")
        out.write(f"{'Variable':<{w}}{'Start':>22}{'End':>22}{'Obs.':>7}{'Gaps':>6}\n", "text")
        for nm, lo, hi, n, ng, gl in summary:
            out.write(f"{nm:<{w}}{fmt(lo):>22}{fmt(hi):>22}{n:>7}{ng:>6}\n", "result")
    pairs.append(("tsfmt", ds.get(info.tvar).fmt))
    kv(s, pairs)


# ---------------------------------------------------------------------------
# rolling
# ---------------------------------------------------------------------------

_ROLL = {"window": 1, "recursive": 3, "rrecursive": 2, "stepsize": 4, "start": 5, "end": 3,
         "keep": 4, "clear": 5, "saving": 6, "nodots": 5, "dots": 4, "noisily": 3, "trace": 5,
         "reject": 3, "force": 5, "double": 6}


def cmd_rolling(s: "Session", args: str) -> None:
    from ..estimation.bootstrap import _evaluate, _subset, parse_exp_list
    from ..estimation.resample import _dot, _dots, _named, _quiet, _run0, _split
    from ..lang.syntax import parse_standard
    from ._util import touse
    ds0 = s.data
    info = require_ts(ds0)
    head, opt_text, command_ = _split(args)
    ifin = ""
    for kw in (" if ", " in "):
        c = find_top(head, kw)
        if c != -1:
            ifin = head[c:]
            head = head[:c]
            break
    o = match_options(opt_text, _ROLL) if opt_text.strip() else {}
    if not isinstance(o.get("window"), str):
        raise StataError(198, "option window() required")
    window = int(float(s.eval(o["window"])))
    if window < 1:
        raise StataError(198, "window() must be positive")
    step = int(float(s.eval(o["stepsize"]))) if isinstance(o.get("stepsize"), str) else 1
    if o.get("recursive") and o.get("rrecursive"):
        raise StataError(198, "options recursive and rrecursive may not be combined")
    items = _named(parse_exp_list(head), "_stat_")
    out = s.output
    noisy = bool(o.get("noisily")) or bool(o.get("trace"))
    cmdname, eclass, e0, est0 = _run0(s, command_, noisy)
    if not items:
        if eclass:
            items = [("", "_b")]
        else:
            items = [(k, f"r({k})") for k, v in s.r.items() if isinstance(v, (int, float))]
            if not items:
                raise StataError(198, "expression list required")
    names, _obs = _evaluate(s, items)
    ds = s.data
    delta = info.delta
    p = parse_standard("x" + ifin)
    base_mask = touse(s, p)
    t = ds.get(info.tvar).data.astype(float)
    tt = t[base_mask & (t < SYS)]
    if not len(tt):
        raise StataError(2000, "no observations")
    first = float(s.eval(o["start"])) if isinstance(o.get("start"), str) else float(tt.min())
    last = float(s.eval(o["end"])) if isinstance(o.get("end"), str) else float(tt.max())
    wins = []
    if o.get("recursive"):
        e = first + (window - 1) * delta
        while e <= last + 1e-9:
            wins.append((first, e))
            e += step * delta
    elif o.get("rrecursive"):
        b = last - (window - 1) * delta
        while b >= first - 1e-9:
            wins.append((b, last))
            b -= step * delta
        wins.reverse()
    else:
        b = first
        while b + (window - 1) * delta <= last + 1e-9:
            wins.append((b, b + (window - 1) * delta))
            b += step * delta
    if not wins:
        raise StataError(198, "window() too large for the sample")      # VERIFY
    keepvar, keep_start = None, False
    if isinstance(o.get("keep"), str):
        kv_, _, kopt = o["keep"].partition(",")
        keepvar, keep_start = kv_.strip(), kopt.strip().startswith("start")
        if not ds.has(keepvar):
            raise StataError(111, f"variable {keepvar} not found")
    out.write(f"(running {cmdname} on estimation sample)\n", "text")
    r0 = dict(s.r)
    nwin = len(wins)
    R = np.full((nwin, len(names)), np.nan)
    kcol = np.full(nwin, np.nan)
    dots = _dots(s, "Rolling replications", nwin, {"nodots": o.get("nodots") or noisy})
    reject = str(o["reject"]) if isinstance(o.get("reject"), str) else None
    try:
        for r, (wb, we) in enumerate(wins):
            sel = np.flatnonzero(base_mask & (t >= wb - 1e-9) & (t <= we + 1e-9))
            s.data = _subset(ds, sel)
            ok = True
            try:
                _quiet(s, command_, noisy)
                _, vals = _evaluate(s, items)
                if len(vals) != len(names):
                    ok = False
                else:
                    R[r] = vals
                    if reject is not None and float(s.eval(reject)) != 0:
                        R[r] = np.nan
            except StataError:
                ok = False
            if keepvar:
                kt = we if not keep_start else wb
                rows = np.flatnonzero(t == kt)
                if len(rows):
                    kcol[r] = to_nan(ds.get(keepvar).data)[rows[0]]
            if dots:
                _dot(s, ok and not np.isnan(R[r]).any(), r, nwin)
    finally:
        s.data = ds
    s.r = r0
    from ..core.dataset import Dataset, Variable

    res = Dataset()
    res.set_obs(nwin)
    tv = ds.get(info.tvar)
    for nm, vals in (("start", np.array([w[0] for w in wins])), ("end", np.array([w[1] for w in wins]))):
        res.add(Variable(nm, tv.vtype if tv.vtype in ("double", "long", "float") else "long", vals, fmt=tv.fmt))
    # _b and _se expand to one variable per coefficient (_b_x, _se_x); the others keep their names
    final = []
    idx = 0
    for nm0, ex in items:
        if ex in ("_b", "_se"):
            pref = "_b_" if ex == "_b" else "_se_"
            n_ex = len(est0.names) if est0 is not None else 0
            n_ex = len(s.e["b"].colnames) if "b" in s.e else n_ex
            for nm in names[idx:idx + n_ex]:
                col = nm.split(":")[-1]
                final.append(pref + (col[1:] if col == "_cons" else col))
            idx += n_ex
        else:
            final.append(nm0 or names[idx])
            idx += 1
    for j, nm in enumerate(final):
        col = np.where(np.isnan(R[:, j]), SYS, R[:, j])
        res.add(Variable(nm[:32], "double" if o.get("double") else "float",
                         col.astype(np.float32).astype(np.float64) if not o.get("double") else col))
    if keepvar:
        res.add(Variable(keepvar, "double", np.where(np.isnan(kcol), SYS, kcol)))
    if isinstance(o.get("saving"), str):
        base = s.data
        s.data = res
        try:
            file, _, opt = o["saving"].partition(",")
            s.interp._execute_expanded(f'save "{file.strip()}"' + (", replace" if "replace" in opt else ""))
        finally:
            s.data = base
    else:
        s.data = res
        s.notify_state()


@command("rolling", "rolling", prefix=True)
def _rolling(s: "Session", args: str) -> None:
    cmd_rolling(s, args)
