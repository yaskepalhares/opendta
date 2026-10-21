"""Panel description commands: xtline, xttab and xtdata ([XT] xtline, [XT] xttab, [XT] xtdata).

xtline is translated into the equivalent twoway command; xttab and xtdata use the functions of
core/panelstats.py. Layouts that the manual does not show carry a VERIFY mark (compat 1905)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import panelstats as ps
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from .tsops import panel_info, ts_info

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _panel_vars(s: "Session", opt_i=None, opt_t=None, need_t: bool = False):
    ds = s.data
    pvar = opt_i or panel_info(ds)
    if not pvar:
        raise StataError(459, "panel variable not set; use xtset varname ...")
    info = ts_info(ds)
    tvar = opt_t or (info.tvar if info is not None and info.tvar != pvar else "")
    if need_t and not tvar:
        raise StataError(459, "time variable not set; use xtset panelvar timevar ...")
    return pvar, tvar


def cmd_xtline(s: "Session", args: str) -> None:
    from ..commands.graphcmd import _twoway
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    own = {"overlay": 3, "i": 1, "t": 1}
    # take our options out of the text and keep the rest for twoway
    keep, mine = [], {}
    from ..lang.syntax import parse_options
    raw = opt.strip()
    parsed = parse_options(raw) if raw else []
    rest_txt = raw
    for name, arg in parsed:
        for full, minimum in own.items():
            if len(name) >= minimum and full.startswith(name):
                mine[full] = arg if arg is not None else True
                txt = name + (f"({arg})" if arg is not None else "")
                rest_txt = rest_txt.replace(txt, "", 1)
                break
    pvar, tvar = _panel_vars(s, str(mine["i"]).strip() if isinstance(mine.get("i"), str) else None,
                             str(mine["t"]).strip() if isinstance(mine.get("t"), str) else None, need_t=True)
    p = parse_standard(body)
    names = s.expand_varlist(p.varlist)
    if not names:
        raise StataError(102, "too few variables specified")
    cond = ""
    if p.if_:
        cond = f" if {p.if_}"
    rest_txt = rest_txt.strip().strip(",").strip()
    if mine.get("overlay"):
        ds = s.data
        ids = ds.get(pvar).data
        vals = sorted(set(ids[ids < SYS].tolist()))
        plots = []
        for v in vals:
            e = f"{pvar}=={v:g}"
            cnd = f" if {e}" + (f" & ({p.if_})" if p.if_ else "")
            for n in names:
                plots.append(f"(line {n} {tvar}{cnd})")
        text = " ".join(plots) + (", " + rest_txt if rest_txt else "")
    else:
        plots = " ".join(f"(line {n} {tvar})" for n in names)
        text = plots + cond + f", by({pvar}, note(\"\"))" + (" " + rest_txt if rest_txt else "")
    _twoway(s, text)


def cmd_xttab(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..commands.tabulate import _categories, _codes, _fmt_freq, _header_lines
    p = parse_standard(args)
    o = match_options(p.options, {"i": 1}) if p.options.strip() else {}
    pvar, _ = _panel_vars(s, str(o["i"]).strip() if isinstance(o.get("i"), str) else None)
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103 if names else 102, "too many variables specified" if names else "varlist required")
    ds = s.data
    v = ds.get(names[0])
    mask = touse(s, p) & (ds.get(pvar).data < SYS)
    cats, labels, ok = _categories(s, v, mask, missing=False, nolabel=False)
    codes = _codes(v, cats)
    sel = ok & (codes >= 0)
    if not sel.any():
        raise StataError(2000, "no observations")
    _, inv = np.unique(ds.get(pvar).data[sel], return_inverse=True)
    r = ps.xttab_counts(codes[sel], inv, len(cats))
    # VERIFY: stub width 9 as in the 1905 log (wider labels widen it)
    title = _header_lines(v.label or v.name, 9)[-1]
    W = min(max([9, len(title)] + [len(lb) for lb in labels]), 30)
    out = s.output
    N, G = r["N"], r["G"]
    out.write("\n", "text")
    out.write(f"{'Overall':>{W + 16}}{'Between':>20}{'Within':>18}\n", "text")
    out.write(f"{title:>{W}} |{'Freq.':>9}{'Percent':>9}{'Freq.':>11}{'Percent':>9}{'Percent':>15}\n", "text")
    rule = "-" * (W + 1) + "+" + "-" * 53 + "\n"
    out.write(rule, "text")
    for k, lab in enumerate(labels):
        out.write(f"{lab[:W]:>{W}} |", "text")
        out.write(f"{_fmt_freq(r['freq'][k]):>8}{100 * r['freq'][k] / N:>10.2f}{_fmt_freq(r['bfreq'][k]):>10}"
                  f"{100 * r['bfreq'][k] / G:>10.2f}{r['within'][k]:>15.2f}\n", "result")
    out.write(rule, "text")
    # VERIFY: the total of the within column is the average of its percentages weighted by the between frequencies
    wtot = float((r["bfreq"] * r["within"]).sum() / r["bfreq"].sum())
    out.write(f"{'Total':>{W}} |", "text")
    out.write(f"{_fmt_freq(N):>8}{100.0:>10.2f}{_fmt_freq(r['bfreq'].sum()):>10}"
              f"{100 * r['bfreq'].sum() / G:>10.2f}{wtot:>15.2f}\n", "result")
    out.write(f"{'':>{W + 23}}(n = {G})\n", "text")
    s.r = {"N": N, "n": float(G), "r": float(len(cats))}


def cmd_xtdata(s: "Session", args: str) -> None:
    from ..commands._util import touse
    p = parse_standard(args)
    o = match_options(p.options, {"fe": 2, "be": 2, "re": 2, "fa": 2, "ratio": 3, "clear": 2, "i": 1}) \
        if p.options.strip() else {}
    modes = [m for m in ("fe", "be", "re", "fa") if o.get(m)]
    if len(modes) != 1:
        raise StataError(198, "one of fe, be, re or fa must be specified")                # VERIFY
    if modes[0] in ("re", "fa"):
        raise StataError(198, f"option {modes[0]} is not available yet in OpenDTA")
    pvar, _ = _panel_vars(s, str(o["i"]).strip() if isinstance(o.get("i"), str) else None)
    ds = s.data
    if ds.changed and not o.get("clear"):
        raise StataError(4, "no; dataset in memory has changed since last saved")     # VERIFY
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else [n for n in ds.names if n != pvar]
    mask = touse(s, p) & (ds.get(pvar).data < SYS)
    _, inv = np.unique(ds.get(pvar).data[mask], return_inverse=True)
    for n in names:
        v = ds.get(n)
        if v.is_string or n == pvar:
            continue
        x = v.data.astype(float)
        m = mask & (x < SYS)
        sub_inv = np.unique(ds.get(pvar).data[m], return_inverse=True)[1]
        new = x.copy()
        new[m] = ps.panel_transform(x[m], sub_inv, modes[0])
        ds.set_numeric(v, new)
    if modes[0] == "be":
        # VERIFY: the between-effects data keep one observation per panel (1905 log)
        ids_ = ds.get(pvar).data
        first = np.zeros(ds.nobs, dtype=bool)
        _, idx = np.unique(np.where(mask, ids_, np.nan), return_index=True)
        first[idx] = True
        ds.keep_obs(first & mask)
    s.notify_state()
