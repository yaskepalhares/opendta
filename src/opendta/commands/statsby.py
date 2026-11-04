"""statsby ([D] statsby): runs a command on each group of the data and collects the results it
returns as a dataset."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Dataset, Variable
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lang.words import strip_outer_quotes
from .datacmds import _need_clear
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


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
        # the command inside is silent even under `capture noisily statsby` (compat 1501)
        noisy_saved = out.noisy_depth
        if not o.get("noisily"):
            out.quiet_depth += 1
            out.noisy_depth = 0
        try:
            s.interp._execute_expanded(cmd)
        except StataError:
            ok = False
        finally:
            if not o.get("noisily"):
                out.quiet_depth -= 1
                out.noisy_depth = noisy_saved
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
