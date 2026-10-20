"""Model building: mkspline, xi, nestreg and stepwise ([R] mkspline, [R] xi, [R] nestreg,
[R] stepwise). The layouts follow the manuals; what they do not show has the VERIFY mark
(compat/do/1903)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import stats as S
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from ..lazy import lazy
from ._util import touse
from .matrix import Matrix
from .registry import command
from .summarize import g9, weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _sample_mask(s: "Session", p, varnames: list[str]):
    ds = s.data
    mask = touse(s, p)
    for v in varnames:
        var = ds.get(v)
        mask = mask & (np.array([x != "" for x in var.raw]) if var.is_string else var.data < SYS)
    return mask


def _add_var(s: "Session", name: str, values: np.ndarray, mask: np.ndarray, vtype: str = "float",
             label: str = "") -> None:
    from ..core.dataset import Variable
    full = np.full(s.data.nobs, SYS)
    full[mask] = values[mask] if len(values) == len(mask) else values
    s.data.add(Variable(name, vtype, full, label=label))


# ---------------------------------------------------------------------------
# mkspline
# ---------------------------------------------------------------------------

_CUBIC_PCT = {3: (10, 50, 90), 4: (5, 35, 65, 95), 5: (5, 27.5, 50, 72.5, 95),
              6: (5, 23, 41, 59, 77, 95), 7: (2.5, 18.33, 34.17, 50, 65.83, 81.67, 97.5)}


def _display_knots(s: "Session", knots: list[float]) -> None:
    out = s.output
    out.write("\n" + f"{'':>12} |" + "".join(f"{'knot' + str(i + 1):>11}" for i in range(len(knots))) + "\n", "text")
    out.write("-" * 13 + "+" + "-" * (11 * len(knots) + 1) + "\n", "text")      # VERIFY: layout
    out.write(f"{'knots':>12} |", "text")
    out.write("".join(f"{g9(k, '%9.0g'):>11}" for k in knots) + "\n", "result")


@command("mkspline", "mkspline")
def cmd_mkspline(s: "Session", args: str) -> None:
    p = parse_standard(args)
    c = find_top(args, "=")
    if c == -1:
        raise StataError(198, "invalid syntax")
    lhs = args[:c].split()
    rest = parse_standard(args[c + 1:])
    o = match_options(rest.options, {"marginal": 3, "displayknots": 1, "pctile": 3, "cubic": 3,
                                     "nknots": 2, "knots": 2}) if rest.options.strip() else {}
    ds = s.data
    xname = s.expand_varlist(rest.varlist)[0]
    x = ds.get(xname)
    if x.is_string:
        raise StataError(109, "type mismatch")
    mask = touse(s, rest) & (x.data < SYS)
    xv = x.data
    xs = np.sort(xv[mask])
    if xs.size == 0:
        raise StataError(2000, "no observations")
    created: list[str] = []
    if o.get("cubic"):
        if len(lhs) != 1:
            raise StataError(198, "invalid syntax")
        stub = lhs[0]
        if isinstance(o.get("knots"), str):
            knots = [float(v) for v in o["knots"].replace(",", " ").split()]
            k = len(knots)
        else:
            k = int(o["nknots"]) if isinstance(o.get("nknots"), str) else 5
            if k not in _CUBIC_PCT:
                raise StataError(198, "nknots() must be between 3 and 7")
            knots = [S.percentile(xs, q) for q in _CUBIC_PCT[k]]
        if k < 3 or k > 7 or sorted(knots) != knots:
            raise StataError(198, "knots() must have between 3 and 7 increasing values")
        tk, tk1, t1 = knots[-1], knots[-2], knots[0]
        pos = lambda v: np.where(v > 0, v, 0.0) ** 3
        cols = [xv.copy()]
        for j in range(k - 2):
            tj = knots[j]
            u = pos(xv - tj) - pos(xv - tk1) * (tk - tj) / (tk - tk1) + pos(xv - tk) * (tk1 - tj) / (tk - tk1)
            cols.append(u / (tk - t1) ** 2)
        names = [f"{stub}{j + 1}" for j in range(len(cols))]
    elif len(lhs) >= 3:
        # newvar1 # newvar2 # ... newvark
        names = lhs[0::2]
        knots = [float(v) for v in lhs[1::2]]
        if len(names) != len(knots) + 1 or sorted(knots) != knots:
            raise StataError(198, "invalid syntax")
        cols = _linear(xv, knots, bool(o.get("marginal")))
    else:
        # stubname # = oldvar: equal-spaced knots (or percentiles)
        if len(lhs) != 2:
            raise StataError(198, "invalid syntax")
        stub, nsp = lhs[0], int(float(lhs[1]))
        if nsp < 2:
            raise StataError(198, "number of spline variables must be at least 2")   # VERIFY
        if o.get("pctile"):
            knots = [S.percentile(xs, 100.0 * i / nsp) for i in range(1, nsp)]
        else:
            lo, hi = float(xs[0]), float(xs[-1])
            knots = [lo + (hi - lo) * i / nsp for i in range(1, nsp)]
        names = [f"{stub}{j + 1}" for j in range(nsp)]
        cols = _linear(xv, knots, bool(o.get("marginal")))
    for nm in names:
        if ds.has(nm):
            raise StataError(110, f"variable {nm} already defined")
    for nm, col in zip(names, cols):
        _add_var(s, nm, col, mask)
        created.append(nm)
    if o.get("displayknots"):
        _display_knots(s, knots)
    s.r = {"N_knots": float(len(knots)), "knots": Matrix(np.array(knots).reshape(1, -1), ["r1"],
                                                         [f"c{i + 1}" for i in range(len(knots))])}


def _linear(x: np.ndarray, knots: list[float], marginal: bool) -> list[np.ndarray]:
    """Linear spline basis: one variable per interval between the knots."""
    if marginal:
        return [x.copy()] + [np.maximum(x - kn, 0.0) for kn in knots]
    cols = [np.minimum(x, knots[0])]
    for j in range(1, len(knots)):
        cols.append(np.clip(x - knots[j - 1], 0.0, knots[j] - knots[j - 1]))
    cols.append(np.maximum(x - knots[-1], 0.0))
    return cols


# ---------------------------------------------------------------------------
# xi
# ---------------------------------------------------------------------------

_XI_TERM = re.compile(r"(?<![\w.])([iI])\.(\w+)(?:([*|])(?:([iI])\.)?(\w+))?")


def _levels_of(s: "Session", name: str):
    v = s.data.get(name)
    if v.is_string:
        raise StataError(109, "type mismatch")      # VERIFY: string variables are not supported yet
    x = v.data
    ok = x < SYS
    return x, ok, sorted(set(x[ok].tolist()))


def _lv(v: float) -> str:
    return f"{v:g}".replace("-", "n")


def _xi_dummies(s: "Session", prefix: str, var: str, omit: bool, made: list[str], out_names: list[str]):
    x, ok, levels = _levels_of(s, var)
    use = levels[1:] if omit else levels
    names = []
    for lv in use:
        nm = f"{prefix}{var}_{_lv(lv)}"
        col = np.where(ok, (x == lv).astype(float), SYS)
        if s.data.has(nm):
            s.data.drop_vars([nm])
        from ..core.dataset import Variable
        s.data.add(Variable(nm, "byte", col, label=f"{var}=={_lv(lv)}"))
        made.append(nm)
        names.append(nm)
    return x, ok, levels, use, names


def _xi_expand(s: "Session", text: str, prefix: str, omit_default: bool, made: list[str]):
    """Replaces the xi terms in `text` by the generated variables and returns (new text, notes)."""
    notes: list[str] = []

    def rep(m):
        flag, a, op, bflag, b = m.groups()
        omit = flag == "i" and omit_default
        x, ok, levels, use, names = _xi_dummies(s, prefix, a, omit, made, [])
        if levels:
            first = f"{prefix}{a}_{_lv(levels[0])}"
            rng = f"{names[0]}-{_lv(levels[-1])}" if names else ""
            notes.append((f"{flag}.{a}", rng, first if omit else None))
        if op is None:
            return " ".join(names)
        # interactions
        out_names = list(names) if op == "*" else []
        if bflag:
            xb, okb, levb, useb, namesb = _xi_dummies(s, prefix, b, omit, made, [])
            if op == "*":
                out_names += namesb
            for la in use:
                for lb in useb:
                    nm = f"{prefix}{a}X{b}_{_lv(la)}_{_lv(lb)}"
                    col = np.where(ok & okb, ((x == la) & (xb == lb)).astype(float), SYS)
                    from ..core.dataset import Variable
                    if s.data.has(nm):
                        s.data.drop_vars([nm])
                    s.data.add(Variable(nm, "byte", col))
                    made.append(nm)
                    out_names.append(nm)
        else:
            cv = s.data.get(b).data
            if op == "*":
                out_names.append(b)
            for la in use:
                nm = f"{prefix}{a}X{b}_{_lv(la)}"
                col = np.where(ok & (cv < SYS), (x == la) * cv, SYS)
                from ..core.dataset import Variable
                if s.data.has(nm):
                    s.data.drop_vars([nm])
                s.data.add(Variable(nm, "float", col))
                made.append(nm)
                out_names.append(nm)
        return " ".join(out_names)
    new = _XI_TERM.sub(rep, text)
    return new, notes


@command("xi", "xi", prefix=True)
def cmd_xi(s: "Session", args: str) -> None:
    colon = find_top(args, ":")
    head, cmdtext = (args, "") if colon == -1 else (args[:colon], args[colon + 1:].strip())
    c = find_top(head, ",")
    body, opts = (head, "") if c == -1 else (head[:c], head[c + 1:])
    o = match_options(opts, {"prefix": 3, "noomit": 3}) if opts.strip() else {}
    prefix = o["prefix"] if isinstance(o.get("prefix"), str) else "_I"
    for nm in getattr(s, "_xi_made", []):
        if s.data.has(nm):
            s.data.drop_vars([nm])
    made: list[str] = []
    text = cmdtext if colon != -1 else body
    new, notes = _xi_expand(s, text, prefix, not o.get("noomit"), made)
    s._xi_made = made
    out = s.output
    for term, rng, omitted in notes:
        tail = f"(naturally coded; {omitted} omitted)" if omitted else "(naturally coded; none omitted)"    # VERIFY
        out.write(f"{term:<18}{rng:<20}{tail}\n", "text")
    if colon != -1:
        s.interp._execute_expanded(new)


# ---------------------------------------------------------------------------
# helpers for nestreg and stepwise
# ---------------------------------------------------------------------------

def _split_prefix(args: str):
    colon = find_top(args, ":")
    if colon == -1:
        raise StataError(198, "invalid syntax")
    head, cmdtext = args[:colon], args[colon + 1:].strip()
    c = find_top(head, ",")
    return (head, "") if c == -1 else (head[:c], head[c + 1:]), cmdtext


def _paren_groups(text: str):
    """Splits `a (b c) d` into [['a'], ['b','c'], ['d']]; returns (groups, rest from the first token
    that is neither a name nor a group)."""
    groups, i, n = [], 0, len(text)
    while i < n:
        if text[i].isspace():
            i += 1
        elif text[i] == "(":
            depth, j = 0, i
            while j < n:
                depth += (text[j] == "(") - (text[j] == ")")
                if depth == 0:
                    break
                j += 1
            groups.append(text[i + 1:j].split())
            i = j + 1
        else:
            m = re.match(r"[^\s(]+", text[i:])
            tok = m.group(0)
            if tok in ("if", "in") or tok.startswith("[") or tok.startswith(","):
                break
            groups.append([tok])
            i += len(tok)
    return groups, text[i:]


class _Spec:
    """command depvar block ... [if] [in] [weight] [, options], ready to be rerun."""

    def __init__(self, s: "Session", cmdtext: str, with_dep: bool = True):
        c = find_top(cmdtext, ",")
        self.options = "" if c == -1 else cmdtext[c + 1:].strip()
        body = cmdtext if c == -1 else cmdtext[:c]
        words = body.split(None, 1)
        self.cmd = words[0]
        rest = words[1] if len(words) > 1 else ""
        self.dep = ""
        if with_dep:
            m = re.match(r"\s*(\S+)\s*(.*)$", rest, re.S)
            self.dep, rest = m.group(1), m.group(2)
        self.groups, tail = _paren_groups(rest)
        self.tail = parse_standard(tail) if tail.strip() else parse_standard("")
        self.s = s
        self.touse_name = ""
        self.mask = None

    def mark(self, names: list[str]) -> None:
        """Common sample: if/in/weight plus no missing values in any variable."""
        from ..commands.syntaxcmd import _temp_name
        from ..core.dataset import Variable
        s = self.s
        mask = _sample_mask(s, self.tail, ([self.dep] if self.dep else []) + names)
        if self.tail.weight is not None:
            w, wt, mask = weights(s, self.tail, mask, ("fweight", "aweight", "pweight", "iweight"))
        self.mask = mask
        self.touse_name = _temp_name(s)
        s.data.add(Variable(self.touse_name, "byte", mask.astype(float)))

    def release(self) -> None:
        if self.touse_name and self.s.data.has(self.touse_name):
            changed = self.s.data.changed
            self.s.data.drop_vars([self.touse_name])
            self.s.data.changed = changed
        self.touse_name = ""

    def text(self, varnames: list[str]) -> str:
        parts = [self.cmd, self.dep] + varnames
        t = self.tail
        if t.weight is not None:
            parts.append(f"[{t.weight[0]}={t.weight[1]}]")
        parts.append(f"if {self.touse_name} == 1")
        out = " ".join(x for x in parts if x)
        return out + (f", {self.options}" if self.options else "")


def _run_quiet(s: "Session", text: str, noisy: bool = False):
    if noisy:
        s.interp._execute_expanded(text)
    else:
        with s.output.quietly():
            s.interp._execute_expanded(text)
    est = getattr(s, "_estimates", None)
    if est is None or "b" not in s.e:
        raise StataError(301, "last estimates not found")
    return est


def _term_idx(est, term: list[str]) -> list[int]:
    idx = []
    for v in term:
        base = v.split(".")[-1]
        for i, nm in enumerate(est.names):
            c = nm.split(":")[-1]
            if c == v or c.endswith("." + base) and "." in c and not nm.startswith("o."):
                idx.append(i)
    return sorted(set(idx))


def _wald(est, idx: list[int]) -> tuple[float, float, float | None]:
    """(statistic, df, p): F with the residual df of the model, or chi2."""
    b = np.asarray(est.b, dtype=float)[idx]
    V = np.asarray(est.V, dtype=float)[np.ix_(idx, idx)]
    k = len(idx)
    if k == 0:
        return 0.0, 0.0, 1.0
    chi2 = float(b @ np.linalg.solve(V, b))
    if est.stat == "t" and est.df_r:
        f = chi2 / k
        return f, float(k), float(st.f.sf(f, k, est.df_r))
    return chi2, float(k), float(chi2_tail(chi2, k))


# ---------------------------------------------------------------------------
# nestreg
# ---------------------------------------------------------------------------

@command("nestreg", "nestreg", prefix=True)
def cmd_nestreg(s: "Session", args: str) -> None:
    (head, opts), cmdtext = _split_prefix(args)
    o = match_options(opts, {"quietly": 1, "lr": 2, "waldtable": 5, "lrtable": 3, "store": 5,
                             "terminator": 4, "forward": 3}) if opts.strip() else {}
    spec = _Spec(s, cmdtext)
    if not spec.groups:
        raise StataError(198, "no blocks specified")   # VERIFY
    names_all = [v for g in spec.groups for v in g]
    spec.mark(s.expand_varlist(" ".join(names_all)) if names_all else [])
    out = s.output
    lr = bool(o.get("lr") or o.get("lrtable"))
    rows = []
    prev = None
    try:
        so_far: list[str] = []
        for b, g in enumerate(spec.groups, start=1):
            so_far += g
            out.write(f"Block {b}: {' '.join(g)}\n", "text")
            est = _run_quiet(s, spec.text(so_far), noisy=not o.get("quietly"))
            e = dict(s.e)
            if isinstance(o.get("store"), str):
                s.interp._execute_expanded(f"estimates store {o['store']}{b}")
            if lr:
                rows.append((b, e.get("ll", SYS), prev, len(g), e))
            else:
                idx = _term_idx(est, s.expand_varlist(" ".join(g)))
                stat, df, pv = _wald(est, idx)
                rows.append((b, stat, df, pv, est.df_r, e))
            prev = e
    finally:
        spec.release()
    _nestreg_table(s, rows, lr, est.stat == "t" and bool(est.df_r))
    s.r = {}


def _nestreg_table(s: "Session", rows, lr: bool, is_f: bool) -> None:
    out = s.output
    out.write("\n", "text")
    if lr:
        hdr1 = "  | Block |        LL        LR     df  Pr > LR       AIC       BIC |"
        out.write("  +" + "-" * (len(hdr1) - 4) + "+\n", "text")
        out.write(hdr1 + "\n", "text")
        out.write("  |-------+" + "-" * (len(hdr1) - 12) + "|\n", "text")
        prev_ll = None
        for b, ll, prev, k, e in rows:
            n = e.get("N", SYS)
            rank = e.get("rank", e.get("k", SYS))
            aic = -2 * ll + 2 * rank
            bic = -2 * ll + rank * np.log(n)
            if prev_ll is None:
                lrs, dfs, pv = "", "", ""
            else:
                lrs_v = 2 * (ll - prev_ll[0])
                df_v = rank - prev_ll[1]
                lrs, dfs, pv = f"{lrs_v:.3f}", f"{int(df_v)}", f"{chi2_tail(lrs_v, df_v):.4f}"
            out.write(f"  | {b:>5} |{ll:>10.3f}{lrs:>10}{dfs:>7}{pv:>10}{aic:>10.3f}{bic:>10.3f} |\n", "result")
            prev_ll = (ll, rank)
        out.write("  +" + "-" * (len(hdr1) - 4) + "+\n", "text")
        return
    stat = "F" if is_f else "Chi2"
    if is_f:
        lines = [f"  | Block |{stat:>9}{'Block':>9}{'Residual':>10}{'':>9}{'':>9}{'Change':>8} |",
                 f"  |       |{'':>9}{'df':>9}{'df':>10}{'Pr > F':>9}{'R2':>9}{'in R2':>8} |"]
    else:
        lines = [f"  | Block |{stat:>9}{'Block':>9}{'':>9}{'Pseudo':>9}{'Change':>9} |",
                 f"  |       |{'':>9}{'df':>9}{'Pr > ' + stat:>9}{'R2':>9}{'in R2':>9} |"]
    width = len(lines[0]) - 4
    out.write("  +" + "-" * width + "+\n", "text")
    for ln in lines:
        out.write(ln + "\n", "text")
    out.write("  |-------+" + "-" * (width - 8) + "|\n", "text")
    prev_r2 = None
    for b, stv, df, pv, dfr, e in rows:
        r2 = e.get("r2", e.get("r2_p", SYS))
        ch = "" if prev_r2 is None else f"{r2 - prev_r2:.4f}"
        if is_f:
            out.write(f"  | {b:>5} |{stv:>9.2f}{int(df):>9}{int(dfr):>10}{pv:>9.4f}{r2:>9.4f}{ch:>8} |\n", "result")
        else:
            out.write(f"  | {b:>5} |{stv:>9.2f}{int(df):>9}{pv:>9.4f}{r2:>9.4f}{ch:>9} |\n", "result")
        prev_r2 = r2
    out.write("  +" + "-" * width + "+\n", "text")


# ---------------------------------------------------------------------------
# stepwise
# ---------------------------------------------------------------------------

@command("stepwise", "stepwise", prefix=True)
def cmd_stepwise(s: "Session", args: str) -> None:
    (head, opts), cmdtext = _split_prefix(args)
    o = match_options(opts, {"pr": 2, "pe": 2, "forward": 1, "hierarchical": 1, "lockterm1": 2}) \
        if opts.strip() else {}
    if not isinstance(o.get("pr"), str) and not isinstance(o.get("pe"), str):
        raise StataError(198, "pr() or pe() required")     # VERIFY
    forward = bool(o.get("forward")) or (isinstance(o.get("pe"), str) and not isinstance(o.get("pr"), str))
    pr = float(o["pr"]) if isinstance(o.get("pr"), str) else None
    pe = float(o["pe"]) if isinstance(o.get("pe"), str) else None
    if forward and pe is None:
        raise StataError(198, "forward stepwise requires pe()")      # VERIFY
    if not forward and pr is None:
        raise StataError(198, "backward stepwise requires pr()")     # VERIFY
    spec = _Spec(s, cmdtext)
    terms = [s.expand_varlist(" ".join(g)) if g else [] for g in spec.groups]
    if not terms:
        raise StataError(198, "no terms")
    spec.mark([v for t in terms for v in t])
    out = s.output
    hier = bool(o.get("hierarchical"))
    lock = bool(o.get("lockterm1"))

    def fit(idx: list[int], noisy=False):
        names = [v for i in idx for v in terms[i]]
        return _run_quiet(s, spec.text(names), noisy)

    def pvals(est, idx: list[int]) -> dict[int, float]:
        res = {}
        for pos, i in enumerate(idx):
            res[i] = _wald(est, _term_idx(est, terms[i]))[2]
        return res
    try:
        incl: list[int] = list(range(len(terms))) if not forward else ([0] if lock else [])
        text_label = "begin with full model" if not forward else "begin with empty model"
        out.write(f"{'':>22}{text_label}\n", "text")      # VERIFY: indentation
        changed = True
        while changed:
            changed = False
            if forward:
                best = None
                cand = [i for i in range(len(terms)) if i not in incl]
                if hier:
                    cand = cand[:1]
                for i in cand:
                    est = fit(sorted(incl + [i]))
                    pv = _wald(est, _term_idx(est, terms[i]))[2]
                    if best is None or pv < best[1]:
                        best = (i, pv)
                if best is not None and best[1] < pe:
                    incl = sorted(incl + [best[0]])
                    out.write(f"p = {best[1]:.4f} <  {pe:.4f}  adding   {' '.join(terms[best[0]])}\n", "text")
                    changed = True
                    if pr is not None and len(incl) > 1:
                        est = fit(incl)
                        pv = pvals(est, incl)
                        cand_r = [i for i in incl if not (lock and i == 0) and i != best[0]]
                        if hier:
                            cand_r = cand_r[-1:]
                        if cand_r:
                            worst = max(cand_r, key=lambda i: pv[i])
                            if pv[worst] >= pr:
                                incl.remove(worst)
                                out.write(f"p = {pv[worst]:.4f} >= {pr:.4f}  removing {' '.join(terms[worst])}\n", "text")
            else:
                if not incl:
                    break
                est = fit(incl)
                pv = pvals(est, incl)
                cand = [i for i in incl if not (lock and i == 0)]
                if hier:
                    cand = cand[-1:]
                if cand:
                    worst = max(cand, key=lambda i: pv[i])
                    if pv[worst] >= pr:
                        incl.remove(worst)
                        out.write(f"p = {pv[worst]:.4f} >= {pr:.4f}  removing {' '.join(terms[worst])}\n", "text")
                        changed = True
        fit(incl, noisy=True)
    finally:
        spec.release()
