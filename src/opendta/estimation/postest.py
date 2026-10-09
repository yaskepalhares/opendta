"""Pós-estimação: replay, predict, test, testparm, lincom e estimates
([R] predict, [R] test, [R] testparm, [R] lincom, [R] estimates)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from . import fvars
from .results import CoefRow, Estimates, coef_table, current, g9

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def replay(s: "Session", cmd: str, args: str) -> None:
    est = current(s)
    if est.cmd != cmd:
        raise StataError(301, "last estimates not found")
    o = match_options(args.lstrip(", "), {"level": 1, "noheader": 6, "notable": 5}) \
        if args.strip(" ,") else {}
    lev = float(o["level"]) if o.get("level") else None
    _display(s, est, level=lev, header=not o.get("noheader"), table=not o.get("notable"))


def _display(s: "Session", est: Estimates, **kw) -> None:
    if est.cmd == "regress":
        from .regress import display_regress
        display_regress(s, est, **kw)
    else:
        disp = getattr(est, "display", None)
        if disp is None:
            raise StataError(301, "last estimates not found")
        disp(s, est, **kw)


# ---------------------------------------------------------------------------
# matriz de delineamento para predict (todas as observações)
# ---------------------------------------------------------------------------

def _columns(s: "Session", est: Estimates):
    """Gera (índice em e(b), coluna) para todas as observações, sem montar X."""
    ds = s.data
    n = ds.nobs
    for j, c in enumerate(est.cols):
        if c.values is None:
            continue
        t = est.terms[c.term]
        y = np.ones(n)
        k = 0
        for comp in t.comps:
            x = fvars.comp_values(ds, comp)
            if comp.factor:
                y = y * (x == c.cell[k])
                k += 1
            else:
                y = y * np.where(x < SYS, x, 0.0)
        yield j, y
    if est.constant:
        yield len(est.cols), np.ones(n)


def design_all(s: "Session", est: Estimates) -> tuple[np.ndarray, np.ndarray]:
    """X para todas as observações (na ordem de e(b)) e máscara de não missing."""
    ds = s.data
    n = ds.nobs
    with fvars.caching():
        ok = ~fvars.term_missing(ds, est.terms)
        k = len(est.b)
        X = np.zeros((n, k))
        for j, col in _columns(s, est):
            X[:, j] = col
    return X, ok


def xb_all(s: "Session", est: Estimates) -> tuple[np.ndarray, np.ndarray]:
    """xb para todas as observações, acumulado coluna a coluna (sem guardar X)."""
    ds = s.data
    with fvars.caching():
        ok = ~fvars.term_missing(ds, est.terms)
        xb = np.zeros(ds.nobs)
        for j, col in _columns(s, est):
            if est.b[j] != 0:
                xb += est.b[j] * col
    return xb, ok


# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

_PRED_OPTS = {"xb": 2, "residuals": 1, "stdp": 4, "stdf": 4, "stdr": 4, "hat": 3, "leverage": 3,
              "rstandard": 3, "rstudent": 4, "cooksd": 4, "pr": 2, "p": 1, "n": 1, "ir": 2,
              "score": 5, "equation": 2, "nooffset": 5, "scores": 5}


def cmd_predict(s: "Session", args: str) -> None:
    est = current(s)
    from ..commands._util import touse
    p = parse_standard(args)
    o = match_options(p.options, _PRED_OPTS) if p.options.strip() else {}
    words = p.varlist.split()
    if not words:
        raise StataError(100, "varlist required")
    vtype = "float"
    if words[0] in ("byte", "int", "long", "float", "double"):
        vtype = words.pop(0)
    if len(words) != 1:
        raise StataError(103, "too many variables specified")
    name = words[0]
    ds = s.data
    if ds.has(name):
        raise StataError(110, f"{name} already defined")
    mask = touse(s, p)
    stats = [k for k in ("xb", "residuals", "stdp", "stdf", "stdr", "hat", "leverage", "rstandard",
                         "rstudent", "cooksd", "pr", "p", "n", "ir") if o.get(k)]
    if len(stats) > 1:
        raise StataError(198, "only one statistic may be specified")   # VERIFICAR
    out = s.output
    pred = getattr(est, "predict", None)
    if pred is not None:
        values, note = pred(s, est, stats[0] if stats else "", mask)
    else:
        values, note = _predict_regress(s, est, stats[0] if stats else "", mask)
    if note:
        out.write(note + "\n", "text")
    values = np.where(mask & (values < SYS) & np.isfinite(values), values, SYS)
    from ..core.dataset import Variable
    var = Variable(name, vtype, np.full(ds.nobs, SYS))
    var.data = values
    ds.add(var)
    nmiss = int(np.sum(var.data >= SYS))
    if nmiss:
        out.write(f"({nmiss:,} missing value{'s' if nmiss != 1 else ''} generated)\n", "text")
    s.notify_state()


def _predict_regress(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    xb, ok = xb_all(s, est)
    xb = np.where(ok, xb, SYS)
    note = ""
    if stat in ("", "xb"):
        if not stat:
            note = "(option xb assumed; fitted values)"
        return xb, note
    ds = s.data
    yv = fvars.comp_values(ds, est.depcomp)
    if stat == "residuals":
        return np.where(ok & (yv < SYS), yv - xb, SYS), ""
    X, _ = design_all(s, est)
    V = est.V
    if stat == "residuals":
        return np.where(ok & (yv < SYS), yv - xb, SYS), ""
    sp = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", X, V, X), 0))
    if stat == "stdp":
        return np.where(ok, sp, SYS), ""
    s2 = est.s2
    if stat == "stdf":
        return np.where(ok, np.sqrt(sp ** 2 + s2), SYS), ""
    # leverage: h = x (X'WX)^-1 x' = stdp^2 / s2 (pesos: VERIFICAR)
    h = np.where(ok, sp ** 2 / s2, SYS)
    if stat in ("hat", "leverage"):
        return h, ""
    e = np.where(ok & (yv < SYS), yv - xb, SYS)
    good = ok & (yv < SYS) & (h < 1)
    stdr = np.where(good, np.sqrt(np.maximum(s2 * (1 - h), 0)), SYS)
    if stat == "stdr":
        return stdr, ""
    rstd = np.where(good, e / np.where(good, stdr, 1), SYS)
    if stat == "rstandard":
        return rstd, ""
    dfr = est.extra["df_r"]
    if stat == "rstudent":
        r = np.where(good, rstd, 0.0)
        val = r * np.sqrt((dfr - 1) / np.maximum(dfr - r ** 2, 1e-300))
        return np.where(good, val, SYS), ""
    if stat == "cooksd":
        k = int(est.extra["rank"])
        val = np.where(good, rstd, 0.0) ** 2 * np.where(good, h, 0) / (k * np.where(good, 1 - h, 1))
        return np.where(good & esample(s), val, SYS), ""
    raise StataError(198, f"option {stat} not allowed")


def esample(s: "Session") -> np.ndarray:
    m = s.data.esample
    if m is None or len(m) != s.data.nobs:
        return np.zeros(s.data.nobs, dtype=bool)
    return m


# ---------------------------------------------------------------------------
# restrições lineares (test, lincom)
# ---------------------------------------------------------------------------

def _coef_index(est: Estimates, s: "Session", name: str, eq: str = "") -> int:
    """Índice em e(b) de um nome como 'mpg', '2.rep78', '_cons', 'L.x'."""
    name = name.strip()
    names = est.names
    eqs = est.eqnames or [""] * len(names)

    def match(target: str) -> int | None:
        for i, n in enumerate(names):
            if (not eq or eqs[i] == eq) and _same(n, target):
                return i
        return None

    i = match(name)
    if i is not None:
        return i
    # abreviação/nome de variável sem operadores
    try:
        from ..core.varlist import resolve_name
        full = resolve_name(s.data, name)
        i = match(full)
        if i is not None:
            return i
    except StataError:
        pass
    m = re.match(r"^(\d+)\.(\w+)$", name)
    if m:
        try:
            from ..core.varlist import resolve_name
            full = resolve_name(s.data, m.group(2))
            i = match(f"{m.group(1)}.{full}")
            if i is not None:
                return i
        except StataError:
            pass
    raise StataError(111, f"{name} not found")


def _strip_marks(n: str) -> str:
    out = []
    for p in n.split("#"):
        p = re.sub(r"^(\d+)(?:b|o|bn)+\.", r"\1.", p)
        p = re.sub(r"^(?:o|co)\.", lambda m: "c." if m.group(0) == "co." else "", p)
        out.append(p)
    return "#".join(out)


def _same(a: str, b: str) -> bool:
    return _strip_marks(a) == _strip_marks(b)


class _Lin:
    """Expressão linear nos coeficientes: {índice: coeficiente} + constante."""

    def __init__(self, s: "Session", est: Estimates, text: str):
        self.s, self.est = s, est
        self.toks = re.findall(r"_b\[[^\]]*\]|_coef\[[^\]]*\]|\[[^\]]*\]\s*_?[\w.#]+|"
                               r"\d+[bno]*\.[A-Za-z_][\w.#]*|"
                               r"\d+\.\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|\d+(?:[eE][-+]?\d+)?|"
                               r"[A-Za-z_][\w.#~]*|[-+*/()]", text)
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else None

    def take(self):
        t = self.peek()
        self.i += 1
        return t

    def parse(self) -> tuple[dict, float]:
        r = self.sum()
        if self.peek() is not None:
            raise StataError(198, "invalid syntax")
        return r

    def sum(self):
        coef, const = self.prod()
        while self.peek() in ("+", "-"):
            op = self.take()
            c2, k2 = self.prod()
            sgn = 1 if op == "+" else -1
            for kk, v in c2.items():
                coef[kk] = coef.get(kk, 0) + sgn * v
            const += sgn * k2
        return coef, const

    def prod(self):
        coef, const = self.unary()
        while self.peek() in ("*", "/"):
            op = self.take()
            c2, k2 = self.unary()
            if op == "*":
                if not coef:
                    coef, const = {kk: v * const for kk, v in c2.items()}, const * k2
                elif not c2:
                    coef, const = {kk: v * k2 for kk, v in coef.items()}, const * k2
                else:
                    raise StataError(131, "not possible with test")   # VERIFICAR
            else:
                if c2:
                    raise StataError(131, "not possible with test")
                coef, const = {kk: v / k2 for kk, v in coef.items()}, const / k2
        return coef, const

    def unary(self):
        t = self.peek()
        if t == "-":
            self.take()
            c, k = self.unary()
            return {kk: -v for kk, v in c.items()}, -k
        if t == "+":
            self.take()
            return self.unary()
        return self.atom()

    def atom(self):
        t = self.take()
        if t is None:
            raise StataError(198, "invalid syntax")
        if t == "(":
            r = self.sum()
            if self.take() != ")":
                raise StataError(198, "invalid syntax")
            return r
        m = re.fullmatch(r"_(?:b|coef)\[([^\]]*)\]", t)
        if m:
            return self._coef(m.group(1)), 0.0
        m = re.fullmatch(r"\[([^\]]*)\]\s*(\S+)", t)
        if m:
            return self._coef(m.group(2), m.group(1)), 0.0
        if re.fullmatch(r"\d+\.\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|\d+(?:[eE][-+]?\d+)?", t):
            return {}, float(t)
        return self._coef(t), 0.0

    def _coef(self, name: str, eq: str = "") -> dict:
        if ":" in name and not eq:
            eq, name = name.split(":", 1)
        return {_coef_index(self.est, self.s, name, eq.strip()): 1.0}


def _constraint(s: "Session", est: Estimates, text: str) -> tuple[dict, float]:
    """'a = b' ou 'a' (= 0) -> R (dict) e r, com R b = r."""
    if "=" in text:
        lhs, rhs = text.split("=", 1)
        c1, k1 = _Lin(s, est, lhs).parse()
        c2, k2 = _Lin(s, est, rhs).parse() if rhs.strip() else ({}, 0.0)
        coef = dict(c1)
        for k, v in c2.items():
            coef[k] = coef.get(k, 0) - v
        return coef, k2 - k1
    c1, k1 = _Lin(s, est, text).parse()
    return c1, -k1


def _ctext(est: Estimates, coef: dict, r: float) -> str:
    """' mpg - weight = 0' como o Stata escreve a restrição."""
    parts = []
    for k in sorted(coef):
        v = coef[k]
        if v == 0:
            continue
        nm = _strip_marks(est.names[k]) if not est.names[k].startswith("o.") else est.names[k]
        if est.eqnames and any(e not in ("", "_") for e in est.eqnames):
            nm = f"[{est.eqnames[k]}]{nm}"   # modelos ml: [y]x1 (compat 0506)
        a = abs(v)
        term = nm if a == 1 else f"{format_value(a, '%9.0g', pad=False).strip()}*{nm}"
        if not parts:
            parts.append(("- " if v < 0 else "") + term)
        else:
            parts.append(("- " if v < 0 else "+ ") + term)
    lhs = " ".join(parts) if parts else "0"
    rhs = format_value(r, "%9.0g", pad=False).strip() if r != 0 else "0"
    return f"{lhs} = {rhs}"


def _wald(s: "Session", est: Estimates, cons: list[tuple[dict, float]], *, header: bool = True,
          mtest: bool = False) -> None:
    out = s.output
    k = len(est.b)
    R = np.zeros((len(cons), k))
    r = np.zeros(len(cons))
    for i, (coef, rv) in enumerate(cons):
        for j, v in coef.items():
            R[i, j] = v
        r[i] = rv
    if header:
        out.write("\n", "text")
        for i, (coef, rv) in enumerate(cons):
            out.write(f" ({i + 1:>2})  {_ctext(est, coef, rv)}\n", "text")
    # restrições dependentes/omitidas
    V = est.V
    d = R @ est.b - r
    RVR = R @ V @ R.T
    notes = []
    keep = []
    for i in range(len(cons)):
        if np.allclose(R[i] @ V @ R[i], 0) and np.allclose(d[i], 0):
            notes.append(i + 1)
        else:
            keep.append(i)
    # dependência linear entre as restrições
    sel = []
    for i in keep:
        trial = sel + [i]
        sub = RVR[np.ix_(trial, trial)]
        if np.linalg.matrix_rank(sub, tol=1e-12 * max(1.0, np.abs(sub).max())) == len(trial):
            sel.append(i)
        else:
            notes.append(i + 1)
    q = len(sel)
    for i in sorted(notes):
        out.write(f"       Constraint {i} dropped\n", "text")   # VERIFICAR
    if q == 0:
        out.write("\n", "text")
        raise StataError(111, "no restrictions")   # VERIFICAR
    W = float(d[sel] @ np.linalg.solve(RVR[np.ix_(sel, sel)], d[sel]))
    out.write("\n", "text")
    s.r = {}
    if est.stat == "t" and est.df_r:
        F = W / q
        dfr = est.df_r
        p = float(st.f.sf(F, q, dfr))
        out.write(f"{'F(' + format(q, '>3') + ',' + format(int(dfr), '>6') + ') =':>22}", "text")
        out.write(f"{F:>8.2f}\n", "result")
        out.write(f"{'Prob > F =':>22}", "text")
        out.write(f"{p:>10.4f}\n", "result")
        s.r = {"drop": 0.0, "df_r": float(dfr), "F": F, "df": float(q), "p": p}
    else:
        p = float(st.chi2.sf(W, q))
        out.write(f"{'chi2(' + format(q, '>3') + ') =':>22}", "text")
        out.write(f"{W:>8.2f}\n", "result")
        out.write(f"{'Prob > chi2 =':>22}", "text")
        out.write(f"{p:>10.4f}\n", "result")
        s.r = {"drop": 0.0, "chi2": W, "df": float(q), "p": p}
    s.r = dict(reversed(list(s.r.items())))


def _split_parens(text: str) -> list[str]:
    out, depth, cur = [], 0, ""
    for ch in text:
        if ch == "(":
            depth += 1
            if depth == 1:
                cur = ""
                continue
        if ch == ")":
            depth -= 1
            if depth == 0:
                out.append(cur)
                continue
        if depth >= 1:
            cur += ch
    return out


def cmd_test(s: "Session", args: str) -> None:
    est = current(s)
    head, comma, opts = args.partition(",")
    o = match_options(opts, {"accumulate": 3, "notest": 5, "mtest": 2, "coef": 4}) \
        if comma and opts.strip() else {}
    text = head.strip()
    if not text:
        raise StataError(198, "invalid syntax")
    if text.startswith("("):
        exprs = _split_parens(text)
    elif "=" in text or re.search(r"[-+*/]|_b\[|\[", text):
        exprs = [text]
    else:
        # test varlist: cada coeficiente = 0 (aceita fatoriais e curingas)
        exprs = _coef_names(s, est, text)
    multi = est.eqnames and len(set(est.eqnames)) > 1
    if multi and not text.startswith("(") and "=" not in text and not re.search(r"[-+*/\[]", text):
        # várias equações: a variável em todas as equações em que é estimada
        # (VERIFICAR a ordem e as equações de base do mlogit)
        cons = []
        for e in exprs:
            for i, (nm, eq) in enumerate(zip(est.names, est.eqnames)):
                if _same(nm, e):
                    cons.append(({i: 1.0}, 0.0))   # a da base sai como "dropped"
        if not cons:
            raise StataError(111, f"{exprs[0]} not found")
    else:
        cons = [_constraint(s, est, e) for e in exprs]
    if o.get("accumulate"):
        cons = getattr(s, "_test_accum", []) + cons
    s._test_accum = cons
    if o.get("notest"):
        out = s.output
        out.write("\n", "text")
        for i, (coef, rv) in enumerate(cons):
            out.write(f" ({i + 1:>2})  {_ctext(est, coef, rv)}\n", "text")
        return
    _wald(s, est, cons)


def _coef_names(s: "Session", est: Estimates, text: str) -> list[str]:
    """Nomes de coeficientes a partir de uma varlist (com fatoriais)."""
    names = []
    for tok in text.split():
        if re.match(r"^\d+\.\w+$", tok):
            names.append(tok)
            continue
        if "#" in tok or re.match(r"^[ibcLFDS]\w*\.", tok, re.I):
            terms = fvars.expand_fv(s.data, tok)
            for c in est.cols:
                t = est.terms[c.term] if c.term < len(est.terms) else None
                if t is not None and any(t.key() == tt.key() for tt in terms) and not c.base:
                    names.append(c.name)
            if not names and re.match(r"^\d+\.", tok):
                names.append(tok)
            continue
        names.append(tok)
    return names


def cmd_testparm(s: "Session", args: str) -> None:
    est = current(s)
    head, comma, opts = args.partition(",")
    o = match_options(opts, {"equal": 2, "equation": 2, "nosvyadjust": 5}) \
        if comma and opts.strip() else {}
    terms = fvars.expand_fv(s.data, head)
    keys = {t.key() for t in terms}
    idx = []
    for i, c in enumerate(est.cols):
        if est.terms[c.term].key() in keys and not c.base and not c.omitted:
            idx.append(i)
    if not idx:
        raise StataError(111, "no such variables")   # VERIFICAR
    if o.get("equal"):
        cons = []
        for j in idx[1:]:
            cons.append(({idx[0]: -1.0, j: 1.0}, 0.0))
    else:
        cons = [({i: 1.0}, 0.0) for i in idx]
    _wald(s, est, cons)


def cmd_lincom(s: "Session", args: str) -> None:
    est = current(s)
    head, comma, opts = args.partition(",")
    o = match_options(opts, {"level": 1, "eform": 2, "or": 2, "hr": 2, "irr": 2, "rrr": 2}) \
        if comma and opts.strip() else {}
    coef, r = _constraint(s, est, head)
    # lincom estima a combinação (sem o lado direito): R b - r
    out = s.output
    out.write("\n", "text")
    out.write(f" ( 1)  {_ctext(est, coef, r)}\n", "text")
    out.write("\n", "text")
    R = np.zeros(len(est.b))
    for j, v in coef.items():
        R[j] = v
    val = float(R @ est.b - r)
    se = float(np.sqrt(max(R @ est.V @ R, 0)))
    lev = float(o["level"]) if o.get("level") else getattr(est, "level", 95.0)
    eform = next((k for k in ("or", "hr", "irr", "rrr", "eform") if o.get(k)), "")
    title = {"or": "Odds Ratio", "hr": "Haz. Ratio", "irr": "IRR", "rrr": "RRR",
             "eform": "exp(b)"}.get(eform, "Coef.")
    if eform:
        from .results import crit, pvalue
        c = crit(est.stat, lev, est.df_r)
        tv = val / se if se > 0 else SYS
        p = pvalue(est.stat, tv, est.df_r)
        ev, ese = np.exp(val), np.exp(val) * se
        _eform_line(s, est, ev, ese, tv, p, np.exp(val - c * se), np.exp(val + c * se), lev, title)
    else:
        coef_table(s, est.depvar, [CoefRow("coef", "(1)", val, se)], stat=est.stat, df=est.df_r,
                   level=lev)
    # r(): só estimate, se e df no Stata 14 (compat 0504); a lista sai invertida
    s.r = {"estimate": val, "se": se}
    if est.df_r:
        s.r["df"] = float(est.df_r)


def _eform_line(s, est, ev, ese, tv, p, lo, hi, lev, title):
    from .results import LINE
    out = s.output
    stat_h = "t    P>|t|" if est.stat == "t" else "z    P>|z|"
    out.write(LINE + "\n", "text")
    out.write(f"{est.depvar[:12]:>12} | {title:>10}   Std. Err.      {stat_h}     "
              f"[{lev:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    out.write(f"{'(1)':>12} |  {g9(ev):>9}  {g9(ese):>9}  {tv:>7.2f}  {p:>6.3f}    {g9(lo):>9}   "
              f"{g9(hi):>9}\n", "result")
    out.write(LINE + "\n", "text")


# ---------------------------------------------------------------------------
# estimates store / restore / dir / drop / clear / replay / table
# ---------------------------------------------------------------------------

def _stored(s: "Session") -> dict:
    if not hasattr(s, "_stored_estimates"):
        s._stored_estimates = {}
    return s._stored_estimates


def _snapshot(s: "Session") -> dict:
    import copy
    return {"e": copy.deepcopy(s.e), "est": copy.deepcopy(getattr(s, "_estimates", None)),
            "sample": None if s.data.esample is None else s.data.esample.copy()}


def _restore(s: "Session", snap: dict) -> None:
    import copy
    s.e = copy.deepcopy(snap["e"])
    s._estimates = copy.deepcopy(snap["est"])
    s.data.esample = None if snap["sample"] is None else snap["sample"].copy()


def cmd_estimates(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    st_ = _stored(s)
    out = s.output
    if sub in ("store", "sto"):
        head, comma, opts = rest.partition(",")
        name = head.strip()
        if not name:
            raise StataError(198, "invalid syntax")
        if not re.fullmatch(r"[A-Za-z_]\w*", name):
            raise StataError(198, f"{name} invalid name")
        if not s.e or "cmd" not in s.e:
            raise StataError(301, "last estimates not found")
        snap = _snapshot(s)
        st_.pop(name, None)
        st_[name] = snap
        return
    if sub in ("restore", "res"):
        name = rest.split(",")[0].strip()
        if name not in st_:
            raise StataError(111, f"estimation result {name} not found")
        _restore(s, st_[name])
        s.e["_estimates_name"] = name   # VERIFICAR: e(_estimates_name)
        s.output.write(f"(results {name} are active now)\n", "text")
        return
    if sub in ("dir",):
        if not st_:
            return
        _est_dir(s, st_)
        return
    if sub in ("drop",):
        names = rest.split()
        if names == ["_all"]:
            st_.clear()
            return
        for n in names:
            if n not in st_:
                raise StataError(111, f"estimation result {n} not found")
            del st_[n]
        return
    if sub in ("clear",):
        st_.clear()
        return
    if sub in ("replay", "rep"):
        names = rest.split(",")[0].split()
        saved = _snapshot(s)
        try:
            for n in names or ["."]:
                if n != ".":
                    if n not in st_:
                        raise StataError(111, f"estimation result {n} not found")
                    _restore(s, st_[n])
                    w = 78
                    lab = f"Model {n}"
                    out.write("\n" + "-" * w + "\n", "text")
                    out.write(lab + "\n", "result")
                    out.write("-" * w + "\n", "text")
                est = current(s)
                _display(s, est)
        finally:
            _restore(s, saved)
        return
    if sub in ("table", "tab"):
        from .esttable import estimates_table
        estimates_table(s, rest)
        return
    raise StataError(198, f"{sub}: unknown subcommand")   # VERIFICAR


def _est_dir(s: "Session", st_: dict) -> None:
    """Layout do Stata 14 (compat 0504): a coluna title é o título dado por
    estimates title, não e(title)."""
    out = s.output
    out.write("\n", "text")
    out.write("-" * 55 + "\n", "text")
    out.write(f"{'name':>12} | {'command':<12} {'depvar':<12} {'npar':>4}  title \n", "text")
    out.write("-" * 13 + "+" + "-" * 41 + "\n", "text")
    for name, snap in st_.items():
        e = snap["e"]
        b = e.get("b")
        npar = b.cols if b is not None else 0
        out.write(f"{name:>12} | ", "text")
        out.write(f"{e.get('cmd', ''):<12} {e.get('depvar', ''):<12} {npar:>4}  "
                  f"{snap.get('title', '')}\n", "result")
    out.write("-" * 55 + "\n", "text")
