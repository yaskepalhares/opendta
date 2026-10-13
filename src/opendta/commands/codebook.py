"""codebook e misstable ([D] codebook, [R] misstable).

Layouts reconstruídos de saídas do Stata 14; o que não pôde ser conferido
tem a marca VERIFICAR e casos em compat/do/0803_codebook_bootstrap.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import stats as S
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.functions import _abbrev
from ..lang.syntax import match_options, parse_standard
from ._util import touse
from .registry import command

if TYPE_CHECKING:
    from ..core.dataset import Dataset, Variable
    from ..session import Session

SYS = M.SYSMISS
def _width(s: "Session") -> int:
    """O codebook ocupa a linha toda: c(linesize) (compat 0803)."""
    return int(s.settings.get("linesize", 80))


def _g(x: float) -> str:
    """Média e desvio-padrão do codebook: 6 algarismos significativos
    (VERIFICAR; reproduz '21.2973' e '5.7855' da documentação)."""
    if x >= SYS or not math.isfinite(x):
        return "."
    t = f"{x:.6g}"
    if "e" in t:
        return format_value(x, "%9.0g", pad=False).strip()
    if t.startswith("0."):
        t = t[1:]
    elif t.startswith("-0."):
        t = "-" + t[2:]
    return t


def _num(x: float) -> str:
    return format_value(float(x), "%9.0g", pad=False).strip()


def _units(vals: np.ndarray, tol: float = 1e-9) -> str:
    """Maior potência de 10 que divide todos os valores (VERIFICAR: o Stata
    também usa unidades como 2 ou 5?)."""
    if len(vals) == 0:
        return "."
    top = float(np.max(np.abs(vals))) or 1.0
    for k in range(int(math.floor(math.log10(top))), -9, -1):
        u = 10.0 ** k
        q = np.round(vals / u)
        if np.all(np.abs(vals - q * u) <= tol * np.maximum(1.0, np.abs(vals))):
            return _num(u)
    return "."


def _header(s: "Session", v: "Variable") -> None:
    out = s.output
    W = _width(s)
    out.write("-" * W + "\n", "text")
    lab = v.label or "(unlabeled)"
    space = W - len(v.name) - len(lab)
    out.write(v.name, "result")
    out.write(" " * max(space, 1) + lab + "\n", "text")
    out.write("-" * W + "\n", "text")
    out.write("\n", "text")


def _line(s: "Session", label: str, value: str, right_label: str = "", right_value: str = "") -> None:
    out = s.output
    out.write(f"{label + ':':>23}  ", "text")
    if right_label:
        out.write(f"{value:<25}", "result")
        out.write(f"{right_label + ':':>10}  ", "text")
        out.write(right_value + "\n", "result")
    else:
        out.write(value + "\n", "result")


def _codebook_var(s: "Session", v: "Variable", mask: np.ndarray, tab_max: int, o: dict) -> None:
    ds = s.data
    out = s.output
    _header(s, v)
    N = int(mask.sum())
    if v.is_string:
        col = np.asarray(v.raw)[mask]
        longest = max((len(x) for x in col), default=0)
        decl = int(v.vtype[3:]) if v.vtype.startswith("str") and v.vtype[3:].isdigit() else 0
        extra = f", but longest is str{max(longest, 1)}" if decl and longest < decl else ""
        _line(s, "type", f"string ({v.vtype}){extra}")
        out.write("\n", "text")
        nonmiss = col[col != ""]
        uniq = sorted(set(nonmiss.tolist()))
        _line(s, "unique values", str(len(uniq)), 'missing ""', f"{N - len(nonmiss)}/{N}")
        out.write("\n", "text")
        if 0 < len(uniq) <= tab_max:
            out.write(f"{'tabulation:':>23}  ", "text")
            out.write("Freq.  Value\n", "text")
            for u in uniq:
                out.write(" " * 25 + f"{int((nonmiss == u).sum()):>5}", "result")
                out.write(f"  \"{u}\"\n", "result")
        elif len(uniq) > 0:
            srt = sorted(nonmiss.tolist())
            n = len(srt)
            # exemplos: posições floor(n·k/5), k = 1..4, na lista ordenada
            # (compat 0803)
            picks = [srt[min(n - 1, n * k // 5)] for k in (1, 2, 3, 4)]
            for k, ex in enumerate(picks):
                lab = "examples:" if k == 0 else ""
                out.write(f"{lab:>23}  ", "text")
                out.write(f"\"{ex}\"\n", "result")
        if any(" " in x for x in uniq):
            out.write("\n", "text")
            out.write(f"{'warning:':>23}  ", "text")
            out.write("variable has embedded blanks\n", "text")
        out.write("\n", "text")
        return
    x = v.data[mask]
    valid = x[x < SYS]
    _line(s, "type", f"numeric ({v.vtype})")
    if v.value_label:
        _line(s, "label", v.value_label + ("" if v.value_label in ds.value_labels else ", but label does not exist"))
    out.write("\n", "text")
    if len(valid):
        rng = f"[{_num(valid.min())},{_num(valid.max())}]"
    else:
        rng = "[.,.]"
    _line(s, "range", rng, "units", _units(valid, 1e-6 if v.vtype == "float" else 1e-12))
    uniq = np.unique(valid)
    _line(s, "unique values", str(len(uniq)), "missing .", f"{N - len(valid)}/{N}")
    # VERIFICAR: com missing estendidos (.a, .b...) o Stata mostra "missing .*:"
    out.write("\n", "text")
    labs = ds.value_labels.get(v.value_label, {}) if v.value_label else {}
    if len(uniq) == 0:
        return
    if len(uniq) <= tab_max:
        out.write(f"{'tabulation:':>23}  ", "text")
        if labs:
            out.write("Freq.   Numeric  Label\n", "text")
        else:
            out.write("Freq.  Value\n", "text")
        for u in uniq:
            n_u = int((valid == u).sum())
            out.write(" " * 25 + f"{n_u:>5}", "result")
            if labs:
                lab = labs.get(int(u), "") if float(u).is_integer() else ""
                out.write(f"{_num(u):>10}  {lab}".rstrip() + "\n", "result")
            else:
                out.write(f"  {_num(u)}\n", "result")
        nm = N - len(valid)
        if nm:
            out.write(" " * 25 + f"{nm:>5}", "result")
            out.write(f"{'.':>10}\n" if labs else "  .\n", "result")
        out.write("\n", "text")
        return
    if labs:
        # VERIFICAR: variável rotulada com muitos valores mostra exemplos?
        pass
    mean = float(valid.mean())
    sd = float(valid.std(ddof=1)) if len(valid) > 1 else SYS
    _line(s, "mean", f"{_g(mean):>8}")
    _line(s, "std. dev", f"{_g(sd):>8}")
    out.write("\n", "text")
    srt = np.sort(valid)
    ps = [S.percentile(srt, q) for q in (10, 25, 50, 75, 90)]
    # primeira coluna com 9 caracteres, as demais com 10 (compat 0803)
    out.write(f"{'percentiles:':>23}  ", "text")
    out.write("".join(f"{str(q) + '%':>{9 if i == 0 else 10}}" for i, q in enumerate((10, 25, 50, 75, 90)))
              + "\n", "text")
    out.write(" " * 25, "text")
    out.write("".join(f"{_num(p):>{9 if i == 0 else 10}}" for i, p in enumerate(ps)) + "\n", "result")
    out.write("\n", "text")


def _compact(s: "Session", names: list[str], mask: np.ndarray) -> None:
    """Uma linha por variável. Larguras pelo conteúdo, como no Stata 14
    (compat 0803: "Variable   Obs Unique      Mean   Min    Max  Label")."""
    ds = s.data
    out = s.output
    rows = []
    for n in names:
        v = ds.get(n)
        if v.is_string:
            col = np.asarray(v.raw)[mask]
            nonmiss = col[col != ""]
            rows.append((n, len(nonmiss), len(set(nonmiss.tolist())), ".", ".", ".", v.label))
            continue
        x = v.data[mask]
        valid = x[x < SYS]
        if len(valid):
            rows.append((n, len(valid), len(np.unique(valid)), _num(float(valid.mean())),
                         _num(valid.min()), _num(valid.max()), v.label))
        else:
            rows.append((n, 0, 0, ".", ".", ".", v.label))
    wn = max([8] + [len(r[0]) for r in rows])
    wm = max([4] + [len(r[3]) for r in rows]) + 2
    wmin = max([3] + [len(r[4]) for r in rows]) + 2
    wmax = max([3] + [len(r[5]) for r in rows]) + 2
    out.write("\n" + f"{'Variable':<{wn}}{'Obs':>6}{'Unique':>7}{'Mean':>{wm}}{'Min':>{wmin}}{'Max':>{wmax}}  Label\n",
              "text")
    out.write("-" * _width(s) + "\n", "text")
    for n, obs, uq, mean, mn, mx, lab in rows:
        out.write(f"{n:<{wn}}{obs:>6}{uq:>7}{mean:>{wm}}{mn:>{wmin}}{mx:>{wmax}}", "result")
        out.write(f"  {lab}".rstrip() + "\n", "text")
    out.write("-" * _width(s) + "\n", "text")


@command("codebook", "codebook")
def cmd_codebook(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"all": 1, "header": 1, "notes": 1, "mv": 2, "tabulate": 3,
                                  "problems": 1, "detail": 1, "compact": 1, "dots": 4, "languages": 3})
    ds = s.data
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else list(ds.names)
    if not names:
        return
    mask = touse(s, p)
    if o.get("compact"):
        _compact(s, names, mask)
        return
    if o.get("problems"):
        _problems(s, names, mask)
        return
    tab_max = int(o["tabulate"]) if isinstance(o.get("tabulate"), str) else 9
    s.output.write("\n", "text")
    if o.get("header") or o.get("all"):
        out = s.output
        out.write(f"{'Dataset:':>17}  {ds.filename or '(unnamed)'}\n", "text")   # VERIFICAR
        out.write(f"{'Last saved:':>17}  {ds.timestamp or 'never'}\n", "text")
        out.write(f"{'Label:':>17}  {ds.label or '(none)'}\n", "text")
        out.write(f"{'Number of vars:':>17}  {ds.nvars:,}\n", "text")
        out.write(f"{'Number of obs:':>17}  {ds.nobs:,}\n", "text")
        out.write(f"{'Size:':>17}  {ds.width() * ds.nobs:,} bytes\n\n", "text")
    for n in names:
        _codebook_var(s, ds.get(n), mask, tab_max, o)


def _problems(s: "Session", names: list[str], mask: np.ndarray) -> None:
    """codebook, problems (compat 0803)."""
    ds = s.data
    out = s.output
    probs: dict[str, list[str]] = {}

    def add(kind, n):
        probs.setdefault(kind, []).append(n)

    for n in names:
        v = ds.get(n)
        if v.is_string:
            col = np.asarray(v.raw)[mask]
            vals = set(col.tolist())
            if len(vals) <= 1:
                add("constant (or all missing) vars", n)
            longest = max((len(x) for x in col), default=0)
            decl = int(v.vtype[3:]) if v.vtype[3:].isdigit() else 0
            if decl and longest < decl:
                add("str# vars that may be compressed", n)
            if any(x != x.lstrip() for x in vals if x):
                add("string vars with leading blanks", n)
            if any(x != x.rstrip() for x in vals if x):
                add("string vars with trailing blanks", n)
            if any(" " in x.strip() for x in vals):
                add("string vars with embedded blanks", n)
        else:
            x = v.data[mask]
            valid = x[x < SYS]
            if len(valid) == 0 or (np.all(valid == valid[0]) and len(valid) == len(x)):
                add("constant (or all missing) vars", n)
            if v.value_label and v.value_label not in ds.value_labels:
                add("vars with nonexistent labels", n)
    out.write("\n", "text")
    if not probs:
        out.write("no potential problems in dataset\n", "text")   # VERIFICAR
        return
    fname = ds.filename or "[unnamed]"
    out.write(f"   Potential problems in dataset   {fname}\n\n", "text")
    out.write(f"{'potential problem':>32}   variables\n", "text")
    out.write("-" * 50 + "\n", "text")
    order = ["constant (or all missing) vars", "vars with nonexistent labels", "incompletely labeled vars",
             "str# vars that may be compressed", "string vars with leading blanks",
             "string vars with trailing blanks", "string vars with embedded blanks"]
    for kind in order:
        if kind in probs:
            out.write(f"{kind:>34}   ", "text")
            out.write(" ".join(probs[kind]) + "\n", "result")
    out.write("-" * 50 + "\n", "text")


# ---------------------------------------------------------------------------
# misstable
# ---------------------------------------------------------------------------

@command("misstable", "misstable")
def cmd_misstable(s: "Session", args: str) -> None:
    t = args.strip()
    sub, _, rest = t.partition(" ")
    if "," in sub:
        sub, _, tail = sub.partition(",")
        rest = "," + tail + " " + rest
    full = None
    for name, m in (("summarize", 3), ("patterns", 3), ("tree", 4), ("nested", 4)):
        if len(sub) >= m and name.startswith(sub):
            full = name
    if full is None:
        raise StataError(198, "invalid syntax")
    p = parse_standard(rest)
    o = match_options(p.options, {"all": 3, "showzeros": 4, "generate": 3, "exok": 4, "frequency": 4,
                                  "asis": 4, "bypatterns": 3, "replace": 7, "clear": 5})
    ds = s.data
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else list(ds.names)
    if full != "summarize":
        names = [n for n in names if not ds.get(n).is_string]
    mask = touse(s, p)
    if full == "summarize":
        _miss_summarize(s, names, mask, o)
    elif full == "patterns":
        _miss_patterns(s, names, mask, o)
    else:
        raise StataError(198, f"misstable {full} not implemented")   # VERIFICAR


def _miss_summarize(s: "Session", names: list[str], mask: np.ndarray, o: dict) -> None:
    ds = s.data
    out = s.output
    rows = []
    for n in names:
        if ds.get(n).is_string:
            if o.get("all"):
                rows.append((n, None, None, None, None, None, None))
            continue
        x = ds.get(n).data[mask]
        eq = int((x == SYS).sum())
        gt = int((x > SYS).sum())
        lt = x[x < SYS]
        if eq + gt == 0 and not o.get("all"):
            continue
        rows.append((n, eq, gt, len(lt), len(np.unique(lt)), lt.min() if len(lt) else SYS,
                     lt.max() if len(lt) else SYS))
    if not rows:
        out.write("(variables nonmissing or string)\n", "text")
        return
    out.write("\n" + " " * 63 + "Obs<.\n", "text")
    out.write(" " * 48 + "+" + "-" * 30 + "\n", "text")
    out.write(f"{'':>15}|{'':>32}| Unique\n", "text")
    out.write(f"{'Variable':>14} |{'Obs=.':>10}{'Obs>.':>10}{'Obs<.':>10}  | values        Min         Max\n", "text")
    out.write("  " + "-" * 13 + "+" + "-" * 32 + "+" + "-" * 30 + "\n", "text")
    for n, eq, gt, lt, uq, mn, mx in rows:
        out.write(f"{_abbrev(n, 12):>14} |", "text")
        if eq is None:
            out.write(f"{' (string variable)':<32}|\n", "text")
            continue
        uq_t = f"{uq:>7}" if uq <= 500 else "   >500"
        out.write(f"{eq if eq or o.get('showzeros') else '':>10}{gt if gt or o.get('showzeros') else '':>10}{lt:>10}"
                  f"  |{uq_t}  {_num(mn):>9}   {_num(mx):>9}\n", "result")
    out.write("  " + "-" * 77 + "\n", "text")
    num = [r for r in rows if r[1] is not None]
    if num:
        n_, eq, gt, lt, uq, mn, mx = num[-1]
        # ordem do return list (compat 0803); valores da última variável numérica
        s.r = {"vartype": "numeric", "N_lt_dot": float(lt), "min": float(mn), "max": float(mx),
               "N_eq_dot": float(eq), "N_gt_dot": float(gt), "K_uniq": float(uq)}


def _miss_patterns(s: "Session", names: list[str], mask: np.ndarray, o: dict) -> None:
    ds = s.data
    out = s.output
    cols = [n for n in names if np.any(ds.get(n).data[mask] >= SYS)]
    if not cols:
        out.write("(no missing values)\n", "text")   # VERIFICAR
        return
    # ordem das variáveis: menos missing primeiro (VERIFICAR)
    cols.sort(key=lambda n: int((ds.get(n).data[mask] >= SYS).sum()))
    M_ = np.column_stack([(ds.get(n).data[mask] < SYS).astype(int) for n in cols])
    pats, counts = np.unique(M_, axis=0, return_counts=True)
    # padrão completo primeiro, depois por frequência (VERIFICAR)
    order = sorted(range(len(pats)), key=lambda j: (0 if pats[j].sum() == len(cols) else 1, -counts[j]))
    N = int(mask.sum())
    freq = bool(o.get("frequency"))
    k = len(cols)
    # layout do Stata 14 (compat 0803)
    out.write("\n" + " " * 3 + "Missing-value patterns\n", "text")
    out.write(" " * 5 + "(1 means complete)\n\n", "text")
    out.write(f"{'':>14}|   Pattern\n", "text")
    head = "    Frequency " if freq else "    Percent   "
    out.write(head + "|  " + "  ".join(str(j + 1) for j in range(k)) + "\n", "text")
    dash = max(13, 3 * k + 7)   # VERIFICAR com mais variáveis
    out.write("  " + "-" * 12 + "+" + "-" * dash + "\n", "text")
    first = True
    for idx in order:
        c = int(counts[idx])
        if freq:
            val = f"{c:>13,} "
        else:
            val = f"{100 * c / N:>9.0f}" + ("%    " if first else "     ")
        out.write(f"{val}|", "text")
        out.write("  " + "  ".join(str(x) for x in pats[idx]) + "\n", "result")
        if first:
            out.write(f"{'':>14}|\n", "text")
            first = False
    out.write("  " + "-" * 12 + "+" + "-" * dash + "\n", "text")
    out.write((f"{N:>13,} " if freq else f"{'100':>9}%    ") + "|\n\n", "text")
    out.write("  Variables are  " + "  ".join(f"({j + 1}) {n}" for j, n in enumerate(cols)) + "\n", "text")
    s.r = {"N_complete": float(int((M_.sum(axis=1) == k).sum())), "N_incomplete": float(int((M_.sum(axis=1) < k).sum())),
           "K": float(len(pats))}   # VERIFICAR nomes
