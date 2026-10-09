"""reshape ([D] reshape): troca entre os formatos largo e longo.

    reshape long stubnames, i(varlist) [j(varname [values]) string]
    reshape wide stubnames, i(varlist) j(varname [values]) [string]
    reshape long | reshape wide       (repete a última especificação)

Nos stubnames, `@` marca onde fica o valor de j (inc@r → inc80r). A
tabela impressa segue o Stata; posições das colunas e ordem das variáveis
no resultado estão marcadas VERIFICAR (caso compat/do/0305_reshape.do).
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import sorting
from ..core.dataset import Dataset, Variable, check_name, smallest_type_for, str_type_for, str_len
from ..core.errors import StataError
from ..core.grouping import group_ids
from ..lang.syntax import parse_options
from ..lang.words import num_str, parse_numlist, split_words
from .combine import common_type
from .registry import command

if TYPE_CHECKING:
    from ..session import Session


def _explain(s: "Session", head: str, paragraphs: list[str]) -> str:
    """Mensagem longa do reshape: parágrafos com recuo de 4, quebrados na
    largura de set linesize (como o Stata)."""
    import textwrap
    width = int(float(s.settings.get("linesize", 80)))
    out = [head]
    for par in paragraphs:
        if par.startswith("@"):            # bloco literal (figura, exemplo)
            out.append(par[1:])
        elif par == "":
            out.append("")
        else:
            out.extend(textwrap.wrap(par, width, initial_indent="    ", subsequent_indent="    ",
                                     break_on_hyphens=False) or ["    "])
    return "\n".join(out)


_PICTURE = """
         long                                wide
        +---------------+                   +------------------+
        | i   j   a   b |                   | i   a1 a2  b1 b2 |
        |---------------| <--- reshape ---> |------------------|
        | 1   1   1   2 |                   | 1   1   3   2  4 |
        | 1   2   3   4 |                   | 2   5   7   6  8 |
        | 2   1   5   6 |                   +------------------+
        | 2   2   7   8 |
        +---------------+"""


def _wide_name(stub: str, j: str) -> str:
    return stub.replace("@", j) if "@" in stub else stub + j


def _long_name(stub: str) -> str:
    return stub.replace("@", "")


def _stub_regex(stub: str, string: bool) -> re.Pattern[str]:
    part = "(.+)" if string else r"(-?\d+)"
    if "@" not in stub:
        stub = stub + "@"
    a, b = stub.split("@", 1)
    return re.compile("^" + re.escape(a) + part + re.escape(b) + "$")


def _jtext(x) -> str:
    return x if isinstance(x, str) else num_str(float(x))


def _table(s: "Session", frm: str, to: str, n0: int, n1: int, k0: int, k1: int,
           jname: str, nj: int, rows: list[tuple[str, str]]) -> None:
    w = s.output.write

    def num_line(label: str, a: int, b: int) -> None:
        w(f"{label}", "text")
        w(f"{a:>{39 - len(label)}}", "result")
        w("   ->", "text")
        w(f"{b:>8}\n", "result")

    def text_line(label: str, a: str, b: str) -> None:
        w(f"{label}{a:>{39 - len(label)}}   ->   {b}\n", "text")

    w("\n", "text")
    text_line("Data", frm, to)
    w("-" * 77 + "\n", "text")
    num_line("Number of obs.", n0, n1)
    num_line("Number of variables", k0, k1)
    jl = f"j variable ({nj} values)"
    if frm == "wide":
        text_line(jl, "", jname)
    else:
        text_line(jl, jname, "(dropped)")
    w("xij variables:\n", "text")
    for left, right in rows:
        text_line("", left, right)
    w("-" * 77 + "\n", "text")


def _parse(s: "Session", args: str):
    text = args.strip()
    comma = text.find(",")
    head, opts_text = (text[:comma], text[comma + 1:]) if comma != -1 else (text, "")
    parts = head.split(None, 1)
    if not parts:
        raise StataError(198, "invalid syntax")
    sub = parts[0].lower()
    stubs = split_words(parts[1]) if len(parts) > 1 else []
    opts = {n.lower(): a for n, a in parse_options(opts_text)} if opts_text.strip() else {}
    for k in opts:
        if k not in ("i", "j", "string", "atwl"):
            raise StataError(198, f"option {k} not allowed")
    return sub, stubs, opts


@command("reshape")
def cmd_reshape(s: "Session", args: str) -> None:
    sub, stubs, opts = _parse(s, args)
    if sub not in ("long", "wide"):
        # VERIFICAR: reshape error, reshape i/j/xij e reshape clear não implementados
        raise StataError(198, "invalid syntax")
    state = getattr(s, "_reshape", None)
    if not stubs:
        if not state:
            raise StataError(198, "reshape:  i(), j(), and stubnames not specified")   # VERIFICAR
        stubs, i_text, j_spec, string = state["stubs"], state["i"], state["j"], state["string"]
    else:
        if not opts.get("i"):
            raise StataError(198, "option i() required")
        i_text = opts["i"]
        j_spec = opts.get("j") or ("_j" if sub == "long" else None)
        if j_spec is None:
            raise StataError(198, "option j() required")
        string = "string" in opts
    j_words = split_words(j_spec)
    jname = j_words[0]
    j_values_spec = " ".join(j_words[1:])
    ivars = s.expand_varlist(i_text)
    s._reshape = {"stubs": stubs, "i": i_text, "j": j_spec, "string": string}
    if sub == "long":
        _to_long(s, stubs, ivars, jname, j_values_spec, string)
    else:
        _to_wide(s, stubs, ivars, jname, j_values_spec, string)
    s.notify_state()


def _to_long(s: "Session", stubs: list[str], ivars: list[str], jname: str,
             jvals_spec: str, string: bool) -> None:
    ds = s.data
    if ds.has(jname):
        raise StataError(110, f"variable {jname} already defined")
    check_name(jname)
    # valores de j e variáveis de cada stub
    found: dict[str, dict] = {}
    all_j: set = set()
    for stub in stubs:
        rx = _stub_regex(stub, string)
        found[stub] = {}
        for name in ds.names:
            m = rx.match(name)
            if m and name not in ivars:
                key = m.group(1) if string else int(m.group(1))
                found[stub][key] = name
        all_j.update(found[stub])
    if jvals_spec:
        jvals = split_words(jvals_spec) if string else [int(x) for x in parse_numlist(jvals_spec)]
    else:
        jvals = sorted(all_j)
    if not jvals:
        raise StataError(111, "no xij variables found")    # VERIFICAR
    if not string:
        bad = [x for x in jvals if x != int(x)]
        if bad:
            raise StataError(198, "j values must be integers")   # VERIFICAR
    # i deve identificar as observações
    s.output.write(f"(note: j = {' '.join(_jtext(x) for x in jvals)})\n", "text")
    for stub in stubs:
        for x in jvals:
            if x not in found[stub]:
                s.output.write(f"(note: {_wide_name(stub, _jtext(x))} not found)\n", "text")
    # i deve identificar as observações (texto observado no Stata 14)
    ids, k = group_ids(ds, ivars)
    if k != ds.nobs:
        il = " ".join(ivars)
        word = "variable" if len(ivars) == 1 else "variables"
        verb = "does" if len(ivars) == 1 else "do"
        raise StataError(9, _explain(s, f"{word} {il} {verb} not uniquely identify the observations", [
            f"Your data are currently wide.  You are performing a reshape long.  You specified i({il}) "
            f"and j({jname}).  In the current wide form, {word} {il} should uniquely identify the "
            "observations.  Remember this picture:",
            "@" + _PICTURE,
            "Type reshape error for a list of the problem observations."]))
    xij_names = {nm for st in stubs for x, nm in found[st].items() if x in set(jvals)}
    for stub in stubs:
        ln = _long_name(stub)
        if ds.has(ln) and ln not in xij_names:
            raise StataError(110, f"variable {ln} already defined")   # VERIFICAR
    n0, k0 = ds.nobs, ds.nvars
    nj = len(jvals)
    rows = np.repeat(np.arange(n0), nj)
    jidx = np.tile(np.arange(nj), n0)
    new = Dataset()
    new.label, new.notes, new.value_labels = ds.label, list(ds.notes), ds.value_labels
    new.nobs = n0 * nj
    for nm in ivars:
        v = ds.get(nm)
        nv = Variable(nm, v.vtype, np.empty(0), fmt=v.fmt, label=v.label, value_label=v.value_label)
        nv.raw = v.raw[rows]
        new.vars.append(nv)
    if string:
        jtype = str_type_for(max((str_len(x) for x in jvals), default=1) or 1)
        jv = Variable(jname, jtype, np.array([jvals[t] for t in jidx], dtype=object))
    else:
        jarr = np.array(jvals, dtype=np.float64)[jidx]
        # tipo mínimo, mas formato %9.0g (byte %9.0g; observado no Stata 14)
        jv = Variable(jname, smallest_type_for(np.array(jvals, dtype=np.float64)), jarr, fmt="%9.0g")
    long_vars: dict[str, Variable] = {}
    last_of: dict[str, str] = {}             # variável larga mais à direita de cada stub
    table_rows = []
    for stub in stubs:
        names = [found[stub].get(x) for x in jvals]
        present = [ds.get(nm) for nm in names if nm]
        if not present:
            raise StataError(111, f"no xij variables found for {stub}")   # VERIFICAR
        if len({v.is_string for v in present}) > 1:
            raise StataError(109, "type mismatch")   # VERIFICAR
        vtype = present[0].vtype
        for v in present[1:]:
            vtype = common_type(vtype, v.vtype)
        is_str = present[0].is_string
        cols = []
        for nm in names:
            src = ds.get(nm) if nm else None
            if is_str:
                cols.append(np.array(src.raw if src else [""] * n0, dtype=object))
            elif src is not None and src.vtype == vtype:
                cols.append(src.raw)
            else:
                cols.append(Variable("_t", vtype, src.data if src else np.full(n0, M.SYSMISS)).raw)
        nv = Variable(_long_name(stub), vtype, np.empty(0), fmt=present[0].fmt,
                      value_label=present[0].value_label)
        # linha a linha: (obs 1, j1), (obs 1, j2), ..., (obs 2, j1), ...
        nv.raw = np.stack(cols, axis=1).reshape(-1).copy()
        long_vars[stub] = nv
        last_of[max((nm for nm in names if nm), key=ds.index)] = stub
        table_rows.append((" ".join(nm for nm in names if nm), _long_name(stub)))
    new.vars.append(jv)
    # cada variável longa fica no lugar da última variável larga do seu stub;
    # as demais mantêm a ordem (observado no Stata 14)
    for v in ds.vars:
        if v.name in ivars:
            continue
        if v.name in xij_names:
            if v.name in last_of:
                new.vars.append(long_vars[last_of[v.name]])
            continue
        nv = Variable(v.name, v.vtype, np.empty(0), fmt=v.fmt, label=v.label, value_label=v.value_label)
        nv.raw = v.raw[rows]
        new.vars.append(nv)
    new.changed = True
    s.data = new
    sorting.sort(new, ivars + [jname])
    _table(s, "wide", "long", n0, new.nobs, k0, new.nvars, jname, nj, table_rows)


def _to_wide(s: "Session", stubs: list[str], ivars: list[str], jname: str,
             jvals_spec: str, string: bool) -> None:
    ds = s.data
    jvar = ds.get(jname)
    if jvar.is_string != string:
        if jvar.is_string:
            raise StataError(109, f"variable {jname} is string; specify string option")   # VERIFICAR
        raise StataError(109, "type mismatch")
    for stub in stubs:
        ds.get(_long_name(stub))
    longs = [_long_name(st) for st in stubs]
    n0, k0 = ds.nobs, ds.nvars
    if string:
        jcol = np.array(jvar.raw, dtype=object)
        jvals = split_words(jvals_spec) if jvals_spec else sorted(set(jcol))
    else:
        jcol = np.asarray(jvar.data, dtype=np.float64)
        if (jcol >= M.SYSMISS).any():
            raise StataError(498, f"variable {jname} contains missing values")   # VERIFICAR
        jvals = [float(x) for x in parse_numlist(jvals_spec)] if jvals_spec else sorted(set(jcol.tolist()))
    jtexts = [_jtext(x) for x in jvals]
    s.output.write(f"(note: j = {' '.join(jtexts)})\n", "text")
    ids, k = group_ids(ds, ivars)
    # (i, j) único
    pair_ids, kp = group_ids(ds, ivars + [jname])
    ilist = " ".join(ivars)
    if kp != n0:
        raise StataError(9, _explain(s, f"values of variable {jname} not unique within {ilist}", [
            f"Your data are currently long.  You are performing a reshape wide.  You specified "
            f"i({ilist}) and j({jname}).  There are observations within i({ilist}) with the same "
            f"value of j({jname}).  In the long data, variables i() and j() together must uniquely "
            "identify the observations.",
            "@" + _PICTURE,
            "Type reshape error for a list of the problem observations."]))   # VERIFICAR
    others = [v for v in ds.vars if v.name not in ivars and v.name != jname and v.name not in longs]
    from ..core.grouping import first_index
    first = first_index(ids, k)
    for v in others:
        col = v.raw
        if v.is_string:
            same = all(col[i] == col[first[ids[i]]] for i in range(n0))
        else:
            same = bool(np.all(col == col[first[ids]]))
        if not same:
            raise StataError(9, _explain(s, f"variable {v.name} not constant within {ilist}", [
                "Your data are currently long.  You are performing a reshape wide.  You typed "
                "something like",
                "",
                f"@        . reshape wide a b, i({ilist}) j({jname})",
                "",
                f"There are variables other than a, b, {ilist}, {jname} in your data.  They must be "
                f"constant within {ilist} because that is the only way they can fit into wide data "
                "without loss of information.",
                "",
                f"The variable or variables listed above are not constant within {ilist}.  Perhaps "
                "the values are in error.  Type reshape error for a list of the problem "
                "observations.",
                "",
                "Either that, or the values vary because they should vary, in which case you must "
                "either add the variables to the list of xij variables to be reshaped, or drop "
                "them."]))
    pos = {x: t for t, x in enumerate(jvals)}
    jpos = np.array([pos.get(x, -1) for x in (jcol.tolist())], dtype=np.int64)
    if (jpos < 0).any():
        keep = jpos >= 0   # valores de j fora da lista dada são descartados (VERIFICAR)
    else:
        keep = np.ones(n0, dtype=bool)
    new = Dataset()
    new.label, new.notes, new.value_labels = ds.label, list(ds.notes), ds.value_labels
    new.nobs = k
    for nm in ivars:
        v = ds.get(nm)
        nv = Variable(nm, v.vtype, np.empty(0), fmt=v.fmt, label=v.label, value_label=v.value_label)
        nv.raw = v.raw[first].copy()
        new.vars.append(nv)
    table_rows = []
    wide_vars: dict[tuple[int, str], Variable] = {}
    for stub, ln in zip(stubs, longs):
        src = ds.get(ln)
        names = []
        for t, jt in enumerate(jtexts):
            wn = _wide_name(stub, jt)
            check_name(wn)
            if wn in ivars or any(v.name == wn for v in others):
                raise StataError(110, f"variable {wn} already defined")
            if src.is_string:
                vals = np.array([""] * k, dtype=object)
            else:
                vals = Variable("_t", src.vtype, np.full(k, M.SYSMISS)).raw.copy()
            sel = keep & (jpos == t)
            vals[ids[sel]] = src.raw[sel]
            # rótulo "80 inc" (valor de j e nome da variável longa; Stata 14)
            nv = Variable(wn, src.vtype, np.empty(0), fmt=src.fmt, value_label=src.value_label,
                          label=f"{jt} {src.label or ln}")
            nv.raw = vals
            wide_vars[(t, stub)] = nv
            names.append(wn)
        table_rows.append((ln, " ".join(names)))
    # VERIFICAR: ordem intercalada por j (inc80 ue80 inc81 ue81 ...)
    for t in range(len(jtexts)):
        for stub in stubs:
            new.vars.append(wide_vars[(t, stub)])
    for v in others:
        nv = Variable(v.name, v.vtype, np.empty(0), fmt=v.fmt, label=v.label, value_label=v.value_label)
        nv.raw = v.raw[first].copy()
        new.vars.append(nv)
    new.changed = True
    s.data = new
    new.sortlist = list(ivars)
    _table(s, "long", "wide", n0, k, k0, new.nvars, jname, len(jvals), table_rows)
