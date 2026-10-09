"""Funções st_* do Mata: leitura e gravação de dados, macros, escalares e
matrizes do Stata ([M-5] st_data(), st_store(), st_local()...).

st_view() devolve uma cópia dos dados: alterar a matriz não altera a base
(no Stata, a view escreve nos dados; VERIFICAR/limitação conhecida).
"""

from __future__ import annotations

import re

import numpy as np

from ..core.dataset import Variable, check_name
from ..core.errors import StataError
from .library import lib
from .values import MV, MataError, SYS, empty, real, string


def _ds(eng):
    return eng.s.data


def _var_names(eng, jv: MV) -> list[str]:
    ds = _ds(eng)
    if jv.t == "string":
        text = " ".join(str(x) for x in jv.a.ravel())
        try:
            return eng.s.expand_varlist(text)
        except StataError:
            raise MataError(3500, "invalid Stata variable name")   # VERIFICAR
    idx = jv.a.ravel()
    if idx.size == 1 and idx[0] >= SYS:
        return list(ds.names)
    names = []
    for k in idx:
        if k >= SYS or not 1 <= k <= ds.nvars:
            raise MataError(3500, "invalid Stata variable name")
        names.append(ds.vars[int(k) - 1].name)
    return names


def _obs_rows(eng, iv: MV | None, sel: MV | None) -> np.ndarray:
    n = _ds(eng).nobs
    if iv is None or (iv.is_scalar and iv.a[0, 0] >= SYS):
        rows = np.arange(n)
    else:
        a = iv.a
        if iv.t != "real":
            raise MataError(3253, "nonreal found where real required")
        flat = a.ravel()
        if np.any(flat >= SYS) or np.any(flat < 1) or np.any(flat > n):
            raise MataError(3301, "subscript invalid")
        rows = flat.astype(np.int64) - 1
    if sel is not None and not (sel.t == "string" and sel.str_scalar() == "") \
            and not (sel.t == "real" and sel.is_scalar and sel.a[0, 0] >= SYS):
        if sel.t == "string":
            var = _ds(eng).get(eng.s.expand_varlist(sel.str_scalar())[0])
        else:
            var = _ds(eng).vars[int(sel.real_scalar()) - 1]
        keep = var.data[rows] != 0
        rows = rows[keep]
    return rows


@lib("st_nobs", 0)
def _st_nobs(eng, v, r):
    return real(float(_ds(eng).nobs))


@lib("st_nvar", 0)
def _st_nvar(eng, v, r):
    return real(float(_ds(eng).nvars))


@lib("st_varindex", 1, 2)
def _st_varindex(eng, v, r):
    ds = _ds(eng)
    out = []
    for nm in " ".join(str(x) for x in v[0].a.ravel()).split():
        try:
            full = eng.s.expand_varlist(nm)[0] if (len(v) > 1 and v[1].real_scalar()) else nm
            out.append(float(ds.index(full) + 1))
        except StataError:
            out.append(SYS)
    return real(np.array(out).reshape(1, -1))


@lib("st_varname", 1)
def _st_varname(eng, v, r):
    ds = _ds(eng)
    idx = v[0].a.ravel()
    arr = np.empty((1, idx.size), dtype=object)
    for k, i in enumerate(idx):
        if i >= SYS or not 1 <= i <= ds.nvars:
            raise MataError(3300, "argument out of range")
        arr[0, k] = ds.vars[int(i) - 1].name
    return MV(arr, "string")


def _data(eng, v, *, strings: bool) -> MV:
    rows = _obs_rows(eng, v[0], v[2] if len(v) > 2 else None)
    names = _var_names(eng, v[1])
    ds = _ds(eng)
    cols = []
    for nm in names:
        var = ds.get(nm)
        if var.is_string != strings:
            if strings:
                raise MataError(3254, "nonstring found where string required")   # VERIFICAR
            cols.append(np.full(len(rows), SYS))
            continue
        cols.append(np.asarray(var.raw if strings else var.data)[rows])
    if strings:
        out = np.empty((len(rows), len(names)), dtype=object)
        for j, c in enumerate(cols):
            out[:, j] = c
        return MV(out, "string")
    return real(np.column_stack(cols) if cols else np.zeros((len(rows), 0)))


@lib("st_data", 2, 3)
def _st_data(eng, v, r):
    return _data(eng, v, strings=False)


@lib("st_sdata", 2, 3)
def _st_sdata(eng, v, r):
    return _data(eng, v, strings=True)


@lib("st_view", 3, 4, outs=(0,))
def _st_view(eng, v, r):
    val = _data(eng, v[1:], strings=False)
    if r[0] is not None:
        r[0].set(val)
    return None


@lib("st_sview", 3, 4, outs=(0,))
def _st_sview(eng, v, r):
    val = _data(eng, v[1:], strings=True)
    if r[0] is not None:
        r[0].set(val)
    return None


def _store(eng, v, *, strings: bool) -> None:
    if len(v) == 4:
        iv, jv, sel, X = v
    else:
        iv, jv, X = v
        sel = None
    rows = _obs_rows(eng, iv, sel)
    names = _var_names(eng, jv)
    ds = _ds(eng)
    data = X.a
    if X.t != ("string" if strings else "real"):
        raise MataError(3250, "type mismatch")
    if data.shape != (len(rows), len(names)):
        if X.is_scalar:
            data = np.broadcast_to(data, (len(rows), len(names)))
        else:
            raise MataError(3200, "conformability error")
    for j, nm in enumerate(names):
        var = ds.get(nm)
        if strings:
            if not var.is_string:
                raise MataError(3250, "type mismatch")
            ds.set_string(var, np.array(data[:, j], dtype=object), rows)
        else:
            if var.is_string:
                raise MataError(3250, "type mismatch")
            ds.set_numeric(var, np.asarray(data[:, j], dtype=np.float64), rows, promote=False)
    eng.s.notify_state()


@lib("st_store", 3, 4)
def _st_store(eng, v, r):
    _store(eng, v, strings=False)
    return None


@lib("st_sstore", 3, 4)
def _st_sstore(eng, v, r):
    _store(eng, v, strings=True)
    return None


@lib("st_addvar", 2, 3)
def _st_addvar(eng, v, r):
    ds = _ds(eng)
    names = " ".join(str(x) for x in v[1].a.ravel()).split()
    types = [str(x) for x in v[0].a.ravel()] if v[0].t == "string" else \
        [f"str{int(x)}" for x in v[0].a.ravel()]
    if len(types) == 1:
        types = types * len(names)
    out = []
    for t, nm in zip(types, names):
        if ds.has(nm):
            raise MataError(3300, f"variable {nm} already defined")   # VERIFICAR
        try:
            check_name(nm)
        except StataError:
            raise MataError(3300, f"{nm} invalid name")
        if t.startswith("str"):
            ds.add(Variable(nm, t, np.array([""] * ds.nobs, dtype=object)))
        else:
            ds.add(Variable(nm, t, np.full(ds.nobs, SYS)))
        out.append(float(ds.index(nm) + 1))
    eng.s.notify_state()
    return real(np.array(out).reshape(1, -1))


@lib("st_addobs", 1, 2)
def _st_addobs(eng, v, r):
    ds = _ds(eng)
    ds.set_obs(ds.nobs + v[0].int_scalar())
    eng.s.notify_state()
    return None


@lib("st_dropvar", 1)
def _st_dropvar(eng, v, r):
    _ds(eng).drop_vars(_var_names(eng, v[0]))
    eng.s.notify_state()
    return None


@lib("st_keepvar", 1)
def _st_keepvar(eng, v, r):
    keep = set(_var_names(eng, v[0]))
    ds = _ds(eng)
    ds.drop_vars([n for n in ds.names if n not in keep])
    eng.s.notify_state()
    return None


# macros, escalares e matrizes --------------------------------------------

_RESULT = re.compile(r"^([res])\((\w+)\)$")


@lib("st_local", 1, 2)
def _st_local(eng, v, r):
    name = v[0].str_scalar()
    if len(v) == 1:
        return string(eng.s.macros.get_local(name))
    eng.s.macros.set_local(name, v[1].str_scalar())
    return None


@lib("st_global", 1, 2)
def _st_global(eng, v, r):
    name = v[0].str_scalar()
    m = _RESULT.match(name)
    if m:
        store = {"r": eng.s.r, "e": eng.s.e, "s": eng.s.sret}[m.group(1)]
        if len(v) == 1:
            val = store.get(m.group(2), "")
            return string(val if isinstance(val, str) else "")
        store[m.group(2)] = v[1].str_scalar()
        return None
    if len(v) == 1:
        return string(eng.s.macros.get_global(name))
    eng.s.macros.set_global(name, v[1].str_scalar())
    return None


@lib("st_numscalar", 1, 2)
def _st_numscalar(eng, v, r):
    name = v[0].str_scalar()
    m = _RESULT.match(name)
    store = {"r": eng.s.r, "e": eng.s.e, "s": eng.s.sret}[m.group(1)] if m else eng.s.scalars
    key = m.group(2) if m else name
    if len(v) == 1:
        val = store.get(key)
        if isinstance(val, (int, float)):
            return real(float(val))
        return empty("real", 0, 0)
    store[key] = v[1].real_scalar()
    return None


@lib("st_strscalar", 1, 2)
def _st_strscalar(eng, v, r):
    name = v[0].str_scalar()
    if len(v) == 1:
        val = eng.s.scalars.get(name)
        return string(val) if isinstance(val, str) else empty("string", 0, 0)
    eng.s.scalars[name] = v[1].str_scalar()
    return None


def _matrix_target(eng, name: str):
    from ..commands.matrix import _store
    m = _RESULT.match(name)
    if m:
        return {"r": eng.s.r, "e": eng.s.e, "s": eng.s.sret}[m.group(1)], m.group(2)
    return _store(eng.s), name


@lib("st_matrix", 1, 2)
def _st_matrix(eng, v, r):
    from ..commands.matrix import Matrix
    name = v[0].str_scalar()
    store, key = _matrix_target(eng, name)
    if len(v) == 1:
        mat = store.get(key)
        if not isinstance(mat, Matrix):
            return empty("real", 0, 0)
        return real(mat.data.copy())
    X = v[1]
    if X.t != "real":
        raise MataError(3253, "nonreal found where real required")
    if X.a.size == 0:
        store.pop(key, None)
        return None
    store[key] = Matrix(X.a.copy())
    return None


def _stripe(eng, v, which: str):
    from ..commands.matrix import Matrix
    name = v[0].str_scalar()
    store, key = _matrix_target(eng, name)
    mat = store.get(key)
    if not isinstance(mat, Matrix):
        raise MataError(3499, f"matrix {name} not found")   # VERIFICAR
    names = mat.rownames if which == "row" else mat.colnames
    if len(v) == 1:
        out = np.empty((len(names), 2), dtype=object)
        for k, nm in enumerate(names):
            eq, _, base = nm.rpartition(":")
            out[k, 0], out[k, 1] = eq, base
        return MV(out, "string")
    S = v[1]
    if S.t != "string" or S.cols != 2:
        raise MataError(3200, "conformability error")
    new = [(f"{a}:{b}" if a else str(b)) for a, b in S.a.tolist()]
    if which == "row":
        mat.rownames = new
    else:
        mat.colnames = new
    return None


@lib("st_matrixrowstripe", 1, 2)
def _st_rowstripe(eng, v, r):
    return _stripe(eng, v, "row")


@lib("st_matrixcolstripe", 1, 2)
def _st_colstripe(eng, v, r):
    return _stripe(eng, v, "col")


@lib("st_rclear", 0)
def _st_rclear(eng, v, r):
    eng.s.r = {}
    return None


@lib("st_eclear", 0)
def _st_eclear(eng, v, r):
    eng.s.e = {}
    return None


# atributos de variáveis ------------------------------------------------------

def _one_var(eng, v):
    names = _var_names(eng, v[0])
    if len(names) != 1:
        raise MataError(3200, "conformability error")
    return _ds(eng).get(names[0])


@lib("st_vartype", 1)
def _st_vartype(eng, v, r):
    return string(_one_var(eng, v).vtype)


@lib("st_varformat", 1, 2)
def _st_varformat(eng, v, r):
    var = _one_var(eng, v)
    if len(v) == 1:
        return string(var.fmt)
    var.fmt = v[1].str_scalar()
    return None


@lib("st_varlabel", 1, 2)
def _st_varlabel(eng, v, r):
    var = _one_var(eng, v)
    if len(v) == 1:
        return string(var.label)
    var.label = v[1].str_scalar()
    return None


@lib("st_varvaluelabel", 1, 2)
def _st_varvaluelabel(eng, v, r):
    var = _one_var(eng, v)
    if len(v) == 1:
        return string(var.value_label)
    var.value_label = v[1].str_scalar()
    return None


@lib("st_isnumvar", 1)
def _st_isnumvar(eng, v, r):
    return real(0.0 if _one_var(eng, v).is_string else 1.0)


@lib("st_isstrvar", 1)
def _st_isstrvar(eng, v, r):
    return real(1.0 if _one_var(eng, v).is_string else 0.0)


@lib("st_tempname", 0, 1)
def _st_tempname(eng, v, r):
    n = v[0].int_scalar() if v else 1
    from ..commands.preserve import _counter, current_scope
    scope = current_scope(eng.s)
    names = [_counter(eng.s) for _ in range(n)]
    scope.tempnames.extend(names)
    arr = np.empty((1, n), dtype=object)
    for k, nm in enumerate(names):
        arr[0, k] = nm
    return MV(arr, "string")


# executar comandos do Stata --------------------------------------------------

@lib("stata", 1, 3)
def _stata(eng, v, r):
    cmd = v[0].str_scalar()
    quiet = len(v) > 1 and v[1] is not None and v[1].real_scalar() != 0
    rc = _run_stata(eng, cmd, quiet)
    if rc:
        raise MataError(rc, "Stata returned error")    # VERIFICAR (3598?)
    return None


@lib("_stata", 1, 3)
def _stata_rc(eng, v, r):
    cmd = v[0].str_scalar()
    quiet = len(v) > 1 and v[1] is not None and v[1].real_scalar() != 0
    return real(float(_run_stata(eng, cmd, quiet)))


def _run_stata(eng, cmd: str, quiet: bool) -> int:
    s = eng.s
    out = s.output
    if quiet:
        out.quiet_depth += 1
    try:
        s.interp.execute(cmd)
        return 0
    except StataError as e:
        if not quiet:
            s.report_error(e)
        return e.rc
    finally:
        if quiet:
            out.quiet_depth -= 1
