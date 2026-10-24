"""Cluster analysis ([MV] cluster, [MV] cluster kmeans and kmedians, [MV] cluster linkage,
[MV] cluster generate, [MV] cluster dendrogram, [MV] cluster stop, [MV] cluster utility,
[MV] matrix dissimilarity).

The computation is in core/cluster.py. A cluster analysis is kept in a dataset characteristic
(JSON) plus the variables Stata creates; layouts that the manual does not show carry a VERIFY
mark and are settled by the reference log of compat case 1906."""

from __future__ import annotations

import json
import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import cluster as K
from ..core import missing as M
from ..core.dataset import Variable
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lang.words import parse_numlist
from .mvcommon import numeric_names, sample
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_LINKAGES = {"singlelinkage": "single", "averagelinkage": "average", "completelinkage": "complete",
             "waveragelinkage": "waverage", "medianlinkage": "median", "centroidlinkage": "centroid",
             "wardslinkage": "ward"}
_SUBS = [("kmeans", 2), ("kmedians", 2), ("singlelinkage", 1), ("averagelinkage", 1), ("completelinkage", 1),
         ("waveragelinkage", 2), ("medianlinkage", 2), ("centroidlinkage", 2), ("wardslinkage", 1),
         ("generate", 1), ("dendrogram", 2), ("list", 1), ("dir", 2), ("drop", 2), ("rename", 3),
         ("stop", 2), ("notes", 1)]
_MEASURE_NAMES = {"l2": "L2", "euclidean": "L2", "l2squared": "L2squared", "l1": "L1", "absolute": "L1",
                  "cityblock": "L1", "manhattan": "L1", "linf": "Linf", "maximum": "Linf", "chebyshev": "Linf",
                  "canberra": "Canberra", "correlation": "correlation", "angular": "angular"}


def _measure(text, default: str) -> tuple[str, float | None]:
    if not isinstance(text, str) or not text.strip():
        return default, None
    t = text.strip()
    low = t.lower()
    if low.startswith("lpower(") or (low.startswith("l(") and low.endswith(")")):
        return "Lp", float(t[t.index("(") + 1:-1])
    if low in _MEASURE_NAMES:
        return _MEASURE_NAMES[low], None
    binary = {m_.lower(): m_ for m_ in K.BINARY}
    if low in binary:
        return binary[low], None
    raise StataError(198, f"invalid measure {t}")


def _store(ds) -> dict:
    return ds.chars.setdefault("_dta", {})


def _info(ds, name: str) -> dict:
    raw = ds.chars.get("_dta", {}).get(f"_clinfo_{name}")
    if raw is None:
        raise StataError(111, f"cluster {name} not found")                                 # VERIFY
    return json.loads(raw)


def _names(ds) -> list[str]:
    return [k[8:] for k in ds.chars.get("_dta", {}) if k.startswith("_clinfo_")]


def _default_name(ds) -> str:
    used = set(_names(ds))
    i = 1
    while f"_clus_{i}" in used:
        i += 1
    return f"_clus_{i}"


def _put_var(ds, name: str, values: np.ndarray, label: str = "") -> None:
    if ds.has(name):
        ds.drop_vars([name])
    var = Variable(name, "double" if name.endswith(("_hgt", "_pht")) else "float", np.full(ds.nobs, SYS))
    var.data = values
    ds.add(var)


def _check_new(ds, names: list[str]) -> None:
    for n in names:
        if ds.has(n):
            raise StataError(110, f"{n} already defined")


@command("cluster", "cluster")
def cmd_cluster(s: "Session", args: str) -> None:
    words = args.strip().split(None, 1)
    if not words:
        raise StataError(198, "invalid syntax")
    sub = words[0].rstrip(",")
    rest = words[1] if len(words) > 1 else ""
    full = next((n for n, m in _SUBS if len(sub) >= m and n.startswith(sub)), None)
    if full is None:
        raise StataError(199, f"unrecognized cluster subcommand {sub}")                  # VERIFY
    if full in _LINKAGES:
        _linkage(s, full, rest)
    elif full in ("kmeans", "kmedians"):
        _kmeans(s, full, rest)
    else:
        {"generate": _generate, "dendrogram": _dendrogram, "list": _list, "dir": _dir, "drop": _drop,
         "rename": _rename, "stop": _stop, "notes": _notes}[full](s, rest)


# ---------------------------------------------------------------------------
# linkage
# ---------------------------------------------------------------------------

def _linkage(s: "Session", sub: str, rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"measure": 3, "name": 2, "generate": 3, "force": 5}) if p.options.strip() else {}
    names = numeric_names(s, p.varlist, 1)
    X, w, wtype, mask = sample(s, p, names, allowed=())
    ds = s.data
    if X.shape[0] < 2:
        raise StataError(2001, "insufficient observations")
    meas, power = _measure(o.get("measure"), "L2squared" if sub in ("medianlinkage", "centroidlinkage") else "L2")
    name = str(o["name"]).strip() if isinstance(o.get("name"), str) else _default_name(ds)
    if name in _names(ds):
        raise StataError(110, f"cluster name {name} already exists")                      # VERIFY
    gen = str(o["generate"]).strip() if isinstance(o.get("generate"), str) else name
    D = K.pairwise(X, None, meas, power)
    _hierarchical(s, D, np.flatnonzero(mask), _LINKAGES[sub], name, gen,
                  meas if power is None else f"L({power:g})", names, f"cluster {sub} {rest.strip()}")


def _hierarchical(s: "Session", D: np.ndarray, rows: np.ndarray, method: str, name: str, gen: str,
                  dissim: str, names: list[str], cmd: str) -> None:
    """Runs the linkage on the dissimilarities D of the observations `rows` and stores the
    cluster variables and the description of the analysis."""
    ds = s.data
    pseudo = method in ("median", "centroid")
    _check_new(ds, [f"{gen}_id", f"{gen}_ord", f"{gen}_hgt"] + ([f"{gen}_pht"] if pseudo else []))
    merges = K.linkage(D, method)
    n = D.shape[0]
    order = K.leaf_order(merges, n)
    pos = np.empty(n, dtype=int)
    pos[order] = np.arange(1, n + 1)
    first_h = np.full(n, np.nan)
    kids = {n + r: (a, b) for r, (a, b, _) in enumerate(merges)}
    # pseudo-heights make the heights of median and centroid linkage monotone (VERIFY)
    pseudo_h = {}
    for r, (a, b, h) in enumerate(merges):
        pseudo_h[n + r] = max([h, pseudo_h.get(a, h), pseudo_h.get(b, h)])

    def leaves(c):
        return [c] if c < n else leaves(kids[c][0]) + leaves(kids[c][1])

    for r, (a, b, h) in enumerate(merges):
        for lf in leaves(a) + leaves(b):
            if np.isnan(first_h[lf]):
                first_h[lf] = h
    ids = np.full(ds.nobs, SYS)
    ords = np.full(ds.nobs, SYS)
    hgts = np.full(ds.nobs, SYS)
    ids[rows] = np.arange(1, n + 1)
    ords[rows] = pos
    hgts[rows] = first_h
    _put_var(ds, f"{gen}_id", ids)
    _put_var(ds, f"{gen}_ord", ords)
    _put_var(ds, f"{gen}_hgt", hgts)
    if pseudo:
        first_p = np.full(n, np.nan)
        for r, (a, b, h) in enumerate(merges):
            for lf in leaves(a) + leaves(b):
                if np.isnan(first_p[lf]):
                    first_p[lf] = pseudo_h[n + r]
        phts = np.full(ds.nobs, SYS)
        phts[rows] = first_p
        _put_var(ds, f"{gen}_pht", phts)
    info = {"type": "hierarchical", "method": method, "dissim": dissim,
            "vars": names, "gen": gen, "pseudo": pseudo, "rows": rows.tolist(), "merges": [[a, b, h] for a, b, h in merges],
            "cmd": cmd, "n": n}
    _store(ds)[f"_clinfo_{name}"] = json.dumps(info)
    s.notify_state()


_SHORT = {"singlelinkage": "single", "averagelinkage": "average", "completelinkage": "complete",
          "waveragelinkage": "waverage", "medianlinkage": "median", "centroidlinkage": "centroid",
          "wardslinkage": "wards"}


@command("clustermat", "clustermat")
def cmd_clustermat(s: "Session", args: str) -> None:
    """clustermat linkage matname [if] [in], shape() add clear labelvar() name() generate()
    ([MV] clustermat): hierarchical clustering of a dissimilarity matrix."""
    from .matrix import get_matrix
    words = args.strip().split(None, 2)
    if len(words) < 2:
        raise StataError(198, "invalid syntax")
    sub = words[0]
    full = next((n for n, m in _SUBS if n in _LINKAGES and (len(sub) >= m and n.startswith(sub) or _SHORT[n] == sub)), None)
    if full is None:
        raise StataError(199, f"unrecognized clustermat subcommand {sub}")              # VERIFY
    rest = words[1] + (" " + words[2] if len(words) > 2 else "")
    p = parse_standard(rest)
    o = match_options(p.options, {"shape": 2, "add": 3, "clear": 5, "labelvar": 5, "name": 2, "generate": 3,
                                  "force": 5}) if p.options.strip() else {}
    mat = get_matrix(s, p.varlist.strip())
    shape = str(o.get("shape", "full")).strip().lower()
    M_ = np.array(mat.data, dtype=float)
    if shape == "full":
        if M_.shape[0] != M_.shape[1]:
            raise StataError(503, "conformability error")
        if not o.get("force") and not np.allclose(M_, M_.T):
            raise StataError(505, f"matrix {p.varlist.strip()} not symmetric")         # VERIFY
        D = (M_ + M_.T) / 2
        labels = list(mat.rownames)
    elif shape in ("lower", "llower", "upper", "uupper"):
        v = M_.ravel()
        diag = shape in ("lower", "upper")
        n = int(round((math.sqrt(8 * v.size + 1) + (-1 if diag else 1)) / 2))
        D = np.zeros((n, n))
        idx = np.tril_indices(n, 0 if diag else -1) if shape.endswith("lower") else np.triu_indices(n, 0 if diag else 1)
        if len(idx[0]) != v.size:
            raise StataError(503, "conformability error")                               # VERIFY
        D[idx] = v
        D = D + D.T - np.diag(np.diag(D))
        labels = [str(i + 1) for i in range(n)]
    else:
        raise StataError(198, f"shape({shape}) invalid")
    np.fill_diagonal(D, 0.0)
    n = D.shape[0]
    ds = s.data
    if o.get("clear") or ds.nobs == 0 and not ds.nvars:
        if o.get("clear"):
            ds.clear()
        ds.set_obs(n)
        mask = np.ones(n, dtype=bool)
    elif o.get("add"):
        from ._util import touse
        mask = touse(s, p)
        if int(mask.sum()) != n:
            raise StataError(198, f"matrix has {n} rows but there are {int(mask.sum())} observations")   # VERIFY
    else:
        raise StataError(4, "no; data in memory would be lost")                         # VERIFY
    if isinstance(o.get("labelvar"), str):
        lv = o["labelvar"].strip()
        _check_new(ds, [lv])
        vals = np.array([""] * ds.nobs, dtype=object)
        vals[np.flatnonzero(mask)] = labels
        ds.add(Variable(lv, f"str{max(len(x) for x in labels)}", vals))
    name = str(o["name"]).strip() if isinstance(o.get("name"), str) else _default_name(ds)
    if name in _names(ds):
        raise StataError(110, f"cluster name {name} already exists")                      # VERIFY
    gen = str(o["generate"]).strip() if isinstance(o.get("generate"), str) else name
    _hierarchical(s, D, np.flatnonzero(mask), _LINKAGES[full], name, gen, f"user matrix {p.varlist.strip()}",
                  [], f"clustermat {args.strip()}")


# ---------------------------------------------------------------------------
# kmeans, kmedians
# ---------------------------------------------------------------------------

def _kmeans(s: "Session", sub: str, rest: str) -> None:
    from ..core import rng
    p = parse_standard(rest)
    o = match_options(p.options, {"k": 1, "name": 2, "measure": 3, "start": 3, "iterate": 4, "generate": 3,
                                  "keepcenters": 4, "force": 5}) if p.options.strip() else {}
    median = sub == "kmedians"
    names = numeric_names(s, p.varlist, 1)
    if not isinstance(o.get("k"), str):
        raise StataError(198, "option k() required")
    k = int(float(o["k"]))
    X, w, wtype, mask = sample(s, p, names, allowed=())
    n = X.shape[0]
    if k < 1 or k > n:
        raise StataError(198, "k() must be between 1 and the number of observations")    # VERIFY
    ds = s.data
    meas, power = _measure(o.get("measure"), "L2")
    name = str(o["name"]).strip() if isinstance(o.get("name"), str) else _default_name(ds)
    if name in _names(ds):
        raise StataError(110, f"cluster name {name} already exists")                      # VERIFY
    gen = str(o["generate"]).strip() if isinstance(o.get("generate"), str) else name
    _check_new(ds, [gen])
    start = str(o.get("start", "krandom")).strip()
    sw = start.split("(")[0].lower()
    arg = start[start.index("(") + 1:start.rindex(")")] if "(" in start else ""
    if sw == "firstk":
        centers = X[:k].copy()
    elif sw == "lastk":
        centers = X[n - k:].copy()
    elif sw in ("krandom", "random"):
        idx = np.argsort(rng.RNG.uniform(n), kind="stable")[:k]
        centers = X[np.sort(idx)].copy()
    else:
        if sw == "everykth":
            grp = np.arange(n) % k
        elif sw == "segments":
            grp = np.minimum((np.arange(n) * k) // n, k - 1)
        elif sw == "prandom":
            grp = np.minimum((rng.RNG.uniform(n) * k).astype(int), k - 1)
        elif sw == "group":
            gv = ds.get(s.expand_varlist(arg)[0]).data[mask]
            lv = sorted(set(gv.tolist()))
            if len(lv) != k:
                raise StataError(198, "group() must define exactly k groups")             # VERIFY
            grp = np.searchsorted(lv, gv)
        else:
            raise StataError(198, f"start({start}) is not allowed")
        centers = np.vstack([np.median(X[grp == g], axis=0) if median else X[grp == g].mean(0) for g in range(k)])
    maxit = int(float(o["iterate"])) if o.get("iterate") else 10000
    labels, centers, it = K.kmeans(X, centers, meas, median, maxit)
    rows = np.flatnonzero(mask)
    vals = np.full(ds.nobs, SYS)
    vals[rows] = labels + 1
    _put_var(ds, gen, vals)
    if o.get("keepcenters"):
        # the k final centers appended as new observations, in the clustered variables (VERIFY
        # the other variables, missing here, and the group variable of the new observations)
        n0 = ds.nobs
        ds.set_obs(n0 + k)
        for j, nm in enumerate(names):
            col = ds.get(nm).data.copy()
            col[n0:] = centers[:, j]
            ds.get(nm).data = col
    info = {"type": "partition", "method": sub, "dissim": meas if power is None else f"L({power:g})",
            "vars": names, "gen": gen, "rows": rows.tolist(), "k": k, "start": start, "iterations": it,
            "cmd": f"cluster {sub} {rest.strip()}", "n": n}
    _store(ds)[f"_clinfo_{name}"] = json.dumps(info)
    s.notify_state()


# ---------------------------------------------------------------------------
# generate
# ---------------------------------------------------------------------------

def _pick(s, o) -> str:
    ds = s.data
    names = _names(ds)
    if isinstance(o.get("name"), str):
        return str(o["name"]).strip()
    if not names:
        raise StataError(198, "no cluster analysis found")                                # VERIFY
    return names[-1]


def _generate(s: "Session", rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"name": 2, "ties": 2}) if p.options.strip() else {}
    if p.exp is None or "(" not in p.exp:
        raise StataError(198, "invalid syntax")
    ds = s.data
    stub = p.varlist.strip()
    what, _, arg = p.exp.partition("(")
    what = what.strip().lower()
    arg = arg.rstrip(") ").strip()
    info = _info(ds, _pick(s, o))
    if info["type"] != "hierarchical":
        raise StataError(198, "cluster generate requires a hierarchical cluster analysis")  # VERIFY
    n = info["n"]
    merges = [tuple(m) for m in info["merges"]]
    rows = np.array(info["rows"])
    if what in ("groups", "group"):
        ks = [int(x) for x in parse_numlist(arg)]
        labs = [K.cut_tree(merges, n, kk) for kk in ks]
    elif what == "cut":
        labs = [K.cut_height(merges, n, float(arg))]
        ks = [len(np.unique(labs[0]))]
    else:
        raise StataError(198, "invalid syntax")
    new = [f"{stub}{kk}" for kk in ks] if len(ks) > 1 else [stub]
    _check_new(ds, new)
    for nm, lab in zip(new, labs):
        vals = np.full(ds.nobs, SYS)
        vals[rows] = lab + 1
        _put_var(ds, nm, vals)
    s.notify_state()


# ---------------------------------------------------------------------------
# dendrogram
# ---------------------------------------------------------------------------

def _dendrogram(s: "Session", rest: str) -> None:
    from ..graphics.model import Graph, Layer, Panel
    from .graphcmd import finish
    c = rest.find(",")
    head, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    o = match_options(opt, {"cutnumber": 5, "cutvalue": 5, "labels": 4, "showcount": 5, "horizontal": 5, "name": 3,
                            "saving": 3, "nodraw": 3, "title": 3, "quick": 3}) if opt.strip() else {}
    ds = s.data
    name = head.strip() or _pick(s, {})
    info = _info(ds, name)
    if info["type"] != "hierarchical":
        raise StataError(198, "cluster dendrogram requires a hierarchical cluster analysis")  # VERIFY
    n = info["n"]
    if n > 100 and not (o.get("cutnumber") or o.get("cutvalue")):
        raise StataError(198, "too many leaves; consider using the cutvalue() or cutnumber() options")   # VERIFY rc
    merges = [tuple(m) for m in info["merges"]]
    order = K.leaf_order(merges, n)
    xpos = {lf: i + 1.0 for i, lf in enumerate(order)}
    ypos = {lf: 0.0 for lf in range(n)}
    xs, ys = [], []
    for r, (a, b, h) in enumerate(merges):
        xa, xb = xpos[a], xpos[b]
        ya, yb = ypos[a], ypos[b]
        xs += [xa, xa, xb, xb, np.nan]
        ys += [ya, h, h, yb, np.nan]
        xpos[n + r] = (xa + xb) / 2.0
        ypos[n + r] = h
    layers = [Layer("line", xs, ys, {}, {}, "")]
    ttl = f"Dendrogram for {name} cluster analysis"
    opts = {"_xtitle": "", "_ytitle": f"{info['dissim']} dissimilarity measure", "title": {"lines": [ttl]}, "legend": {"off": True}}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), {k: v for k, v in o.items() if k in ("name", "saving", "nodraw")})


# ---------------------------------------------------------------------------
# list, dir, drop, rename, notes
# ---------------------------------------------------------------------------

def _list(s: "Session", rest: str) -> None:
    c = rest.find(",")
    head, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    o = match_options(opt, {"type": 1, "state": 1, "dissimilarity": 1, "vars": 1, "all": 1}) if opt.strip() else {}
    ds = s.data
    names = head.split() or list(reversed(_names(ds)))
    out = s.output
    if not names:
        raise StataError(198, "no cluster analysis found")                                # VERIFY
    for nm in names:
        info = _info(ds, nm)
        method = "wards" if info["method"] == "ward" else info["method"]
        out.write(f"{nm}  (type: {info['type']},  method: {method},  dissimilarity: {info['dissim']})\n", "text")
        if info["type"] == "hierarchical":
            g = info["gen"]
            hgt = "real_height" if info.get("pseudo") else "height"
            out.write(f"      vars: {g}_id (id variable)\n            {g}_ord (order variable)\n"
                      f"            {g}_hgt ({hgt} variable)\n", "text")
            if info.get("pseudo"):
                out.write(f"            {g}_pht (pseudo_height variable)\n", "text")
        else:
            out.write(f"      vars: {info['gen']} (group variable)\n", "text")
        out.write(f"     other: cmd: {info['cmd']}\n            varlist: {' '.join(info['vars'])}\n", "text")
        if info["type"] == "partition":
            out.write(f"            k: {info['k']}\n            start: {info['start']}\n", "text")
        out.write("            range: 0 .\n\n", "text")


def _dir(s: "Session", rest: str) -> None:
    ds = s.data
    for nm in reversed(_names(ds)):
        s.output.write(f"{nm}\n", "text")


def _drop(s: "Session", rest: str) -> None:
    ds = s.data
    names = rest.split()
    if names == ["_all"] or not names:
        names = _names(ds)
    for nm in names:
        info = _info(ds, nm)
        g = info["gen"]
        drops = [f"{g}_id", f"{g}_ord", f"{g}_hgt", f"{g}_pht"] if info["type"] == "hierarchical" else [g]
        ds.drop_vars([d for d in drops if ds.has(d)])
        del _store(ds)[f"_clinfo_{nm}"]
    s.notify_state()


def _rename(s: "Session", rest: str) -> None:
    ds = s.data
    w = rest.split()
    if len(w) != 2:
        raise StataError(198, "invalid syntax")
    info = _info(ds, w[0])
    if w[1] in _names(ds):
        raise StataError(110, f"cluster name {w[1]} already exists")                      # VERIFY
    _store(ds)[f"_clinfo_{w[1]}"] = _store(ds).pop(f"_clinfo_{w[0]}")


def _notes(s: "Session", rest: str) -> None:
    raise StataError(198, "cluster notes is not available yet in OpenDTA")


# ---------------------------------------------------------------------------
# stop
# ---------------------------------------------------------------------------

def _stop(s: "Session", rest: str) -> None:
    c = rest.find(",")
    head, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    o = match_options(opt, {"rule": 2, "groups": 2, "matrix": 3, "variables": 3, "force": 5}) if opt.strip() else {}
    ds = s.data
    info = _info(ds, head.strip() or _pick(s, {}))
    rule = str(o.get("rule", "calinski")).strip().lower()
    if info["type"] != "hierarchical":
        raise StataError(198, "cluster stop requires a hierarchical cluster analysis")      # VERIFY
    names = s.expand_varlist(str(o["variables"])) if isinstance(o.get("variables"), str) else info["vars"]
    rows = np.array(info["rows"])
    X = np.column_stack([ds.get(v).data[rows] for v in names]).astype(float)
    n = info["n"]
    merges = [tuple(m) for m in info["merges"]]
    ks = [int(x) for x in parse_numlist(str(o["groups"]))] if isinstance(o.get("groups"), str) else list(range(1, 16))
    ks = [k for k in ks if 1 <= k <= n]
    out = s.output
    out.write("\n", "text")
    res = []

    def cell(text: str, left: int | None = None) -> str:
        l_ = (13 - len(text)) // 2 if left is None else left
        return " " * l_ + text + " " * (13 - len(text) - l_)

    if "calinski".startswith(rule) and len(rule) >= 1:
        out.write("+" + "-" * 27 + "+\n", "text")
        out.write(f"|{cell('')}|{cell('Calinski/')}|\n|{cell('Number of')}|{cell('Harabasz')}|\n", "text")
        out.write(f"|{cell('clusters')}|{cell('pseudo-F')}|\n|{'-' * 13}+{'-' * 13}|\n", "text")
        for k in ks:
            if k < 2:
                res.append([k, float("nan")])
                out.write(f"|{cell(str(k))}|{cell('.')}|\n", "text")
                continue
            f_ = K.calinski_harabasz(X, K.cut_tree(merges, n, k))
            res.append([k, f_])
            out.write(f"|{cell(str(k))}|{cell(f'{f_:.2f}')}|\n", "text")
        out.write("+" + "-" * 27 + "+\n", "text")
    elif "duda".startswith(rule) and len(rule) >= 1:
        out.write("+" + "-" * 41 + "+\n", "text")
        out.write(f"|{cell('')}|{cell('Duda/Hart', 9).rstrip() + ' ' * 0:<27}|\n", "text")
        out.write(f"|{cell('Number of')}|{cell('')}|{cell('pseudo', 2)}|\n", "text")
        out.write(f"|{cell('clusters')}|{cell('Je(2)/Je(1)')}|{cell('T-squared')}|\n", "text")
        out.write(f"|{'-' * 13}+{'-' * 13}+{'-' * 13}|\n", "text")
        # the row labelled k is the split of the k-th cluster, that is, the step from k to k + 1 clusters
        for k in ks:
            if k + 1 > n:
                continue
            ratio, t2 = K.duda_hart(X, merges, k + 1)
            res.append([k, ratio, t2])
            out.write(f"|{cell(str(k))}|{cell(f'{ratio:.4f}')}|{cell(f'{t2:.2f}')}|\n", "text")
        out.write("+" + "-" * 41 + "+\n", "text")
    else:
        raise StataError(198, f"rule({rule}) is not allowed")
    if isinstance(o.get("matrix"), str):
        from .matrix import Matrix, _store as mstore
        arr = np.array(res, dtype=float)
        cols = ["k", "calinski"] if arr.shape[1] == 2 else ["k", "Je2/Je1", "pseudoT2"]
        mstore(s)[str(o["matrix"]).strip()] = Matrix(arr, [f"r{i + 1}" for i in range(len(arr))], cols)
    s.r = {}


# ---------------------------------------------------------------------------
# matrix dissimilarity
# ---------------------------------------------------------------------------

def matrix_dissimilarity(s: "Session", rest: str) -> None:
    """matrix dissimilarity matname = [varlist] [if] [in] [, measure observations variables]"""
    from .matrix import Matrix, _store as mstore
    p = parse_standard(rest)
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    mat = p.varlist.strip()
    o = {}
    meas, power = "L2", None
    if p.options.strip():
        o = match_options(p.options.lower(), {"observations": 3, "variables": 3, "names": 5, "l2": 2, "euclidean": 3, "l2squared": 3,
                                      "l1": 2, "absolute": 3, "linf": 2, "maximum": 3, "canberra": 3,
                                      "correlation": 4, "angular": 3, "lpower": 2, "dissim": 4,
                                      **{m_.lower(): len(m_) for m_ in K.BINARY}})
        for key, val in (("l2", "L2"), ("euclidean", "L2"), ("l2squared", "L2squared"), ("l1", "L1"), ("absolute", "L1"),
                         ("linf", "Linf"), ("maximum", "Linf"), ("canberra", "Canberra"), ("correlation", "correlation"),
                         ("angular", "angular")):
            if o.get(key):
                meas = val
        for m_ in K.BINARY:
            if o.get(m_.lower()):
                meas = m_
        if isinstance(o.get("lpower"), str):
            meas, power = "Lp", float(o["lpower"])
    q = parse_standard(p.exp + (" if " + p.if_ if p.if_ else "") + (" in " + p.in_ if p.in_ else ""))
    names = numeric_names(s, q.varlist, 1) if q.varlist.strip() else [v.name for v in s.data.vars if not v.is_string]
    X, w, wtype, mask = sample(s, q, names, allowed=())
    similarity = meas in ("correlation", "angular") or meas in K.BINARY    # the matrix holds the similarity itself (VERIFY dissim())

    def dissimilarities(Z):
        D = K.pairwise(Z, None, meas, power)
        if similarity:
            D = 1.0 - D
            np.fill_diagonal(D, 1.0)
        return np.tril(D) + np.tril(D, -1).T        # exactly symmetric: BLAS rounding may break it

    if o.get("variables"):
        mstore(s)[mat] = Matrix(dissimilarities(X.T), names, names)
    else:
        D = dissimilarities(X)
        labs = [f"obs{i + 1}" for i in range(len(D))]                                   # VERIFY
        mstore(s)[mat] = Matrix(D, labs, labs)
