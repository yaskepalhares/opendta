"""irt ([IRT] irt): item response theory models, estat report and irtgraph.

Every model is a gsem model with one latent variable Theta ~ N(0, 1) (variance fixed at 1) and
one response per item; the coefficients are stored in slope-intercept form, as gsem does, and
the table shows the discrimination a and the difficulties b by the delta method:

    1pl, 2pl   P(y = 1 | theta) = invlogit(a (theta - b)),  a common (1pl) or by item
    3pl        P(y = 1 | theta) = g + (1 - g) invlogit(a (theta - b)), g common to the items
               unless sepguessing
    grm        P(y >= k | theta) = invlogit(a (theta - b_k))
    nrm        P(y = k | theta) proportional to exp(a_k theta + c_k), a_0 = c_0 = 0
    pcm, gpcm  P(y = k | theta) proportional to exp(sum_{j<=k} a (theta - b_j)), a common (pcm)
    rsm        as pcm with b_j = b_i + kappa_j shared by the items, kappa_1 = 0

1pl, 2pl and grm are written as gsem paths; the other models are "custom" responses of gsem
(analytic derivatives with respect to theta, used by its adaptive quadrature). The names of the
rows of the nominal and rating-scale models and the layout of byparm are VERIFY (compat case
1909)."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from ..lazy import lazy
from .gsem import Eq, GModel, _g, _post, fit, lse
from .results import Estimates, crit

st = lazy("scipy.stats")
special = lazy("scipy.special")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_TITLES = {"1pl": "One-parameter logistic model", "2pl": "Two-parameter logistic model",
           "3pl": "Three-parameter logistic model", "grm": "Graded response model",
           "nrm": "Nominal response model", "pcm": "Partial credit model",
           "gpcm": "Generalized partial credit model", "rsm": "Rating scale model",
           "hybrid": "Hybrid model"}
_KINDS = ("1pl", "2pl", "3pl", "grm", "nrm", "pcm", "gpcm", "rsm")


# ---------------------------------------------------------------------------
# custom responses
# ---------------------------------------------------------------------------

class ThreePL:
    """P(1) = g + (1 - g) invlogit(a theta + c), g = invlogit(lg)."""

    def __init__(self, P, name, y, shared=True):
        self.y = (y != 0)
        p = min(max(float(np.mean(self.y)), .05), .95)
        self.ia = P.add(name, "Theta", 1.0)
        self.ic = P.add(name, "_cons", math.log(p / (1 - p)))
        # one pseudoguessing parameter for all the items unless sepguessing
        self.ig = P.add(name, "Guess", -2.0, "aux", label="_g_3pl" if shared else None)

    def logp(self, th, z):
        a, c, lg = th[self.ia], th[self.ic], th[self.ig]
        u = a * z + c
        y = self.y[:, None]
        # in logs, so that a guessing parameter near 0 or 1 does not overflow
        lgu, lgd = float(special.log_expit(lg)), float(special.log_expit(-lg))
        l1 = np.logaddexp(lgu, lgd + special.log_expit(u))
        l0 = lgd + special.log_expit(-u)
        return np.where(y, l1, l0)

    def d12(self, th, z):
        a, c, g = th[self.ia], th[self.ic], float(special.expit(th[self.ig]))
        u = a * z + c
        s = special.expit(u)
        s1 = s * (1 - s)
        s2 = s1 * (1 - 2 * s)
        P = g + (1 - g) * s
        d1 = (1 - g) * a * s1
        d2 = (1 - g) * a * a * s2
        l1 = d1 / P
        l2 = d2 / P - l1 * l1
        return (np.where(self.y, l1, -a * s), np.where(self.y, l2, -a * a * s1))


class Cat:
    """Categorical item: P(y = k) proportional to exp(A_k theta + C_k)."""

    def __init__(self, P, name, y, kind):
        self.kind = kind
        self.levels = np.unique(y)
        J = self.J = len(self.levels)
        if J < 3 and kind != "nrm" or J < 2:
            raise StataError(198, f"{name} has fewer than three categories")   # VERIFY
        self.yi = np.searchsorted(self.levels, y).astype(int)
        n = np.maximum(np.bincount(self.yi, minlength=J).astype(float), 1.0)
        eta0 = np.log(n / n[0])
        d0 = np.diff(eta0)
        if kind == "nrm":
            self.ia = [P.add(name, f"{k}.Theta", .5 * k) for k in range(1, J)]
            self.ic = [P.add(name, f"{k}._cons", eta0[k]) for k in range(1, J)]
        else:
            lab = {"pcm": "_a_pcm", "rsm": "_a_rsm"}.get(kind)
            self.ia = [P.add(name, "Theta", 1.0, "coef", lab)]
            if kind == "rsm":
                self.ib = P.add(name, "_cons", float(np.mean(d0)))
                self.ie = [P.add("kappa", f"kappa{j}", 0.0, "aux", f"_rsm_kappa{j}") for j in range(2, J)]
            else:
                self.ic = [P.add(name, f"{k}._cons", d0[k - 1]) for k in range(1, J)]

    def AC(self, th):
        J = self.J
        if self.kind == "nrm":
            return (np.concatenate([[0.0], th[self.ia]]), np.concatenate([[0.0], th[self.ic]]))
        A = th[self.ia[0]] * np.arange(J)
        if self.kind == "rsm":
            e = np.concatenate([[0.0, 0.0], th[self.ie]])
            C = th[self.ib] * np.arange(J) + np.cumsum(e)
        else:
            C = np.concatenate([[0.0], np.cumsum(th[self.ic])])
        return A, C

    def logp(self, th, z):
        A, C = self.AC(th)
        eta = z[:, :, None] * A + C
        norm = lse(eta, axis=2)
        pick = np.take_along_axis(eta, self.yi[:, None, None], axis=2)[..., 0]
        return pick - norm

    def d12(self, th, z):
        A, C = self.AC(th)
        eta = z[:, None] * A + C
        p = np.exp(eta - lse(eta, axis=1)[:, None])
        m1 = p @ A
        m2 = p @ (A * A)
        return A[self.yi] - m1, -(m2 - m1 * m1)


# ---------------------------------------------------------------------------
# items: parameters as reported and probabilities
# ---------------------------------------------------------------------------

class Item:
    def __init__(self, kind: str, name: str, r, levels=None):
        self.kind, self.name, self.r = kind, name, r
        self.levels = levels
        c = r.custom
        if kind in ("1pl", "2pl"):
            self.ia = r.cref[0][0]
            self.ic = r.xref[-1][0]
        elif kind == "grm":
            self.ia = r.cref[0][0]
            self.cuts = list(r.aux)
            self.levels = r.levels
        elif kind == "3pl":
            self.ia, self.ic, self.ig = c.ia, c.ic, c.ig
        else:
            self.c = c
            self.levels = c.levels

    def shared_a(self) -> bool:
        return self.kind in ("1pl", "pcm", "rsm")

    def rows(self, th):
        """[(kind, label, value)]: kind "row" has a value, "head" is a sub-header."""
        k = self.kind
        if k in ("1pl", "2pl"):
            a, c = th[self.ia], th[self.ic]
            return [("row", "Discrim", a), ("row", "Diff", -c / a)]
        if k == "3pl":
            a, c = th[self.ia], th[self.ic]
            return [("row", "Discrim", a), ("row", "Diff", -c / a),
                    ("row", "Guess", float(special.expit(th[self.ig])))]
        if k == "grm":
            a = th[self.ia]
            out = [("row", "Discrim", a), ("head", "Diff", 0.0)]
            for j, i in enumerate(self.cuts):
                out.append(("row", f">={self.levels[j + 1]:g}", th[i] / a))
            return out
        c = self.c
        lv = self.levels
        if k == "nrm":
            out = [("head", "Discrim", 0.0)]
            out += [("row", f"{lv[j]:g}", th[c.ia[j - 1]]) for j in range(1, len(lv))]
            out += [("head", "Diff", 0.0)]
            out += [("row", f"{lv[j]:g}", -th[c.ic[j - 1]] / th[c.ia[j - 1]]) for j in range(1, len(lv))]
            return out
        a = th[c.ia[0]]
        out = [("row", "Discrim", a), ("head", "Diff", 0.0)]
        if k == "rsm":
            out = [("row", "Discrim", a), ("row", "Diff", -th[c.ib] / a)]
            return out
        for j in range(1, len(lv)):
            out.append(("row", f"{lv[j - 1]:g}-{lv[j]:g}", -th[c.ic[j - 1]] / a))
        return out

    def probs(self, th, t):
        """Probability of each category (rows) on the grid t."""
        k = self.kind
        if k in ("1pl", "2pl"):
            p = special.expit(th[self.ia] * t + th[self.ic])
            return np.vstack([1 - p, p])
        if k == "3pl":
            g = float(special.expit(th[self.ig]))
            p = g + (1 - g) * special.expit(th[self.ia] * t + th[self.ic])
            return np.vstack([1 - p, p])
        if k == "grm":
            a = th[self.ia]
            cum = [np.ones_like(t)] + [special.expit(a * t - th[i]) for i in self.cuts] + [np.zeros_like(t)]
            return np.vstack([cum[j] - cum[j + 1] for j in range(len(cum) - 1)])
        A, C = self.c.AC(th)
        eta = A[:, None] * t[None, :] + C[:, None]
        return np.exp(eta - lse(eta, axis=0)[None, :])

    def info(self, th, t):
        h = 1e-5
        p = self.probs(th, t)
        d = (self.probs(th, t + h) - self.probs(th, t - h)) / (2 * h)
        return np.sum(d * d / np.maximum(p, 1e-300), axis=0)


def _expected(item: Item, th, t):
    p = item.probs(th, t)
    sc = np.arange(p.shape[0], dtype=float)
    return sc @ p


# ---------------------------------------------------------------------------
# command
# ---------------------------------------------------------------------------

_OPTS = {"sepguessing": 5, "listwise": 4, "intmethod": 4, "intpoints": 4, "level": 1, "notable": 5, "noheader": 6,
         "nolog": 4, "log": 3, "iterate": 4, "vce": 3, "group": 5, "cns": 3, "difficult": 4,
         "technique": 4, "from": 4, "coeflegend": 5, "noci": 4, "nopvalues": 5, "tolerance": 3,
         "ltolerance": 4, "nrtolerance": 4}


def _items_eq(s, kind: str, text: str, taken: set, sepguess: bool = False) -> list[Eq]:
    names = s.expand_varlist(text)
    eqs = []
    for n in names:
        if s.data.get(n).is_string:
            raise StataError(109, "type mismatch")
        if n in taken:
            raise StataError(198, f"{n} is specified more than once")   # VERIFY
        taken.add(n)
        e = Eq(n, "resp")
        if kind in ("1pl", "2pl"):
            e.family, e.link = "bernoulli", "logit"
            e.terms = [("lat", "Theta", None, "_a_1pl" if kind == "1pl" else None)]
        elif kind == "grm":
            e.family, e.link = "ordinal", "logit"
            e.terms = [("lat", "Theta", None, None)]
        elif kind == "3pl":
            e.family, e.link = "bernoulli", "logit"
            e.irt = lambda P, name, y: ThreePL(P, name, y, shared=not sepguess)
        else:
            e.family, e.link = "nominal" if kind == "nrm" else "ordinal", "mlogit"
            e.irt = (lambda k: lambda P, name, y: Cat(P, name, y, k))(kind)
        eqs.append(e)
    return eqs


def _parse(s, args: str):
    c = find_top(args, ",")
    head, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    m = re.match(r"\s*(\w+)\s*(.*)$", head, re.S)
    if not m:
        raise StataError(198, "irt requires a model")   # VERIFY
    word = m.group(1).lower()
    kinds = [k for k in _KINDS + ("hybrid",) if k.startswith(word)]
    if word in _KINDS + ("hybrid",):
        kinds = [word]
    if len(kinds) != 1:
        raise StataError(198, f"irt {word}: invalid model")   # VERIFY
    kind, rest = kinds[0], m.group(2)
    ifm = re.search(r"\s((?:if|in)\b.*)$", " " + rest, re.S)
    body, tail = (rest[:ifm.start()], ifm.group(1)) if ifm else (rest, "")
    p = parse_standard(tail + ("," + opts if opts.strip() else ""))
    return kind, body, p, match_options(p.options, _OPTS)


def cmd_irt(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "gsem", args)
        return
    from ..commands._util import touse
    kind, body, p, o = _parse(s, args)
    for bad in ("group", "cns"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() is not available in OpenDTA yet")   # OpenDTA
    vce = str(o.get("vce") or "oim").strip().lower()
    if vce not in ("oim", "oi", "o"):
        raise StataError(198, f"vce({vce}) is not available in OpenDTA yet")   # OpenDTA
    taken: set = set()
    eqs: list[Eq] = []
    spec: list[tuple[str, list[str]]] = []
    if kind == "hybrid":
        groups = re.findall(r"\(\s*(\w+)\s+([^()]*)\)", body)
        if not groups:
            raise StataError(198, "irt hybrid requires (model varlist) specifications")   # VERIFY
        for mk, vl in groups:
            mk = mk.lower()
            if mk not in _KINDS:
                raise StataError(198, f"irt hybrid: invalid model {mk}")   # VERIFY
            new = _items_eq(s, mk, vl, taken, sepguess=bool(o.get("sepguessing")))
            eqs += new
            spec.append((mk, [e.name for e in new]))
    else:
        if not body.strip():
            raise StataError(100, "varlist required")
        eqs = _items_eq(s, kind, body, taken, sepguess=bool(o.get("sepguessing")))
        spec.append((kind, [e.name for e in eqs]))
    o2 = dict(o)
    o2["_covs"], o2["_varfix"] = [], {"Theta": 1.0}
    mask = touse(s, p)
    m = GModel(s, eqs, ["Theta"], [], [], o2, mask)
    log = not o.get("nolog")
    maxiter = int(o["iterate"]) if o.get("iterate") else 300
    th, V, ll, conv, iters = fit(m, s, log, maxiter)
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", 95))
    kinds = {n: k for k, ns in spec for n in ns}
    items = [Item(kinds[r.eq.name], r.eq.name, r) for r in m.resps]
    title = _TITLES[kind]

    def hook(est, scalars, macros):
        est.display = display_irt
        est.predict = None
        est.estat = {"report": (3, estat_report)}
        est.irt = dict(items=items, title=title, kind=kind, spec=spec,
                       notable=bool(o.get("notable")), noheader=bool(o.get("noheader")))
        # e() in Stata 14's order (case 1909); VERIFY the items and models of hybrid
        from .results import datasignature
        groups = spec if kind == "hybrid" else [(kind, [r.eq.name for r in m.resps])]
        names = [r.eq.name for r in m.resps]
        scalars.append(("irt_k_eq", float(len(groups))))
        for g, (_k, its) in enumerate(groups, 1):
            scalars.append((f"k_items{g}", float(len(its))))
        drop = {"cmdline", "cmd", "title", "predict", "estat_cmd"}
        rest = [x for x in macros if x[0] not in drop]
        lead = [("cmdline", f"irt {args.strip()}"), ("predict", "irt_p"), ("cmd2", "irt"),
                ("estat_cmd", "irt_estat")]
        for g, (k_, its) in enumerate(groups, 1):
            lead += [(f"items{g}", " ".join(its)), (f"model{g}", k_)]
        lead += [("title", title), ("cmd", "gsem")]
        meth = [x for x in rest if x[0] == "method"]
        rest = [x for x in rest if x[0] != "method"]
        sig = [("datasignaturevars", " ".join(names)), ("datasignature", datasignature(s, names, m.sample))]
        macros[:] = lead + meth + sig + rest
    _post(s, args, m, th, V, ll, conv, iters, level, hook=hook, quiet=True)
    est = s._estimates
    display_irt(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _table_rows(est, items=None, byparm=False):
    """Rows with their values, standard errors: [(item or None, kind, label, b, se)]."""
    x = est.extra
    th, V = x["rep"], x["Vr"]
    its = est.irt["items"] if items is None else items

    def vals(t):
        out = []
        for it in its:
            out.extend(v for _k, _l, v in it.rows(t))
        return np.array(out)
    v0 = vals(th)
    J = np.zeros((len(v0), len(th)))
    for i in range(len(th)):
        h = 1e-6 * (abs(th[i]) + 1)
        e = np.zeros(len(th))
        e[i] = h
        J[:, i] = (vals(th + e) - vals(th - e)) / (2 * h)
    se = np.sqrt(np.maximum(np.diag(J @ V @ J.T), 0))
    rows, k = [], 0
    for it in its:
        for kind, lab, _v in it.rows(th):
            rows.append((it, kind, lab, v0[k], se[k]))
            k += 1
    return rows


def _line(out, label, b, s_, c, W=12, left=False):
    out.write((f"{label:<{W}} |" if left else f"{label:>{W}} |"), "text")
    if s_ > 0 and np.isfinite(s_):
        z = b / s_
        pv = 2 * float(st.norm.sf(abs(z)))
        out.write(f"  {_g(b):>9}  {_g(s_):>9}  {z:>7.2f}  {pv:>6.3f}    {_g(b - c * s_):>9}"
                  f"   {_g(b + c * s_):>9}\n", "result")
    else:
        out.write(f"  {_g(b):>9}  {'.':>9}\n", "result")


def display_irt(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    info = est.irt
    lv = level or est.level
    if header:
        out.write("\n" + f"{info['title']:<48}{'Number of obs':<18}= ", "text")
        out.write(f"{x['N']:>10,}\n", "result")
        out.write("Log likelihood = ", "text")
        out.write(f"{_g(x['ll'], '%10.0g')}\n\n", "result")
    if not table:
        return
    _print(s, est, _table_rows(est), lv)


def _print(s, est, rows, lv, byparm=False):
    out = s.output
    c = crit("z", lv, None)
    W = 12
    dash = "-" * (W + 1) + "+" + "-" * 64 + "\n"
    out.write("-" * (W + 66) + "\n", "text")
    out.write(f"{'':>{W}} |      Coef.   Std. Err.      z    P>|z|     [{lv:g}% Conf. Interval]\n", "text")
    out.write(dash, "text")
    its = est.irt["items"]
    if byparm:
        order, seen = [], []
        for it, kind, lab, b, se in rows:
            key = lab if kind == "row" and not lab.startswith(">=") and "-" not in lab else None
        labels = []
        for it, kind, lab, b, se in rows:
            if kind == "head":
                continue
            lbl = lab if lab in ("Discrim", "Diff", "Guess") else "Diff"
            if lbl not in labels:
                labels.append(lbl)
        first = True
        for lbl in labels:
            if not first:
                out.write(dash, "text")
            first = False
            out.write(f"{lbl:<{W}} |\n", "text")
            for it, kind, lab, b, se in rows:
                if kind == "head":
                    continue
                cur = lab if lab in ("Discrim", "Diff", "Guess") else "Diff"
                if cur == lbl:
                    _line(out, it.name if lab == lbl else f"{it.name} {lab}", b, se, c)
        out.write("-" * (W + 66) + "\n", "text")
        return
    shared = [it for it in its if it.shared_a()]
    done_shared: set = set()
    last = None
    first = True
    for it, kind, lab, b, se in rows:
        if it.shared_a() and lab == "Discrim":
            key = it.kind
            if key in done_shared:
                continue
            done_shared.add(key)
            if not first:
                out.write(dash, "text")
            first = False
            _line(out, "Discrim", b, se, c)
            last = None
            continue
        if last is not it:
            if not first:
                out.write(dash, "text")
            first = False
            out.write(f"{it.name:<{W}} |\n", "text")
            last = it
        if kind == "head":
            out.write(f"{lab:>{W}} |\n", "text")
        else:
            _line(out, lab, b, se, c)
    # rating scale: the thresholds shared by the items
    rs = next((it for it in its if it.kind == "rsm"), None)
    if rs is not None:
        th, V = est.extra["rep"], est.extra["Vr"]
        ia = rs.c.ia[0]
        out.write(dash, "text")
        out.write(f"{'Kappa':<{W}} |\n", "text")
        for j, ie in enumerate(rs.c.ie):
            def f(t, ie=ie):
                return -t[ie] / t[ia]
            b = f(th)
            g = np.zeros(len(th))
            for i in (ie, ia):
                e = np.zeros(len(th))
                h = 1e-6 * (abs(th[i]) + 1)
                e[i] = h
                g[i] = (f(th + e) - f(th - e)) / (2 * h)
            _line(out, f"{rs.c.levels[j + 2]:g}", b, float(np.sqrt(max(g @ V @ g, 0))), c)
    out.write("-" * (W + 66) + "\n", "text")


def estat_report(s: "Session", est: Estimates, args: str) -> None:
    """estat report [varlist] [, byparm sort(a|b|c|#) descending level(#)]"""
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    o = match_options(opt, {"byparm": 2, "sort": 3, "descending": 4, "level": 1, "alabel": 2,
                            "blabel": 2, "cblabel": 2, "seqlabel": 3, "post": 4}) if opt.strip() else {}
    its = est.irt["items"]
    if body.strip():
        names = s.expand_varlist(body)
        sel = []
        for n in names:
            hit = [it for it in its if it.name == n]
            if not hit:
                raise StataError(111, f"{n} not found among the items")   # VERIFY
            sel.append(hit[0])
        its = sel
    th = est.extra["rep"]
    key = str(o.get("sort") or "").strip().lower()
    if key:
        def keyf(it):
            rows = it.rows(th)
            vals = {l: v for k, l, v in rows if k == "row"}
            if key == "a":
                return vals.get("Discrim", 0.0)
            if key == "b":
                ds = [v for k, l, v in rows if k == "row" and l != "Discrim" and l != "Guess"]
                return float(np.mean(ds)) if ds else 0.0
            if key == "c":
                return vals.get("Guess", 0.0)
            raise StataError(198, f"sort({key}) invalid")   # VERIFY
        its = sorted(its, key=keyf, reverse=bool(o.get("descending")))
    lv = float(o["level"]) if o.get("level") else est.level
    s.output.write("\n", "text")
    _print(s, est, _table_rows(est, its), lv, byparm=bool(o.get("byparm")))


# ---------------------------------------------------------------------------
# irtgraph
# ---------------------------------------------------------------------------

def cmd_irtgraph(s: "Session", args: str) -> None:
    from .results import current
    est = current(s)
    if est.cmd != "gsem" or not hasattr(est, "irt"):
        raise StataError(301, "last estimates not found")
    m = re.match(r"\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "irtgraph requires icc, iif, tcc or tif")   # VERIFY
    sub = m.group(1).lower()
    names = [w for w in ("icc", "iif", "tcc", "tif") if w.startswith(sub) and len(sub) >= 2]
    if sub in ("icc", "iif", "tcc", "tif"):
        names = [sub]
    if len(names) != 1:
        raise StataError(198, f"irtgraph {sub}: invalid subcommand")   # VERIFY
    sub = names[0]
    rest = m.group(2)
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    rng = re.search(r"\bxrange\(([^)]*)\)|\brange\(([^)]*)\)", opt)
    lo, hi = -4.0, 4.0
    if rng:
        vals = (rng.group(1) or rng.group(2)).split()
        lo, hi = float(vals[0]), float(vals[1])
        opt = opt[:rng.start()] + opt[rng.end():]
    its_all = est.irt["items"]
    if body.strip():
        want = s.expand_varlist(body)
        its = []
        for n in want:
            hit = [it for it in its_all if it.name == n]
            if not hit:
                raise StataError(111, f"{n} not found among the items")   # VERIFY
            its.append(hit[0])
    else:
        its = list(its_all)
    th = est.extra["rep"]
    t = np.linspace(lo, hi, 401)
    curves: list[tuple[str, np.ndarray]] = []
    if sub == "icc":
        for it in its:
            p = it.probs(th, t)
            if p.shape[0] == 2:
                curves.append((it.name, p[1]))
            else:
                lv = it.levels if it.levels is not None else np.arange(p.shape[0])
                for k in range(p.shape[0]):
                    curves.append((f"{it.name} {lv[k]:g}", p[k]))
        title, yt = "Item Characteristic Curves", "Probability"
    elif sub == "iif":
        curves = [(it.name, it.info(th, t)) for it in its]
        title, yt = "Item Information Functions", "Information"
    elif sub == "tcc":
        curves = [("Test characteristic curve", sum(_expected(it, th, t) for it in its))]
        title, yt = "Test Characteristic Curve", "Expected score"
    else:
        curves = [("Test information function", sum(it.info(th, t) for it in its))]
        title, yt = "Test Information Function", "Information"
    _plot(s, t, curves, title, yt, opt)


def _plot(s, t, curves, title, ytitle, opt: str) -> None:
    from ..core.dataset import Variable
    ds = s.data
    s.run("preserve")
    try:
        s.run("quietly clear")
        s.run(f"quietly set obs {len(t)}")
        ds = s.data
        ds.add(Variable("theta", "double", np.asarray(t, float)))
        plots, labels = [], []
        for k, (lab, y) in enumerate(curves, 1):
            ds.add(Variable(f"c{k}", "double", np.asarray(y, float)))
            plots.append(f"(line c{k} theta)")
            labels.append(f'label({k} "{lab}")')
        user = opt or ""
        defaults = []
        if not re.search(r"\btitle\(", user):
            defaults.append(f'title("{title}")')
        if not re.search(r"\bxtitle\(", user):
            defaults.append('xtitle("Theta")')
        if not re.search(r"\bytitle\(", user):
            defaults.append(f'ytitle("{ytitle}")')
        if not re.search(r"\blegend\(", user):
            defaults.append("legend(" + " ".join(labels) + ")" if len(curves) > 1 else "legend(off)")
        sep = "," if not user.strip().startswith(",") else ""
        s.run(f"twoway {' '.join(plots)}, {' '.join(defaults)} {user}".replace(", ,", ","))
    finally:
        s.run("restore")
