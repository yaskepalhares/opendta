"""Nonlinear estimation with substitutable expressions: nl, nlsur and gmm ([R] nl, [R] nlsur,
[R] gmm). A parameter is written {name}, {name=init} or {lc: varlist [, noconstant]} inside an
expression; the expression is evaluated with Stata's own evaluator at the current parameter
values and the derivatives are numerical.

Layouts come from the manuals; the lines marked VERIFY wait for the Stata 14 log of case 1908."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.varlist import resolve_name
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from .mleq import aux_row
from .models2 import g
from .regress import cmdline, comma
from .results import Estimates, coef_table, post

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_BRACE = re.compile(r"\{([^{}]*)\}")


def _lit(v) -> str:
    """A number as text the expression parser reads back exactly."""
    return format(float(v), ".17g")


class SubExpr:
    """An expression with parameters. `names` are the parameters in order of appearance (a linear
    combination {xb: x1 x2} contributes xb_x1, xb_x2 and xb__cons)."""

    def __init__(self, s: "Session", text: str):
        self.s = s
        self.text = text
        self.names: list[str] = []
        self.init: dict[str, float | None] = {}
        self.lcs: dict[str, list[tuple[str, str]]] = {}     # label -> [(param, variable or "")]
        for m in _BRACE.finditer(text):
            self._add(m.group(1))

    def _add(self, body: str) -> None:
        body = body.strip()
        if ":" in body:
            label, rest = body.split(":", 1)
            label = label.strip() or "xb"
            noconstant = False
            if "," in rest:
                rest, opt = rest.split(",", 1)
                noconstant = opt.strip().startswith("nocons")
            if label in self.lcs:
                return
            parts = []
            for v in rest.split():
                v = resolve_name(self.s.data, v)
                parts.append((f"{label}_{v}", v))
            if not noconstant:
                parts.append((f"{label}__cons", ""))
            self.lcs[label] = parts
            for nm, _ in parts:
                self._register(nm, None)
            return
        init = None
        if "=" in body:
            name, val = body.split("=", 1)
            name = name.strip()
            try:
                init = float(val.strip())
            except ValueError:
                raise StataError(198, f"invalid initial value for {name}")   # VERIFY
        else:
            name = body
        if not re.fullmatch(r"[A-Za-z_]\w*", name):
            raise StataError(198, f"invalid parameter name {{{body}}}")   # VERIFY
        self._register(name, init)

    def _register(self, name: str, init: float | None) -> None:
        if name not in self.names:
            self.names.append(name)
            self.init[name] = init
        elif init is not None:
            self.init[name] = init

    def variables(self) -> list[str]:
        return [v for parts in self.lcs.values() for _, v in parts if v]

    def render(self, theta: dict) -> str:
        def sub(m):
            body = m.group(1).strip()
            if ":" in body:
                label = body.split(":", 1)[0].strip() or "xb"
                terms = []
                for nm, v in self.lcs[label]:
                    terms.append(f"{_lit(theta[nm])}*{v}" if v else f"{_lit(theta[nm])}")
                return "(" + " + ".join(terms) + ")"
            name = body.split("=", 1)[0].strip()
            return f"({_lit(theta[name])})"
        return _BRACE.sub(sub, self.text)

    def evaluate(self, theta: dict) -> np.ndarray:
        from ..commands._util import eval_vector
        v = eval_vector(self.s, self.render(theta))
        v = np.asarray(v, dtype=float)
        n = self.s.data.nobs
        if v.ndim == 0:
            v = np.full(n, float(v))
        return v


def _split_equations(text: str) -> tuple[list[str], str]:
    """'(a) (b) if x, opts' -> ['a', 'b'], 'if x, opts'."""
    t = text.strip()
    eqs, depth, start, end = [], 0, None, 0
    for i, ch in enumerate(t):
        if ch == "(":
            if depth == 0:
                start = i
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0 and start is not None:
                eqs.append(t[start + 1:i])
                end = i + 1
        elif depth == 0 and ch not in " \t":
            break
    return eqs, t[end:].strip()


def _lhs_rhs(body: str) -> tuple[str, str]:
    depth = 0
    for i, ch in enumerate(body):
        if ch in "({[":
            depth += 1
        elif ch in ")}]":
            depth -= 1
        elif ch == "=" and depth == 0 and body[i - 1:i + 1] != "==" and body[i:i + 2] != "==":
            return body[:i].strip(), body[i + 1:].strip()
    raise StataError(198, "equation must be specified as (depvar = expression)")   # VERIFY


def _num_jac(fun, theta: np.ndarray, h: float = 1e-6) -> np.ndarray:
    """Jacobian of fun(theta) (n values) by central differences."""
    cols = []
    for j in range(len(theta)):
        step = h * max(abs(theta[j]), 1e-2)
        up, dn = theta.copy(), theta.copy()
        up[j] += step
        dn[j] -= step
        cols.append((fun(up) - fun(dn)) / (2 * step))
    return np.column_stack(cols)


def gauss_newton(resid, theta0: np.ndarray, W: np.ndarray, *, eps: float = 1e-5, maxit: int = 100, log=None):
    """Minimizes r' W r with r = resid(theta) (W a vector of weights or a matrix) by Gauss-Newton with
    step halving. `log(it, value)` is called at every iteration. Returns (theta, J, r, iterations)."""
    def quad(r):
        return float(r @ (W * r)) if W.ndim == 1 else float(r @ W @ r)

    def wmul(A):
        return A * W[:, None] if W.ndim == 1 else W @ A

    theta = theta0.copy()
    f = lambda t: -resid(t)                         # the model: r = y - f
    r = resid(theta)
    q = quad(r)
    it = 0
    while True:
        J = _num_jac(f, theta)
        if log:
            log(it, q)
        WJ = wmul(J)
        A = J.T @ WJ
        d = np.linalg.solve(A + 1e-12 * np.trace(A) / max(len(theta), 1) * np.eye(len(theta)), WJ.T @ r)
        crit = float(r @ WJ @ np.linalg.solve(A + 1e-14 * np.eye(len(theta)), WJ.T @ r)) / max(q, 1e-300)
        if crit < eps * eps or it >= maxit:
            break
        lam = 1.0
        while lam > 1e-10:
            tn = theta + lam * d
            rn = resid(tn)
            qn = quad(rn)
            if np.isfinite(qn) and qn <= q:
                break
            lam /= 2
        else:
            break
        converged_rel = abs(q - qn) <= 1e-12 * max(q, 1e-300)
        theta, r, q = tn, rn, qn
        it += 1
        if converged_rel:
            if log:
                log(it, q)
            J = _num_jac(f, theta)
            break
    return theta, J, r, it


# ---------------------------------------------------------------------------
# nl
# ---------------------------------------------------------------------------

class ProgFunc:
    """The function evaluator program of nl, nlsur and gmm ([R] nl, [R] nlsur, [R] gmm), called
    as `prog newvars varlist if touse, at(matrix)`: it fills the m new variables (the fitted
    values of nl and nlsur, the residuals of gmm). The program runs on every observation and the
    caller keeps those of its sample. `names` are the parameters (parameters(), or b1, b2, ...
    with nparameters(#)); `eq(q)` is equation q with the interface of SubExpr."""

    def __init__(self, s: "Session", prog: str, m: int, names: list[str], varlist: list[str],
                 noisy_first: bool = False):
        self.s = s
        self.noisy_first = noisy_first      # gmm shows the output of the first call (compat 1923)
        self.prog = prog
        self.m = m
        self.names = names
        self.init: dict[str, float | None] = {nm: None for nm in names}
        self.vars = varlist
        self._key = None
        self._val = None

    def evaluate_all(self, theta: dict) -> np.ndarray:
        key = tuple(float(theta[nm]) for nm in self.names)
        if key == self._key:
            return self._val
        from ..commands.matrix import Matrix, _store
        from ..commands.syntaxcmd import _temp_name
        from ..core.dataset import Variable
        s, ds = self.s, self.s.data
        tmp = []
        for i in range(self.m + 1):
            nm = _temp_name(s)
            ds.add(Variable(nm, "double", np.ones(ds.nobs) if i == 0 else np.full(ds.nobs, SYS)))
            tmp.append(nm)
        mats = _store(s)
        mats["__nlat"] = Matrix(np.array([key]), ["r1"], list(self.names))
        try:
            call = f"{self.prog} {' '.join(tmp[1:] + self.vars)} if {tmp[0]}, at(__nlat)"
            if self.noisy_first:
                self.noisy_first = False
                s.interp._execute_expanded(call)
            else:
                with s.output.quietly():
                    s.interp._execute_expanded(call)
            val = np.column_stack([ds.get(nm).data.astype(float) for nm in tmp[1:]])
        finally:
            ds.drop_vars([nm for nm in tmp if ds.has(nm)])
            mats.pop("__nlat", None)
        self._key, self._val = key, val
        return val

    def eq(self, q: int) -> "_ProgEq":
        return _ProgEq(self, q)


class _ProgEq:
    def __init__(self, pf: ProgFunc, q: int):
        self.pf, self.q = pf, q
        self.names = pf.names
        self.init = pf.init
        self.lcs: dict = {}
        self.text = ""

    def variables(self) -> list[str]:
        return list(self.pf.vars)

    def evaluate(self, theta: dict) -> np.ndarray:
        return self.pf.evaluate_all(theta)[:, self.q]


_PROG_FORM = re.compile(r"^\s*([A-Za-z_]\w*)\s*([@:])(.*)$", re.S)


def _program(s: "Session", name: str) -> str:
    if name not in s.programs:
        raise StataError(199, f"unrecognized command:  {name}")
    return name


def _split_varlist(s: "Session", text: str) -> tuple[list[str], str]:
    """'y x1 x2 if z, opts' -> (['y', 'x1', 'x2'], 'if z, opts')."""
    text = text.strip()
    cut = len(text)
    for mm in re.finditer(r",|\[|(?<!\S)(?:if|in)(?!\S)", text):
        cut = mm.start()
        break
    head = text[:cut].split()
    tail = text[cut:].strip()
    from ..core.varlist import expand
    return (expand(s.data, " ".join(head)) if head else []), tail


def _param_names(opts_text: str) -> list[str]:
    """parameters(namelist) or nparameters(#) (b1, b2, ...) of a function evaluator program."""
    m = re.search(r"\bparam\w*\s*\(([^)]*)\)", opts_text)
    if m and not re.search(r"\bnparam\w*\s*\(", opts_text[:m.start() + 1]):
        return [nm.split(":")[-1] for nm in m.group(1).split()]
    m = re.search(r"\bnparam\w*\s*\(\s*(\d+)\s*\)", opts_text)
    if m:
        return [f"b{i}" for i in range(1, int(m.group(1)) + 1)]
    raise StataError(198, "you must specify parameters() or nparameters()")   # VERIFY


def _int_option(opts_text: str, stem: str) -> int | None:
    m = re.search(rf"\b{stem}\w*\s*\(\s*(\d+)\s*\)", opts_text)
    return int(m.group(1)) if m else None


def _sexp_program(s: "Session", prog: str, varlist: list[str], tail: str) -> None:
    """Runs a substitutable expression program (`nl name : varlist`), whose r() holds the
    equations."""
    from ..lang.syntax import parse_standard as _ps
    p = _ps(tail)
    cond = ""
    if p.if_:
        cond += f" if {p.if_}"
    if p.in_:
        cond += f" in {p.in_}"
    with s.output.quietly():
        s.interp._execute_expanded(f"{prog} {' '.join(varlist)}{cond}")


def _options_tail(tail: str) -> str:
    i = tail.find(",")
    return tail[i + 1:] if i >= 0 else ""


def cmd_nl(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..commands.summarize import weights as get_weights
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "nl", args)
        return
    ds = s.data
    orig_args = args
    pm = _PROG_FORM.match(args)
    funcprog = ""
    nl_type = "3" if pm and pm.group(2) == "@" else ("2" if pm else "1")
    if pm and pm.group(2) == ":":
        # substitutable expression program: nl<name> leaves the equation in r(eq) and a title in
        # r(title)
        prog = _program(s, "nl" + pm.group(1))
        varlist, tail = _split_varlist(s, pm.group(3))
        _sexp_program(s, prog, varlist, tail)
        eq = str(s.r.get("eq", ""))
        if not eq:
            raise StataError(198, f"{prog} did not return r(eq)")   # VERIFY
        if "=" not in eq.split("{")[0]:
            eq = f"{varlist[0]} = {eq}"
        title = str(s.r.get("title", ""))
        if title and not re.search(r"\btitle\s*\(", tail):
            tail = (tail + ("" if "," in tail else ",") + f' title("{title}")')
        args = f"({eq}) {tail}"
    if pm and pm.group(2) == "@":
        # function evaluator program: nl<name> fills the dependent variable (a temporary one)
        funcprog = _program(s, "nl" + pm.group(1))
        varlist, tail = _split_varlist(s, pm.group(3))
        if not varlist:
            raise StataError(102, "too few variables specified")
        depname, rhs = varlist[0], ""
        se = ProgFunc(s, funcprog, 1, _param_names(_options_tail(tail)), varlist[1:]).eq(0)
    else:
        eqs, tail = _split_equations(args)
        if len(eqs) != 1:
            raise StataError(198, "nl requires one equation: (depvar = expression)")   # VERIFY
        dep, rhs = _lhs_rhs(eqs[0])
        depname = resolve_name(ds, dep)
        se = SubExpr(s, rhs)
        if not se.names:
            raise StataError(198, "no parameters in the expression")   # VERIFY
    p = parse_standard(tail)
    o = match_options(p.options, {"initial": 4, "nolog": 4, "log": 3, "eps": 3, "iterate": 4, "level": 1,
                                  "variables": 3, "lnlsq": 3, "title": 5, "robust": 1, "vce": 3, "cluster": 2,
                                  "noconstant": 4, "hasconstant": 3, "nocnsreport": 5, "noheader": 6,
                                  "notable": 5, "delta": 3, "parameters": 5, "nparameters": 6}) \
        if p.options.strip() else {}
    mask = touse(s, p)
    y_all = ds.get(depname).data.astype(float)
    mask &= y_all < SYS
    for v in se.variables():
        mask &= ds.get(v).data < SYS
    w, wtype, mask = get_weights(s, p, mask, ("aweight", "fweight", "iweight", "pweight"))
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    init = dict(se.init)
    if o.get("initial"):
        for m in re.finditer(r"(\w+)\s+([-+0-9.eE]+)", str(o["initial"])):
            init[m.group(1)] = float(m.group(2))
    theta0 = np.array([0.0 if init.get(nm) is None else init[nm] for nm in se.names])
    for label, parts in se.lcs.items():
        if all(init.get(nm) is None for nm, _ in parts):
            # linear combinations start at the regression of the dependent variable on the variables
            X = np.column_stack([ds.get(v).data[mask] if v else np.ones(int(mask.sum())) for _, v in parts])
            bb = np.linalg.lstsq(X, y_all[mask], rcond=None)[0]
            for (nm, _), v in zip(parts, bb):
                theta0[se.names.index(nm)] = float(v)
    with np.errstate(all="ignore"):
        f0 = se.evaluate(dict(zip(se.names, theta0)))
    mask &= f0 < SYS
    mask &= np.isfinite(f0)
    n = int(mask.sum())
    if n == 0:
        raise StataError(2000, "no observations")
    wv = w[mask].astype(float) if w is not None else np.ones(n)
    if wtype in ("aweight", "pweight"):
        wv = wv * n / wv.sum()
    N = float(wv.sum()) if wtype in ("fweight", "iweight") else float(n)
    y = y_all[mask]
    out = s.output
    k = len(se.names)
    out.write(f"(obs = {comma(N)})\n", "text")

    def f_theta(t):
        with np.errstate(all="ignore"):
            return se.evaluate(dict(zip(se.names, t)))[mask]

    def resid(t):
        return y - f_theta(t)

    def log(it, q):
        if not o.get("nolog"):
            if it == 0:
                out.write("\n", "text")
            out.write(f"Iteration {it}:  residual SS = {g(q, '%9.0g'):>9}\n", "text")

    theta, J, r, iters = gauss_newton(resid, theta0, wv, eps=float(o["eps"]) if o.get("eps") else 1e-5,
                                      maxit=int(o["iterate"]) if o.get("iterate") else 100, log=log)
    if not o.get("nolog"):
        out.write("\n", "text")
    rss = float(r @ (wv * r))
    df_r = N - k
    s2 = rss / df_r
    A = J.T @ (J * wv[:, None])
    Ainv = np.linalg.inv(A)
    V = s2 * Ainv
    vcetype = ""
    if o.get("robust") or str(o.get("vce", "")).startswith("rob") or wtype == "pweight":
        sc = J * (wv * r)[:, None]
        V = (N / (N - k)) * Ainv @ (sc.T @ sc) @ Ainv
        vcetype = "Robust"
    V = (V + V.T) / 2
    # a parameter that enters additively is the constant term of the model
    const_idx = None
    for j in range(k):
        if np.allclose(J[:, j], 1.0, atol=1e-6):
            const_idx = j
            break
    if o.get("noconstant"):
        const_idx = None
    if const_idx is not None or o.get("hasconstant"):
        ybar = float(np.sum(wv * y) / wv.sum())
        tss = float(np.sum(wv * (y - ybar) ** 2))
        df_m = k - 1
        df_tot = N - 1
    else:
        tss = float(np.sum(wv * y ** 2))
        df_m = k
        df_tot = N
    mss = tss - rss
    r2 = mss / tss if tss > 0 else SYS
    r2_a = 1 - (1 - r2) * df_tot / df_r if tss > 0 else SYS
    rmse = float(np.sqrt(s2))
    dev = N * (1 + np.log(2 * np.pi) + np.log(rss / N))
    est = Estimates("nl", depname, list(se.names), theta, V, int(N), "t", df_r)
    est.level, est.vcetype, est.clustvar, est.n_clust = level, vcetype, "", None
    est.eqnames = ["/"] * k
    est.rows = [aux_row("/" + nm, float(theta[i]), float(np.sqrt(max(V[i, i], 0))), level) for i, nm in enumerate(se.names)]
    for r_ in est.rows:
        r_.noeform = True
    est.extra = {"N": N, "mss": mss, "rss": rss, "tss": tss, "df_m": df_m, "df_r": df_r, "df_t": df_tot,
                 "r2": r2, "r2_a": r2_a, "rmse": rmse, "dev": dev, "const": se.names[const_idx] if const_idx is not None else ""}
    est.sexp = se
    est.title = str(o.get("title") or "").strip().strip('"')
    est.display = _display_nl

    def predict(s_, e, stat, mask_):
        t = dict(zip(e.sexp.names, e.b))
        v = e.sexp.evaluate(t)
        if stat in ("", "yhat"):
            return v, ("(option yhat assumed; fitted values)" if not stat else "")   # compat 1923
        if stat == "residuals":
            return ds.get(e.depvar).data.astype(float) - v, ""
        raise StataError(198, f"option {stat} not allowed")   # VERIFY
    est.predict = predict
    est.predict_opts = {"yhat": 1, "residuals": 3}
    # e() as Stata lists it (compat 1923)
    scalars = [("df_m", df_m), ("df_r", df_r), ("mss", mss), ("mms", mss / df_m if df_m else SYS),
               ("msr", rss / df_r), ("rmse", rmse), ("ll", -dev / 2), ("r2", r2), ("r2_a", r2_a), ("dev", dev),
               ("converge", 1.0), ("N", N), ("rss", rss), ("df_t", df_tot), ("tss", tss),
               ("cj", float(const_idx + 1) if const_idx is not None else 0.0), ("gm_2", 1.0),   # VERIFY gm_2
               ("delta", float(o["delta"]) if o.get("delta") else 4e-7), ("k", float(k)), ("k_eq_model", 0.0),
               ("rank", float(np.linalg.matrix_rank(V))), ("log_t", 1.0 if o.get("lnlsq") else 0.0), ("ic", iters)]
    covs = " ".join(se.variables()) if nl_type == "3" else " ".join(dict.fromkeys(
        w for w in re.findall(r"[A-Za-z_]\w*", _BRACE.sub(" ", rhs)) if ds.has(w)))
    macros = [("cmdline", cmdline("nl", orig_args)), ("cmd", "nl"), ("predict", "nl_p"),
              ("marginsnotok", "SCores Residuals"), ("marginsok", "default Yhat"), ("vce", "gnr")]
    if funcprog:
        macros.append(("funcprog", funcprog))
    macros += [("covariates", covs), ("rhs", covs), ("type", nl_type)]
    if est.title:
        macros.append(("title", est.title))
    if nl_type != "3":
        # the expression as Stata stores it: "{a} *xp^ {b}" (compat 1923; VERIFY the general rule)
        def _sx(mm_):
            return "{" + mm_.group(1).split("=")[0].strip() + "}" + (" " if mm_.start() == 0 else "")
        sx = re.sub(r"\s+", "", rhs)
        sx = re.sub(r"\{([^}]*)\}", lambda mm_: _sx(mm_) if mm_.start() == 0 else " " + _sx(mm_), sx)
        macros.append(("sexp", sx))
    macros += [("params", " ".join(se.names)), ("depvar", depname), ("properties", "b V")]
    from ..commands.matrix import Matrix
    post(s, est, scalars, macros, mask, extra_mats=[("init", Matrix(theta0.reshape(1, -1), ["y1"], list(se.names)))])
    _display_nl(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _display_nl(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    out = s.output
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        out.write("\n", "text")
        msm = x["mss"] / x["df_m"] if x["df_m"] > 0 else SYS

        def ss_row(nm, ss, df, ms):
            return f"{nm:>12} |{g(ss, '%10.0g'):>11}{int(df):>11}{g(ms, '%11.0g'):>12}"

        left = [f"{'Source':>12} |      SS            df       MS",
                "-" * 13 + "+" + "-" * 34,
                ss_row("Model", x["mss"], x["df_m"], msm),
                ss_row("Residual", x["rss"], x["df_r"], x["rss"] / x["df_r"]),
                "-" * 13 + "+" + "-" * 34,
                ss_row("Total", x["tss"], x["df_t"], x["tss"] / x["df_t"])]
        rights = {1: ("Number of obs", comma(x["N"])), 2: ("R-squared", g(x["r2"], "%9.4f")),
                  3: ("Adj R-squared", g(x["r2_a"], "%9.4f")), 4: ("Root MSE", g(x["rmse"], "%9.0g")),
                  5: ("Res. dev.", g(x["dev"], "%9.0g"))}
        for i, lft in enumerate(left):
            out.write(lft, "text")
            if i in rights:
                lab, val = rights[i]
                out.write(f"{'':<{48 - len(lft)}}    {lab:<14}=", "text")
                out.write(f"{val:>11}", "result")
            out.write("\n", "text")
        if est.title:
            out.write(f"{est.title}\n", "text")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=lev, vcetype=est.vcetype)
        if x["const"]:
            out.write(f"  Parameter {x['const']} taken as constant term in model & ANOVA table\n", "text")


# ---------------------------------------------------------------------------
# nlsur
# ---------------------------------------------------------------------------

def cmd_nlsur(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "nlsur", args)
        return
    ds = s.data
    orig_args = args
    pm = _PROG_FORM.match(args)
    funcprog = ""
    if pm and pm.group(2) == ":":
        # substitutable expression program: nlsur<name> returns r(n_eq) and the equations in
        # r(eq_1), r(eq_2), ... (VERIFY)
        prog = _program(s, "nlsur" + pm.group(1))
        varlist, tail = _split_varlist(s, pm.group(3))
        _sexp_program(s, prog, varlist, tail)
        try:
            neq = int(float(s.r.get("n_eq", 0)))
        except (TypeError, ValueError):
            neq = 0
        if neq < 1:
            raise StataError(198, f"{prog} did not return r(n_eq)")   # VERIFY
        args = " ".join(f"({s.r.get(f'eq_{q}', '')})" for q in range(1, neq + 1)) + " " + tail
    deps, ses = [], []
    if pm and pm.group(2) == "@":
        # function evaluator program: nlsur<name> fills the dependent variables (temporary ones)
        funcprog = _program(s, "nlsur" + pm.group(1))
        varlist, tail = _split_varlist(s, pm.group(3))
        neq = _int_option(_options_tail(tail), "neq")
        if not neq:
            raise StataError(198, "option nequations() required")   # VERIFY
        if len(varlist) < neq:
            raise StataError(102, "too few variables specified")
        deps = varlist[:neq]
        pf = ProgFunc(s, funcprog, neq, _param_names(_options_tail(tail)), varlist[neq:])
        ses = [pf.eq(q) for q in range(neq)]
        eqs = ["" for _ in range(neq)]
    else:
        eqs, tail = _split_equations(args)
        if not eqs:
            raise StataError(198, "equation required")   # VERIFY
        for body in eqs:
            dep, rhs = _lhs_rhs(body)
            deps.append(resolve_name(ds, dep))
            ses.append(SubExpr(s, rhs))
    names: list[str] = []
    init: dict = {}
    for se in ses:
        for nm in se.names:
            if nm not in names:
                names.append(nm)
            if se.init.get(nm) is not None:
                init[nm] = se.init[nm]
    p = parse_standard(tail)
    o = match_options(p.options, {"fgnls": 2, "nls": 3, "initial": 4, "nolog": 4, "log": 3, "eps": 3,
                                  "iterate": 4, "level": 1, "variables": 3, "noconstant": 4, "title": 5,
                                  "noheader": 6, "notable": 5, "delta": 3, "robust": 1, "nequations": 3,
                                  "parameters": 5, "nparameters": 6, "lnlsq": 3, "ifgnls": 3}) \
        if p.options.strip() else {}
    mask = touse(s, p)
    for dn in deps:
        mask &= ds.get(dn).data < SYS
    for se in ses:
        for v in se.variables():
            mask &= ds.get(v).data < SYS
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    theta0 = np.array([init.get(nm, 0.0) for nm in names])
    for se in ses:
        for label, parts in se.lcs.items():
            if all(nm not in init for nm, _ in parts):
                X = np.column_stack([ds.get(v).data[mask] if v else np.ones(int(mask.sum())) for _, v in parts])
                bb = np.linalg.lstsq(X, ds.get(deps[ses.index(se)]).data[mask], rcond=None)[0]
                for (nm, _), v in zip(parts, bb):
                    theta0[names.index(nm)] = float(v)
    for se in ses:
        with np.errstate(all="ignore"):
            f0 = se.evaluate({nm: theta0[names.index(nm)] for nm in se.names})
        mask &= (f0 < SYS) & np.isfinite(f0)
    n = int(mask.sum())
    m = len(eqs)
    Y = np.column_stack([ds.get(dn).data[mask].astype(float) for dn in deps])
    out = s.output
    out.write(f"(obs = {comma(n)})\n\n", "text")

    def fit_vals(t):
        cols = []
        with np.errstate(all="ignore"):
            for se in ses:
                cols.append(se.evaluate({nm: t[names.index(nm)] for nm in se.names})[mask])
        return np.column_stack(cols)

    def resid_vec(t):
        return (Y - fit_vals(t)).T.reshape(-1)            # equation by equation

    def logger(label):
        def log(it, q):
            if not o.get("nolog"):
                lab = "Residual SS" if label == "nls" else "Scaled RSS"
                out.write(f"Iteration {it}:  {lab} = {g(q, '%9.0g'):>9}\n", "text")
        return log

    if not o.get("nolog"):
        out.write("Calculating NLS estimates...\n", "text")     # VERIFY
    I = np.ones(n * m)
    theta, J, r, it = gauss_newton(resid_vec, theta0, I, log=logger("nls"))
    method = "nls" if o.get("nls") else ("ifgnls" if o.get("ifgnls") else "fgnls")
    sigma = None
    if method != "nls":
        for rounds in range(1 if method == "fgnls" else 50):
            E = (Y - fit_vals(theta))
            sigma = E.T @ E / n
            Wm = np.kron(np.linalg.inv(sigma), np.eye(n))
            if not o.get("nolog"):
                out.write("Calculating FGNLS estimates...\n" if rounds == 0 else "", "text")   # VERIFY
            tnew, J, r, it = gauss_newton(resid_vec, theta, Wm, log=logger("fgnls"))
            done = np.max(np.abs(tnew - theta)) < 1e-8
            theta = tnew
            if done:
                break
        E = (Y - fit_vals(theta))
        if method == "ifgnls":
            sigma = E.T @ E / n
        Wm = np.kron(np.linalg.inv(sigma), np.eye(n))
        J = _num_jac(lambda t: -resid_vec(t), theta)
        A = J.T @ Wm @ J
        V = np.linalg.inv(A)
    else:
        # Stata's nls: the equations weighted by the inverse of their residual variances from the
        # unweighted fit (no covariances), one more Gauss-Newton fit, V = (J' W J)^-1 (compat 1923)
        E = (Y - fit_vals(theta))
        sigma = np.diag(np.diag(E.T @ E / n))
        Wm = np.kron(np.linalg.inv(sigma), np.eye(n))
        theta, J, r, it = gauss_newton(resid_vec, theta, Wm, log=None)
        J = _num_jac(lambda t: -resid_vec(t), theta)
        V = np.linalg.inv(J.T @ Wm @ J)
    V = (V + V.T) / 2
    k = len(names)
    E = Y - fit_vals(theta)
    info = []
    Jf = _num_jac(lambda t: -resid_vec(t), theta)
    for i in range(m):
        rss = float(E[:, i] @ E[:, i])
        yy = Y[:, i]
        tss = float(np.sum((yy - yy.mean()) ** 2))
        kp = len([nm for nm in names if nm in ses[i].names])
        const = ""
        for j, nm in enumerate(names):
            if nm in ses[i].names and np.allclose(np.abs(Jf[i * n:(i + 1) * n, j]), 1.0, atol=1e-6):
                const = nm
                break
        if funcprog:
            # a function evaluator program: Stata cannot tell the parameters or the constant of each
            # equation, and reports the uncentered R-squared (compat 1923)
            tss_u = float(yy @ yy)
            info.append({"obs": n, "parms": ".", "rmse": float(np.sqrt(rss / n)), "r2": 1 - rss / tss_u,
                         "rss": rss, "mss": tss_u - rss, "const": "(none)", "uncentered": True})
        else:
            info.append({"obs": n, "parms": kp, "rmse": float(np.sqrt(rss / n)), "r2": 1 - rss / tss if tss > 0 else SYS,
                         "rss": rss, "mss": tss - rss, "const": const})
    est = Estimates("nlsur", deps[0], list(names), theta, V, n, "z", None)
    est.level, est.vcetype, est.clustvar, est.n_clust = level, "", "", None
    est.eqnames = ["/"] * k
    est.rows = [aux_row("/" + nm, float(theta[i]), float(np.sqrt(max(V[i, i], 0))), level) for i, nm in enumerate(names)]
    est.extra = {"N": n, "info": info, "deps": deps, "method": method, "sigma": sigma}
    est.sexps, est.pnames = ses, names
    est.display = _display_nlsur

    def predict(s_, e, stat, mask_):
        t = dict(zip(e.pnames, e.b))
        v = e.sexps[0].evaluate(t)
        return v, ""
    est.predict = predict
    # e() as Stata lists it (compat 1923, function evaluator program; VERIFY the other forms)
    scalars = [("N", n)]
    for i, e in enumerate(info, 1):
        scalars += [(f"rss_{i}", e["rss"]), (f"mss_{i}", e["mss"]), (f"r2_{i}", e["r2"]), (f"rmse_{i}", e["rmse"])]
    Sfin = sigma if sigma is not None else (E.T @ E / n)
    # the log likelihood with the covariance of the final residuals (compat 1923)
    ll = -0.5 * n * (m * (1 + np.log(2 * np.pi)) + np.linalg.slogdet(E.T @ E / n)[1])
    scalars += [("ll", ll), ("n_eq", float(m)), ("k", float(k)), ("k_eq_model", 0.0), ("converge", 1.0),
                ("rank", float(np.linalg.matrix_rank(V))), ("k_eq", float(m))]
    covs = " ".join(ses[0].variables()) if funcprog else " ".join(dict.fromkeys(
        v for se in ses for v in re.findall(r"[A-Za-z_]\w*", _BRACE.sub(" ", se.text)) if ds.has(v) and v not in deps))
    macros = [("cmdline", cmdline("nlsur", orig_args)), ("cmd", "nlsur"), ("predict", "nlsur_p"),
              ("estat_cmd", "nlsur_estat"), ("covariates", covs), ("rhs", covs)]
    if funcprog:
        macros.append(("funcprog", funcprog))
    macros += [("vce", "gnr"), ("method", method),
               ("constants", " ".join(e["const"] if e["const"] not in ("", "(none)") else "." for e in info)),
               ("type", "3" if funcprog else ("2" if pm else "1")), ("params", " ".join(names)),
               ("depvar", " ".join(deps))]
    macros += [(f"depvar_{i}", deps[i - 1]) for i in range(m, 0, -1)]
    macros.append(("properties", "b V"))
    from ..commands.matrix import Matrix
    nlsur_mats = [("init", Matrix(theta0.reshape(1, -1), ["y1"], list(names))),
                  ("Sigma", Matrix(Sfin, list(deps), list(deps)))]
    post(s, est, scalars, macros, mask, extra_mats=nlsur_mats)
    _display_nlsur(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _display_nlsur(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    out = s.output
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        title = {"fgnls": "FGNLS regression", "ifgnls": "Iterated FGNLS regression",
                 "nls": "Nonlinear least-squares regression"}[x["method"]]   # compat 1923
        out.write(f"\n{title} \n", "text")
        out.write("-" * 71 + "\n", "text")
        out.write(f"{'Equation':>15} |{'Obs':>11}{'Parms':>8}{'RMSE':>11}{'R-sq':>10}{'Constant':>13}\n", "text")
        out.write("-" * 16 + "+" + "-" * 54 + "\n", "text")
        for i, (dn, e) in enumerate(zip(x["deps"], x["info"]), 1):
            # R-sq is written from column 4 of its field, so a negative value shifts the rest
            r2 = g(e['r2'], '%9.4f').strip() + ("*" if e.get("uncentered") else "")
            out.write(f"{i:>2}{dn:>13} |{e['obs']:>11}{e['parms']:>8}{g(e['rmse'], '%9.0g'):>11}"
                      f"    {r2:<7}{e['const'] or '.':>12}\n", "text")
        out.write("-" * 71 + "\n", "text")
        if any(e.get("uncentered") for e in x["info"]):
            out.write("* Uncentered R-sq\n", "text")
        out.write("\n", "text")
    if table:
        coef_table(s, "", est.rows, stat="z", df=None, level=lev)
