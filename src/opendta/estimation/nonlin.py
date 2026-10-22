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
from ..core.stats import chi2_tail
from ..core.varlist import resolve_name
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from . import fvars
from .mleq import aux_row
from .models2 import _header, g
from .regress import cmdline, comma, sig5
from .results import CoefRow, Estimates, coef_table, post

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

def cmd_nl(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..commands.summarize import weights as get_weights
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "nl", args)
        return
    eqs, tail = _split_equations(args)
    if len(eqs) != 1:
        raise StataError(198, "nl requires one equation: (depvar = expression)")   # VERIFY
    dep, rhs = _lhs_rhs(eqs[0])
    ds = s.data
    depname = resolve_name(ds, dep)
    se = SubExpr(s, rhs)
    if not se.names:
        raise StataError(198, "no parameters in the expression")   # VERIFY
    p = parse_standard(tail)
    o = match_options(p.options, {"initial": 4, "nolog": 4, "log": 3, "eps": 3, "iterate": 4, "level": 1,
                                  "variables": 3, "lnlsq": 3, "title": 5, "robust": 1, "vce": 3, "cluster": 2,
                                  "noconstant": 4, "hasconstant": 3, "nocnsreport": 5, "noheader": 6,
                                  "notable": 5, "delta": 3}) if p.options.strip() else {}
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
    est.title = str(o.get("title") or "")
    est.display = _display_nl

    def predict(s_, e, stat, mask_):
        t = dict(zip(e.sexp.names, e.b))
        v = e.sexp.evaluate(t)
        if stat in ("", "yhat"):
            return v, ""
        if stat == "residuals":
            return ds.get(e.depvar).data.astype(float) - v, ""
        raise StataError(198, f"option {stat} not allowed")   # VERIFY
    est.predict = predict
    est.predict_opts = {"yhat": 1, "residuals": 3}
    scalars = [("N", N), ("df_m", df_m), ("df_r", df_r), ("mss", mss), ("rss", rss), ("r2", r2),
               ("r2_a", r2_a), ("rmse", rmse), ("dev", dev), ("ic", iters), ("converged", 1.0)]
    macros = [("cmdline", cmdline("nl", args)), ("title", est.title or "Nonlinear regression"),
              ("depvar", depname), ("cmd", "nl"), ("properties", "b V"), ("predict", "nl_p"),
              ("expression", rhs)]
    post(s, est, scalars, macros, mask)
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
    eqs, tail = _split_equations(args)
    if not eqs:
        raise StataError(198, "equation required")   # VERIFY
    ds = s.data
    deps, ses = [], []
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
    o = match_options(p.options, {"fgnls": 2, "ifgnls": 2, "nls": 3, "initial": 4, "nolog": 4, "log": 3, "eps": 3,
                                  "iterate": 4, "level": 1, "variables": 3, "noconstant": 4, "title": 5,
                                  "noheader": 6, "notable": 5, "delta": 3, "robust": 1, "nequations": 3,
                                  "parameters": 5, "lnlsq": 3, "ifgnls": 3}) if p.options.strip() else {}
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
        E = (Y - fit_vals(theta))
        sigma = E.T @ E / n
        A = J.T @ J
        Ai = np.linalg.inv(A)
        # equation-by-equation residual variances
        Wm = np.kron(np.linalg.inv(sigma), np.eye(n))
        V = Ai @ (J.T @ np.kron(sigma, np.eye(n)) @ J) @ Ai
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
        info.append({"obs": n, "parms": kp, "rmse": float(np.sqrt(rss / n)), "r2": 1 - rss / tss if tss > 0 else SYS,
                     "rss": rss, "const": const})
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
    scalars = [("N", n), ("k", k), ("k_eq", m), ("rank", np.linalg.matrix_rank(V)), ("converged", 1.0)]
    macros = [("cmdline", cmdline("nlsur", args)), ("title", {"fgnls": "FGNLS regression", "ifgnls": "Iterated FGNLS regression",
                                                             "nls": "NLS regression"}[method]),
              ("depvar", " ".join(deps)), ("cmd", "nlsur"), ("properties", "b V"), ("method", method)]
    post(s, est, scalars, macros, mask)
    _display_nlsur(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _display_nlsur(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    out = s.output
    x = est.extra
    lev = level if level is not None else est.level
    if header:
        title = {"fgnls": "FGNLS regression", "ifgnls": "Iterated FGNLS regression", "nls": "NLS regression"}[x["method"]]
        out.write(f"\n{title} \n", "text")
        out.write("-" * 71 + "\n", "text")
        out.write(f"{'Equation':>15} |{'Obs':>11}{'Parms':>8}{'RMSE':>11}{'R-sq':>10}{'Constant':>13}\n", "text")
        out.write("-" * 16 + "+" + "-" * 54 + "\n", "text")
        for i, (dn, e) in enumerate(zip(x["deps"], x["info"]), 1):
            out.write(f"{i:>2}{dn:>13} |{e['obs']:>11}{e['parms']:>8}{g(e['rmse'], '%9.0g'):>11}"
                      f"{g(e['r2'], '%10.4f'):>10}{e['const'] or '.':>13}\n", "text")
        out.write("-" * 71 + "\n\n", "text")
    if table:
        coef_table(s, "", est.rows, stat="z", df=None, level=lev)


# ---------------------------------------------------------------------------
# gmm
# ---------------------------------------------------------------------------

def cmd_gmm(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..commands.summarize import weights as get_weights
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "gmm", args)
        return
    eqs, tail = _split_equations(args)
    if not eqs:
        raise StataError(198, "equation required")   # VERIFY
    ds = s.data
    ses = [SubExpr(s, body) for body in eqs]
    names: list[str] = []
    init: dict = {}
    for se in ses:
        for nm in se.names:
            if nm not in names:
                names.append(nm)
            if se.init.get(nm) is not None:
                init[nm] = se.init[nm]
    p = parse_standard(tail)
    # instruments() may repeat, so it is cut from the option text before matching the others
    opt_text = p.options
    inst_specs = re.findall(r"instruments\s*\(((?:[^()]|\([^()]*\))*)\)", opt_text)
    opt_text = re.sub(r"instruments\s*\(((?:[^()]|\([^()]*\))*)\)", " ", opt_text)
    o = match_options(opt_text, {"onestep": 3, "twostep": 3, "igmm": 4, "winitial": 4, "wmatrix": 2, "vce": 3,
                                 "level": 1, "nolog": 4, "log": 3, "eps": 3, "iterate": 4, "from": 3, "robust": 1,
                                 "cluster": 2, "noheader": 6, "notable": 5, "nocommonesample": 5, "quickderivatives": 5,
                                 "derivatives": 3, "variables": 3, "hasconstant": 3, "parameters": 5, "nequations": 3,
                                 "haslfderivatives": 5, "conv_maxiter": 5, "technique": 4, "nparameters": 5}) \
        if opt_text.strip() else {}
    if o.get("derivatives") or o.get("parameters") or o.get("nparameters") or o.get("igmm"):
        raise StataError(198, "that option of gmm is not available yet")   # VERIFY
    mask = touse(s, p)
    for se in ses:
        for v in se.variables():
            mask &= ds.get(v).data < SYS
    # instruments per equation (eqno: list), or one list for all
    inst_by_eq: dict[int, tuple[list[str], bool]] = {}
    for spec in inst_specs:
        mm = re.match(r"^\s*(?:(\d+)\s*:)?\s*(.*)$", spec, re.S)
        eqno = int(mm.group(1)) - 1 if mm.group(1) else None
        body = mm.group(2)
        noc = False
        if "," in body:
            body, opt = body.split(",", 1)
            noc = "noconstant".startswith(opt.strip()) and bool(opt.strip())
        vs = [resolve_name(ds, v) for v in body.split()]
        for q in range(len(ses)) if eqno is None else [eqno]:
            inst_by_eq[q] = (vs, noc)
    for q in range(len(ses)):
        if q not in inst_by_eq:
            raise StataError(198, "option instruments() required")   # VERIFY
    for vs, _ in inst_by_eq.values():
        for v in vs:
            mask &= ds.get(v).data < SYS
    theta0 = np.array([init.get(nm, 0.0) for nm in names])
    if o.get("from"):
        for mm in re.finditer(r"(\w+)\s+([-+0-9.eE]+)", str(o["from"])):
            if mm.group(1) in names:
                theta0[names.index(mm.group(1))] = float(mm.group(2))
    for se in ses:
        with np.errstate(all="ignore"):
            f0 = se.evaluate({nm: theta0[names.index(nm)] for nm in se.names})
        mask &= (f0 < SYS) & np.isfinite(f0)
    n = int(mask.sum())
    if n == 0:
        raise StataError(2000, "no observations")
    m = len(ses)
    Zs, znames = [], []
    for q in range(m):
        vs, noc = inst_by_eq[q]
        cols = [ds.get(v).data[mask].astype(float) for v in vs] + ([] if noc else [np.ones(n)])
        Zs.append(np.column_stack(cols))
        znames.append(" ".join(vs + ([] if noc else ["_cons"])))
    offs = np.concatenate([[0], np.cumsum([Z.shape[1] for Z in Zs])]).astype(int)
    n_mom = int(offs[-1])
    k = len(names)
    if n_mom < k:
        raise StataError(481, "equation not identified; number of moments is less than the number of parameters")   # VERIFY
    out = s.output
    wm = str(o.get("wmatrix", "robust")).split()
    wkind = wm[0].lower()[:3]
    cl = None
    clustvar = ""
    if wkind == "clu":
        clustvar = resolve_name(ds, wm[1])
        cv = ds.get(clustvar).data[mask]
        _, cl = np.unique(cv, return_inverse=True)
    if o.get("cluster"):
        clustvar = resolve_name(ds, str(o["cluster"]))
        _, cl = np.unique(ds.get(clustvar).data[mask], return_inverse=True)
        wkind = "clu"

    def resid(t):
        cols = []
        with np.errstate(all="ignore"):
            for se in ses:
                cols.append(se.evaluate({nm: t[names.index(nm)] for nm in se.names})[mask])
        return np.column_stack(cols)

    def moments(t):
        e = resid(t)
        return np.column_stack([Zs[q] * e[:, [q]] for q in range(m)])

    def mbar(t):
        return moments(t).mean(axis=0)

    def meat(t, kind):
        G_ = moments(t)
        if kind == "clu" and cl is not None:
            u = np.zeros((cl.max() + 1, G_.shape[1]))
            np.add.at(u, cl, G_)
            return u.T @ u / n
        return G_.T @ G_ / n

    def step(t, W, label):
        out.write(f"\n{label}\n" if label else "", "text")
        q_old = None
        for it in range(100):
            gb = mbar(t)
            G = _num_jac(mbar, t)
            q = float(gb @ W @ gb)
            if not o.get("nolog"):
                out.write(f"Iteration {it}:   GMM criterion Q(b) = {g(q, '%10.0g'):>10}\n", "text")
            if (q_old is not None and abs(q_old - q) < 1e-10 * (1 + abs(q))) or q < 1e-30:
                break
            d = np.linalg.solve(G.T @ W @ G + 1e-12 * np.eye(k), G.T @ W @ gb)
            lam, tn = 1.0, t
            while lam > 1e-8:
                tn = t - lam * d
                gn = mbar(tn)
                qn = float(gn @ W @ gn)
                if np.isfinite(qn) and qn <= q + 1e-12:
                    break
                lam /= 2
            q_old, t = q, tn
        return t

    # initial weight matrix: unadjusted, Z'Z/n blocks scaled by the residual covariance
    e0 = resid(theta0)
    S0 = np.zeros((n_mom, n_mom))
    sig0 = e0.T @ e0 / n
    for a in range(m):
        for b in range(m):
            S0[offs[a]:offs[a + 1], offs[b]:offs[b + 1]] = sig0[a, b] * (Zs[a].T @ Zs[b]) / n
    winit = str(o.get("winitial", "unadjusted")).strip().lower()
    if winit.startswith("iden"):
        W1 = np.eye(n_mom)
        winit_label = "Identity"
    else:
        W1 = np.linalg.inv(S0 + 1e-12 * np.eye(n_mom))
        winit_label = "Unadjusted"
    theta = step(theta0, W1, "Step 1")
    if o.get("onestep"):
        W = W1
        wm_label = winit_label
    else:
        S = meat(theta, wkind) if wkind in ("rob", "clu") else None
        if S is None:
            e1 = resid(theta)
            sg = e1.T @ e1 / n
            S = np.zeros((n_mom, n_mom))
            for a in range(m):
                for b in range(m):
                    S[offs[a]:offs[a + 1], offs[b]:offs[b + 1]] = sg[a, b] * (Zs[a].T @ Zs[b]) / n
        W = np.linalg.inv(S)
        theta = step(theta, W, "Step 2")
        wm_label = {"rob": "Robust", "una": "Unadjusted", "clu": f"Cluster ({clustvar})"}[wkind]
    gb = mbar(theta)
    G = _num_jac(mbar, theta)
    Q = float(gb @ W @ gb)
    Ainv = np.linalg.inv(G.T @ W @ G)
    vkind = "clu" if (cl is not None) else ("rob" if wkind == "rob" or str(o.get("vce", "")).startswith("rob")
                                              or o.get("robust") else "una")
    if vkind == "una":
        V = Ainv / n
    else:
        V = Ainv @ (G.T @ W @ meat(theta, vkind) @ W @ G) @ Ainv / n
    V = (V + V.T) / 2
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    est = Estimates("gmm", "", list(names), theta, V, n, "z", None)
    est.level = level
    est.vcetype = "Robust" if vkind in ("rob", "clu") else ""
    est.clustvar, est.n_clust = clustvar, (int(cl.max() + 1) if cl is not None else None)
    est.eqnames = ["/"] * k
    est.rows = [aux_row("/" + nm, float(theta[i]), float(np.sqrt(max(V[i, i], 0))), level) for i, nm in enumerate(names)]
    est.extra = dict(N=n, n_par=k, n_mom=n_mom, wm_label=wm_label, winit=winit_label, Q=Q, znames=znames,
                     onestep=bool(o.get("onestep")))
    est.sexps, est.pnames = ses, names

    def display(s_, e, *, header=True, table=True, level=None, **_kw):
        o_ = s_.output
        x = e.extra
        if header:
            o_.write("\nGMM estimation \n\n", "text")
            o_.write(f"Number of parameters = {x['n_par']:>3}\n", "text")
            o_.write(f"Number of moments    = {x['n_mom']:>3}\n", "text")
            o_.write(f"{'Initial weight matrix: ' + x['winit']:<50}{'Number of obs':<16}=", "text")
            o_.write(f"{comma(x['N']):>11}\n", "result")
            if not x["onestep"]:
                o_.write(f"GMM weight matrix:     {x['wm_label']}\n", "text")
            o_.write("\n", "text")
        if not table:
            return
        coef_table(s_, "", e.rows, stat="z", df=None, level=level or e.level, vcetype=e.vcetype)
        for i, zn in enumerate(x["znames"], 1):
            o_.write(f"Instruments for equation {i}: {zn}\n", "text")
    est.display = display
    scalars = [("N", n), ("k", k), ("k_eq", m), ("n_moments", n_mom), ("Q", Q), ("J", n * Q), ("J_df", n_mom - k),
               ("rank", np.linalg.matrix_rank(V)), ("converged", 1.0)]
    macros = [("cmdline", cmdline("gmm", args)), ("cmd", "gmm"), ("title", "GMM estimation"),
              ("wmatrix", wm_label.lower().split()[0]), ("estat_cmd", "gmm_estat"), ("properties", "b V")]
    post(s, est, scalars, macros, mask)
    display(s, est, header=not o.get("noheader"), table=not o.get("notable"))
