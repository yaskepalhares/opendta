"""glm ([R] glm): modelos lineares generalizados por máxima verossimilhança
(Newton-Raphson, erros-padrão OIM), famílias gaussian, binomial, poisson,
gamma, igaussian e nbinomial; ligações identity, log, logit, probit,
cloglog, loglog, logc, power #, opower #."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import special, stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from .models2 import _Prep, _bmap, _cluster_note, _expand, _finish, _log, _semap, num_hess
from .results import CoefRow, Estimates, coef_table, rows_from_columns

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def g(x: float, fmt: str) -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


# ---------------------------------------------------------------------------
# ligações: (g(μ), μ = g⁻¹(η), dμ/dη, texto, nome)
# ---------------------------------------------------------------------------

def _link(name: str, power: float | None = None):
    n = name.lower()
    if n in ("i", "identity", "id", "ide", "iden"):
        return (lambda m: m, lambda e: e, lambda e: np.ones_like(e), "g(u) = u", "Identity")
    if n in ("log",):
        return (np.log, lambda e: np.exp(np.clip(e, -700, 700)), lambda e: np.exp(np.clip(e, -700, 700)),
                "g(u) = ln(u)", "Log")
    if n in ("logit", "l"):
        return (special.logit, special.expit, lambda e: special.expit(e) * (1 - special.expit(e)),
                "g(u) = ln(u/(1-u))", "Logit")
    if n in ("probit", "p"):
        return (st.norm.ppf, st.norm.cdf, st.norm.pdf, "g(u) = invnorm(u)", "Probit")
    if n in ("cloglog", "c"):
        return (lambda m: np.log(-np.log(1 - m)), lambda e: 1 - np.exp(-np.exp(e)),
                lambda e: np.exp(e - np.exp(e)), "g(u) = ln(-ln(1-u))", "Complementary log-log")
    if n in ("loglog", "ll"):
        return (lambda m: -np.log(-np.log(m)), lambda e: np.exp(-np.exp(-e)),
                lambda e: np.exp(-e - np.exp(-e)), "g(u) = -ln(-ln(u))", "Log-log")
    if n in ("logc",):
        return (lambda m: np.log(1 - m), lambda e: 1 - np.exp(e), lambda e: -np.exp(e),
                "g(u) = ln(1-u)", "Log-complement")
    if n in ("power", "pow"):
        a = 1.0 if power is None else power
        if a == 0:
            return _link("log")
        return (lambda m: m ** a, lambda e: np.sign(e) * np.abs(e) ** (1 / a),
                lambda e: (1 / a) * np.abs(e) ** (1 / a - 1), f"g(u) = u^({a:g})", f"Power({a:g})")
    raise StataError(198, f"unknown link {name}")   # VERIFICAR


_FAMILY_DEFAULT_LINK = {"gaussian": ("identity", None), "binomial": ("logit", None), "poisson": ("log", None),
                        "gamma": ("power", -1.0), "igaussian": ("power", -2.0), "nbinomial": ("log", None)}

_FAMILY_ALIASES = {"gau": "gaussian", "gaussian": "gaussian", "normal": "gaussian", "b": "binomial",
                   "bin": "binomial", "binomial": "binomial", "p": "poisson", "poi": "poisson",
                   "poisson": "poisson", "gam": "gamma", "gamma": "gamma", "ig": "igaussian",
                   "igaussian": "igaussian", "inv": "igaussian", "nb": "nbinomial", "nbinomial": "nbinomial"}


def cmd_glm(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "glm", args)
        return
    pr = _Prep(s, args, {"family": 1, "link": 1, "eform": 2, "irls": 4, "offset": 3, "exposure": 3,
                         "scale": 3})
    o = pr.o
    out = s.output
    pr.write_notes(out)
    fam_words = str(o.get("family") or "gaussian").split()
    fam = _FAMILY_ALIASES.get(fam_words[0].lower())
    if fam is None:
        raise StataError(198, f"unknown family {fam_words[0]}")   # VERIFICAR
    ds = s.data
    y = pr.y
    m_trials = np.ones(len(y))
    if fam == "binomial" and len(fam_words) > 1:
        tok = fam_words[1]
        try:
            m_trials = np.full(len(y), float(tok))
        except ValueError:
            m_trials = ds.get(tok).data[pr.smp.mask].astype(float)
    nb_k = float(fam_words[1]) if fam == "nbinomial" and len(fam_words) > 1 else 1.0
    lk_words = str(o.get("link") or "").split()
    if lk_words:
        lname = lk_words[0]
        lpow = float(lk_words[1]) if len(lk_words) > 1 else None
    else:
        lname, lpow = _FAMILY_DEFAULT_LINK[fam]
    gfun, ginv, dmu, link_text, link_name = _link(lname, lpow)
    off = np.zeros(len(y))
    if o.get("exposure"):
        off = np.log(ds.get(str(o["exposure"]).strip()).data[pr.smp.mask].astype(float))
    elif o.get("offset"):
        off = ds.get(str(o["offset"]).strip()).data[pr.smp.mask].astype(float)
    X = np.column_stack([pr.X, np.ones(len(y))]) if pr.constant else pr.X
    k = X.shape[1]
    w = pr.w
    if fam == "binomial" and (np.any(y < 0) or np.any(y > m_trials)):
        raise StataError(499, f"{pr.smp.depname} has values outside [0, n]")   # VERIFICAR

    def var(mu):
        if fam == "gaussian":
            return np.ones_like(mu)
        if fam == "binomial":
            return mu * (1 - mu / m_trials)
        if fam == "poisson":
            return mu
        if fam == "gamma":
            return mu ** 2
        if fam == "igaussian":
            return mu ** 3
        return mu + nb_k * mu ** 2

    def deviance_i(mu):
        with np.errstate(divide="ignore", invalid="ignore"):
            if fam == "gaussian":
                return (y - mu) ** 2
            if fam == "binomial":
                t1 = np.where(y > 0, y * np.log(y / mu), 0.0)
                t2 = np.where(m_trials - y > 0, (m_trials - y) * np.log((m_trials - y) / (m_trials - mu)), 0.0)
                return 2 * (t1 + t2)
            if fam == "poisson":
                return 2 * (np.where(y > 0, y * np.log(y / mu), 0.0) - (y - mu))
            if fam == "gamma":
                return -2 * (np.log(y / mu) - (y - mu) / mu)
            if fam == "igaussian":
                return (y - mu) ** 2 / (mu ** 2 * y)
            return 2 * (np.where(y > 0, y * np.log(y / mu), 0.0)
                        - (1 + nb_k * y) / nb_k * np.log((1 + nb_k * y) / (1 + nb_k * mu)))

    def mu_of(b):
        return ginv(X @ b + off)

    def loglik(mu, phi):
        with np.errstate(divide="ignore", invalid="ignore"):
            if fam == "gaussian":
                return -0.5 * ((y - mu) ** 2 / phi + np.log(2 * np.pi * phi))
            if fam == "binomial":
                p = np.clip(mu / m_trials, 1e-300, 1 - 1e-16)
                return (special.gammaln(m_trials + 1) - special.gammaln(y + 1) - special.gammaln(m_trials - y + 1)
                        + y * np.log(p) + (m_trials - y) * np.log1p(-p))
            if fam == "poisson":
                return -mu + y * np.log(mu) - special.gammaln(y + 1)
            if fam == "gamma":
                a = 1 / phi
                return a * np.log(a * y / mu) - a * y / mu - np.log(y) - special.gammaln(a)
            if fam == "igaussian":
                return -0.5 * ((y - mu) ** 2 / (y * mu ** 2 * phi) + np.log(2 * np.pi * phi * y ** 3))
            mm = 1 / nb_k
            return (special.gammaln(y + mm) - special.gammaln(mm) - special.gammaln(y + 1)
                    + mm * np.log(mm / (mm + mu)) + y * np.log(mu / (mm + mu)))

    # quase-verossimilhança para o NR (dispersão 1): soma de (y-μ)/V dμ/dη x
    def grad(b):
        eta = X @ b + off
        mu = ginv(eta)
        return X.T @ (w * (y - mu) / var(mu) * dmu(eta))

    def fun(b, todo):
        mu = mu_of(b)
        if fam == "binomial" and (np.any(mu <= 0) or np.any(mu >= m_trials)):
            return -np.inf, None, None
        if fam in ("poisson", "gamma", "igaussian", "nbinomial") and np.any(mu <= 0):
            return -np.inf, None, None
        # NR na verossimilhança com dispersão fixa em 1 (gaussian: equivale a MQO)
        ll = float(w @ loglik(mu, 1.0)) if fam not in ("gaussian",) else float(-0.5 * w @ (y - mu) ** 2)
        if fam in ("gamma", "igaussian"):
            ll = float(-0.5 * w @ deviance_i(mu))
        if todo == 0:
            return ll, None, None
        return ll, grad(b), num_hess(grad, b)

    # início: IRLS com μ = (y + ȳ)/2 (VERIFICAR: valores iniciais do Stata)
    ybar = float(np.sum(w * y) / np.sum(w))
    mu0 = (y + ybar) / 2 if fam != "binomial" else (y + 0.5 * m_trials) / 2
    eta0 = gfun(mu0)
    z = eta0 - off
    sw = w * dmu(eta0) ** 2 / var(mu0)
    b0 = np.linalg.lstsq(X * np.sqrt(sw)[:, None], z * np.sqrt(sw), rcond=None)[0]
    res = pr.maximize(fun, b0)
    mu = mu_of(res.b)
    N = pr.N
    df = N - k
    dev = float(w @ deviance_i(mu))
    pearson = float(w @ ((y - mu) ** 2 / var(mu)))
    if fam in ("binomial", "poisson", "nbinomial"):
        phi = 1.0
    else:
        phi = pearson / df
    ll_rep = float(w @ loglik(mu, phi if fam != "gaussian" else dev / N))
    # log de iterações: nas famílias sem escala, o objetivo do NR é a própria
    # verossimilhança. VERIFICAR o log do Stata para gaussian, gamma e igaussian
    if fam in ("binomial", "poisson", "nbinomial"):
        _log(s, pr, res)
    H = num_hess(grad, res.b)
    V, n_clust = pr.vce(H, None)
    V = V * phi
    if pr.robust:
        eta = X @ res.b + off
        sc = X * ((y - mu) / var(mu) * dmu(eta))[:, None]
        V2, n_clust = pr.vce(H, sc)
        V = V2
    aic = (-2 * ll_rep + 2 * k) / N
    bic = dev - df * np.log(N)
    names = []
    for c in pr.smp.cols:
        from . import fvars
        names.append(c.name if not c.omitted else fvars.omitted_name(c.name))
    rows = rows_from_columns(pr.smp.cols, _bmap(pr, res.b, k), _semap(pr, V, k), ds, constant=pr.constant)
    layout = [(nm, cc.values is not None) for nm, cc in zip(names, pr.smp.cols)]
    if pr.constant:
        layout.append(("_cons", True))
        names = names + ["_cons"]
    b_full, V_full = _expand(pr, res.b, V, layout)
    fam_text = {"gaussian": ("V(u) = 1", "Gaussian"), "binomial": ("V(u) = u*(1-u)", "Bernoulli")
                if np.all(m_trials == 1) else ("V(u) = u*(1-u/N)", "Binomial"),
                "poisson": ("V(u) = u", "Poisson"), "gamma": ("V(u) = u^2", "Gamma"),
                "igaussian": ("V(u) = u^3", "Inverse Gaussian"),
                "nbinomial": (f"V(u) = u+({nb_k:g})u^2", "Neg. Binomial")}[fam]
    est = _finish(s, pr, "glm", args, title="Generalized linear models", names=names,
                  eqnames=[pr.smp.depname] * len(names), b=b_full, V=V_full, rows=rows, res=res, ll0=SYS,
                  df_m=k - (1 if pr.constant else 0), chi2kind="Wald", chi2=SYS, p=SYS, r2p=None,
                  display=display_glm, n_clust=n_clust,
                  extra_scalars=[("phi", phi), ("aic", aic), ("bic", bic), ("deviance", dev),
                                 ("deviance_p", pearson), ("dispers", dev / df), ("dispers_p", pearson / df),
                                 ("df", df)],
                  extra_macros=[("varfunc", fam_text[0]), ("varfunct", fam_text[1]),
                                ("link", link_text), ("linkt", link_name)])
    s.e["ll"] = ll_rep
    est.extra.update({"ll": ll_rep, "dev": dev, "pearson": pearson, "df": df, "phi": phi, "aic": aic,
                      "bic": bic, "fam_text": fam_text, "link_text": link_text, "link_name": link_name})
    est.glm_ginv, est.glm_fam = ginv, fam
    est.predict = _predict_glm
    if o.get("eform"):
        est.coef_title, est.eform_on = "exp(b)", True   # VERIFICAR título por família/ligação
    display_glm(s, est)
    est.coef_title, est.eform_on = "Coef.", False


def display_glm(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        L = 50   # VERIFICAR colunas
        rows = [
            (f"{est.title}", ("No. of obs", g(x["N"], "%10.0fc") if False else f"{int(x['N']):,}")),
            ("Optimization     : ML", ("Residual df", f"{int(x['df']):,}")),
            ("", ("Scale parameter", g(x["phi"], "%10.0g"))),
            (f"Deviance         = {g(x['dev'], '%12.0g'):>12}", ("(1/df) Deviance", g(x["dev"] / x["df"], "%10.0g"))),
            (f"Pearson          = {g(x['pearson'], '%12.0g'):>12}", ("(1/df) Pearson", g(x["pearson"] / x["df"], "%10.0g"))),
        ]
        for lft, (rl, rv) in rows:
            out.write(f"{lft:<{L}}{rl:<16}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
        vf, vn = x["fam_text"]
        out.write(f"{'Variance function: ' + vf:<{L}}[{vn}]\n", "text")
        out.write(f"{'Link function    : ' + x['link_text']:<{L}}[{x['link_name']}]\n", "text")
        out.write("\n", "text")
        out.write(f"{'':<{L}}{'AIC':<16}= ", "text")
        out.write(f"{g(x['aic'], '%10.0g'):>10}\n", "result")
        lab = "Log pseudolikelihood" if est.robust else "Log likelihood"
        out.write(f"{(lab + '   = ' + g(x['ll'], '%12.0g')) if not est.robust else (lab + ' = ' + g(x['ll'], '%12.0g')):<{L}}"
                  f"{'BIC':<16}= ", "text")
        out.write(f"{g(x['bic'], '%10.0g'):>10}\n", "result")
        out.write("\n", "text")
    if table:
        _cluster_note(s, est)
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level,
                   vcetype=est.vcetype or "OIM", coef_title=getattr(est, "coef_title", "Coef."),
                   eform=getattr(est, "eform_on", False))


def _predict_glm(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat in ("", "mu"):
        note = "(option mu assumed; predicted mean " + est.depvar + ")" if not stat else ""
        return np.where(ok, est.glm_ginv(xb), SYS), note
    raise StataError(198, f"option {stat} not allowed")
