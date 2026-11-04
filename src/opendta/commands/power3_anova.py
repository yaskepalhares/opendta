"""power twoway and power repeated ([PSS] power twoway, power repeated): two-way and repeated-measures
ANOVA designs. The formulas are in core/power_suites.py."""

from __future__ import annotations

import math
import re

import numpy as np

from ..core import power_suites as C
from ..core.errors import StataError
from . import power as P
from .matrix import Matrix, get_matrix
from .power3_base import _I, _M, _MAT, _SOLVER, _T, _fail, _grid, _solver_vals, _store


class TwoWay(_M):
    method = "twoway"
    allowed = {"npercell", "cellweights", "nrows", "ncols", "factor", "vareffect", "varrow",
               "varcolumn", "varrowcolumn", "varerror", "showcellsizes"}
    _NAMES = {"row": "row", "column": "column", "rowcol": "row-by-column"}
    _VAR = {"row": "Var_r", "column": "Var_c", "rowcol": "Var_rc"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        b = body.strip()
        self.means = None
        if b:
            if re.fullmatch(r"[A-Za-z_]\w*", b):
                self.means = np.array(get_matrix(s, b).data, dtype=float)
            else:
                self.means = np.array(_grid(b), dtype=float)
            if min(self.means.shape) < 2:
                _fail("at least two rows and two columns of means are required")
            self.J, self.K = self.means.shape
        else:
            nr, nc = self.one("nrows"), self.one("ncols")
            if nr is None or nc is None:
                _fail("nrows() and ncols() are required without cell means")
            self.J, self.K = int(nr), int(nc)
            if min(self.J, self.K) < 2:
                _fail("nrows() and ncols() must be at least 2")
        fac = o.get("factor")
        which = None
        if fac not in (None, True):
            f = str(fac).strip().lower()
            which = ("rowcol" if f.startswith("rowc") else "row" if f in ("r", "ro", "row")
                     else "column" if f and "column".startswith(f) else None)
            if which is None:
                _fail("option factor() incorrectly specified")
        named = [n for n in ("varrow", "varcolumn", "varrowcolumn") if o.get(n)]
        if len(named) > 1 or (named and which is not None):
            _fail("only one of factor(), varrow(), varcolumn() or varrowcolumn() may be specified")
        if named:
            which = {"varrow": "row", "varcolumn": "column", "varrowcolumn": "rowcol"}[named[0]]
            self.veff = self.vals(named[0])
        else:
            self.veff = self.vals("vareffect")
        self.which = which or "row"
        self.cw = None
        if o.get("cellweights") not in (None, True):
            t = str(o["cellweights"]).strip()
            arr = (np.array(get_matrix(s, t).data, dtype=float) if re.fullmatch(r"[A-Za-z_]\w*", t)
                   else np.array(_grid(t), dtype=float))
            if arr.shape != (self.J, self.K) or (arr <= 0).any():
                _fail("option cellweights() incorrectly specified")
            if not self.flag("nfractional") and not np.allclose(arr, np.round(arr)):
                _fail("cellweights() must contain integers unless nfractional is specified")
            self.cw = arr
            self.we = np.array(C.standard_weights(list(arr.reshape(-1)))).reshape(arr.shape)
        if self.cw is not None and o.get("npercell"):
            _fail("npercell() and cellweights() may not be combined")

    def base(self):
        b = self.common()
        v = self.vals("npercell")
        if v is not None:
            b["npercell"] = v
        if self.veff is not None:
            b["veff"] = self.veff
        b["var_e"] = self.vals("varerror") or [1.0]
        return b

    def title(self, what):
        if what == "n":
            return "Estimated sample size for two-way ANOVA"
        if what == "power":
            return "Estimated power for two-way ANOVA"
        return f"Estimated {self._NAMES[self.which]} variance for two-way ANOVA"

    def hyp_lines(self, c, what):
        return [f"F test for {self._NAMES[self.which]} effect", "Ho: delta = 0  versus  Ha: delta != 0"]

    def compute(self, c):
        self.fix(c)
        J, K, which, a = self.J, self.K, self.which, c["alpha"]
        ve = c["var_e"]
        if self.means is not None:
            veff = C.twoway_variances(self.means)[("row", "column", "rowcol").index(which)]
        else:
            veff = c.get("veff")
        have_n = c.get("n") is not None or c.get("npercell") is not None
        if self.means is None and veff is None:
            if not have_n:
                _fail("cell means or vareffect() required")
            what = "effect"
            if c.get("power") is None:
                _fail("power() is required to compute the effect size")          # VERIFY
        else:
            # VERIFY: vareffect() given together with cell means is ignored
            if have_n:
                if c.get("power") is not None:
                    _fail("only two of n(), power() and the effect size may be specified")
                what = "power"
            else:
                c.setdefault("power", 0.8)
                what = "n"
        sw_bal = np.full((J, K), 1.0)
        w = self.we if self.cw is not None else sw_bal
        frac = w / w.sum()

        def lam(N, fr):
            if self.means is not None:
                return C.twoway_noncentrality(self.means, fr, N, which, ve)
            return N * veff / ve

        if what == "power":
            N = float(c["n"]) if c.get("n") is not None else J * K * c["npercell"]
            fr = frac if self.cw is not None else sw_bal / (J * K)
            Ne = N
            if self.cw is None and not self.flag("nfractional"):       # whole subjects per cell
                Ne = J * K * math.floor(N / (J * K) + 1e-9)
            c["power_out"] = C.twoway_power(lam(Ne, fr), J, K, Ne, which, a)
            c["N"], c["delta2"], c["Ne"] = N, lam(Ne, fr) / Ne, Ne
        elif what == "n":
            sw = w.sum()

            def pw(k):
                return C.twoway_power(lam(k * sw, frac), J, K, k * sw, which, a)
            k = self.size(pw, math.floor(J * K / sw) + 1, c["power"])
            c["N"], c["delta2"] = k * sw, lam(k * sw, frac) / (k * sw)
            c["power_out"] = pw(k)
            c["_init"] = float(J * K)
            if self.exact is not None:
                c["N_a"] = self.exact * sw
        else:
            if self.cw is not None:
                _fail("option cellweights() is not available yet in OpenDTA for the effect size")
            N = float(c["n"]) if c.get("n") is not None else J * K * c["npercell"]
            c["N"] = N
            lm = C.lambda_for_power(C.twoway_df1(J, K, which), N - J * K, a, c["power"])
            c["delta2"] = lm / N
            c["power_out"] = c["power"]
            c["_init"] = 1.0
        c["veff"] = c["delta2"] * ve if (what == "effect" or self.cw is not None) else veff
        c["balanced"] = self.cw is None
        c["cells"] = (w / w.sum() * c["N"]) if self.cw is not None else None
        c["delta"] = math.sqrt(c["delta2"])
        return what

    def sizes(self, c) -> list:
        N = c["N"]
        if c["balanced"]:
            return [("N", _I(N)), ("N per cell", _I(c.get("Ne", N) / (self.J * self.K)))]
        return [("N", _I(N)), ("Average N", P._f4(N / (self.J * self.K))), ("Nij", _MAT())]

    def show(self, c, what):
        vn = self._VAR[self.which]
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "n":
            if what == "effect":
                items.append(("power", P._f4(c["power"])))
            items += self.sizes(c)
        else:
            items.append(("power", P._f4(c["power"])))
        if what != "effect":
            items.append(("delta", P._f4(c["delta"])))
        items += [("N_r", _I(self.J)), ("N_c", _I(self.K))]
        if self.means is not None:
            items.append(("means", _MAT()))
        if what != "effect":
            items.append((vn, P._f4(c["veff"])))
        items.append(("Var_e", P._f4(c["var_e"])))
        self.put(items)
        if what == "n":
            self.results("Estimated sample sizes:", self.sizes(c))
            if self.flag("showcellsizes") and not c["balanced"]:
                self.cellsizes(c)
        elif what == "power":
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results(f"Estimated effect size and {self._NAMES[self.which]} variance:",
                         [("delta", P._f4(c["delta"])), (vn, P._f4(c["veff"]))])
        self.returns(c, what)

    def cellsizes(self, c):                       # VERIFY: layout of the cell sizes
        out = self.s.output
        out.write("\nCell sample sizes\n", "text")
        out.write("columns\n" + " ".join(f"{j + 1:>6}" for j in range(self.K)) + "\nrows\n", "text")
        for i in range(self.J):
            out.write(f"{i + 1} " + " ".join(f"{P._nfmt(v):>6}" for v in c["cells"][i]) + "\n", "result")

    def columns(self, c, what):
        cols = ["alpha", "power", "N", "N_per_cell" if c["balanced"] else "N_avg", "delta", "N_r", "N_c"]
        return cols + [self._VAR[self.which], "Var_e"]

    def row(self, c, what):
        N = c["N"]
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": N, "N_per_cell": c.get("Ne", N) / (self.J * self.K), "N_avg": N / (self.J * self.K),
                "delta": c["delta"], "N_r": self.J, "N_c": self.K, self._VAR[self.which]: c["veff"],
                "Var_e": c["var_e"]}

    _ORDER = (["Var_e", "balanced", "N_per_cell", "N_avg", "N_rc", "N_r", "N_c", "N_a", "m1_1", "m1_2",
               "delta"] + _SOLVER + ["onesided", "nfractional", "alpha", "beta", "power", "N", "m2_3",
                                     "m2_2", "m2_1", "m1_3", "type", "method", "means", "Nij"])

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        N, JK = c["N"], self.J * self.K
        v = {"N": float(N), "N_avg": N / JK, "N_rc": float(JK), "N_a": c.get("N_a", float(N)),
             "power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": 0.0,
             "balanced": float(c["balanced"]), "N_r": float(self.J), "N_c": float(self.K),
             self._VAR[self.which]: float(c["veff"]), "Var_e": float(c["var_e"])}
        if c["balanced"]:
            v["N_per_cell"] = c.get("Ne", N) / JK
        v.update(_solver_vals(c))
        v.update({"type": "test", "method": "twoway"})
        cells = np.full((self.J, self.K), c.get("Ne", N) / JK) if c["balanced"] else c["cells"]
        v["Nij"] = Matrix(cells)
        if self.means is not None:
            v["means"] = Matrix(self.means)
            for i in range(self.J):
                for j in range(self.K):
                    v[f"m{i + 1}_{j + 1}"] = float(self.means[i, j])
        _store(self.s, v, [self._VAR[self.which]] + self._ORDER)


class Repeated(_M):
    method = "repeated"
    allowed = {"npergroup", "grweights", "ngroups", "nrepeated", "corr", "covmatrix", "factor",
               "vareffect", "varbetween", "varwithin", "varbwithin", "varerror"}
    _VARN = {"between": ("Var_b", "Var_be"), "within": ("Var_w", "Var_we"),
             "bwithin": ("Var_bw", "Var_bwe")}
    _TEST = {"between": "F test for between subjects",
             "within": "F test for within subject with Greenhouse-Geisser correction",
             "bwithin": "F test for between-within subjects with Greenhouse-Geisser correction"}
    _EFF = {"between": "between-subjects", "within": "within-subject",
            "bwithin": "between-within-subjects"}          # VERIFY: names in the effect-size title

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        b = body.strip()
        self.means = None
        self.cov = None
        if o.get("covmatrix") not in (None, True):
            self.cov = np.array(get_matrix(s, str(o["covmatrix"])).data, dtype=float)
        elif not o.get("corr"):
            _fail("corr() or covmatrix() is required")
        if b:
            self.means = (np.array(get_matrix(s, b).data, dtype=float)
                          if re.fullmatch(r"[A-Za-z_]\w*", b) else np.array(_grid(b), dtype=float))
            self.J, self.K = self.means.shape
            if self.K < 2:
                _fail("at least two repeated measurements are required")
        else:
            ng, nr = self.one("ngroups"), self.one("nrepeated")
            if ng is None:
                _fail("ngroups() is required without cell means")
            self.J = int(ng)
            self.K = int(nr) if nr is not None else (self.cov.shape[0] if self.cov is not None else 0)
            if self.K < 2:
                _fail("nrepeated() is required")
        if self.cov is not None and self.cov.shape != (self.K, self.K):
            _fail("covmatrix() does not conform to the number of repeated measurements")
        fac = o.get("factor")
        which = "between"
        if fac not in (None, True):
            f = str(fac).strip().lower()
            which = ("bwithin" if f.startswith("bw") else "between" if f and "between".startswith(f)
                     else "within" if f and "within".startswith(f) else None)
            if which is None:
                _fail("option factor() incorrectly specified")
        named = [n for n in ("varbetween", "varwithin", "varbwithin") if o.get(n)]
        if len(named) > 1 or (named and fac not in (None, True)):
            _fail("only one of factor(), varbetween(), varwithin() or varbwithin() may be specified")
        if named:
            which = {"varbetween": "between", "varwithin": "within", "varbwithin": "bwithin"}[named[0]]
            self.veff = self.vals(named[0])
        else:
            self.veff = self.vals("vareffect")
        self.which = which
        if self.cov is not None and o.get("varerror"):
            _fail("varerror() may not be combined with covmatrix()")
        self.gw = None
        if o.get("grweights") not in (None, True):
            w = P._values(str(o["grweights"]))
            if len(w) != self.J or min(w) <= 0:
                _fail("option grweights() incorrectly specified")
            self.gw = w
        if self.gw is not None and o.get("npergroup"):
            _fail("npergroup() and grweights() may not be combined")

    def base(self):
        b = self.common()
        v = self.vals("npergroup")
        if v is not None:
            b["npergroup"] = v
        if self.veff is not None:
            b["veff"] = self.veff
        if self.cov is None:
            b["corr"] = self.vals("corr")
            b["var_e"] = self.vals("varerror") or [1.0]
        return b

    def title(self, what):
        if what == "n":
            return "Estimated sample size for repeated-measures ANOVA"
        if what == "power":
            return "Estimated power for repeated-measures ANOVA"
        return f"Estimated {self._EFF[self.which]} variance for repeated-measures ANOVA"

    def hyp_lines(self, c, what):
        return [self._TEST[self.which], "Ho: delta = 0  versus  Ha: delta != 0"]

    def compute(self, c):
        self.fix(c)
        J, K, a = self.J, self.K, c["alpha"]
        cov = self.cov if self.cov is not None else C.cs_covariance(K, c["var_e"], c["corr"])
        ws = C.standard_weights(self.gw) if self.gw is not None else [1.0] * J
        wf = np.array(ws) / sum(ws)
        have_n = c.get("n") is not None or c.get("npergroup") is not None
        have_eff = self.means is not None or c.get("veff") is not None
        try:
            ve_arg = c.get("veff") if c.get("veff") is not None else (1.0 if self.means is None else None)
            sm = C.repeated_summary(J, K, cov, self.which, wf, self.means, ve_arg)
        except ValueError as e:
            _fail(str(e))
        if not have_eff:
            if not have_n:
                _fail("cell means or a variance option required")
            what = "effect"
            if c.get("power") is None:
                _fail("power() is required to compute the effect size")          # VERIFY
        elif have_n:
            if c.get("power") is not None:
                _fail("only two of n(), power() and the effect size may be specified")
            what = "power"
        else:
            c.setdefault("power", 0.8)
            what = "n"
        pw = lambda N, d2: C.repeated_power(d2, sm["dc"], sm["du"], J, N, sm["eps"], a)
        sw = sum(ws)
        if what == "power":
            N = float(c["n"]) if c.get("n") is not None else float(J * c["npergroup"])
            c["N"], c["power_out"], c["delta2"] = N, pw(N, sm["delta2"]), sm["delta2"]
            c["beta_out"] = C.repeated_power(sm["delta2"], sm["dc"], sm["du"], J, N, sm["eps"], a, lower=True)
        elif what == "n":
            k = self.size(lambda k: pw(k * sw, sm["delta2"]), math.floor(J / sw) + 1, c["power"])
            c["N"], c["delta2"] = k * sw, sm["delta2"]
            c["power_out"] = pw(k * sw, sm["delta2"])
            c["_init"] = float(J)
            if self.exact is not None:
                c["N_a"] = self.exact * sw
        else:
            N = float(c["n"]) if c.get("n") is not None else float(J * c["npergroup"])
            try:
                d2 = C.solve_up(lambda d: pw(N, d) - c["power"], 1e-9, 1.0)
            except ValueError:
                raise StataError(498, "failed to find a solution")
            c["N"], c["delta2"], c["power_out"] = N, d2, c["power"]
            c["_init"] = 1.0
        c["sm"], c["cov"] = sm, cov
        c["delta"] = math.sqrt(c["delta2"])
        c["var_x"] = c["delta2"] * sm["var_xe"] if what == "effect" else sm["var_x"]
        c["balanced"] = self.gw is None or len(set(ws)) == 1
        c["groups"] = [c["N"] * x for x in wf]
        return what

    def sizes(self, c):
        if c["balanced"]:
            return [("N", _I(c["N"])), ("N per group", _I(c["N"] / self.J))]
        return [("N", _I(c["N"])), ("Average N", P._f4(c["N"] / self.J))] + [
            (f"N{j + 1}", _I(g)) for j, g in enumerate(c["groups"])]

    def show(self, c, what):
        vx, vxe = self._VARN[self.which]
        items = [("alpha", P._f4(c["alpha"]))]
        if what == "effect":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items += self.sizes(c)
        else:
            items.append(("power", P._f4(c["power"])))
        if what != "effect":
            items.append(("delta", P._f4(c["delta"])))
        items += [("N_g", _I(self.J)), ("N_rep", _I(self.K))]
        if self.means is not None:
            items.append(("means", _MAT()))
        if what != "effect":
            items.append((vx, P._f4(c["var_x"])))
        items.append((vxe, P._f4(c["sm"]["var_xe"])))
        if self.cov is not None:
            items += [("Cov", _MAT()), ("spherical", _T("true" if c["sm"]["spherical"] else "false"))]
        else:
            items += [("Var_e", P._f4(c["var_e"])), ("rho", P._f4(c["corr"]))]
        self.put(items)
        if what == "n":
            self.results("Estimated sample sizes:", self.sizes(c))
        elif what == "power":
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results(f"Estimated effect size and {self._EFF[self.which]} variance:",
                         [("delta", P._f4(c["delta"])), (vx, P._f4(c["var_x"]))])
        self.returns(c, what)

    def columns(self, c, what):
        vx, vxe = self._VARN[self.which]
        cols = ["alpha", "power", "N", "N_per_group" if c["balanced"] else "N_avg", "delta", "N_g", "N_rep", vx]
        return cols + ([vxe] if self.cov is not None else [vxe, "rho"])

    def row(self, c, what):
        vx, vxe = self._VARN[self.which]
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": c["N"], "N_per_group": c["N"] / self.J, "N_avg": c["N"] / self.J,
                "delta": c["delta"], "N_g": self.J, "N_rep": self.K, vx: c["var_x"],
                vxe: c["sm"]["var_xe"], "rho": c.get("corr")}

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        vx, vxe = self._VARN[self.which]
        N = c["N"]
        v = {"N": float(N), "N_per_group": N / self.J, "N_avg": N / self.J, "N_a": c.get("N_a", float(N)),
             "power": pw, "beta": c.get("beta_out", 1 - pw), "alpha": c["alpha"],
             "delta": c["delta"], "nfractional": float(self.flag("nfractional")), "onesided": 0.0,
             "balanced": float(c["balanced"]), "N_g": float(self.J), "N_rep": float(self.K),
             vx: float(c["var_x"]), vxe: float(c["sm"]["var_xe"]), "epsilon": float(c["sm"]["eps"]),
             "epsilon_m": float(c["sm"]["eps"]),                  # VERIFY: epsilon_m when not spherical
             "spherical": float(c["sm"]["spherical"])}
        if self.cov is None:
            v.update({"Var_e": float(c["var_e"]), "corr": float(c["corr"])})
        v.update(_solver_vals(c))
        v.update({"factor": self.which, "type": "test", "method": "repeated"})
        v["Cov"] = Matrix(c["cov"])
        if self.means is not None:
            v["means"] = Matrix(self.means)
        order = (["balanced", "N_per_group", "N_avg", "N_g", "N_rep", "N_a", "N", "power", "beta", "alpha"]
                 + _SOLVER + ["nfractional", "onesided", "Var_e", "corr", vx, "delta", "epsilon_m",
                              "epsilon", "spherical", vxe, "factor", "type", "method", "Cov", "means"])
        _store(self.s, v, order)
