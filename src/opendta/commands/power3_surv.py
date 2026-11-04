"""power cox, power exponential and power logrank ([PSS] power cox, power exponential, power logrank):
survival designs. The formulas are in core/power_suites.py."""

from __future__ import annotations

import math
import re

import numpy as np

from ..core import power_suites as C
from ..core.errors import StataError
from . import power as P
from .matrix import get_matrix
from .power3_base import _I, _M, _SOLVER, _fail, _solver_vals, _store


class Cox(_M):
    method = "cox"
    iterates = False                       # closed forms: no iteration message (checked, case 1919)
    allowed = {"hratio", "sd", "r2", "eventprob", "failprob", "wdprob", "effect", "direction", "onesided"}
    _EFFN = {"coefficient": "coefficient", "hratio": "hazard ratio", "lnhratio": "log hazard-ratio"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        t = body.strip()
        self.b1 = P._values(t) if t else None
        if self.b1 is not None and o.get("hratio"):
            _fail("b1 and hratio() may not be combined")
        eff = str(o.get("effect") or "coefficient").lower()
        self.eff = next((e for e in self._EFFN if e.startswith(eff)), None)
        if self.eff is None:
            _fail("option effect() incorrectly specified")
        self.up_dir = str(o.get("direction") or "lower").lower().startswith("u")
        if self.vals("eventprob") and self.vals("failprob"):
            _fail("only one of eventprob() or failprob() may be specified")

    def base(self):
        b = self.common()
        if self.b1 is not None:
            b["b1"] = self.b1
        hr = self.vals("hratio")
        if hr is not None:
            b["hratio"] = hr
        for k, dflt in (("sd", [0.5]), ("r2", [0.0]), ("wdprob", [0.0])):
            b[k] = self.vals(k) or dflt
        b["pe"] = self.vals("eventprob") or self.vals("failprob") or [1.0]
        return b

    def title(self, what):
        w = {"n": "sample size", "power": "power", "effect": "target coefficient"}[what]
        return f"Estimated {w} for Cox PH regression"

    def hyp_lines(self, c, what):
        up = self.up_dir if what == "effect" else c["b1"] > 0
        extra = f"; beta1 {'>' if self.up_dir else '<'} 0" if (what == "effect" and self.sides() == 2) else ""
        return ["Wald test", f"Ho: beta1 = 0  versus  Ha: beta1 {self.hyp_op(c, up)} 0{extra}"]

    def compute(self, c):
        self.fix(c)
        a, sides = c["alpha"], self.sides()
        sd, r2, pe, wd = c["sd"], c["r2"], c["pe"], c["wdprob"]
        if sd <= 0 or not 0 <= r2 < 1 or not 0 < pe <= 1 or not 0 <= wd < 1:
            _fail("invalid argument value")                          # VERIFY
        b1 = c.get("b1")
        if b1 is None and c.get("hratio") is not None:
            b1 = math.log(c["hratio"])
        have_n = c.get("n") is not None
        if b1 is None:
            if have_n:
                what = "effect"
                c.setdefault("power", 0.8)
            else:
                b1 = math.log(0.5)                 # hratio() defaults to 0.5
                what = "n"
        if b1 is not None:
            if have_n:
                if c.get("power") is not None:
                    _fail("only two of n(), power() and the effect size may be specified")
                what = "power"
            else:
                what = "n"
        if what == "n":
            c.setdefault("power", 0.8)
            ev = C.cox_events(b1, sd, r2, a, sides, c["power"])
            n = ev / pe                          # N from the events before rounding (Stata 14: 164 with wdprob(.2))
            if wd > 0:
                n = n / (1 - wd)
            if not self.flag("nfractional"):
                n = math.ceil(n - 1e-9)
            ev = n * (1 - wd) * pe                  # events shown: from the final N (Stata 14)
            c["E"], c["N"] = float(ev if self.flag("nfractional") else math.ceil(ev - 1e-9)), float(n)
            c["power_out"] = C.cox_power(b1, sd, r2, n * (1 - wd), pe, a, sides)
        elif what == "power":
            c["N"] = float(c["n"])
            c["power_out"] = C.cox_power(b1, sd, r2, c["N"], pe, a, sides)
            c["E"] = c["N"] * pe if self.flag("nfractional") else math.ceil(c["N"] * pe - 1e-9)
        else:
            c["N"] = float(c["n"])
            mag = C.cox_effect(sd, r2, c["N"], pe, a, sides, c["power"])
            b1 = mag if self.up_dir else -mag
            c["power_out"] = c["power"]
            c["E"] = c["N"] * pe if self.flag("nfractional") else math.ceil(c["N"] * pe - 1e-9)
        c["b1"], c["hr"] = b1, math.exp(b1)
        c["delta"] = {"coefficient": b1, "hratio": math.exp(b1), "lnhratio": b1}[self.eff]
        return what

    def show(self, c, what):
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "power":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items.append(("N", _I(c["N"])))
        if what != "effect":
            items.append(("delta", P._f4(c["delta"]) + f"  ({self._EFFN[self.eff]})"))
            items.append(("hratio" if self.eff == "coefficient" else "b1",
                          P._f4(c["hr"] if self.eff == "coefficient" else c["b1"])))
        items.append(("sd", P._f4(c["sd"])))
        if self.o.get("r2"):
            items.append(("R2", P._f4(c["r2"])))
        self.put(items)
        cens = [("Pr_E", P._f4(c["pe"]))] + ([("Pr_w", P._f4(c["wdprob"]))] if self.o.get("wdprob") else [])
        if what == "n":
            self.results("Censoring and withdrawal:" if self.o.get("wdprob") else "Censoring:", cens)
            self.results("Estimated number of events and sample size:", [("E", _I(c["E"])), ("N", _I(c["N"]))])
        elif what == "power":
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results(f"Estimated effect size and target {self._EFFN[self.eff]}:",
                         [("delta", P._f4(c["delta"]) + f"  ({self._EFFN[self.eff]})"), ("b1", P._f4(c["b1"]))])
        self.returns(c, what)

    def columns(self, c, what):
        return ["alpha", "power", "N", "E", "delta", "hratio" if self.eff == "coefficient" else "b1", "sd"] + (
            ["R2"] if self.o.get("r2") else []) + ["Pr_E"]

    def row(self, c, what):
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": c["N"], "E": c["E"], "delta": c["delta"], "b1": c["b1"], "hratio": c["hr"],
                "sd": c["sd"], "R2": c["r2"], "Pr_E": c["pe"]}

    _ORDER = ["N", "power", "beta", "alpha", "delta", "nfractional", "onesided", "E", "hratio", "b1", "sd",
              "R2", "Pr_E", "Pr_w", "type", "method", "effect"]

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        v = {"N": c["N"], "power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": float(self.sides() == 1),
             "E": c["E"], "hratio": c["hr"], "b1": c["b1"], "sd": c["sd"], "R2": c["r2"], "Pr_E": c["pe"],
             "Pr_w": c["wdprob"], "type": "test", "method": "cox", "effect": self.eff}
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)


class Exponential(_M):
    method = "exponential"
    iterates = False
    allowed = {"time", "nratio", "hratio", "lnhratio", "hdifference", "difference", "diff",
               "loghazard", "unconditional", "onesided", "studytime", "fperiod", "aperiod",
               "aprob", "aptime", "atime", "ashape", "lossprob", "lossprob1", "lossprob2",
               "losstime", "losshaz", "losshaz1", "losshaz2", "effect"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        if "effect" in o:
            _fail("option effect() is not available yet in OpenDTA")
        toks = P._split_args(body)
        if len(toks) > 2:
            _fail("too many values specified")
        self.args = [P._values(t) for t in toks]
        self.surv = bool(o.get("time"))
        self.logtest = self.flag("loghazard") or not self.args
        self.unc = self.flag("unconditional")
        # accrual and follow-up (scalars)
        T, f, r = self.one("studytime"), self.one("fperiod"), self.one("aperiod")
        if r is not None and f is None and T is None:
            _fail("two of studytime(), fperiod() and aperiod() are required")    # VERIFY
        self.r = r or 0.0
        if T is not None and f is not None:
            self.r = T - f if r is None else r
            if abs(self.r + f - T) > 1e-9:
                _fail("studytime(), fperiod() and aperiod() are inconsistent")
        elif T is not None:
            f = T - self.r
        elif f is not None:
            T = self.r + f
        self.T, self.f = T, f
        if self.r < 0 or (T is not None and T <= 0):
            _fail("invalid accrual or follow-up period")
        self.gamma = self.one("ashape", 0.0)
        self.pa = self.ta = None
        if any(o.get(k) for k in ("aprob", "aptime", "atime")):
            if self.r <= 0 or o.get("ashape"):
                _fail("aprob(), aptime() and atime() require aperiod() and may not be combined with ashape()")
            self.pa = self.one("aprob", 0.5)
            self.ta = self.one("atime") if o.get("atime") else self.one("aptime", 0.5) * self.r
            if not 0 < self.ta <= self.r or not 0 < self.pa < 1:
                _fail("option atime(), aptime() or aprob() incorrectly specified")
            self.gamma = C.exp_gamma(self.r, self.ta, self.pa)
        # losses
        self.eta1 = self.eta2 = 0.0
        self.loss = None
        tl = self.one("losstime", 1.0)
        lp, lp1, lp2 = self.one("lossprob"), self.one("lossprob1"), self.one("lossprob2")
        lh, lh1, lh2 = self.one("losshaz"), self.one("losshaz1"), self.one("losshaz2")
        if (lp is not None or lp1 is not None or lp2 is not None) and (lh is not None or lh1 is not None
                                                                         or lh2 is not None):
            _fail("lossprob() and losshaz() may not be combined")
        if any(v is not None for v in (lp, lp1, lp2, lh, lh1, lh2)):
            if T is None:
                _fail("loss options require aperiod() or fperiod()")
            def eta_of(lpg, lpc, lhg, lhc):
                if lpg is not None or lpc is not None:
                    return -math.log(1 - (lpg if lpg is not None else lpc)) / tl
                return lhg if lhg is not None else (lhc or 0.0)
            self.eta1 = eta_of(lp1, lp, lh1, lh)
            self.eta2 = eta_of(lp2, lp, lh2, lh)
            if any(v is not None for v in (lp, lp1, lp2)):
                self.loss = [("lpr1", 1 - math.exp(-self.eta1 * tl)), ("lpr2", 1 - math.exp(-self.eta2 * tl)),
                             ("lt", tl)]
            else:
                self.loss = [("lh1", self.eta1), ("lh2", self.eta2)]                # VERIFY: names with losshaz()
        self.up_dir = None

    def base(self):
        b = self.common()
        names = ("s1", "s2") if self.surv else ("h1", "h2")
        for nm, v in zip(names, self.args):
            b[nm] = v
        for k in ("hratio", "lnhratio", "hdifference", "time"):
            v = self.vals(k)
            if v is not None:
                b[k] = v
        d = self.vals("diff") or self.vals("difference")
        if d is not None:
            b["hdifference"] = d
        b["nratio"] = self.vals("nratio") or [1.0]
        return b

    def title(self, what):
        return ("Estimated sample sizes for two-sample comparison of survivor functions" if what == "n"
                else "Estimated power for two-sample comparison of survivor functions")

    def hyp_lines(self, c, what):
        lt = c["logtest"]
        up = c["hr"] > 1
        if lt:
            return [c["test"], f"Ho: ln(HR) = 0  versus  Ha: ln(HR) {self.hyp_op(c, up)} 0"]
        return [c["test"], f"Ho: h2 = h1  versus  Ha: h2 {self.hyp_op(c, up)} h1"]

    def compute(self, c):
        self.fix(c)
        a, sides, R = c["alpha"], self.sides(), c["nratio"]
        if R <= 0:
            _fail("nratio() must be positive")
        l1 = l2 = None
        if self.surv:
            t = c.get("time")
            if t is None or t <= 0:
                _fail("time() is required with survival probabilities")
            if c.get("s1") is not None:
                l1 = -math.log(c["s1"]) / t
            if c.get("s2") is not None:
                l2 = -math.log(c["s2"]) / t
        else:
            l1, l2 = c.get("h1"), c.get("h2")
        eff = [k for k in ("hratio", "lnhratio", "hdifference") if c.get(k) is not None]
        if l1 is not None and l2 is not None:
            if eff:
                _fail("hratio(), lnhratio() and hdifference() may not be combined with both groups given")
        elif len(eff) > 1:
            _fail("only one of hratio(), lnhratio() and hdifference() may be specified")
        elif l1 is not None:
            if eff == ["hdifference"]:
                l2 = l1 + c["hdifference"]
            elif eff == ["lnhratio"]:
                l2 = l1 * math.exp(c["lnhratio"])
            else:
                l2 = l1 * (c["hratio"] if eff else 0.5)
        else:
            if eff == ["hdifference"]:
                _fail("hdifference() requires the control-group hazard or survival")
            hr0 = c["hratio"] if eff == ["hratio"] else (math.exp(c["lnhratio"]) if eff else 0.5)
            if self.T is not None:
                _fail("the control-group hazard or survival is required with accrual, follow-up or losses")
            l1, l2 = 1.0, hr0                   # VERIFY: without censoring only the ratio matters
            c["hz_unknown"] = True
        if min(l1, l2) <= 0:
            _fail("hazards must be positive")
        hr = l2 / l1
        logtest = self.logtest
        c.update({"l1": l1, "l2": l2, "hr": hr, "logtest": logtest})
        delta = math.log(hr) if logtest else l2 - l1
        c["delta"] = delta
        c["test"] = ("Exponential test, log hazard-ratio, " if logtest else "Exponential test, hazard difference, ") + \
            ("unconditional" if self.unc else "conditional")
        xo, xa = C.exp_xi(l1, l2, self.eta1, self.eta2, R, self.T, self.r, self.gamma, logtest, self.unc)
        have_n = c.get("n") is not None
        if have_n:
            if c.get("power") is not None:
                _fail("only two of n(), power() and the effect size may be specified")
            c["power_out"] = C.exp_power(delta, xo, xa, float(c["n"]), a, sides, self.unc)
            N = float(c["n"])
            n1 = N / (1 + R)
            n2 = N - n1
            what = "power"
        else:
            c.setdefault("power", 0.8)
            Nc = C.exp_size(delta, xo, xa, a, sides, c["power"], self.unc)
            if self.flag("nfractional"):
                n1 = Nc / (1 + R)
                n2 = R * n1
            else:
                n1 = float(math.ceil(Nc / (1 + R) - 1e-9))        # VERIFY: rounding of N1 and N2
                n2 = float(math.ceil(R * n1 - 1e-9))
            N = n1 + n2
            c["power_out"] = C.exp_power(delta, xo, xa, N, a, sides, self.unc)
            what = "n"
        c["N"], c["N1"], c["N2"] = N, n1, n2
        c["R_a"] = n2 / n1
        return what

    def show(self, c, what):
        out = self.s.output
        items = [("alpha", P._f4(c["alpha"]))]
        if what == "n":
            items.append(("power", P._f4(c["power"])))
        else:
            items += [("N", _I(c["N"])), ("N per group", _I(c["N1"])) if c["N1"] == c["N2"] else ("N1", _I(c["N1"]))]
            if c["N1"] != c["N2"]:
                items.append(("N2", _I(c["N2"])))
        lab = "log hazard-ratio" if c["logtest"] else "hazard difference"
        items.append(("delta", P._f4(c["delta"]) + f"  ({lab})"))
        if c["nratio"] != 1:
            items.append(("N2/N1", P._f4(c["nratio"])))
        self.put(items)
        if self.T is not None:
            out.write("\nAccrual and follow-up information:\n\n", "text")
            self.line("duration", P._f4(self.T))
            self.line("follow-up", P._f4(self.f))
            self.line("accrual", P._f4(self.r) + ("  (uniform)" if abs(self.gamma) < 1e-6 else "  (exponential)"))
            if self.pa is not None:
                self.line("accrual(%)", P._f4(100 * self.pa))
                self.line("by time", P._f4(self.ta))
        out.write("\nSurvival information:\n\n", "text")
        if c.get("hz_unknown"):
            self.line("hratio", P._f4(c["hr"]))
        elif self.surv:
            for h, nm in ((c["l1"], "1"), (c["l2"], "2")):
                sv = c.get("s" + nm) if c.get("s" + nm) is not None else math.exp(-h * c["time"])
                out.write(f"{'h' + nm:>{self.w}} = ", "text")
                out.write(P._f4(h), "result")
                out.write(f"{'s' + nm:>{self.w}} = ", "text")
                out.write(P._f4(sv) + "\n", "result")
            out.write(" " * (self.w + 12) + f"{'t':>{self.w}} = ", "text")
            out.write(P._f4(c["time"]) + "\n", "result")
        else:
            self.line("h1", P._f4(c["l1"]))
            self.line("h2", P._f4(c["l2"]))
        if self.loss is not None:
            out.write("\nLoss-to-follow-up information:\n\n", "text")
            for nm, v in self.loss:
                self.line(nm, P._f4(v))
        if what == "n":
            res = [("N", _I(c["N"]))]
            if c["N1"] == c["N2"]:
                res.append(("N per group", _I(c["N1"])))
            else:
                res += [("N1", _I(c["N1"])), ("N2", _I(c["N2"]))]
            self.results("Estimated sample sizes:", res)
        else:
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        if c.get("hz_unknown"):
            out.write("Note: Individual hazard rates are not reported because a\n"
                      "      control-group hazard rate is not specified.\n", "text")
        self.returns(c, what)

    def header(self, c, what, lines):
        if self.args:
            self.s.output.write("note: input parameters are " + ("survival probabilities" if self.surv else "hazard rates")
                                + "\n", "text")
        super().header(c, what, lines)

    def columns(self, c, what):
        return ["alpha", "power", "N", "delta", "h1", "h2", "hratio"]

    def row(self, c, what):
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"], "N": c["N"],
                "delta": c["delta"], "h1": c["l1"], "h2": c["l2"], "hratio": c["hr"]}

    _ORDER = ["power", "beta", "alpha", "delta", "nfractional", "onesided", "N", "N1", "N2", "nratio",
              "hratio", "lnhratio", "h1", "h2", "hdiff", "unconditional", "test", "effect", "type", "method"]

    def returns(self, c, what):
        # r(E0_vec), r(Ea_vec) and r(Pr_vec) of Stata are not stored: the manual does not define them
        pw = c["power_out"] if what == "power" else c["power"]
        v = {"power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": float(self.sides() == 1),
             "N": c["N"], "N1": c["N1"], "N2": c["N2"], "nratio": c["nratio"], "hratio": c["hr"],
             "lnhratio": math.log(c["hr"]), "h1": c["l1"], "h2": c["l2"], "hdiff": c["l2"] - c["l1"],
             "unconditional": float(self.unc),
             "test": "log hazard-ratio" if c["logtest"] else "hazard difference",
             "effect": "lnhratio" if c["logtest"] else "difference",          # VERIFY: effect with loghazard
             "type": "test", "method": "exponential"}
        if self.surv:
            v.update({"s1": math.exp(-c["l1"] * c["time"]), "s2": math.exp(-c["l2"] * c["time"]),
                      "time": c["time"]})
        if self.T is not None:
            v.update({"aperiod": self.r, "fperiod": self.f, "studytime": self.T, "ashape": self.gamma})
        if self.eta1 or self.eta2:
            v.update({"losshaz1": self.eta1, "losshaz2": self.eta2})
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)


class Logrank(_M):
    method = "logrank"
    iterates = False
    allowed = {"nratio", "hratio", "lnhratio", "schoenfeld", "direction", "onesided", "simpson", "wdprob"}

    def __init__(self, s, body, o):
        super().__init__(s, body, o)
        toks = P._split_args(body)
        if len(toks) > 2:
            _fail("too many values specified")
        self.args = [P._values(t) for t in toks]
        self.schoen = self.flag("schoenfeld")
        self.up_dir = str(o.get("direction") or "lower").lower().startswith("u")
        self.simp = None
        if o.get("simpson") not in (None, True):
            t = str(o["simpson"]).strip()
            v = (list(np.array(get_matrix(s, t).data, dtype=float).reshape(-1))
                 if re.fullmatch(r"[A-Za-z_]\w*", t) else P._values(t))
            if len(v) != 3 or any(not 0 < x < 1 for x in v):
                _fail("option simpson() incorrectly specified")
            self.simp = tuple(v)
            if self.args:
                _fail("simpson() may not be combined with surv1 and surv2")

    def base(self):
        b = self.common()
        for nm, v in zip(("s1", "s2"), self.args):
            b[nm] = v
        for k in ("hratio", "lnhratio", "wdprob"):
            v = self.vals(k)
            if v is not None:
                b[k] = v
        b["nratio"] = self.vals("nratio") or [1.0]
        return b

    def title(self, what):
        w = {"n": "sample sizes", "power": "power", "effect": "hazard ratio"}[what]
        return f"Estimated {w} for two-sample comparison of survivor functions"

    def hyp_lines(self, c, what):
        up = self.up_dir if what == "effect" else c["hr"] > 1
        op = self.hyp_op(c, up)
        sub = "Log-rank test, Schoenfeld method" if self.schoen else "Log-rank test, Freedman method"
        hyp = f"Ho: ln(HR) = 0  versus  Ha: ln(HR) {op} 0" if self.schoen else f"Ho: HR = 1  versus  Ha: HR {op} 1"
        if what == "effect" and self.sides() == 2:
            hyp += f"; {'ln(HR)' if self.schoen else 'HR'} {'>' if self.up_dir else '<'} {'0' if self.schoen else '1'}"
        return [sub, hyp]

    def compute(self, c):
        self.fix(c)
        a, sides, R = c["alpha"], self.sides(), c["nratio"]
        s1, s2 = c.get("s1"), c.get("s2")
        wd = c.get("wdprob", 0.0)
        eff = [k for k in ("hratio", "lnhratio") if c.get(k) is not None]
        if len(eff) > 1 or (eff and s1 is not None and s2 is not None):
            _fail("hratio() and lnhratio() may not be combined with each other or with surv2")
        have_n = c.get("n") is not None
        if have_n and c.get("power") is not None:
            what = "effect"
            if eff or s2 is not None:
                _fail("the hazard ratio may not be specified when computing the effect size")
        elif have_n:
            what = "power"
        else:
            what = "n"
            c.setdefault("power", 0.8)
        sch = self.schoen
        simp = self.simp

        def pe_of(hr):
            if simp is not None:
                return C.logrank_pe(hr, R, None, None, simp)
            return C.logrank_pe(hr, R, s1, None)

        if what == "effect":
            N = float(c["n"])
            sgn = 1.0 if self.up_dir else -1.0
            try:
                t = C.solve_up(lambda t: C.logrank_power(math.exp(sgn * t), R, N, pe_of(math.exp(sgn * t)),
                                                         a, sides, sch) - c["power"], 1e-9, 1.0)
            except ValueError:
                raise StataError(498, "failed to find a solution")
            hr = math.exp(sgn * t)
            c["_init"] = 1.0
        else:
            if s1 is not None and s2 is not None:
                hr = math.log(s2) / math.log(s1)
            elif eff:
                hr = c["hratio"] if eff == ["hratio"] else math.exp(c["lnhratio"])
            else:
                hr = 0.5
        if hr <= 0 or hr == 1:
            _fail("the hazard ratio must be positive and different from 1")
        pe = pe_of(hr)
        c.update({"hr": hr, "pe": pe, "wd": wd})
        c["delta"] = math.log(hr) if sch else hr
        if what == "n":
            ev = C.logrank_events(hr, R, a, sides, c["power"], sch)
            n_c = ev / pe / (1 - wd)                # N from the events before rounding (Stata 14, wdprob(.1))
            if self.flag("nfractional"):
                n1 = n_c / (1 + R)
                n2 = R * n1
            else:
                n1 = float(math.ceil(n_c / (1 + R) - 1e-9))
                n2 = float(math.ceil(R * n1 - 1e-9))
            N = n1 + n2
            ev = N * (1 - wd) * pe                  # events shown: from the final N (Stata 14)
            c["E"] = float(ev if self.flag("nfractional") else math.ceil(ev - 1e-9))
            c["power_out"] = C.logrank_power(hr, R, N * (1 - wd), pe, a, sides, sch)
        else:
            N = float(c["n"])
            n1 = N / (1 + R)
            n2 = N - n1
            ev = N * (1 - wd) * pe
            c["E"] = ev if self.flag("nfractional") else math.ceil(ev - 1e-9)             # Stata 14 shows 81 for 80.46
            c["power_out"] = c["power"] if what == "effect" else C.logrank_power(hr, R, N, pe, a, sides, sch)
        c["N"], c["N1"], c["N2"] = N, n1, n2
        return what

    def header(self, c, what, lines):
        out = self.s.output
        if self.simp is not None:
            s1 = self.simp
            s2 = [x ** c["hr"] for x in s1]
            out.write("note: probability of an event is computed using Simpson's rule with\n", "text")
            for g, v in (("1", s1), ("2", s2)):
                out.write(f"      S{g}(f) = {v[0]:.2f}, S{g}(f+r/2) = {v[1]:.2f}, S{g}(T) = {v[2]:.2f}\n", "text")
        self.iterates = what == "effect"              # the other cases are closed forms
        super().header(c, what, lines)

    def show(self, c, what):
        items = [("alpha", P._f4(c["alpha"]))]
        if what != "power":
            items.append(("power", P._f4(c["power"])))
        if what != "n":
            items += [("N", _I(c["N"])), ("N per group", _I(c["N1"]))] if c["N1"] == c["N2"] else [
                ("N", _I(c["N"])), ("N1", _I(c["N1"])), ("N2", _I(c["N2"]))]
        lab = "log hazard-ratio" if self.schoen else "hazard ratio"
        if what != "effect":
            items.append(("delta", P._f4(c["delta"]) + f"  ({lab})"))
            items.append(("hratio", P._f4(c["hr"])))
        if c["nratio"] != 1:
            items.append(("N2/N1", P._f4(c["nratio"])))
        self.put(items)
        cens = []
        if c.get("s1") is not None:
            cens += [("s1", P._f4(c["s1"])), ("s2", P._f4(c["s1"] ** c["hr"]))]
        cens.append(("Pr_E", P._f4(c["pe"])))
        if c["wd"]:
            cens.append(("Pr_w", P._f4(c["wd"])))
        if what == "n":
            self.results("Censoring and withdrawal:" if c["wd"] else "Censoring:", cens)
            res = [("E", _I(c["E"])), ("N", _I(c["N"]))]
            res += ([("N per group", _I(c["N1"]))] if c["N1"] == c["N2"] else
                    [("N1", _I(c["N1"])), ("N2", _I(c["N2"]))])
            self.results("Estimated number of events and sample sizes:", res)
        elif what == "power":
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results("Estimated power:", [("power", P._f4(c["power_out"]))])
        else:
            self.results("Number of events and censoring:", [("E", _I(c["E"]))] + cens)
            self.results("Estimated effect size and hazard ratio:",
                         [("delta", P._f4(c["delta"]) + f"  ({lab})"), ("hratio", P._f4(c["hr"]))])
        self.returns(c, what)

    def columns(self, c, what):
        return ["alpha", "power", "N", "N1", "N2", "E", "delta", "hratio", "Pr_E"]

    def row(self, c, what):
        return {"alpha": c["alpha"], "power": c["power_out"] if what == "power" else c["power"],
                "N": c["N"], "N1": c["N1"], "N2": c["N2"], "E": c["E"], "delta": c["delta"],
                "hratio": c["hr"], "Pr_E": c["pe"]}

    _ORDER = (["power", "beta", "alpha"] + _SOLVER + ["delta", "nfractional", "onesided", "N", "N1", "N2",
              "nratio", "hratio", "E", "s1", "s2", "Pr_E", "Pr_w", "type", "method", "test", "effect"])

    def returns(self, c, what):
        pw = c["power_out"] if what == "power" else c["power"]
        v = {"power": pw, "beta": 1 - pw, "alpha": c["alpha"], "delta": c["delta"],
             "nfractional": float(self.flag("nfractional")), "onesided": float(self.sides() == 1),
             "N": c["N"], "N1": c["N1"], "N2": c["N2"], "nratio": c["nratio"], "hratio": c["hr"],
             "E": c["E"], "Pr_E": c["pe"], "Pr_w": c["wd"]}
        if c.get("s1") is not None:
            v.update({"s1": c["s1"], "s2": c["s1"] ** c["hr"]})
        v.update(_solver_vals(c))
        v.update({"type": "test", "method": "logrank", "test": "Schoenfeld" if self.schoen else "Freedman",
                  "effect": "hratio"})                                  # VERIFY: effect with schoenfeld
        _store(self.s, {k: (float(x) if not isinstance(x, str) else x) for k, x in v.items()}, self._ORDER)
