"""Dialog boxes of the Statistics and Graphics menus.

Each dialog is described by a Spec (command, fields and how to build the line)
and drawn by a single generic CommandDialog. The dialog shows the command
that will run while the user fills in the fields and offers four buttons,
as in Stata: OK (runs and closes), Submit (runs and stays open), Copy (writes
the command to the Command window, without running it) and Cancel.

The drawing and layout are OpenDTA's own; they do not reproduce Stata's dialogs.
The generated command follows the syntax documented in the manuals.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from ..i18n import tr


# ---------------------------------------------------------------------------
# dialog descriptions (without Qt, to test how the command is built)
# ---------------------------------------------------------------------------

@dataclass
class Field:
    key: str
    label: str
    kind: str                  # var, vars, exp, check, opt, choice, weight
    opt: str = ""              # option template: "detail", "by({})", "level({})"
    choices: list = field(default_factory=list)   # [(label, value)]
    default: object = ""
    required: bool = False
    group: str = "main"        # main or options


@dataclass
class Spec:
    ident: str
    title: str
    cmd: str
    fields: list[Field]
    templates: tuple[str, ...] = ()          # menu items replaced
    ifin: bool = True
    weights: tuple[str, ...] = ()
    build: Callable[["Spec", dict], str] | None = None

    def command(self, values: dict) -> str:
        line = self.build(self, values) if self.build is not None else default_build(self, values)
        return " ".join(line.split()) if '"' not in line else line.strip()


def _opt(f: Field, v) -> str:
    if f.kind == "check":
        return f.opt if v else ""
    if f.kind == "choice":
        return f.opt.format(v) if (v and f.opt) else (v or "")
    v = str(v).strip()
    return f.opt.format(v) if v else ""


def positional(spec: Spec, values: dict) -> list[str]:
    parts = []
    for f in spec.fields:
        if f.group != "main" or f.kind not in ("var", "vars", "exp"):
            continue
        v = str(values.get(f.key, "") or "").strip()
        if v:
            parts.append(f.opt.format(v) if f.opt else v)
    return parts


def options(spec: Spec, values: dict, skip: tuple[str, ...] = ()) -> list[str]:
    out = []
    for f in spec.fields:
        if f.key in skip or (f.group == "main" and f.kind in ("var", "vars", "exp")):
            continue
        o = _opt(f, values.get(f.key, f.default))
        if o:
            out.append(o)
    return out


def qualifiers(spec: Spec, values: dict) -> list[str]:
    out = []
    w = values.get("_weight")
    wv = str(values.get("_weightvar", "") or "").strip()
    if w and wv:
        out.append(f"[{w}={wv}]")
    if spec.ifin:
        e = str(values.get("_if", "") or "").strip()
        r = str(values.get("_in", "") or "").strip()
        if e:
            out.append(f"if {e}")
        if r:
            out.append(f"in {r}")
    return out


def default_build(spec: Spec, values: dict, *, skip: tuple[str, ...] = ()) -> str:
    parts = [spec.cmd] + positional(spec, values) + qualifiers(spec, values)
    line = " ".join(parts)
    opts = options(spec, values, skip)
    if opts:
        line += ", " + " ".join(opts)
    return line


def missing_required(spec: Spec, values: dict) -> list[str]:
    return [f.label for f in spec.fields if f.required and not str(values.get(f.key, "") or "").strip()]


# -- montagens especiais -------------------------------------------------------

def _build_test_by(spec: Spec, values: dict) -> str:
    """ttest/prtest: one group (var == value) or two groups (by())."""
    var = str(values.get("var", "")).strip()
    by = str(values.get("by", "")).strip()
    val = str(values.get("value", "")).strip()
    line = f"{spec.cmd} {var}"
    if not by and val:
        line += f" == {val}"
    q = qualifiers(spec, values)
    if q:
        line += " " + " ".join(q)
    opts = [f"by({by})"] if by else []
    opts += options(spec, values, skip=("by", "value"))
    return line + (", " + " ".join(opts) if opts else "")


def _build_scatter(spec: Spec, values: dict) -> str:
    y = str(values.get("y", "")).strip()
    x = str(values.get("x", "")).strip()
    q = " ".join(qualifiers(spec, values))
    q = (" " + q) if q else ""
    opts = options(spec, values, skip=("lfit",))
    if values.get("lfit"):
        line = f"twoway (scatter {y} {x}{q}) (lfit {y} {x}{q})"
    else:
        line = f"twoway scatter {y} {x}{q}"
    return line + (", " + " ".join(opts) if opts else "")


def _build_svy(spec: Spec, values: dict) -> str:
    sub = str(values.get("subpop", "")).strip()
    prefix = f"svy, subpop({sub}): " if sub else "svy: "
    return prefix + default_build(spec, values, skip=("subpop",))


def _build_svyset(spec: Spec, values: dict) -> str:
    psu = str(values.get("psu", "")).strip() or "_n"
    w = str(values.get("pw", "")).strip()
    line = f"svyset {psu}"
    if w:
        line += f" [pweight={w}]"
    opts = options(spec, values, skip=("psu", "pw"))
    return line + (", " + " ".join(opts) if opts else "")


def _build_tset(spec: Spec, values: dict) -> str:
    return " ".join([spec.cmd] + positional(spec, values))


# -- campos reutilizados -----------------------------------------------------------

def _robust() -> Field:
    return Field("robust", "Robust standard errors", "check", "vce(robust)", group="options")


def _cluster() -> Field:
    return Field("cluster", "Cluster variable:", "var", "vce(cluster {})", group="options")


def _level() -> Field:
    return Field("level", "Confidence level:", "opt", "level({})", group="options")


def _title() -> Field:
    return Field("title", "Graph title:", "opt", 'title("{}")', group="options")


def _model(ident, title, cmd, extra=(), templates=None, weights=("fweight", "pweight", "iweight")) -> Spec:
    return Spec(ident, title, cmd,
                [Field("dep", "Dependent variable:", "var", required=True),
                 Field("indep", "Independent variables:", "vars"),
                 *extra, _robust(), _cluster(), _level()],
                templates=templates or (cmd + " ",), weights=weights)


_FAMILIES = [("Gaussian", "gaussian"), ("Binomial", "binomial"), ("Poisson", "poisson"),
             ("Negative binomial", "nbinomial"), ("Gamma", "gamma"), ("Inverse Gaussian", "igaussian")]
_LINKS = [("(default)", ""), ("Identity", "identity"), ("Log", "log"), ("Logit", "logit"),
          ("Probit", "probit"), ("Complementary log-log", "cloglog"), ("Power -1", "power -1")]


SPECS: list[Spec] = [
    Spec("summarize", "summarize - Summary statistics", "summarize",
         [Field("vars", "Variables:", "vars"),
          Field("detail", "Detail", "check", "detail", group="options")],
         templates=("summarize ",), weights=("aweight", "fweight")),
    Spec("tab1", "tabulate - One-way table", "tabulate",
         [Field("var", "Variable:", "var", required=True),
          Field("missing", "Treat missing values like other values", "check", "missing", group="options"),
          Field("nolabel", "Display numeric codes instead of value labels", "check", "nolabel",
                group="options")],
         templates=("tabulate ",), weights=("fweight", "aweight")),
    Spec("tab2", "tabulate - Two-way table", "tabulate",
         [Field("row", "Row variable:", "var", required=True),
          Field("col", "Column variable:", "var", required=True),
          Field("chi2", "Chi-squared test", "check", "chi2", default=True, group="options"),
          Field("exact", "Fisher's exact test", "check", "exact", group="options"),
          Field("rowp", "Row percentages", "check", "row", group="options"),
          Field("colp", "Column percentages", "check", "column", group="options"),
          Field("missing", "Treat missing values like other values", "check", "missing", group="options")],
         templates=("tabulate  , chi2",), weights=("fweight", "aweight")),
    Spec("tabstat", "tabstat - Table of statistics", "tabstat",
         [Field("vars", "Variables:", "vars", required=True),
          Field("stats", "Statistic:", "opt", "statistics({})", default="mean sd n", group="options"),
          Field("by", "Group variable:", "var", "by({})", group="options")],
         templates=("tabstat , statistics(mean sd n)",), weights=("aweight", "fweight")),
    Spec("correlate", "correlate - Correlations", "correlate",
         [Field("vars", "Variables:", "vars"),
          Field("cov", "Covariances instead of correlations", "check", "covariance", group="options")],
         templates=("correlate ",), weights=("aweight", "fweight")),
    Spec("ttest", "ttest - t test", "ttest",
         [Field("var", "Variable:", "var", required=True),
          Field("value", "Hypothesized mean:", "opt"),
          Field("by", "Two-sample test by group:", "var", group="options"),
          Field("unequal", "Unequal variances", "check", "unequal", group="options"),
          _level()],
         templates=("ttest ",), build=_build_test_by),
    Spec("prtest", "prtest - Proportion test", "prtest",
         [Field("var", "Variable:", "var", required=True),
          Field("value", "Hypothesized proportion:", "opt"),
          Field("by", "Two-sample test by group:", "var", group="options"),
          _level()],
         templates=("prtest ",), build=_build_test_by),
    Spec("ci", "ci - Confidence intervals for means", "ci means",
         [Field("vars", "Variables:", "vars"), _level()],
         templates=("ci means ",), weights=("aweight", "fweight")),
    Spec("anova", "oneway - One-way analysis of variance", "oneway",
         [Field("var", "Response variable:", "var", required=True),
          Field("group", "Group variable:", "var", required=True),
          Field("tab", "Summary table", "check", "tabulate", group="options"),
          Field("bonf", "Bonferroni comparisons", "check", "bonferroni", group="options")],
         templates=("oneway ",), weights=("aweight", "fweight")),
    Spec("swilk", "swilk - Shapiro-Wilk normality test", "swilk",
         [Field("vars", "Variables:", "vars", required=True)],
         templates=("swilk ",)),
    _model("regress", "regress - Linear regression", "regress",
           [Field("beta", "Standardized beta coefficients", "check", "beta", group="options")],
           weights=("aweight", "fweight", "pweight", "iweight")),
    _model("logistic", "logistic - Logistic regression, reporting odds ratios", "logistic"),
    _model("logit", "logit - Logistic regression", "logit",
           [Field("or", "Odds ratios", "check", "or", group="options")]),
    _model("probit", "probit - Probit regression", "probit"),
    _model("ologit", "ologit - Ordered logistic regression", "ologit",
           [Field("or", "Odds ratios", "check", "or", group="options")]),
    _model("mlogit", "mlogit - Multinomial logistic regression", "mlogit",
           [Field("base", "Base outcome:", "opt", "baseoutcome({})", group="options"),
            Field("rrr", "Relative-risk ratios", "check", "rrr", group="options")]),
    _model("poisson", "poisson - Poisson regression", "poisson",
           [Field("exposure", "Exposure variable:", "var", "exposure({})", group="options"),
            Field("irr", "Incidence-rate ratios", "check", "irr", group="options")]),
    _model("nbreg", "nbreg - Negative binomial regression", "nbreg",
           [Field("exposure", "Exposure variable:", "var", "exposure({})", group="options"),
            Field("irr", "Incidence-rate ratios", "check", "irr", group="options")]),
    _model("glm", "glm - Generalized linear models", "glm",
           [Field("family", "Family:", "choice", "family({})", _FAMILIES, default="gaussian", group="options"),
            Field("link", "Link:", "choice", "link({})", _LINKS, group="options"),
            Field("eform", "Exponentiated coefficients", "check", "eform", group="options")],
           templates=("glm , family() link()",)),
    _model("qreg", "qreg - Quantile regression", "qreg",
           [Field("q", "Quantile:", "opt", "quantile({})", default="0.5", group="options")],
           weights=("aweight", "fweight", "iweight")),
    Spec("xtreg", "xtreg - Panel-data linear models", "xtreg",
         [Field("dep", "Dependent variable:", "var", required=True),
          Field("indep", "Independent variables:", "vars"),
          Field("est", "Estimator:", "choice", "",
                [("Fixed effects", "fe"), ("Random effects", "re"), ("Between effects", "be"),
                 ("Maximum likelihood", "mle")], default="fe", group="options"),
          _robust(), _level()],
         templates=("xtreg , fe", "xtreg , re"), weights=()),
    Spec("stset", "stset - Declare survival-time data", "stset",
         [Field("time", "Time variable:", "var", required=True),
          Field("fail", "Failure variable:", "opt", "failure({})"),
          Field("id", "Subject ID variable:", "var", "id({})", group="options"),
          Field("origin", "Origin (time expression):", "opt", "origin(time {})", group="options"),
          Field("scale", "Scale:", "opt", "scale({})", group="options")],
         templates=("stset , failure()",), ifin=False),
    Spec("stcox", "stcox - Cox proportional hazards model", "stcox",
         [Field("indep", "Independent variables:", "vars", required=True),
          Field("efron", "Efron method for ties", "check", "efron", group="options"),
          Field("nohr", "Report coefficients, not hazard ratios", "check", "nohr", group="options"),
          Field("strata", "Strata variables:", "opt", "strata({})", group="options"),
          _robust(), _level()],
         templates=("stcox ",)),
    Spec("streg", "streg - Parametric survival models", "streg",
         [Field("indep", "Independent variables:", "vars"),
          Field("dist", "Distribution:", "choice", "distribution({})",
                [("Exponential", "exponential"), ("Weibull", "weibull"), ("Gompertz", "gompertz"),
                 ("Lognormal", "lognormal"), ("Loglogistic", "loglogistic")], default="weibull"),
          Field("time", "Accelerated failure-time metric", "check", "time", group="options"),
          _robust(), _level()],
         templates=("streg , distribution(weibull)",)),
    Spec("stsgraph", "sts graph - Kaplan-Meier survivor function", "sts graph",
         [Field("by", "Group variable:", "var", "by({})", group="options"),
          Field("failure", "Failure function (1 - S)", "check", "failure", group="options"),
          Field("ci", "Confidence interval", "check", "ci", group="options"),
          _title()],
         templates=("sts graph",)),
    Spec("ststest", "sts test - Test equality of survivor functions", "sts test",
         [Field("vars", "Group variable:", "vars", required=True),
          Field("test", "Test:", "choice", "",
                [("Log-rank", "logrank"), ("Wilcoxon (Breslow)", "wilcoxon"),
                 ("Tarone-Ware", "tware"), ("Peto-Peto", "peto")], default="logrank", group="options")],
         templates=("sts test ",)),
    Spec("xtset", "xtset - Declare panel data", "xtset",
         [Field("panel", "Panel variable:", "var", required=True),
          Field("time", "Time variable:", "var")],
         templates=("xtset ",), ifin=False, build=_build_tset),
    Spec("tsset", "tsset - Declare time-series data", "tsset",
         [Field("time", "Time variable:", "var", required=True)],
         templates=("tsset ",), ifin=False, build=_build_tset),
    Spec("svyset", "svyset - Declare survey design", "svyset",
         [Field("psu", "PSU variable:", "var"),
          Field("pw", "Sampling weight:", "var"),
          Field("strata", "Strata variable:", "var", "strata({})", group="options"),
          Field("fpc", "FPC variable:", "var", "fpc({})", group="options"),
          Field("vce", "Variance estimation:", "choice", "vce({})",
                [("Linearized", ""), ("Jackknife", "jackknife"), ("BRR", "brr"),
                 ("Bootstrap", "bootstrap"), ("SDR", "sdr")], group="options"),
          Field("single", "Strata with a single PSU:", "choice", "singleunit({})",
                [("(missing)", ""), ("Certainty", "certainty"), ("Scaled", "scaled"),
                 ("Centered", "centered")], group="options")],
         templates=("svyset ",), ifin=False, build=_build_svyset),
    Spec("svymean", "svy: mean - Survey means", "mean",
         [Field("vars", "Variables:", "vars", required=True),
          Field("over", "Over:", "vars", "over({})", group="options"),
          Field("subpop", "Subpopulation (variable):", "var", group="options"), _level()],
         templates=("svy: mean ",), build=_build_svy),
    Spec("svyprop", "svy: proportion - Survey proportions", "proportion",
         [Field("vars", "Variables:", "vars", required=True),
          Field("over", "Over:", "vars", "over({})", group="options"),
          Field("subpop", "Subpopulation (variable):", "var", group="options"), _level()],
         templates=("svy: proportion ",), build=_build_svy),
    Spec("svyreg", "svy: regress - Survey linear regression", "regress",
         [Field("dep", "Dependent variable:", "var", required=True),
          Field("indep", "Independent variables:", "vars"),
          Field("subpop", "Subpopulation (variable):", "var", group="options"), _level()],
         templates=("svy: regress ",), build=_build_svy),
    Spec("svytab", "svy: tabulate - Survey tables", "tabulate",
         [Field("row", "Row variable:", "var", required=True),
          Field("col", "Column variable:", "var"),
          Field("rowp", "Row percentages", "check", "row", group="options"),
          Field("colp", "Column percentages", "check", "column", group="options"),
          Field("ci", "Confidence intervals", "check", "ci", group="options"),
          Field("percent", "Percentages instead of proportions", "check", "percent", group="options"),
          Field("subpop", "Subpopulation (variable):", "var", group="options")],
         templates=("svy: tabulate ",), build=_build_svy),
    Spec("histogram", "histogram - Histogram", "histogram",
         [Field("var", "Variable:", "var", required=True),
          Field("bins", "Number of bins:", "opt", "bin({})", group="options"),
          Field("scale", "Scale:", "choice", "",
                [("Density", ""), ("Fraction", "fraction"), ("Frequency", "frequency"),
                 ("Percent", "percent")], group="options"),
          Field("normal", "Normal density", "check", "normal", group="options"),
          _title()],
         templates=("histogram ",), weights=("fweight",)),
    Spec("scatter", "twoway scatter - Scatterplot", "twoway scatter",
         [Field("y", "Y variable:", "var", required=True),
          Field("x", "X variable:", "var", required=True),
          Field("lfit", "Add linear fit", "check", group="options"),
          _title()],
         templates=("twoway scatter ",), build=_build_scatter),
    Spec("box", "graph box - Box plot", "graph box",
         [Field("vars", "Variables:", "vars", required=True),
          Field("over", "Over:", "var", "over({})", group="options"), _title()],
         templates=("graph box ",)),
    Spec("bar", "graph bar - Bar chart of means", "graph bar (mean)",
         [Field("vars", "Variables:", "vars", required=True),
          Field("over", "Over:", "var", "over({})", group="options"), _title()],
         templates=("graph bar (mean) ",)),
    Spec("kdensity", "kdensity - Kernel density", "kdensity",
         [Field("var", "Variable:", "var", required=True),
          Field("normal", "Normal density", "check", "normal", group="options"), _title()],
         templates=("kdensity ",)),
]

_BY_TEMPLATE = {t: sp for sp in SPECS for t in sp.templates}
_BY_ID = {sp.ident: sp for sp in SPECS}


def spec_for(template: str) -> Spec | None:
    return _BY_TEMPLATE.get(template)


def spec(ident: str) -> Spec:
    return _BY_ID[ident]


def initial_values(sp: Spec, template: str = "") -> dict:
    vals = {f.key: f.default for f in sp.fields}
    if sp.ident == "xtreg" and template.endswith("re"):
        vals["est"] = "re"
    return vals


# ---------------------------------------------------------------------------
# a janela
# ---------------------------------------------------------------------------

def _qt():
    from PySide6 import QtCore, QtGui, QtWidgets
    return QtCore, QtGui, QtWidgets


def CommandDialog(sp: Spec, session, parent=None, *, run=None, to_command=None,  # noqa: N802
                  template: str = ""):
    """Creates the dialog (factory function: Qt is only imported when needed)."""
    return _dialog_class()(sp, session, parent, run=run, to_command=to_command, template=template)


_CLASS = None


def _dialog_class():
    global _CLASS
    if _CLASS is not None:
        return _CLASS
    QtCore, QtGui, QtWidgets = _qt()
    Qt = QtCore.Qt

    class _VarList(QtWidgets.QWidget):
        """Variable-list field: free text and a button that adds
        names from the dataset in memory."""
        changed = QtCore.Signal()

        def __init__(self, names: list[str], single: bool):
            super().__init__()
            lay = QtWidgets.QHBoxLayout(self)
            lay.setContentsMargins(0, 0, 0, 0)
            lay.setSpacing(2)
            self.single = single
            self.edit = QtWidgets.QLineEdit()
            comp = QtWidgets.QCompleter(names, self)
            comp.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
            if single:
                self.edit.setCompleter(comp)
            self.edit.textChanged.connect(self.changed)
            btn = QtWidgets.QToolButton()
            btn.setText("▾")
            btn.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
            menu = QtWidgets.QMenu(btn)
            for n in names[:400]:
                menu.addAction(n, lambda n=n: self._add(n))
            btn.setMenu(menu)
            btn.setEnabled(bool(names))
            lay.addWidget(self.edit, 1)
            lay.addWidget(btn)

        def _add(self, name: str) -> None:
            if self.single:
                self.edit.setText(name)
            else:
                cur = self.edit.text().rstrip()
                self.edit.setText((cur + " " if cur else "") + name)

        def value(self) -> str:
            return self.edit.text()

        def set_value(self, v) -> None:
            self.edit.setText(str(v or ""))

    class _Dialog(QtWidgets.QDialog):
        def __init__(self, sp: Spec, session, parent=None, *, run=None, to_command=None,
                     template: str = ""):
            super().__init__(parent)
            self.spec = sp
            self.session = session
            self._run = run
            self._to_command = to_command
            self.setWindowTitle(sp.title)
            self.widgets: dict[str, object] = {}
            names = list(session.data.names) if session is not None else []

            lay = QtWidgets.QVBoxLayout(self)
            main = QtWidgets.QFormLayout()
            opts_box = QtWidgets.QGroupBox(tr("Options:").rstrip(":"))
            opts = QtWidgets.QFormLayout(opts_box)
            vals = initial_values(sp, template)
            for f in sp.fields:
                w = self._make(f, names)
                self.widgets[f.key] = w
                self._set(f, w, vals.get(f.key))
                target = main if f.group == "main" else opts
                if f.kind == "check":
                    target.addRow("", w)
                else:
                    target.addRow(tr(f.label), w)
            lay.addLayout(main)
            if opts.rowCount():
                lay.addWidget(opts_box)

            if sp.ifin or sp.weights:
                box = QtWidgets.QGroupBox("if/in/" + tr("Weights:").rstrip(":").lower())
                form = QtWidgets.QFormLayout(box)
                if sp.ifin:
                    self.if_edit = QtWidgets.QLineEdit()
                    self.in_edit = QtWidgets.QLineEdit()
                    self.if_edit.textChanged.connect(self._update)
                    self.in_edit.textChanged.connect(self._update)
                    form.addRow(tr("If (expression):"), self.if_edit)
                    form.addRow(tr("In (range):"), self.in_edit)
                if sp.weights:
                    row = QtWidgets.QHBoxLayout()
                    self.wtype = QtWidgets.QComboBox()
                    self.wtype.addItem("", "")
                    for w in sp.weights:
                        self.wtype.addItem(w, w)
                    self.wvar = _VarList(names, True)
                    self.wtype.currentIndexChanged.connect(self._update)
                    self.wvar.changed.connect(self._update)
                    row.addWidget(self.wtype)
                    row.addWidget(self.wvar, 1)
                    form.addRow(tr("Weights:"), row)
                lay.addWidget(box)

            lay.addWidget(QtWidgets.QLabel(tr("The command will be:")))
            self.preview = QtWidgets.QLineEdit()
            self.preview.setReadOnly(True)
            mono = QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont)
            self.preview.setFont(mono)
            lay.addWidget(self.preview)

            buttons = QtWidgets.QHBoxLayout()
            self.btn_copy = QtWidgets.QPushButton(tr("Copy to Command window"))
            self.btn_ok = QtWidgets.QPushButton(tr("OK"))
            self.btn_submit = QtWidgets.QPushButton(tr("Submit"))
            self.btn_cancel = QtWidgets.QPushButton(tr("Cancel"))
            self.btn_ok.setDefault(True)
            self.btn_copy.clicked.connect(self._copy)
            self.btn_ok.clicked.connect(self._ok)
            self.btn_submit.clicked.connect(self._submit)
            self.btn_cancel.clicked.connect(self.reject)
            buttons.addWidget(self.btn_copy)
            buttons.addStretch(1)
            for b in (self.btn_ok, self.btn_submit, self.btn_cancel):
                buttons.addWidget(b)
            lay.addLayout(buttons)
            self.setMinimumWidth(520)
            self._update()

        def _make(self, f: Field, names):
            if f.kind in ("var", "vars"):
                w = _VarList(names, f.kind == "var")
                w.changed.connect(self._update)
                return w
            if f.kind == "check":
                w = QtWidgets.QCheckBox(tr(f.label))
                w.toggled.connect(self._update)
                return w
            if f.kind == "choice":
                w = QtWidgets.QComboBox()
                for lab, val in f.choices:
                    w.addItem(tr(lab), val)
                w.currentIndexChanged.connect(self._update)
                return w
            w = QtWidgets.QLineEdit()
            w.textChanged.connect(self._update)
            return w

        @staticmethod
        def _set(f: Field, w, v) -> None:
            if f.kind == "check":
                w.setChecked(bool(v))
            elif f.kind == "choice":
                k = w.findData(v)
                w.setCurrentIndex(max(k, 0))
            elif f.kind in ("var", "vars"):
                w.set_value(v)
            else:
                w.setText(str(v or ""))

        def values(self) -> dict:
            out = {}
            for f in self.spec.fields:
                w = self.widgets[f.key]
                if f.kind == "check":
                    out[f.key] = w.isChecked()
                elif f.kind == "choice":
                    out[f.key] = w.currentData() or ""
                elif f.kind in ("var", "vars"):
                    out[f.key] = w.value()
                else:
                    out[f.key] = w.text()
            if self.spec.ifin:
                out["_if"] = self.if_edit.text()
                out["_in"] = self.in_edit.text()
            if self.spec.weights:
                out["_weight"] = self.wtype.currentData() or ""
                out["_weightvar"] = self.wvar.value()
            return out

        def command(self) -> str:
            return self.spec.command(self.values())

        def _update(self) -> None:
            if not hasattr(self, "preview"):
                return
            self.preview.setText(self.command())
            ok = not missing_required(self.spec, self.values())
            for b in (self.btn_ok, self.btn_submit):
                b.setEnabled(ok)

        def _copy(self) -> None:
            if self._to_command is not None:
                self._to_command(self.command())

        def _submit(self) -> None:
            if self._run is not None:
                self._run(self.command())

        def _ok(self) -> None:
            self._submit()
            self.accept()

    _CLASS = _Dialog
    return _CLASS
