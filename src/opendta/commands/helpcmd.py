"""help, view and search.

help looks for the topic first in the adopath folders (.sthlp and
.hlp files of third-party packages) and then in OpenDTA's own pages
(src/opendta/help). With the graphical interface, the page opens in the Viewer;
without it, the text goes to Results.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from ..core.errors import StataError
from ..core.smcl import render_text
from ..lang.adopath import find_file
from ..lang.words import strip_outer_quotes
from .registry import REGISTRY, command

if TYPE_CHECKING:
    from ..session import Session

HELP_DIR = Path(__file__).resolve().parents[1] / "help"

ALIASES = {
    "cloglog": "binmodels", "scobit": "binmodels", "hetprobit": "binmodels", "fracreg": "binmodels", "glogit": "binmodels", "binreg": "binmodels", "zip": "binmodels", "zinb": "binmodels", "tpoisson": "binmodels", "tnbreg": "binmodels", "gnbreg": "binmodels",
    "intreg": "censmodels", "truncreg": "censmodels", "heckman": "censmodels", "heckprobit": "censmodels", "heckoprobit": "censmodels", "biprobit": "censmodels",
    "ivprobit": "ivmodels", "ivtobit": "ivmodels", "ivpoisson": "ivmodels",
    "cnsreg": "sysmodels", "constraint": "sysmodels", "sureg": "sysmodels", "reg3": "sysmodels", "eivreg": "sysmodels", "vwls": "sysmodels", "rreg": "sysmodels",
    "clogit": "choicemodels", "rologit": "choicemodels", "asclogit": "choicemodels", "slogit": "choicemodels", "mprobit": "choicemodels", "asmprobit": "choicemodels", "asroprobit": "choicemodels",
    "exlogistic": "exactmodels", "expoisson": "exactmodels",
    "nl": "nlmodels", "nlsur": "nlmodels", "gmm": "nlmodels", "frontier": "nlmodels",
    "irtgraph": "irt", "bayes": "bayesmh", "bayesstats": "bayesmh", "bayestest": "bayesmh", "bayesgraph": "bayesmh",
    "teffects": "treatment", "eteffects": "treatment", "stteffects": "treatment", "etregress": "treatment", "etpoisson": "treatment",
    "xtprobit": "xtquad", "xtcloglog": "xtquad", "xtologit": "xtquad", "xtoprobit": "xtquad", "xttobit": "xtquad", "xtintreg": "xtquad", "xtnbreg": "xtquad", "quadchk": "xtquad", "xtgls": "xtlinear", "xtpcse": "xtlinear", "xtregar": "xtlinear", "xtrc": "xtlinear", "xthtaylor": "xtlinear", "xtivreg": "xtlinear", "xtabond": "xtdpd", "xtdpdsys": "xtdpd", "estat_abond": "xtdpd", "estat_sargan": "xtdpd", "xtline": "xtdesc", "xttab": "xtdesc", "xtdata": "xtdesc",
    "jackknife": "resample", "permute": "resample", "simulate": "resample",
    "mlexp": "ml", "mleval": "ml", "mlsum": "ml", "mlvecsum": "ml", "mlmatsum": "ml",
    "mkspline": "modelbuild", "xi": "modelbuild", "nestreg": "modelbuild", "stepwise": "modelbuild", "fp": "modelbuild", "mfp": "modelbuild",
    "testnl": "postest2", "predictnl": "postest2", "linktest": "postest2", "contrast": "postest2", "pwcompare": "postest2", "pwmean": "postest2", "suest": "postest2",
    "bitest": "basetests", "bitesti": "basetests", "ksmirnov": "basetests", "median": "basetests", "nptrend": "basetests", "runtest": "basetests", "tetrachoric": "basetests", "ktau": "basetests", "pcorr": "basetests",
    "ameans": "basestats", "icc": "basestats", "loneway": "basestats", "esize": "basestats", "esizei": "basestats", "lnskew0": "basestats", "ladder": "basestats", "boxcox": "basestats",
    "quantile": "basegraphs", "symplot": "basegraphs", "serrbar": "basegraphs", "dotplot": "basegraphs", "sunflower": "basegraphs", "lpoly": "basegraphs", "cusum": "basegraphs", "gladder": "basegraphs",
    "cumul": "basedisp", "stem": "basedisp", "lv": "basedisp", "smooth": "basedisp",
    "stjoin": "stsplit", "stfill": "stsplit", "stvary": "stsplit", "stgen": "stsplit", "stbase": "stsplit", "snapspan": "stsplit", "stptime": "strate", "stci": "strate", "stir": "strate", "stmh": "strate", "stmc": "strate", "ltable": "strate", "ctset": "ct", "cttost": "ct", "sttoct": "ct", "sttocc": "ct",
    "moptimize": "mata", "optimize": "mata", "m2_class": "mata", "class": "mata",
    "xtmixed": "mixed", "meprobit": "melogit", "mepoisson": "melogit", "mecloglog": "meglm", "menbreg": "meglm",
    "meologit": "meglm", "meoprobit": "meglm", "meqrlogit": "meglm", "meqrpoisson": "meglm", "estat_icc": "mixed",
    "varsoc": "var", "vargranger": "var",
    "tsappend": "tsdata", "tsfill": "tsdata", "tsreport": "tsdata", "rolling": "tsdata",
    "tssmooth": "tsfilter", "pergram": "tsspectral", "cumsp": "tsspectral", "xcorr": "tsspectral",
    "wntestb": "tsspectral", "wntestq": "tsspectral", "psdensity": "tsspectral",
    "dfgls": "tsunit", "pperron": "tsunit", "varlmar": "varpost", "varnorm": "varpost",
    "varstable": "varpost", "varwle": "varpost", "varbasic": "varpost", "fcast": "varpost",
    "vecrank": "vec", "svar": "vec", "arch": "archmodels", "arfima": "archmodels",
    "ucm": "ssm", "sspace": "ssm", "dfactor": "ssm",
    "factor": "pca", "rotate": "pca", "alpha": "pca", "screeplot": "pca", "loadingplot": "pca", "scoreplot": "pca",
    "mvreg": "hotelling", "manova": "hotelling", "manovatest": "hotelling", "canon": "hotelling", "mvtest": "hotelling",
    "candisc": "discrim", "mca": "ca", "mds": "ca", "mdsmat": "ca", "procrustes": "ca", "biplot": "ca", "camat": "ca", "cabiplot": "ca", "caprojection": "ca", "mcaplot": "ca", "mcaprojection": "ca", "mdslong": "ca", "mdsconfig": "ca", "mdsshepard": "ca", "procoverlay": "ca", "clustermat": "cluster",
    "st": "stset", "streset": "stset", "stdes": "stdescribe", "estat_phtest": "stcox",
    "d": "describe", "de": "describe", "des": "describe", "g": "generate", "gen": "generate",
    "replace": "generate", "l": "list", "li": "list", "saveold": "save", "export": "import",
    "insheet": "import", "outsheet": "import", "import_delimited": "import",
    "import_excel": "import", "export_excel": "import", "export_delimited": "import",
    "infix": "infile", "restore": "preserve", "tempvar": "preserve", "tempfile": "preserve",
    "tempname": "preserve", "return": "program", "ereturn": "program", "sreturn": "program",
    "args": "syntax", "gettoken": "syntax", "marksample": "syntax", "markout": "syntax",
    "bysort": "by", "bys": "by", "missing_values": "missing", "mi_values": "missing",
    "mat": "matrix", "cmdlog": "log", "set_hints": "hints", "hint": "hints",
    "set_superscript": "superscript", "set_numerics": "numerics",
    "su": "summarize", "sum": "summarize", "summ": "summarize", "tabstat": "summarize",
    "centile": "summarize", "pctile": "summarize", "_pctile": "summarize", "xtile": "summarize",
    "correlate": "summarize", "corr": "summarize", "pwcorr": "summarize",
    "prtest": "ttest", "ci": "ttest", "ta": "tabulate", "tab": "tabulate", "tab1": "tabulate",
    "tab2": "tabulate", "table": "tabulate", "contract": "collapse", "expand": "collapse",
    "fillin": "collapse", "append": "merge", "joinby": "merge", "cross": "merge",
    "encode": "recode", "decode": "recode", "destring": "recode", "tostring": "recode",
    "split": "recode", "mvencode": "recode", "mvdecode": "recode", "isid": "duplicates",
    "levelsof": "duplicates", "seed": "set_seed", "rng": "set_seed", "set_rng": "set_seed",
    "runiform": "set_seed", "rnormal": "set_seed",
    "gr": "graph", "graph_export": "graph", "graph_save": "graph", "graph_use": "graph",
    "graph_dir": "graph", "graph_drop": "graph", "graph_rename": "graph", "graph_combine": "graph",
    "graph_display": "graph", "graph_copy": "graph", "set_scheme": "graph", "scheme": "graph",
    "set_graphics": "graph", "graph_twoway": "twoway", "tw": "twoway", "scatter": "twoway",
    "sc": "twoway", "line": "twoway", "connected": "twoway", "tsline": "twoway", "lfit": "twoway",
    "qfit": "twoway", "lfitci": "twoway", "qfitci": "twoway", "twoway_function": "twoway",
    "function": "twoway", "rcap": "twoway", "rarea": "twoway", "hist": "histogram",
    "kdensity": "histogram", "graph_hbar": "graph_bar", "graph_box": "graph_bar",
    "graph_hbox": "graph_bar", "graph_pie": "graph_bar", "graph_matrix": "graph_bar",
    "graph_dot": "graph_bar", "doedit": "doedit",
    # pages that cover several commands (phase 9)
    "areg": "regress", "ivregress": "regress", "reg": "regress", "logistic": "logit",
    "probit": "logit", "oprobit": "ologit", "mlogit": "ologit", "nbreg": "poisson",
    "xtset": "xtreg", "tsset": "xtreg", "testparm": "test", "lincom": "test", "nlcom": "test",
    "est": "estimates", "estimates_table": "estimates", "estimates_store": "estimates",
    "proportion": "mean", "ratio": "mean", "total": "mean", "ranksum": "nptests",
    "signrank": "nptests", "kwallis": "nptests", "spearman": "nptests", "swilk": "nptests",
    "sktest": "nptests", "sdtest": "nptests", "cii": "immediate", "tabi": "immediate",
    "ttesti": "immediate", "prtesti": "immediate", "drop": "datamgmt", "keep": "datamgmt",
    "sort": "datamgmt", "gsort": "datamgmt", "rename": "datamgmt", "order": "datamgmt",
    "count": "datamgmt", "notes": "label", "char": "label", "recast": "format",
    "compress": "format", "di": "display", "local": "macro", "global": "macro",
    "scalar": "macro", "tokenize": "macro", "set_language": "language",
    "set_messages": "language", "messages": "language", "cd": "files", "pwd": "files",
    "erase": "files", "rm": "files", "type": "files", "which": "files", "findfile": "files",
    "sysdir": "files", "adopath": "files", "discard": "files", "dir": "files", "copy": "files",
    "run": "do", "include": "do", "capture": "do", "exit": "do", "continue": "do",
    "error": "do", "confirm": "do", "version": "do", "timer": "do", "quietly": "do",
    "browse": "interface", "edit": "interface", "varmanage": "interface", "help": "interface",
    "search": "interface", "view": "interface", "about": "interface", "creturn": "interface",
    "clear": "interface", "file": "filecmd", "misstable": "codebook", "svyset": "svy",
    "svydescribe": "svy",
    "cs": "epitab", "csi": "epitab", "cc": "epitab", "cci": "epitab", "ir": "epitab", "iri": "epitab",
    "mcc": "epitab", "mcci": "epitab", "tabodds": "epitab", "mhodds": "epitab", "istdize": "dstdize",
    "kap": "kappa", "roctab": "roc", "roccomp": "roc", "rocfit": "roc", "lroc": "roc", "lsens": "roc",
    "copy": "filesys", "mkdir": "filesys", "rmdir": "filesys", "dir": "filesys", "ls": "filesys",
    "shell": "filesys", "filefilter": "filesys", "zipfile": "filesys", "unzipfile": "filesys",
    "checksum": "filesys", "clonevar": "datautil", "insobs": "datautil", "range": "datautil",
    "separate": "datautil", "ipolate": "datautil", "expandcl": "datautil", "stack": "datautil",
    "xpose": "datautil", "lookfor": "ds", "inspect": "ds", "labelbook": "ds", "cf": "compare",
    "datasignature": "compare", "snapshot": "compare", "corr2data": "drawnorm", "webuse": "sysuse",
    "icd9": "icd", "icd9p": "icd", "icd10": "icd",
    "unab": "varlistcmds", "tsunab": "varlistcmds", "fvunab": "varlistcmds",
    "unabcmd": "varlistcmds", "fvexpand": "varlistcmds", "fvrevar": "varlistcmds",
    "tsrevar": "varlistcmds", "_rmcoll": "varlistcmds", "rmcoll": "varlistcmds",
    "post": "postfile", "postclose": "postfile", "postutil": "postfile",
    "pause": "progutil", "sleep": "progutil", "window": "progutil", "plugin": "progutil",
    "viewsource": "progutil", "signestimationsample": "estprog", "checkestimationsample": "estprog",
    "makecns": "estprog", "matcproc": "estprog", "mark": "syntax", "markin": "syntax",
    "sortpreserve": "program",
}


def find_help(s: "Session", topic: str) -> Path | None:
    t = "_".join(topic.split()) or "opendta"
    t = ALIASES.get(t, t)
    for ext in (".sthlp", ".hlp"):
        hits = find_file(s, t + ext)
        if hits:
            return hits[0]
    own = HELP_DIR / f"{t}.sthlp"
    if own.exists():
        return own
    full = REGISTRY.get(t)
    if full is not None and (HELP_DIR / f"{full.name}.sthlp").exists():
        return HELP_DIR / f"{full.name}.sthlp"
    return None


def available_topics() -> list[str]:
    return sorted(p.stem for p in HELP_DIR.glob("*.sthlp"))


def _read(path: Path) -> str:
    raw = path.read_bytes()
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


@command("help", "h")
def cmd_help(s: "Session", args: str) -> None:
    topic = args.strip().split(",")[0].strip()
    path = find_help(s, topic)
    if path is None:
        # VERIFY: Stata 14 message for a nonexistent topic
        raise StataError(111, f"help for {topic} not found")
    hook = s.ui_hooks.get("help")
    if hook is not None:
        hook(path, topic or "opendta")
        return
    s.output.write(render_text(_read(path), int(float(s.settings.get("linesize", 80))) - 1), "text")


@command("view")
def cmd_view(s: "Session", args: str) -> None:
    t = args.strip()
    if t.startswith("help "):
        cmd_help(s, t[5:])
        return
    if t.startswith("file "):
        t = t[5:]
    path = Path(strip_outer_quotes(t)).expanduser()
    if not path.exists():
        raise StataError(601, f"file {path} not found")
    hook = s.ui_hooks.get("view")
    if hook is not None:
        hook(path, path.name)
        return
    text = _read(path)
    if path.suffix in (".smcl", ".sthlp", ".hlp") or text.lstrip().startswith("{smcl}"):
        text = render_text(text)
    s.output.write(text if text.endswith("\n") else text + "\n", "text")


@command("search")
def cmd_search(s: "Session", args: str) -> None:
    """Searches for words in the titles and text of the help pages."""
    words = [w.lower() for w in args.split() if not w.startswith(",")]
    if not words:
        raise StataError(198, "nothing to search for")   # VERIFY
    out = s.output
    hits = []
    for path in sorted(HELP_DIR.glob("*.sthlp")):
        text = render_text(_read(path)).lower()
        if all(w in text for w in words):
            title = next((ln.strip() for ln in render_text(_read(path)).splitlines()
                          if "--" in ln), path.stem)
            hits.append((path.stem, title))
    if not hits:
        out.write(f"no entries found for search on \"{' '.join(words)}\"\n", "text")
        return
    for topic, title in hits:
        out.write(f"  [{topic}]  ", "result")
        out.write(title + "\n", "text")
