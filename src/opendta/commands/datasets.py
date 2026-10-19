"""Example datasets, medical codes, ODBC and business calendars: sysuse, webuse
([D] sysuse, webuse), icd9, icd9p and icd10 ([D] icd), odbc ([D] odbc) and bcal
([D] bcal).

sysuse loads OpenDTA's own example datasets (synthetic data built by
tools/make_datasets.py; Stata's example datasets are StataCorp files and are not
included) and datasets in the adopath. webuse downloads a dataset from a URL
prefix (webuse set), by default the public address of Stata's example datasets.
icd9 and icd10 check codes against the public-domain ICD-9-CM (diagnoses) and
ICD-10-CM tables of the U.S. CMS (core/icd.py) and give their descriptions;
icd9p checks procedure codes by their format only. odbc works through the
pyodbc package when it is installed. Messages not shown in the manual are marked
VERIFY and checked in compat/do/1501_datamgmt.do."""

from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Dataset, Variable
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import touse
from .data import _note
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
DATA_DIR = Path(__file__).resolve().parent.parent / "data"
# the public web address of the example datasets of the target version ([D] webuse); the files are
# downloaded by the user on demand, not distributed with OpenDTA
WEBUSE_DEFAULT = "https://www.stata-press.com/data/r14/"
STATA_EXAMPLES = {"auto", "auto2", "census", "lifeexp", "nlsw88", "uslifeexp", "bplong", "bpwide",
                  "cancer", "citytemp", "voter", "pop2000", "sp500", "xtline1", "nhanes2"}


def _split(args: str) -> tuple[list[str], str]:
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    return [strip_outer_quotes(w) for w in re.findall(r'"[^"]*"|\S+', body)], opt


def _use(s: "Session", path: Path, clear: bool) -> None:
    s.interp._execute_expanded(f'use "{path}"' + (", clear" if clear else ""))


# ---------------------------------------------------------------------------
# sysuse, webuse
# ---------------------------------------------------------------------------

def _sys_dirs(s: "Session") -> list[Path]:
    from ..lang.adopath import adopath, resolve_entry
    return [DATA_DIR] + [resolve_entry(s, e) for e in adopath(s)]


def _find_sys(s: "Session", name: str) -> Path | None:
    from ..lang.adopath import find_file
    fname = name if name.endswith(".dta") else name + ".dta"
    if (DATA_DIR / fname).exists():
        return DATA_DIR / fname
    hits = find_file(s, fname)
    return hits[0] if hits else None


@command("sysuse", "sysuse")
def cmd_sysuse(s: "Session", args: str) -> None:
    words, opt = _split(args)
    if not words:
        raise StataError(198, "invalid syntax")
    out = s.output
    if words[0] == "dir":
        match_options(opt, {"all": 3})
        names = sorted({f.name for d in _sys_dirs(s) if d.is_dir() for f in d.glob("*.dta")})
        width = max((len(n) for n in names), default=0) + 3
        per = max(1, 76 // width)
        for i in range(0, len(names), per):
            out.write("  " + "".join(f"{n:<{width}}" for n in names[i:i + per]).rstrip() + "\n", "text")
        return
    if words[0] == "describe":
        if len(words) != 2:
            raise StataError(198, "invalid syntax")
        f = _find_sys(s, words[1])
        if f is None:
            raise _not_found(words[1])
        from ..io.dta import read_dta
        from .inspect import cmd_describe
        saved = s.data
        try:
            s.data = read_dta(f)
            s.data.filename = f.name
            cmd_describe(s, "")
        finally:
            s.data = saved
        return
    o = match_options(opt, {"clear": 5})
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    f = _find_sys(s, words[0])
    if f is None:
        raise _not_found(words[0])
    _use(s, f, bool(o.get("clear")))


def _not_found(name: str) -> StataError:
    base = name[:-4] if name.endswith(".dta") else name
    if base in STATA_EXAMPLES:
        # Stata's example datasets are StataCorp files; OpenDTA has its own (sysuse dir)
        return StataError(601, f"file {base}.dta not found\n{base}.dta is one of Stata's example datasets, "
                               f"which OpenDTA does not include; webuse {base} downloads it from Stata's "
                               "website, and sysuse dir lists OpenDTA's own")
    return StataError(601, f"file {base}.dta not found")


@command("webuse", "webuse")
def cmd_webuse(s: "Session", args: str) -> None:
    words, opt = _split(args)
    out = s.output
    if words and words[0] == "set":
        url = words[1] if len(words) > 1 else WEBUSE_DEFAULT
        if not url.endswith("/"):
            url += "/"
        s.settings["webuse_url"] = url
        out.write(f'(prefix now "{url}")\n', "text")
        return
    if words and words[0] == "query":
        out.write(f'(prefix is "{s.settings.get("webuse_url", WEBUSE_DEFAULT)}")\n', "text")
        return
    o = match_options(opt, {"clear": 5})
    if len(words) != 1:
        raise StataError(198, "invalid syntax")
    name = words[0] if words[0].endswith(".dta") else words[0] + ".dta"
    url = s.settings.get("webuse_url", WEBUSE_DEFAULT) + name
    import urllib.request
    from ..updater import ssl_context
    try:
        with urllib.request.urlopen(url, timeout=30, context=ssl_context()) as r:
            data = r.read()
    except Exception:   # noqa: BLE001
        raise StataError(601, f"file {url} not found")
    tmp = Path(tempfile.mkdtemp(prefix="opendta-webuse-")) / name
    tmp.write_bytes(data)
    _use(s, tmp, bool(o.get("clear")))
    s.data.filename = url


# ---------------------------------------------------------------------------
# icd9, icd9p, icd10
# ---------------------------------------------------------------------------

_PROBLEMS = {
    # icd10 follows ICD-10-CM: up to 7 characters, letters allowed after the 2nd (VERIFY: Stata's
    # icd10 follows the WHO ICD-10, whose codes have at most 4)
    "icd10": ["Invalid placement of period", "Too many periods", "Code too short", "Code too long",
              "Invalid 1st char (not A-Z)", "Invalid 2nd char (not 0-9)", "Invalid 3rd char (not 0-9 or A-Z)",
              "Invalid 4th-7th char (not 0-9 or A-Z)", "Code not defined"],
    "icd9": ["Invalid placement of period", "Too many periods", "Code too short", "Code too long",
             "Invalid 1st char (not 0-9, E, or V)", "Invalid 2nd char (not 0-9)", "Invalid 3rd char (not 0-9)",
             "Invalid 4th char (not 0-9)", "Invalid 5th char (not 0-9)", "Code not defined"],
    "icd9p": ["Invalid placement of period", "Too many periods", "Code too short", "Code too long",
              "Invalid 1st char (not 0-9)", "Invalid 2nd char (not 0-9)", "Invalid 3rd char (not 0-9)",
              "Invalid 4th char (not 0-9)"],
}


def icd_problem(kind: str, code: str, fmtonly: bool = False) -> int:
    """0 if the code is valid, otherwise the number of the problem (1-based, as listed in
    _PROBLEMS). icd9 and icd10 also check that the code is defined in the public-domain
    ICD-9-CM and ICD-10-CM tables (core/icd.py), unless fmtonly; icd9p checks the format only."""
    p = _format_problem(kind, code)
    if p or fmtonly or kind == "icd9p":
        return p
    from ..core import icd
    return 0 if icd.defined(kind, code) else len(_PROBLEMS[kind])


def _format_problem(kind: str, code: str) -> int:
    c = code.strip().upper()
    if c.count(".") > 1:
        return 2
    body = c.replace(".", "")
    head = {"icd10": 3, "icd9p": 2}.get(kind, 4 if c.startswith("E") else 3)
    if "." in c and c.index(".") != head:
        return 1
    maxlen = {"icd10": 7, "icd9": 5, "icd9p": 4}[kind]
    minlen = 3 if kind != "icd9p" else 2
    if len(body) < minlen:
        return 3
    if len(body) > maxlen:
        return 4
    if kind == "icd10":
        if not ("A" <= body[0] <= "Z"):
            return 5
        if not body[1].isdigit():
            return 6
        if not body[2].isalnum():
            return 7
        if any(not ch.isalnum() for ch in body[3:]):
            return 8
        return 0
    if kind == "icd9":
        if not (body[0].isdigit() or body[0] in "EV"):
            return 5
        for i, ch in enumerate(body[1:], start=2):
            if not ch.isdigit():
                return 4 + i
        return 0
    for i, ch in enumerate(body, start=1):
        if not ch.isdigit():
            return 4 + i
    return 0


def _icd_info(s: "Session", kind: str, sub: str, rest: str) -> None:
    """icd9/icd10 lookup (codes, ranges a/b and wildcards), search (words in the descriptions)
    and query (the code list in use)."""
    from ..core import icd
    out = s.output
    if sub == "query":
        out.write(f"\n    {kind} codes: {icd.VERSIONS[kind]}, {len(icd.table(kind)):,} codes\n", "text")
        return
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    o = match_options(opt, {"or": 2, "long": 1, "version": 1})
    t = icd.table(kind)
    codes = sorted(t)
    if sub == "lookup":
        hits = []
        for item in body.split():
            it = icd.norm(item)
            if "/" in it:
                a, b = it.split("/", 1)
                hits += [k for k in codes if a <= k[:len(a)] and k[:len(b)] <= b]
            elif it.endswith("*"):
                hits += [k for k in codes if k.startswith(it[:-1])]
            else:
                hits += [k for k in codes if k == it] or ([it] if icd.defined(kind, it) else [])
        if not hits:
            out.write("(no matches found)\n", "text")   # VERIFY
    else:
        words = [strip_outer_quotes(w).lower() for w in re.findall(r'"[^"]*"|\S+', body)]
        if not words:
            raise StataError(198, "invalid syntax")
        test = any if o.get("or") else all
        hits = [k for k in codes if test(w in t[k][2].lower() for w in words)]
        if not hits:
            out.write("(no matches found)\n", "text")
    for k in dict.fromkeys(hits):
        out.write(f"{icd.dotted(kind, k):>10}  ", "result")
        out.write(icd.describe(kind, k, long=bool(o.get("long")) or sub == "search") + "\n", "text")
    s.r = {"N_codes": float(len(dict.fromkeys(hits)))}


def _icd(s: "Session", kind: str, args: str) -> None:
    words = args.split(None, 1)
    if not words:
        raise StataError(198, "invalid syntax")
    sub, rest = words[0], (words[1] if len(words) > 1 else "")
    if sub in ("lookup", "search", "query"):
        if kind == "icd9p":
            raise StataError(198, "icd9p lookup, search and query need the ICD-9-CM procedure code list, which "
                                  "OpenDTA does not include; check, clean and generate work by the format")
        _icd_info(s, kind, sub, rest)
        return
    p = parse_standard(rest)
    ds = s.data
    out = s.output
    if sub == "check":
        o = match_options(p.options, {"any": 3, "list": 1, "generate": 3, "summary": 3, "version": 1,
                                      "fmtonly": 3, "system": 3})
        fmtonly = bool(o.get("any") or o.get("fmtonly"))   # VERIFY: meaning of any in Stata 14
        name = s.expand_varlist(p.varlist)[0]
        v = ds.get(name)
        if not v.is_string:
            raise StataError(109, f"{name} must be a string variable")   # VERIFY
        mask = touse(s, p)
        codes = [str(x) for x in v.raw]
        probs = np.array([icd_problem(kind, c, fmtonly) if m and c.strip() else 0 for c, m in zip(codes, mask)])
        nmiss = sum(1 for c, m in zip(codes, mask) if m and not c.strip())
        bad = int((probs > 0).sum())
        # layout of Stata 14 (compat 1501; VERIFY: the wording with missing values and for
        # icd10/icd9p)
        tail = "no missing values" if not nmiss else f"{nmiss:,} missing value{'s' if nmiss != 1 else ''}"
        out.write(f"({name} contains {tail})\n", "text")
        if bad == 0:
            out.write(f"\n{name} contains valid codes\n", "text")   # VERIFY
        else:
            out.write(f"\n{name} contains invalid codes:\n\n", "text")
            for i, lab in enumerate(_PROBLEMS[kind], start=1):
                out.write(f"{i:>5}.  {lab:<40}", "text")
                out.write(f"{int((probs == i).sum()):>11,}\n", "result")
            out.write(" " * 48 + "-" * 11 + "\n", "text")
            out.write(f"{'':8}{'Total':<40}", "text")
            out.write(f"{bad:>11,}\n", "result")
        if isinstance(o.get("generate"), str):
            ds.add(Variable(str(o["generate"]), "byte", probs.astype(float)))
            s.notify_state()
        if o.get("list") and bad:
            # a list of the code and the problem without its parenthesis (compat 1501)
            from .inspect import cmd_list
            out.write("\nListing of invalid and undefined codes\n", "text")
            short = np.array([_PROBLEMS[kind][k - 1].split(" (")[0] if k else "" for k in probs], dtype=object)
            w = max(len(x) for x in short)
            ds.add(Variable("__prob", f"str{w}", short, fmt=f"%-{w}s"))
            try:
                cmd_list(s, f'{name} __prob if __prob != ""')
            finally:
                ds.drop_vars(["__prob"])
        if o.get("summary"):
            vals, counts = np.unique([c.strip().upper() for c, m in zip(codes, mask) if m and c.strip()],
                                     return_counts=True)
            out.write(f"\n{name:>12} |      Freq.\n" + "-" * 13 + "+" + "-" * 11 + "\n", "text")
            for c, n in zip(vals, counts):
                out.write(f"{c:>12} |", "text")
                out.write(f"{n:>11,}\n", "result")
        s.r = {"esum": float(bad)}
        for i in range(len(_PROBLEMS[kind])):
            s.r[f"e{i + 1}"] = float((probs == i + 1).sum())
        return
    if sub == "clean":
        # icd9 and icd9p clean take only dots and pad in Stata 14 (compat 1501)
        allowed = {"generate": 3, "replace": 7, "dots": 3, "pad": 3, "check": 5} if kind == "icd10" \
            else {"dots": 3, "pad": 3}
        o = match_options(p.options, allowed)
        name = s.expand_varlist(p.varlist)[0]
        v = ds.get(name)
        if not v.is_string:
            raise StataError(109, f"{name} must be a string variable")   # VERIFY
        if kind != "icd10" and any(str(c).strip() and _format_problem(kind, str(c)) for c in v.raw):
            # clean refuses badly formed codes (compat 1501; VERIFY: icd9p wording)
            raise StataError(459, f"{name} contains invalid ICD-9 codes")
        head = {"icd10": 3, "icd9p": 2}.get(kind, 3)
        new_vals, changes = [], 0
        for c in v.raw:
            t = str(c).strip().upper().replace(" ", "").replace(".", "")
            if o.get("dots") and len(t) > head:
                h = 4 if kind == "icd9" and t.startswith("E") else head
                t = t[:h] + "." + t[h:]
            changes += t != c
            new_vals.append(t)
        if isinstance(o.get("generate"), str):
            ds.add(Variable(str(o["generate"]), v.vtype, np.array(new_vals, dtype=object)))
        else:
            v.data = new_vals
            if kind != "icd10":
                # Stata 14 leaves the cleaned variable left-justified in list (compat 1501; VERIFY: the width)
                v.fmt = "%-" + v.fmt.lstrip("%-")
        _note(s, f"({changes:,} change{'s' if changes != 1 else ''} made)")   # VERIFY
        s.notify_state()
        return
    if sub in ("generate", "gen"):
        o = match_options(p.options, {"category": 3, "short": 2, "range": 2, "description": 4,
                                      "addcode": 3, "nodots": 3, "long": 2, "end": 2, "version": 1})
        if not p.exp:
            raise StataError(198, "invalid syntax")
        newname = p.varlist.strip()
        src = ds.get(s.expand_varlist(p.exp.strip())[0])
        codes = [str(c).strip().upper().replace(".", "") for c in src.raw]
        if o.get("description"):
            if kind == "icd9p":
                raise StataError(198, "icd9p generate, description needs the ICD-9-CM procedure code list, "
                                      "which OpenDTA does not include")
            from ..core import icd
            where = str(o["addcode"]) if isinstance(o.get("addcode"), str) else ""
            vals = []
            for c in codes:
                d = icd.describe(kind, c, long=bool(o.get("long"))) if c and icd.defined(kind, c) else ""
                if d and where.startswith("b"):
                    d = f"{icd.dotted(kind, c) if not o.get('nodots') else c} {d}"
                elif d and where.startswith("e"):
                    d = f"{d} {icd.dotted(kind, c) if not o.get('nodots') else c}"
                vals.append(d)
            ds.add(Variable(newname, f"str{max([len(x) for x in vals] + [1])}", np.array(vals, dtype=object)))
        elif o.get("category"):
            k = {"icd10": 3, "icd9p": 2}.get(kind, 3)
            vals = [c[:4] if kind == "icd9" and c.startswith("E") else c[:k] for c in codes]
            ds.add(Variable(newname, f"str{max([len(x) for x in vals] + [1])}", np.array(vals, dtype=object)))
        elif o.get("short"):
            vals = [c[:4] for c in codes]
            ds.add(Variable(newname, f"str{max([len(x) for x in vals] + [1])}", np.array(vals, dtype=object)))
        elif isinstance(o.get("range"), str):
            ranges = []
            for part in str(o["range"]).split():
                a, _, b = part.upper().replace(".", "").partition("/")
                ranges.append((a, b or a))

            def inside(c: str) -> bool:
                return any(a <= c[:len(a)] <= b or (c.startswith(a) and not b) or (a <= c <= b + "~")
                           for a, b in ranges) if c else False
            vals = np.array([1.0 if inside(c) else (SYS if not c else 0.0) for c in codes])
            ds.add(Variable(newname, "byte", vals))
        else:
            raise StataError(198, "category, short or range() required")   # VERIFY
        s.notify_state()
        return
    raise StataError(198, f"{kind} {sub}: invalid subcommand")


@command("icd10", "icd10")
def cmd_icd10(s: "Session", args: str) -> None:
    _icd(s, "icd10", args)


@command("icd9", "icd9")
def cmd_icd9(s: "Session", args: str) -> None:
    _icd(s, "icd9", args)


@command("icd9p", "icd9p")
def cmd_icd9p(s: "Session", args: str) -> None:
    _icd(s, "icd9p", args)


# ---------------------------------------------------------------------------
# odbc (through pyodbc)
# ---------------------------------------------------------------------------

def _pyodbc():
    try:
        import pyodbc
        return pyodbc
    except ImportError:
        raise StataError(198, "odbc needs the Python package pyodbc and an ODBC driver for your database "
                              "(pip install pyodbc)")


def _connect(s: "Session", o: dict):
    pyodbc = _pyodbc()
    dsn = strip_outer_quotes(str(o["dsn"])) if isinstance(o.get("dsn"), str) else s.settings.get("odbc_dsn")
    conn = strip_outer_quotes(str(o["connectionstring"])) if isinstance(o.get("connectionstring"), str) else None
    if not dsn and not conn:
        raise StataError(198, "dsn() or connectionstring() required")
    try:
        return pyodbc.connect(conn or f"DSN={dsn}", autocommit=True)
    except Exception as e:   # noqa: BLE001
        raise StataError(682, f"could not connect to {dsn or conn} ({e})")


_ODBC_OPTS = {"dsn": 3, "connectionstring": 4, "table": 1, "exec": 4, "clear": 5, "noquote": 5,
              "lowercase": 5, "sqlshow": 4, "allstring": 4, "datestring": 5, "user": 1, "password": 4}


@command("odbc", "odbc")
def cmd_odbc(s: "Session", args: str) -> None:
    words = args.strip().split(None, 1)
    if not words:
        raise StataError(198, "invalid syntax")
    sub, rest = words[0], (words[1] if len(words) > 1 else "")
    out = s.output
    c = find_top(rest, ",")
    body, opt = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    if sub == "list":
        pyodbc = _pyodbc()
        out.write("\nData Source Name                   Driver\n" + "-" * 79 + "\n", "text")
        for name, drv in sorted(pyodbc.dataSources().items()):
            out.write(f"{name:<35}{drv}\n", "text")
        out.write("-" * 79 + "\n", "text")
        return
    if sub == "query":
        o = match_options(opt, _ODBC_OPTS)
        if body.strip():
            o["dsn"] = body.strip()
        cn = _connect(s, o)
        out.write(f"\nDataSource: {strip_outer_quotes(str(o.get('dsn', '')))}\n" + "-" * 79 + "\n", "text")
        for row in cn.cursor().tables(tableType="TABLE"):
            out.write(f"{row.table_name}\n", "text")
        out.write("-" * 79 + "\n", "text")
        return
    if sub == "describe":
        o = match_options(opt, _ODBC_OPTS)
        table = strip_outer_quotes(body.strip())
        cn = _connect(s, o)
        out.write(f"\nTable:  {table}\n" + "-" * 79 + "\n" + f"{'Variable Name':<35}Variable Type\n"
                  + "-" * 79 + "\n", "text")
        for col in cn.cursor().columns(table=table):
            out.write(f"{col.column_name:<35}{col.type_name}\n", "text")
        out.write("-" * 79 + "\n", "text")
        return
    if sub == "load":
        o = match_options(opt, _ODBC_OPTS)
        if s.data.changed and s.data.nvars and not o.get("clear"):
            raise StataError(4, "no; data in memory would be lost")
        if isinstance(o.get("exec"), str):
            sql = strip_outer_quotes(str(o["exec"]))
        elif isinstance(o.get("table"), str):
            cols = ", ".join(body.split()) if body.strip() else "*"
            sql = f"SELECT {cols} FROM {strip_outer_quotes(str(o['table']))}"
        else:
            raise StataError(198, "table() or exec() required")
        if o.get("sqlshow"):
            out.write(f"SQL: {sql}\n", "text")
        cn = _connect(s, o)
        cur = cn.cursor()
        cur.execute(sql)
        names = [d[0].lower() if o.get("lowercase") else d[0] for d in cur.description]
        rows = cur.fetchall()
        ds = Dataset()
        ds.set_obs(len(rows))
        for j, n in enumerate(names):
            vals = [r[j] for r in rows]
            vname = re.sub(r"\W", "_", n)[:32]
            if not o.get("allstring") and all(v is None or isinstance(v, (int, float)) for v in vals):
                ds.add(Variable(vname, "double", np.array([SYS if v is None else float(v) for v in vals])))
            else:
                texts = ["" if v is None else str(v) for v in vals]
                ds.add(Variable(vname, f"str{min(max([len(t) for t in texts] + [1]), 2045)}",
                                np.array(texts, dtype=object)))
        ds.changed = True
        s.data = ds
        s.notify_state()
        return
    if sub == "exec":
        m = re.match(r'\s*\((.*)\)\s*$', body) or re.match(r"(.*)", body)
        o = match_options(opt, _ODBC_OPTS)
        sql = strip_outer_quotes(m.group(1).strip())
        cn = _connect(s, o)
        cn.cursor().execute(sql)
        return
    raise StataError(198, f"odbc {sub} is not available in OpenDTA")


# ---------------------------------------------------------------------------
# bcal
# ---------------------------------------------------------------------------

@command("bcal", "bcal")
def cmd_bcal(s: "Session", args: str) -> None:
    from ..core import bcal
    from ..core.formats import format_value
    words, opt = _split(args)
    if not words:
        raise StataError(198, "invalid syntax")
    sub = words[0]
    out = s.output
    from ..lang.adopath import adopath, resolve_entry
    for e in adopath(s):
        bcal.add_path(resolve_entry(s, e))

    def td(x):
        return format_value(float(x), "%td", pad=False).strip()
    if sub == "check":
        o = match_options(opt, {"rc0": 3})
        names = words[1:] or sorted({re.match(r"%-?tb(\w+)", v.fmt).group(1) for v in s.data.vars
                                     if re.match(r"%-?tb\w", v.fmt)})
        ok = True
        for n in names:
            try:
                cal = bcal.get(n)
                users = [v.name for v in s.data.vars if re.match(rf"%-?tb{n}\b", v.fmt)]
                out.write(f"{'%tb' + n:>12}:  defined" + (", used by variable" + ("s" if len(users) > 1 else "")
                                                         if users else "") + "\n", "text")
                if users:
                    out.write(" " * 16 + " ".join(users) + "\n", "text")
            except StataError as e:
                ok = False
                out.write(f"{'%tb' + n:>12}:  {e.message}\n", "text")
        if not ok and not o.get("rc0"):
            raise StataError(459, "")   # VERIFY rc
        return
    if sub == "load":
        if len(words) != 2:
            raise StataError(198, "invalid syntax")
        bcal.clear_cache()
        cal = bcal.get(words[1])
        out.write(f"(calendar {cal.path} loaded)\n", "text")   # VERIFY
        return
    if sub == "describe":
        if len(words) != 2:
            raise StataError(198, "invalid syntax")
        cal = bcal.get(words[1])
        n_in = len(cal.days)
        span = cal.hi - cal.lo + 1
        years = span / 365.25
        b_lo, b_hi = -cal.c_index, n_in - 1 - cal.c_index
        out.write(f"\n  Business calendar {cal.name} (format %tb{cal.name}):\n\n", "text")
        if cal.purpose:
            out.write(f"    purpose:  {cal.purpose}\n\n", "text")
        out.write(f"      range:  {td(cal.lo):<11}{td(cal.hi)}\n", "text")
        out.write(f"{cal.lo:>18}{cal.hi:>11}   in %td units\n", "text")
        out.write(f"{b_lo:>18}{b_hi:>11}   in %tb{cal.name} units\n\n", "text")
        out.write(f"     center:  {td(cal.center)}\n{cal.center:>18}   in %td units\n{0:>18}   in %tb{cal.name} "
                  f"units\n\n", "text")
        out.write(f"    omitted:{span - n_in:>11,}   days\n{(span - n_in) / years:>21.1f}   approx. days/year\n\n",
                  "text")
        out.write(f"   included:{n_in:>11,}   days\n{n_in / years:>21.1f}   approx. days/year\n", "text")
        s.r = {"name": cal.name, "purpose": cal.purpose, "min_date_td": float(cal.lo), "max_date_td": float(cal.hi),
               "ctr_date_td": float(cal.center), "omitted": float(span - n_in), "included": float(n_in)}
        return
    if sub == "create":
        o = match_options(opt, {"from": 4, "generate": 3, "purpose": 4, "dateformat": 5, "range": 3,
                                "centerdate": 6, "excludemissing": 4, "maxgap": 4, "replace": 7, "personal": 4})
        if len(words) != 2 or not isinstance(o.get("from"), str):
            raise StataError(198, "invalid syntax")
        path = Path(os.path.expanduser(words[1]))
        if path.suffix != ".stbcal":
            path = path.with_suffix(".stbcal")
        if path.exists() and not o.get("replace"):
            raise StataError(602, f"file {path} already exists")
        v = s.data.get(s.expand_varlist(str(o["from"]))[0])
        d = v.data
        days = np.unique(d[d < SYS]).astype(int)
        if not len(days):
            raise StataError(2000, "no observations")
        lo, hi = int(days.min()), int(days.max())
        if isinstance(o.get("range"), str):
            a, b = str(o["range"]).split()
            lo, hi = _td(a), _td(b)
        present = set(days.tolist())
        import datetime as _dt
        epoch = _dt.date(1960, 1, 1)
        all_days = range(lo, hi + 1)
        dows = {(epoch + _dt.timedelta(days=t)).weekday() for t in present}
        drop_dow = [k for k in range(7) if k not in dows]
        names = ["Mo", "Tu", "We", "Th", "Fr", "Sa", "Su"]
        center = _td(str(o["centerdate"])) if isinstance(o.get("centerdate"), str) else min(present)

        def ymd(t):
            return (epoch + _dt.timedelta(days=t)).strftime("%Y%b%d").lower()
        lines = ["* Business calendar created by bcal create (OpenDTA)", "version 14",
                 f'purpose "{strip_outer_quotes(str(o.get("purpose")))}"' if isinstance(o.get("purpose"), str) else "",
                 "dateformat ymd", f"range {ymd(lo)} {ymd(hi)}", f"centerdate {ymd(center)}"]
        if drop_dow:
            lines.append("omit dayofweek (" + " ".join(names[k] for k in drop_dow) + ")")
        for t in all_days:
            if t not in present and (epoch + _dt.timedelta(days=t)).weekday() not in drop_dow:
                lines.append(f"omit date {ymd(t)}")
        path.write_text("\n".join(x for x in lines if x) + "\n", encoding="utf-8")
        bcal.clear_cache()
        out.write(f"\n  Business calendar {path.stem} (format %tb{path.stem}):\n", "text")
        out.write(f"    file {path} saved\n", "text")   # VERIFY
        if isinstance(o.get("generate"), str):
            cal = bcal.get(path.stem) if path.parent.resolve() == Path.cwd().resolve() else \
                bcal.parse(path.stem, path.read_text(encoding="utf-8"), str(path))
            vals = np.array([cal.bofd(x) if x < SYS else SYS for x in d])
            vals = np.where(np.isnan(vals), SYS, vals)
            s.data.add(Variable(str(o["generate"]), "long", vals, fmt=f"%tb{path.stem}"))
            s.notify_state()
        return
    raise StataError(198, f"bcal {sub}: invalid subcommand")


def _td(text: str) -> int:
    """A date such as 01nov2011 or 2011nov01 in %td units."""
    from ..core import bcal
    for order in ("dmy", "ymd", "mdy"):
        try:
            return (bcal._date(strip_outer_quotes(text), order) - bcal.EPOCH).days
        except (StataError, ValueError):
            continue
    raise StataError(198, f"{text}: invalid date")
