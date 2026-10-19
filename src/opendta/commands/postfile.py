"""postfile, post, postclose and postutil ([P] postfile).

    postfile postname newvarlist using filename [, every(#) replace]
    post postname (exp) (exp) ... (exp)
    postclose postname

The rows posted are kept in memory and written to the .dta file by postclose; the file is created
(empty) by postfile, so that `replace` and an existing file are checked as in Stata. The data in
memory are not touched. A postfile named by a tempname is forgotten when the program ends.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..core.dataset import NAME_RE, Dataset, Variable
from ..core.errors import StataError
from ..io.dta import write_dta
from ..lang.syntax import match_options, parse_standard
from .files import dta_path
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

_TYPES = ("byte", "int", "long", "float", "double")


@dataclass
class _Post:
    path: Path
    names: list[str]
    types: list[str]
    rows: list[list] = field(default_factory=list)


def _posts(s: "Session") -> dict[str, _Post]:
    if not hasattr(s, "posts"):
        s.posts = {}
    return s.posts


def _is_type(w: str) -> bool:
    return w in _TYPES or re.fullmatch(r"str\d+|strL", w) is not None


def _declare(s: "Session", text: str) -> tuple[list[str], list[str]]:
    """'a b str20 c double(d e)' -> names and types (the default type is `set type`)."""
    names: list[str] = []
    types: list[str] = []
    pending = None
    for m in re.finditer(r"(\w+)\(([^)]*)\)|(\S+)", text):
        if m.group(1):
            if not _is_type(m.group(1)):
                raise StataError(198, f"{m.group(1)} invalid type")
            group = [(w, m.group(1)) for w in m.group(2).split()]
        elif _is_type(m.group(3)):
            pending = m.group(3)
            continue
        else:
            group = [(m.group(3), pending or s.default_type())]
            pending = None
        for name, vtype in group:
            if not NAME_RE.match(name):
                raise StataError(198, f"{name} invalid name")
            if name in names:
                raise StataError(110, f"postfile: variable {name} already exists")   # VERIFY: rc
            names.append(name)
            types.append(vtype)
    return names, types


def _dataset(post: _Post) -> Dataset:
    ds = Dataset()
    ds.nobs = len(post.rows)
    for k, (name, vtype) in enumerate(zip(post.names, post.types)):
        col = [row[k] for row in post.rows]
        if vtype.startswith("str"):
            data = np.array(col, dtype=object)
        else:
            data = np.array(col, dtype=np.float64)
        ds.vars.append(Variable(name, vtype, data))
    return ds


@command("postfile")
def cmd_postfile(s: "Session", args: str) -> None:
    p = parse_standard(args)
    words = p.varlist.split(None, 1)
    if len(words) < 2 or p.using is None:
        raise StataError(198, "invalid syntax")
    postname, decl = words
    if not NAME_RE.match(postname):
        raise StataError(198, f"{postname} invalid name")
    o = match_options(p.options, {"every": 5, "replace": 7}) if p.options.strip() else {}
    if o.get("every") not in (None, True):
        int(float(s.eval(str(o["every"]))))
    posts = _posts(s)
    if postname in posts:
        raise StataError(110, f"postfile {postname} already exists")      # VERIFY: return code
    if re.search(r"(^|\s)_all(\s|$)", decl):
        raise StataError(101, "_all not allowed")
    names, types = _declare(s, decl)
    path = dta_path(p.using)
    if path.exists() and not o.get("replace"):
        raise StataError(602, f"file {path} already exists")
    if not path.exists() and o.get("replace"):
        s.output.write(f"(note: file {path} not found)\n", "text")
    post = _Post(path, names, types)
    try:
        write_dta(_dataset(post), path)
    except OSError as e:
        raise StataError(603, f"file {path} could not be opened ({e.strerror})")
    posts[postname] = post


def _groups(text: str) -> list[str]:
    """The contents of the parenthesized expressions '(a) (b(c))' -> ['a', 'b(c)']."""
    out, i, n = [], 0, len(text)
    while i < n:
        if text[i].isspace():
            i += 1
            continue
        if text[i] != "(":
            raise StataError(198, "invalid syntax")
        depth, j, start = 0, i, i + 1
        while j < n:
            c = text[j]
            if text.startswith('`"', j):
                k = text.find("\"'", j)
                j = n if k == -1 else k + 1
            elif c == '"':
                k = text.find('"', j + 1)
                j = n if k == -1 else k
            elif c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if depth != 0:
            raise StataError(198, "invalid syntax")
        out.append(text[start:j])
        i = j + 1
    return out


@command("post")
def cmd_post(s: "Session", args: str) -> None:
    postname, _, rest = args.strip().partition(" ")
    post = _posts(s).get(postname)
    if post is None:
        raise StataError(111, f"post {postname} not found")              # VERIFY: rc
    try:
        exps = _groups(rest)
        if len(exps) != len(post.names):
            found = (f"only {len(exps)}" if len(exps) < len(post.names)
                     else f"more than {len(post.names)}")
            raise StataError(198, f"post:  {len(post.names)} expressions expected and "
                                  f"{found} found")                      
        row = []
        for k, (name, vtype, text) in enumerate(zip(post.names, post.types, exps), 1):
            try:
                v = s.eval(text)
                if vtype.startswith("str"):
                    if not isinstance(v, str):
                        raise StataError(109, "type mismatch")
                    row.append(v if vtype == "strL" else v[:int(vtype[3:])])
                else:
                    if isinstance(v, str):
                        raise StataError(109, "type mismatch")
                    row.append(float(v))
            except StataError as e:
                raise StataError(e.rc, f"{e.message}\npost:  above message corresponds to "
                                       f"expression {k}, variable {name}")
    except StataError:
        # an error in post closes the postfile (the rows posted so far are saved)
        _close(s, postname)
        raise
    post.rows.append(row)


def _close(s: "Session", name: str) -> None:
    post = _posts(s).pop(name)
    try:
        write_dta(_dataset(post), post.path)
    except OSError as e:
        raise StataError(603, f"file {post.path} could not be opened ({e.strerror})")


@command("postclose")
def cmd_postclose(s: "Session", args: str) -> None:
    name = args.strip()
    if name not in _posts(s):
        raise StataError(111, f"post {name} not found")                  # VERIFY: rc
    _close(s, name)


@command("postutil")
def cmd_postutil(s: "Session", args: str) -> None:
    sub = args.strip()
    posts = _posts(s)
    if sub == "dir":
        for name in posts:
            s.output.write(f"{name}\n", "text")
    elif sub == "clear":
        for name in list(posts):
            _close(s, name)
    else:
        raise StataError(198, "invalid syntax")
