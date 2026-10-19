"""ICD code tables for icd9 and icd10: ICD-9-CM diagnosis codes (version 29) and
ICD-10-CM (2018), from the public-domain files of the U.S. Centers for Medicare &
Medicaid Services (see src/opendta/data/icd/README.md and
tools/make_icd_tables.py). Codes are kept without the period."""

from __future__ import annotations

import gzip
from functools import lru_cache
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data" / "icd"
FILES = {"icd9": "icd9cm.tsv.gz", "icd10": "icd10cm.tsv.gz"}
VERSIONS = {"icd9": "ICD-9-CM diagnosis codes, version 29 (CMS)", "icd10": "ICD-10-CM, 2018 (CMS)"}


@lru_cache(maxsize=None)
def table(kind: str) -> dict[str, tuple[bool, str, str]]:
    """{code: (billable, short description, long description)}."""
    out = {}
    with gzip.open(DATA / FILES[kind], "rt", encoding="utf-8") as f:
        for line in f:
            code, bill, short, long_ = line.rstrip("\n").split("\t")
            out[code] = (bill == "1", short, long_)
    return out


@lru_cache(maxsize=None)
def _prefixes(kind: str) -> frozenset[str]:
    """Every leading part of a code: categories and subcategories of the codes in the table
    (the ICD-9-CM file lists only the most detailed codes)."""
    pre = set()
    for code in table(kind):
        for k in range(3, len(code) + 1):
            pre.add(code[:k])
    return frozenset(pre)


def norm(code: str) -> str:
    return code.strip().upper().replace(".", "")


def defined(kind: str, code: str) -> bool:
    return norm(code) in _prefixes(kind)


def describe(kind: str, code: str, long: bool = False) -> str:
    """Description of a code in the table; empty for a category the table does not list (the
    ICD-9-CM file has no titles for its categories)."""
    entry = table(kind).get(norm(code))
    return "" if entry is None else entry[2 if long else 1]


def dotted(kind: str, code: str) -> str:
    """The code with its period: after the 3rd character (4th for ICD-9 E codes)."""
    c = norm(code)
    head = 4 if kind == "icd9" and c.startswith("E") else 3
    return c if len(c) <= head else c[:head] + "." + c[head:]
