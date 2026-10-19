"""Builds the ICD code tables used by icd9 and icd10 from the files published by the
U.S. Centers for Medicare & Medicaid Services (CMS), which are in the public
domain:

- ICD-9-CM diagnosis codes, version 29 master descriptions (tab-separated:
  code, long description, short description), e.g. the copy of
  cmsv29_master_descriptions in https://github.com/drobbins/ICD9 (icd9.txt);
- ICD-10-CM 2018 order file (icd10cm_order_2018.txt, fixed columns: order,
  code, header flag, short description, long description), e.g. the copy in
  https://github.com/kamillamagna/ICD-10-CSV (codes.txt).

    python tools/make_icd_tables.py ICD9_FILE ICD10_ORDER_FILE

writes src/opendta/data/icd/icd9cm.tsv.gz and icd10cm.tsv.gz: one line per code,
code (no period) <TAB> 1 if billable else 0 <TAB> short description <TAB> long
description.
"""

from __future__ import annotations

import csv
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "src" / "opendta" / "data" / "icd"


def icd9(path: Path) -> list[tuple[str, str, str, str]]:
    rows = []
    with open(path, encoding="latin-1", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)                                  # header
        for rec in reader:
            if not rec or not rec[0].strip():
                continue
            code, long_, short = rec[0].strip(), rec[1].strip(), rec[2].strip()
            rows.append((code, "1", short, long_))
    return rows


def icd10(path: Path) -> list[tuple[str, str, str, str]]:
    rows = []
    with open(path, encoding="latin-1") as f:
        for line in f:
            if len(line) < 16:
                continue
            code = line[6:13].strip()
            billable = line[14:15]
            short = line[16:76].strip()
            long_ = line[77:].strip()
            rows.append((code, billable, short, long_))
    return rows


def write(rows, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT / name, "wt", encoding="utf-8", compresslevel=9) as f:
        for r in rows:
            f.write("\t".join(x.replace("\t", " ") for x in r) + "\n")
    print(f"{name}: {len(rows):,} codes")


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    write(icd9(Path(sys.argv[1])), "icd9cm.tsv.gz")
    write(icd10(Path(sys.argv[2])), "icd10cm.tsv.gz")


if __name__ == "__main__":
    main()
