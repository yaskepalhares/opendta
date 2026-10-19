# ICD code tables

Used by the `icd9` and `icd10` commands (see `src/opendta/core/icd.py`).

| File | Content | Source |
|---|---|---|
| `icd9cm.tsv.gz` | ICD-9-CM diagnosis codes, version 29 (14,567 codes) | CMS, *ICD-9-CM Diagnosis and Procedure Codes: Abbreviated and Full Code Titles* (cmsv29_master_descriptions) |
| `icd10cm.tsv.gz` | ICD-10-CM, fiscal year 2018 (94,127 codes, with categories) | CMS, *2018 ICD-10-CM* code descriptions in tabular order (icd10cm_order_2018.txt) |

Both lists were published by the U.S. Centers for Medicare & Medicaid Services and are in the public domain. Each line holds the code without the period, 1 if the code is billable (0 for a category header), the short description and the long description, separated by tabs. `tools/make_icd_tables.py` rebuilds the files from the CMS originals.

The ICD-9-CM procedure codes (`icd9p`) are not included: `icd9p` checks codes by their format only.
