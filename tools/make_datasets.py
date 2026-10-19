"""Builds the example datasets shipped with OpenDTA (sysuse, webuse).

They are synthetic data made by this script with a fixed seed, so anyone can
rebuild them and see where every number comes from; no third-party dataset is
copied. Run from the repository root:

    python tools/make_datasets.py

writes src/opendta/data/*.dta.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from opendta.core.dataset import Dataset, Variable   # noqa: E402
from opendta.io.dta import write_dta                  # noqa: E402

OUT = ROOT / "src" / "opendta" / "data"


def _ds(label: str, cols, value_labels=None) -> Dataset:
    ds = Dataset()
    ds.set_obs(len(cols[0][2]))
    for name, vtype, values, lab, *rest in cols:
        fmt = rest[0] if rest else ""
        vl = rest[1] if len(rest) > 1 else ""
        ds.add(Variable(name, vtype, np.asarray(values, dtype=object if vtype.startswith("str") else float),
                        fmt=fmt, label=lab, value_label=vl))
    ds.label = label
    ds.value_labels = value_labels or {}
    return ds


def diet(rng) -> Dataset:
    """Adults' diet and body measures (cross-sectional survey with strata and PSUs)."""
    n = 400
    sex = rng.integers(0, 2, n)
    age = rng.integers(18, 80, n)
    region = rng.integers(1, 5, n)
    upf = np.clip(rng.normal(30 + 5 * (region == 1) - 0.1 * (age - 45), 10, n), 2, 80)
    fv = np.clip(rng.poisson(3 + (sex == 1) + 0.02 * (age - 18), n), 0, 12)
    bmi = 22 + 0.08 * upf - 0.25 * fv + 0.04 * (age - 18) + rng.normal(0, 3.5, n)
    obese = (bmi >= 30).astype(float)
    strata = region
    psu = rng.integers(1, 9, n) + 10 * region
    w = np.round(rng.uniform(500, 3000, n))
    return _ds("Diet and body measures (synthetic, OpenDTA)", [
        ("id", "int", np.arange(1, n + 1), "Respondent"),
        ("sex", "byte", sex, "Sex", "%8.0g", "sexlbl"),
        ("age", "byte", age, "Age (years)"),
        ("region", "byte", region, "Region", "%9.0g", "reglbl"),
        ("upf", "float", np.round(upf, 1), "Ultra-processed foods (% of energy)", "%5.1f"),
        ("fruitveg", "byte", fv, "Servings of fruit and vegetables per day"),
        ("bmi", "float", np.round(bmi, 1), "Body mass index (kg/m2)", "%5.1f"),
        ("obese", "byte", obese, "BMI >= 30", "%8.0g", "yesno"),
        ("strata", "byte", strata, "Sampling stratum"),
        ("psu", "int", psu, "Primary sampling unit"),
        ("weight", "float", w, "Sampling weight"),
    ], {"sexlbl": {0: "Male", 1: "Female"}, "reglbl": {1: "North", 2: "South", 3: "East", 4: "West"},
        "yesno": {0: "No", 1: "Yes"}})


def cohort(rng) -> Dataset:
    """Follow-up of a cohort: exposure, person-years and events."""
    n = 300
    exposed = rng.integers(0, 2, n)
    agegrp = rng.integers(1, 4, n)
    rate = 0.02 * np.exp(0.6 * exposed + 0.5 * (agegrp - 1))
    t = rng.exponential(1 / rate)
    cens = rng.uniform(2, 10, n)
    time = np.minimum(t, cens)
    died = (t <= cens).astype(float)
    return _ds("Cohort follow-up (synthetic, OpenDTA)", [
        ("id", "int", np.arange(1, n + 1), "Subject"),
        ("exposed", "byte", exposed, "Exposed", "%8.0g", "yesno"),
        ("agegrp", "byte", agegrp, "Age group", "%9.0g", "agelbl"),
        ("time", "float", np.round(time, 2), "Years of follow-up", "%6.2f"),
        ("died", "byte", died, "Died", "%8.0g", "yesno"),
    ], {"yesno": {0: "No", 1: "Yes"}, "agelbl": {1: "40-49", 2: "50-59", 3: "60-69"}})


def growth(rng) -> Dataset:
    """Height of children measured every year (panel)."""
    kids, waves = 50, 6
    id_ = np.repeat(np.arange(1, kids + 1), waves)
    age = np.tile(np.arange(2, 2 + waves), kids)
    girl = np.repeat(rng.integers(0, 2, kids), waves)
    base = np.repeat(rng.normal(86, 3, kids), waves)
    slope = np.repeat(rng.normal(6.5, 0.6, kids), waves)
    height = base + slope * (age - 2) - 1.0 * girl + rng.normal(0, 0.8, kids * waves)
    return _ds("Children's growth (synthetic, OpenDTA)", [
        ("id", "int", id_, "Child"),
        ("age", "byte", age, "Age (years)"),
        ("girl", "byte", girl, "Girl", "%8.0g", "yesno"),
        ("height", "float", np.round(height, 1), "Height (cm)", "%5.1f"),
    ], {"yesno": {0: "No", 1: "Yes"}})


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for i, (name, fn) in enumerate((("diet", diet), ("cohort", cohort), ("growth", growth))):
        ds = fn(np.random.default_rng(2026 + i))
        write_dta(ds, OUT / f"{name}.dta")
        print(f"{name}.dta: {ds.nobs} obs, {ds.nvars} vars")


if __name__ == "__main__":
    main()
