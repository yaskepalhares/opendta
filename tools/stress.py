"""File stress test: how far OpenDTA can go.

Generates synthetic datasets of increasing size and measures, at each step, time,
peak memory (RSS) and file size for `save` and `use` (.dta) and for
`export delimited` / `import delimited` (CSV). Each measurement runs in a
separate subprocess, watched by a memory cap: if it goes over the cap, the
subprocess is killed and the series stops there. The files live in a
temporary folder deleted at the end.

Usage:
    python tools/stress.py --cap-gb 5 --output results.json
    python tools/stress.py --quick

Dataset profiles:
    survey    100 variables per block: 60 byte, 15 int, 5 long, 12 float,
              4 double and 4 str16 (similar to a survey dataset)
    numeric   the same without the strings
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PROFILE = [("byte", 60), ("int", 15), ("long", 5), ("float", 12), ("double", 4), ("str16", 4)]
WORDS = ["Belo Horizonte", "Contagem", "Betim", "Sabará", "Nova Lima", "Ibirité",
         "Santa Luzia", "Vespasiano", "Lagoa Santa", ""]


# ---------------------------------------------------------------------------
# generation (in the subprocess)
# ---------------------------------------------------------------------------

def build(nobs: int, nvars: int, profile: str, seed: int = 1):
    import numpy as np

    from opendta.core import missing as M
    from opendta.core.dataset import Dataset, Variable

    rng = np.random.default_rng(seed)
    ds = Dataset()
    ds.nobs = nobs
    spec = [(t, n) for t, n in PROFILE if not (profile == "numeric" and t.startswith("str"))]
    total = sum(n for _, n in spec)
    counts = [int(nvars * n / total) for _, n in spec]
    counts[0] += nvars - sum(counts)              # the rest goes to byte
    # interleaves the types, as in a real dataset
    order = []
    left = counts[:]
    while len(order) < nvars:
        for j, (vtype, _) in enumerate(spec):
            if left[j]:
                order.append(vtype)
                left[j] -= 1
    for k, vtype in enumerate(order):
        name = f"v{k + 1}"
        if vtype == "byte":
            data = rng.integers(1, 10, nobs).astype(np.float64)
        elif vtype == "int":
            data = rng.integers(-30000, 30000, nobs).astype(np.float64)
        elif vtype == "long":
            data = rng.integers(-2_000_000_000, 2_000_000_000, nobs).astype(np.float64)
        elif vtype == "float":
            data = rng.normal(25, 5, nobs).astype(np.float32).astype(np.float64)
        elif vtype == "double":
            data = rng.normal(0, 1, nobs)
        else:
            data = np.array(WORDS, dtype=object)[rng.integers(0, len(WORDS), nobs)]
        if not vtype.startswith("str") and nobs:
            data[rng.random(nobs) < 0.05] = M.SYSMISS
            data[0] = M.EXTENDED["a"]
        ds.vars.append(Variable(name, vtype, data))
    ds.label = f"estresse {nobs}x{nvars}"
    ds.vars[0].label = "first variable"
    ds.value_labels = {"sim": {1: "Sim", 2: "Não"}}
    ds.vars[0].value_label = "sim"
    return ds


def digest(ds) -> str:
    h = hashlib.sha256()
    for v in ds.vars:
        h.update(v.name.encode() + v.vtype.encode())
        if v.is_string:
            h.update("\x1f".join(v.data.tolist()).encode())
        else:
            h.update(v.data.tobytes())
    return h.hexdigest()


def rss_mb() -> float:
    import psutil
    return psutil.Process().memory_info().rss / 2**20


def peak_mb() -> float:
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return kb / 1024 if sys.platform != "darwin" else kb / 2**20


def child(args) -> dict:
    """One measurement. Prints a JSON on the last line."""
    out: dict = {}
    path = Path(args.path)
    if args.phase == "save":
        ds = build(args.nobs, args.nvars, args.profile)
        out["gen_rss_mb"] = rss_mb()
        out["digest"] = digest(ds)
        from opendta.io.dta import write_dta
        t = time.perf_counter()
        write_dta(ds, path)
        out["secs"] = time.perf_counter() - t
        out["file_mb"] = path.stat().st_size / 2**20
        out["release"] = path.read_bytes()[28:31].decode() if path.stat().st_size < 2**31 else "?"
    elif args.phase == "use":
        from opendta.io.dta import read_dta
        base = rss_mb()
        t = time.perf_counter()
        ds = read_dta(path)
        out["secs"] = time.perf_counter() - t
        out["data_rss_mb"] = rss_mb() - base
        out["digest"] = digest(ds)
        out["nobs"], out["nvars"] = ds.nobs, ds.nvars
    elif args.phase == "readstat":
        import pyreadstat
        t = time.perf_counter()
        df, meta = pyreadstat.read_dta(str(path), output_format="dict")
        out["secs"] = time.perf_counter() - t
        out["ok"] = (meta.number_rows == args.nobs and len(df) == args.nvars)
    elif args.phase == "export":
        import numpy as np

        from opendta.io.delimited import write_delimited
        ds = build(args.nobs, args.nvars, args.profile)
        t = time.perf_counter()
        write_delimited(ds, path, ds.names, np.arange(ds.nobs))
        out["secs"] = time.perf_counter() - t
        out["file_mb"] = path.stat().st_size / 2**20
    elif args.phase == "import":
        from opendta.io.delimited import read_delimited
        t = time.perf_counter()
        ds = read_delimited(path)
        out["secs"] = time.perf_counter() - t
        out["nobs"], out["nvars"] = ds.nobs, ds.nvars
    out["peak_mb"] = peak_mb()
    print(json.dumps(out))
    return out


# ---------------------------------------------------------------------------
# orchestration (in the main process)
# ---------------------------------------------------------------------------

def run_step(phase: str, path: Path, nobs: int, nvars: int, profile: str,
             ceiling_mb: float, timeout: float) -> dict:
    import psutil

    cmd = [sys.executable, __file__, "--child", phase, "--path", str(path),
           "--nobs", str(nobs), "--nvars", str(nvars), "--profile", profile]
    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    p = psutil.Process(proc.pid)
    watched_peak = 0.0
    reason = ""
    while proc.poll() is None:
        try:
            rss = p.memory_info().rss / 2**20
        except psutil.NoSuchProcess:
            break
        watched_peak = max(watched_peak, rss)
        if rss > ceiling_mb:
            reason = f"went over the cap of {ceiling_mb:.0f} MB"
            proc.kill()
            break
        if time.perf_counter() - t0 > timeout:
            reason = f"took more than {timeout:.0f} s"
            proc.kill()
            break
        time.sleep(0.02)
    stdout, stderr = proc.communicate()
    if reason:
        return {"phase": phase, "fail": reason, "watched_peak_mb": watched_peak}
    if proc.returncode != 0:
        last = (stderr.strip().splitlines() or ["?"])[-1]
        return {"phase": phase, "fail": f"error: {last}", "watched_peak_mb": watched_peak}
    res = json.loads(stdout.strip().splitlines()[-1])
    res["phase"] = phase
    res["wall_secs"] = time.perf_counter() - t0
    return res


def series(name, steps, profile, phases, tmp, ceiling_mb, timeout, log):
    results = []
    for nobs, nvars in steps:
        row = {"series": name, "nobs": nobs, "nvars": nvars, "profile": profile}
        ext = ".csv" if "export" in phases else ".dta"
        path = tmp / f"s_{nobs}_{nvars}{ext}"
        stop = False
        for phase in phases:
            r = run_step(phase, path, nobs, nvars, profile, ceiling_mb, timeout)
            row[phase] = r
            if "fail" in r:
                stop = True
                break
        if "save" in row and "use" in row and "digest" in row["use"]:
            row["roundtrip_ok"] = row["save"]["digest"] == row["use"]["digest"]
        log(row)
        results.append(row)
        if path.exists():
            path.unlink()
        if stop:
            break
    return results


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--child")
    ap.add_argument("--path")
    ap.add_argument("--nobs", type=int, default=0)
    ap.add_argument("--nvars", type=int, default=0)
    ap.add_argument("--profile", default="survey")
    ap.add_argument("--cap-gb", type=float, default=5.0)
    ap.add_argument("--timeout", type=float, default=900)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--series", default="long,numeric,wide,csv")
    ap.add_argument("--output", default="")
    args = ap.parse_args()
    if args.child:
        args.phase = args.child
        child(args)
        return

    ceiling = args.cap_gb * 1024
    if args.quick:
        long_steps = [(10_000, 100), (100_000, 100)]
        wide_steps = [(1_000, 1_000)]
        csv_steps = [(10_000, 100)]
    else:
        long_steps = [(n, 100) for n in (10_000, 100_000, 250_000, 500_000, 1_000_000,
                                         2_000_000, 3_000_000, 4_000_000, 6_000_000, 8_000_000,
                                         12_000_000, 16_000_000, 20_000_000, 24_000_000)]
        wide_steps = [(1_000, k) for k in (1_000, 5_000, 10_000, 32_767, 40_000)]
        csv_steps = [(n, 100) for n in (10_000, 100_000, 500_000, 1_000_000, 2_000_000, 4_000_000)]

    all_rows = []

    def log(row):
        parts = [f"{row['series']:<9} {row['nobs']:>10,} obs x {row['nvars']:>6,} vars"]
        for ph in ("save", "use", "readstat", "export", "import"):
            r = row.get(ph)
            if not r:
                continue
            if "fail" in r:
                parts.append(f"{ph}: FAILED ({r['fail']})")
            else:
                extra = f" {r['file_mb']:.0f}MB" if "file_mb" in r else ""
                parts.append(f"{ph} {r['secs']:.1f}s peak {r['peak_mb']:.0f}MB{extra}")
        if "roundtrip_ok" in row:
            parts.append("identical" if row["roundtrip_ok"] else "DIFFERENT")
        print(" | ".join(parts), flush=True)
        all_rows.append(row)

    with tempfile.TemporaryDirectory(prefix="opendta-stress-") as td:
        tmp = Path(td)
        wanted = args.series.split(",")
        if "long" in wanted:
            series("long", long_steps, "survey", ["save", "use"], tmp, ceiling, args.timeout, log)
        if "numeric" in wanted:
            series("numeric", long_steps, "numeric", ["save", "use"], tmp, ceiling, args.timeout, log)
        if "wide" in wanted:
            series("wide", wide_steps, "survey", ["save", "use", "readstat"], tmp, ceiling,
                   args.timeout, log)
        if "csv" in wanted:
            series("csv", csv_steps, "survey", ["export", "import"], tmp, ceiling, args.timeout, log)

    if args.output:
        Path(args.output).write_text(json.dumps(all_rows, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
