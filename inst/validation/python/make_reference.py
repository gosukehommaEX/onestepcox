"""Datasets and reference values for inst/validation/validate_python.R.

Generates 203 datasets with numpy (seed 20261028): 200 random two-group
trials over sample sizes from 30 to 2,000, hazard ratios from 0.1 to 10,
allocation ratios 1:1, 2:1, and 1:2, with and without dropout, continuous
times or times rounded up or down to a unit (ties and zero times), and 1, 3,
or 10 strata with their own baseline hazards; and 3 small datasets without an
estimate (no event in the control group; all events in the treatment group
occur when no control is at risk; no event at all). Times are written with
six significant digits, and both implementations read the written values.
reference.py computes all estimators for both handlings of ties.

Usage (from the root of the package):
    python3 inst/validation/python/make_reference.py
Writes inst/validation/data/python_datasets.csv (one row per subject),
python_datasets_info.csv (one row per dataset), and python_results.csv.
"""
import csv
import os
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import reference  # noqa: E402

out_dir = os.path.join(here, "..", "data")
rng = np.random.default_rng(20261028)
lam_C = np.log(2) / 18


def trial(n, theta, ratio, dropout, unit, rounding, n_strata, e_frac):
    n_T = int(round(n * ratio / (1 + ratio)))
    g = np.r_[np.ones(n_T, bool), np.zeros(n - n_T, bool)]
    stratum = rng.integers(1, n_strata + 1, n)
    mult = np.exp(rng.normal(0, 0.7, n_strata)) if n_strata > 1 else np.ones(1)
    haz = lam_C * mult[stratum - 1] * np.where(g, theta, 1.0)
    acc = rng.uniform(0, 12, n)
    t = rng.exponential(1 / haz)
    d = rng.exponential(1 / dropout, n) if dropout > 0 else np.full(n, np.inf)
    tte = np.minimum(t, d)
    ev_cal = np.sort((acc + t)[t <= d])
    E = min(int(round(e_frac * n)), len(ev_cal))
    cut = ev_cal[E - 1]
    keep = acc <= cut
    time = (np.minimum(acc + tte, cut) - acc)[keep]
    event = ((t <= d) & (acc + t <= cut))[keep].astype(int)
    if unit > 0:
        time = (np.ceil(time / unit) if rounding == "up" else np.floor(time / unit)) * unit
    return time, event, g[keep], stratum[keep]


datasets = []
sizes = [30] * 60 + [60] * 50 + [150] * 40 + [400] * 30 + [1000] * 15 + [2000] * 5
for i, n in enumerate(sizes, start=1):
    theta = float(rng.choice([0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 3.0, 10.0]))
    ratio = float(rng.choice([1.0, 2.0, 0.5]))
    dropout = float(rng.choice([0.0, 0.02]))
    unit = float(rng.choice([0.0, 1.0, 3.0]))
    rounding = str(rng.choice(["up", "down"])) if unit > 0 else "none"
    n_strata = int(rng.choice([1, 1, 3, 10]))
    e_frac = float(rng.choice([0.5, 0.8]))
    time, event, g, s = trial(n, theta, ratio, dropout, unit, rounding, n_strata, e_frac)
    info = dict(dataset=i, n=len(time), theta=theta, ratio=ratio, dropout=dropout, unit=unit,
                rounding=rounding, n_strata=n_strata, e_frac=e_frac, kind="random")
    datasets.append((info, time, event, g, s))

# Datasets without an estimate
special = [
    ("no control event", [1, 2, 3, 4, 5, 6], [1, 0, 1, 0, 1, 0], [1, 0, 1, 0, 1, 0]),
    ("treated events with no control at risk", [1, 2, 5, 6], [1, 1, 1, 1], [0, 0, 1, 1]),
    ("no event", [1, 2, 3, 4], [0, 0, 0, 0], [1, 0, 1, 0]),
]
for j, (kind, t, e, g) in enumerate(special, start=len(datasets) + 1):
    info = dict(dataset=j, n=len(t), theta=float("nan"), ratio=float("nan"), dropout=0.0, unit=0.0,
                rounding="none", n_strata=1, e_frac=float("nan"), kind=kind)
    datasets.append((info, np.array(t, float), np.array(e), np.array(g, bool), np.ones(len(t), int)))

fmt = lambda x: "NA" if not np.isfinite(x) else repr(float(x))
with open(os.path.join(out_dir, "python_datasets.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["dataset", "time", "event", "group", "stratum"])
    for info, time, event, g, s in datasets:
        for k in range(len(time)):
            w.writerow([info["dataset"], f"{time[k]:.6g}", int(event[k]), "T" if g[k] else "C", int(s[k])])
with open(os.path.join(out_dir, "python_datasets_info.csv"), "w", newline="") as f:
    keys = ["dataset", "n", "theta", "ratio", "dropout", "unit", "rounding", "n_strata", "e_frac", "kind"]
    w = csv.writer(f, lineterminator="\n")
    w.writerow(keys)
    for info, *_ in datasets:
        w.writerow([info[k] if not isinstance(info[k], float) else fmt(info[k]) for k in keys])

# Reference values from the written (rounded) data
rows = list(csv.DictReader(open(os.path.join(out_dir, "python_datasets.csv"))))
by_id = {}
for r in rows:
    by_id.setdefault(int(r["dataset"]), []).append(r)
with open(os.path.join(out_dir, "python_results.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["dataset", "ties", "method", "coef", "se"])
    for i in sorted(by_id):
        rr = by_id[i]
        time = np.array([float(r["time"]) for r in rr])
        event = np.array([int(r["event"]) for r in rr])
        trt = np.array([r["group"] == "T" for r in rr])
        stratum = np.array([int(r["stratum"]) for r in rr])
        for ties in ("efron", "breslow"):
            est = reference.estimates(time, event, trt, stratum, ties)
            for m in reference.METHODS:
                w.writerow([i, ties, m, fmt(est[m][0]), fmt(est[m][1])])
print(len(datasets), "datasets,", len(rows), "subjects")
