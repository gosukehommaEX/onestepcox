"""Export the columns of the R datasets used in the tests to CSV.

Sources (read with rda.py, a minimal reader of R's .rda format):
  survival 3.8-12, data/cancer.rda (objects ovarian and gbsg)
    https://github.com/cran/survival, commit 40b954baeb1b8ca9c67faaca9f9fdbcd242ee903
  asaur 0.50, data/pharmacoSmoking.rda
    https://github.com/cran/asaur, commit 0d782825331f15a729be9f999faec7e13180bc06

Usage: python3 export_data.py <survival source dir> <asaur source dir>
"""
import csv
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from rda import read_rda, as_columns  # noqa: E402

surv_dir, asaur_dir = sys.argv[1], sys.argv[2]
cancer = read_rda(os.path.join(surv_dir, "data", "cancer.rda"))
smoke = read_rda(os.path.join(asaur_dir, "data", "pharmacoSmoking.rda"))

sets = {
    "ovarian": (as_columns(cancer["ovarian"]), ["futime", "fustat", "rx"]),
    "gbsg": (as_columns(cancer["gbsg"]), ["rfstime", "status", "hormon", "grade"]),
    "pharmacoSmoking": (as_columns(smoke["pharmacoSmoking"]), ["ttr", "relapse", "grp", "ageGroup2"]),
}
for name, (cols, keep) in sets.items():
    path = os.path.join(here, "data", name + ".csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(keep)
        for i in range(len(cols[keep[0]])):
            row = []
            for k in keep:
                v = cols[k][i]
                row.append(v if isinstance(v, str) else repr(float(v)))
            w.writerow(row)
    print(name, len(cols[keep[0]]), "rows ->", os.path.relpath(path, here))
