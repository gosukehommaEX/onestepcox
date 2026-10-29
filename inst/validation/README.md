# Validation

The scripts compare `onestepcox` with independent routes to the same
quantities and write a table with a status per check to `results/`: PASS,
EXPLAINED (a failure whose reason is registered in `explained.csv`), FAIL,
or INFO (reported without a judgement). Each script first checks that its
judgement fails for a deliberately wrong value. Run the scripts from the
root of the package after `devtools::load_all()` (or with the installed
package), in this order:

| Script | Content | Time |
|---|---|---|
| `validate_python.R` | All estimators of `EstHR()` (coef and se, Efron and Breslow) against an independent Python implementation on 203 datasets: random trials with 30 to 2,000 subjects, hazard ratios 0.1 to 10, allocation 1:1, 2:1, and 1:2, dropout, ties and zero times from rounded times, and 1, 3, or 10 strata, and 3 datasets without an estimate | about a minute |
| `validate_LimitHR.R` | `LimitHR()` against large-sample simulations: the mean difference between each estimator and the Cox estimate at n = 100,000 against the limiting bias, for hazard ratios 0.2, 0.5, 2, 5, and 10 (n = 1,000 and 10,000 as INFO) | a few minutes |
| `summarize_validation.R` | Counts of each status for every results table of `inst/reproduce/` and `inst/validation/`, written to `results/summary.csv` | seconds |

## Python implementation

`python/reference.py` computes the estimators from the definitions in the
documentation of `EstHR()`: the risk sets are counted by comparing every
subject with every event time, and the Cox estimate is the root of the score
found by Brent's method. `python/make_reference.py` generates the datasets
and the reference values in `data/` (numpy and scipy are required; run
`python3 inst/validation/python/make_reference.py` from the root of the
package). The times are written with six significant digits, and both
implementations use the written values.

## Judgement

- `validate_python.R`: PASS if both values are missing, or if
  `|R - Python| <= tol * max(1, |Python|)` with `tol = 1e-10` (`1e-9` for
  `"Cox_NR"` and `1e-8` for `"Cox"`, which uses `survival::coxph()`). The
  comparison of `"Cox"` is INFO when the Cox estimate does not exist.
- `validate_LimitHR.R`: PASS if `|mean(d) - bias| <= 4 SE + 0.01 |bias|`,
  where `d` is the difference between the estimator and the Cox estimate in
  each of 50 trials of 100,000 subjects and `SE` is the standard error of its
  mean. This is checked for MH, Pike, Peto, and the one-step corrections from
  the Pike anchor (30 checks). The one-step corrections from the MH anchor
  are INFO: their limit equals that of MH, and their deviation from the Cox
  estimate at a finite n is a higher-order term of the anchor error (for
  Newton of fixed sign and of order 1/n), which a test of a zero mean at a
  fixed n would reject although the limit is 0.
