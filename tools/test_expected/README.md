# Expected values of the tests

The expected values of the closed-form estimators, the Newton-Raphson Cox
estimate, the score derivatives, the risk-set totals, and the large-sample
limits in
`tests/testthat/helper-expected.R` are computed by an implementation in Python
that does not share code with the package. This folder is excluded from the
package build (`.Rbuildignore`).

| File | Content |
|---|---|
| `ref.py` | Independent implementation: risk sets counted by brute force (one comparison per subject and event time), Breslow and Efron derivatives written per event time, Cox estimate by Brent's method on the score |
| `limit_ref.py` | Independent implementation of the large-sample limits of `LimitHR()` with scipy: analysis time from the event fraction integrated numerically over the accrual time, each limiting sum a separate integral, Cox limit by Brent's method; the designs of the tests are defined here (`DESIGNS`) |
| `rda.py` | Minimal reader of R's `.rda` format, used once to export the datasets |
| `export_data.py` | Writes `data/*.csv` from `ovarian` and `gbsg` of 'survival' 3.8-12 and `pharmacoSmoking` of 'asaur' 0.50 (source commits in the script) |
| `make_expected.py` | Reads `data/*.csv`, runs `ref.py` and `limit_ref.py`, and writes `tests/testthat/helper-expected.R` and `expected_values.txt` |
| `data/` | The columns of the three datasets used by the tests |

Regenerate the helper file (numpy and scipy are required):

```sh
python3 tools/test_expected/make_expected.py .
```

The agreement of the Python Cox estimates with `survival::coxph()` is checked
in `test-EstHR.R` ("Cox and Cox_NR agree with the Python root of the score"),
and `test-RiskTable.R` checks that the datasets installed with 'survival' and
'asaur' have the column sums recorded in the helper file.
