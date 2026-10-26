# onestepcox

`onestepcox` computes closed-form estimators of the hazard ratio of a
two-group comparison of survival data that approximate the Cox partial
likelihood estimator without iteration. The proposed estimator starts from the
Mantel-Haenszel estimator and takes one Halley step on the Cox score equation.
Its standard error is updated in the same pass over the risk sets, and a Wald
confidence interval is returned. The Efron and Breslow handling of tied event
times and stratification are supported. The package is written in R only.

## Installation

```r
# install.packages("remotes")
remotes::install_github("gosukehommaEX/onestepcox")
```

## Methods of `EstHR()`

| `method` | Estimator of the log hazard ratio |
|---|---|
| `"MH_Halley"` | Proposed: Mantel-Haenszel anchor and one Halley step |
| `"MH_Newton"`, `"MH_Chebyshev"` | Mantel-Haenszel anchor and one Newton or Chebyshev step |
| `"Pike_Halley"`, `"Pike_Newton"`, `"Pike_Chebyshev"` | Pike anchor and one Halley, Newton, or Chebyshev step |
| `"MH"`, `"Pike"`, `"Peto"` | Closed-form estimators |
| `"Cox"` | Cox estimate by `survival::coxph()` with a tight tolerance |
| `"Cox_NR"` | Cox estimate by Newton-Raphson iteration on the same risk-set summaries as the closed-form estimators |

Every method returns `coef`, `exp(coef)`, `se(coef)`, and the confidence
limits, named as in `summary()` of a `coxph` fit. `RiskTable()` returns the
risk-set summaries (numbers at risk and events at each distinct event time)
that all methods except `"Cox"` share.

## Example

```r
library(onestepcox)
d <- survival::gbsg
EstHR(d$rfstime, d$status, d$hormon, control = 0, method = "MH_Halley")
EstHR(d$rfstime, d$status, d$hormon, control = 0, method = "Cox")

# Breslow ties and strata by tumour grade
EstHR(d$rfstime, d$status, d$hormon, control = 0, method = "MH_Halley",
      ties = "breslow", strata = d$grade)

# Risk-set summaries
head(RiskTable(d$rfstime, d$status, d$hormon, control = 0))
```

## Simulation

`GenSimData()` generates two-group trials with exponential survival times,
uniform accrual, optional dropout, and an analysis at a fixed number of
events. `CompareOC()` runs `GenSimData()` and `EstHR()` over a set of
scenarios and summarizes the estimates.

## Tests

The tests compare the closed-form estimators, the score derivatives, and the
risk-set summaries with values computed beforehand by an independent Python
implementation (`tools/test_expected/`, not part of the package build), and
compare `"Cox"` and `"Cox_NR"` with `survival::coxph()`.

## License

GPL (>= 2)
