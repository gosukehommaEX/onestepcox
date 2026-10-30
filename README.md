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

The vignettes describe the estimators, their agreement with the `survival`
package, simulations, and the large-sample limits. To build them, install with
`remotes::install_github("gosukehommaEX/onestepcox", build_vignettes = TRUE)`
and see `browseVignettes("onestepcox")`.

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
limits, named as in `summary()` of a `coxph` fit; several methods at once give
a matrix with one row per method. `RiskTable()` returns the risk-set summaries
(numbers at risk and events at each distinct event time) that all methods
except `"Cox"` share.

The data are given by a formula `Surv(time, event) ~ group + strata(s)` and a
data frame, as in `survival::coxph()`, or as vectors. With a formula, the
control group is the first level of the group variable (or its smallest
value) unless `control` is given, and rows with missing values are removed.

## Example

```r
library(onestepcox)
est <- EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
             method = c("MH_Halley", "Pike", "Peto", "Cox"))
est        # print(): groups, estimates, z, p-values, and intervals
plot(est)  # forest plot of the hazard ratios

# Breslow ties and strata by tumour grade
EstHR(Surv(rfstime, status) ~ hormon + strata(grade), data = survival::gbsg,
      method = "MH_Halley", ties = "breslow")

# The same estimate from vectors
d <- survival::gbsg
EstHR(d$rfstime, d$status, d$hormon, control = 0, method = "MH_Halley",
      ties = "breslow", strata = d$grade)

# Risk-set summaries
head(RiskTable(Surv(rfstime, status) ~ hormon, data = survival::gbsg))
```

The `plot()` methods of `EstHR()`, `CompareOC()`, and `LimitHR()` results
draw with ggplot2 and return a `ggplot` object, which can be changed further
with ggplot2 (for example `plot(est) + ggplot2::theme(...)`) or saved with
`ggplot2::ggsave()`. They take the labels of the methods (`labels`; by
default the method names with hyphens, such as `MH-Halley`), the colours
(`col`), point shapes and sizes (`shape`, `size`), line types and widths
(`linetype`, `linewidth`), the base font size in points (`base_size`), and
the position of the legend (`legend_pos`: `"bottom"`, `"top"`, `"right"`, or
`"left"` outside the panel, coordinates `c(x, y)` inside it, or `"none"`).
The legend is below the panel by default, except in the forest plot of
`EstHR()`, where the methods label the vertical axis.

## Simulation and large-sample limits

`GenSimData()` generates two-group trials with exponential survival times,
uniform accrual, optional dropout, and an analysis at a fixed number of
events, optionally with strata that have their own baseline hazards and with
times rounded up to a unit (which produces ties). `CompareOC()` runs
`GenSimData()` and `EstHR()` over a set of scenarios and reports, for each
method, the distribution of the deviation from a reference (by default the
Cox estimate), also relative to the standard error of the reference, the
proportion of replicates in which the hazard ratio with two decimals differs
from that of the reference, the ratio of the standard errors, the bias,
standard deviation, RMSE, mean standard error, and coverage relative to the
true hazard ratio, and the computing time.

`LimitHR()` computes the limits of the estimators as the sample size grows
under the design of `GenSimData()`. The Mantel-Haenszel and Cox estimators
converge to the true log hazard ratio, whereas the Pike and Peto estimators,
and the one-step corrections from the Pike anchor, converge to other values.

```r
lim <- LimitHR(HR = c(0.5, 2, 10), hazard_C = log(2) / 18,
               n_T = 100, n_C = 100, a = 12, E = 160)
lim        # bias by method and hazard ratio
plot(lim)  # absolute bias against the hazard ratio
```

## Tests, reproduction, and validation

The tests compare the closed-form estimators, the score derivatives, the
risk-set summaries, and the large-sample limits with values computed
beforehand by an independent Python implementation (`tools/test_expected/`,
not part of the package build), and compare `"Cox"` and `"Cox_NR"` with
`survival::coxph()`.

`inst/reproduce/` reproduces published results: Figure 1 of Berry, Kitchin
and Mock (1991) (Peto and Pike estimates in simulated trials) and the
closed-form results for test data 1 of Therneau and Grambsch (2000).
`inst/validation/` compares all estimators with an independent Python
implementation on 203 datasets and `LimitHR()` with simulations of 100,000
subjects. Each script writes its checks to a `results/` folder, and
`inst/validation/summarize_validation.R` counts them.

## License

GPL (>= 2)
