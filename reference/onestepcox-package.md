# onestepcox: One-Step Closed-Form Hazard Ratio Estimators

Closed-form estimators of the hazard ratio of a two-group comparison
that approximate the Cox partial likelihood estimator without iteration.
The proposed estimator takes one Halley step on the Cox score equation
from the Mantel-Haenszel estimator, and its standard error is updated in
the same pass over the risk sets. The Efron and Breslow handling of tied
event times and stratification are supported.

## Details

The main functions are

- [`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md):

  Hazard ratio, standard error, and confidence interval by the proposed
  estimator, other one-step corrections, the Mantel-Haenszel, Pike, and
  Peto estimators, and the Cox estimator.

- [`RiskTable`](https://gosukehommaex.github.io/onestepcox/reference/RiskTable.md):

  Risk-set summaries shared by the estimators.

- [`GenSimData`](https://gosukehommaex.github.io/onestepcox/reference/GenSimData.md):

  Simulated two-group trial data.

- [`CompareOC`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md):

  Operating characteristics of the estimators across simulation
  scenarios.

- [`LimitHR`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md):

  Large-sample limits of the estimators under the design of
  [`GenSimData()`](https://gosukehommaex.github.io/onestepcox/reference/GenSimData.md).

[`EstHR()`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
and
[`RiskTable()`](https://gosukehommaex.github.io/onestepcox/reference/RiskTable.md)
accept vectors or a formula `Surv(time, event) ~ group + strata(s)`, as
in [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).
The results of
[`EstHR()`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md),
[`CompareOC()`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md),
and
[`LimitHR()`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md)
have methods of [`print()`](https://rdrr.io/r/base/print.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html); the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods draw
with `ggplot2` and return `ggplot` objects.

## See also

Useful links:

- <https://github.com/gosukehommaEX/onestepcox>

- <https://gosukehommaex.github.io/onestepcox/>

- Report bugs at <https://github.com/gosukehommaEX/onestepcox/issues>

## Author

**Maintainer**: Gosuke Homma <my.name.is.gosuke@gmail.com>
([ORCID](https://orcid.org/0000-0002-6854-5627))

Authors:

- Gosuke Homma <my.name.is.gosuke@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-6854-5627))
