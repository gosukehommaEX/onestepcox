# Changelog

## onestepcox 0.1.0

- First version, accompanying the article on the one-step closed-form
  hazard ratio estimator anchored at the Mantel-Haenszel estimator.
- [`EstHR()`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
  returns the log hazard ratio, the hazard ratio, the standard error,
  and a Wald confidence interval by eleven methods: the proposed
  `"MH_Halley"`, the one-step corrections `"MH_Newton"`,
  `"Pike_Halley"`, `"Pike_Newton"`, `"MH_Chebyshev"`, and
  `"Pike_Chebyshev"`, the closed-form `"MH"`, `"Pike"`, and `"Peto"`,
  and the Cox estimate by
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  (`"Cox"`) or by Newton-Raphson iteration on the same risk-set
  summaries (`"Cox_NR"`). Tied event times are handled by the Efron
  (default) or the Breslow method, and stratification is supported.
- [`EstHR()`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
  and
  [`RiskTable()`](https://gosukehommaex.github.io/onestepcox/reference/RiskTable.md)
  accept a formula `Surv(time, event) ~ group + strata(s)` with a data
  frame, as in
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html), as
  well as vectors.
  [`EstHR()`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
  takes several methods at once and returns an object of class `"EstHR"`
  with [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.
- [`RiskTable()`](https://gosukehommaex.github.io/onestepcox/reference/RiskTable.md)
  returns the risk-set summaries shared by the estimators.
- [`GenSimData()`](https://gosukehommaex.github.io/onestepcox/reference/GenSimData.md),
  carried over from the package accompanying an earlier version of the
  article, generates stratified data (`strata_prop`, `strata_hazard`)
  and tied times (`time_unit`). Without these arguments it gives the
  same data as before for a given seed.
- [`CompareOC()`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md)
  reports, for each method, the median, 90th and 99th percentiles, and
  maximum of the deviation from a reference estimator, the deviation
  divided by the standard error of the reference, the proportion of
  replicates whose hazard ratio with two decimals differs from that of
  the reference, the ratio of the standard errors, the bias, standard
  deviation, RMSE, mean standard error, and coverage relative to the
  true hazard ratio, and the computing time. Scenarios may have strata
  and rounded times.
- [`LimitHR()`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md)
  computes the large-sample limits of the estimators under the design of
  [`GenSimData()`](https://gosukehommaex.github.io/onestepcox/reference/GenSimData.md).
- The results of
  [`CompareOC()`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md)
  and
  [`LimitHR()`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md)
  have classes `"CompareOC"` and `"LimitHR"` with
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. The
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods of
  the three classes draw with ‘ggplot2’ and return `ggplot` objects.
  They take labels of the methods (by default with hyphens, such as
  `MH-Halley`), colors, point shapes and sizes, line types and widths,
  the base font size, and the position of the legend (outside the panel
  at the bottom, top, right, or left, or inside it). The legend is below
  the panel by default for the results of
  [`CompareOC()`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md)
  and
  [`LimitHR()`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md).
- `inst/reproduce/` reproduces Figure 1 of Berry, Kitchin and
  Mock (1991) and test data 1 of Therneau and Grambsch (2000);
  `inst/validation/` compares the estimators with an independent Python
  implementation and
  [`LimitHR()`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md)
  with large-sample simulations.
- Vignettes (`browseVignettes("onestepcox")`): an introduction, the
  agreement with
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) and
  [`survival::survdiff()`](https://rdrr.io/pkg/survival/man/survdiff.html),
  simulations with
  [`GenSimData()`](https://gosukehommaex.github.io/onestepcox/reference/GenSimData.md)
  and
  [`CompareOC()`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md),
  and the large-sample limits of
  [`LimitHR()`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md).
