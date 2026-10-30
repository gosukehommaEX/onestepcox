# onestepcox 0.1.0

* First version, accompanying the article on the one-step closed-form hazard
  ratio estimator anchored at the Mantel-Haenszel estimator.
* `EstHR()` returns the log hazard ratio, the hazard ratio, the standard
  error, and a Wald confidence interval by eleven methods: the proposed
  `"MH_Halley"`, the one-step corrections `"MH_Newton"`, `"Pike_Halley"`,
  `"Pike_Newton"`, `"MH_Chebyshev"`, and `"Pike_Chebyshev"`, the closed-form
  `"MH"`, `"Pike"`, and `"Peto"`, and the Cox estimate by
  `survival::coxph()` (`"Cox"`) or by Newton-Raphson iteration on the same
  risk-set summaries (`"Cox_NR"`). Tied event times are handled by the Efron
  (default) or the Breslow method, and stratification is supported.
* `EstHR()` and `RiskTable()` accept a formula
  `Surv(time, event) ~ group + strata(s)` with a data frame, as in
  `survival::coxph()`, as well as vectors. `EstHR()` takes several methods at
  once and returns an object of class `"EstHR"` with `print()` and `plot()`
  methods.
* `RiskTable()` returns the risk-set summaries shared by the estimators.
* `GenSimData()`, carried over from the package accompanying an earlier
  version of the article, generates stratified data (`strata_prop`,
  `strata_hazard`) and tied times (`time_unit`). Without these arguments it
  gives the same data as before for a given seed.
* `CompareOC()` reports, for each method, the median, 90th and 99th
  percentiles, and maximum of the deviation from a reference estimator, the
  deviation divided by the standard error of the reference, the proportion
  of replicates whose hazard ratio with two decimals differs from that of
  the reference, the ratio of the standard errors, the bias, standard
  deviation, RMSE, mean standard error, and coverage relative to the true
  hazard ratio, and the computing time. Scenarios may have strata and
  rounded times.
* `LimitHR()` computes the large-sample limits of the estimators under the
  design of `GenSimData()`.
* The results of `CompareOC()` and `LimitHR()` have classes `"CompareOC"`
  and `"LimitHR"` with `print()` and `plot()` methods. The `plot()` methods
  of the three classes draw with 'ggplot2' and return `ggplot` objects. They
  take labels of the methods (by default with hyphens, such as
  `MH-Halley`), colours, point shapes and sizes, line types and widths, the
  base font size, and the position of the legend (outside the panel at the
  bottom, top, right, or left, or inside it). The legend is below the panel
  by default for the results of `CompareOC()` and `LimitHR()`.
* `inst/reproduce/` reproduces Figure 1 of Berry, Kitchin and Mock (1991) and
  test data 1 of Therneau and Grambsch (2000); `inst/validation/` compares the
  estimators with an independent Python implementation and `LimitHR()` with
  large-sample simulations.
* A vignette (`vignette("onestepcox")`) shows the estimators, the simulation
  functions, and the large-sample limits.
