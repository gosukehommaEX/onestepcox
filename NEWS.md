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
* `RiskTable()` returns the risk-set summaries shared by the estimators.
* `GenSimData()`, carried over from the package accompanying an earlier
  version of the article, generates stratified data (`strata_prop`,
  `strata_hazard`) and tied times (`time_unit`). Without these arguments it
  gives the same data as before for a given seed.
* `CompareOC()` reports, for each method, the median, 90th and 99th
  percentiles, and maximum of the deviation from a reference estimator, the
  ratio of the standard errors, the bias, RMSE, and coverage relative to the
  true hazard ratio, and the computing time. Scenarios may have strata and
  rounded times.
* `LimitHR()` computes the large-sample limits of the estimators under the
  design of `GenSimData()`.
* `inst/reproduce/` reproduces Figure 1 of Berry, Kitchin and Mock (1991) and
  test data 1 of Therneau and Grambsch (2000); `inst/validation/` compares the
  estimators with an independent Python implementation and `LimitHR()` with
  large-sample simulations.
