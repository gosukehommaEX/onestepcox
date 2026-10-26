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
* `GenSimData()` and `CompareOC()` are carried over from the package
  accompanying an earlier version of the article; `CompareOC()` uses the
  hazard ratio element of the new `EstHR()` output.
