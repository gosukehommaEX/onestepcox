# Print operating characteristics of hazard ratio estimators

Prints, for each scenario of
[`CompareOC`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md),
a table with one row per method of the number of replicates without an
estimate, the median and maximum deviation from the reference estimator,
the median ratio of the standard errors, the bias, RMSE, and coverage,
and the computing time per replicate. All columns are in the component
`summary`.

## Usage

``` r
# S3 method for class 'CompareOC'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  An object of class `"CompareOC"` returned by
  [`CompareOC()`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md).

- digits:

  Number of significant digits. Default 3.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`CompareOC`](https://gosukehommaex.github.io/onestepcox/reference/CompareOC.md),
[`plot.CompareOC`](https://gosukehommaex.github.io/onestepcox/reference/plot.CompareOC.md).

## Examples

``` r
scenarios <- data.frame(label = "demo", n_T = 100, n_C = 100,
                        hazard_C = log(2) / 18, HR = 0.7, a = 12,
                        E = 110)
res <- CompareOC(scenarios, c("Cox", "Pike", "MH_Halley"), nsim = 20,
                 seed = 1)
print(res, digits = 2)
#> Operating characteristics of hazard ratio estimators
#> 1 scenario; 20 replicates each; reference: Cox; ties: efron; conf.level: 0.95
#> 
#> Scenario demo (true HR = 0.7)
#>           n_NA dev_med dev_max se_ratio   bias rmse coverage    time
#> Cox          0 0.0e+00 0.0e+00     1.00 -0.088 0.21     0.95 0.00185
#> Pike         0 1.4e-03 1.0e-02     0.99 -0.085 0.20     0.95 0.00025
#> MH_Halley    0 1.4e-09 4.2e-08     1.00 -0.088 0.21     0.95 0.00030
#> 
#> dev: |coef - coef of the reference| (median and maximum); se_ratio: median of
#> se / se of the reference; time: seconds per replicate. All summaries are in x$summary.
```
