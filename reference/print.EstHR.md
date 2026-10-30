# Print hazard ratio estimates

Prints the estimates of
[`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
in the layout of [`summary()`](https://rdrr.io/r/base/summary.html) of a
`coxph` fit: the numbers of subjects and events by group, a table of
`coef`, `exp(coef)`, `se(coef)`, the Wald statistic
`z = coef / se(coef)` and its two-sided p-value, a table of the hazard
ratio and its confidence limits, and a description of each method.

## Usage

``` r
# S3 method for class 'EstHR'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  An object of class `"EstHR"` returned by
  [`EstHR()`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md).

- digits:

  Number of significant digits. Default
  `max(3, getOption("digits") - 3)`.

- ...:

  Not used.

## Value

`x`, invisibly.

## See also

[`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md),
[`plot.EstHR`](https://gosukehommaex.github.io/onestepcox/reference/plot.EstHR.md).

## Examples

``` r
est <- EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
             method = c("MH_Halley", "Cox"))
print(est)
#> Hazard ratio of the treatment group relative to the control group
#> 
#>   Treatment (1):  n = 246, events = 94
#>   Control (0):    n = 440, events = 205
#>   Ties: efron
#> 
#>              coef exp(coef) se(coef)      z Pr(>|z|)
#> MH_Halley -0.3640    0.6949   0.1250 -2.911   0.0036
#> Cox       -0.3640    0.6949   0.1250 -2.911   0.0036
#> 
#>           exp(coef) lower .95 upper .95
#> MH_Halley    0.6949    0.5438    0.8879
#> Cox          0.6949    0.5438    0.8879
#> 
#> Methods:
#>   MH_Halley  Halley step from the Mantel-Haenszel estimate (proposed)
#>   Cox        Cox estimate by survival::coxph()
print(est, digits = 6)
#> Hazard ratio of the treatment group relative to the control group
#> 
#>   Treatment (1):  n = 246, events = 94
#>   Control (0):    n = 440, events = 205
#>   Ties: efron
#> 
#>                coef exp(coef)  se(coef)        z  Pr(>|z|)
#> MH_Halley -0.364010  0.694884  0.125045 -2.91104 0.0036023
#> Cox       -0.364010  0.694884  0.125045 -2.91104 0.0036023
#> 
#>           exp(coef) lower .95 upper .95
#> MH_Halley  0.694884  0.543844  0.887873
#> Cox        0.694884  0.543844  0.887873
#> 
#> Methods:
#>   MH_Halley  Halley step from the Mantel-Haenszel estimate (proposed)
#>   Cox        Cox estimate by survival::coxph()
```
