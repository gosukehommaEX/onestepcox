# Print the large-sample limits of hazard ratio estimators

Prints the design and, for each true hazard ratio, the calendar time of
the analysis in the limit and the large-sample bias of each estimator of
the log hazard ratio (`limit - log(HR)`), as a table with one row per
method and one column per hazard ratio.

## Usage

``` r
# S3 method for class 'LimitHR'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  An object of class `"LimitHR"` returned by
  [`LimitHR()`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md).

- digits:

  Number of significant digits. Default 3.

- ...:

  Not used.

## Value

`x`, invisibly. `as.data.frame(x)` gives the data in long format.

## See also

[`LimitHR`](https://gosukehommaex.github.io/onestepcox/reference/LimitHR.md),
[`plot.LimitHR`](https://gosukehommaex.github.io/onestepcox/reference/plot.LimitHR.md).

## Examples

``` r
lim <- LimitHR(HR = c(0.5, 2, 10), hazard_C = log(2) / 18,
               n_T = 100, n_C = 100, a = 12, E = 160)
print(lim)
#> Large-sample limits of the log hazard ratio estimators
#> Design: hazard_C = 0.0385, n_T : n_C = 100 : 100, a = 12, E = 160, dropout hazard = 0
#> 
#> Analysis time:
#> HR = 0.5   HR = 2  HR = 10 
#>     67.6       37       30 
#> 
#> Bias (limit - log HR):
#>                 HR = 0.5    HR = 2  HR = 10
#> MH_Halley              0         0        0
#> MH_Newton              0         0        0
#> Pike_Halley     3.58e-07 -3.85e-07  -0.0101
#> Pike_Newton     1.71e-05 -1.79e-05  -0.0494
#> MH_Chebyshev           0         0        0
#> Pike_Chebyshev  3.75e-07 -4.03e-07  -0.0125
#> MH             -2.22e-16         0        0
#> Pike              0.0167   -0.0171   -0.654
#> Peto             -0.0208    0.0209    0.167
#> Cox                    0 -1.11e-16 -1.6e-14
print(lim, digits = 6)
#> Large-sample limits of the log hazard ratio estimators
#> Design: hazard_C = 0.0385082, n_T : n_C = 100 : 100, a = 12, E = 160, dropout hazard = 0
#> 
#> Analysis time:
#> HR = 0.5   HR = 2  HR = 10 
#>  67.6232  37.0497  30.0387 
#> 
#> Bias (limit - log HR):
#>                    HR = 0.5       HR = 2      HR = 10
#> MH_Halley                 0            0            0
#> MH_Newton                 0            0            0
#> Pike_Halley     3.57758e-07 -3.85178e-07   -0.0101153
#> Pike_Newton     1.71161e-05 -1.78753e-05   -0.0493771
#> MH_Chebyshev              0            0            0
#> Pike_Chebyshev  3.74611e-07 -4.03077e-07   -0.0125088
#> MH             -2.22045e-16            0            0
#> Pike              0.0166641    -0.017091    -0.654148
#> Peto             -0.0207536    0.0208575     0.166844
#> Cox                       0 -1.11022e-16 -1.59872e-14
```
