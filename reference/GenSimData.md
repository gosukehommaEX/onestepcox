# Generate simulated time-to-event data for two-group trials

Generates Monte Carlo replicates of a two-group parallel-group trial
with exponential survival times, uniform accrual, optional independent
dropout, and analysis triggered when the E-th overall event is observed.
Unequal allocation between treatment and control groups is supported via
separate `n_T` and `n_C` arguments. Optionally, subjects are assigned to
strata with their own baseline hazards, and the recorded times are
rounded up to a time unit, which produces tied times.

## Usage

``` r
GenSimData(
  nsim = 1000,
  n_T,
  n_C,
  hazard_C,
  HR,
  a,
  E,
  dropout_hazard = 0,
  strata_prop = NULL,
  strata_hazard = NULL,
  time_unit = 0,
  seed = NULL
)
```

## Arguments

- nsim:

  Integer. Number of Monte Carlo replicates.

- n_T:

  Integer. Number of subjects on the treatment group (group label
  `"T"`).

- n_C:

  Integer. Number of subjects on the control group (group label `"C"`).

- hazard_C:

  Positive numeric. Hazard rate for the control group.

- HR:

  Positive numeric. True hazard ratio, `hazard_T = HR * hazard_C`.

- a:

  Positive numeric. Length of the accrual period; subject accrual times
  are drawn from `Uniform(0, a)`.

- E:

  Positive integer. Required total number of observed events (pooled
  across both groups) at which the analysis is triggered.

- dropout_hazard:

  Non-negative numeric. Hazard rate of independent dropout. Set to 0
  (default) to disable dropout.

- strata_prop:

  Optional numeric vector of stratum probabilities (positive, summing to
  1). Each subject is assigned to stratum `s` with probability
  `strata_prop[s]`, independently of the group. Default `NULL` (no
  strata).

- strata_hazard:

  Optional positive numeric vector of the same length as `strata_prop`:
  the hazard of a subject in stratum `s` is multiplied by
  `strata_hazard[s]` in both groups, so that the hazard ratio is common
  to all strata. Default `NULL`, which means 1 for every stratum.

- time_unit:

  Non-negative numeric. If positive, the recorded times are rounded up
  to multiples of `time_unit`, which produces tied times; whether a
  subject has an event and the analysis time are determined from the
  exact times. Default 0 (no rounding).

- seed:

  Optional integer seed for reproducibility.

## Value

A `data.table` with one row per subject per simulation, carrying the
columns `simID`, `group` (factor with levels `"T"` and `"C"`), `time`,
and `event`, and, when `strata_prop` is given, `stratum` (integer
stratum number). Subjects accrued after the E-th event are dropped, and
subjects whose total calendar time exceeds the analysis time are
administratively censored. The rows are sorted by `simID` and `time`, or
by `simID`, `stratum`, and `time` when there are strata, so that the
rows of each simulation can be passed to
[`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
with `is_sorted = TRUE`.

## Details

The data-generating process is:

- Accrual times `u_{i,j} ~ Uniform(0, a)`;

- Survival times `t_{i,j} ~ Exp(lambda_j)` with
  `lambda_T = HR * lambda_C`;

- Optional independent dropout `d_{i,j} ~ Exp(lambda_D)`;

- Analysis triggered when the E-th overall event is observed.

With strata, `lambda_C` is multiplied by `strata_hazard[s]` for the
subjects of stratum `s`, and `lambda_T` by the same factor. The `time`
and `event` columns (and `stratum`) are the variables to be passed to
[`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md).
Without `strata_prop` and with `time_unit = 0`, the random numbers are
drawn in the same order as in earlier versions, so that a given `seed`
gives the same data.

## Examples

``` r
d <- GenSimData(nsim = 100, n_T = 100, n_C = 100,
                hazard_C = log(2) / 18, HR = 0.7,
                a = 12, E = 110, seed = 1)
head(d)
#> Key: <simID, time>
#>    simID  group        time event
#>    <int> <fctr>       <num> <int>
#> 1:     1      C 0.002155094     1
#> 2:     1      T 0.199681083     1
#> 3:     1      C 0.442931665     1
#> 4:     1      C 1.244341219     1
#> 5:     1      C 1.743443761     1
#> 6:     1      T 1.804099459     1

# Two strata with different baseline hazards, times in whole months
d2 <- GenSimData(nsim = 10, n_T = 100, n_C = 100,
                 hazard_C = log(2) / 18, HR = 0.7, a = 12, E = 110,
                 strata_prop = c(0.4, 0.6), strata_hazard = c(1, 2),
                 time_unit = 1, seed = 1)
head(d2)
#> Key: <simID, stratum, time>
#>    simID stratum  group  time event
#>    <int>   <int> <fctr> <num> <int>
#> 1:     1       1      T     2     1
#> 2:     1       1      C     2     1
#> 3:     1       1      C     2     1
#> 4:     1       1      C     3     1
#> 5:     1       1      C     3     1
#> 6:     1       1      C     3     1
```
