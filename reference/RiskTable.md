# Risk-set summaries at the distinct event times

Tabulates, for every distinct event time (within strata when `strata` is
given), the numbers of treated and control subjects at risk just before
that time and the numbers of treated and total events at that time. All
estimators of
[`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
except `"Cox"` are computed from these summaries, and this function
returns them from the same internal code.

## Usage

``` r
RiskTable(
  time,
  event,
  group,
  control,
  strata = NULL,
  is_sorted = FALSE,
  data = NULL
)
```

## Arguments

- time:

  Numeric vector of follow-up times, or a formula
  `Surv(time, event) ~ group`, optionally with `+ strata(s)` (see
  Formula interface in Details).

- event:

  Vector of event indicators, 1 (event) or 0 (censored), or logical,
  aligned with `time`. With a formula, the data frame may be given here
  instead of in `data`, as in `coxph(formula, data)`.

- group:

  Vector of group labels aligned with `time`. Subjects with
  `group != control` form the treatment group. Not used with a formula.

- control:

  Single value of `group` that identifies the control group. With a
  formula, the default is the first level of the group variable if it is
  a factor and its smallest value otherwise, which is the reference
  level of `coxph()`.

- strata:

  Optional vector of stratum labels aligned with `time`. Each stratum
  has its own baseline hazard and the HR is common to all strata.
  Default `NULL` (no stratification). With a formula, strata are given
  by `strata()` in the formula instead.

- is_sorted:

  Logical. If `TRUE`, the rows are taken to be sorted in ascending order
  of `time` within strata, with the rows of each stratum contiguous, and
  the internal sorting is skipped. An error is raised if they are not.
  Default `FALSE`.

- data:

  A data frame (or list) in which the variables of a formula are
  evaluated. Used only with a formula. Default `NULL` (the environment
  of the formula).

## Value

A data frame with one row per distinct event time within stratum,
ordered by stratum and time, and the columns

- stratum:

  Stratum label (only when `strata` is given).

- time:

  Event time.

- n_T, n_C:

  Numbers at risk in the treatment and control groups just before
  `time`.

- O_T:

  Number of events in the treatment group at `time`.

- O:

  Total number of events at `time`.

The counts are stored as doubles.

## Details

Subjects censored at an event time are counted as at risk at that time,
as in
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).
Times are tied only when they are exactly equal. The data are given as
vectors or by a formula `Surv(time, event) ~ group + strata(s)`, as in
[`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md)
(see Formula interface in its Details); rows with missing values are
then removed.

## See also

[`EstHR`](https://gosukehommaex.github.io/onestepcox/reference/EstHR.md).

## Examples

``` r
time  <- c(1, 1, 2, 2, 3, 3, 4, 5, 5, 6)
event <- c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1)
group <- c("C", "T", "C", "T", "C", "T", "T", "C", "C", "T")
RiskTable(time, event, group, control = "C")
#>   time n_T n_C O_T O
#> 1    1   5   5   1 2
#> 2    2   4   4   0 1
#> 3    3   3   3   1 2
#> 4    4   2   2   1 1
#> 5    5   1   2   0 1
#> 6    6   1   0   1 1
RiskTable(time, event, group, control = "C",
          strata = rep(c("a", "b"), each = 5))
#>   stratum time n_T n_C O_T O
#> 1       a    1   2   3   1 2
#> 2       a    2   1   2   0 1
#> 3       a    3   0   1   0 1
#> 4       b    3   3   2   1 1
#> 5       b    4   2   2   1 1
#> 6       b    5   1   2   0 1
#> 7       b    6   1   0   1 1

# Formula interface
head(RiskTable(Surv(rfstime, status) ~ hormon + strata(grade),
               data = survival::gbsg))
#>   stratum time n_T n_C O_T O
#> 1       1  476  28  45   0 1
#> 2       1  559  28  43   1 1
#> 3       1  637  27  43   0 1
#> 4       1  662  27  42   1 1
#> 5       1  707  25  42   0 1
#> 6       1  795  23  39   0 1
```
