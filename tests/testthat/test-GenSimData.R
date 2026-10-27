# Tests for GenSimData()

test_that("Output has the documented structure", {
  d <- GenSimData(nsim = 5, n_T = 50, n_C = 50,
                  hazard_C = log(2) / 18, HR = 0.7,
                  a = 12, E = 60, seed = 1)
  expect_s3_class(d, "data.table")
  expect_named(d, c("simID", "group", "time", "event"),
               ignore.order = TRUE)
  expect_s3_class(d$group, "factor")
  expect_setequal(levels(d$group), c("T", "C"))
  expect_true(all(d$event %in% c(0L, 1L)))
  expect_true(all(d$time >= 0))
})

test_that("Same seed yields identical output", {
  args <- list(nsim = 10, n_T = 50, n_C = 50,
               hazard_C = log(2) / 18, HR = 0.8,
               a = 12, E = 70, seed = 42)
  d1 <- do.call(GenSimData, args)
  d2 <- do.call(GenSimData, args)
  expect_equal(d1, d2)
})

test_that("Each simulation contains exactly E events", {
  E <- 60
  d <- GenSimData(nsim = 5, n_T = 50, n_C = 50,
                  hazard_C = log(2) / 18, HR = 0.7,
                  a = 12, E = E, seed = 3)
  events_per_sim <- d[, sum(event), by = simID]$V1
  expect_true(all(events_per_sim == E))
})

test_that("Unequal allocation is respected", {
  d <- GenSimData(nsim = 3, n_T = 30, n_C = 60,
                  hazard_C = log(2) / 18, HR = 0.7,
                  a = 12, E = 50, seed = 4)
  tab <- table(d[simID == 1L, group])
  # All accrued subjects from sim 1 should still be present at most up to
  # the original allocation; the relative ordering n_C > n_T must hold.
  expect_gte(tab[["C"]], tab[["T"]])
})

test_that("E larger than n_T + n_C errors out", {
  expect_error(
    GenSimData(nsim = 2, n_T = 10, n_C = 10,
               hazard_C = 0.1, HR = 0.7, a = 12, E = 100),
    "larger than the total sample size"
  )
})

test_that("Dropout introduces censored subjects beyond the no-dropout case", {
  # Use a small dropout hazard relative to the event hazard so both
  # configurations can hit E events. With dropout, some subjects are
  # censored due to dropout (event = 0), in addition to the
  # administrative censoring induced by the analysis cutoff.
  args_no <- list(nsim = 20, n_T = 50, n_C = 50,
                  hazard_C = log(2) / 18, HR = 0.7,
                  a = 12, E = 40, seed = 5)
  d_no <- do.call(GenSimData, args_no)
  args_yes <- modifyList(args_no, list(dropout_hazard = 0.01))
  d_yes <- do.call(GenSimData, args_yes)
  n_no_event_no  <- d_no[event == 0L, .N]
  n_no_event_yes <- d_yes[event == 0L, .N]
  expect_gte(n_no_event_yes, n_no_event_no)
})

# ---- Strata and rounding of times ---------------------------------------

gen_args <- list(nsim = 20, n_T = 100, n_C = 100, hazard_C = log(2) / 18,
                 HR = 0.7, a = 12, E = 120)

test_that("Default strata and rounding arguments give the earlier output", {
  d1 <- do.call(GenSimData, c(gen_args, list(seed = 1)))
  d2 <- do.call(GenSimData, c(gen_args, list(strata_prop = NULL,
                                             strata_hazard = NULL,
                                             time_unit = 0, seed = 1)))
  expect_identical(d1, d2)
  expect_named(d1, c("simID", "group", "time", "event"))
})

test_that("Strata are assigned with the given probabilities", {
  d <- do.call(GenSimData, c(gen_args, list(strata_prop = c(0.2, 0.3, 0.5),
                                            seed = 6)))
  expect_named(d, c("simID", "stratum", "group", "time", "event"))
  expect_true(all(d$stratum %in% 1:3))
  p <- as.numeric(table(factor(d$stratum, levels = 1:3))) / nrow(d)
  # About 4,000 subjects: the standard error of each proportion is below 0.01
  expect_true(all(abs(p - c(0.2, 0.3, 0.5)) < 0.04))
  expect_true(all(d[, sum(event), by = simID]$V1 == 120))
})

test_that("A larger stratum hazard gives a larger event fraction", {
  d <- do.call(GenSimData, c(gen_args, list(strata_prop = c(0.5, 0.5),
                                            strata_hazard = c(1, 4),
                                            seed = 7)))
  expect_gt(mean(d$event[d$stratum == 2]), mean(d$event[d$stratum == 1]))
})

test_that("Stratified output is sorted for is_sorted = TRUE", {
  d <- do.call(GenSimData, c(gen_args, list(strata_prop = c(0.5, 0.5),
                                            strata_hazard = c(1, 2),
                                            seed = 8)))
  d1 <- d[simID == 3L]
  a <- EstHR(d1$time, d1$event, d1$group, "C", "MH_Halley",
             strata = d1$stratum)
  b <- EstHR(d1$time, d1$event, d1$group, "C", "MH_Halley",
             strata = d1$stratum, is_sorted = TRUE)
  expect_identical(a, b)
})

test_that("time_unit rounds the recorded times up and keeps the events", {
  u <- 1
  d0 <- do.call(GenSimData, c(gen_args, list(seed = 9)))
  d1 <- do.call(GenSimData, c(gen_args, list(time_unit = u, seed = 9)))
  expect_true(all(abs(d1$time / u - round(d1$time / u)) < 1e-9))
  expect_equal(sum(d1$event), sum(d0$event))
  for (s in 1:3) {
    expect_equal(sort(d1$time[d1$simID == s]),
                 sort(ceiling(d0$time[d0$simID == s] / u) * u))
  }
  ev_times <- d1$time[d1$simID == 1L & d1$event == 1L]
  expect_gt(anyDuplicated(ev_times), 0L)
})

test_that("Invalid strata and rounding arguments give errors", {
  expect_error(do.call(GenSimData, c(gen_args, list(strata_prop = c(0.5, 0.6)))),
               "strata_prop")
  expect_error(do.call(GenSimData, c(gen_args, list(strata_prop = c(0.5, 0.5),
                                                    strata_hazard = c(1, 2, 3)))),
               "strata_hazard")
  expect_error(do.call(GenSimData, c(gen_args, list(strata_hazard = c(1, 2)))),
               "requires")
  expect_error(do.call(GenSimData, c(gen_args, list(time_unit = -1))),
               "time_unit")
})
