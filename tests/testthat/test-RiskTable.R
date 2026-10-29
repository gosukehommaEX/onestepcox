# Tests for RiskTable()

# Ties between groups, and a subject censored at an event time (time 5)
# placed before the event in the input
hand_time  <- c(1, 1, 2, 2, 3, 3, 4, 5, 5, 6)
hand_event <- c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1)
hand_group <- c("C", "T", "C", "T", "C", "T", "T", "C", "C", "T")

test_that("RiskTable() reproduces a hand-computed table", {
  rt <- RiskTable(hand_time, hand_event, hand_group, "C")
  expect_equal(rt, data.frame(time = c(1, 2, 3, 4, 5, 6),
                              n_T  = c(5, 4, 3, 2, 1, 1),
                              n_C  = c(5, 4, 3, 2, 2, 0),
                              O_T  = c(1, 0, 1, 1, 0, 1),
                              O    = c(2, 1, 2, 1, 1, 1)))
})

test_that("RiskTable() reproduces a hand-computed stratified table", {
  rt <- RiskTable(hand_time, hand_event, hand_group, "C",
                  strata = rep(c("a", "b"), each = 5))
  expect_equal(rt, data.frame(stratum = c("a", "a", "a", "b", "b", "b", "b"),
                              time = c(1, 2, 3, 3, 4, 5, 6),
                              n_T  = c(2, 1, 0, 3, 2, 1, 1),
                              n_C  = c(3, 2, 1, 2, 2, 2, 0),
                              O_T  = c(1, 0, 0, 1, 1, 0, 1),
                              O    = c(2, 1, 1, 1, 1, 1, 1)))
})

# Compare the row count and column totals with the Python values
check_totals <- function(name) {
  d <- test_data(name)
  rt <- RiskTable(d$time, d$event, d$group, d$control, strata = d$strata)
  ref <- expected_risk_totals[expected_risk_totals$data == name, ]
  expect_equal(nrow(ref), 1L, info = name)
  expect_equal(nrow(rt), ref$rows, info = name)
  expect_equal(c(sum(rt$n_T), sum(rt$n_C), sum(rt$O_T), sum(rt$O)),
               c(ref$n_T, ref$n_C, ref$O_T, ref$O), info = name)
}

test_that("Totals match the Python values (survival data)", {
  for (nm in c("ovarian", "gbsg", "gbsg_grade")) {
    check_totals(nm)
  }
})

test_that("Totals match the Python values (pharmacoSmoking)", {
  skip_if_not_installed("asaur")
  for (nm in c("pharmacoSmoking", "pharmacoSmoking_ageGroup2")) {
    check_totals(nm)
  }
})

test_that("The survival datasets are those used for the Python values", {
  ov <- survival::ovarian
  expect_equal(c(futime = sum(ov$futime), fustat = sum(ov$fustat),
                 rx = sum(ov$rx)),
               expected_data_sums$ovarian)
  gb <- survival::gbsg
  expect_equal(c(rfstime = sum(gb$rfstime), status = sum(gb$status),
                 hormon = sum(gb$hormon), grade = sum(gb$grade)),
               expected_data_sums$gbsg)
})

test_that("The asaur dataset is that used for the Python values", {
  skip_if_not_installed("asaur")
  ps <- asaur::pharmacoSmoking
  expect_equal(c(ttr = sum(ps$ttr), relapse = sum(ps$relapse),
                 grp_patchOnly = sum(ps$grp == "patchOnly"),
                 ageGroup2_50plus = sum(ps$ageGroup2 == "50+")),
               expected_data_sums$pharmacoSmoking)
})

test_that("Row order does not matter and is_sorted = TRUE gives the same table", {
  d <- test_data("gbsg_grade")
  rt <- RiskTable(d$time, d$event, d$group, d$control, strata = d$strata)
  set.seed(5)
  p <- sample(length(d$time))
  expect_equal(RiskTable(d$time[p], d$event[p], d$group[p], d$control,
                         strata = d$strata[p]), rt)
  o <- order(d$strata, d$time)
  expect_equal(RiskTable(d$time[o], d$event[o], d$group[o], d$control,
                         strata = d$strata[o], is_sorted = TRUE), rt)
})

test_that("The formula interface gives the table of the vector interface", {
  gb <- survival::gbsg
  a <- RiskTable(gb$rfstime, gb$status, gb$hormon, 0, strata = gb$grade)
  b <- RiskTable(Surv(rfstime, status) ~ hormon + strata(grade), data = gb)
  expect_identical(b, a)
  a <- RiskTable(gb$rfstime, gb$status, gb$hormon, 0)
  b <- RiskTable(Surv(rfstime, status) ~ hormon, gb)
  expect_identical(b, a)
  expect_error(RiskTable(hand_time, hand_event, hand_group, "C", data = gb),
               "only with a formula")
})

test_that("Empty input gives a table without rows", {
  rt <- RiskTable(numeric(0), numeric(0), character(0), "C")
  expect_equal(nrow(rt), 0L)
  expect_named(rt, c("time", "n_T", "n_C", "O_T", "O"))
})

test_that("Invalid input gives informative errors", {
  expect_error(RiskTable(hand_time, hand_event[-1], hand_group, "C"),
               "same length")
  expect_error(RiskTable(hand_time, hand_event, hand_group, "X"), "control")
  expect_error(RiskTable(hand_time, hand_event, hand_group, "C",
                         strata = 1:3), "strata")
  expect_error(RiskTable(rev(hand_time), hand_event, hand_group, "C",
                         is_sorted = TRUE), "is_sorted")
})
