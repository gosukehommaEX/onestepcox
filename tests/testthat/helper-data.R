# Datasets of the tests, in the form of the arguments of EstHR(). The
# expected values of helper-expected.R refer to these names.
test_data <- function(name) {
  if (startsWith(name, "gbsg")) {
    d <- survival::gbsg
    out <- list(time = d$rfstime, event = d$status, group = d$hormon,
                control = 0,
                strata = if (name == "gbsg_grade") d$grade else NULL)
  } else if (startsWith(name, "pharmacoSmoking")) {
    d <- asaur::pharmacoSmoking
    out <- list(time = d$ttr, event = d$relapse, group = d$grp,
                control = "combination",
                strata = if (name == "pharmacoSmoking_ageGroup2") {
                  d$ageGroup2
                } else {
                  NULL
                })
  } else {
    d <- survival::ovarian
    out <- list(time = d$futime, event = d$fustat, group = d$rx,
                control = 1, strata = NULL)
  }
  out
}

# The eleven methods of EstHR()
all_methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
                 "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto",
                 "Cox", "Cox_NR")
