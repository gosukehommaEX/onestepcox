#' Compare operating characteristics of HR estimators across scenarios
#'
#' Runs Monte Carlo simulations across one or more two-group trial scenarios
#' and compares the operating characteristics of the supplied HR
#' estimation methods. Each scenario specifies its own design parameters
#' (sample sizes per group, hazard rate, true HR, accrual length, target
#' number of events, optional dropout). The function uses
#' \code{\link{GenSimData}} to generate replicate datasets and
#' \code{\link{EstHR}} to compute estimates.
#'
#' @param scenarios A data frame with one row per scenario. Required
#'   columns: \code{label} (character or factor),
#'   \code{n_T} (treatment-group size), \code{n_C} (control-group size),
#'   \code{hazard_C}, \code{HR}, \code{a}, \code{E}.
#'   Optional column: \code{dropout_hazard} (defaults to 0).
#' @param methods Character vector of methods to compare. Each must be a
#'   value accepted by \code{\link{EstHR}}, e.g. one of \code{"Cox"},
#'   \code{"Peto"}, \code{"Pike"}, \code{"MH_Halley"}. The estimates are
#'   computed with the default settings of \code{EstHR()} (Efron handling of
#'   ties, no strata).
#' @param nsim Integer. Number of Monte Carlo replicates per scenario.
#' @param seed Optional integer. Base seed; scenario \code{i} uses
#'   \code{seed + 100 * (i - 1)}.
#' @param reference Character string naming the reference method whose
#'   per-replicate estimates are used as the target for bias and RMSE
#'   calculations. Default \code{"Cox"}. The reference method itself
#'   appears in \code{summary} with \code{bias = 0} and \code{rmse = 0}.
#'
#' @return A list with two components.
#'   \describe{
#'     \item{estimates}{A data frame in long format with columns
#'       \code{label}, \code{simID}, \code{method}, \code{HR}.}
#'     \item{summary}{A data frame with one row per (scenario, method)
#'       combination, containing \code{label}, \code{HR_true},
#'       \code{method}, \code{mean}, \code{bias}, \code{rmse},
#'       \code{time_per_sim} (seconds per simulation, averaged across
#'       replicates).}
#'   }
#'
#' @details
#' The reference method is run first for each scenario; its per-replicate
#' estimates form the target against which other methods' bias and RMSE
#' are computed. This matches the article convention where the Cox MLE
#' plays the role of "ground truth".
#'
#' @examples
#' scenarios <- data.frame(
#'   label    = "demo",
#'   n_T      = 100,
#'   n_C      = 100,
#'   hazard_C = log(2) / 18,
#'   HR       = 0.7,
#'   a        = 12,
#'   E        = 110
#' )
#' res <- CompareOC(
#'   scenarios = scenarios,
#'   methods   = c("Cox", "Peto", "Pike", "MH_Halley"),
#'   nsim      = 50,
#'   seed      = 1
#' )
#' res$summary
#'
#' @import data.table
#' @export
CompareOC <- function(scenarios, methods, nsim,
                      seed = NULL, reference = "Cox") {

  required_cols <- c("label", "n_T", "n_C", "hazard_C", "HR", "a", "E")
  missing_cols  <- setdiff(required_cols, names(scenarios))
  if (length(missing_cols) > 0L) {
    stop("scenarios is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!"dropout_hazard" %in% names(scenarios)) {
    scenarios$dropout_hazard <- 0
  }
  if (!reference %in% methods) {
    stop("'reference' (", reference, ") must be one of 'methods'")
  }

  # Move the reference method to the front so it is computed first
  methods_ord <- c(reference, setdiff(methods, reference))

  base_seed <- if (is.null(seed)) NULL else as.integer(seed)

  est_rows  <- vector("list", nrow(scenarios))
  summ_rows <- vector("list", nrow(scenarios) * length(methods_ord))
  summ_idx  <- 1L

  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, ]
    seed_i <- if (is.null(base_seed)) NULL else base_seed + 100L * (i - 1L)

    # Generate one big replicate set for this scenario
    dat <- GenSimData(
      nsim           = nsim,
      n_T            = as.integer(sc$n_T),
      n_C            = as.integer(sc$n_C),
      hazard_C       = sc$hazard_C,
      HR             = sc$HR,
      a              = sc$a,
      E              = as.integer(sc$E),
      dropout_hazard = sc$dropout_hazard,
      seed           = seed_i
    )
    data.table::setkey(dat, simID, time)

    # Run each method; reference goes first
    HR_by_method <- vector("list", length(methods_ord))
    names(HR_by_method) <- methods_ord
    time_by_method <- numeric(length(methods_ord))
    names(time_by_method) <- methods_ord

    for (m in methods_ord) {
      tm <- system.time({
        ests <- dat[, .(HR = EstHR(time = time,
                                   event = event,
                                   group = group,
                                   control = "C",
                                   method = m,
                                   is_sorted = TRUE)[["exp(coef)"]]),
                    by = simID]
      })[["elapsed"]]
      HR_by_method[[m]]   <- ests$HR
      time_by_method[m]   <- tm / nsim
    }

    sim_ids <- sort(unique(dat$simID))

    # Estimate-level long table
    est_i <- do.call(rbind, lapply(methods_ord, function(m) {
      data.frame(label  = as.character(sc$label),
                 simID  = sim_ids,
                 method = m,
                 HR     = HR_by_method[[m]],
                 stringsAsFactors = FALSE)
    }))
    est_rows[[i]] <- est_i

    # Summary table (one row per (scenario, method))
    ref_HR <- HR_by_method[[reference]]
    for (m in methods_ord) {
      hr_m <- HR_by_method[[m]]
      bias <- mean(hr_m - ref_HR, na.rm = TRUE)
      rmse <- if (m == reference) 0
      else sqrt(mean((hr_m - ref_HR)^2, na.rm = TRUE))
      summ_rows[[summ_idx]] <- data.frame(
        label        = as.character(sc$label),
        HR_true      = sc$HR,
        method       = m,
        mean         = mean(hr_m, na.rm = TRUE),
        bias         = bias,
        rmse         = rmse,
        time_per_sim = unname(time_by_method[m]),
        stringsAsFactors = FALSE
      )
      summ_idx <- summ_idx + 1L
    }
  }

  list(
    estimates = do.call(rbind, est_rows),
    summary   = do.call(rbind, summ_rows)
  )
}
