#' Risk-set summaries at the distinct event times
#'
#' Internal workhorse of \code{RiskTable()} and \code{EstHR()}. For every
#' stratum and every distinct event time in it, counts the treated and control
#' subjects at risk just before that time and the treated and total events at
#' that time. Subjects censored at an event time are counted as at risk at
#' that time, as in \code{survival::coxph()}. Times are tied only when they
#' are exactly equal.
#'
#' @param time Numeric vector of follow-up times.
#' @param event Vector of event indicators coded 1 (event) and 0 (censored),
#'   or logical.
#' @param trt Logical vector, \code{TRUE} for the treatment group.
#' @param strata Optional vector of stratum labels, or \code{NULL}.
#' @param is_sorted Logical. If \code{TRUE}, the rows are taken to be sorted in
#'   ascending order of \code{time} within strata, with the rows of each
#'   stratum contiguous.
#'
#' @return A list with elements \code{stratum} (\code{NULL} when
#'   \code{strata} is \code{NULL}), \code{time}, \code{n_T}, \code{n_C},
#'   \code{O_T}, and \code{O}, one entry per distinct event time within
#'   stratum, ordered by stratum and time. The counts are stored as doubles
#'   so that products of counts cannot overflow.
#'
#' @noRd
risk_sets <- function(time, event, trt, strata = NULL, is_sorted = FALSE) {
  n <- length(time)
  if (n == 0L) {
    return(list(stratum = if (is.null(strata)) NULL else strata[0L],
                time = numeric(0), n_T = numeric(0), n_C = numeric(0),
                O_T = numeric(0), O = numeric(0)))
  }

  # Sort by stratum and time
  if (!is_sorted) {
    ord <- if (is.null(strata)) order(time) else order(strata, time)
    time <- time[ord]
    event <- event[ord]
    trt <- trt[ord]
    if (!is.null(strata)) strata <- strata[ord]
  }

  # First row of each stratum and of each distinct time within stratum
  new_time <- c(TRUE, time[-1L] != time[-n])
  if (is.null(strata)) {
    str_start <- 1L
    if (is_sorted && is.unsorted(time)) {
      stop("'is_sorted = TRUE' but 'time' is not in ascending order")
    }
  } else {
    new_str <- c(TRUE, strata[-1L] != strata[-n])
    str_start <- which(new_str)
    if (is_sorted) {
      if (anyDuplicated(strata[str_start]) > 0L) {
        stop("'is_sorted = TRUE' but the rows of a stratum are not contiguous")
      }
      if (any(time[-1L] < time[-n] & !new_str[-1L])) {
        stop("'is_sorted = TRUE' but 'time' is not in ascending order ",
             "within strata")
      }
    }
    new_time <- new_time | new_str
  }
  str_end <- c(str_start[-1L] - 1L, n)
  str_id <- rep.int(seq_along(str_start), diff(c(str_start, n + 1L)))

  first <- which(new_time)
  blk <- cumsum(new_time)
  n_blk <- length(first)

  # Events per block (a block is one distinct time within a stratum)
  ev <- event == 1
  O <- tabulate(blk[ev], nbins = n_blk)
  O_T <- tabulate(blk[ev & trt], nbins = n_blk)

  # At risk at the first row of a block: the rows from there to the end of
  # its stratum. cum_T[i] is the number of treated rows from row i to row n.
  end_blk <- str_end[str_id[first]]
  n_at <- end_blk - first + 1L
  cum_T <- c(rev(cumsum(rev(as.integer(trt)))), 0L)
  n_T <- cum_T[first] - cum_T[end_blk + 1L]
  n_C <- n_at - n_T

  keep <- O > 0L
  pos <- first[keep]
  list(stratum = if (is.null(strata)) NULL else strata[pos],
       time = time[pos],
       n_T = as.numeric(n_T[keep]),
       n_C = as.numeric(n_C[keep]),
       O_T = as.numeric(O_T[keep]),
       O = as.numeric(O[keep]))
}
