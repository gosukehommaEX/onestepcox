#' Large-sample limits of the hazard ratio estimators
#'
#' Computes the limits in probability of the log hazard ratio estimators of
#' \code{\link{EstHR}} under the design of \code{\link{GenSimData}}, as the
#' sample size grows with the ratios \code{n_T : n_C : E} held fixed. The
#' limit of the Mantel-Haenszel and Cox estimators is the true log hazard
#' ratio, whereas the Pike and Peto estimators, and the one-step corrections
#' from the Pike anchor, converge to other values. The difference from the
#' true log hazard ratio is the large-sample bias of each estimator.
#'
#' @param HR Positive numeric vector of true hazard ratios.
#' @param hazard_C Positive numeric. Hazard rate of the control group.
#' @param n_T,n_C Positive numbers of subjects in the treatment and control
#'   groups. Only their ratio matters.
#' @param a Positive numeric. Length of the accrual period.
#' @param E Positive number of events at which the analysis is triggered.
#'   Only the ratio \code{E / (n_T + n_C)} matters, and it must be below the
#'   fraction of subjects who would ever have an event.
#' @param dropout_hazard Non-negative numeric. Hazard rate of independent
#'   dropout. Default 0.
#'
#' @return An object of class \code{c("LimitHR", "data.frame")}, a data
#'   frame with one row per hazard ratio and method and the columns
#'   \code{HR}, \code{analysis_time} (calendar time of the analysis in the
#'   limit), \code{method} (as in \code{EstHR()}), \code{limit} (limit of the
#'   estimator of the log hazard ratio), and \code{bias}
#'   (\code{limit - log(HR)}). The methods are \code{"MH_Halley"},
#'   \code{"MH_Newton"}, \code{"Pike_Halley"}, \code{"Pike_Newton"},
#'   \code{"MH_Chebyshev"}, \code{"Pike_Chebyshev"}, \code{"MH"},
#'   \code{"Pike"}, \code{"Peto"}, and \code{"Cox"}. The attribute
#'   \code{design} holds the arguments other than \code{HR}.
#'   \code{as.data.frame()} gives an object of class \code{"data.frame"}. See
#'   \code{\link{print.LimitHR}} and \code{\link{plot.LimitHR}} for the
#'   methods of \code{print()} and \code{plot()}.
#'
#' @details
#' Let \code{pi_T = n_T / (n_T + n_C)} and \code{pi_C = 1 - pi_T},
#' \code{lambda_C = hazard_C}, \code{lambda_T = HR * hazard_C}, and
#' \code{lambda_D = dropout_hazard}. In the limit, the analysis takes place at
#' the calendar time \code{tau} at which the expected fraction of subjects
#' with an event reaches \code{E / (n_T + n_C)}. Divided by the sample size,
#' the number at risk in group j at time \code{t} since entry tends to
#' \code{y_j(t) = pi_j G(t) exp(-(lambda_j + lambda_D) t)}, where
#' \code{G(t) = min(1, (tau - t) / a)} for \code{0 <= t <= tau} is the
#' fraction of subjects followed for at least \code{t}, and the sums over the
#' event times in \code{EstHR()} tend to integrals over \code{(0, tau)} with
#' respect to \code{dN(t) = (lambda_T y_T(t) + lambda_C y_C(t)) dt}. For
#' example, the score \code{U(eta)} divided by the sample size tends to
#' \code{u(eta) = integral of (lambda_T y_T q - lambda_C y_C p) dt} with
#' \code{p = y_T exp(eta) / (y_C + y_T exp(eta))} and \code{q = 1 - p}, which
#' vanishes at \code{eta = log(HR)}. The limits of the Mantel-Haenszel,
#' Pike, and Peto estimators are obtained in the same way, and the limits of
#' the one-step corrections apply the steps of \code{EstHR()} to the limits
#' of the score and its derivatives at the limit of the anchor. The Cox limit
#' is the root of \code{u(eta)}, found numerically as a check.
#'
#' The integrals are computed by \code{stats::integrate()} with a relative
#' tolerance of \code{1e-10}, split at \code{tau - a} when \code{tau > a}.
#' Times are continuous, so the Efron and Breslow handling of ties give the
#' same limits, and strata are not considered.
#'
#' @seealso \code{\link{GenSimData}}, \code{\link{EstHR}}.
#'
#' @examples
#' lim <- LimitHR(HR = c(0.5, 2, 10), hazard_C = log(2) / 18,
#'                n_T = 100, n_C = 100, a = 12, E = 160)
#' lim
#' plot(lim)
#' # Bias of the Pike anchor and of the Halley steps from each anchor
#' subset(as.data.frame(lim), method == "Pike" | method == "Pike_Halley" |
#'          method == "MH_Halley")
#'
#' @export
LimitHR <- function(HR, hazard_C, n_T, n_C, a, E, dropout_hazard = 0) {
  if (!is.numeric(HR) || length(HR) == 0L || anyNA(HR) || any(HR <= 0)) {
    stop("'HR' must be positive numbers")
  }
  stopifnot(length(hazard_C) == 1L, hazard_C > 0, length(n_T) == 1L,
            n_T > 0, length(n_C) == 1L, n_C > 0, length(a) == 1L, a > 0,
            length(E) == 1L, E > 0, length(dropout_hazard) == 1L,
            dropout_hazard >= 0)
  pi_T <- n_T / (n_T + n_C)
  pi_C <- 1 - pi_T
  e_frac <- E / (n_T + n_C)
  methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
               "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto", "Cox")

  rows <- lapply(HR, function(theta) {
    lam_T <- theta * hazard_C
    lam_C <- hazard_C
    tau <- limit_analysis_time(lam_T, lam_C, pi_T, a, e_frac, dropout_hazard)
    breaks <- if (tau > a) c(0, tau - a, tau) else c(0, tau)

    # Limits of the numbers at risk per subject, and the event intensity
    G <- function(t) pmin(1, (tau - t) / a)
    y_T <- function(t) pi_T * G(t) * exp(-(lam_T + dropout_hazard) * t)
    y_C <- function(t) pi_C * G(t) * exp(-(lam_C + dropout_hazard) * t)
    # Probability that an event at time t is in the treatment group under
    # the model with log hazard ratio eta
    p_eta <- function(t, eta) {
      stats::plogis(log(pi_T / pi_C) - (lam_T - lam_C) * t + eta)
    }
    integ <- function(f) {
      total <- 0
      for (k in seq_len(length(breaks) - 1L)) {
        total <- total + stats::integrate(f, breaks[k], breaks[k + 1L],
                                          rel.tol = 1e-10, abs.tol = 1e-15,
                                          subdivisions = 1000L)$value
      }
      total
    }
    # Score, written as one integrand so that it is accurate near its root,
    # and its first two derivatives (I = -u', J = -u'')
    u_score <- function(eta) {
      integ(function(t) {
        p <- p_eta(t, eta)
        lam_T * y_T(t) * (1 - p) - lam_C * y_C(t) * p
      })
    }
    u_derivs <- function(eta) {
      u <- u_score(eta)
      I <- integ(function(t) {
        p <- p_eta(t, eta)
        p * (1 - p) * (lam_T * y_T(t) + lam_C * y_C(t))
      })
      J <- integ(function(t) {
        p <- p_eta(t, eta)
        p * (1 - p) * (1 - 2 * p) * (lam_T * y_T(t) + lam_C * y_C(t))
      })
      c(U = u, I = I, J = J)
    }

    # Closed-form estimators
    O_T <- integ(function(t) lam_T * y_T(t))
    O_C <- integ(function(t) lam_C * y_C(t))
    E_T <- integ(function(t) {
      p_eta(t, 0) * (lam_T * y_T(t) + lam_C * y_C(t))
    })
    E_C <- integ(function(t) {
      (1 - p_eta(t, 0)) * (lam_T * y_T(t) + lam_C * y_C(t))
    })
    R <- integ(function(t) lam_T * y_T(t) * (1 - p_eta(t, 0)))
    S <- integ(function(t) lam_C * y_C(t) * p_eta(t, 0))
    d_null <- u_derivs(0)
    lim <- c(MH = log(R / S),
             Pike = log((O_T * E_C) / (O_C * E_T)),
             Peto = d_null[["U"]] / d_null[["I"]])

    # One-step corrections from the limits of the anchors
    for (anchor in c("MH", "Pike")) {
      eta0 <- lim[[anchor]]
      d0 <- u_derivs(eta0)
      U <- d0[["U"]]
      I <- d0[["I"]]
      J <- d0[["J"]]
      lim[paste0(anchor, "_Newton")] <- eta0 + U / I
      lim[paste0(anchor, "_Halley")] <- eta0 + 2 * U * I / (2 * I * I + U * J)
      lim[paste0(anchor, "_Chebyshev")] <- eta0 + U / I - J * U * U /
        (2 * I * I * I)
    }

    # Cox: root of the limiting score
    lim["Cox"] <- stats::uniroot(u_score, log(theta) + c(-1, 1),
                                 tol = 1e-12, maxiter = 1000L)$root

    data.frame(HR = theta, analysis_time = tau, method = methods,
               limit = unname(lim[methods]),
               bias = unname(lim[methods]) - log(theta),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  class(out) <- c("LimitHR", "data.frame")
  attr(out, "design") <- c(hazard_C = hazard_C, n_T = n_T, n_C = n_C,
                           a = a, E = E, dropout_hazard = dropout_hazard)
  out
}
