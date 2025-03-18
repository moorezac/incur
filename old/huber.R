#' @title IRWLS via Huber loss function.
#' @description Perform iterative-reweighted least squares regression via Huber loss function.
#' @param .fit A fitted `nlsModel` object
#' @param .iter_max The number of iterations to perform
#' @param .k The tuning constant to be used.
#' @param .tol Tolerance level required to be achieved before stopping.
#' @return A fitted `nlsModel` object with updated coefficents.
#' @importFrom dplyr case_when
#' 
huber <- function(.fit, .arguments, .iter_max = 100, .k = 1.345, .tol = 1e-6) {
  iter <- 0
  converged <- FALSE
  fit_old <- .fit

  # loop
  while (iter < .iter_max && !converged) {
    coef_old <- coef(fit_old)
    resids <- residuals(fit_old)
    # mad
    s <- median(abs(resids - median(resids))) * 1.4826
    u <- resids / (s + 1e-10)
    weights_vec <- dplyr::case_when(
      abs(u) <= .k ~ 1,
      .default = .k / abs(u)
    )

    arguments_new <- append(
      .arguments,
      list(weights = weights_vec)
    )

    fit_new <- try_fit(.data, .curve_func, arguments_new)
    coef_new <- coef(fit_new)

    # check convergence
    coef_change <- max(abs(coef_new - coef_old) / (abs(coef_old) + 1e-6))
    converged <- coef_change < .tol

    # update for next iteration
    fit_old <- fit_new
    iter <- iter + 1
  }
  # message(str_glue("achieved {tol} tolerance in {iter} iterations"))
  fit <- fit_new

  return(fit)
}
