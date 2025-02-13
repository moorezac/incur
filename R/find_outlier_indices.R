#' An implementation of the ROUT method described in Motulsky and Brown (2006).
#'
#' @param .fit A fitted `nlsLM` object.
#' @param .q The maximum desired FDR.
#' @param .scale_method The method in which to scale the residuals. The original paper uses the `quantile` method. The `MAD` implementation is used within the `dr4pl` package.
#' @return A numeric vector of indices individual data points as outliers. This corresponds to the row number of the data.frame used in the original fit.
#'
#' @export

find_outlier_indices <- function(.fit, .q = 1e-5, .scale_method = "quantile") {
  residuals_vec <- resid(.fit)
  n <- length(residuals_vec)

  # dr4pl implementation
  if (.scale_method == "mad") {
    scale <- mad(residuals_vec)
  }
  # original paper
  if (.scale_method == "quantile") {
    scale <- quantile(
      x = residuals_vec,
      probs = pnorm(1, 0, 1) - pnorm(-1, 0, 1)
    ) * (n / (n - length(coef(.fit))))
  }

  abs_res_sorted <- sort(abs(residuals_vec), index.return = TRUE)$x
  indices_sorted <- sort(abs(residuals_vec), index.return = TRUE)$ix

  alphas <- .q * seq(from = n, to = 1, by = -1) / n
  p_values <- 2 * pt(
    q = abs_res_sorted / scale,
    df = n - length(coef(.fit)),
    lower.tail = FALSE
  )

  indices_fdr <- which(p_values < alphas)

  if (length(indices_fdr) == 0) {
    indices_outlier <- NULL
  } else {
    indices_outlier <-
      indices_sorted[seq(from = min(indices_fdr), to = n, by = 1)]
  }

  indices_outlier
}
