#' @title Calculate AUC for a fitted model.
#' @description A function that finds area under the curve for a fitted `nls` object across a lower and upper bound.
#' @param .fit A fitted `nlsLM` object.
#' @param .lower The lower bound of where to calculate AUC.
#' @param .upper The upper bound of where to calculate AUC.
#' @return A list of class `integrate` as produced by `stats::integrate`.
#' @importFrom stats integrate predict
#' @export

calculate_auc <- function(.fit, .lower, .upper) {
  stats::integrate(
    f = function(x) stats::predict(.fit, data.frame(x)),
    lower = .lower,
    upper = .upper
  )
}
