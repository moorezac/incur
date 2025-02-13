#' A function that finds the value of `x` for a fitted `nlsModel` for a targeted `y` value
#'
#' @param .fit A fitted `nlsLM` object.
#' @param .x_vals A vector of `x` values in which to search across.
#' @param .target The target `y` value to search for.
#' @return A list produced by `uniroot`. See page for precise details.
#'
#' @export

find_root_fit <- function(.fit, .x_vals, .target) {
  uniroot(
    f = function(x) predict(fit, tibble(x)) - target,
    interval = c(min(x_vals), max(x_vals))
  )
}
