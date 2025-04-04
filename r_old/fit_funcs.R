#' @title Attempt to fit a curve via `minpack.lm`
#' @description Try to fit a curve within a new environment. See help page on `rlang::exec` for background on this approach.
#' @param .data A A `data.frame` or `data.frame` extension (tibble) in long format with columns `x` and `y`. 
#' @param .curve_func A function that describes a curve/model in terms of `x`.
#' @param .arguments A named list of arguments to be injected into the `minpack.lm::nlsLM` function.
#' @return Either a fitted `nls` object if expression is evaluated without error, or an invisible object with class `try-error`.
#' @importFrom rlang env expr
#' @importFrom minpack.lm nlsLM
#' 
try_fit <- function(.data, .curve_func, .arguments) {
  f <- .curve_func
  fit <- try({
    new_env <- rlang::env(dat = .data)
    eval(
      rlang::expr(minpack.lm::nlsLM(data = dat, !!!.arguments)),
      envir = new_env
    )
  })
  return(fit)
}


#' A function that finds the value of `x` for a fitted `nlsModel` for a targeted `y` value
#'
#' @param .fit A fitted `nlsLM` object.
#' @param .x_vals A vector of `x` values in which to search across.
#' @param .target The target `y` value to search for.
#' @return A list produced by `uniroot`. See page for precise details.
#' @export
#' 
find_root_fit <- function(.fit, .x_vals, .target) {
  uniroot(
    f = function(x) predict(.fit, tibble(x)) - .target,
    interval = c(min(.x_vals), max(.x_vals))
  )
}

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

predict_data <- function(.fit, .lower_x, .upper_x, .group = NULL, .num_points = 1e3) {
  # full range of x values
  x_vals <- seq(.lower_x, .upper_x, length.out = .num_points)
  if (!is_null(.group)) {
    # unique group values
    .data <- .data |> mutate(group = !!ensym(.group))
    group_vals <- unique(.data$group)
    # go through each and collate
    predicted <- map(group_vals, function(i) {
      tibble(
        x = x_vals,
        y = predict(.fit, newdata = tibble(x = x_vals, group = i)),
        group = i
      )
    }) |>
      bind_rows()
  } else {
    predicted <- tibble(
      x = x_vals,
      y = predict(.fit, newdata = tibble(x = x_vals))
    )
  }
  return(predicted)
}
