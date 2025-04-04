#' @title Fit a curve with `incur`.
#' @description Given a set of data, fit an curve via `minpack.lm` with options for shared parameters, bounds, and automated outlier detection and removal.
#' @param .data A `data.frame` or `data.frame` extension (tibble) in long format.
#' @param .x_var An quoted/unquoted argument that refers to the the `x` value within `.data`.
#' @param .y_var An quoted/unquoted argument that refers to the the `y` value within `.data`.
#' @param .curve_func A function that describes a curve/model in terms of `x`. See `incur::incur_models` for common examples. Alternatively, see below for examples.
#' @param .start_func A function with arguments `x` and `y` that returns a named list of starting values for all arguments/parameters within `.curv_func`. See `incur::incur_models` for common examples. Alternatively, see below for example.
#' @param .start_vals A named list of starting values for all arguments/parameters within `curv_func`. This cannot be used in conjunction with `.detect_outliers`.
#' @param .huber Perform iterative reweighted least squares non-linear regression using the Huber loss function.
#' @param .detect_outliers Boolean to indicate whether to filter outliers as described in Motulsky and Brown (2006). It is highly recommended that the `.huber` argument is set to true in order to start with "robust" regression.
#' @param .shared_group For nested fits (shared parameters), the quoted/unquoted argument for the column in `.data` that refers to individual groups to share fitted parameters.
#' @param .shared_params For nested fits (shared parameters), a character vector that indicates which parameter(s) in `.curve_func` are to be shared across groups.
#' @param .lower_bounds A named list that contains the lower bounds for specified parameters. All other lower bounds will be set at `-Inf`.
#' @param .upper_bounds A named list that contains the upper bounds for specified parameters. All other upper bounds will be set at `Inf`.
#' @param .return_func For nested fits (shared parameters), a boolean to indicate whether to return a modified `.curve_func` used in this process.
#' @param ... Other arguments to be passed to `minpack.lm::nlsLM`, such as `control` or `weights`.
#' @return A named list containing:
#'  \itemize{
#'    \item `fit`: the fitted `nlsModel` object.
#'    \item `data`: the original data used in the fit. If `.detect_outliers` is true then column added in the format `outlier_{y_var}` will be added to indicate which points are detected as outliers.
#'  }
#'  @importFrom checkmate assert_character assert_data_frame assert_function assert_list assert_logical
#'  @importFrom dplyr mutate relocate
#'  @importFrom rlang is_null
#'  @importFrom tidyr drop_na
#' @export
#' @examples
#' # exponential plateau
# func <- function(x, y0, ym, k) ym - (ym - y0) * exp(-k * x)
# func_start <- function(x, y) list(ym = max(y), y0 = min(y), k = 1)
#'
#' # fit the curve to the data
#' # share y0 and k across state and force y0 > 0
# fit <- fit_model(
#   .data = Puromycin,
#   .x_var = "conc",
#   .y_var = "rate",
#   .curve_func = func,
#   .start_func = func_start,
#   .detect_outliers = TRUE,
#   .lower_bounds = list(y0 = 0),
#   .shared_group = "state",
#   .shared_params = c("y0", "k"),
#   control = minpack.lm::nls.lm.control(maxiter = 1e3)
# )
#' # create data from each group
#' predicted <- map(unique(Puromycin$state), function(i) {
#'   # more data points = smoother curve
#'   # this is easier than a geom_function approach
#'   x_vals <- seq(min(fit$data$x), max(fit$data$x), length.out = 1e4)
#'   # create the data and return
#'   tibble(
#'     x = x_vals,
#'     y = predict(res$fit, newdata = tibble(x = x_vals, group = i)),
#'     group = i
#'   )
#' })
#' # collate
#' predicted <- bind_rows(predicted)
#' # plot
#' ggplot(mapping = aes(x, y, colour = group)) +
#'   geom_point(data = fit$data) +
#'   geom_line(data = predicted)
fit_model <- function(.data, .x_var, .y_var, .curve_func, .start_func = NULL, .start_vals = NULL, .huber = FALSE, .detect_outliers = FALSE, .shared_group = NULL, .shared_params = NULL, .lower_bounds = NULL, .upper_bounds = NULL, .return_func = FALSE, ...) {
  # dots dots dots
  .dots <- list(...)
  
  # asserts
  checkmate::assert_data_frame(as.data.frame(.data))
  checkmate::assert_character(.x_var)
  checkmate::assert_character(.y_var)
  checkmate::assert_function(.curve_func)
  if (!rlang::is_null(.start_func)) checkmate::assert_function(.start_func)
  if (!rlang::is_null(.start_vals)) checkmate::assert_list(.start_vals)
  if (!rlang::is_null(.huber)) checkmate::assert_logical(.huber)
  if (!rlang::is_null(.detect_outliers)) checkmate::assert_logical(.detect_outliers)
  if (!rlang::is_null(.shared_group)) checkmate::assert_character(.shared_group)
  if (!rlang::is_null(.shared_params)) checkmate::assert_character(.shared_params)
  if (!rlang::is_null(.lower_bounds)) checkmate::assert_list(.lower_bounds)
  if (!rlang::is_null(.upper_bounds)) checkmate::assert_list(.upper_bounds)
  if (!rlang::is_null(.return_func)) checkmate::assert_logical(.return_func)
  
  # checks
  if (rlang::is_null(.start_func) & !rlang::is_null(.detect_outliers)) {
    stop("if detecting outliers a function needs to be provided to generate starting values")
  }
  if (xor(rlang::is_null(.shared_group), rlang::is_null(.shared_params))) {
    stop("'.shared_group' and '.shared_params' need to be provided together")
  }
  
  # mutate data to x and y
  .data <- dplyr::mutate(.data, x = !!rlang::ensym(.x_var), y = !!rlang::ensym(.y_var))
  .data <- tidyr::drop_na(.data, x, y)
  .data <- dplyr::mutate(.data, x = as.numeric(x), y = as.numeric(y))
  .data <- dplyr::relocate(.data, x, y)
  
  # start values
  if (rlang::is_null(.start_vals)) {
    .start_vals <- .start_func(x = .data$x, y = .data$y)
  }

  # shared curve parameters
  if (rlang::is_null(.shared_group) && rlang::is_null(.shared_params)) {
    shared_parameter_list <- make_shared_parameters(.data, .shared_params, .shared_group, .curve_func, .start_func)
    .data <- shared_parameter_list[[1]]
    .curve_func <- shared_parameter_list[[2]]
    .start_vals <- shared_parameter_list[[3]]
  }

  # final_arguments
  final_arguments <- make_final_arguments(.curve_func, .start_vals, .lower_bounds, .upper_bounds, .dots)

  fit <- try_fit(.data, .curve_func, final_arguments)

  if (inherits(fit, "try-error")) {
    stop("error in initial model fit")
  }

  if (.huber) {
    fit <- huber(fit, final_arguments, .iter_max = 100, .k = 1.345, .tol = 1e-6)
  }
  
  if (!.detect_outliers) {
    final_fit <- fit
  } else {
    outlier_list <- detect_outliers(.data, .x_var, .y_var, fit, .curve_func, .start_func, .upper_bounds, .lower_bounds, .dots)
    final_fit <- outlier_list[[1]]
    .data <- outlier_list[[2]]
  }
  
  # return
  if (.return_func) {
    return(list(fit = final_fit, data = .data, func = .curve_func))
  } else {
    return(list(fit = final_fit, data = .data))
  }
}
