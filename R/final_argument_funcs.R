#' Create Formula Expression from Model Function
#'
#' @description
#' Generate a formula object in the form `y ~ f(arguments)` from a provided
#' model function, suitable for use with nonlinear least squares fitting.
#'
#' @param model_func A function that describes a mathematical curve in terms of `x`
#'   and individual curve parameters.
#'
#' @return A language object representing the formula `y ~ f(parameters)`.
#'
#' @details
#' This is an internal helper function used to prepare arguments for `minpack.lm::nlsLM()`.
#' The function extracts the argument structure from the provided function
#' and creates an appropriate formula object.
#'
#' @family model_helpers
#' @importFrom stringr str_c str_flatten str_replace
#' @keywords internal
make_formula <- function(func) {
  function_out <- capture.output(args(func))
  function_out <- function_out[function_out != "NULL"]
  function_out[-1] <- trimws(function_out[-1])
  function_out <- stringr::str_flatten(function_out)
  function_out <- stringr::str_replace(function_out, "function ", "f")
  function_out <- stringr::str_c("y ~", function_out)

  str2lang(function_out)
}

#' Create Parameter Bounds List for Model Fitting
#'
#' @description
#' Generate named lists of upper and lower parameter bounds for nonlinear
#' model fitting, setting specified bounds and defaults for unspecified parameters.
#'
#' @param model_func A function defining the mathematical model.
#' @param lower_bounds A named list specifying lower bounds for specific parameters.
#'   Unspecified parameters default to `-Inf`.
#' @param upper_bounds A named list specifying upper bounds for specific parameters.
#'   Unspecified parameters default to `Inf`.
#'
#' @return A named list containing `lower` and `upper` elements with bounds
#'   for all parameters in `model_func`.
#'
#' @details
#' This function ensures all parameters in the model function have defined bounds,
#' either from user specification or from defaults. Parameter names must match
#' those in the model function arguments.
#'
#' @family model_helpers
#' @importFrom dplyr case_when
#' @importFrom purrr discard iwalk list_flatten
#' @importFrom rlang is_null
#' @importFrom stringr str_flatten_comma str_c
#' @keywords internal
make_bounds <- function(func, lower_x = NULL, upper_x = NULL) {
  formal_arguments <- names(formals(func))
  formal_arguments <- formal_arguments[!formal_arguments %in% c("x", "group")]

  # if only one is supplied
  final_list <- list(lower = lower_x, upper = upper_x) |>
    purrr::discard(rlang::is_null)

  # check bound args are there
  purrr::iwalk(final_list, function(x, i) {
    if (!any(names(x) %in% formal_arguments)) {
      missing <- names(x)[!names(x) %in% formal_arguments]
      stop(str_c(
        i,
        " arg(s) not found in : ",
        stringr::str_flatten_comma(missing)
      ))
    }
  })

  # the bounds need to be in order
  purrr::imap(final_list, function(x, i) {
    full <- dplyr::case_when(
      i == "lower" ~ rep(-Inf, length(formal_arguments)),
      i == "upper" ~ rep(Inf, length(formal_arguments))
    )
    names(full) <- formal_arguments

    replace <- unlist(purrr::list_flatten(x))
    full[names(full) %in% names(replace)] <- replace

    full
  })
}

#' Prepare Final Arguments for Model Fitting
#'
#' @description
#' Assemble all necessary arguments for `minpack.lm::nlsLM()` including formula,
#' starting values, and optional parameter bounds.
#'
#' @param model_func A function defining the mathematical model.
#' @param start_values A named list of starting parameter values.
#' @param lower_bounds Optional named list of lower parameter bounds.
#' @param upper_bounds Optional named list of upper parameter bounds.
#'
#' @return A named list ready for injection into `minpack.lm::nlsLM()` containing:
#' \describe{
#'   \item{formula}{Model formula object}
#'   \item{start}{Starting parameter values}
#'   \item{lower}{Lower parameter bounds (if specified)}
#'   \item{upper}{Upper parameter bounds (if specified)}
#' }
#'
#' @family model_helpers
#' @importFrom purrr discard list_flatten
#' @importFrom rlang is_null
#' @keywords internal
prepare_final_arguments <- function(
  model_func,
  start_values,
  lower_bounds = NULL,
  upper_bounds = NULL
) {
  # formula
  function_formula <- make_formula(model_func)

  final_arguments <- list(
    formula = function_formula,
    start = purrr::list_flatten(start_values)
  )
  # bounds
  if (!rlang::is_null(lower_bounds) | !rlang::is_null(upper_bounds)) {
    bounds <- make_bounds(model_func, lower_bounds, upper_bounds)
    final_arguments <- append(final_arguments, bounds)
  }

  # final final
  final_arguments <- purrr::discard(final_arguments, rlang::is_null)

  final_arguments
}
