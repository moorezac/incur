#' @title Make a formula expression from a function
#' @description From a provided function, create a formula in terms of `y ~ f({function arguments})`.
#' @param func A function that describes a curve in terms of `x` and individual curve parameters.
#' @return A language object in the form of `y ~ f()`
#' @importFrom stringr str_c str_flatten str_replace

make_formula <- function(func) {
  function_out <- capture.output(args(func))
  function_out <- function_out[function_out != "NULL"]
  function_out[-1] <- trimws(function_out[-1])
  function_out <- stringr::str_flatten(function_out)
  function_out <- stringr::str_replace(function_out, "function ", "f")
  function_out <- stringr::str_c("y ~", function_out)
  function_out <- str2lang(function_out)
  
  return(function_out)
}

#' @title Make a list of upper/lower bounds.
#' @description From a function and a list of lower and/or upper bounds to be set, create a named list with bounds set for all arguments in `func` as provided, with all others set to Inf/-Inf for upper/lower bounds respectively.
#' @param func A function for which to create a list of bounds from.
#' @param lower A named list for which arguments in `func` to set the lower bounds for.
#' @param upper A named list for which arguments in `func` to set the upper bounds for.
#' @param lower_bounds A named list that contains the lower bounds for specified parameters. All other lower bounds will be set at `-Inf`.
#' @param upper_bounds A named list that contains the upper bounds for specified parameters. All other upper bounds will be set at `Inf`.
#' @return A named list of all arguments as seen in `func` with provided bounds set, and Inf/-Inf for all other upper/lower bounds respectively.
#' @importFrom dplyr case_when
#' @importFrom purrr discard iwalk list_flatten
#' @importFrom rlang is_null
#' @importFrom stringr str_flatten_comma

make_bounds <- function(func, lower = NULL, upper = NULL) {
  formal_arguments <- names(formals(func))
  formal_arguments <- formal_arguments[!formal_arguments %in% c("x", "group")]
  
  # if only one is supplied
  final_list <- list(lower = lower, upper = upper) |> purrr::discard(rlang::is_null)
  
  # check bound args are there
  purrr::iwalk(final_list, function(x, i) {
    if (!any(names(x) %in% formal_arguments)) {
      missing <- names(x)[!names(x) %in% formal_arguments]
      stop(str_c(i, " arg(s) not found in : ", stringr::str_flatten_comma(missing)))
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
    
    return(full)
  })
}

#' @title Run the complete final argument pipeline.
#' @description Make a list of all arguments to be injected into `minpack.lm::nlsLM` to fit a model.
#' @param curve_func A function that describes a curve/model in terms of `x`.
#' @param start_vals A named list of starting values for all arguments/parameters within `curv_func`.
#' @return A named list to be injected into `minpack.lm::nlsLM`.
#' @importFrom purrr discard list_flatten
#' @importFrom rlang is_null

make_final_arguments <- function(curve_func, start_vals, lower_bounds = NULL, upper_bounds = NULL, dots) {
  # formula
  function_formula <- make_formula(curve_func)
  
  final_arguments <- list(
    formula = function_formula,
    start = purrr::list_flatten(start_vals)
  )
  # bounds
  if (!rlang::is_null(lower_bounds) | !rlang::is_null(upper_bounds)) {
    bounds <- make_bounds(curve_func, lower_bounds, upper_bounds)
    final_arguments <- append(final_arguments, bounds)
  }
  # dots dots dots
  final_arguments <- append(final_arguments, dots)
  
  # final final
  final_arguments <- purrr::discard(final_arguments, is_null)
  
  return(final_arguments)
}
