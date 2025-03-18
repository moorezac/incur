#' @title Make a list of upper/lower bounds.
#' @description From a function and a list of lower and/or upper bounds to be set, create a named list with bounds set for all arguments in `.func` as provided, with all others set to Inf/-Inf for upper/lower bounds respectively.
#' @param .func A function for which to create a list of bounds from.
#' @param .lower A named list for which arguments in `.func` to set the lower bounds for.
#' @param .upper A named list for which arguments in `.func` to set the upper bounds for.
#' @return A named list of all arguments as seen in `.func` with provided bounds set, and Inf/-Inf for all other upper/lower bounds respectively.
#' @importFrom dplyr case_when
#' @importFrom purrr discard iwalk list_flatten
#' @importFrom rlang is_null
#' @importFrom stringr str_flatten_comma

make_bounds <- function(.func, .lower = NULL, .upper = NULL) {
  formal_arguments <- names(formals(.func))
  formal_arguments <- formal_arguments[!formal_arguments %in% c("x", "group")]

  # if only one is supplied
  final_list <- list(lower = .lower, upper = .upper) |> purrr::discard(rlang::is_null)

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
