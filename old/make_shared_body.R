#' @title Make a shared function body.
#' @description Create a shared function body across arguments and individual groups.
#' @param .func A function for which the body is to be modified.
#' @param .group A character vector that contains all unique grouping values.
#' @param .params A character vector of arguments in `.func` that are to be shared across `.groups`.
#' @return An unevaluated `call` that contains a vectorised `case_when` to modify shared parameters based on `.group` values. This can be used to then create a new function.
#' @importFrom dplyr case_when
#' @importFrom purrr map pmap
#' @importFrom stringr str_c str_glue

make_shared_body <- function(.func, .group, .params) {
  formal_arguments <- names(formals(.func))
  unique_arguments <- formal_arguments[!formal_arguments %in% c(.params, "x")]
  append_arguments <- make_append_arguments(unique_arguments, .group)

  individual_group_list <- rep(list(.group), length(unique_arguments))

  # this is not recommended
  expression_list <- purrr::pmap(
    .l = list(a = unique_arguments, b = individual_group_list, c = append_arguments),
    .f = function(a, b, c) {
      # this creates the case_when
      # has to be vectorized - cannot use base if
      stringr::str_c(
        stringr::str_glue("{a} <- dplyr::case_when("),
        stringr::str_glue("group == '{b}' ~ {c}") |> stringr::str_c(collapse = ", "),
        ")"
      )
    }
  )
  expression_list <- purrr::map(expression_list, str2lang)

  return(as.call(c(as.name("{"), expression_list, body(.func))))
}
