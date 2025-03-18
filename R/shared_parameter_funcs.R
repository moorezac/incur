#' @title Make a list of shared starting values.
#' @description If modelling a curve with shared parameters across individual groups, create a named list of starting values for both sets of groups to improve the fitting process.,
#' @param .data A dataframe or dataframe-like object.
#' @param .group A character vector that contains all unique grouping values as seen in `.data`.
#' @param .start_func  function with arguments `x` and `y` that returns a named list of starting values.
#' @param .params A character vector that indicates which parameter(s) are to be shared across groups.
#' @return A named list of starting values for provided arguments and group combinations.
#' @importFrom dplyr case_when filter pull
#' @importFrom purrr map imap list_flatten
#' @importFrom stringr str_c
#' 
make_shared_start_vals <- function(.data, .start_func, .group, .params) {
  # use all data for shared params
  start_vals_global <- .start_func(
    x = .data |> dplyr::pull(x),
    y = .data |> dplyr::pull(y)
  )
  
  start_vals_group <- purrr::map(.group, function(i) {
    .start_func(
      x = .data |> dplyr::filter(group == i) |> dplyr::pull(x),
      y = .data |> dplyr::filter(group == i) |> dplyr::pull(y)
    )
  })
  names(start_vals_group) <- .group
  
  # make names
  start_vals_group <- purrr::imap(start_vals_group, function(x, i) {
    names(x) <- dplyr::case_when(
      !names(x) %in% .params ~ stringr::str_c(names(x), "_", i),
      .default = names(x)
    )
    # filter to group values
    x[!names(x) %in% names(start_vals_global)]
  })
  
  # append
  start_vals_global <- start_vals_global[names(start_vals_global) %in% .params]
  
  return(append(start_vals_global, purrr::list_flatten(start_vals_group, name_spec = "{inner}")))
}

#' @title Make shared function formals across arguments and individual groups.
#' @description From a provided function, modify arguments so that shared parameters are replaced in the form `{argument}_{group}`.
#' @param .func A function for which arguments are to be modified.
#' @param .group A character vector that contains all unique grouping values.
#' @param .params A character vector of arguments in `.func` that are to be shared across `.groups`.
#' @return An `pairlist` object that that can be used to create a new function
#' @importFrom rlang expr
#' @importFrom purrr map
#' @importFrom stringr str_flatten_comma
#' 
make_shared_formals <- function(.func, .group, .params) {
  # easier to make this a function
  make_append_arguments <- function(.unique, .group) {
    len <- length(.group)
    purrr::map(.unique, str_c, "_", .group)
  }
  
  formal_arguments <- names(formals(.func))
  unique_arguments <- formal_arguments[!formal_arguments %in% c(.params, "x")]
  append_arguments <- make_append_arguments(unique_arguments, .group)
  
  if (!any(.params %in% formal_arguments)) {
    missing <- .params[!.params %in% formal_arguments]
    stop(str_c("model params not found: ", stringr::str_flatten_comma(missing)))
  }
  
  final_arguments <- c(
    "x",
    formal_arguments[formal_arguments %in% .params],
    unlist(append_arguments),
    "group"
  )
  
  # use this to create alist
  # TODO is this jank?
  final_formals <- rep(list(rlang::expr()), length(final_arguments))
  names(final_formals) <- final_arguments
  
  return(as.pairlist(final_formals))
}

#' @title Make a shared function body.
#' @description Create a shared function body across arguments and individual groups.
#' @param .func A function for which the body is to be modified.
#' @param .group A character vector that contains all unique grouping values.
#' @param .params A character vector of arguments in `.func` that are to be shared across `.groups`.
#' @return An unevaluated `call` that contains a vectorised `case_when` to modify shared parameters based on `.group` values. This can be used to then create a new function.
#' @importFrom dplyr case_when
#' @importFrom purrr map pmap
#' @importFrom stringr str_c str_glue
#' 
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

#' @title Run the complete shared parameters pipeline .
#' @description Create a new list of starting values, and custom function that allows for shared parameters across groups.
#' @param .func A function for which arguments are to be modified.
#' @param .group A character vector that contains all unique grouping values.
#' @param .params A character vector of arguments in `.func` that are to be shared across `.groups`.
#' @return An list data contains in order: 
#' \itemize{
#'    \item Original `.data` with modified columns to represent groups.
#'    \item A custom `.curve_func` with modified formals and body.
#'    \item A named list of new starting values for each group.
#'  } 
#' @importFrom dplyr mutate relocate
#' @importFrom rlang ensym new_function
#' 
make_shared_parameters <- function(.data, .shared_params, .shared_group, .curve_func, .start_func) {
  # asserts
  if ("group" %in% colnames(.data)) {
    message("detected `group` as a col. name in provided data")
    message("`group` is a reserved col. name in incur")
    message("renaming `group` to `group_original`")
    .data <- dplyr::mutate(.data, group_orignal = group)
  }
  
  # mutate in place
  .data <- dplyr::mutate(.data, group = !!rlang::ensym(.shared_group))
  .data <- dplyr::mutate(.data, group = as.character(group))
  .data <- dplyr::relocate(.data, group, .after = y)
  
  group_vec <- unique(.data$group)
  
  .start_vals <- make_shared_start_vals(.data, .start_func, group_vec, .shared_params)
  
  .shared_arguments <- make_shared_formals(.curve_func, group_vec, .shared_params)
  .shared_body <- make_shared_body(.curve_func, group_vec, .shared_params)
  .curve_func <- rlang::new_function(as.pairlist(.shared_arguments), .shared_body)
  
  return(list(.data, .curve_func, .start_vals))
}
