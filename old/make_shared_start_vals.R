#' @title Make a list of shared starting values.
#' @description If modelling a curve with shared parameters across individual groups, create a named list of starting values for both sets of groups to improve the fitting process.,
#' @param .data A dataframe or dataframe-like object.
#' @param .group A character vector that contains all unique grouping values as seen in `.data`.
#' @param .start_func  function with arguments `x` and `y` that returns a named list of starting values.
#' @param .params A character vector that indicates which parameter(s) are to be shared across groups.
#' @importFrom dplyr case_when filter pull
#' @importFrom purrr map imap list_flatten
#' @importFrom stringr str_c

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
