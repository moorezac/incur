#' @title Run the complete final argument pipeline.
#' @description From a `nls` object, detect outliers via the ROUT method, and then re-fit with outliers removed.
#' @param .data A `data.frame` or `data.frame` extension (tibble) in long format.
#' @param .x_var An quoted/unquoted argument that refers to the the `x` value within `.data`.
#' @param .y_var An quoted/unquoted argument that refers to the the `y` value within `.data`.
#' @param .fit A fitted `nls` object.
#' @param .curve_func A function that describes a curve/model in terms of `x`. See `incur::incur_models` for common examples. Alternatively, see below for examples.
#' @param .start_func A function with arguments `x` and `y` that returns a named list of starting values for all arguments/parameters within `.curv_func`. See `incur::incur_models` for common examples. Alternatively, see below for example.
#' @param .lower_bounds A named list that contains the lower bounds for specified parameters. All other lower bounds will be set at `-Inf`.
#' @param .upper_bounds A named list that contains the upper bounds for specified parameters. All other upper bounds will be set at `Inf`.
#' @param .dots Other arguments to be passed to `minpack.lm::nlsLM`, such as `control` or `weights`.
#' @return An list data contains in order: 
#' \itemize{
#'    \item A fitted `nls` object.
#'    \item Original `.data` with modified columns to represent outliers
#'  } 
#'  @importFrom dplyr filter if_else mutate pull relocate row_number
#'  @importFrom purrr map_vec
#'  @importFrom rlang as_name enquo ensym sym
#'  @importFrom stringr str_c
#'  
detect_outliers <- function(.data, .x_var, .y_var, .fit, .curve_func, .start_func, .upper_bounds, .lower_bounds, .dots) {
  # which points are outliers within the fit
  outlier_indices <- find_outlier_indices(.fit, .scale_method = "mad")

  # name a new column based on y_var
  outlier_column <- stringr::str_c("outlier", rlang::as_name(rlang::enquo(.y_var)), sep = "_")
  # add in
  .data <- dplyr::mutate(.data, !!outlier_column := purrr::map_vec(dplyr::row_number(), function(x) {
    dplyr::if_else(x %in% outlier_indices, TRUE, FALSE)
  }))
  .data <- dplyr::relocate(.data, !!outlier_column, .after = y)

  # if we have identified no outliers we can skip to the end
  if (all(.data[[outlier_column]] == FALSE)) {
    final_fit <- .fit
  } else {
    # this time we need to filter
    .start_vals_filtered <- .start_func(
      .data |> dplyr::filter(!!rlang::sym(outlier_column) == FALSE) |> dplyr::pull(x),
      .data |> dplyr::filter(!!rlang::sym(outlier_column) == FALSE) |> dplyr::pull(y)
    )
    
    # new arguments with outliers removed
    final_arguments <- create_final_arguments(.curve_func, .start_vals_filtered, .lower_bounds, .upper_bounds, .dots)

    # new fit with outliers removed
    fit_outlier_removed <- try_fit(dplyr::filter(.data, !!rlang::ensym(outlier_column) == FALSE), .curve_func, final_arguments)

    # if this didn't fit then return original
    if (inherits(fit_outlier_removed, "try-error")) {
      message("error in the model fit with outliers removed")
      message("returning the model fit with outliers included")
      final_fit <- .fit
    } else {
      final_fit <- fit_outlier_removed
    }
  }
  
  return(list(fit = final_fit, data = .data))
}
