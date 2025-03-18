#' @title Find outliers within a fitted `nls` object.
#' @description An implementation of the ROUT method described in Motulsky and Brown (2006).
#' @param .fit A fitted `nls` object.
#' @param .q The maximum desired FDR.
#' @param .scale_method The method in which to scale the residuals. The original paper uses the `quantile` method. The `MAD` implementation is used within the `dr4pl` package.
#' @return A numeric vector of indices individual data points as outliers. This corresponds to the row number of the data.frame used in the original fit.
#' @export
#' 
find_outlier_indices <- function(.fit, .q = 1e-5, .scale_method = "quantile") {
  residuals_vec <- resid(.fit)
  n <- length(residuals_vec)

  # dr4pl implementation
  if (.scale_method == "mad") {
    scale <- mad(residuals_vec)
  }
  # original paper
  if (.scale_method == "quantile") {
    scale <- quantile(
      x = residuals_vec,
      probs = pnorm(1, 0, 1) - pnorm(-1, 0, 1)
    ) * (n / (n - length(coef(.fit))))
  }

  abs_res_sorted <- sort(abs(residuals_vec), index.return = TRUE)$x
  indices_sorted <- sort(abs(residuals_vec), index.return = TRUE)$ix

  alphas <- .q * seq(from = n, to = 1, by = -1) / n
  p_values <- 2 * pt(
    q = abs_res_sorted / scale,
    df = n - length(coef(.fit)),
    lower.tail = FALSE
  )

  indices_fdr <- which(p_values < alphas)

  if (length(indices_fdr) == 0) {
    indices_outlier <- NULL
  } else {
    indices_outlier <-
      indices_sorted[seq(from = min(indices_fdr), to = n, by = 1)]
  }

  return(indices_outlier)
}


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
