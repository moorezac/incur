#' @title IRWLS via Huber loss function.
#' @description Perform iterative-reweighted least squares regression via Huber loss function.
#' @param fit A fitted `nlsModel` object
#' @param .iter_max The number of iterations to perform
#' @param .k The tuning constant to be used.
#' @param .tol Tolerance level required to be achieved before stopping.
#' @return A fitted `nlsModel` object with updated coefficents.
#' @importFrom dplyr case_when
#' 
huber <- function(fit, arguments, data, curve_func, iter_max = 100, k = 1.345, tol = 1e-6) {
  iter <- 0
  converged <- FALSE
  fit_old <- fit
  
  # loop
  while (iter < iter_max && !converged) {
    coef_old <- coef(fit_old)
    resids <- residuals(fit_old)
    # mad
    s <- median(abs(resids - median(resids))) * 1.4826
    u <- resids / (s + 1e-10)
    weights_vec <- dplyr::case_when(
      abs(u) <= k ~ 1,
      .default = k / abs(u)
    )
    
    arguments_new <- append(arguments, list(weights = weights_vec))
    
    fit_new <- try_fit(data, curve_func, arguments_new)
    
    if (inherits(fit_new, "try-error")) {
      message("error in huber")
      return(fit_old)
    }
    
    coef_new <- coef(fit_new)
    
    # check convergence
    coef_change <- max(abs(coef_new - coef_old) / (abs(coef_old) + 1e-6))
    converged <- coef_change < tol
    
    # update for next iteration
    fit_old <- fit_new
    iter <- iter + 1
  }
  # message(str_glue("achieved {tol} tolerance in {iter} iterations"))
  fit <- fit_new
  
  return(fit)
}

#' @title Find outliers within a fitted `nls` object.
#' @description An implementation of the ROUT method described in Motulsky and Brown (2006).
#' @param fit A fitted `nls` object.
#' @param .q The maximum desired FDR.
#' @param scale_method The method in which to scale the residuals. The original paper uses the `quantile` method. The `MAD` implementation is used within the `dr4pl` package.
#' @return A numeric vector of indices individual data points as outliers. This corresponds to the row number of the data.frame used in the original fit.
#' @export
#' 
find_outlier_indices <- function(fit, q = 1e-5, scale_method = "quantile") {
  residuals_vec <- resid(fit)
  n <- length(residuals_vec)
  
  # dr4pl implementation
  if (scale_method == "mad") {
    scale <- mad(residuals_vec)
  }
  # original paper
  if (scale_method == "quantile") {
    scale <- quantile(
      x = residuals_vec,
      probs = pnorm(1, 0, 1) - pnorm(-1, 0, 1)
    ) * (n / (n - length(coef(fit))))
  }
  
  abs_res_sorted <- sort(abs(residuals_vec), index.return = TRUE)$x
  indices_sorted <- sort(abs(residuals_vec), index.return = TRUE)$ix
  
  alphas <- q * seq(from = n, to = 1, by = -1) / n
  p_values <- 2 * pt(
    q = abs_res_sorted / scale,
    df = n - length(coef(fit)),
    lower.tail = FALSE
  )
  
  indices_fdr <- which(p_values < alphas)
  
  if (length(indices_fdr) == 0) {
    indices_outlier <- NULL
  } else {
    indices_outlier <- indices_sorted[seq(from = min(indices_fdr), to = n, by = 1)]
  }
  
  return(indices_outlier)
}

#' @title Run the complete final argument pipeline.
#' @description From a `nls` object, detect outliers via the ROUT method, and then re-fit with outliers removed.
#' @param data A `data.frame` or `data.frame` extension (tibble) in long format.
#' @param x_var An quoted/unquoted argument that refers to the the `x` value within `data`.
#' @param y_var An quoted/unquoted argument that refers to the the `y` value within `data`.
#' @param fit A fitted `nls` object.
#' @param curve_func A function that describes a curve/model in terms of `x`. See `incur::incur_models` for common examples. Alternatively, see below for examples.
#' @param start_func A function with arguments `x` and `y` that returns a named list of starting values for all arguments/parameters within `.curv_func`. See `incur::incur_models` for common examples. Alternatively, see below for example.
#' @param lower_bounds A named list that contains the lower bounds for specified parameters. All other lower bounds will be set at `-Inf`.
#' @param upper_bounds A named list that contains the upper bounds for specified parameters. All other upper bounds will be set at `Inf`.
#' @param dots Other arguments to be passed to `minpack.lm::nlsLM`, such as `control` or `weights`.
#' @return An list data contains in order: 
#' \itemize{
#'    \item A fitted `nls` object.
#'    \item Original `data` with modified columns to represent outliers
#'  } 
#' @importFrom dplyr filter if_else mutate pull relocate row_number
#' @importFrom purrr map_vec
#' @importFrom rlang as_name enquo ensym sym
#' @importFrom stringr str_c
#'  
detect_outliers <- function(data, x_var, y_var, fit, curve_func, start_func, upper_bounds, lower_bounds, dots) {
  # which points are outliers within the fit
  outlier_indices <- find_outlier_indices(fit, scale_method = "mad")
  
  # name a new column based on y_var
  # outlier_column <- stringr::str_c("outlier", rlang::as_name(rlang::enquo(y_var)), sep = "_")
  outlier_column <- stringr::str_c("outlier", y_var, sep = "_")
  # add in
  data <- dplyr::mutate(data, !!outlier_column := purrr::map_vec(dplyr::row_number(), function(x) {
    dplyr::if_else(x %in% outlier_indices, TRUE, FALSE)
  }))
  data <- dplyr::relocate(data, !!outlier_column, .after = y)
  
  # if we have identified no outliers we can skip to the end
  if (all(data[[outlier_column]] == FALSE)) {
    final_fit <- fit
  } else {
    # this time we need to filter
    start_vals_filtered <- start_func(
      data |> dplyr::filter(!!rlang::sym(outlier_column) == FALSE) |> dplyr::pull(x),
      data |> dplyr::filter(!!rlang::sym(outlier_column) == FALSE) |> dplyr::pull(y)
    )
    
    # new arguments with outliers removed
    final_arguments <- make_final_arguments(curve_func, start_vals_filtered, lower_bounds, upper_bounds, dots)
    
    # new fit with outliers removed
    fit_outlier_removed <- try_fit(dplyr::filter(data, !!rlang::ensym(outlier_column) == FALSE), curve_func, final_arguments)
    
    # if this didn't fit then return original
    if (inherits(fit_outlier_removed, "try-error")) {
      message("error in the model fit with outliers removed")
      message("returning the model fit with outliers included")
      final_fit <- fit
    } else {
      final_fit <- fit_outlier_removed
    }
  }
  
  return(list(fit = final_fit, data = data))
}
