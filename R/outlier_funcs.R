#' Apply Huber Robust Iteratively Reweighted Least Squares
#'
#' @description
#' Apply Huber robust regression to reduce the influence of outliers on model fits
#' using iteratively reweighted least squares (IRWLS) algorithm.
#'
#' @param fit Initial nls model fit object.
#' @param arguments Arguments list for nlsLM function.
#' @param data Data frame containing the fitting data.
#' @param model_func Function defining the mathematical model.
#' @param huber_opts List of Huber regression options containing:
#'   \itemize{
#'     \item `iter_max`: Maximum number of iterations (default: 100)
#'     \item `k`: Huber constant for weight calculation (default: 1.345)
#'     \item `tol`: Convergence tolerance (default: 1e-6)
#'   }
#'
#' @return Updated nls fit object with robust parameter estimates.
#'
#' @details
#' The Huber robust regression method reduces the influence of outliers by
#' applying weights to data points based on their residuals. Points with
#' large residuals receive lower weights in subsequent fitting iterations.
#'
#' @family outlier_detection
#' @importFrom dplyr case_when
#' @importFrom rlang is_null
#' @keywords internal
huber_irwls <- function(
  fit,
  arguments,
  data,
  model_func,
  huber_opts = NULL
) {
  if (rlang::is_null(huber_opts)) {
    iter_max <- 100
    k <- 1.345
    tol <- 1e-6
  } else {
    iter_max <- huber_opts$iter_max
    k <- huber_opts$k
    iter_tol <- huber_opts$tol
  }

  iter <- 0
  converged <- FALSE
  current_fit <- fit

  while (iter < iter_max && !converged) {
    old_coef <- coef(current_fit)
    residuals_vec <- residuals(current_fit)

    # Calculate robust weights
    s <- median(abs(residuals_vec - median(residuals_vec))) * 1.4826
    u <- residuals_vec / (s + 1e-10)
    weights <- dplyr::case_when(
      abs(u) <= k ~ 1,
      .default = k / abs(u)
    )

    # Fit with weights
    weighted_args <- append(arguments, list(weights = weights))
    new_fit <- try_fit(data, model_func, weighted_args)

    if (inherits(new_fit, "try-error")) {
      warning("Error in Huber iteration, returning previous fit")
      return(current_fit)
    }

    # Check convergence
    new_coef <- coef(new_fit)
    coef_change <- max(abs(new_coef - old_coef) / (abs(old_coef) + 1e-6))
    converged <- coef_change < tol

    current_fit <- new_fit
    iter <- iter + 1
  }

  current_fit
}

#' Find Outlier Indices Using ROUT Method
#'
#' @description
#' Identify outliers in model residuals using the ROUT (Robust regression and
#' Outlier removal) method, which controls the False Discovery Rate (FDR).
#'
#' @param fit A fitted nls model object.
#' @param rout_opts List of ROUT options containing:
#'   \itemize{
#'     \item `q`: Maximum desired False Discovery Rate (default: 1e-3)
#'     \item `scale_method`: Method for scaling residuals, either "mad" or "quantile" (default: "mad")
#'   }
#'
#' @return Numeric vector of outlier indices, or NULL if no outliers detected.
#'
#' @details
#' The ROUT method uses a sequential testing procedure that controls the FDR
#' at the specified level `q`. It's more appropriate than traditional methods
#' when testing multiple points for outlier status.
#'
#' @family outlier_detection
#' @importFrom stats mad quantile pt
#' @export
#'
#' @examples
#' \dontrun{
#' # Identify outliers in fitted model
#' outlier_indices <- find_rout_indices(
#'   fit = fitted_model,
#'   rout_opts = list(q = 0.01, scale_method = "mad")
#' )
#' }
find_rout_indices <- function(
  fit,
  rout_opts = NULL
) {
  if (rlang::is_null(rout_opts)) {
    q <- 1e-3
    scale_method <- "mad"
  } else {
    q <- rout_opts$q
    scale_method <- rout_opts$scale_method
  }

  residuals_vec <- resid(fit)
  n <- length(residuals_vec)

  # Calculate scale factor
  if (scale_method == "mad") {
    scale_factor <- mad(residuals_vec)
  } else if (scale_method == "quantile") {
    scale_factor <- quantile(
      residuals_vec,
      probs = pnorm(1, 0, 1) - pnorm(-1, 0, 1)
    ) *
      (n / (n - length(coef(fit))))
  } else {
    stop("scale_method must be 'mad' or 'quantile'")
  }

  # Sort residuals and get indices
  abs_res_sorted <- sort(abs(residuals_vec), index.return = TRUE)$x
  indices_sorted <- sort(abs(residuals_vec), index.return = TRUE)$ix

  # Calculate FDR-adjusted p-values
  alphas <- q * seq(from = n, to = 1, by = -1) / n
  p_values <- 2 *
    pt(
      q = abs_res_sorted / scale_factor,
      df = n - length(coef(fit)),
      lower.tail = FALSE
    )

  # Find outliers
  indices_fdr <- which(p_values < alphas)

  if (length(indices_fdr) == 0) {
    NULL
  } else {
    indices_sorted[seq(from = min(indices_fdr), to = n, by = 1)]
  }
}

#' Detect Outliers and Refit Model
#'
#' @description
#' Internal function that identifies outliers using the ROUT method and refits
#' the model with outliers removed. Adds outlier flags to the data.
#'
#' @param data Input data frame.
#' @param x_var Character string for x variable name.
#' @param y_var Character string for y variable name.
#' @param model_func Function defining the mathematical model.
#' @param start_func Starting values function.
#' @param rout_opts ROUT method options.
#' @param share_group Grouping variable for shared parameters.
#' @param share_params Vector of shared parameter names.
#' @param upper_bounds Upper parameter bounds.
#' @param lower_bounds Lower parameter bounds.
#' @param fit Initial fitted model.
#'
#' @return List containing updated fit object and data with outlier flags.
#'
#' @family outlier_detection
#' @importFrom dplyr mutate filter relocate row_number if_else pull
#' @importFrom purrr map_vec
#' @importFrom stringr str_c
#' @importFrom rlang sym is_null
#' @keywords internal
detect_outliers_and_refit <- function(
  data,
  x_var,
  y_var,
  model_func,
  start_func,
  rout_opts,
  share_group,
  share_params,
  upper_bounds,
  lower_bounds,
  fit
) {
  # Identify outliers
  outlier_indices <- find_rout_indices(fit, rout_opts)

  # Add outlier column to data
  outlier_column <- stringr::str_c("outlier_", y_var)
  data <- dplyr::mutate(
    data,
    !!outlier_column := purrr::map_vec(dplyr::row_number(), function(x) {
      dplyr::if_else(x %in% outlier_indices, TRUE, FALSE)
    })
  )
  data <- dplyr::relocate(data, !!outlier_column, .after = y)

  # If no outliers, return original fit
  if (all(data[[outlier_column]] == FALSE)) {
    return(list(fit, data))
  }

  # Refit without outliers
  data_filtered <- dplyr::filter(data, !!rlang::sym(outlier_column) == FALSE)

  # Generate new starting values
  if (!rlang::is_null(share_group) && !rlang::is_null(share_params)) {
    group_vec <- unique(data$group)

    start_values_filtered <- make_shared_start_values(
      data,
      model_func,
      start_func,
      group_vec,
      share_params
    )
  } else {
    start_values_filtered <- start_func(
      x = dplyr::pull(data_filtered, x),
      y = dplyr::pull(data_filtered, y)
    )
  }

  # Create final arguments
  final_arguments <- prepare_final_arguments(
    model_func,
    start_values_filtered,
    lower_bounds,
    upper_bounds
  )

  fit_outlier_removed <- try_fit(data_filtered, model_func, final_arguments)

  if (inherits(fit_outlier_removed, "try-error")) {
    message("Error fitting model without outliers, returning original fit")
    list(fit, data)
  } else {
    list(fit_outlier_removed, data)
  }
}
