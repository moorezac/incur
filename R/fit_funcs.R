#' Find X Value for Target Y Value from Fitted Model
#'
#' @description
#' Use root-finding to determine the x value that produces a specified y value
#' from a fitted nonlinear model. Commonly used to find IC50, EC50, or other
#' benchmark concentrations.
#'
#' @param fit A fitted `nls` object from `fit_model()` or similar.
#' @param x_values A numeric vector defining the search range for x values.
#' @param target The target y value to solve for.
#'
#' @return A list object from `uniroot()` containing:
#' \describe{
#'   \item{root}{The x value that produces the target y value}
#'   \item{f.root}{The function value at the root (should be near zero)}
#'   \item{iter}{Number of iterations used}
#'   \item{init.it}{Number of initial iterations}
#'   \item{estim.prec}{Estimated precision}
#' }
#'
#' @family model_analysis
#' @importFrom stats uniroot predict
#' @export
#'
#' @examples
#' \dontrun{
#' # Find IC50 (concentration giving 50% inhibition)
#' ic50_result <- find_x_for_y(
#'   fit = dose_response_fit,
#'   x_values = c(-9, -4),  # Log concentration range
#'   target = 50  # 50% viability
#' )
#' ic50_value <- ic50_result$root
#' }
find_x_for_y <- function(fit, x_values, target) {
  uniroot(
    f = function(x) {
      predict(fit, data.frame(x = x)) - target
    },
    interval = c(min(x_values), max(x_values))
  )
}

#' Calculate Area Under Curve for Fitted Model
#'
#' @description
#' Calculate the area under the curve (AUC) for a fitted model between specified
#' bounds using numerical integration. Useful for pharmacokinetic analysis or
#' measuring total response over time.
#'
#' @param fit A fitted `nls` object from `fit_model()` or similar.
#' @param lower_x The lower bound for integration.
#' @param upper_x The upper bound for integration.
#'
#' @return A list of class `integrate` containing:
#' \describe{
#'   \item{value}{The calculated area under the curve}
#'   \item{abs.error}{Estimate of absolute error}
#'   \item{subdivisions}{Number of subdivisions used}
#'   \item{message}{Status message}
#' }
#'
#' @family model_analysis
#' @importFrom stats integrate predict
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate total growth over time
#' auc_result <- calculate_auc(
#'   fit = growth_fit,
#'   lower_x = 0,
#'   upper_x = 72  # 72 hours
#' )
#' total_growth <- auc_result$value
#' }
calculate_auc <- function(fit, lower_x, upper_x) {
  integrate(
    f = function(x) {
      predict(fit, tibble(x = x))
    },
    lower = lower_x,
    upper = upper_x
  )
}

#' Generate Predictions from Fitted Model
#'
#' @description
#' Generate smooth prediction curves from fitted models over a specified range
#' of x values. Supports both single models and grouped models with shared parameters.
#'
#' @param fit A fitted `nls` object from `fit_model()`.
#' @param lower_x The minimum x value for predictions.
#' @param upper_x The maximum x value for predictions.
#' @param group Optional character vector of group values for shared parameter models.
#' @param num_points Number of prediction points between lower_x and upper_x (default: 1000).
#'
#' @return A tibble containing:
#' \describe{
#'   \item{x}{Prediction x values}
#'   \item{y}{Predicted y values}
#'   \item{group}{Group identifier (if applicable)}
#' }
#'
#' @family model_analysis
#' @importFrom dplyr bind_rows
#' @importFrom purrr map
#' @importFrom rlang is_null
#' @importFrom tibble tibble
#' @importFrom stats predict
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate smooth prediction curve
#' predictions <- predict_data(
#'   fit = growth_fit,
#'   lower_x = 0,
#'   upper_x = 100,
#'   num_points = 500
#' )
#'
#' # For grouped model
#' predictions <- predict_data(
#'   fit = grouped_fit,
#'   lower_x = 0,
#'   upper_x = 100,
#'   group = c("cell_line_A", "cell_line_B")
#' )
#' }
predict_data <- function(
  fit,
  lower_x,
  upper_x,
  group = NULL,
  num_points = 1e3
) {
  x_values <- seq(lower_x, upper_x, length.out = num_points)

  if (!rlang::is_null(group)) {
    purrr::map(group, function(i) {
      tibble(
        x = x_values,
        y = predict(fit, newdata = tibble(x = x_values, group = i)),
        group = i
      )
    }) |>
      dplyr::bind_rows()
  } else {
    tibble::tibble(
      x = x_values,
      y = predict(fit, newdata = tibble(x = x_values))
    )
  }
}

#' Fit Models to Concentration Series Data
#'
#' @description
#' Apply model fitting across multiple concentration levels in concentration-response
#' experiments. Handles parallel processing and error recovery for failed fits.
#'
#' @param data A data frame containing concentration-response data with a `concentration` column.
#' @param x_var Character string specifying the time variable column name.
#' @param y_var Character string specifying the response variable column name.
#' @param model Character string specifying the built-in model to use.
#' @param model_func Optional custom model function.
#' @param start_func Optional starting values function.
#' @param start_values Optional named list of starting values.
#' @param huber Logical for Huber robust regression (default: FALSE).
#' @param huber_opts Options for Huber regression.
#' @param rout Logical for ROUT outlier detection (default: FALSE).
#' @param rout_opts Options for ROUT method.
#' @param share_group Character string for grouping variable.
#' @param share_params Character vector of shared parameters.
#' @param lower_bounds Named list of lower bounds.
#' @param upper_bounds Named list of upper bounds.
#' @param return_func Logical for returning model function.
#'
#' @return A list where each element corresponds to a concentration level,
#'   containing the fit results and processed data.
#'
#' @details
#' This function automatically handles:
#' \itemize{
#'   \item Iteration over unique concentration levels
#'   \item Progress tracking with progress bars
#'   \item Error handling for failed fits
#'   \item Exclusion of completely flagged datasets
#'   \item Warning suppression for cleaner output
#' }
#'
#' @family model_fitting
#' @importFrom purrr map
#' @importFrom dplyr filter
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit growth curves across concentration series
#' conc_fits <- fit_model_conc(
#'   data = concentration_data,
#'   x_var = "time_hours",
#'   y_var = "cell_area",
#'   model = "logistic_growth"
#' )
#' }
fit_model_conc <- function(
  data,
  x_var,
  y_var,
  model,
  model_func = NULL,
  start_func = NULL,
  start_values = NULL,
  huber = FALSE,
  huber_opts = NULL,
  rout = FALSE,
  rout_opts = NULL,
  share_group = NULL,
  share_params = NULL,
  lower_bounds = NULL,
  upper_bounds = NULL,
  return_func = FALSE
) {
  purrr::map(.progress = TRUE, unique(data$concentration), function(a) {
    data_filt <- dplyr::filter(data, concentration == !!a)

    if ("exclude" %in% colnames(data_filt)) {
      if (all(data_filt$exclude)) {
        return(list(
          fit = NULL,
          data = dplyr::filter(data, concentration == !!a)
        ))
      }
    }

    res <- tryCatch(
      withCallingHandlers(
        fit_model(
          dplyr::filter(data, concentration == !!a),
          x_var,
          y_var,
          model,
          model_func = NULL,
          start_func = NULL,
          start_values = NULL,
          huber_opts = NULL,
          rout = TRUE,
          rout_opts = NULL,
          share_group = NULL,
          share_params = NULL,
          lower_bounds = NULL,
          upper_bounds = NULL,
          return_func = FALSE
        ),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )

    if (inherits(res, "try-error")) {
      res <- list(
        fit = NULL,
        data = dplyr::filter(data, concentration == !!a)
      )
    }
    res
  })
}
