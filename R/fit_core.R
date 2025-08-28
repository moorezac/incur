#' Attempt Model Fitting with Error Handling
#'
#' @description
#' Safely attempt to fit a nonlinear model using `minpack.lm::nlsLM()` with
#' proper error handling and environment management.
#'
#' @param data A data frame containing the data to fit, with columns `x` and `y`.
#' @param model_func A function defining the mathematical model.
#' @param arguments A named list of arguments to pass to `minpack.lm::nlsLM()`.
#'
#' @return Either a fitted `nls` object if successful, or an object of class
#'   `try-error` if fitting fails.
#'
#' @family model_fitting
#' @importFrom rlang env expr
#' @importFrom minpack.lm nlsLM
#' @keywords internal
try_fit <- function(data, model_func, arguments) {
  f <- model_func
  try(
    {
      new_env <- rlang::env(data = data)
      eval(
        expr(minpack.lm::nlsLM(data = data, !!!arguments)),
        envir = new_env
      )
    },
    silent = TRUE
  )
}

#' Execute Complete Model Fitting Pipeline
#'
#' @description
#' Internal function that executes the full model fitting process including
#' initial fitting, optional robust regression, and outlier detection.
#'
#' @param data A data frame with standardized `x` and `y` columns.
#' @param x_var Character string specifying the original x variable name.
#' @param y_var Character string specifying the original y variable name.
#' @param model_func A function defining the mathematical model.
#' @param start_func A function to generate starting parameter values.
#' @param start_values Optional named list of starting values (overrides start_func).
#' @param huber Logical indicating whether to apply Huber robust regression.
#' @param huber_opts List of options for Huber regression (iter_max, k, tol).
#' @param rout Logical indicating whether to apply ROUT outlier detection.
#' @param rout_opts List of options for ROUT method (q, scale_method).
#' @param share_group Character string for shared parameter grouping variable.
#' @param share_params Character vector of parameters to share across groups.
#' @param upper_bounds Named list of upper parameter bounds.
#' @param lower_bounds Named list of lower parameter bounds.
#' @param final_arguments Prepared arguments list for nlsLM.
#'
#' @return A list containing the fitted model object and processed data.
#'
#' @family model_fitting
#' @keywords internal
execute_fit <- function(
  data,
  x_var,
  y_var,
  model_func,
  start_func,
  start_values,
  huber,
  huber_opts,
  rout,
  rout_opts,
  share_group,
  share_params,
  upper_bounds,
  lower_bounds,
  final_arguments
) {
  # Initial fit
  fit <- try_fit(data, model_func, final_arguments)
  if (inherits(fit, "try-error")) {
    stop("Error in initial model fit: ", attr(fit, "condition")$message)
  }

  # Apply Huber IRWLS
  if (huber) {
    fit <- huber_irwls(fit, final_arguments, data, model_func, huber_opts)
  }

  # rout
  if (rout) {
    if (rlang::is_null(huber)) {
      huber <- TRUE
    }

    detect_outliers_and_refit(
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
    )
  } else {
    list(fit = fit, data = data)
  }
}

#' Fit Mathematical Models to Biological Data
#'
#' @description
#' Comprehensive function for fitting nonlinear mathematical models to biological
#' data with support for built-in models, custom functions, robust regression,
#' outlier detection, and shared parameters across groups.
#'
#' @param data A data frame or tibble containing the experimental data.
#' @param x_var Character string specifying the independent variable column name.
#' @param y_var Character string specifying the dependent variable column name.
#' @param model Character string specifying a built-in model name from `incur_models`.
#' @param model_func Optional custom function defining the mathematical model.
#' @param start_func Optional function to generate starting parameter values.
#' @param start_values Optional named list of starting parameter values.
#' @param huber Logical indicating whether to apply Huber robust regression (default: FALSE).
#' @param huber_opts List of Huber regression options (iter_max, k, tol).
#' @param rout Logical indicating whether to apply ROUT outlier detection (default: FALSE).
#' @param rout_opts List of ROUT options (q, scale_method).
#' @param share_group Character string specifying grouping variable for shared parameters.
#' @param share_params Character vector of parameter names to share across groups.
#' @param lower_bounds Named list of lower parameter bounds.
#' @param upper_bounds Named list of upper parameter bounds.
#' @param return_func Logical indicating whether to return the model function (default: FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{fit}{Fitted nls model object}
#'   \item{data}{Processed data with any outlier flags}
#'   \item{func}{Model function (if return_func = TRUE)}
#' }
#'
#' @details
#' This function provides a comprehensive interface for fitting mathematical models
#' to biological data. It supports:
#' \itemize{
#'   \item Built-in models via the `model` parameter
#'   \item Custom models via `model_func` and `start_func`
#'   \item Robust regression to reduce outlier influence
#'   \item Automatic outlier detection and removal
#'   \item Shared parameters across experimental groups
#'   \item Parameter bounds for constrained optimization
#' }
#'
#' @family model_fitting
#' @importFrom dplyr mutate relocate
#' @importFrom rlang is_null sym
#' @importFrom stringr str_glue
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit exponential growth model to cell growth data
#' result <- fit_model(
#'   data = growth_data,
#'   x_var = "time_hours",
#'   y_var = "cell_count",
#'   model = "exponential_growth"
#' )
#'
#' # Fit dose-response curve with outlier detection
#' result <- fit_model(
#'   data = dose_response,
#'   x_var = "concentration",
#'   y_var = "viability",
#'   model = "four_param_sigmoid_log",
#'   rout = TRUE
#' )
#'
#' # Fit with shared parameters across cell lines
#' result <- fit_model(
#'   data = multi_line_data,
#'   x_var = "time",
#'   y_var = "growth",
#'   model = "logistic_growth",
#'   share_group = "cell_line",
#'   share_params = c("k")  # Share growth rate parameter
#' )
#' }
fit_model <- function(
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
  # Check for built-in
  if (!rlang::is_null(model)) {
    if (!model %in% names(incur_models)) {
      stop(stringr::str_glue("Not a built-in model: {model}"))
    }
    model_func <- incur_models[[model]]$model_func
    start_func <- incur_models[[model]]$start_func
  }

  # Validate inputs
  validate_inputs(
    data,
    x_var,
    y_var,
    model_func,
    start_func,
    start_values,
    share_group,
    share_params
  )

  # Prepare data
  data <- prepare_data(data, x_var, y_var)

  # Prepare starting values
  if (rlang::is_null(start_values)) {
    start_values <- start_func(x = data$x, y = data$y)
  }

  # Apply shared parameters if needed
  if (!rlang::is_null(share_group) && !rlang::is_null(share_params)) {
    shared_parameters <- prepare_shared_parameters(
      data,
      share_params,
      share_group,
      model_func,
      start_func
    )
    data <- shared_parameters[[1]]
    model_func <- shared_parameters[[2]]
    start_values <- shared_parameters[[3]]
  }

  # Create final arguments for nlsLM
  final_arguments <- prepare_final_arguments(
    model_func,
    start_values,
    lower_bounds,
    upper_bounds
  )

  # Execute fit with post-processing
  result <- execute_fit(
    data,
    x_var,
    y_var,
    model_func,
    start_func,
    start_values,
    huber,
    huber_opts,
    rout,
    rout_opts,
    share_group,
    share_params,
    upper_bounds,
    lower_bounds,
    final_arguments
  )

  # Return results
  if (return_func) {
    list(fit = result[[1]], data = result[[2]], func = curve_func)
  } else {
    list(fit = result[[1]], data = result[[2]])
  }
}
