#' This the function used to fit
#'
#' @param .data A `data.frame` or `data.frame` extension (tibble) in long format.
#' @param .x_var An quoted/unquoted argument that refers to the the `x` value within `.data`.
#' @param .y_var An quoted/unquoted argument that refers to the the `y` value within `.data`.
#' @param .curve_func A function that describes a curve/model in terms of `x`. See below for examples.
#' @param .start_func A function with arguments `x` and `y` that returns a named list of starting values for all arguments/parameters within `.curv_func`. See below for examples.
#' @param .start_vals A named list of starting values for all arguments/parameters within `curv_func`. This cannot be used in conjunction with `.detect_outliers`.
#' @param .huber Perform iterative reweighted least squares non-linear regression using the Huber loss function.
#' @param .detect_outliers Boolean to indicate whether to filter outliers as described in Motulsky and Brown (2006). It is highly recommended that the `.huber` argument is set to true in order to start with "robust" regression.
#' @param .shared_group For nested fits (shared parameters), the quoted/unquoted argument for the column in `.data` that refers to individual groups to share fitted parameters.
#' @param .shared_params For nested fits (shared parameters), a character vector that indicates which parameter(s) in `.curve_func` are to be shared across groups.
#' @param .lower_bounds A named list that contains the lower bounds for specified parameters. All other lower bounds will be set at `-Inf`.
#' @param .upper_bounds A named list that contains the upper bounds for specified parameters. All other upper bounds will be set at `Inf`.
#' @param .return_func For nested fits (shared parameters), a boolean to indicate whether to return a modified `.curve_func` used in this process.
#' @param ... Other arguments to be passed to `minpack.lm::nlsLM`, such as `control` or `weights`.
#' @return A named list containing:
#'  \itemize{
#'    \item `fit`: the fitted `nlsModel` object 
#'    \item `data`: the original data used in the fit. If `.detect_outliers` is true then column added in the format `outlier_{y_var}` will be added to indicate which points are detected as outliers.
#'  }
#'  
#' @export
#' @examples
#' # exponential plateau
#' func <- function(x, y0, ym, k) ym - (ym - y0) * exp(-k * x)
#' func_start <- function(x, y) list(ym = max(y), y0 = min(y), k = 1)
#'
#' # fit the curve to the data
#' # share y0 and k across state and force y0 > 0
#' fit <- fit_model(
#'   .data = Puromycin,
#'   .x_var = conc,
#'   .y_var = rate,
#'   .curve_func = func,
#'   .start_func = func_start,
#'   .detect_outliers = TRUE,
#'   .lower_bounds = list(y0 = 0),
#'   .shared_group = state,
#'   .shared_params = c("y0", "k"),
#'   control = minpack.lm::nls.lm.control(maxiter = 1e3)
#' )
#' # create data from each group
#' predicted <- map(unique(Puromycin$state), function(i) {
#'   # more data points = smoother curve
#'   # this is easier than a geom_function approach
#'   x_vals <- seq(min(fit$data$x), max(fit$data$x), length.out = 1e4)
#'   # create the data and return
#'   tibble(
#'     x = x_vals,
#'     y = predict(res$fit, newdata = tibble(x = x_vals, group = i)),
#'     group = i
#'   )
#' })
#' # collate
#' predicted <- bind_rows(predicted)
#' # plot
#' ggplot(mapping = aes(x, y, colour = group)) +
#'   geom_point(data = fit$data) +
#'   geom_line(data = predicted)
fit_model <- function(
    .data,
    .x_var,
    .y_var,
    .curve_func,
    .start_func = NULL,
    .start_vals = NULL,
    .huber = FALSE,
    .detect_outliers = FALSE,
    .shared_group = NULL,
    .shared_params = NULL,
    .lower_bounds = NULL,
    .upper_bounds = NULL,
    .return_func = FALSE,
    ...) {
  # modify data
  .data <- mutate(.data, x = !!ensym(.x_var), y = !!ensym(.y_var))
  .data <- mutate(.data, x = as.numeric(x), y = as.numeric(y))
  .data <- mutate(.data, x = as.numeric(x), y = as.numeric(y))
  .data <- relocate(.data, x, y)

  .data <- drop_na(.data, x, y)
  
  # check input
  if (is_null(.start_func) & is_null(.start_vals)) {
    stop("both start_func and start_vals provided")
  }
  if (is_null(.start_func) & !is_null(.detect_outliers)) {
    stop("if detecting outliers a function needs to be provided to generate starting values")
  }

  # start values
  if (is_null(.start_vals)) {
    .start_vals <- .start_func(x = .data$x, y = .data$y)
  }

  # shared curve parameters
  if (!is_null(.shared_params)) {
    .data <- mutate(.data, group = !!ensym(.shared_group))
    .data <- mutate(.data, group = as.character(group))
    .data <- relocate(.data, group, .after = y)

    group_vec <- unique(.data$group)

    .shared_arguments <- make_shared_formals(
      .func = .curve_func,
      .group = group_vec,
      .params = .shared_params
    )
    .shared_body <- make_shared_body(
      .func = .curve_func,
      .group = group_vec,
      .params = .shared_params
    )
    .curve_func <- rlang::new_function(
      args = as.pairlist(.shared_arguments),
      body = .shared_body
    )

    .start_vals <- make_shared_start_vals(
      .data = .data,
      .start_func = .start_func,
      .group = group_vec,
      .params = .shared_params
    )
  }

  # formula
  function_formula <- make_formula(.curve_func)
  f <- .curve_func

  # create a list sans data to inject into minpack.lm::nlsLM
  final_arguments <- list(
    formula = function_formula,
    start = list_flatten(.start_vals)
  )

  # bounds
  if (!is_null(.lower_bounds) | !is_null(.upper_bounds)) {
    bounds <- make_bounds(.curve_func, .lower_bounds, .upper_bounds)
    final_arguments <- append(final_arguments, bounds)
  }
  # dots dots dots
  dots <- list(...)
  final_arguments <- append(final_arguments, dots)

  # final final
  final_arguments <- discard(final_arguments, is_null)

  # see help page on rlang::exec for background on this approach
  fit <- try({
    new_env <- rlang::env(dat = .data, f = .curve_func)
    eval(
      # expr(robustbase::nlrob(data = data, !!!final_arguments)),
      expr(minpack.lm::nlsLM(data = dat, !!!final_arguments)),
      envir = new_env
    )
  })

  if (inherits(fit, "try-error")) {
    stop("error in initial model fit")
  }

  if (.huber) {
    iter <- 0
    iter_max <- 100
    converged <- FALSE
    k <- 1.345
    tol <- 1e-6
    fit_old <- fit

    # loop
    while (iter < iter_max && !converged) {
      coef_old <- coef(fit_old)
      resids <- residuals(fit_old)
      # mad
      s <- median(abs(resids - median(resids))) * 1.4826
      u <- resids / (s + 1e-10)
      weights_vec <- case_when(
        abs(u) <= k ~ 1,
        .default = k / abs(u)
      )

      arguments_new <- append(
        final_arguments,
        list(weights = weights_vec)
      )

      fit_new <- try({
        new_env <- rlang::env(dat = .data, f = .curve_func)
        eval(
          expr(minpack.lm::nlsLM(data = dat, !!!arguments_new)),
          envir = new_env
        )
      })
      coef_new <- coef(fit_new)

      # check convergence
      coef_change <- max(abs(coef_new - coef_old) / (abs(coef_old) + 1e-6))
      converged <- coef_change < tol

      # update for next iteration
      fit_old <- fit_new
      iter <- iter + 1
    }
    # message(str_glue("huber achieved {tol} tolerance in {iter} iterations"))
    fit <- fit_new
  }

  # outlier detection
  if (!.detect_outliers) {
    # return first fit
    final_fit <- fit
  } else {
    # which points are outliers within the fit
    outlier_indices <- find_outlier_indices(fit, .scale_method = "mad")

    # name a new column based on y_var and add in
    outlier_column <- str_c("outlier", rlang::as_name(enquo(.y_var)), sep = "_")
    .data <- mutate(
      .data,
      !!outlier_column := map_vec(row_number(), function(x) {
        if_else(x %in% outlier_indices, TRUE, FALSE)
      })
    )
    .data <- relocate(.data, !!outlier_column, .after = y)

    # if we have identified no outliers we can skip to the end
    if (all(.data[[outlier_column]] == FALSE)) {
      final_fit <- fit
    } else {
      # this time we need to filter
      .start_vals_filtered <- .start_func(
        .data |> filter(!!sym(outlier_column) == FALSE) |> pull(x),
        .data |> filter(!!sym(outlier_column) == FALSE) |> pull(y)
      )

      # create a list sans data to inject into minpack.lm::nlsLM
      final_arguments <- list(
        formula = function_formula,
        start = list_flatten(.start_vals)
      )
      # bounds
      if (!is_null(.lower_bounds) | !is_null(.upper_bounds)) {
        bounds <- make_bounds(.curve_func, .lower_bounds, .upper_bounds)
        final_arguments <- append(final_arguments, bounds)
      }
      # dots dots dots
      dots <- list(...)
      final_arguments <- append(final_arguments, dots)

      # final
      final_arguments <- discard(final_arguments, is_null)

      fit_outlier_removed <- try({
        new_env <- rlang::env(data = filter(.data, !!ensym(outlier_column) == FALSE))
        eval(
          # expr(robustbase::nlrob(data = data, !!!final_arguments)),
          expr(minpack.lm::nlsLM(data = data, !!!final_arguments)),
          envir = new_env
        )
      })

      # if this didn't fit then return original
      if (inherits(fit_outlier_removed, "try-error")) {
        message("error in the model fit with outliers removed")
        message("returning the model fit with outliers included")
        final_fit <- fit
      } else {
        final_fit <- fit_outlier_removed
      }
    }
    # end outlier detection
  }

  # return
  if (.return_func) {
    list(fit = final_fit, data = .data, func = .curve_func)
  } else {
    list(fit = final_fit, data = .data)
  }
}
