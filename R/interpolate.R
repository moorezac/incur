#' Fit Multiple Models and Select Best by BIC
#' @description
#' Attempts to fit each specified model to the data and returns the best
#' fit according to BIC, along with all successful fits.
#' @param data A data frame with x and y columns.
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @param models Character vector of model names from `incur_models`.
#' @param outlier_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `huber`: Logical; use Huber robust regression? (default: FALSE).
#'     \item `huber_iter`: Maximum iterations for Huber (default: 100).
#'     \item `huber_k`: Huber tuning constant (default: 1.345).
#'     \item `huber_tol`: Convergence tolerance for Huber (default: 1e-6).
#'     \item `rout`: Logical; use ROUT outlier detection? (default: FALSE).
#'     \item `rout_q`: False discovery rate for ROUT (default: 1e-3).
#'     \item `rout_scale`: Scale estimator for ROUT, either "mad" or "quantile" (default: "mad").
#'   }
#' @return
#' A list containing:
#'  \itemize{
#'    \item `selected`: The best-fitting model object.
#'    \item `all_fits`: Named list of all successful fit objects.
#'    \item `chosen_model`: Character string naming the selected model.
#'    \item `bic_values`: Named numeric vector of BIC values.
#'  }
#'  Returns NA if no models fit successfully.
#' @keywords internal
select_best_model <- function(
  data,
  x_var,
  y_var,
  models,
  outlier_opts = NA
) {
  successful_fits <- list()

  for (model in models) {
    curve_opts <- list(model = model)

    # This feels a little jank
    fit_result <- tryCatch(
      {
        expr <- suppressMessages(suppressWarnings(fit_curve(
          data = data,
          x_var = x_var,
          y_var = y_var,
          curve_opts = curve_opts,
          outlier_opts = outlier_opts
        )))
      },
      error = function(e) {
        invisible()
      }
    )

    if (!inherits(fit_result, "try-error")) {
      successful_fits[[model]] <- fit_result
    }
  }

  if (length(successful_fits) == 0) {
    return(NA)
  }

  # Extract fit objects and calculate BIC
  obj_list <- lapply(successful_fits, function(x) {
    x$fit$obj
  })
  bic_vals <- sapply(obj_list, BIC)

  # Select best model
  chosen <- names(which.min(bic_vals))

  list(
    selected = successful_fits[[chosen]],
    all_fits = successful_fits,
    chosen_model = chosen,
    bic_values = bic_vals
  )
}

#' Fit and Select Best Models Across Concentrations
#' @description
#' Fits multiple growth models to time-course data for each treatment-concentration
#' combination, selects the best model by BIC, and generates a diagnostic plot
#' comparing all fitted models.
#' @param data A data frame containing time-course measurements.
#' @param loess Logical; if TRUE, fits LOESS curves instead of parametric models.
#'   Default is FALSE.
#' @param loess_span Span to use for LOESS.
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#' @param outlier_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `huber`: Logical; use Huber robust regression? (default: FALSE).
#'     \item `huber_iter`: Maximum iterations for Huber (default: 100).
#'     \item `huber_k`: Huber tuning constant (default: 1.345).
#'     \item `huber_tol`: Convergence tolerance for Huber (default: 1e-6).
#'     \item `rout`: Logical; use ROUT outlier detection? (default: FALSE).
#'     \item `rout_q`: False discovery rate for ROUT (default: 1e-3).
#'     \item `rout_scale`: Scale estimator for ROUT, either "mad" or "quantile" (default: "mad").
#'   }
#' @param models Character vector of model names from `incur_models`.
#' @return
#' A list containing:
#'  \itemize{
#'    \item `model_list`: A named list (one element per treatment-concentration
#'      combination) where each element contains:
#'      \itemize{
#'        \item `selected`: The best-fitting model object (or NA if excluded).
#'        \item `data`: The subset of data for that combination.
#'      }
#'    \item: A \code{ggplot2} object showing all fitted models for each
#'      concentration, with the selected model highlighted.
#'  }
#' @note Only a single treatment compound (plus controls) should be present in
#'   the data. The function will error if multiple non-control treatments are
#'   detected.
#' @seealso \code{\link{fit_curve}}, \code{\link{calc_inhibition_metrics}},
#'   \code{\link{predict_data}}
#' @importFrom ggplot2 ggplot aes geom_point geom_line guides guide_legend labs
#'   theme element_text facet_wrap
#' @export
interpolate_curve_concentration <- function(
  data,
  x_var,
  y_var,
  treatment_column,
  concentration_column,
  negative_control_name,
  positive_control_name = NA,
  loess = FALSE,
  loess_span = 0.3,
  outlier_opts = NA,
  models = names(incur_models)
) {
  data <- prep_data(data, x_var, y_var)

  data <- assign_condition(
    data,
    treatment_column,
    negative_control_name,
    positive_control_name
  )

  # This avoids a nested loop
  combinations <- lapply(
    split(data[[concentration_column]], data[[treatment_column]]),
    unique
  )
  combinations <- stack(combinations)
  names(combinations) <- c("concentration", "treatment")

  shared_opts <- NA
  outlier_opts_original <- outlier_opts

  if (loess) {
    chosen_model <- "loess"
  }

  model_list <- mapply(
    combinations$concentration,
    combinations$treatment,
    SIMPLIFY = FALSE,
    FUN = function(conc, treat) {
      data_filtered <- data[
        data[[concentration_column]] == conc &
          data[[treatment_column]] == treat,
      ]

      # Check if this concentration should be excluded
      if (
        "exclude" %in%
          names(data_filtered) &&
          any(data_filtered$exclude == TRUE)
      ) {
        return(list(
          selected = NA,
          data = data_filtered,
          prediction_data = NA
        ))
      }

      if (loess) {
        chosen_model <- "loess"
        obj <- loess(formula = y ~ x, data = data_filtered, span = loess_span)
        model_selection <- list()
        model_selection$all_fits$loess$fit$obj <- obj
        model_selection$all_fits$loess$data <- data_filtered
      } else {
        # Fit parametric models and select best
        model_selection <- select_best_model(
          data = data_filtered,
          x_var = x_var,
          y_var = y_var,
          models = models,
          outlier_opts = outlier_opts_original
        )
        chosen_model <- model_selection$chosen_model

        if (!length(model_selection$all_fits)) {
          stop(sprintf(
            "No models converged for %s at concentration %s",
            treat,
            conc
          ))
        }
      }

      obj_list <- lapply(model_selection$all_fits, function(x) {
        x$fit$obj
      })

      # Generate predictions
      predictions <- mapply(
        obj_list,
        names(obj_list),
        SIMPLIFY = FALSE,
        FUN = function(y, i) {
          pred <- predict_data(
            obj = y,
            lower_x = min(data_filtered$x),
            upper_x = max(data_filtered$x)
          )
          pred$model <- i
          pred[[concentration_column]] <- conc
          pred[[treatment_column]] <- treat
          pred$label <- ifelse(pred$model == chosen_model, chosen_model, NA)

          return(pred)
        }
      )
      prediction_data <- do.call(rbind, predictions)

      return(list(
        fit = model_selection$all_fits[[chosen_model]]$fit,
        data = model_selection$all_fits[[chosen_model]]$data,
        prediction_data = prediction_data
      ))
    }
  )
  names(model_list) <- paste(
    combinations$treatment,
    combinations$concentration,
    sep = "_"
  )

  # Extract predictions
  prediction_list <- lapply(model_list, function(x) {
    x$prediction_data
  })
  prediction_data <- do.call(rbind, prediction_list)
  rownames(prediction_data) <- NULL
  prediction_data <- prediction_data[!is.na(prediction_data$x), ]

  # This is easier for facet_wrap
  gg_data <- data
  gg_data$facet <- paste(
    gg_data[[treatment_column]],
    gg_data[[concentration_column]],
    sep = "_"
  )
  prediction_data$facet <- paste(
    prediction_data[[treatment_column]],
    prediction_data[[concentration_column]],
    sep = "_"
  )

  # Plot
  gg_all <- ggplot2::ggplot(
    mapping = ggplot2::aes(x = x, y = y, colour = model)
  ) +
    ggplot2::geom_point(
      data = gg_data,
      colour = "black",
      alpha = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = prediction_data,
      alpha = 0.5
    ) +
    ggplot2::geom_line(
      data = prediction_data[!is.na(prediction_data$label), ],
      linewidth = 2,
      alpha = 0.75,
      show.legend = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        title = "Model",
        override.aes = list(alpha = 1)
      )
    ) +
    ggplot2::labs(x = x_var, y = y_var) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 8)) +
    theme_incur() +
    ggplot2::facet_wrap(
      ~facet,
      scales = "free_y"
    )

  return(list(
    model_list = lapply(model_list, function(x) {
      x[!names(x) == "prediction_data"]
    }),
    plot = gg_all
  ))
}
