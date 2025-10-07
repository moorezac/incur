#' Plot Fitted Models with Data Points
#'
#' @description
#' Create publication-quality plots showing fitted model curves overlaid on
#' experimental data points, with support for grouped models and outlier flagging.
#'
#' @param data A data frame containing the experimental data.
#' #' @param fit A fitted model object from `fit_model()`.
#' @param x_var Character string specifying the x-axis variable name.
#' @param y_var Character string specifying the y-axis variable name.
#' @param return_data Logical indicating whether to return plot data instead of
#'   the plot object (default: FALSE).
#'
#' @return A ggplot object (default) or a list containing original and predicted
#'   data for external plotting (if return_data = TRUE).
#'
#' @details
#' The function automatically detects and handles:
#' \itemize{
#'   \item Grouped models (shared parameters): adds color coding by group
#'   \item Outlier flags: uses different point shapes for flagged outliers
#'   \item Legend formatting: provides informative legends for groups and outliers
#'   \item Smooth curves: generates high-resolution prediction curves
#' }
#'
#' @family visualization
#' @importFrom dplyr mutate
#' @importFrom ggplot2 geom_line geom_point ggplot labs aes guides guide_legend
#' @importFrom rlang as_name enquo ensym
#' @importFrom stringr str_c
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot single model fit
#' model_plot <- plot_model(
#'   data = growth_data,
#'   x_var = "time_hours",
#'   y_var = "cell_count",
#'   fit = fitted_model
#' )
#' print(model_plot)
#'
#' # Return data for custom plotting
#' plot_data <- plot_model(
#'   data = growth_data,
#'   x_var = "time_hours",
#'   y_var = "cell_count",
#'   fit = fitted_model,
#'   return_data = TRUE
#' )
#' }
plot_model <- function(data, fit, x_var, y_var, return_data = FALSE) {
  data <- prepare_data(data, x_var, y_var)

  # Check for grouping
  formula_vars <- all.vars(formula(fit))
  if ("group" %in% formula_vars) {
    group_vec <- unique(data$group)
    predicted <- predict_data(fit, min(data$x), max(data$x), group = group_vec)
    use_colour <- TRUE

    original_group <- find_identical_to_column(
      data = data,
      target = "group",
      coerce_to_char = TRUE
    )
  } else {
    predicted <- predict_data(fit, min(data$x), max(data$x))
    # gg <- ggplot2::ggplot(mapping = aes(x, y))
    use_colour <- FALSE
  }

  # Check for outliers
  outlier_column <- stringr::str_c(
    "outlier",
    rlang::as_name(rlang::enquo(y_var)),
    sep = "_"
  )
  if (outlier_column %in% colnames(data)) {
    use_shape <- TRUE
  } else {
    use_shape <- FALSE
  }

  # Aesthetics
  map_point <- map_line <- ggplot2::aes(x = x, y = y)
  if (use_colour) {
    map_point$colour <- as.name("group")
    map_line$colour <- as.name("group")
  }
  if (use_shape) {
    map_point$shape <- as.name(outlier_column)
  }

  # If wanted to plot elsewhere
  if (return_data) {
    return(list(data = data, predicted = predicted))
  }

  # Plot
  ggplot2::ggplot() +
    ggplot2::geom_point(data = data, mapping = map_point, alpha = 0.25) +
    ggplot2::geom_line(data = predicted, mapping = map_line) +
    ggplot2::labs(x = rlang::enquo(x_var), y = rlang::enquo(y_var)) +
    {
      if ("group" %in% formula_vars) {
        ggplot2::guides(colour = ggplot2::guide_legend(original_group))
      }
    } +
    {
      if (outlier_column %in% colnames(data)) {
        ggplot2::guides(shape = ggplot2::guide_legend("Outlier"))
      }
    } +
    theme_incur()
}
