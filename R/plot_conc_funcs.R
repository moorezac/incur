#' Plot Concentration-Response Model Fits
#'
#' @description
#' Create comprehensive visualization of fitted models across multiple concentration
#' levels, with automatic color coding and concentration labeling.
#'
#' @param model_list List of fitted model results from `fit_model_conc()`.
#' @param x_var Character string specifying the x-axis variable name.
#' @param y_var Character string specifying the y-axis variable name.
#' @param start_colour Starting color for the concentration color palette (default: "#0085ca").
#'
#' @return A ggplot object showing data points and fitted curves for each concentration,
#'   with properly formatted concentration labels and color-coded legend.
#'
#' @details
#' The function automatically:
#' \itemize{
#'   \item Generates concentration-appropriate color palettes
#'   \item Formats concentration labels with appropriate units (M, mM, µM, nM, etc.)
#'   \item Handles excluded data points by coloring them red
#'   \item Places vehicle/control conditions appropriately in the legend
#'   \item Creates smooth fitted curves overlaid on data points
#' }
#'
#' @family visualization
#' @importFrom dplyr bind_rows filter pull mutate case_when
#' @importFrom purrr map map_lgl
#' @importFrom ggplot2 ggplot aes geom_point geom_line scale_colour_manual scale_fill_manual labs guides theme
#' @importFrom scales hue_pal
#' @importFrom tinter tinter
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot concentration-response curves
#' conc_plot <- plot_model_conc(
#'   model_list = fitted_concentrations,
#'   x_var = "time_hours",
#'   y_var = "cell_viability",
#'   start_colour = "#2E86AB"
#' )
#' print(conc_plot)
#' }
plot_model_conc <- function(
  model_list,
  x_var,
  y_var,
  start_colour = "#0085ca"
) {
  data_list <- purrr::map(model_list, 2)
  data <- dplyr::bind_rows(data_list)
  data <- prepare_data(data, x_var, y_var)

  fit_list <- purrr::map(model_list, 1)
  fitted_list <- purrr::map2(data_list, fit_list, function(a, b) {
    # If there is no fit object
    if (rlang::is_null(b)) {
      return(NULL)
    }

    data <- prepare_data(a, x_var, y_var)
    fitted <- predict_data(
      fit = b,
      lower_x = min(data$x),
      upper_x = max(data$x)
    )
    fitted$concentration <- unique(data$concentration)

    fitted
  })
  fitted <- dplyr::bind_rows(fitted_list)

  # Concentrations
  uni_conc <- sort(unique(data$concentration))
  # Vehicle/0 will be at the wrong end
  uni_conc <- c(uni_conc[length(uni_conc)], uni_conc[-length(uni_conc)])
  num_conc <- length(uni_conc)

  # Labels
  labels <- log_molar_to_string(uni_conc)

  # Tints
  if (num_conc != 1) {
    steps <- floor(num_conc / 2) + 1
  } else {
    steps <- num_conc
  }
  tints <- tinter::tinter(
    x = start_colour,
    steps = steps,
    adjust = -0.1
  )[c(1:num_conc)]
  tints <- rev(tints)
  names(tints) <- uni_conc

  tints[names(tints) == 0] <- "#000000"

  if ("exclude" %in% colnames(data)) {
    keep <- dplyr::filter(data, exclude == FALSE) |>
      dplyr::pull(concentration) |>
      unique()
    drop <- dplyr::filter(data, exclude == TRUE) |>
      dplyr::pull(concentration) |>
      unique()

    # If there are still values left ie. don't exclude all replicates
    drop_idx <- purrr::map_lgl(drop, function(a) {
      data_filt <- dplyr::filter(data, concentration == a)
      all(data_filt$exclude)
    })
    drop <- drop[drop_idx]

    # Adjust colours
    tints[names(tints) %in% drop] <- "#ff0000"
  } else {
    # Just add this in
    data$exclude <- FALSE
  }

  data_filt <- dplyr::filter(data, exclude == FALSE)

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data_filt,
      mapping = ggplot2::aes(
        x,
        y,
        fill = factor(concentration),
        colour = factor(concentration)
      ),
      alpha = 0.25,
      key_glyph = draw_key_rect,
      show.legend = TRUE
    ) +
    ggplot2::geom_line(
      fitted,
      mapping = ggplot2::aes(
        x,
        y,
        colour = factor(concentration)
      ),
      key_glyph = draw_key_rect,
      show.legend = TRUE,
      linewidth = 1
    ) +
    ggplot2::scale_colour_manual(
      values = tints,
      limits = factor(uni_conc),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend("Concentration")
    ) +
    ggplot2::scale_fill_manual(
      values = tints,
      limits = factor(uni_conc),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend("Concentration")
    ) +
    ggplot2::labs(x = x_var, y = y_var) +
    ggplot2::guides(colour = NULL) +
    theme_incur()
}

#' Plot GR Values Over Time for Concentration Series
#'
#' @description
#' Generate time-course plots of Growth Rate (GR) values calculated from
#' concentration-response model fits, showing growth inhibition dynamics.
#'
#' @param model_list List of fitted model results from `fit_model_conc()`.
#' @param x_var Character string specifying the time variable name.
#' @param y_var Character string specifying the response variable name.
#' @param start_colour Starting color for the concentration palette (default: "#0085ca").
#'
#' @return A list containing:
#' \describe{
#'   \item{data}{Data frame with calculated GR values}
#'   \item{plot}{ggplot object showing GR vs time for each concentration}
#' }
#'
#' @details
#' GR values provide normalized growth inhibition metrics:
#' \itemize{
#'   \item GR = 0: no growth inhibition
#'   \item GR = 0.5: partial growth inhibition
#'   \item GR = 1: complete growth inhibition
#'   \item GR = -1: cell death (for excluded concentrations)
#' }
#'
#' @family visualization
#' @importFrom dplyr bind_rows filter
#' @importFrom purrr map
#' @importFrom ggplot2 geom_hline
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot GR over time
#' gr_time_result <- plot_model_conc_gr_time(
#'   model_list = fitted_concentrations,
#'   x_var = "time_hours",
#'   y_var = "cell_area"
#' )
#'
#' # View the plot
#' gr_time_result$plot
#'
#' # Access the data
#' gr_data <- gr_time_result$data
#' }
plot_model_conc_gr_time <- function(
  model_list,
  x_var,
  y_var,
  start_colour = "#0085ca"
) {
  data_list <- purrr::map(model_list, 2)
  data <- dplyr::bind_rows(data_list)
  data <- prepare_data(data, x_var, y_var)

  fit_list <- purrr::map(model_list, 1)
  fitted_list <- purrr::map2(data_list, fit_list, function(a, b) {
    # If there is no fit object
    if (rlang::is_null(b)) {
      return(NULL)
    }

    data <- prepare_data(a, x_var, y_var)
    fitted <- predict_data(
      fit = b,
      lower_x = min(data$x),
      upper_x = max(data$x)
    )
    fitted$concentration <- unique(data$concentration)

    fitted
  })
  fitted <- dplyr::bind_rows(fitted_list)

  fitted <- calculate_gr(
    data = fitted,
    x_var = "x",
    y_var = "y",
    share_group = "concentration",
    control_name = 0
  )

  # Concentrations
  uni_conc <- sort(unique(data$concentration))
  # Vehicle/0 will be at the wrong end
  uni_conc <- c(uni_conc[length(uni_conc)], uni_conc[-length(uni_conc)])
  num_conc <- length(uni_conc)

  # Labels
  labels <- log_molar_to_string(uni_conc)

  # Tints
  if (num_conc != 1) {
    steps <- floor(num_conc / 2) + 1
  } else {
    steps <- num_conc
  }
  tints <- tinter::tinter(
    x = start_colour,
    steps = steps,
    adjust = -0.1
  )[c(1:num_conc)]
  tints <- rev(tints)
  names(tints) <- uni_conc

  tints[names(tints) == 0] <- "#000000"

  if ("exclude" %in% colnames(data)) {
    keep <- dplyr::filter(data, exclude == FALSE) |>
      dplyr::pull(concentration) |>
      unique()
    drop <- dplyr::filter(data, exclude == TRUE) |>
      dplyr::pull(concentration) |>
      unique()

    # If there are still values left ie. don't exclude all replicates
    drop_idx <- purrr::map_lgl(drop, function(a) {
      data_filt <- dplyr::filter(data, concentration == a)
      all(data_filt$exclude)
    })
    drop <- drop[drop_idx]

    # Adjust colours
    tints[names(tints) %in% drop] <- "#ff0000"

    # Add in 0 values
    x_vals <- unique(fitted$x)
    zero_vals <- purrr::map(drop, function(a) {
      tibble(x = x_vals, concentration = a, gr = -1)
    })
    fitted <- dplyr::bind_rows(fitted, zero_vals)
  }

  # The first GR values is NaN
  fitted <- dplyr::filter(fitted, is.finite(gr))

  gg <- ggplot2::ggplot() +
    ggplot2::geom_line(
      fitted,
      mapping = ggplot2::aes(
        x,
        gr,
        colour = factor(concentration)
      ),
      key_glyph = draw_key_rect,
      show.legend = TRUE,
      linewidth = 1
    ) +
    ggplot2::scale_colour_manual(
      values = tints,
      limits = factor(uni_conc),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend("Concentration")
    ) +
    ggplot2::scale_fill_manual(
      values = tints,
      limits = factor(uni_conc),
      breaks = names(tints),
      labels = labels,
      drop = FALSE,
      guide = ggplot2::guide_legend("Concentration")
    ) +
    ggplot2::geom_hline(
      yintercept = 0
    ) +
    ggplot2::guides(colour = NULL) +
    ggplot2::labs(x = x_var, y = "GR") +
    theme_incur()

  list(data = fitted, plot = gg)
}

#' Plot GR Dose-Response Relationship with Model Fit
#'
#' @description
#' Create dose-response plots using GR values, with optional model fitting
#' to estimate GR50 (concentration for 50% growth inhibition).
#'
#' @param model_list List of fitted model results from `fit_model_conc()`.
#' @param x_var Character string for the time variable name.
#' @param y_var Character string for the response variable name.
#' @param model Character string specifying model for GR dose-response fitting
#'   (default: "five_param_sigmoid_log").
#' @param start_colour Color for the fitted dose-response curve (default: "#0085ca").
#'
#' @return A ggplot object showing GR values vs concentration with fitted dose-response
#'   curve and GR50 annotation if successfully calculated.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Calculates GR values for each concentration and time point
#'   \item Fits a dose-response model to GR vs concentration data
#'   \item Estimates GR50 (concentration giving 50% growth inhibition)
#'   \item Creates publication-ready plots with log-scale concentration axis
#'   \item Colors points by time to show temporal progression
#' }
#'
#' @family visualization
#' @importFrom dplyr bind_rows filter rename
#' @importFrom ggplot2 scale_colour_viridis_c scale_x_continuous geom_hline
#' @importFrom geomtextpath geom_textvline
#' @export
#'
#' @examples
#' \dontrun{
#' # Create GR dose-response plot with GR50 estimation
#' gr_dose_plot <- plot_model_conc_gr_dose(
#'   model_list = fitted_concentrations,
#'   x_var = "time_hours",
#'   y_var = "cell_viability",
#'   model = "four_param_sigmoid_log"
#' )
#' print(gr_dose_plot)
#' }
plot_model_conc_gr_dose <- function(
  model_list,
  x_var,
  y_var,
  model = "five_param_sigmoid_log",
  start_colour = "#0085ca"
) {
  data_list <- purrr::map(model_list, 2)
  data <- dplyr::bind_rows(data_list)
  data <- prepare_data(data, x_var, y_var)

  fit_list <- purrr::map(model_list, 1)
  fitted_list <- purrr::map2(data_list, fit_list, function(a, b) {
    # If there is no fit object
    if (rlang::is_null(b)) {
      return(NULL)
    }

    data <- prepare_data(a, x_var, y_var)
    fitted <- predict_data(
      fit = b,
      lower_x = min(data$x),
      upper_x = max(data$x)
    )
    fitted$concentration <- unique(data$concentration)

    fitted
  })
  fitted <- dplyr::bind_rows(fitted_list)

  fitted <- calculate_gr(
    data = fitted,
    x_var = "x",
    y_var = "y",
    share_group = "concentration",
    control_name = 0
  )

  if ("exclude" %in% colnames(data)) {
    keep <- dplyr::filter(data, exclude == FALSE) |>
      dplyr::pull(concentration) |>
      unique()
    drop <- dplyr::filter(data, exclude == TRUE) |>
      dplyr::pull(concentration) |>
      unique()

    # If there are still values left ie. don't exclude all replicates
    drop_idx <- purrr::map_lgl(drop, function(a) {
      data_filt <- dplyr::filter(data, concentration == a)
      all(data_filt$exclude)
    })
    drop <- drop[drop_idx]

    # Add in 0 values
    x_vals <- unique(fitted$x)
    zero_vals <- purrr::map(drop, function(a) {
      tibble(x = x_vals, gr = -1, concentration = a)
    })
    fitted <- dplyr::bind_rows(fitted, zero_vals)
  }

  fitted <- dplyr::filter(fitted, is.finite(gr))

  fitted_gr_model <- try(
    {
      fit_model(
        data = dplyr::filter(fitted, !concentration == 0),
        x_var = "concentration",
        y_var = "gr",
        model = model,
        upper_bounds = list(top = 1, bottom = -1),
        lower_bounds = list(top = 1, bottom = -1)
      )
    },
    silent = TRUE
  )
  # If I can't fit a model
  if (inherits(fitted_gr_model, "try-error")) {
    GR50 <- NULL
    fitted_gr <- NULL
  } else {
    # Get the predictions
    fitted_gr <- predict_data(
      fit = fitted_gr_model$fit,
      lower_x = min(fitted_gr_model$data$x),
      upper_x = max(fitted_gr_model$data$x)
    )

    # If I can't find the GR50
    GR50 <- try({
      find_x_for_y(
        fit = fitted_gr_model$fit,
        x_values = fitted$concentration,
        target = 0.5
      )
    })
    if (inherits(GR50, "try-error")) {
      GR50 <- NULL
    }
  }

  # For ease in plots
  fitted <- dplyr::rename(
    fitted,
    x_val = x,
    y_val = y,
    x = concentration,
    y = gr
  )
  fitted_filt <- dplyr::filter(fitted, !x == 0)

  # Plot
  ggplot2::ggplot() +
    ggplot2::geom_point(
      fitted_filt,
      mapping = ggplot2::aes(
        x,
        y,
        colour = x_val
      ),
      show.legend = TRUE
    ) +
    ggplot2::scale_colour_viridis_c() +
    {
      if (!rlang::is_null(fitted_gr)) {
        ggplot2::geom_line(
          fitted_gr,
          mapping = ggplot2::aes(
            x,
            y
          ),
          show.legend = TRUE,
          linewidth = 1,
          colour = start_colour
        )
      }
    } +
    {
      if (!rlang::is_null(GR50)) {
        geomtextpath::geom_textvline(
          xintercept = GR50$root,
          label = log_molar_to_string(GR50$root),
          colour = start_colour,
          linewidth = 0.5,
          hjust = 0.1,
          size = 10 / .pt
        )
      }
    } +
    ggplot2::guides(colour = ggplot2::guide_colourbar(x_var)) +
    ggplot2::scale_x_continuous(
      labels = \(a) log_molar_to_string(a, digits = 0),
      limits = c(
        min(fitted_filt$x),
        max(fitted_filt$x)
      ),
      minor_breaks = log10(
        rep(1:9, 21) * (10^rep(-10:10, each = 9))
      )
    ) +
    ggplot2::geom_hline(
      yintercept = 0.5,
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = 0
    ) +
    ggplot2::labs(x = "Concentration", y = "GR") +
    theme_incur()
}
