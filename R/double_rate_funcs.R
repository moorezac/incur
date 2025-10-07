#' Calculate Doubling Time from Passage Data
#'
#' @description
#' Calculate cell doubling time from known cell counts at two different time points,
#' typically from passage data or manual cell counting.
#'
#' @param start_date Character string or Date object representing the first measurement date in format `YYYY-MM-DD`.
#' @param end_date Character string or Date object representing the second measurement date in format `YYYY-MM-DD`.
#' @param start_n Numeric value for the number of cells at the first measurement.
#' @param end_n Numeric value for the number of cells at the second measurement.
#'
#' @return Numeric value representing the doubling time in hours.
#'
#' @details
#' Uses the exponential growth formula to calculate doubling time:
#' t_double = (t2 - t1) * ln(2) / ln(N2/N1)
#'
#' @family growth_analysis
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate doubling time from passage data
#' doubling_hours <- calc_double_rate_passage(
#'   start_date = "2024-01-01",
#'   end_date = "2024-01-03",
#'   start_n = 1e5,
#'   end_n = 4e5
#' )
#' }
calc_double_rate_passage <- function(start_date, end_date, start_n, end_n) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  days_between <- difftime(end_date, start_date) |> as.numeric()

  days_between * log(2) / (log(end_n / start_n)) * 24
}

#' Calculate Doubling Time from Image-Based Growth Data
#'
#' @description
#' Calculate cell doubling time from time-series imaging data by identifying
#' the exponential growth phase and fitting a linear model to log-transformed data.
#'
#' @param data A data frame containing time-series growth measurements.
#' @param x_var Character string specifying the time variable column name.
#' @param y_var Character string specifying the growth measurement column name (e.g., cell area, count).
#' @param cell_column Optional character string specifying a column for grouping multiple known cell populations.
#' @param dimension_factor Numeric multiplier to account for dimensional readout (default = 1 for linear measurements, use 1.5 for area measurements).
#'
#' @return A list containing:
#' \describe{
#'   \item{double_time}{Numeric value of calculated doubling time}
#'   \item{plot}{ggplot object showing linear and log-scale growth curves with fitted lines}
#' }
#'
#' @details
#' The function automatically identifies the exponential growth phase using
#' inflection point detection, then fits a linear model to log-transformed data
#' in this region. The dimension factor adjusts for different measurement types.
#'
#' @family growth_analysis
#' @importFrom dplyr filter group_by mutate ungroup case_when first pull bind_rows
#' @importFrom purrr map
#' @importFrom rlang is_null sym
#' @importFrom inflection ese
#' @importFrom ggplot2 ggplot aes geom_line geom_vline labs ggtitle theme guides guide_legend
#' @importFrom scales hue_pal
#' @importFrom stringr str_glue str_c
#' @importFrom geomtextpath geom_textabline
#' @importFrom patchwork wrap_plots
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate doubling time from imaging data
#' result <- calc_double_rate_data(
#'   data = growth_data,
#'   x_var = "time_hours",
#'   y_var = "total_area",
#'   cell_column = "cell_line"
#' )
#'
#' # View doubling time
#' result$double_time
#'
#' # Display plot
#' result$plot
#' }
calc_double_rate_data <- function(
  data,
  x_var,
  y_var,
  cell_column = NULL,
  dimension_factor = 1
) {
  data <- prepare_data(data, x_var, y_var)

  if (!rlang::is_null(cell_column)) {
    data <- dplyr::mutate(data, cells = !!rlang::sym(cell_column))

    # Check
    if (length(unique(data$cells)) == 1) {
      stop(stringr::str_glue("Only 1 unique value found for `{cell_column}`"))
    }
    if (!is.numeric(data$cells)) {
      stop("Provided `cell_column` is not numeric")
    }

    data_filt <- purrr::map(unique(data$cells), function(a) {
      data_filt <- dplyr::filter(data, cells == a)
      # Find where section is linear, this is exponential growth
      inflection_points <- inflection::ese(data_filt$x, log(data_filt$y), 0)
      linear_indexes <- inflection_points[1]:inflection_points[2]

      data_filt[linear_indexes, ]
    })
    data_filt <- dplyr::bind_rows(data_filt)

    # Fit model with offset
    data_filt <- dplyr::group_by(data_filt, cells)
    data_filt <- dplyr::mutate(
      data_filt,
      n0 = dplyr::first(y),
      log_n0 = log(n0)
    )
    data_filt <- dplyr::ungroup(data_filt, cells)

    model <- lm(
      log(y) ~ x + offset(log_n0),
      data_filt
    )
    model_double <- lm(
      dimension_factor * log(y) ~ x + offset(log_n0),
      data_filt
    )
  } else {
    # Single sample
    # Find where section is linear, this is exponential growth
    inflection_points <- inflection::ese(data$x, log(data$y), 0)
    linear_indexes <- inflection_points[1]:inflection_points[2]

    data_filt <- data[linear_indexes, ]

    # Fit model on single
    model <- lm(log(y) ~ x, data_filt)
    model_double <- lm(
      dimension_factor * log(y) ~ x,
      data_filt
    )
  }

  # Double calculation
  double_time <- log(2) / coef(model_double)[["x"]]

  ## Plots
  # Add in linear section
  if (!rlang::is_null(cell_column)) {
    data <- purrr::map(unique(data_filt$cells), function(a) {
      data_filt_cell <- dplyr::filter(data_filt, cells == a)
      data_cell <- dplyr::filter(data, cells == a)
      dplyr::mutate(
        data_cell,
        linear = dplyr::case_when(
          x %in% data_filt_cell$x ~ TRUE,
          .default = FALSE
        )
      )
    })
    data <- dplyr::bind_rows(data)
  } else {
    data <- dplyr::mutate(
      data,
      linear = dplyr::case_when(
        x %in% data_filt$x ~ TRUE,
        .default = FALSE
      )
    )
  }

  # Colours
  if (!rlang::is_null(cell_column)) {
    colours <- scales::hue_pal()(length(unique(data$cells)))
    names(colours) <- unique(data$cells)
  }

  # Create aes()
  if (!rlang::is_null(cell_column)) {
    aes_lin <- ggplot2::aes(x, y, colour = factor(cells))
    aes_log <- ggplot2::aes(x, log(y), colour = factor(cells))
  } else {
    aes_lin <- ggplot2::aes(x, y)
    aes_log <- ggplot2::aes(x, log(y))
  }

  # Linear
  gg_lin <- ggplot2::ggplot(data, aes_lin) +
    ggplot2::geom_line() +
    ggplot2::labs(x = x_var, y = y_var) +
    {
      if (!rlang::is_null(cell_column)) {
        ggplot2::guides(colour = guide_legend(title = cell_column))
      }
    } +
    purrr::map(unique(data$cells), function(a) {
      if (!rlang::is_null(cell_column)) {
        data_filt_cell <- dplyr::filter(data_filt, cells == a)
        ggplot2::geom_vline(
          xintercept = c(
            min(data_filt_cell$x),
            max(data_filt_cell$x)
          ),
          colour = colours[names(colours) == a],
          linetype = "dashed"
        )
      } else {
        ggplot2::geom_vline(
          xintercept = c(
            min(data_filt$x),
            max(data_filt$x)
          ),
          linetype = "dashed"
        )
      }
    })

  # Logarthmic
  gg_log <- ggplot2::ggplot(data, aes_log) +
    ggplot2::geom_line() +
    ggplot2::labs(x = x_var, y = stringr::str_glue("log({y_var})")) +
    {
      if (!rlang::is_null(cell_column)) {
        ggplot2::guides(colour = guide_legend(title = cell_column))
      }
    } +
    purrr::map(unique(data$cells), function(a) {
      if (!rlang::is_null(cell_column)) {
        data_filt_cell <- dplyr::filter(data_filt, cells == a)
        ggplot2::geom_vline(
          xintercept = c(
            min(data_filt_cell$x),
            max(data_filt_cell$x)
          ),
          colour = colours[names(colours) == a],
          linetype = "dashed"
        )
      } else {
        ggplot2::geom_vline(
          xintercept = c(
            min(data_filt$x),
            max(data_filt$x)
          ),
          linetype = "dashed"
        )
      }
    }) +
    purrr::map(unique(data$cells), function(a) {
      if (!rlang::is_null(cell_column)) {
        data_filt_cell <- dplyr::filter(data_filt, cells == a)
        intercept <- data_filt_cell$log_n0[1] + coef(model)[1]
      } else {
        intercept <- coef(model)[[1]]
      }

      # Return
      geomtextpath::geom_textabline(
        intercept = intercept,
        slope = coef(model)[[2]],
        label = (stringr::str_c(
          "y = ",
          signif(coef(model)[[2]], digits = 3),
          "\u00d7",
          "x + ",
          signif(intercept, digits = 3)
        )),
        size = 6 / ggplot2::.pt,
        vjust = 1.5
      )
    })

  gg_final <- patchwork::wrap_plots(
    gg_lin +
      ggplot2::ggtitle(stringr::str_c(
        "Time to Double: ",
        signif(double_time, 3),
        " ",
        x_var
      )) +
      ggplot2::theme(aspect.ratio = 1) +
      theme_incur(),
    gg_log +
      ggplot2::theme(aspect.ratio = 1) +
      theme_incur(),
    guides = "collect",
    axes = "collect"
  )

  list(double_time = double_time, plot = gg_final)
}
