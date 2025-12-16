#' Calculate Growth Rate Inhibition Metrics
#' @description
#' Computes GR (growth rate inhibition) values and optional NDR (normalised
#'    drug response) values from time-course data. Supports input as either raw
#'    data or pre-fitted model objects.
#' @param data A data frame containing time-course measurements. Required if
#'   \code{model_list} is NULL.
#' @param cap Logical; if TRUE (default), GR and NDR values are capped at 1.
#' @param return_all_columns Logical; if TRUE, returns all intermediate columns
#'   used in calculations. If FALSE (default), returns only original columns
#'   plus \code{gr} (and \code{ndr} if applicable).
#' @param model_list A named list of model results, typically output from
#'   \code{\link{interpolate_curve_concentration}}. Each element should contain
#'   \code{selected} (the fitted model) and \code{data} (the subset of data).
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
#' @return
#' A data frame containing the original data plus:
#'  \itemize{
#'    \item `gr`: Growth rate inhibition values, ranging from -1 (complete cell
#'      death) through 0 (cytostasis) to 1 (no effect).
#'    \item `ndr`: Normalised drug response values (only if
#'      \code{positive_control_name} is provided), ranging from -1 to 1.
#'  }
#' If \code{return_all_columns = TRUE}, additional columns include:
#'  \code{time_zero_values}, \code{negative_control_zero}, \code{normalised},
#'  \code{negative_control_norm}, and (if applicable) \code{positive_control_zero}
#'   and \code{positive_control_norm}.
#' @details
#' GR values are calculated using the formula from Hafner et al. (2016):
#' \deqn{GR = 2^{log_2(x(c)/x_0) / log_2(x_{ctrl}/x_0)} - 1}
#' where \eqn{x(c)} is the measured value at concentration \eqn{c}, \eqn{x_0}
#' is the value at time zero, and \eqn{x_{ctrl}} is the negative control value.
#' NDR values provide an additional normalisation against a positive control,
#' useful when a maximal effect reference is available.
#' Groups marked with \code{exclude = TRUE} in the input data will have their
#' GR and NDR values set to -1.
#' @references
#' Hafner, M., Niepel, M., Chung, M., & Sorger, P. K. (2016). Growth rate
#' inhibition metrics correct for confounders in measuring sensitivity to
#' cancer drugs. \emph{Nature Methods}, 13(6), 521-527.
#' @export
calc_inhibition_metrics <- function(
  data = NA,
  model_list = NA,
  x_var,
  y_var,
  treatment_column,
  concentration_column,
  negative_control_name,
  positive_control_name = NA,
  cap = TRUE,
  return_all_columns = FALSE
) {
  if (!all(is.na(model_list))) {
    obj_list <- lapply(model_list, function(x) {
      x$fit$obj
    })
    data_list <- lapply(model_list, function(x) {
      x$data
    })
    all_cols <- Reduce(union, lapply(data_list, names))
    data_list <- lapply(data_list, function(d) {
      d[all_cols] <- lapply(all_cols, function(x) d[[x]])
      d
    })
    data <- do.call(rbind, data_list)
    rownames(data) <- NULL

    # If some scans are missing
    data_all <- do.call(rbind, data_list)

    predictions <- mapply(
      obj_list,
      data_list,
      SIMPLIFY = FALSE,
      FUN = function(obj, data) {
        if ("exclude" %in% names(data) && all(data$exclude == TRUE)) {
          pred <- data.frame(
            x = seq(min(data_all$x), max(data_all$x), length.out = 1e3),
            y = 0
          )
          pred[[concentration_column]] <- unique(data[[concentration_column]])
          pred[[treatment_column]] <- unique(data[[treatment_column]])
          pred$exclude <- TRUE

          return(pred)
        }

        pred <- predict_data(
          obj = obj,
          lower_x = min(data_all$x),
          upper_x = max(data_all$x)
        )
        pred[[concentration_column]] <- unique(data[[concentration_column]])
        pred[[treatment_column]] <- unique(data[[treatment_column]])

        if ("exclude" %in% names(data) && all(data$exclude == FALSE)) {
          pred$exclude <- FALSE
        }
        return(pred)
      }
    )
    data <- do.call(rbind, predictions)
  } else {
    data <- prep_data(data, x_var, y_var)
  }
  rownames(data) <- NULL

  # Store original columns for later selection
  original_columns <- colnames(data)

  data <- data[order(data$x), ]

  data <- assign_condition(
    data,
    treatment_column,
    negative_control_name,
    positive_control_name
  )

  combinations <- lapply(
    split(data[[concentration_column]], data[[treatment_column]]),
    unique
  )
  combinations <- stack(combinations)
  names(combinations) <- c("concentration", "treatment")

  # Calculate time zero values for each treatment-concentration group
  group_id <- paste(
    data[[treatment_column]],
    data[[concentration_column]],
    sep = "_"
  )
  data$time_zero_values <- ave(
    data$y,
    group_id,
    FUN = function(y) rep(y[1], length(y))
  )

  # Extract control time zero value
  negative_control_mask <- data[[treatment_column]] == negative_control_name
  negative_control_zero <- unique(data$time_zero_values[negative_control_mask])

  if (length(negative_control_zero) != 1) {
    stop("Control group must have a unique time zero value")
  }

  # Add control time zero to all rows
  data$negative_control_zero <- negative_control_zero

  # Handle positive control if provided
  if (!is.na(positive_control_name)) {
    positive_control_mask <- data[[treatment_column]] == positive_control_name
    positive_control_zero <- unique(data$time_zero_values[
      positive_control_mask
    ])

    if (length(positive_control_zero) != 1) {
      stop("Positive control group must have a unique time zero value")
    }

    data$positive_control_zero <- positive_control_zero
  }

  # Calculate normalized values
  data$normalised <- data$y / data$time_zero_values

  # Get control normalized values for each time point
  unique_times <- unique(data$x)
  negative_control_at_time_norm <- numeric(length(
    unique_times
  ))
  for (i in seq_along(unique_times)) {
    t <- unique_times[i]
    negative_control_at_time <- data$normalised[
      negative_control_mask & data$x == t
    ]
    negative_control_at_time_norm[i] <- mean(
      negative_control_at_time,
      na.rm = TRUE
    )
  }

  # Create lookup and join
  negative_control_norm <- numeric(nrow(data))
  for (i in seq_along(unique_times)) {
    time_mask <- data$x == unique_times[i]
    negative_control_norm[time_mask] <- negative_control_at_time_norm[i]
  }
  data$negative_control_norm <- negative_control_norm

  # Handle positive control normalized values if provided
  if (!is.na(positive_control_name)) {
    positive_control_at_time_norm <- numeric(length(
      unique_times
    ))

    for (i in seq_along(unique_times)) {
      t <- unique_times[i]
      positive_control_at_time <- data$normalised[
        positive_control_mask & data$x == t
      ]
      positive_control_at_time_norm[i] <- mean(
        positive_control_at_time,
        na.rm = TRUE
      )
    }

    # Create lookup and join
    positive_control_norm <- numeric(nrow(data))
    for (i in seq_along(unique_times)) {
      time_mask <- data$x == unique_times[i]
      positive_control_norm[time_mask] <- positive_control_at_time_norm[i]
    }
    data$positive_control_norm <- positive_control_norm
  }

  all_treatment_names <- unique(data[[treatment_column]])

  # Remove last time point per group
  keep_rows <- logical(nrow(data))
  for (treat in all_treatment_names) {
    treat_rows <- which(data[[treatment_column]] == treat)
    if (length(treat_rows) > 1) {
      keep_rows[treat_rows[-length(treat_rows)]] <- TRUE
    }
  }
  data <- data[keep_rows, ]

  # Calculate GR values
  data$gr <- 2^(log2(data$normalised) / log2(data$negative_control_norm)) - 1

  # Apply cap to GR values if requested
  if (cap) {
    data$gr[data$gr > 1] <- 1
  }

  # Calculate NDR values if positive control provided
  if (!is.na(positive_control_name)) {
    data$ndr <- (1 -
      2^(log2(data$normalised) / log2(data$positive_control_norm))) /
      (1 -
        2^(log2(data$negative_control_norm) / log2(data$positive_control_norm)))

    # Clamp NDR values
    data$ndr[data$ndr > 1] <- 1
    data$ndr[data$ndr < -1] <- -1
  }

  rownames(data) <- NULL

  if ("exclude" %in% colnames(data)) {
    mask_exclude <- data$exclude
    data[mask_exclude, ]$gr <- -1
    if ("ndr" %in% colnames(data)) {
      data[mask_exclude, ]$ndr <- -1
    }
  }

  # Select columns to return
  if (return_all_columns) {
    return(data)
  } else {
    # Return original columns plus GR (and NDR if calculated)
    cols_to_select <- c(original_columns, "gr")
    if (!is.na(positive_control_name)) {
      cols_to_select <- c(cols_to_select, "ndr")
    }
    return(data[, cols_to_select])
  }
}
