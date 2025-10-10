#' Find Index of Closest Value in Vector
#'
#' @description
#' Locate the index position of the value in a numeric vector that is closest
#' to a specified target value.
#'
#' @param vector A numeric vector to search within.
#' @param target The target value to find the closest match for.
#'
#' @return An integer index indicating the position of the closest value in `vector`.
#'
#' @family utility_functions
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' values <- c(1.2, 3.7, 5.1, 8.9, 12.3)
#' closest_idx <- index_of_closest_value(values, 5.0)  # Returns 3
#' }
index_of_closest_value <- function(vector, target) {
  which(abs(vector - target) == min(abs(vector - target)))
}

#' Find X Value at Y Midpoint
#'
#' @description
#' For paired x and y vectors, find the x value where y reaches its midpoint
#' between minimum and maximum values. Useful for estimating EC50/IC50 starting values.
#'
#' @param x A numeric vector of x values.
#' @param y A numeric vector of y values (must be same length as x).
#'
#' @return The x value corresponding to the y midpoint.
#'
#' @family utility_functions
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' x_vals <- c(1, 2, 3, 4, 5)
#' y_vals <- c(10, 25, 50, 75, 90)
#' midpoint_x <- x_at_y_mid(x_vals, y_vals)  # Returns x where y ≈ 50
#' }
x_at_y_mid <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have equal lengths")
  }

  mid <- min(y) + (max(y) - min(y)) / 2
  index <- which.min(abs(y - mid))

  return(x[index])
}

#' Prepare Data with Standardized Column Names
#'
#' @description
#' Internal function to standardize data by creating `x` and `y` columns from
#' user-specified variable names, with data cleaning and type conversion.
#'
#' @param data Input data frame.
#' @param x_var Character string specifying the x variable column name.
#' @param y_var Character string specifying the y variable column name.
#'
#' @return A data frame with standardized `x` and `y` columns, missing values removed,
#'   and variables converted to numeric type.
#'
#' @family utility_functions
#' @importFrom dplyr mutate relocate
#' @importFrom tidyr drop_na
#' @importFrom rlang sym
#' @keywords internal
prepare_data <- function(data, x_var, y_var) {
  data |>
    dplyr::mutate(
      x = !!sym(x_var),
      y = !!sym(y_var)
    ) |>
    tidyr::drop_na(x, y) |>
    dplyr::mutate(x = as.numeric(x), y = as.numeric(y)) |>
    dplyr::relocate(x, y)
}

#' Find Columns Identical to Target Column
#'
#' @description
#' Identify columns in a data frame that contain identical values to a specified
#' target column. Useful for detecting redundant variables or validating data structure.
#'
#' @param data A data frame to examine.
#' @param target Character string specifying the target column name.
#' @param coerce_to_char Logical indicating whether to coerce columns to character
#'   for comparison (default: TRUE).
#'
#' @return Character vector of column names that are identical to the target column.
#'
#' @family utility_functions
#' @importFrom purrr map map_lgl
#' @keywords internal
find_identical_to_column <- function(data, target, coerce_to_char = TRUE) {
  # Data to match for
  data_target <- data[[target]]

  if (is.character(target) && coerce_to_char) {
    data <- data.frame(
      purrr::map(data, as.character),
      stringsAsFactors = FALSE
    )
  }

  identical_cols <- purrr::map_lgl(colnames(data), function(x) {
    if (x != target) {
      # Exclude the target column itself
      identical(data[[x]], data_target)
    } else {
      FALSE
    }
  })

  colnames(data)[identical_cols]
}

#' Convert Log Molar Concentrations to Human-Readable Strings
#'
#' @description
#' Convert log10 molar concentration values to formatted strings with appropriate
#' unit prefixes (M, mM, µM, nM, pM, fM) for display in plots and tables.
#'
#' @param log_m Numeric vector of log10 molar concentration values.
#' @param digits Number of decimal places for formatting (default: 2).
#'
#' @return Character vector of formatted concentration strings.
#'
#' @family utility_functions
#' @importFrom purrr map_chr
#' @importFrom stringr str_c
#' @export
#'
#' @examples
#' # Convert log concentrations to readable format
#' log_concs <- c(-9, -6, -3, 0)
#' formatted <- log_molar_to_string(log_concs)
#' # Returns: c("1.00 nM", "1.00 µM", "1.00 mM", "1.00 M")
log_molar_to_string <- function(log_m, digits = 2) {
  purrr::map_chr(log_m, function(x) {
    if (x == 0) {
      return("Untreated")
    }
    # Define prefixes and scale factors
    units <- c("M", "mM", "µM", "nM", "pM", "fM")
    scale <- c(0, -3, -6, -9, -12, -15)

    # Find the best unit
    idx <- min(which(x >= scale))

    # Scale concentration and format
    value <- 10^x / 10^scale[idx]

    sprintf(stringr::str_c("%.", digits, "f %s"), value, units[idx])
  })
}
