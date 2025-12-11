#' NA Coalescing Operator
#' @param a TODO: description.
#' @param b TODO: description.
#' @keywords internal
`%??%` <- function(a, b) {
  if (!all(is.na(a))) {
    a
  } else {
    b
  }
}


#' Prepare Data
#' @param data A dataframe.
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#' @description
#' Create `x` and `y` columns from specified variable names.
#' @keywords internal
prep_data <- function(data, x_var, y_var) {
  # If already done
  if (all(c("x", "y", "x_original", "y_original") %in% colnames(data))) {
    return(data)
  }
  
  # Handle existing columns
  cols_to_rename <- c("x", "y")[c("x", "y") %in% colnames(data)]
  if (length(cols_to_rename) > 0) {
    names(data)[names(data) %in% cols_to_rename] <- paste(
      cols_to_rename,
      "_original",
      sep = ""
    )
  }
  
  # Extract columns by name
  x <- data[[x_var]]
  y <- data[[y_var]]
  
  # Subset data and coerce to numeric
  keep <- !(is.na(x) | is.na(y))
  
  data <- data[keep, , drop = FALSE]
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])
  
  # Reorder columns to have x, y first
  result <- cbind(data.frame(x = x, y = y), data)
  rownames(result) <- NULL
  
  return(result)
}


#' Convert Log Molar Concentrations to Human-Readable Strings
#' @description
#' Convert log10 molar concentration values to strings with appropriate
#'    unit prefixes (M, mM, µM, nM, pM, fM) for display in plots and tables.
#' @param log_m Numeric vector of log10 molar concentration values.
#' @param digits Number of decimal places for formatting (default: 2).
#' @return 
#' Character vector of formatted concentration strings.
#' @examples
#' # Convert log concentrations to readable format
#' log_concs <- c(-9, -6, -3, 0)
#' formatted <- log_molar_to_string(log_concs)
#' # Returns: c("1.00 nM", "1.00 µM", "1.00 mM", "1.00 M")
#' @export
log_m_to_str <- function(log_m, digits = 2) {
  sapply(
    log_m,
    function(x) {
      if (is.na(x)) {
        return(NA_character_)
      }
      if (x == 0) {
        return("Negative Control")
      }
      units <- c("M", "mM", "µM", "nM", "pM", "fM")
      scale <- c(0, -3, -6, -9, -12, -15)
      idx <- min(which(x >= scale))
      value <- 10^x / 10^scale[idx]
      sprintf(paste0("%.", digits, "f %s"), value, units[idx])
    },
    USE.NAMES = FALSE
  )
}
