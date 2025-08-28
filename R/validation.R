#' Validate Data Frame Structure and Required Columns
#'
#' @description
#' Check that input data is a valid data frame and contains all required columns,
#' providing clear error messages for missing components.
#'
#' @param data A data frame to validate.
#' @param required_cols Character vector of required column names.
#'
#' @return Invisible. Throws error if validation fails.
#'
#' @family validation
#' @importFrom stringr str_flatten_comma
#' @keywords internal
validate_data_frame_columns <- function(data, required_cols) {
  checkmate::assert_data_frame(data)
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Missing columns in data: ", stringr::str_flatten_comma(missing))
  }
}

#' Validate Model Fitting Function Inputs
#'
#' @description
#' Comprehensive validation of inputs to `fit_model()` and related functions,
#' checking data structure, variable names, function validity, and parameter consistency.
#'
#' @param data Input data frame.
#' @param x_var Character string for x variable name.
#' @param y_var Character string for y variable name.
#' @param model_func Mathematical model function.
#' @param start_func Starting values function.
#' @param start_values Starting values list.
#' @param shared_group Shared parameter group variable.
#' @param shared_params Shared parameter names.
#'
#' @return Invisible. Throws informative errors if validation fails.
#'
#' @family validation
#' @importFrom checkmate assert_data_frame assert_character assert_function assert_list
#' @importFrom rlang is_null
#' @keywords internal
validate_inputs <- function(
  data,
  x_var,
  y_var,
  model_func,
  start_func,
  start_values,
  shared_group,
  shared_params
) {
  # Basic input validation
  checkmate::assert_data_frame(as.data.frame(data))
  checkmate::assert_character(x_var, len = 1)
  checkmate::assert_character(y_var, len = 1)
  checkmate::assert_function(model_func)

  # Check data columns exist
  validate_data_frame_columns(data, c(x_var, y_var))

  # Validate start function or values
  if (!rlang::is_null(start_func)) {
    checkmate::assert_function(start_func)
  }
  if (!rlang::is_null(start_values)) {
    checkmate::assert_list(start_values)
  }

  # Shared parameters validation
  if (!rlang::is_null(shared_group) || !rlang::is_null(shared_params)) {
    validate_shared_parameters(data, shared_group, shared_params, model_func)
  }
}

#' Validate Shared Parameter Configuration
#'
#' @description
#' Ensure that shared parameter setup is correctly configured, including
#' proper pairing of group and parameter specifications and parameter name validity.
#'
#' @param data Input data frame.
#' @param shared_group Grouping variable name.
#' @param shared_params Parameter names to share.
#' @param model_func Mathematical model function.
#'
#' @return Invisible. Throws error if configuration is invalid.
#'
#' @family validation
#' @importFrom checkmate assert_character
#' @importFrom rlang is_null
#' @importFrom stringr str_c str_flatten_comma
#' @keywords internal
validate_shared_parameters <- function(
  data,
  shared_group,
  shared_params,
  model_func
) {
  # Both must be provided together
  if (xor(rlang::is_null(shared_group), rlang::is_null(shared_params))) {
    stop("'shared_group' and 'shared_params' must be provided together")
  }

  if (!rlang::is_null(shared_group)) {
    checkmate::assert_character(shared_group, len = 1)
    if (!shared_group %in% colnames(data)) {
      stop(stringr::str_c("Not found in data columns: ", shared_group))
    }
  }

  if (!rlang::is_null(shared_params)) {
    checkmate::assert_character(shared_params)
    formal_args <- names(formals(model_func))
    invalid <- setdiff(shared_params, formal_args)
    if (length(invalid) > 0) {
      stop(
        "Not found in model_func parameters: ",
        stringr::str_flatten_comma(invalid)
      )
    }
  }
}

#' Check Character to Numeric Conversion Feasibility
#'
#' @description
#' Test whether a character vector can be successfully converted to numeric
#' values, handling common formatting issues like comma separators.
#'
#' @param x Character vector to test.
#'
#' @return Logical indicating whether conversion is possible.
#'
#' @family validation
#' @importFrom stringr str_trim str_replace_all
#' @keywords internal
can_convert_to_numeric <- function(x) {
  cleaned <- stringr::str_trim(stringr::str_replace_all(x, ",", ""))
  tryCatch(
    {
      as.numeric(cleaned)
      TRUE
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
}

#' Check Character to Integer Conversion Feasibility
#'
#' @description
#' Test whether a character vector represents integer values and can be
#' converted to integer type without loss of information.
#'
#' @param x Character vector to test.
#'
#' @return Logical indicating whether integer conversion is appropriate.
#'
#' @family validation
#' @importFrom stringr str_trim str_replace_all
#' @keywords internal
can_convert_to_integer <- function(x) {
  cleaned <- stringr::str_trim(stringr::str_replace_all(x, ",", ""))
  tryCatch(
    {
      converted <- as.numeric(cleaned)
      all(converted == as.integer(converted), na.rm = TRUE)
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )
}

#' Convert Character Columns to Appropriate Numeric Types
#'
#' @description
#' Automatically convert character columns to integer or numeric types
#' based on their content, with preference for integer conversion when possible.
#'
#' @param data A data frame with potentially convertible character columns.
#'
#' @return Data frame with character columns converted to appropriate numeric types.
#'
#' @details
#' The conversion hierarchy is:
#' \enumerate{
#'   \item Integer conversion (if all values represent whole numbers)
#'   \item Numeric conversion (if values represent decimal numbers)
#'   \item Character retention (if values cannot be converted)
#' }
#'
#' @family validation
#' @importFrom dplyr mutate across
#' @importFrom tidyselect where
#' @keywords internal
convert_character_columns <- function(data) {
  dplyr::mutate(
    data,
    dplyr::across(tidyselect::where(is.character), function(x) {
      # Order - integer, numeric, character
      if (can_convert_to_integer(x)) {
        as.integer(x)
      } else if (can_convert_to_numeric(x)) {
        as.numeric(x)
      } else {
        x
      }
    })
  )
}
