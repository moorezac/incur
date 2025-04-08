#' @title Find the index of a vector closest to a target value.
#' @description From a vector of `numeric`, find which the index of which value is closest to a desired value.
#' @param vector A numeric vector to search within.
#' @param target The target value to find the closest value to.
#' @return An index of `vector` that contains the closest value to `target`
#'
index_of_closest_value <- function(vector, target) {
  which(abs(vector - target) == min(abs(vector - target)))
}

#' @title Find the closest x value at the midpoint of y.
#' @description For two vectors x and y of equal length, this will find the x value where y is at it's mid point.
#' @param x A numeric vector
#' @param y A numeric vector
#' @return The closest value for `x` that represents the midpoint of `y`.
#'
x_at_y_mid <- function(x, y) {
  if (length(x) != length(y)) {
    stop("x and y must have equal lengths")
  }

  mid <- min(y) + (max(y) - min(y)) / 2
  index <- which.min(abs(y - mid))

  return(x[index])
}

can_convert_to_numeric <- function(x) {
  # remove commas and trim whitespace
  cleaned <- str_trim(str_replace_all(x, ",", ""))

  # check if all values can be converted to numeric
  result <- tryCatch(
    expr = {
      as.numeric(cleaned)
      TRUE
    },
    warning = function(w) FALSE,
    error = function(e) FALSE
  )

  return(result)
}

format_concentrations <- function(vec, signif = 2) {
  vec <- as.numeric(vec)
  map_vec(
    x = vec,
    .f = function(x) {
      x_molar <- 10^x
      # convert to µM
      x_micro_molar <- x_molar * 1e6

      if (x_micro_molar >= 1) {
        # for concentrations ≥ 1 µM, display in µM
        # return(sprintf("%.2f µM", x_micro_molar))
        if (signif == 0) {
          return(str_c(round(x_micro_molar), " \u03bcM"))
        } else {
          return(str_c(format(x_micro_molar, digits = signif, nsmall = signif), " \u03bcM"))
        }
      } else {
        # for concentrations < 1 µM, convert to nM
        # return(sprintf("%.2f nM", x_nano_molar))
        x_nano_molar <- x_micro_molar * 1000
        if (signif == 0) {
          return(str_c(round(x_nano_molar), " nM"))
        } else {
          return(str_c(format(x_nano_molar, digits = signif, nsmall = signif), " nM"))
        }
      }
    }
  )
}

#' Automatically converted character vectors to numerics across a dataframe/tibble.
#'
#' @param df A dataframe or dataframe-like object.
#' @return The original `df` with character columns converted to integer/numeric.
#' @importFrom stringr str_replace_all str_trim
#' @importFrom dplyr across mutate
#' @importFrom tidyselect where

auto_convert_numerics <- function(df) {
  # fun to check if a character vector can be converted to numeric
  can_convert_to_numeric <- function(x) {
    # remove commas and trim whitespace
    cleaned <- stringr::str_trim(stringr::str_replace_all(x, ",", ""))

    # check if all non-NA values can be converted to numeric
    tryCatch(
      {
        as.numeric(cleaned)
        TRUE
      },
      warning = function(w) FALSE,
      error = function(e) FALSE
    )
  }

  # function to check if a character vector can be converted to int
  can_convert_to_integer <- function(x) {
    # remove commas and trim whitespace
    cleaned <- stringr::str_trim(stringr::str_replace_all(x, ",", ""))

    # check if all non-NA values can be converted to int without loss
    tryCatch(
      {
        converted <- as.numeric(cleaned)
        all(converted == as.integer(converted), na.rm = TRUE)
      },
      warning = function(w) FALSE,
      error = function(e) FALSE
    )
  }

  # iterate through columns and convert
  df <- dplyr::mutate(df, dplyr::across(tidyselect::where(is.character), function(x) {
    # int first, numeric second, else kep
    if (can_convert_to_integer(x)) {
      as.integer(x)
      # as.integer(stringr::str_trim(stringr::str_replace_all(x, ",", "")))
    } else if (can_convert_to_numeric(x)) {
      as.numeric(x)
      # as.numeric(stringr::str_trim(stringr::str_replace_all(x, ",", "")))
    } else {
      x
    }
  }))

  return(df)
}

