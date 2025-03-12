#' Automatically converted character vectors to numerics across a dataframe/tibble.
#'
#' @param .df A dataframe or dataframe-like object.
#' @return The original `.df` with character columns converted to integer/numeric.

auto_convert_numerics <- function(.df) {
  # fun to check if a character vector can be converted to numeric
  can_convert_to_numeric <- function(x) {
    # remove commas and trim whitespace
    cleaned <- str_trim(str_replace_all(x, ",", ""))

    # check if all non-NA values can be converted to numeric
    tryCatch(
      expr = {
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
    cleaned <- str_trim(str_replace_all(x, ",", ""))

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
  .df |>
    mutate(across(
      where(is.character),
      ~ {
        # try to convert to int first
        if (can_convert_to_integer(.x)) {
          as.integer(str_trim(str_replace_all(.x, ",", "")))
          # then try to convert to numeric
        } else if (can_convert_to_numeric(.x)) {
          as.numeric(str_trim(str_replace_all(.x, ",", "")))
        }
        # if conversion not possible, keep as is
        else {
          .x
        }
      }
    ))
}
