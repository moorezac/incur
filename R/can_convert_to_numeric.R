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
