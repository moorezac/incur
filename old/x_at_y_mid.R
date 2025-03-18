#' @title Find the closest x value at the midpoint of y.
#' @description For two vectors x and y of equal length, this will find the x value where y is at it's mid point.
#' @param .x A numeric vector
#' @param .y A numeric vector
#' @return The closest value for `x` that represents the midpoint of `y`.
#'
x_at_y_mid <- function(.x, .y) {
  if (length(.x) != length(.y)) {
    stop("x and y must have equal lengths")
  }

  mid <- min(.y) + (max(.y) - min(.y)) / 2
  index <- which.min(abs(.y - mid))

  return(.x[index])
}
