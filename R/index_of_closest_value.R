#' Find the index of a vector closest to a target value
#'
#' @param vector A numeric vector to search within.
#' @param target The target value to find the closest value to.
#' @return An index of `vector` that contains the closest value to `target`
#'
#' @export

index_of_closest_value <- function(vector, target) {
  which(abs(vector - target) == min(abs(vector - target)))
}
