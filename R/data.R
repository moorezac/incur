#' A collection of models with associated functions for both curves and starting value.
#'
#' @format A named list of length 10:
#' \describe{
#'   \item{curve_func}{A function that describes a curve in terms of `x`}
#'   \item{start_func}{A function that takes in vectors `x` and `y` of equal length and produces a named list corresponding to parameters in `curve_func`}
#'   }
"incur_models"
