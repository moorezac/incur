#' @title Make a formula expression from a function
#' @description From a provided function, create a formula in terms of `y ~ f({function arguments})`.
#' @param .func A function that describes a curve in terms of `x` and individual curve parameters.
#' @return A language object in the form of `y ~ f()`
#' @importFrom stringr str_c str_flatten str_replace

make_formula <- function(.func) {
  function_out <- capture.output(args(.func))
  function_out <- function_out[function_out != "NULL"]
  function_out[-1] <- trimws(function_out[-1])
  function_out <- stringr::str_flatten(function_out)
  function_out <- stringr::str_replace(function_out, "function ", "f")
  function_out <- stringr::str_c("y ~", function_out)
  function_out <- str2lang(function_out)

  return(function_out)
}
