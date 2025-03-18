#' @title Attempt to fit a curve via `minpack.lm`
#' @description Try to fit a curve within a new environment. See help page on `rlang::exec` for background on this approach.
#' @param .data A A `data.frame` or `data.frame` extension (tibble) in long format with columns `x` and `y`. 
#' @param .curve_func A function that describes a curve/model in terms of `x`.
#' @param .arguments A named list of arguments to be injected into the `minpack.lm::nlsLM` function.
#' @return Either a fitted `nls` object if expression is evaluated without error, or an invisible object from class `try-error`.
#' @importFrom rlang env expr
#' @importFrom minpack.lm nlsLM
#' 
try_fit <- function(.data, .curve_func, .arguments) {
  f <- .curve_func
  fit <- try({
    new_env <- rlang::env(dat = .data)
    eval(
      rlang::expr(minpack.lm::nlsLM(data = dat, !!!.arguments)),
      envir = new_env
    )
  })
  
  return(fit)
}
