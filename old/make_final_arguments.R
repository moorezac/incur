make_final_arguments <- function(.curve_func, .start_vals, .lower_bounds = NULL, .upper_bounds = NULL, .dots) {
  # formula
  function_formula <- make_formula(.curve_func)
  
  final_arguments <- list(
    formula = function_formula,
    start = purrr::list_flatten(.start_vals)
  )
  # bounds
  if (!purrr::is_null(.lower_bounds) | !purrr::is_null(.upper_bounds)) {
    bounds <- make_bounds(.curve_func, .lower_bounds, .upper_bounds)
    final_arguments <- append(final_arguments, bounds)
  }
  # dots dots dots
  final_arguments <- append(final_arguments, .dots)
  
  # final final
  final_arguments <- purrr::discard(final_arguments, is_null)
  
  return(final_arguments)
}
