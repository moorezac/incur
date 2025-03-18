
make_shared_parameters <- function(.data, .shared_params, .shared_group, .curve_func, .start_func) {
  
  if ("group" %in% colnames(.data)) {
    message("detected `group` as a col. name in provided data")
    message("`group` is a reserved col. name in incur")
    message("renaming `group` to `group_original`")
    .data <- dplyr::mutate(.data, group_orignal = group)
  }
  
  .data <- dplyr::mutate(.data, group = !!rlang::ensym(.shared_group))
  .data <- dplyr::mutate(.data, group = as.character(group))
  .data <- dplyr::relocate(.data, group, .after = y)
  
  group_vec <- unique(.data$group)
  
  .shared_arguments <- make_shared_formals(
    .func = .curve_func,
    .group = group_vec,
    .params = .shared_params
  )
  .shared_body <- make_shared_body(
    .func = .curve_func,
    .group = group_vec,
    .params = .shared_params
  )
  .curve_func <- rlang::new_function(
    args = as.pairlist(.shared_arguments),
    body = .shared_body
  )
  
  .start_vals <- make_shared_start_vals(
    .data = .data,
    .start_func = .start_func,
    .group = group_vec,
    .params = .shared_params
  )
  
  return(list(.data, .curve_func, .start_vals))
}
