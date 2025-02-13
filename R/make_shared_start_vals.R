make_shared_start_vals <- function(.data, .start_func, .group, .params) {
  # use all data for shared params
  start_vals_global <- .start_func(
    x = .data |> pull(x),
    y = .data |> pull(y)
  )

  start_vals_group <- map(.group, function(i) {
    .start_func(
      x = .data |> filter(group == i) |> pull(x),
      y = .data |> filter(group == i) |> pull(y)
    )
  })
  names(start_vals_group) <- .group

  # make names
  start_vals_group <- imap(start_vals_group, function(x, i) {
    names(x) <- case_when(
      !names(x) %in% .params ~ str_c(names(x), "_", i),
      .default = names(x)
    )
    # filter to group values
    x[!names(x) %in% names(start_vals_global)]
  })

  # append
  start_vals_global <- start_vals_global[names(start_vals_global) %in% .params]

  append(start_vals_global, list_flatten(start_vals_group, name_spec = "{inner}"))
}
