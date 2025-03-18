predict_data <- function(.fit, .lower_x, .upper_x, .group = NULL, .num_points = 1e3) {
  # full range of x values
  x_vals <- seq(.lower_x, .upper_x, length.out = .num_points)
  if (!is_null(.group)) {
    # unique group values
    .data <- .data |> mutate(group = !!ensym(.group))
    group_vals <- unique(.data$group)
    # go through each and collate
    predicted <- map(group_vals, function(i) {
      tibble(
        x = x_vals,
        y = predict(.fit, newdata = tibble(x = x_vals, group = i)),
        group = i
      )
    }) |>
      bind_rows()
  } else {
    predicted <- tibble(
      x = x_vals,
      y = predict(.fit, newdata = tibble(x = x_vals))
    )
  }
  return(predicted)
}
