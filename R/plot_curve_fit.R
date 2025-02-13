plot_curve_fit <- function(.data, .x_var, .y_var, .fit) {
  # raw
  data_raw <- mutate(.data, x = !!ensym(.x_var), y = !!ensym(.y_var))
  # predicted
  x_vals <- unique(.data$x)

  x_vals_full <- seq(min(.data$x), max(.data$x), length.out = 1e4)
  data_predicted <- tibble(
    x = a$data$x,
    y = predict(.fit)
  )

  ggplot() +
    geom_point(aes(x = x, y = y), data_raw) +
    geom_line(aes(x = x, y = y), data_predicted)
}
