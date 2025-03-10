plot_curve_fit <- function(.data, .x_var, .y_var, .fit, .colour = NULL) {
  # raw
  data_raw <- mutate(.data, x = !!ensym(.x_var), y = !!ensym(.y_var))
  
  if(is.numeric(data_raw$x))
  
  # predicted
  x_vals <- unique(.data$x)

  x_vals_full <- seq(min(data_raw$x), max(data_raw$x), length.out = 1e4)
  data_predicted <- tibble(
    x = x_vals_full,
    y = predict(.fit, newdata = tibble(x = x_vals_full))
  )
  
  ggplot() +
    geom_point(aes(x = x, y = y, colour = !!ensym(.colour)), data_raw |> auto_convert_numerics()) +
    geom_line(aes(x = x, y = y), data_predicted)
}
