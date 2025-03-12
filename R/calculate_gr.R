calculate_gr <- function(.best_fit_values, .group, .control_name, .cap = TRUE) {
  # .best_fit_values <- predicted_all
  
  # groups
  .best_fit_values <- group_by(.best_fit_values, !!ensym(.group))
  .best_fit_values <- mutate(.best_fit_values, time_zero = first(y))
  # time zeroes
  control_time_zero <- filter(.best_fit_values, !!ensym(.group) == !!enquo(.control_name)) |> pull(time_zero) |> unique()
  .best_fit_values <- mutate(.best_fit_values, control_time_zero = control_time_zero)
  
  # normalised
  .best_fit_values <- mutate(.best_fit_values, normalised = y / time_zero)
  control_normalised <- filter(.best_fit_values, !!ensym(.group) == !!enquo(.control_name)) |> pull(normalised)
  .best_fit_values <- mutate(.best_fit_values, control_normalised = control_normalised)
  
  # calc
  .best_fit_values <- slice_head(.best_fit_values, n = -1)
  .best_fit_values <- mutate(.best_fit_values, gr = 2^((log2(normalised) / log2(control_normalised))) - 1)
  
  if (.cap) {
    .best_fit_values <- mutate(.best_fit_values, gr = case_when(gr > 1 ~ 1, .default = gr))
  }
  
  .best_fit_values <- ungroup(.best_fit_values)
  
  return(.best_fit_values)
}
