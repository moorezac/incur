calculate_gr <- function(best_fit_values, cap = TRUE) {
  # time zeroes
  dat <- best_fit_values |>
    group_by(treatment_name, concentration) |>
    mutate(time_zero = dplyr::first(y)) %>%
    mutate(
      vehicle_time_zero = {
        . |>
          filter(treatment_name == "vehicle") |>
          pull(time_zero)
      }
    )

  # normalised
  dat <- dat |>
    group_by(treatment_name, concentration) |>
    mutate(normalised = y / time_zero) %>%
    mutate(
      vehicle_normalised = {
        . |>
          filter(treatment_name == "vehicle") |>
          pull(normalised)
      }
    ) |>
    ungroup()

  dat <- dat |>
    group_by(treatment_name, concentration) |>
    slice_head(n = -1) |>
    # slice_tail(n = -1) |>
    mutate(gr = 2^((log(x = normalised, base = 2) / log(x = vehicle_normalised, base = 2))) - 1) |>
    ungroup() |>
    relocate(gr)

  if (cap) {
    dat <- dat |>
      mutate(
        gr = case_when(
          gr > 1 ~ 1,
          .default = gr
        )
      )
  }

  dat
}
