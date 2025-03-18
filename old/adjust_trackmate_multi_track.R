adjust_trackmate_multi_track <- function(
    list,
    time_col = "datetime",
    multi_spot_fun = list(max = max)) {
  map(
    .x = list,
    .f = function(x) {
      x |>
        group_by(!!ensym(time_col)) |>
        summarise(
          across(
            .cols = everything(),
            .fns = multi_spot_fun,
            .names = "{.col}"
          )
        ) |>
        ungroup()
    }
  )
}
