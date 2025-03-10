extract_best_fit_values <- function(data_nest, y_var) {
  x_vals <- data_nest$data |>
    bind_rows() |>
    pull(x) |>
    unique() |>
    as_tibble_col(column_name = "x")

  pred_list <- pmap(
    .l = list(
      treat = data_nest$treatment_name,
      conc = data_nest$concentration,
      fit = data_nest[[y_var]] |> map(1)
    ),
    .f = function(treat, conc, fit) {
      if (is_null(fit)) {
        return()
      }
      broom::augment(
        x = fit,
        newdata = x_vals
      ) |>
        # select(1, 4) |>
        rename(x = 1, y = 2) |>
        distinct() |>
        mutate(
          treatment_name = treat,
          concentration = conc
          # test_utc_norm = unique(data$test_utc_norm)
        )
    }
  ) |>
    discard(.p = is_null)

  left_join(
    x = bind_rows(pred_list),
    y = plate_params$plate_setup |> select(-well),
    by = join_by(treatment_name, concentration),
    relationship = "many-to-many"
  ) |>
    distinct()
}
