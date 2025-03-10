library(minpack.lm)
library(rlang)
library(tidyverse)

roxygen2::roxygenise()

treatment_names <- c(
  "temozolomide",
  "doramapimod",
  "pp121",
  "vorinostat",
  "maramistat",
  "gabapentin",
  "selinexor",
  "s63845"
)
plate_params <- list(
  cell_line = "GL0101",
  start_date_passage = "2024/06/08",
  end_date_passage = "2024/06/19",
  start_n_passage = 1000000,
  end_n_passage = 4160000,
  microns_per_pixel = 1.24,
  time_between_scans = 4,
  plate_setup = bind_rows(
    tibble(
      treatment_name = treatment_names[1],
      concentration = rep(x = seq(from = -8, to = -4, by = 0.5), each = 3),
      well = outer(LETTERS[3:5], 3:11, paste0) |> as.vector()
    ),
    tibble(
      treatment_name = treatment_names[2],
      concentration = rep(x = seq(from = -8, to = -4, by = 0.5), each = 3),
      well = outer(LETTERS[6:8], 3:11, paste0) |> as.vector(),
    ),
    tibble(
      treatment_name = treatment_names[3],
      concentration = rep(x = seq(from = -8, to = -4, by = 0.5), each = 3),
      well = outer(LETTERS[9:11], 3:11, paste0) |> as.vector(),
    ),
    tibble(
      treatment_name = treatment_names[4],
      concentration = rep(x = seq(from = -8, to = -4, by = 0.5), each = 3),
      well = outer(LETTERS[12:14], 3:11, paste0) |> as.vector(),
    ),
    tibble(
      treatment_name = treatment_names[5],
      concentration = rep(x = seq(from = -8, to = -4, by = 0.5), each = 3),
      well = outer(LETTERS[3:5], 14:22, paste0) |> as.vector(),
    ),
    tibble(
      treatment_name = treatment_names[6],
      concentration = rep(x = seq(from = -8, to = -4, by = 0.5), each = 3),
      well = outer(LETTERS[6:8], 14:22, paste0) |> as.vector(),
    ),
    tibble(
      treatment_name = treatment_names[7],
      concentration = rep(x = seq(from = -8, to = -4, by = 0.5), each = 3),
      well = outer(LETTERS[9:11], 14:22, paste0) |> as.vector(),
    ),
    tibble(
      treatment_name = treatment_names[8],
      concentration = rep(x = seq(from = -9, to = -5, by = 0.5), each = 3),
      well = outer(LETTERS[12:14], 14:22, paste0) |> as.vector(),
    ),
    tibble(
      treatment_name = "media_n",
      concentration = rep(x = NA, each = 6),
      well = outer(LETTERS[3:5], 12:13, paste0) |> as.vector()
    ),
    tibble(
      treatment_name = "media_two_n",
      concentration = rep(x = NA, each = 6),
      well = outer(LETTERS[6:8], 12:13, paste0) |> as.vector()
    ),
    tibble(
      treatment_name = "vehicle",
      concentration = rep(x = NA, each = 6),
      well = outer(LETTERS[9:11], 12:13, paste0) |> as.vector()
    ),
    tibble(
      treatment_name = "staurosporine",
      concentration = rep(x = NA, each = 6),
      well = outer(LETTERS[12:14], 12:13, paste0) |> as.vector()
    )
  ) |>
    mutate(
      concentration = as.character(concentration)
    ) |>
    mutate(
      concentration = case_when(
        is.na(concentration) ~ treatment_name,
        .default = concentration
      )
    )
)

# import
spot_data <- import_trackmate_spot_data(
  paths = "/stornext/Bioinf/data/lab_brain_cancer/projects/drug_screen/incu/data/new/20240619b/data/test/composite/"
)

# filter to what we want
pp121_veh_wells <- plate_params$plate_setup |>
  filter(treatment_name %in% c("pp121", "vehicle")) |>
  pull(well)
spot_data <- spot_data[names(spot_data) %in% pp121_veh_wells]

# plate params
spot_data <- map(
  .x = spot_data,
  .f = function(x) {
    left_join(
      x = x,
      y = plate_params$plate_setup,
      by = join_by(well)
    )
  }
)

# area correction
spot_data <- map(
  .x = spot_data,
  .f = function(x) {
    x |> mutate(area = area * plate_params$microns_per_pixel^2)
  }
)

# wells
spot_data <- map(
  .x = spot_data,
  .f = function(x) {
    x |> mutate(
      row = str_extract(well, "[A-Za-z]+"),
      column = str_extract(well, "\\d+") |> as.numeric()
    )
  }
)

# bring together in long format
spot_collated <- bind_rows(spot_data)

spot_collated <- spot_collated |> filter(area > 2.5e4)

# time
spot_collated <- spot_collated |>
  mutate(
    hours_within = difftime(
      time1 = datetime,
      time2 = min(datetime),
      units = "hours"
    ) |>
      as.double.difftime()
  )

spot_collated

# go to models

vehicle_test <- spot_collated |> filter(treatment_name == "vehicle")

vehicle_results <- fit_model(
  .data = test,
  .x_var = hours_within,
  .y_var = area,
  .curve_func = incur_models$five_param_sigmoid$curve_func,
  .start_func = incur_models$five_param_sigmoid$start_func,
  .huber = TRUE,
  .detect_outliers = TRUE,
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)

# vehicle fit
plot_curve_fit(
  .data = vehicle_results$data,
  .fit = vehicle_results$fit,
  .x_var = hours_within,
  .y_var = area
)

# nested data
data_nest <- spot_collated |>
  nest(.by = c(concentration, treatment_name)) |>
  arrange(treatment_name, concentration)

data_nest <- data_nest |>
  mutate(
    area = map(
      .progress = TRUE,
      .x = data,
      .f = function(x) {
        # this returns a named list of data, formula and fit
        fit_model(
          .data = x,
          .x_var = hours_within,
          .y_var = area,
          .curve_func = incur_models$five_param_sigmoid$curve_func,
          .start_func = incur_models$five_param_sigmoid$start_func,
          .huber = TRUE,
          .detect_outliers = TRUE
        )
      }
    ),
    # this overwrites data column with added outlier column
    data = map(area, 2)
  )

plot_readout(
  data_nest,
  "hours_from",
  "area"
)

best_fit_values <- extract_best_fit_values(
  data_nest = data_nest,
  y_var = "area"
)

best_fit_values <- calculate_gr(best_fit_values = best_fit_values)

plot_gr_over_time(best_fit_values = best_fit_values)

gr_fit <- fit_model(
  .data = best_fit_values,
  .x_var = "concentration", 
  .y_var = "gr",
  .curve_func = incur_models$five_param_sigmoid_log$curve_func,
  .start_func = incur_models$five_param_sigmoid_log$start_func,
  .huber = TRUE,
  .detect_outliers = FALSE
)

plot_curve_fit(
  .data = gr_fit$data,
  .fit = gr_fit$fit,
  .x_var = concentration,
  .y_var = gr
)

find_root_fit(.fit = gr_fit$fit, .x_vals = gr_fit$data$x, .target = 0.5)

10^-6.895572 * 10e6

data_nest <- data_nest |> mutate(
  area_auc = map(
  .x = data_nest$area |> map(1),
  .f = function(x) {
    calculate_auc(
      .fit = x,
      .lower = min(spot_collated$hours_within),
      .upper = max(spot_collated$hours_within)
    )$value
  }
) |> 
  unlist(1)
)

auc_fit <- fit_model(
  .data = data_nest,
  .x_var = "concentration", 
  .y_var = "auc_res",
  .curve_func = incur_models$five_param_sigmoid_log$curve_func,
  .start_func = incur_models$five_param_sigmoid_log$start_func,
  .huber = TRUE,
  .detect_outliers = FALSE
)

plot_curve_fit(
  .data = auc_fit$data,
  .fit = auc_fit$fit,
  .x_var = concentration,
  .y_var = auc_res
)


