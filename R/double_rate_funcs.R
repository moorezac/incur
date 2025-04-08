#' @title Calculate the double rate for area measurements of spheroids over time.
#' @description From a set of matched values that represent area/growth readouts where one group is a known double of another, calculate a doubling time. This function was designed to measure spheroid area from IncuCyte images were cells were initially seeded at 1,000 and 2,000 cells separately. It is recommended that data be generated from the `fit_model` function with shared parameters, and then `predict_data` used to generate a large number of points to aid in calculations.
#' @param data A `data.frame` or `data.frame` extension (tibble) in long format.
#' @param x_var A character string that refers to the the `x` value within `data`.
#' @param y_var A character string that refers to the the `y` value within `data`.
#' @param shared_group A string for the column in `data` that refers to individual groups to share fitted parameters.
#' @param n A character string in `shared_group` that refers to samples which represent samples composed of baseline numbers of cells.
#' @param two_n A character string in `shared_group` that refers to samples which represent samples composed of two times baseline numbers of cells.
#'
calc_double_rate_group <- function(data, x_var = "x", y_var = "y", shared_group = "group", n, two_n, x_lab = "x", y_lab = "y") {
  # data <- double_fit_shared$data
  
  # mutate
  data <- dplyr::mutate(data, x = !!rlang::ensym(x_var), y = !!rlang::ensym(y_var))
  data <- dplyr::mutate(data, group = !!rlang::sym(shared_group))
  
  # checks
  group_vec <- unique(data$group)
  if (!n %in% group_vec && !two_n %in% group_vec) {
    stop("provided 'n' and/or 'two_n' not found in provided 'shared_group'")
  }
  
  # extract data
  data_x <- dplyr::filter(data, group == n) |> dplyr::pull(x)
  data_n <- dplyr::ilter(data, group == n) |> dplyr::pull(y)
  data_two_n <- dplyr::filter(data, group == two_n) |> dplyr::pull(y)
  
  # ratio approach
  data_ratio <- log(data_two_n / data_n)
  
  # if assume sigmoid growth
  # when two_n > n the ratio max is around two_n inflection point
  max_index <- which(ratio == max(ratio))
  data_x_filt <- data_x[0:max_index]
  data_ratio_filt <- data_ratio[0:max_index]
  # plot(data_ratio_filt ~ data_x_filt)
  
  # find where section is linear
  # this represents exponential growth of ratio
  inflection_points <- inflection::ese(data_x_filt, data_ratio_filt, 0)
  linear_indexes <- inflection_points[1]:inflection_points[2]
  # plot(data_ratio_filt[linear_indexes] ~ data_x_filt[linear_indexes])
  
  # fit a curve
  linear_fit <- lm(data_ratio_filt[linear_indexes] ~ data_x_filt[linear_indexes])
  
  # to double spherical volume radius must increase 2^(1/3)
  # to double spherical volume area must increase 2^(2/3)
  double_rate <- log(2)/((2^(2/3)) * coef(linear_fit)[[2]])
  
  # plot
  df_overall <- dplyr::mutate(data, linear = dplyr::case_when(x %in% data_x[linear_indexes] ~ TRUE, .default = FALSE))
  df_overall <- dplyr::mutate(df_overall, group = dplyr::case_when(linear == TRUE ~ "linear", .default = group))
  gg_overall <- ggplot2::ggplot(df_overall, mapping = aes(x, y, colour = group)) + 
    ggplot2::geom_point() + 
    ggplot2::scale_colour_manual(values = c("#FFC20A", "#0C7BDC", "#40B0A6") |> rlang::set_names(c(n, two_n, "linear")), breaks = c(n, two_n)) + 
    ggplot2::geom_vline(xintercept = c(min(data_x_filt[linear_indexes]), max(data_x_filt[linear_indexes])), linetype = "dashed") + 
    ggplot2::labs(x = x_lab, y = y_lab) + 
    theme_incur
  
  df_ratio <- data.frame(x = data_x, y = data_ratio)
  df_ratio <- dplyr::mutate(df_ratio, linear = dplyr::case_when(x %in% data_x[linear_indexes] ~ TRUE, .default = FALSE))
  gg_ratio <- ggplot2::ggplot(df_ratio, mapping = aes(x, y, colour = linear)) + 
    ggplot2::geom_point() + 
    ggplot2::scale_colour_manual(values = c("#40B0A6", "#E1BE6A") |> rlang::set_names(TRUE, FALSE)) + 
    # geom_abline(intercept = coef(linear_fit)[[1]], slope = coef(linear_fit)[[2]]) +
    geomtextpath::geom_textabline(
      intercept = coef(linear_fit)[[1]], slope = coef(linear_fit)[[2]],
      label = (stringr::str_c(
        signif(coef(linear_fit)[[1]], digits = 3), "\u00d7", "x + ", signif(coef(linear_fit)[[2]], digits = 3)
        )),
      size = 6/ggplot2::.pt,
      vjust = 1.5
    ) +
    ggplot2::labs(x = x_lab, y = str2expression(str_c("ln(", two_n, "/", n, ")")))  + 
    theme_incur
  
  gg <- patchwork::wrap_plots(gg_overall, gg_ratio, guides = "collect", axis_titles = "collect") + 
    ggplot2::theme(aspect.ratio = 1)
  
  return(list(double_rate = double_rate, double_plot = gg))
  
  # bad approach
  # # where is the inflection point aka steepest part
  # # this is NOT the ec50!
  # curve_index <- inflection::check_curve(data_x, data_n)$index
  # inflection_n <- inflection::ese(data_x, data_n, index = curve_index)
  # inflection_two_n <- inflection::ese(data_x, data_two_n, index = curve_index)
  # inflection_diff <- inflection_n[3] - inflection_two_n[3]
}

calc_double_rate_passage <- function(start_date, end_date, start_n, end_n) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  days_between <-
    difftime(
      time1 = end_date,
      time2 = start_date
    ) |>
    as.numeric()
  
  return(days_between * log(2) / (log(end_n / start_n)) * 24)
}
