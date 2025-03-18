plot_model <- function(
    .data,
    .x_var,
    .y_var,
    .fit,
    .return_data) {
  .data <- mutate(.data, x = !!ensym(.x_var), y = !!ensym(.y_var))

  all.vars(formula(.fit))

  # check for group
  formula_vars <- all.vars(formula(.fit))
  if ("group" %in% formula_vars) {
    predicted <- predict_data(.fit, .lower_x = min(.data$x), max(.data$x), .group = "group")
    gg <- ggplot(mapping = aes(x, y, colour = group))
  } else {
    predicted <- predict_data(.fit, .lower_x = min(.data$x), max(.data$x))
    gg <- ggplot(mapping = aes(x, y))
  }

  if (return_data) {
    return(list(collated = data_collated, predicted = data_predicted))
  }

  # plot
  gg <- gg +
    geom_point(data = .data) +
    geom_line(data = predicted) +
    labs(x = enquo(.x_var), y = enquo(.y_var))

  return(gg)
}

plot_models <- function(
    .data_list,
    .x_var,
    .y_var,
    .fit_list,
    .nest_vec,
    .nest_vec_name,
    .colour_vector = NULL,
    .label_vector = NULL,
    .return_data = FALSE) {
  # .data_list <- data_nest$data
  # .x_var <- "hours_within"
  # .y_var <- "area"
  # .fit_list <- data_nest$area |> map(1)
  # .nest_vec <- data_nest$concentration
  # .nest_vec_name <- "treatment"

  all_same <- function(x) {
    if (!is.numeric(x)) {
      return(length(unique(x)) == 1)
    }
    var(x) == 0
  }

  if (!all_same(c(length(.data_list), length(.fit_list), length(.nest_vec)))) {
    stop("mismatch in lengths")
  }

  .data_list <- map2(.data_list, .nest_vec, function(a, b) {
    # add in x and y
    a <- mutate(a, x = !!ensym(.x_var), y = !!ensym(.y_var))
    # add in names
    a <- mutate(a, !!ensym(.nest_vec_name) := b)
    a <- relocate(a, !!ensym(.nest_vec_name), .after = y)
    return(a)
  })
  data_collated <- bind_rows(.data_list)

  data_predicted <- map2(.data_list, .fit_list, function(a, b) {
    # check for group
    formula_vars <- all.vars(formula(b))
    if ("group" %in% formula_vars) {
      stop("incur cannot plot more than one model that contains shared parameters")
    }

    # predicted data
    predicted <- predict_data(b, .lower_x = min(a$x), max(a$x))
    # extract the nested_vec value and add in
    nested_vec <- pull(a, !!ensym(.nest_vec_name))
    predicted <- mutate(predicted, !!ensym(.nest_vec_name) := unique(nested_vec))

    return(predicted)
  })
  data_predicted <- bind_rows(data_predicted)

  if (.return_data) {
    return(list(collated = data_collated, predicted = data_predicted))
  }

  set_1 <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")
  set_3 <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")


  # deal with outliers
  outlier_column <- str_c("outlier", rlang::as_name(enquo(.y_var)), sep = "_")
  # colour_vector <- NULL

  if (outlier_column %in% colnames(data_collated)) {
    .shape_vector <- c(1, 16)
    names(.shape_vector) <- c("True", "False")
    .shape_labels <- c("True", "False")
    names(.shape_labels) <- c("True", "False")

    # add in - not sure we can set values on scale_shape with TRUE/FALSE?
    data_collated <- mutate(data_collated, !!ensym(outlier_column) := str_to_title(!!ensym(outlier_column)))

    gg <- ggplot(mapping = aes(x, y, colour = !!ensym(.nest_vec_name), shape = !!ensym(outlier_column))) +
      guides(shape = guide_legend(title = "Outlier", override.aes = list(alpha = 1))) +
      scale_shape_manual(values = .shape_vector, labels = .shape_labels)
    data_predicted <- mutate(data_predicted, !!ensym(outlier_column) := "False")
  } else {
    gg <- gg <- ggplot(mapping = aes(x, y, colour = !!ensym(.nest_vec_name)))
  }

  # plot
  gg <- gg +
    geom_point(data = data_collated, alpha = 0.25) +
    geom_line(mapping = aes(group = !!ensym(.nest_vec_name)), data = data_predicted, linewidth = 1.1, colour = "black") +
    geom_line(data = data_predicted, linewidth = 1) +
    labs(x = enquo(.x_var), y = enquo(.y_var))

  if (!is_null(.label_vector) & !is_null(.colour_vector)) {
    gg <- gg + scale_colour_manual(values = .colour_vector, labels = .label_vector)
  } else {
    if (!is_null(.colour_vector)) {
      gg <- gg + scale_colour_manual(values = .colour_vector)
    } else {
      number_models <- length(.fit_list)
      .colour_vector <- case_when(
        between(number_models, 1, 9) ~ set_1[1:number_models],
        between(number_models, 10, 12) ~ set_3[1:number_models],
        .default = scales::hue_pal()(number_models)
      )
      names(.colour_vector) <- .nest_vec

      if (!is_null(.label_vector)) {
        gg <- gg + scale_colour_manual(values = .colour_vector, labels = .label_vector)
      } else {
        gg <- gg + scale_colour_manual(values = .colour_vector)
      }
    }
  }

  gg <- gg + theme_incur

  return(gg)
}
