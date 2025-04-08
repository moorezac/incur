#' @title Plot a model fitted with `incur`.
#' @description Plot a model fitted with `incur`.
#' @param data A `data.frame` or `data.frame` extension (tibble) in long format.
#' @param x_var A string that refers to the the `x` value within `data`.
#' @param y_var A string that refers to the the `y` value within `data`.
#' @param fit The fitted `nlsModel` object.
#' @param return_data Whether to instead return the original and predicted data to be plot elsewhere.
#' @return A `ggplot` object.
#' @importFrom dplyr mutate
#' @importFrom ggplot2 geom_line geom_point ggplot labs
#' @importFrom rlang as_name enquo ensym
#' @export
#' 
plot_model <- function(data, x_var, y_var, fit, return_data = FALSE) {
  # justin
  data <- dplyr::mutate(data, x = !!rlang::ensym(x_var), y = !!rlang::ensym(y_var))
  
  # check for group
  formula_vars <- all.vars(formula(fit))
  if ("group" %in% formula_vars) {
    group_vec <- unique(data$group)
    predicted <- predict_data(fit, min(data$x), max(data$x), group = group_vec)
    # gg <- ggplot2::ggplot(mapping = aes(x, y, colour = group))
    use_colour <- TRUE
  } else {
    predicted <- predict_data(fit, min(data$x), max(data$x))
    # gg <- ggplot2::ggplot(mapping = aes(x, y))
    use_colour <- FALSE
  }
  
  # check for outliers
  outlier_column <- stringr::str_c("outlier", rlang::as_name(rlang::enquo(y_var)), sep = "_")
  if (outlier_column %in% colnames(data)) {
    use_shape <- TRUE
  } else {
    use_shape <- FALSE
  }
  
  # aes
  map_point <- map_line <- aes(x = x, y = y)
  if (use_colour) {
    map_point$colour <- as.name("group")
    map_line$colour <- as.name("group")
  }
  if (use_shape) {
    map_point$shape <- as.name(outlier_column)
  }
  
  # if wanted to plot elsewhere
  if (return_data) {
    return(list(data = data, predicted = predicted))
  }
  
  # plot
  gg <- ggplot2::ggplot() +
    ggplot2::geom_point(data = data, mapping = map_point, alpha = 0.25) +
    ggplot2::geom_line(data = predicted, mapping = map_line) +
    ggplot2::labs(x = rlang::enquo(x_var), y = rlang::enquo(y_var))

  return(gg)
}

plot_models <- function(data_list, x_var, y_var, fit_list, nest_vec, nest_vec_name, colour_vector = NULL, label_vector = NULL, return_data = FALSE) {
  # data_list <- data_nest$data
  # x_var <- "hours_within"
  # y_var <- "area"
  # fit_list <- data_nest$area |> map(1)
  # nest_vec <- data_nest$concentration
  # nest_vec_name <- "treatment"
  
  all_same <- function(x) {
    if (!is.numeric(x)) {
      return(length(unique(x)) == 1)
    }
    return(var(x) == 0)
  }
  
  if (!all_same(c(length(data_list), length(fit_list), length(nest_vec)))) {
    stop("mismatch in lengths")
  }
  
  data_list <- map2(data_list, nest_vec, function(a, b) {
    # add in x and y
    a <- mutate(a, x = !!ensym(x_var), y = !!ensym(y_var))
    # add in names
    a <- mutate(a, !!ensym(nest_vec_name) := b)
    a <- relocate(a, !!ensym(nest_vec_name), .after = y)
    return(a)
  })
  data_collated <- bind_rows(data_list)
  
  data_predicted <- map2(data_list, fit_list, function(a, b) {
    # check for group
    formula_vars <- all.vars(formula(b))
    if ("group" %in% formula_vars) {
      stop("incur cannot plot more than one model that contains shared parameters")
    }
    
    # predicted data
    predicted <- predict_data(b, lower_x = min(a$x), max(a$x))
    # extract the nested_vec value and add in
    nested_vec <- pull(a, !!ensym(nest_vec_name))
    predicted <- mutate(predicted, !!ensym(nest_vec_name) := unique(nested_vec))
    
    return(predicted)
  })
  data_predicted <- bind_rows(data_predicted)
  
  if (return_data) {
    return(list(collated = data_collated, predicted = data_predicted))
  }
  
  set_1 <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999")
  set_3 <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
  
  
  # deal with outliers
  outlier_column <- str_c("outlier", rlang::as_name(enquo(y_var)), sep = "_")
  # colour_vector <- NULL
  
  if (outlier_column %in% colnames(data_collated)) {
    .shape_vector <- c(1, 16)
    names(.shape_vector) <- c("True", "False")
    .shape_labels <- c("True", "False")
    names(.shape_labels) <- c("True", "False")
    
    # add in - not sure we can set values on scale_shape with TRUE/FALSE?
    data_collated <- mutate(data_collated, !!ensym(outlier_column) := str_to_title(!!ensym(outlier_column)))
    
    gg <- ggplot(mapping = aes(x, y, colour = !!ensym(nest_vec_name), shape = !!ensym(outlier_column))) +
      guides(shape = guide_legend(title = "Outlier", override.aes = list(alpha = 1))) +
      scale_shape_manual(values = .shape_vector, labels = .shape_labels)
    data_predicted <- mutate(data_predicted, !!ensym(outlier_column) := "False")
  } else {
    gg <- gg <- ggplot(mapping = aes(x, y, colour = !!ensym(nest_vec_name)))
  }
  
  # plot
  gg <- gg +
    geom_point(data = data_collated, alpha = 0.25) +
    geom_line(mapping = aes(group = !!ensym(nest_vec_name)), data = data_predicted, linewidth = 1.1, colour = "black") +
    geom_line(data = data_predicted, linewidth = 1) +
    labs(x = enquo(x_var), y = enquo(y_var))
  
  if (!is_null(label_vector) & !is_null(colour_vector)) {
    gg <- gg + scale_colour_manual(values = colour_vector, labels = label_vector)
  } else {
    if (!is_null(colour_vector)) {
      gg <- gg + scale_colour_manual(values = colour_vector)
    } else {
      number_models <- length(fit_list)
      colour_vector <- case_when(
        between(number_models, 1, 9) ~ set_1[1:number_models],
        between(number_models, 10, 12) ~ set_3[1:number_models],
        .default = scales::hue_pal()(number_models)
      )
      names(colour_vector) <- nest_vec
      
      if (!is_null(label_vector)) {
        gg <- gg + scale_colour_manual(values = colour_vector, labels = label_vector)
      } else {
        gg <- gg + scale_colour_manual(values = colour_vector)
      }
    }
  }
  
  gg <- gg + theme_incur
  
  return(gg)
}

plot_readout <- function(data_nest, x_var, y_var, .conc_var) {
  # raw data
  dat_raw <- map2(
    x = data_nest$treatment_name,
    y = data_nest$concentration,
    .f = function(x, y) {
      data_nest |>
        filter(treatment_name == x & concentration == y) |>
        select(data) |>
        unnest(cols = data) |>
        mutate(treatment_name = !!x, concentration = !!y)
    }
  ) |>
    bind_rows()
  # there should be a mutate(x = !!enysm(x_var)) here
  # but it doens't work??
  # regardless x should be set
  
  # curve data
  dat_pred <- map2(
    x = data_nest$treatment_name,
    y = data_nest$concentration,
    .f = function(x, y) {
      fit <- data_nest |>
        filter(treatment_name == x & concentration == y) |>
        pull(!!ensym(y_var)) |>
        map(1) |>
        pluck(1)
      
      if (is_null(fit)) {
        return()
      }
      
      dat <- data_nest |>
        filter(treatment_name == x & concentration == y) |>
        pull(!!ensym(y_var)) |>
        map(2) |>
        pluck(1)
      
      pred <- tibble(
        x = seq(min(dat$x), max(dat$x), length.out = 1000)
      )
      pred <- pred |>
        mutate(pred = predict(fit, newdata = pred)) |>
        mutate(pred = as.numeric(pred)) |>
        mutate(treatment_name = !!x, concentration = !!y)
      
      return(pred)
    }
  ) |>
    discard(.p = is_null) |>
    bind_rows()
  # add in others
  # dat_pred <- left_join(
  #   dat_pred,
  #   dat_raw |> select(x, treatment_name, concentration, tint, test_divisions),
  #   join_by(x, treatment_name, concentration)
  # ) |>
  #   fill(tint)
  
  # there may or may not be outlier_well
  # there may or may not be outlier_col
  outlier_readout <- str_c("outlier_", y_var)
  outlier_well <- "outlier_well"
  
  dat_raw <- dat_raw |>
    mutate(
      outlier_total = case_when(
        !!ensym(outlier_readout) == TRUE ~ TRUE,
        !!ensym(outlier_well) == TRUE ~ TRUE,
        .default = FALSE
      )
    )
  
  if (!exists("colour_pal")) {
    colour_pal <- scales::pal_hue()(length(unique(data_nest$concentration))) |>
      set_names(unique(data_nest$concentration))
  }
  
  map(
    x = unique(dat_raw$treatment_name) |> discard(.p = \(x) x == "vehicle"),
    .f = function(treat) {
      # use pred to construct the legend
      concs_total <- data_nest |>
        filter(treatment_name == treat) |>
        pull(concentration) |>
        unique()
      concs_actual <- dat_pred |>
        filter(treatment_name == treat) |>
        pull(concentration) |>
        unique()
      concs_miss <- concs_total[!concs_total %in% concs_actual]
      
      # val <- with(
      #   colour_pal |>
      #     filter(treatment_name %in% c("vehicle", treat)) |>
      #     distinct(treatment_name, concentration, tint),
      #   set_names(tint, concentration)
      # )
      # val <- val[str_sort(names(val), decreasing = TRUE)]
      # val <- c(
      #   val,
      #   rep("red3", times = length(concs_miss)) |> set_names(concs_miss)
      # )
      # lab <- c("Vehicle", format_conc(concs_total))
      
      ggplot() +
        geom_point(
          mapping = aes(
            x = x,
            y = y,
            colour = concentration,
            shape = outlier_total
          ),
          data = dat_raw |>
            mutate(y = !!ensym(y_var)) |>
            filter(treatment_name %in% c("vehicle", treat)),
          alpha = 0,
          # size = 1
        ) +
        geom_point(
          mapping = aes(
            x = x,
            y = y,
            colour = concentration,
            fill = concentration,
            shape = outlier_total,
          ),
          data = dat_raw |>
            filter(treatment_name %in% c("vehicle", treat)) |>
            bind_rows(
              tibble(
                concentration = concs_miss,
                treatment_name = treat,
                outlier_total = TRUE
              )
            ) |>
            mutate(y = !!ensym(y_var)) |>
            mutate(outlier_total = factor(outlier_total, levels = c(TRUE, FALSE))),
          alpha = 1,
          # size = 1,
          key_glyph = draw_key_rect
        ) +
        geom_line(
          aes(
            x = x,
            y = y,
            colour = concentration
          ),
          data = dat_pred |>
            filter(treatment_name %in% c("vehicle", treat)) |>
            mutate(y = pred),
          linewidth = 2 / .pt,
          show.legend = FALSE
        ) +
        # this is a lot easier than other approaches
        # scale_colour_manual(
        #   values = val,
        #   labels = lab,
        #   limits = names(val)
        # ) +
        # scale_fill_manual(
        #   values = val,
        #   labels = lab,
        #   limits = names(val)
        # ) +
        # scale_shape_manual(values = c(1, 16)) +
        # scale_x_continuous(
        #   name = "Divisions",
        #   position = "bottom",
        #   sec.axis = dup_axis(
        #     name = "Hours",
        #     trans = ~ . * double_rate_incu,
        #     guide = guide_axis(position = "top")
        #   )
        # ) +
        ggtitle(treat) +
        # ylab(str_c(to_sentence_case(y_var), " (\u03bcm\u00B2)")) +
        guides(
          colour = "none",
          # fill = guide_legend(
          #   title = "Concentration",
          #   override.aes = list(alpha = 1),
          #   order = 1
          # ),
          shape = guide_legend(
            title = "outlier",
            override.aes = list(alpha = 1),
            order = 2
          )
        )
      # theme(
      #   axis.titlex = element_text(size = 8),
      #   axis.titley = element_text(size = 8),
      #   legend.box = "horizontal",
      #   legend.location = "plot",
      #   legend.justification = c(1, 0.85),
      #   legend.key = element_rect(
      #     colour = "black",
      #     linewidth = 1 / .pt
      #   ),
      #   legend.key.size = unit(0.5, "cm"),
      #   legend.title = element_text(size = 10),
      #   panel.grid = element_line(
      #     colour = "grey90",
      #     linewidth = 1 / .pt
      #   ),
      #   panel.background = element_rect(fill = "white"),
      #   panel.border = element_rect(
      #     colour = "black",
      #     fill = NA,
      #     linewidth = 1 / .pt
      #   ),
      #   plot.title = element_text(size = 10),
      #   text = element_text(size = 10)
      # )
    }
  ) |>
    set_names(unique(dat_raw$treatment_name) |> discard(.p = \(x) x == "vehicle"))
}

plot_gr_over_time <- function(best_fit_values) {
  # best_fit_values <- best_fit_values_area
  map(
    x = unique(best_fit_values$treatment_name) |> discard(.p = \(x) x == "vehicle"),
    .f = function(treat) {
      # use pred to construct the legend
      concs_total <- data_nest |>
        filter(treatment_name == treat) |>
        pull(concentration) |>
        unique()
      concs_actual <- best_fit_values |>
        filter(treatment_name == treat) |>
        pull(concentration) |>
        unique()
      concs_miss <- concs_total[!concs_total %in% concs_actual]
      
      # val <- with(
      #   colour_pal |>
      #     filter(treatment_name %in% c("vehicle", treat)) |>
      #     distinct(treatment_name, concentration, tint),
      #   set_names(tint, concentration)
      # )
      # val <- val[str_sort(names(val), decreasing = TRUE)]
      # val[concs_miss] <- "red3"
      # lab <- c("Vehicle", format_conc(concs_total))
      x_vals <- best_fit_values$x |> unique()
      
      ggplot() +
        geom_line(
          mapping = aes(
            x = x,
            y = gr,
            colour = concentration
          ),
          data = best_fit_values |>
            filter(treatment_name %in% c("vehicle", treat)) |>
            bind_rows(
              tibble(
                gr = -1,
                x = rep(x_vals, length(concs_miss)),
                concentration = rep(concs_miss, length(x_vals)),
                treatment_name = treat
              )
            ),
          linewidth = 1,
          key_glyph = draw_key_rect
        ) +
        # this is a lot easier than other approaches
        # scale_colour_manual(
        #   values = val,
        #   labels = lab,
        #   limits = names(val)
        # ) +
        # scale_x_continuous(
        #   name = "Divisions",
        #   position = "bottom",
        #   sec.axis = dup_axis(
        #     name = "Hours",
        #     trans = ~ . * double_rate_incu,
        #     guide = guide_axis(position = "top")
        #   )
        # ) +
        ylim(c(-1, 1)) +
        ggtitle(treat) +
        ylab("GR")
      # guides(
      #   colour = guide_legend(
      #     title = "Concentration",
      #     override.aes = list(alpha = 1),
      #     order = 1
      #   )
      # ) +
      # theme(
      #   axis.titlex = element_text(size = 8),
      #   axis.titley = element_text(size = 8),
      #   legend.box = "horizontal",
      #   legend.location = "plot",
      #   legend.justification = c(1, 0.85),
      #   legend.key = element_rect(
      #     colour = "black",
      #     linewidth = 1 / .pt
      #   ),
      #   legend.key.size = unit(0.5, "cm"),
      #   legend.title = element_text(size = 10),
      #   panel.grid = element_line(
      #     colour = "grey90",
      #     linewidth = 1 / .pt
      #   ),
      #   panel.background = element_rect(fill = "white"),
      #   panel.border = element_rect(
      #     colour = "black",
      #     fill = NA,
      #     linewidth = 1 / .pt
      #   ),
      #   plot.title = element_text(size = 10),
      #   text = element_text(size = 10)
      # )
      # ggh4x::force_panelsizes(unit(5, "cm"), unit(5, "cm"))
    }
  ) |>
    set_names(unique(best_fit_values$treatment_name) |> discard(.p = \(x) x == "vehicle"))
}
