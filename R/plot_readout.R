plot_readout <- function(.data_nest, .x_var, .y_var) {
  # raw data
  dat_raw <- map2(
    .x = data_nest$treatment_name,
    .y = data_nest$concentration,
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
    .x = data_nest$treatment_name,
    .y = data_nest$concentration,
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
    .x = unique(dat_raw$treatment_name) |> discard(.p = \(x) x == "vehicle"),
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
        #   axis.title.x = element_text(size = 8),
        #   axis.title.y = element_text(size = 8),
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

format_conc <- function(vec, signif = 2) {
  vec <- as.numeric(vec)
  map_vec(
    .x = vec,
    .f = function(x) {
      x_molar <- 10^x
      # convert to µM
      x_micro_molar <- x_molar * 1e6

      if (x_micro_molar >= 1) {
        # for concentrations ≥ 1 µM, display in µM
        # return(sprintf("%.2f µM", x_micro_molar))
        if (signif == 0) {
          return(str_c(round(x_micro_molar), " \u03bcM"))
        } else {
          return(str_c(format(x_micro_molar, digits = signif, nsmall = signif), " \u03bcM"))
        }
      } else {
        # for concentrations < 1 µM, convert to nM
        # return(sprintf("%.2f nM", x_nano_molar))
        x_nano_molar <- x_micro_molar * 1000
        if (signif == 0) {
          return(str_c(round(x_nano_molar), " nM"))
        } else {
          return(str_c(format(x_nano_molar, digits = signif, nsmall = signif), " nM"))
        }
      }
    }
  )
}
