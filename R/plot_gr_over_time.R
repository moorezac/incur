plot_gr_over_time <- function(best_fit_values) {
  # best_fit_values <- best_fit_values_area
  map(
    .x = unique(best_fit_values$treatment_name) |> discard(.p = \(x) x == "vehicle"),
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
        )  + 
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
        # ggh4x::force_panelsizes(unit(5, "cm"), unit(5, "cm"))
    }
  ) |>
    set_names(unique(best_fit_values$treatment_name) |> discard(.p = \(x) x == "vehicle"))
  
}
