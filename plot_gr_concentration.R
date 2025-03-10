plot_gr_concentration <- function(best_fit_values, fits) {
  map(
    .x = unique(best_fit_values$treatment_name) |> discard(.p = \(x) x == "vehicle"),
    .f = function(treat) {
     fit <- fits[[treat]]$fit
     
     concs_total <- data_nest |>
       filter(treatment_name == treat) |>
       pull(concentration) |>
       unique()
     concs_actual <- best_fit_values |>
       filter(treatment_name == treat) |>
       pull(concentration) |>
       unique()
     concs_miss <- concs_total[!concs_total %in% concs_actual]
     
     x_vals <- best_fit_values |>
       pull(x) |>
       unique()
     
     dat <- best_fit_values |>
       filter(treatment_name == treat) |>
       filter(!is.nan(gr)) |>
       bind_rows(
         tibble(
           gr = -1,
           x = rep(x_vals, length(concs_miss)),
           concentration = rep(concs_miss, length(x_vals)),
           treatment_name = treat
         )
       ) |>
       mutate(concentration = as.numeric(concentration))
     
     dat_curve <- tibble(
       x = seq(
         min(dat$concentration),
         max(dat$concentration),
         length.out = 1000
       )
     )
     dat_curve <- dat_curve |> 
       mutate(gr = predict(fit, newdata = dat_curve)) |> 
       mutate(gr = as.numeric(gr), concentration = x)
     
     val <- with(
       colour_pal |>
         filter(treatment_name %in% c(treat)) |>
         distinct(treatment_name, concentration, tint),
       set_names(tint, concentration)
     )
     val <- val[str_sort(names(val), decreasing = TRUE)]
     lab <- c(format_conc(concs_total))
     
     intercept_50 <- try({
       find_root_fit(fit = fit, x_vals = dat$concentration, target = 0.5)$root
     })
     
     gg <- ggplot() +
       geom_point(
         mapping = aes(
           x = concentration,
           y = gr,
           colour = x
         ),
         data = dat,
         size = 1
       ) + 
       geom_line(
         aes(
           x = concentration,
           y = gr
         ),
         data = dat_curve,
         linewidth = 2 / .pt,
         show.legend = FALSE,
         colour = colour_pal |> 
           filter(treatment_name == treat) |> 
           pull(colour) |> 
           unique()
       ) + 
       scale_colour_viridis_c(
         guide = guide_colorbar(
           barheight = unit(5, "cm"),
           barwidth = unit(0.5, "cm"),
           title.position = "right",
           title = element_text("Divisions", size = 8 / .pt, angle = -90),
           ticks.colour = "black",
           frame.colour = "black",
           frame.linewidth = 0.2,
           reverse = TRUE
         )
         # direction = 1
       ) + 
       scale_x_continuous(
         labels = \(x) format_conc(x, signif = 0),
         breaks =  \(x) seq(ceiling(x[1]), floor(x[2]), by = 1),
         guide = guide_axis(angle = 45)
       ) + 
       ggtitle(treat) +
       xlab("Concentration") + 
       ylab("GR") +
       ylim(c(-1, 1)) + 
       theme(
         axis.title.x = element_text(size = 8),
         axis.title.y = element_text(size = 8),
         legend.box = "horizontal",
         legend.location = "plot",
         # legend.justification = c(1, 1.2),
         legend.key = element_rect(
           colour = "black",
           linewidth = 1 / .pt
         ),
         # legend.key.height = unit(1, "cm"),
         # legend.key.width = unit(0.5, "cm"),
         legend.title = element_text(size = 8, angle = -90, hjust = 0.5),
         panel.grid = element_line(
           colour = "grey90",
           linewidth = 1 / .pt
         ),
         panel.background = element_rect(fill = "white"),
         panel.border = element_rect(
           colour = "black",
           fill = NA,
           linewidth = 1 / .pt
         ),
         plot.title = element_text(size = 10),
         text = element_text(size = 10)
       ) +
       ggh4x::force_panelsizes(unit(5, "cm"), unit(5, "cm"))
     
     if( class(intercept_50) != "try-error") {
       gg <- gg + 
         geomtextpath::geom_textvline(
           xintercept = intercept_50,
           label = str_c(
             "GR\u2085\u2080 =	", format_conc(intercept_50, signif = 2)
           ),
           offset = -0.1,
           size = 8/.pt,
           colour = colour_pal |> 
             filter(treatment_name == treat) |> 
             pull(colour) |> 
             unique()
         )
     }
     return(gg)
    }
  ) |> 
    set_names(unique(best_fit_values$treatment_name) |> discard(.p = \(x) x == "vehicle"))
}
