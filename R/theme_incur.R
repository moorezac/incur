# theme_incur
font_size <- 10
font_small <- 8
font_large <- 12
linewidth <- 1

theme_incur <- theme(
  axis.line = element_line(linewidth = linewidth / .pt),
  axis.text = element_text(size = font_small),
  axis.title.x = element_text(size = font_size),
  axis.title.y = element_text(size = font_size),
  # legend.box = "horizontal",
  # legend.location = "plot",
  # legend.justification = c(1, 1),
  legend.box.margin = margin(0, 0, 0, 0),
  legend.key = element_rect(colour = "black", linewidth = linewidth / .pt),
  legend.key.size = unit(1.25 * font_size, "pt"),
  legend.title = element_text(size = font_size),
  panel.grid = element_line(
    colour = "grey90",
    linewidth = linewidth / .pt
  ),
  panel.background = element_rect(fill = "white"),
  # panel.border = element_rect(
  #   colour = "black",
  #   fill = NA,
  #   linewidth = linewidth / .pt
  # ),
  panel.border = element_blank(),
  plot.title = element_text(size = font_size),
  text = element_text(size = font_size)
)
