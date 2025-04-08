# theme_incur
font_size <- 10
font_small <- 8
font_large <- 12
linewidth <- 1

theme_incur <- ggplot2::theme(
  axis.line = ggplot2::element_line(linewidth = linewidth / ggplot2::.pt),
  axis.text = ggplot2::element_text(size = font_small),
  axis.title.x = ggplot2::element_text(size = font_size),
  axis.title.y = ggplot2::element_text(size = font_size),
  # legend.box = "horizontal",
  # legend.location = "plot",
  # legend.justification = c(1, 1),
  legend.box.margin = ggplot2::margin(0, 0, 0, 0),
  legend.key = ggplot2::element_rect(colour = "black", linewidth = linewidth / ggplot2::.pt),
  legend.key.size = ggplot2::unit(1.25 * font_size, "pt"),
  legend.title = ggplot2::element_text(size = font_size),
  panel.grid = ggplot2::element_line(
    colour = "grey90",
    linewidth = linewidth / ggplot2::.pt
  ),
  panel.background = ggplot2::element_rect(fill = "white"),
  # panel.border = element_rect(
  #   colour = "black",
  #   fill = NA,
  #   linewidth = linewidth / .pt
  # ),
  panel.border = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(size = font_size),
  text = ggplot2::element_text(size = font_size)
)
