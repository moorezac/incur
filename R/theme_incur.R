#' ggplot2 Theme for incur Package
#'
#' @description
#' A clean, publication-ready ggplot2 theme designed for scientific plots
#' created by the incur package. Provides consistent styling across all
#' plotting functions.
#'
#' @details
#' The theme features:
#' \itemize{
#'   \item Clean white background with subtle grid lines
#'   \item Consistent font sizing (10pt base, 8pt small, 12pt large)
#'   \item Black borders on legend keys for clarity
#'   \item Minimal margins and spacing for efficient plot real estate
#'   \item Axis lines instead of full plot borders
#'   \item Grey grid lines for subtle reference
#' }
#'
#' @family visualization
#' @importFrom ggplot2 .pt theme element_line element_text element_rect margin unit
#' @export
#'
#' @examples
#' \dontrun{
#' # Apply incur theme to any ggplot
#' library(ggplot2)
#' ggplot(data, aes(x, y)) +
#'   geom_point() +
#'   theme_incur()
#' }
theme_incur <- function(
  font_size = 10,
  font_small = 8,
  font_large = 12,
  linewidth = 1
) {
  ggplot2::theme(
    axis.line = ggplot2::element_line(linewidth = linewidth / ggplot2::.pt),
    axis.text = ggplot2::element_text(size = font_small),
    axis.title.x = ggplot2::element_text(size = font_size),
    axis.title.y = ggplot2::element_text(size = font_size),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.key = ggplot2::element_rect(
      colour = "black",
      linewidth = linewidth / ggplot2::.pt
    ),
    legend.key.size = ggplot2::unit(1.25 * font_size, "pt"),
    legend.title = ggplot2::element_text(size = font_size),
    panel.grid = ggplot2::element_line(
      colour = "grey90",
      linewidth = linewidth / ggplot2::.pt
    ),
    panel.background = ggplot2::element_rect(fill = "white"),
    panel.border = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(size = font_size),
    text = ggplot2::element_text(size = font_size)
  )
}
