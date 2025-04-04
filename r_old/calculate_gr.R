#' @title Calculate GR values.
#' @description A function that finds area under the curve for a fitted `nls` object across a lower and upper bound.
#' @param .best_fit_values A fitted `nlsLM` object.
#' @param .lower The lower bound of where to calculate AUC.
#' @param .upper The upper bound of where to calculate AUC.
#' @return A list of class `integrate` as produced by `stats::integrate`.
#' @importFrom dplyr case_when filter first group_by mutate pull slice_head ungroup
#' @importFrom rlang enquo ensym
#' @export

calculate_gr <- function(.best_fit_values, .group, .control_name, .cap = TRUE) {
  # .best_fit_values <- predicted_all

  # groups
  .best_fit_values <- dplyr::group_by(.best_fit_values, !!rlang::ensym(.group))
  .best_fit_values <- dplyr::mutate(.best_fit_values, time_zero = dplyr::first(y))
  # time zeroes
  control_time_zero <- dplyr::filter(.best_fit_values, !!rlang::ensym(.group) == !!rlang::enquo(.control_name)) |>
    dplyr::pull(time_zero) |>
    unique()
  .best_fit_values <- dplyr::mutate(.best_fit_values, control_time_zero = control_time_zero)

  # normalised
  .best_fit_values <- dplyr::mutate(.best_fit_values, normalised = y / time_zero)
  control_normalised <- dplyr::filter(.best_fit_values, !!rlang::ensym(.group) == !!rlang::enquo(.control_name)) |> dplyr::pull(normalised)
  .best_fit_values <- dplyr::mutate(.best_fit_values, control_normalised = control_normalised)

  # calc
  .best_fit_values <- dplyr::slice_head(.best_fit_values, n = -1)
  .best_fit_values <- dplyr::mutate(.best_fit_values, gr = 2^((log2(normalised) / log2(control_normalised))) - 1)

  if (.cap) {
    .best_fit_values <- dplyr::mutate(.best_fit_values, gr = dplyr::case_when(gr > 1 ~ 1, .default = gr))
  }

  .best_fit_values <- dplyr::ungroup(.best_fit_values)

  return(.best_fit_values)
}
