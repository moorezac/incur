#' Calculate Growth Rate (GR) Values
#'
#' @description
#' Calculate growth rate (GR) values from time-series data, comparing treatment
#' groups to a control group. GR values provide a normalized measure of growth
#' inhibition that accounts for differences in growth rates between cell lines.
#'
#' @param data A data frame containing the time-series data.
#' @param x_var Character string specifying the name of the time variable column.
#' @param y_var Character string specifying the name of the response variable column.
#' @param share_group Character string specifying the column name for grouping (e.g., treatment groups).
#' @param control_name The value in `share_group` that represents the control condition.
#' @param cap Logical. If `TRUE` (default), cap GR values at 1 to prevent values > 1.
#' @param return_all_columns Logical. If `TRUE`, return all original columns plus GR values.
#'   If `FALSE` (default), return only original columns plus GR values.
#'
#' @return A data frame with GR values added. GR = 0 indicates no growth inhibition,
#'   GR = 1 indicates complete growth inhibition, and GR = -1 indicates cell death.
#'
#' @details
#' The GR metric normalizes growth inhibition relative to a control condition,
#' making it suitable for comparing effects across different cell lines or
#' experimental conditions. The calculation uses the formula:
#' GR = 2^((log2(normalized_treatment) / log2(normalized_control))) - 1
#'
#' @family growth_analysis
#' @importFrom dplyr case_when filter first group_by mutate pull slice_head ungroup select
#' @importFrom rlang enquo ensym
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate GR values comparing treatments to vehicle control
#' gr_data <- calculate_gr(
#'   data = cell_data,
#'   x_var = "time_hours",
#'   y_var = "cell_area",
#'   share_group = "treatment",
#'   control_name = "vehicle"
#' )
#' }
calculate_gr <- function(
  data,
  x_var,
  y_var,
  share_group,
  control_name,
  cap = TRUE,
  return_all_columns = FALSE
) {
  # Original columns
  original_columns <- colnames(data)

  # Prepare data
  data <- prepare_data(data, x_var, y_var)

  # Group data and calculate time zero values
  data <- dplyr::group_by(data, !!rlang::ensym(share_group))
  data <- dplyr::mutate(data, time_zero = dplyr::first(y))

  # Extract time zero and add to all groups
  control_time_zero <- dplyr::filter(
    data,
    !!rlang::ensym(share_group) == !!rlang::enquo(control_name)
  ) |>
    dplyr::pull(time_zero) |>
    unique()
  data <- dplyr::mutate(
    data,
    control_time_zero = control_time_zero
  )

  # Extract normalised and add to all groups
  data <- dplyr::mutate(data, normalised = y / time_zero)
  control_normalised <- dplyr::filter(
    data,
    !!rlang::ensym(share_group) == !!rlang::enquo(control_name)
  ) |>
    dplyr::pull(normalised)
  data <- dplyr::mutate(data, control_normalised = control_normalised)

  # Calculate GR values
  data <- dplyr::slice_head(data, n = -1)
  data <- dplyr::mutate(
    data,
    gr = 2^((log2(normalised) / log2(control_normalised))) - 1
  )
  # Apply cap if requested
  if (cap) {
    data <- dplyr::mutate(
      data,
      gr = dplyr::case_when(gr > 1 ~ 1, .default = gr)
    )
  }

  # Remove grouping
  data <- dplyr::ungroup(data)

  # Select relevant colums
  if (return_all_columns) {
    data
  } else {
    dplyr::select(data, !!original_columns, gr)
  }
}
