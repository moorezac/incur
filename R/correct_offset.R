#' @title Correct Offset Discontinuities in Time Series Data
#'
#' @description Detects and corrects for sudden jumps or offsets in time series
#' data by identifying outliers in slope changes and applying offset corrections
#' to realign displaced segments. The function uses robust smoothing and outlier
#' detection to identify discontinuities, then applies sequential corrections to
#' maintain trend continuity.
#'
#' @param x Numeric vector of x-axis values (typically time points).
#' @param y Numeric vector of y-axis values corresponding to x.
#'
#' @return Numeric vector of corrected y values with offset discontinuities
#'   removed. If no discontinuities are detected, returns the smoothed trend.
#'
#' @details The function works by:
#' \enumerate{
#'   \item Applying robust smoothing (3RS3R) to the input data
#'   \item Calculating slopes and second differences to identify sudden changes
#'   \item Using robust outlier detection to find discontinuity points
#'   \item Applying offset corrections to realign displaced segments
#' }
#'
#' @importFrom purrr map
#' @importFrom rlang is_null
#' @importFrom tibble tibble
#' @export
#'
#' @examples
#' \dontrun{
#' # Create data with an artificial offset
#' x <- 1:100
#' y <- sin(x / 10) + c(rep(0, 50), rep(2, 50))  # Jump at point 50
#'
#' # Correct the offset
#' y_corrected <- correct_offset(x, y)
#'
#' # Plot comparison
#' plot(x, y, col = "red")
#' lines(x, y_corrected, col = "blue")
#' }
correct_offset <- function(x, y) {
  smooth_trend <- smooth(y, kind = "3RS3R")
  smooth_trend <- as.numeric(smooth_trend)
  #smooth_trend <- loess(y ~ x)
  #smooth_trend <- smooth_trend$fitted |> plot()

  # Get slopes
  slopes <- diff(smooth_trend) / diff(x)
  dslopes <- diff(slopes)

  fit <- lm(
    y ~ x,
    tibble::tibble(x = seq_len(length(dslopes)), y = dslopes)
  )

  # Where are the outliers?
  indices <- incur::find_rout_indices(
    fit = fit,
    rout_opts = list(q = 1e-5, scale_method = "mad")
  )
  # plot(resid(fit))

  # Simple
  if (rlang::is_null(indices) || all(indices == 1)) {
    return(y)
  }

  indices <- indices + 1

  indices <- sort(indices)

  indices_filt <- indices[!indices %in% 1]
  # This should never be duplicated because of the diff funcs.
  indices_filt <- c(indices_filt, length(smooth_trend))

  # Get required pairs of indices
  pairs <- purrr::map(
    seq(from = 2, to = length(indices_filt), by = 2),
    function(i) {
      if (i == length(indices)) {
        c(indices_filt[i], length(y))
      } else {
        c(indices_filt[i - 1], indices_filt[i])
      }
    }
  )

  # If just the one
  if (length(pairs) == 1 & diff(pairs[[1]]) == 1) {
    pairs <- list(c(pairs[[1]][2], length(smooth_trend)))
  }

  # Process sequentially
  corrected <- y

  for (i in seq_len(length(pairs))) {
    # Get the end of previous segment and start of current segment
    prev_end <- pairs[[i]][1] - 1
    curr_start <- pairs[[i]][1]
    curr_end <- pairs[[i]][2]

    # Calculate offset (difference between end of previous and start of current)
    offset <- abs(corrected[prev_end] - corrected[curr_start])

    offset <- offset * sign(slopes[curr_start + 2])
    # Otherwise will be equal
    # offset <- offset + mean(slopes[prev_end - 1], slopes[curr_start + 1])

    # Apply correction to current segment
    corrected[curr_start:curr_end] <- corrected[curr_start:curr_end] + offset
  }

  # Last value NA
  if (is.na(corrected[length(corrected)])) {
    corrected <- corrected[-length(corrected)]
  }

  corrected
}
