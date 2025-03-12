format_concentrations <- function(vec, signif = 2) {
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
