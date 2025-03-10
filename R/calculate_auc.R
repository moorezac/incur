calculate_auc <- function(.fit, .lower, .upper) {
  integrate(
    f = function(x) predict(.fit, data.frame(x)),
    lower = .lower,
    upper = .upper
  )
}
