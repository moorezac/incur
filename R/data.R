#' Built-in Model Functions and Starting Value Functions
#'
#' A collection of mathematical models commonly used in biological curve fitting,
#' along with their associated starting value estimation functions.
#'
#' @format A named list with 10 model entries, each containing:
#' \describe{
#'   \item{model_func}{A function that describes a curve in terms of `x` and model parameters}
#'   \item{start_func}{A function that takes vectors `x` and `y` of equal length and
#'     produces a named list of starting parameter values for `model_func`}
#' }
#'
#' @details
#' Available models include:
#' - `exponential_growth`: Simple exponential growth model
#' - `exponential_plateau`: Exponential approach to plateau
#' - `four_param_sigmoid`: Four-parameter sigmoidal dose-response curve
#' - `four_param_sigmoid_log`: Four-parameter sigmoid with log-transformed EC50
#' - `biphasic_four_param_sigmoid`: Biphasic dose-response curve
#' - `five_param_sigmoid`: Five-parameter sigmoid with asymmetry parameter
#' - `five_param_sigmoid_log`: Five-parameter sigmoid with log EC50
#' - `logistic_growth`: Logistic growth model
#' - `gompertz_growth`: Gompertz growth model
#' - `beta_growth_decay`: Beta growth-decay model
#'
#' @family model_functions
"incur_models"

#' Example Concentration Response Dataset
#'
#' A sample dataset demonstrating concentration-response data structure
#' for use with incur modeling functions.
#'
#' @format A data frame with concentration-response experimental data
#' @family example_data
"conc_data"

#' Example Cell Doubling Dataset
#'
#' A sample dataset for demonstrating cell doubling rate calculations
#' and growth analysis functions.
#'
#' @format A data frame with time-series cell growth data
#' @family example_data
"double_data"

#' Example SNP Dataset
#'
#' A sample dataset for demonstrating cell line authentication
#' functions.
#'
#' @format A named list with SNP array data
#' @family example_data
"snp_data"
