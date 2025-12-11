#' @param args A named list of arguments to pass to `minpack.lm::nlsLM()`.
#'
#' @param concentration_column Character string specifying the column containing
#'   drug concentrations in log10 molar format.
#'
#' @param curve_opts A named list of curve fitting options. Elements include:
#'   \itemize{
#'     \item `model`: Character string specifying a built-in model from \code{incur_models} or "loess".
#'     \item `model_func`: A function that describes a curve in terms of x (used if \code{model} is NA).
#'     \item `start_func`: A function that generates a named list of starting values for `model_func` (used if \code{model} is NA).
#'     \item `start_values`: A named list of starting values for `model_func`.
#'     \item `lower_bounds`: A named list that specifies the lower bounds for parameters in `model_func`.
#'     \item `upper_bounds`:A named list that specifies the upper bounds for parameters in `model_func`.
#'   }
#'
#' @param exclude_concs Numeric vector of concentrations to mark as excluded.
#'
#' @param func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#'
#' @param group_vec Character vector of unique group identifiers.
#'
#' @param metric Character string specifying the response metric to use,
#'   either \code{"gr"} for growth rate inhibition or \code{"ndr"} for
#'   normalised drug response. Default is \code{"gr"}.
#'
#' @param models Character vector of model names from `incur_models`.
#'
#' @param model_func A function that describes a curve in terms of `x`
#'   and individual curve parameters.
#'
#' @param model_list A named list of model results, typically output from
#'   \code{\link{interpolate_curve_concentration}}. Each element should contain
#'   \code{selected} (the fitted model) and \code{data} (the subset of data).
#'
#' @param negative_control_name Character string identifying the negative
#'   control (e.g., "DMSO" or "Vehicle") in \code{treatment_column}.
#'
#' @param obj A fitted `nls` object from `fit_curve()` or similar.
#'
#' @param outlier_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `huber`: Logical; use Huber robust regression? (default: FALSE).
#'     \item `huber_iter`: Maximum iterations for Huber (default: 100).
#'     \item `huber_k`: Huber tuning constant (default: 1.345).
#'     \item `huber_tol`: Convergence tolerance for Huber (default: 1e-6).
#'     \item `rout`: Logical; use ROUT outlier detection? (default: FALSE).
#'     \item `rout_q`: False discovery rate for ROUT (default: 1e-3).
#'     \item `rout_scale`: Scale estimator for ROUT, either "mad" or "quantile" (default: "mad").
#'   }
#'
#' @param positive_control_name Character string identifying the positive
#'   control (e.g., a cytotoxic agent) in \code{treatment_column}. If provided,
#'   NDR values are calculated. Default is NULL.
#'
#' @param share_params Character vector of parameter names to be shared.
#'
#' @param shared_opts A named list of shared parameter options. Elements include:
#'  \itemize{
#'     \item `share_group`: Character string specifying the column name for the grouping variable (default: NA).
#'     \item `share_params`: Character string specifying which parameters in the model are to be shared (default: NA).
#'     \item `return_func`: Logical; return fitted function? (default: FALSE).
#'   }
#'
#' @param start_colour Hex colour string for the highest concentration. Lower
#'   concentrations are displayed as progressively lighter tints. Default is
#'   \code{"#0085ca"}.
#'   
#' @param start_func A function for generating a maed list of starting parameter values.
#'
#' @param start_values A named list of starting parameter values.
#'
#' @param treatment_column Character string specifying the column containing
#'   treatment identifiers.
#'   
#' @param unique_concs Numeric vector of unique concentrations.
#'
#' @param x_var Character string specifying the column name for the independent
#'   variable (typically time).
#'
#' @param y_var Character string specifying the column name for the dependent
#'   variable (e.g., cell count, confluence, or a calculated metric).
#'   
NULL
