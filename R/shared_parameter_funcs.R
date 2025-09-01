#' Create Starting Values for Shared Parameter Models
#'
#' @description
#' Generate appropriate starting parameter values when fitting models with
#' shared parameters across experimental groups.
#'
#' @param data A data frame containing the experimental data with group information.
#' @param model_func A function defining the mathematical model.
#' @param start_func A function for generating starting parameter values.
#' @param group_vec Character vector of unique group identifiers.
#' @param share_params Character vector of parameter names to be shared across groups.
#'
#' @return A named list of starting values suitable for shared parameter model fitting.
#'
#' @details
#' This function creates starting values by:
#' \itemize{
#'   \item Calculating global starting values for shared parameters using all data
#'   \item Calculating group-specific starting values for non-shared parameters
#'   \item Properly naming group-specific parameters with group suffixes
#'   \item Combining shared and group-specific values into a single list
#' }
#'
#' @family shared_parameters
#' @importFrom dplyr case_when filter pull
#' @importFrom purrr map imap list_flatten
#' @importFrom stringr str_c
#' @keywords internal
make_shared_start_values <- function(
  data,
  model_func,
  start_func,
  group_vec,
  share_params
) {
  # Global starting values for shared parameters
  start_values_global <- start_func(
    x = data |> dplyr::pull(x),
    y = data |> dplyr::pull(y)
  )

  # Group-specific starting values
  start_values_group <- purrr::map(group_vec, function(i) {
    start_func(
      x = data |> dplyr::filter(group == i) |> dplyr::pull(x),
      y = data |> dplyr::filter(group == i) |> dplyr::pull(y)
    )
  })
  names(start_values_group) <- group_vec

  # Modify names for group-specific parameters
  start_values_group <- purrr::imap(start_values_group, function(x, i) {
    names(x) <- dplyr::case_when(
      !names(x) %in% share_params ~ stringr::str_c(names(x), "_", i),
      .default = names(x)
    )
    # Keep only group-specific parameters
    x[!names(x) %in% names(start_values_global)]
  })

  # Combine shared and group-specific values
  start_values_global <- start_values_global[
    names(start_values_global) %in% share_params
  ]

  append(
    start_values_global,
    purrr::list_flatten(start_values_group, name_spec = "{inner}")
  )
}

#' Create Function Arguments for Shared Parameter Models
#'
#' @description
#' Modify function formal arguments to support shared parameters across groups
#' by creating group-specific parameter names for non-shared parameters.
#'
#' @param model_func A function defining the mathematical model.
#' @param group_vec Character vector of unique group identifiers.
#' @param share_params Character vector of parameter names to share across groups.
#'
#' @return A pairlist object suitable for creating a new function with
#'   group-specific parameters.
#'
#' @family shared_parameters
#' @importFrom rlang expr
#' @importFrom purrr map
#' @importFrom stringr str_c str_flatten_comma
#' @keywords internal
make_shared_formals <- function(func, group_vec, share_params) {
  formal_arguments <- names(formals(func))
  unique_arguments <- formal_arguments[
    !formal_arguments %in% c(share_params, "x")
  ]

  # Create group-specific parameter names
  append_arguments <- purrr::map(
    unique_arguments,
    stringr::str_c,
    "_",
    group_vec
  )

  # Validate shared parameters exist
  missing_params <- share_params[!share_params %in% formal_arguments]
  if (length(missing_params) > 0) {
    stop("Model params not found: ", str_flatten_comma(missing_params))
  }

  # Combine all arguments
  final_arguments <- c(
    "x",
    formal_arguments[formal_arguments %in% share_params],
    unlist(append_arguments),
    "group"
  )

  # Create formals list
  # TODO Is this jank??
  final_formals <- rep(list(rlang::expr()), length(final_arguments))
  names(final_formals) <- final_arguments

  as.pairlist(final_formals)
}

#' Create Function Body for Shared Parameter Models
#'
#' @description
#' Generate a function body that implements shared parameter logic using
#' `dplyr::case_when()` to assign appropriate parameter values based on group membership.
#'
#' @param model_func A function defining the mathematical model.
#' @param group_vec Character vector of unique group identifiers.
#' @param share_params Character vector of parameter names to share across groups.
#'
#' @return A call object representing the modified function body.
#'
#' @family shared_parameters
#' @importFrom dplyr case_when
#' @importFrom purrr map pmap
#' @importFrom stringr str_c str_glue
#' @keywords internal
make_shared_body <- function(func, group_vec, share_params) {
  formal_arguments <- names(formals(func))
  unique_arguments <- formal_arguments[
    !formal_arguments %in% c(share_params, "x")
  ]

  # Create parameter mapping expressions
  append_arguments <- purrr::map(
    unique_arguments,
    stringr::str_c,
    "_",
    group_vec
  )
  individual_group_list <- rep(list(group_vec), length(unique_arguments))

  # Generate case_when expressions for each parameter
  # Has to be vectorized - cannot use base if
  # TODO: Is there a better method to achieve this?
  expression_list <- purrr::pmap(
    .l = list(
      a = unique_arguments,
      b = individual_group_list,
      c = append_arguments
    ),
    .f = function(a, b, c) {
      stringr::str_c(
        stringr::str_glue("{a} <- dplyr::case_when("),
        stringr::str_c(
          stringr::str_glue("group == '{b}' ~ {c}"),
          collapse = ", "
        ),
        ")"
      )
    }
  )
  expression_list <- purrr::map(expression_list, str2lang)

  # Combine with original function body
  as.call(c(as.name("{"), expression_list, body(func)))
}

#' Prepare Complete Shared Parameter Model Setup
#'
#' @description
#' Execute the complete pipeline for creating shared parameter models including
#' data preparation, function modification, and starting value generation.
#'
#' @param data Input data frame.
#' @param share_params Character vector of parameters to share across groups.
#' @param share_group Character string specifying the grouping variable column name.
#' @param model_func Function defining the mathematical model.
#' @param start_func Function for generating starting parameter values.
#'
#' @return A list containing:
#' \describe{
#'   \item{[[1]]}{Modified data frame with standardized group column}
#'   \item{[[2]]}{Modified model function with shared parameter support}
#'   \item{[[3]]}{Starting parameter values for the shared parameter model}
#' }
#'
#' @details
#' This function orchestrates the complete shared parameter workflow:
#' \itemize{
#'   \item Handles existing 'group' columns by renaming to avoid conflicts
#'   \item Creates standardized group variable from specified column
#'   \item Generates appropriate starting values for shared and group-specific parameters
#'   \item Creates modified model function with group-aware parameter assignment
#'   \item Returns all components needed for shared parameter model fitting
#' }
#'
#' @family shared_parameters
#' @importFrom dplyr mutate relocate
#' @importFrom rlang ensym sym new_function
#' @keywords internal
prepare_shared_parameters <- function(
  data,
  share_params,
  share_group,
  model_func,
  start_func
) {
  # Handle existing 'group' column
  if ("group" %in% colnames(data)) {
    message("Detected 'group' column in data, renaming to 'group_original'")
    data <- dplyr::mutate(data, group_original = group)
  }

  # Add group column
  data <- dplyr::mutate(data, group = !!rlang::sym(share_group))
  data <- dplyr::mutate(data, group = as.character(group))
  data <- dplyr::relocate(data, group, .after = y)

  group_vec <- unique(data$group)

  start_values <- make_shared_start_values(
    data,
    model_func,
    start_func,
    group_vec,
    share_params
  )
  shared_arguments <- make_shared_formals(
    model_func,
    group_vec,
    share_params
  )
  shared_body <- make_shared_body(model_func, group_vec, share_params)

  # Create modified function
  model_func <- rlang::new_function(as.pairlist(shared_arguments), shared_body)

  list(data, model_func, start_values)
}
