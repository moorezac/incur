#' Make shared function formals across arguments and individual groups.
#'
#' @param .func A function for which arguments are to be modified.
#' @param .group A character vector that contains all unique grouping values.
#' @param .params A character vector of arguments in `.func` that are to be shared across `.groups`.
#'
#' @return An `pairlist` object that that can be used to create a new function
#'

make_shared_formals <- function(.func, .group, .params) {
  formal_arguments <- names(formals(.func))
  unique_arguments <- formal_arguments[!formal_arguments %in% c(.params, "x")]
  append_arguments <- make_append_arguments(unique_arguments, .group)

  if (!any(.params %in% formal_arguments)) {
    missing <- .params[!.params %in% formal_arguments]
    stop(str_c("model params not found: ", str_flatten_comma(missing)))
  }

  final_arguments <- c(
    "x",
    formal_arguments[formal_arguments %in% .params],
    unlist(append_arguments),
    "group"
  )

  # use this to create alist
  # TODO is this jank?
  final_formals <- rep(list(rlang::expr()), length(final_arguments))
  names(final_formals) <- final_arguments

  as.pairlist(final_formals)
}

make_append_arguments <- function(.unique, .group) {
  len <- length(.group)
  map(.unique, str_c, "_", .group)
}

# create_formals_shared(func = f, shared_args = shared_args)
