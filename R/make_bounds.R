make_bounds <- function(.func, .lower = NULL, .upper = NULL) {
  formal_arguments <- names(formals(.func))
  formal_arguments <- formal_arguments[!formal_arguments %in% c("x", "group")]

  # if only one is supplied
  final_list <- list(lower = .lower, upper = .upper) |> discard(is_null)

  # check bound args are there
  iwalk(final_list, function(x, i) {
    if (!any(names(x) %in% formal_arguments)) {
      missing <- names(x)[!names(x) %in% formal_arguments]
      stop(str_c(i, " arg(s) not found in : ", str_flatten_comma(missing)))
    }
  })

  # the bounds need to be in order
  imap(final_list, function(x, i) {
    full <- case_when(
      i == "lower" ~ rep(-Inf, length(formal_arguments)),
      i == "upper" ~ rep(Inf, length(formal_arguments))
    )
    names(full) <- formal_arguments

    replace <- unlist(list_flatten(x))
    full[names(full) %in% names(replace)] <- replace

    full
  })
}

# make_bounds(.func, lower = list(bottom = 0, top = 0), upper = list(top = 1))
# make_bounds(.func, lower = NULL, upper = list(top = 1))
