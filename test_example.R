func <- function(x, y0, k) y0 * exp(k * x)
func_start <- function(x, y) list(y0 = min(y), k = 0)



func_start(x = Puromycin$conc, y = Puromycin$rate)

plot(Puromycin$rate ~ Puromycin$conc, col = Puromycin$state)


library(minpack.lm)
# fit the curve


res$data
res$data$outlier_density |> table()


ggplot() +
  geom_jitter(aes(x, y), res$data) +
  geom_line(aes(x, y), tibble(x = res$data |> pull(x), y = predict(res$fit)))

# exponential plateau
func <- function(x, y0, ym, k) ym - (ym - y0) * exp(-k * x)
func_start <- function(x, y) list(ym = max(y), y0 = min(y), k = 1)

# fit the curve to the data
fit <- fit_model(
  .data = Puromycin,
  .x_var = conc,
  .y_var = rate,
  .curve_func = func,
  .start_func = func_start,
  .detect_outliers = TRUE,
  .lower_bounds = list(y0 = 50),
  .shared_group = state,
  .shared_params = c("k", "y0"),
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)

# create data from each group
predicted <- map(unique(Puromycin$state), function(i) {
  # more data points = smoother curve
  # this is easier than a geom_function approach
  x_vals <- seq(min(fit$data$x), max(fit$data$x), length.out = 1e4)
  # create the data and return
  tibble(
    x = x_vals,
    y = predict(fit$fit, newdata = tibble(x = x_vals, group = i)),
    group = i
  )
})
# collate predicted data
predicted <- bind_rows(predicted)

# plot
ggplot(mapping = aes(x, y, colour = group)) +
  geom_point(data = fit$data) +
  geom_line(data = predicted)
