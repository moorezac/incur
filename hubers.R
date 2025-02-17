func <- function(x, y0, ym, k) ym - (ym - y0) * exp(-k * x)
func_start <- function(x, y) list(ym = max(y), y0 = min(y), k = 1)

fit <- fit_model(
  .data = Puromycin,
  .x_var = conc,
  .y_var = rate,
  .curve_func = func,
  .start_func = func_start,
  .detect_outliers = TRUE,
  .lower_bounds = list(y0 = 0),
  .shared_group = state,
  .shared_params = c("y0", "k"),
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
# collate
predicted <- bind_rows(predicted)
# plot
ggplot(mapping = aes(x, y, colour = group)) +
  geom_point(data = fit$data) +
  geom_line(data = predicted)


dat <- datasets::DNase

dat <- rbind(
  dat,
  data.frame(Run = c(1, 2, 3), conc = c(12.5, 12.5, 12.5), density = c(0.5, 0.75, 1))
)

plot(dat$density ~ dat$conc)

fit <- fit_model(
  .data = dat,
  .x_var = conc,
  .y_var = density,
  .curve_func = func,
  .start_func = func_start,
  .huber = FALSE,
  .detect_outliers = FALSE,
  .lower_bounds = list(y0 = 0),
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)
# more data points = smoother curve
# this is easier than a geom_function approach
x_vals <- seq(min(fit$data$x), max(fit$data$x), length.out = 1e4)
# create the data and return
predicted <- tibble(
  x = x_vals,
  y = predict(fit$fit, newdata = tibble(x = x_vals))
)
# plot
ggplot(mapping = aes(x, y)) +
  geom_jitter(data = fit$data) +
  geom_line(data = predicted)

fit_a <- fit_model(
  .data = dat,
  .x_var = conc,
  .y_var = density,
  .curve_func = func,
  .start_func = func_start,
  .huber = TRUE,
  .detect_outliers = FALSE,
  .lower_bounds = list(y0 = 0),
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)
# more data points = smoother curve
# this is easier than a geom_function approach
x_vals <- seq(min(fit_a$data$x), max(fit_a$data$x), length.out = 1e4)
# create the data and return
predicted <- tibble(
  x = x_vals,
  y = predict(fit_a$fit, newdata = tibble(x = x_vals))
)
# plot
ggplot(mapping = aes(x, y)) +
  geom_jitter(data = fit$data) +
  geom_line(data = predicted)
weights(fit_a$fit)

fit_b <- fit_model(
  .data = dat,
  .x_var = conc,
  .y_var = density,
  .curve_func = func,
  .start_func = func_start,
  .huber = TRUE,
  .detect_outliers = TRUE,
  .lower_bounds = list(y0 = 0),
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)
# more data points = smoother curve
# this is easier than a geom_function approach
x_vals <- seq(min(fit_b$data$x), max(fit_b$data$x), length.out = 1e4)
# create the data and return
predicted <- tibble(
  x = x_vals,
  y = predict(fit_b$fit, newdata = tibble(x = x_vals))
)
# plot
ggplot(mapping = aes(x, y)) +
  geom_jitter(aes(shape = outlier_density), data = fit_b$data) +
  geom_line(data = predicted)
