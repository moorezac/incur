library(tidyverse)
library(minpack.lm)

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

fit

dat <- datasets::DNase
dat <- rbind(
  dat,
  data.frame(Run = c(1, 2, 3), conc = c(12.5, 12.5, 12.5), density = c(0.5, 0.75, 1))
)

library(tictoc)
{
  tic()
  model_a <- fit_model(
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
  toc()
}
model_a$fit

# more data points = smoother curve
# this is easier than a geom_function approach
x_vals <- seq(min(model_a$data$x), max(model_a$data$x), length.out = 1e4)
# create the data and return
predicted_a <- tibble(
  x = x_vals,
  y = predict(model_a$fit, newdata = tibble(x = x_vals))
)
# plot
ggplot(mapping = aes(x, y)) +
  geom_jitter(data = model_a$data) +
  geom_line(data = predicted_a) + 
  ggtitle("least squares")

model_b <- fit_model(
  .data = dat,
  .x_var = conc,
  .y_var = density,
  .curve_func = func,
  .start_func = func_start,
  .huber = TRUE,
  .lower_bounds = list(y0 = 0),
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)
# more data points = smoother curve
# this is easier than a geom_function approach
x_vals <- seq(min(model_b$data$x), max(model_b$data$x), length.out = 1e4)
# create the data and return
predicted_b <- tibble(
  x = x_vals,
  y = predict(model_b$fit, newdata = tibble(x = x_vals))
)
# plot
ggplot(mapping = aes(x, y)) +
  geom_jitter(data = model_b$data) +
  geom_line(data = predicted_b) + 
  ggtitle("iterative reweighted least squares")
weights(model_b$fit)

model_c <- fit_model(
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
x_vals <- seq(min(model_c$data$x), max(model_c$data$x), length.out = 1e4)
# create the data and return
predicted_c <- tibble(
  x = x_vals,
  y = predict(model_c$fit, newdata = tibble(x = x_vals))
)
# plot
ggplot(mapping = aes(x, y)) +
  geom_jitter(aes(shape = outlier_density), data = model_c$data) +
  geom_line(data = predicted_c) + 
  ggtitle("irwls + outlier removal")
