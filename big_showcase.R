library(rlang)
library(tidyverse)
library(minpack.lm)

roxygen2::roxygenise()

# incur needs a function that describes a curve:

# exponential plateau
curve_func <- function(x, y0, ym, k) {
  ym - (ym - y0) * exp(-k * x)
}
# y0 = minimum
# ym = plateau
# k = rate of growth
curve_func(x = 0:10, y0 = 0, k = 1, ym = 1) |> plot()
curve_func(x = 0:10, y0 = 0, k = 0.5, ym = 1) |> plot()

# and a function that gives me starting values
start_func <- function(x, y) {
  list(
    ym = max(y),
    y0 = min(y),
    k = 1
  )
}

source("R/models.R")

Puromycin |> View()

ggplot(Puromycin, aes(conc, rate, colour = state)) + 
  geom_point()

initial_fit <- incur::fit_model(
  .data = Puromycin,
  .x_var = conc,
  .y_var = rate,
  .curve_func = curve_func,
  .start_func = start_func,
  .lower_bounds = list(y0 = 0),
  # dots dots dots
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)

# more data points = smoother curve
# this is easier than a geom_function approach
x_vals <- seq(min(initial_fit$data$x), max(initial_fit$data$x), length.out = 1e4)
# create the data and return
predicted <- tibble(
  x = x_vals,
  y = predict(initial_fit$fit, newdata = tibble(x = x_vals))
)
# plot
ggplot(mapping = aes(x, y)) +
  geom_point(data = initial_fit$data, mapping = aes(x, y, colour = state)) +
  geom_line(data = predicted)

Puromycin

shared_fit <- incur::fit_model(
  .data = Puromycin,
  .x_var = conc,
  .y_var = rate,
  .curve_func = curve_func,
  .start_func = start_func,
  .lower_bounds = list(y0 = 0),
  # shared coefs:
  .shared_group = state,
  .shared_params = c("y0", "k"),
  .return_func = TRUE,
  control = minpack.lm::nls.lm.control(maxiter = 1e3)
)

# look at the coefs!!
initial_fit$fit
shared_fit$fit

# the function modifications:
curve_func
shared_fit$func

# this is generalisable!

# generics work
class(shared_fit$fit)

coef(shared_fit$fit)
residuals(shared_fit$fit)
predict(shared_fit$fit)
logLik(shared_fit$fit)

# create data from each group
predicted <- map(unique(Puromycin$state), function(i) {
  # more data points = smoother curve
  # this is easier than a geom_function approach
  x_vals <- seq(min(shared_fit$data$x), max(shared_fit$data$x), length.out = 1e4)
  # create the data and return
  tibble(
    x = x_vals,
    y = predict(shared_fit$fit, newdata = tibble(x = x_vals, group = i)),
    group = i
  )
})
# collate
predicted <- bind_rows(predicted)
# plot
ggplot(mapping = aes(x, y, colour = group)) +
  geom_point(data = shared_fit$data) +
  geom_line(data = predicted)

# outliers
dat <- datasets::DNase
dat

ggplot(dat, aes(conc, density)) + 
  geom_jitter()

# add in
dat <- rbind(
  dat,
  data.frame(Run = c(1, 2, 3), conc = c(12.5, 12.5, 12.5), density = c(0.5, 0.75, 1))
)

ggplot(dat, aes(conc, density)) + 
  geom_jitter()

library(tictoc)
{
  tic()
  model_a <- incur::fit_model(
    .data = dat,
    .x_var = conc,
    .y_var = density,
    .curve_func = curve_func,
    .start_func = start_func,
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

# add in hubers
model_b <- incur::fit_model(
  .data = dat,
  .x_var = conc,
  .y_var = density,
  .curve_func = curve_func,
  .start_func = start_func,
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

# outlier detection via ROUT
model_c <- incur::fit_model(
  .data = dat,
  .x_var = conc,
  .y_var = density,
  .curve_func = curve_func,
  .start_func = start_func,
  .huber = TRUE,
  .detect_outliers = TRUE,
  .lower_bounds = list(y0 = 0)
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

# how much does it change fitted coefs?
model_b$fit
model_c$fit

plot_model(
  .data = shared_fit$data,
  .fit = shared_fit$fit,
  .x_var = conc,
  .y_var = rate
)
plot_model(
  .data = model_c$data,
  .fit = model_c$fit,
  .x_var = conc,
  .y_var = density
)

incur::find_root_fit(
  .fit = model_c$fit,
  .x_vals = model_c$data$x,
  .target = 1
)
incur::calculate_auc(
  .fit = model_c$fit,
  .lower = min(model_c$data$x),
  .upper = max(model_c$data$x)
)
