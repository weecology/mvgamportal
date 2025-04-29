library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)
library(portalr)

data_all <- readRDS("data_pb_regime.rds")

split_train_test <- function(data_all, train_start, train_end, test_start, test_end) {
  train_inds <- which(data_all$newmoonnumber >= train_start &
    data_all$newmoonnumber <= train_end)
  data_train <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][train_inds, ]
    } else {
      data_all[[x]][train_inds]
    }
  })

  test_inds <- which(data_all$newmoonnumber >= test_start &
    data_all$newmoonnumber <= test_end)
  data_test <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][test_inds, ]
    } else {
      data_all[[x]][test_inds]
    }
  })

  names(data_train) <- names(data_test) <- names(data_all)

  return(list(train = data_train, test = data_test))
}

data_split <- split_train_test(data_all, train_start = 279, train_end = 394, test_start = 395, test_end = 395 + 12)
data_train <- data_split$train
data_test <- data_split$test

# Base priors

priors <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)

# Update the prior for the NDVI random slopes
priors <- c(
  priors,
  prior(inv_gamma(2.3693353, 0.7311319), class = sigma_raw_trend)
)

# The observation-level intercept is a nuisance parameter in this model that we
# don't really need (and cannot adequately estimate anyway).
# At present there is no way to drop this term using mvgam, but we can
# heavily regularize it to zero with a strong prior
priors <- c(
  priors,
  prior(normal(0, 0.001), class = Intercept)
)

# Fit the model
model_gam_var <- mvgam(
  formula = y ~ -1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
    te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr")),
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = "VAR1cor",
  priors = priors,
  samples = 1600
)

model_gam_ar <- mvgam(
  formula = y ~ -1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
    te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr")),
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  priors = priors,
  samples = 1600
)

# Prior adjustment from AR model in Clark et al. 2025
# https://github.com/nicholasjclark/portal_VAR/blob/main/2.%20models.R
priors <- c(priors, prior(std_normal(), class = b))

model_ar <- mvgam(
  formula = y ~ series,
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = AR(),
  priors = priors,
  samples = 1600
)

saveRDS(model_gam_var, "gam_var_one_year_into_transition_output_y_minus_1.rds")
saveRDS(model_gam_ar, "gam_ar_one_year_into_transition_output_y_minus_1.rds")
saveRDS(model_ar, "ar_one_year_into_transition_output_y_minus_1.rds")

model_gam_var <- readRDS("gam_var_one_year_into_transition_output_y_minus_1.rds")
model_gam_ar <- readRDS("gam_ar_one_year_into_transition_output_y_minus_1.rds")
model_ar <- readRDS("ar_one_year_into_transition_output_y_minus_1.rds")

par(mfrow = c(3, 3))

plot(model_gam_var, "forecast", series = 1)
plot(model_gam_var, "forecast", series = 2)
plot(model_gam_var, "forecast", series = 3)
plot(model_gam_var, "forecast", series = 4)
plot(model_gam_var, "forecast", series = 5)
plot(model_gam_var, "forecast", series = 6)
plot(model_gam_var, "forecast", series = 7)


dev.print(pdf, "gam_var_one_year_into_transition_y_minus_1.pdf")

par(mfrow = c(3, 3))

plot(model_gam_ar, "forecast", series = 1)
plot(model_gam_ar, "forecast", series = 2)
plot(model_gam_ar, "forecast", series = 3)
plot(model_gam_ar, "forecast", series = 4)
plot(model_gam_ar, "forecast", series = 5)
plot(model_gam_ar, "forecast", series = 6)
plot(model_gam_ar, "forecast", series = 7)

dev.print(pdf, "gam_ar_one_year_into_transition.pdf")

par(mfrow = c(3, 3))

plot(model_ar, "forecast", series = 1)
plot(model_ar, "forecast", series = 2)
plot(model_ar, "forecast", series = 3)
plot(model_ar, "forecast", series = 4)
plot(model_ar, "forecast", series = 5)
plot(model_ar, "forecast", series = 6)
plot(model_ar, "forecast", series = 7)

dev.print(pdf, "ar_one_year_into_transition_y_just_series.pdf")
