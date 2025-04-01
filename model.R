library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)
library(portalr)
library(portalr)

data_all <- readRDS("data_all.rds")

# Test/train split to mirror Clarke et al.
# Time 135 is equivalent to newmoon 230, which is the start of 1996
# Times 419-431 are equivalent to newmoons 514-526 which is 2019
train_inds <- which(data_all$time >= 135 & data_all$time < 419)
data_train <- lapply(seq_along(data_all), function(x) {
  if (is.matrix(data_all[[x]])) {
    data_all[[x]][train_inds, ]
  } else {
    data_all[[x]][train_inds]
  }
})

test_inds <- which(data_all$time >= 419 & data_all$time <= 431)
data_test <- lapply(seq_along(data_all), function(x) {
  if (is.matrix(data_all[[x]])) {
    data_all[[x]][test_inds, ]
  } else {
    data_all[[x]][test_inds]
  }
})

names(data_train) <- names(data_test) <- names(data_all)

priors <- get_mvgam_priors(
  formula = y ~ 1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
    te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr")),
  data = data_train,
  family = poisson(),
  trend_model = "VAR1cor"
)

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
modvar <- mvgam(
  formula = y ~ 1,
  trend_formula = ~ s(ndvi_ma12, trend, bs = "re") +
    te(mintemp, lag, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_dm, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_do, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_ot, k = c(3, 4), bs = c("tp", "cr")) +
    te(mintemp, lag, by = weights_pp, k = c(3, 4), bs = c("tp", "cr")),
  data = data_train,
  newdata = data_test,
  family = poisson(),
  trend_model = "VAR1cor",
  priors = priors,
  samples = 1600
)

saveRDS(modvar, "modvar_output.rds")

modvar <- readRDS("modvar_output.rds")

par(mfrow = c(3, 3))

plot(modvar, "forecast", series = 1)
plot(modvar, "forecast", series = 2)
plot(modvar, "forecast", series = 3)
plot(modvar, "forecast", series = 4)
plot(modvar, "forecast", series = 5)
plot(modvar, "forecast", series = 6)
plot(modvar, "forecast", series = 7)
plot(modvar, "forecast", series = 8)
plot(modvar, "forecast", series = 9)
