library(dplyr)
library(cmdstanr)
library(forecast)
library(mvgam)
library(portalr)
# source("R/get_regime.R")

data_all <- readRDS("data_pb_regime.rds")

split_train_test <- function(data_all, gap, train_start, train_end, test_start, test_end) {
  train_inds <- which(data_all$newmoonnumber >= train_start &
    data_all$newmoonnumber <= (train_end + gap + 1))
  data_train <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][train_inds, ]
    } else {
      data_all[[x]][train_inds]
    }
  })

  names(data_train) <- names(data_all)
  gap_inds <- which(data_train$newmoonnumber >= (train_end + 1) &
    data_train$newmoonnumber <= (train_end + gap))
  data_train$y[gap_inds] <- NA

  test_inds <- which(data_all$newmoonnumber >= (test_start + 1) &
    data_all$newmoonnumber <= test_end)
  data_test <- lapply(seq_along(data_all), function(x) {
    if (is.matrix(data_all[[x]])) {
      data_all[[x]][test_inds, ]
    } else {
      data_all[[x]][test_inds]
    }
  })

  names(data_test) <- names(data_all)

  return(list(train = data_train, test = data_test))
}
### Conduct within regime forecasts

within_train_starts = c(287,299,311) # starts each train set in August - same as the
                              # out- training set
without_train_start = c(336)

fit_models = function(data,train_starts,type){
  train_win_width = 60
  gap = 16
  horizon = 7 # shortened from 13 to allow 3 train sets within regime
  for (i in seq_along(train_starts)){
    train_start <- train_starts[i]
    train_stop <- train_start + 60 - 1
    test_start <- train_stop + gap + 1
    test_stop = test_start + horizon - 1
    data_split <- split_train_test(
      data_all,
      gap = gap,
      train_start = train_start,
      train_end = train_stop,
      test_start = test_start,
      test_end = test_stop
  )
  data_train <- data_split$train
  data_test <- data_split$test
  
  ## Set base priors
  
  sigma_prior <- prior(beta(10, 10), class = sigma, lb = 0.2, ub = 1)
  
  # The observation-level intercept is a nuisance parameter in this model that we
  # don't really need (and cannot adequately estimate anyway).
  # At present there is no way to drop this term using mvgam, but we can
  # heavily regularize it to zero with a strong prior
  intercept_prior <- prior(normal(0, 0.001), class = Intercept)
  
  ndvi_random_slopes_prior <- prior(
    inv_gamma(2.3693353, 0.7311319),
    class = sigma_raw_trend
  )
  
  # AR model prior in Clark et al. 2025 https://github.com/nicholasjclark/portal_VAR/blob/main/2.%20models.R
  ar_sp_intercept_prior <- prior(std_normal(), class = b)
  ar_priors <- c(sigma_prior, intercept_prior, ar_sp_intercept_prior)
  
  gam_ar_priors <- c(sigma_prior, intercept_prior, ndvi_random_slopes_prior)
  gam_var_priors <- c(sigma_prior, intercept_prior, ndvi_random_slopes_prior)
  
  
  ## Fit models
  print("start baseline model")
  baseline_model <- mvgam(
    formula = y ~ series,
    data = data_train,
    newdata = data_test,
    family = poisson(),
    priors = ar_priors,
    samples = 1600
  )
  print("start gam_var")
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
    trend_model = "VAR1",
    priors = gam_var_priors,
    samples = 1600
  )
  print("start gam_ar")
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
    priors = gam_ar_priors,
    samples = 1600
  )
  print("start ar_model")
  model_ar <- mvgam(
    formula = y ~ series,
    data = data_train,
    newdata = data_test,
    family = poisson(),
    trend_model = AR(),
    priors = ar_priors,
    samples = 1600
  )
  
  saveRDS(model_gam_var, paste("gam_var_",type,"regime_output",test_start,".rds", sep=""))
  saveRDS(model_gam_ar, paste("gam_ar_",type,"regime_output",test_start,".rds", sep=""))
  saveRDS(model_ar, paste("ar_",type,"regime_output",test_start,".rds", sep=""))
  saveRDS(baseline_model, paste("baseline_",type,"regime_output",test_start,".rds", sep=""))

  }
}

fit_models(data_all, within_train_starts, "in")


